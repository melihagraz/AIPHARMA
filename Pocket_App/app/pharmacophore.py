"""3D pharmacophore extraction from a bound ligand.

Strategy: parse the bound ligand from a co-crystal PDB, build an RDKit
molecule with explicit 3D coordinates, then run RDKit's BaseFeatures
fdef on it. Each feature carries a 3D position (averaged across the
matching atoms) and a type — Donor, Acceptor, Aromatic, Hydrophobe,
PosIonizable, NegIonizable.

The returned pharmacophore is the "ground truth" pattern for a known
EGFR/JAK2/etc. inhibitor and is used downstream to:

  - filter REINVENT-generated molecules (must hit ≥ 4 of N features)
  - score candidates by feature-overlap quality
  - guide RL fine-tuning (Phase 4)

We deliberately use RDKit's stock ``BaseFeatures.fdef`` (ships with
RDKit, MIT-equivalent licence). No proprietary feature definitions.
"""

from __future__ import annotations

import logging
from io import StringIO
from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures
from rdkit.Geometry import Point3D
from rdkit import RDConfig

from .pocket import _fetch_pdb, _iter_hetatm_records, NON_LIGAND_RESNAMES, list_bundled_targets
from .schemas import PharmacophoreFeature

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Feature definition factory (cached)
# ---------------------------------------------------------------------------

_FEATURE_FACTORY: Optional[ChemicalFeatures.MolChemicalFeatureFactory] = None


def _feature_factory() -> ChemicalFeatures.MolChemicalFeatureFactory:
    """Lazily build the BaseFeatures factory shipped with RDKit."""
    global _FEATURE_FACTORY
    if _FEATURE_FACTORY is None:
        fdef_path = Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"
        if not fdef_path.is_file():
            raise RuntimeError(
                f"RDKit BaseFeatures.fdef not found at {fdef_path}. "
                "RDKit installation may be incomplete."
            )
        _FEATURE_FACTORY = ChemicalFeatures.BuildFeatureFactory(str(fdef_path))
    return _FEATURE_FACTORY


# Map RDKit feature family names to the typed enum we expose on the wire.
_FEATURE_TYPE_MAP = {
    "Donor": "Donor",
    "Acceptor": "Acceptor",
    "Aromatic": "Aromatic",
    "Hydrophobe": "Hydrophobe",
    "LumpedHydrophobe": "Hydrophobe",
    "PosIonizable": "PosIonizable",
    "NegIonizable": "NegIonizable",
}


# ---------------------------------------------------------------------------
# Bundled PDB lookup — for a target name, find a co-crystal PDB we ship
# ---------------------------------------------------------------------------

_BUNDLED_COCRYSTALS: dict[str, tuple[str, Optional[str]]] = {
    # target → (pdb_id, ligand_resname). Ligand resname optional; if None
    # the largest non-trivial HETATM is used.
    "EGFR":  ("1m17", "AQ4"),  # 1M17 has erlotinib (AQ4)
    "JAK2":  ("3lpb", None),   # JAK2 with inhibitor
    "BACE1": ("3ind", None),
    "KDR":   ("3vhe", None),   # VEGFR2
    "MET":   ("3lq8", None),
    "AKT1":  ("3o96", None),
    "LCK":   ("2of4", None),
    "KIT":   ("4u0i", None),
}


def _bundled_pdb_for_target(target: str) -> Optional[tuple[str, Optional[str]]]:
    if not target:
        return None
    return _BUNDLED_COCRYSTALS.get(target.strip().upper())


# ---------------------------------------------------------------------------
# Ligand → RDKit Mol with 3D coords
# ---------------------------------------------------------------------------


def _extract_ligand_pdb_block(
    pdb_text: str, ligand_resname: str
) -> str:
    """Return a PDB block containing only the requested HETATM ligand."""
    target = ligand_resname.strip().upper()
    out: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith("HETATM") and len(line) >= 20:
            if line[17:20].strip().upper() == target:
                out.append(line)
        elif line.startswith("CONECT"):
            # CONECT records don't reference resname directly. Keep them all
            # — RDKit will ignore atoms it doesn't have.
            out.append(line)
    if not out:
        raise ValueError(f"Ligand {ligand_resname!r} not found in PDB.")
    return "\n".join(out) + "\nEND\n"


# ---------------------------------------------------------------------------
# RCSB Chemical Component Dictionary — SMILES template lookup
# ---------------------------------------------------------------------------
# PDB files carry atom coords + CONECT connectivity but NO bond orders. The
# robust way to recover bond orders for a HETATM ligand is to grab the
# correct Kekulé SMILES from RCSB's Chemical Component Dictionary (CCD) and
# overlay it onto the parsed mol via ``AssignBondOrdersFromTemplate``. The
# CCD is RCSB's curated database of every HETATM 3-letter code (AQ4 =
# erlotinib, YY3 = osimertinib, 66T = AZD4547, etc.), exposed as a free
# REST endpoint with no auth.
#
# Cached in-memory so we don't hammer RCSB across requests.
_CCD_SMILES_CACHE: dict[str, Optional[str]] = {}
_CCD_REST_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{resname}"


def _fetch_ccd_smiles(resname: str) -> Optional[str]:
    """Fetch the canonical SMILES for a HETATM residue from RCSB CCD.

    Returns ``None`` on network failure / missing record so the caller can
    fall back to perception-only parsing. We prefer the OpenEye descriptor
    when present (matches RDKit's atom indexing better) and fall back to
    the CACTVS or InChI-derived descriptors.
    """
    if not resname:
        return None
    key = resname.strip().upper()
    if key in _CCD_SMILES_CACHE:
        return _CCD_SMILES_CACHE[key]
    try:
        import urllib.request
        import json as _json
        url = _CCD_REST_URL.format(resname=key)
        with urllib.request.urlopen(url, timeout=5.0) as resp:  # noqa: S310
            payload = _json.loads(resp.read().decode("utf-8"))
    except Exception as exc:  # noqa: BLE001
        logger.info("RCSB CCD lookup failed for %r: %s", key, exc)
        _CCD_SMILES_CACHE[key] = None
        return None
    # RCSB shape: payload["pdbx_chem_comp_descriptor"] is a LIST of
    # {"descriptor": "...", "type": "SMILES"|"InChI"|"SMILES_CANONICAL", "program": "..."}.
    # Prefer OpenEye OEToolkits (best aromatic perception), then anything
    # labelled SMILES. Skip InChI strings.
    descriptors = payload.get("pdbx_chem_comp_descriptor") or []
    if not isinstance(descriptors, list):
        descriptors = []

    def _pick(filter_fn) -> Optional[str]:
        for d in descriptors:
            if not isinstance(d, dict):
                continue
            if filter_fn(d):
                s = d.get("descriptor")
                if isinstance(s, str) and s:
                    return s
        return None

    smiles = (
        _pick(lambda d: d.get("type") in ("SMILES_CANONICAL", "SMILES") and "OpenEye" in str(d.get("program") or ""))
        or _pick(lambda d: d.get("type") == "SMILES_CANONICAL")
        or _pick(lambda d: d.get("type") == "SMILES")
    )
    _CCD_SMILES_CACHE[key] = smiles
    return smiles


def _build_ligand_mol(pdb_block: str, resname: Optional[str] = None) -> Chem.Mol:
    """Parse a PDB block into an RDKit Mol with proper bond orders.

    PDB ligands lack bond orders, so naive parsing yields a fully-single-
    bonded skeleton with no aromaticity — every aromatic-ring feature in
    RDKit's BaseFeatures.fdef misses. (This is why YY3 / osimertinib was
    reporting 0 Aromatic features despite having three aromatic rings:
    indole, pyrimidine, anilino phenyl.)

    Strategy:
      1. Parse PDB with ``sanitize=False`` → raw 3D connectivity.
      2. If we know the HETATM resname (e.g. "YY3"), fetch its canonical
         SMILES from RCSB's Chemical Component Dictionary and overlay it
         with ``AssignBondOrdersFromTemplate`` — this is the standard
         recipe used by every serious cheminformatics pipeline.
      3. Full ``SanitizeMol`` to perceive aromaticity now that bond orders
         are correct.
      4. Lenient fallback chain so unusual ligands still yield *something*.

    Earlier attempt used ``rdDetermineBonds.DetermineBonds`` which infers
    bond orders from coordinates alone — but it fails on most drug-sized
    molecules (charge mismatch) and degrades feature counts. Removed.
    """
    mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError("RDKit could not parse the ligand PDB block.")

    # Some PDB entries pack multiple copies of the ligand into one HETATM
    # residue (e.g. 4V05's "66T" is a dimer of AZD4547: 67 atoms in the
    # crystal asymmetric unit, but the CCD template is a single 34-atom
    # molecule). AssignBondOrdersFromTemplate requires the atom count to
    # match — for a dimer it would silently fail and we'd lose aromaticity
    # perception. Split the parsed mol into disconnected fragments and
    # pick the *largest* one for downstream pharmacophore extraction.
    # A pharmacophore is a property of the molecule, not crystal packing,
    # so one representative copy is the scientifically correct choice.
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    except Exception:  # noqa: BLE001
        frags = (mol,)
    if len(frags) > 1:
        # Pick the heaviest fragment (most heavy atoms). Falls back to the
        # first if heavy-atom counts tie.
        frags_sorted = sorted(frags, key=lambda m: m.GetNumHeavyAtoms(), reverse=True)
        mol = frags_sorted[0]
        logger.info(
            "Ligand %r split into %d fragments; using the largest "
            "(%d heavy atoms) for pharmacophore extraction.",
            resname, len(frags), mol.GetNumHeavyAtoms(),
        )

    # Try to assign bond orders from the CCD template SMILES.
    if resname:
        template_smiles = _fetch_ccd_smiles(resname)
        if template_smiles:
            try:
                from rdkit.Chem import AllChem
                template = Chem.MolFromSmiles(template_smiles)
                if template is not None and template.GetNumAtoms() > 0:
                    mol_with_orders = AllChem.AssignBondOrdersFromTemplate(template, mol)
                    if mol_with_orders is not None:
                        mol = mol_with_orders
            except Exception as exc:  # noqa: BLE001
                logger.info(
                    "AssignBondOrdersFromTemplate failed for %r (%s); "
                    "falling back to perception-only parsing.",
                    resname, exc,
                )

    # Full sanitise — this is what actually triggers aromaticity perception
    # against the now-correct Kekulé bond orders.
    try:
        Chem.SanitizeMol(mol)
    except (Chem.AtomValenceException, Chem.KekulizeException):
        # Lenient round — drop kekulisation + property normalisation but
        # keep aromaticity perception itself.
        try:
            Chem.SanitizeMol(
                mol,
                sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                & ~Chem.SanitizeFlags.SANITIZE_KEKULIZE
                & ~Chem.SanitizeFlags.SANITIZE_PROPERTIES,
            )
        except Exception:  # noqa: BLE001
            logger.warning("Lenient sanitize also failed; using mol as-is")
    return mol


# ---------------------------------------------------------------------------
# Pharmacophore feature extraction
# ---------------------------------------------------------------------------


def _features_from_mol(mol: Chem.Mol) -> list[PharmacophoreFeature]:
    """Extract pharmacophore features with their 3D positions.

    RDKit's ``BaseFeatures.fdef`` ships multiple SMARTS patterns per family
    (e.g. ``SingleAtomDonor``, ``ChainTwoWayAttach`` and ``H2N_NoCharge`` all
    match a primary amine). For a typical drug-sized ligand this produces
    redundant features: the same atom set, matched by the same family,
    appears 2–4 times. Pharma DD expects 4–8 donors for V6G-sized molecules,
    not 11.

    We dedupe by ``(family, atom_set)``: identical atoms in the same family
    collapse to one feature. Cross-family overlap is preserved on purpose —
    a hydroxyl really IS both donor and acceptor, and the docking layer
    needs that to reason about H-bond geometry.
    """
    factory = _feature_factory()
    raw = factory.GetFeaturesForMol(mol)
    out: list[PharmacophoreFeature] = []
    conf = mol.GetConformer() if mol.GetNumConformers() else None
    # (ftype, sorted_atom_ids) → already emitted? If so, skip.
    seen: set[tuple] = set()
    for f in raw:
        family = f.GetFamily()
        ftype = _FEATURE_TYPE_MAP.get(family)
        if ftype is None:
            continue
        atom_ids_tuple = tuple(sorted(f.GetAtomIds()))
        dedup_key = (ftype, atom_ids_tuple)
        if dedup_key in seen:
            continue
        seen.add(dedup_key)
        # RDKit's GetPos() gives the centroid of the matching atom set
        try:
            pos: Point3D = f.GetPos()
        except RuntimeError:
            # Some features may fail to compute pos if the conformer is missing
            if conf is None:
                continue
            if not atom_ids_tuple:
                continue
            xs, ys, zs = 0.0, 0.0, 0.0
            for aid in atom_ids_tuple:
                p = conf.GetAtomPosition(aid)
                xs += p.x
                ys += p.y
                zs += p.z
            n = len(atom_ids_tuple)
            pos = Point3D(xs / n, ys / n, zs / n)
        out.append(
            PharmacophoreFeature(
                type=ftype,  # type: ignore[arg-type]
                x=round(pos.x, 3),
                y=round(pos.y, 3),
                z=round(pos.z, 3),
                tolerance=1.5,
                label=f.GetType(),
            )
        )
    return out


# ---------------------------------------------------------------------------
# High-level entrypoints
# ---------------------------------------------------------------------------


def pharmacophore_from_pdb(
    pdb_id: str,
    ligand_resname: Optional[str] = None,
    prepared_pdb_path: Optional[str] = None,
) -> tuple[list[PharmacophoreFeature], Optional[str], Optional[str]]:
    """Extract a pharmacophore from a co-crystal PDB.

    Returns ``(features, chosen_ligand_resname, ligand_smiles)``.

    When ``prepared_pdb_path`` is provided, the file is read directly and the
    RCSB fetch is skipped. The caller is responsible for upstream preparation
    (PDBFixer cleanup, protonation, water handling).
    """
    from pathlib import Path

    if prepared_pdb_path:
        prep_path = Path(prepared_pdb_path)
        if not prep_path.is_file():
            raise ValueError(
                f"prepared_pdb_path does not exist or is not a file: {prepared_pdb_path}"
            )
        pdb_path = prep_path
    else:
        pdb_path = _fetch_pdb(pdb_id)
    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")

    # Pick the ligand to use
    chosen: Optional[str] = None
    if ligand_resname:
        chosen = ligand_resname.strip().upper()
    else:
        candidates = list(_iter_hetatm_records(pdb_text))
        best_size = 0
        for resname, coords in candidates:
            if resname.upper() in NON_LIGAND_RESNAMES:
                continue
            if len(coords) > best_size and len(coords) >= 6:
                best_size = len(coords)
                chosen = resname
        if chosen is None:
            raise ValueError(
                f"No drug-like HETATM in PDB {pdb_id}. Provide a co-crystal PDB."
            )

    pdb_block = _extract_ligand_pdb_block(pdb_text, chosen)
    # Pass the HETATM resname so _build_ligand_mol can look up the correct
    # Kekulé SMILES in RCSB's Chemical Component Dictionary and assign bond
    # orders — required for aromaticity perception on aromatic ligands like
    # YY3 (osimertinib), AQ4 (erlotinib), V6G (LY3502970), etc.
    mol = _build_ligand_mol(pdb_block, resname=chosen)
    features = _features_from_mol(mol)

    smiles: Optional[str] = None
    try:
        smiles = Chem.MolToSmiles(mol)
    except Exception:  # noqa: BLE001
        smiles = None

    return features, chosen, smiles


def pharmacophore_from_target(
    target: str,
) -> tuple[list[PharmacophoreFeature], Optional[str], Optional[str], Optional[str]]:
    """Extract a pharmacophore for a bundled target via its known co-crystal.

    Returns ``(features, pdb_id, ligand_resname, ligand_smiles)``.
    """
    spec = _bundled_pdb_for_target(target)
    if spec is None:
        raise ValueError(
            f"No bundled co-crystal PDB known for target {target!r}. "
            f"Provide pdb_id explicitly. Available: {list(_BUNDLED_COCRYSTALS)}."
        )
    pdb_id, ligand_resname = spec
    feats, chosen, smiles = pharmacophore_from_pdb(pdb_id, ligand_resname)
    return feats, pdb_id, chosen, smiles
