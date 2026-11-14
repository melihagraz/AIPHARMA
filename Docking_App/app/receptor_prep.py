"""On-the-fly receptor PDB → PDBQT preparation using meeko.

When the De Novo pipeline runs against a target that has *no* pre-bundled
PDBQT under ``receptors/`` — typically because the user uploaded their own
PDB via ``/api/pocket/upload-pdb`` for a non-standard target, allosteric
pocket, or mutant — we have to prepare a receptor PDBQT on the fly before
Vina can do anything with it.

Pipeline:
    PDB text → meeko Polymer.from_pdb_string (handles altlocs, missing
    residues with allow_bad_res, charge assignment) → PDBQTWriterLegacy
    → cached on disk keyed by content hash so repeat runs against the
    same uploaded PDB don't pay the (~10-20 s) prep cost again.

Failures are *not* silent: callers receive a clear error message instead
of a corrupt PDBQT, so the caller can fall back to the bundled receptor
or surface the issue to the user.
"""

from __future__ import annotations

import hashlib
import logging
import warnings
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Cache prepared receptors next to the bundled ones so the existing
# resolve_target logic can find them without changes.
_RECEPTORS_DIR = Path(__file__).resolve().parent.parent / "receptors"
_CUSTOM_CACHE_DIR = _RECEPTORS_DIR / "custom_cache"


def _cache_key(pdb_content: str) -> str:
    """SHA-256 of the canonicalised PDB content (whitespace normalised)."""
    canonical = "\n".join(line.rstrip() for line in pdb_content.splitlines() if line.strip())
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:16]


def prepare_custom_receptor(
    pdb_content: str,
    *,
    pdb_id: Optional[str] = None,
    default_altloc: str = "A",
    allow_bad_res: bool = True,
) -> tuple[Path, str]:
    """Prepare a receptor PDBQT from raw PDB content.

    Returns ``(pdbqt_path, cache_key)``. ``pdbqt_path`` points at the cached
    .pdbqt file inside ``Docking_App/receptors/custom_cache/`` — Vina's
    ``--receptor`` flag can read it directly. ``cache_key`` is the
    16-char hash so callers can advertise it in logs / debug surfaces.

    The ``pdb_id`` hint (when supplied) is used purely for the filename
    suffix so the cached file is human-recognisable on disk. It does not
    affect the cache key — that's content-derived so identical PDB
    content under different IDs still de-duplicates.

    Raises ``ValueError`` when meeko can't build a valid Polymer from the
    PDB content; the message includes the meeko log so the caller can
    surface it to the user.
    """
    if not pdb_content or not pdb_content.strip():
        raise ValueError("Empty PDB content")
    key = _cache_key(pdb_content)
    label = (pdb_id or "custom").upper()[:8]
    _CUSTOM_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cached = _CUSTOM_CACHE_DIR / f"{label}_{key}.pdbqt"

    if cached.exists() and cached.stat().st_size > 0:
        logger.info("Custom receptor cache hit: %s (%d bytes)", cached.name, cached.stat().st_size)
        return cached, key

    # Lazy meeko import — only paid when we actually need to prep a
    # custom receptor, not on every Docking_App startup.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from meeko import (
            MoleculePreparation,
            PDBQTWriterLegacy,
            Polymer,
            ResidueChemTemplates,
        )

        templates = ResidueChemTemplates.create_from_defaults()
        mk_prep = MoleculePreparation()
        try:
            poly = Polymer.from_pdb_string(
                pdb_content,
                chem_templates=templates,
                mk_prep=mk_prep,
                set_template=None,
                residues_to_delete=None,
                allow_bad_res=allow_bad_res,
                default_altloc=default_altloc,
            )
        except Exception as exc:  # noqa: BLE001
            raise ValueError(
                f"meeko Polymer build failed for {label}: {exc.__class__.__name__}: {exc}"
            ) from exc

        if poly is None or not getattr(poly, "monomers", None):
            raise ValueError(
                f"meeko returned an empty Polymer for {label} — no usable monomers"
            )

        result = PDBQTWriterLegacy.write_string_from_polymer(poly)
        if isinstance(result, tuple):
            pdbqt_text, meeko_log = result[0], result[1] if len(result) > 1 else ""
        else:
            pdbqt_text, meeko_log = result, ""

        if not pdbqt_text or len(pdbqt_text.strip()) < 100:
            raise ValueError(
                f"meeko PDBQT output for {label} is empty / too short. "
                f"meeko log: {str(meeko_log)[:300]}"
            )

    cached.write_text(pdbqt_text, encoding="utf-8")
    logger.info(
        "Prepared custom receptor: %s (%d bytes, %d monomers)",
        cached.name, cached.stat().st_size, len(poly.monomers),
    )
    return cached, key
