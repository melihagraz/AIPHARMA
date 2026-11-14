"""Extended receptor catalog for on-demand RCSB → meeko → Vina docking.

The 8 entries under Docking_App/receptors/ are *bundled* — their PDBQT
files and Vina conf.txt grid boxes ship with the install. Everything in
this file is *on-demand*: the user picks a name, the backend fetches
the canonical PDB from RCSB, runs Meeko to produce a PDBQT, and feeds
the rest of the pipeline. Pocket geometry comes from Pocket_App (the
chemist runs Step 1 on the auto-fetched PDB) — we don't ship a
pre-computed grid for these because that would require curating ~60
co-crystal ligand centroids.

Each entry pins:
    pdb       — canonical RCSB id (Homo sapiens, ideally co-crystal with
                a reference inhibitor at the orthosteric/ATP site)
    gene      — HGNC gene symbol; used by rcsb_search title-token
                filtering and the ChEMBL actives fetch
    family    — coarse class used by the UI to group dropdown options;
                see ``FAMILIES`` below
    aliases   — alternative names the user might type; matched
                case-insensitively in the lookup endpoint

PDB choices: when possible the canonical inhibitor-bound crystal
(e.g. 1M17 EGFR + erlotinib, 2RH1 ADRB2 + carazolol, 1T46 KIT +
imatinib). For receptors without a famous co-crystal we pick the
highest-resolution Homo sapiens entry surfaced by the gene-name
strict RCSB filter — chemists routinely re-pick a different PDB
per project, so this is a starting point, not a commitment.

Verified against RCSB Search API v2 on 2026-05-31: every gene name
below returns ≥1 Homo sapiens hit. Two entries (VEGFR3, N-Cadherin)
needed a UniProt-id fallback because their gene_name annotations are
inconsistent across RCSB entries — see comments.
"""

from __future__ import annotations

from dataclasses import dataclass


FAMILIES = (
    "rtk",           # receptor tyrosine kinases
    "hormone",       # peptide / metabolic hormone receptors
    "cytokine",      # interleukin / interferon / EPO receptors
    "checkpoint",    # immune checkpoint receptors
    "tlr",           # toll-like receptors
    "gpcr_amine",    # monoamine + chemokine GPCRs
    "lgic",          # ligand-gated ion channels
    "adhesion",      # integrins + cadherins
    "death",         # death receptors
    "antigen",       # TCR + BCR subunits
)


@dataclass(frozen=True)
class ReceptorEntry:
    name: str          # canonical display name
    pdb: str           # 4-char RCSB id
    gene: str          # HGNC gene symbol
    family: str        # one of FAMILIES
    aliases: tuple[str, ...] = ()  # other names users might type


# Ordered for stable JSON output and predictable UI rendering. Insertion
# order is intentionally by family so the frontend dropdown groups
# without an extra sort step.
EXTENDED_CATALOG: tuple[ReceptorEntry, ...] = (
    # ---- Receptor Tyrosine Kinases ----
    ReceptorEntry("EGFR",   "1M17", "EGFR",  "rtk", ("HER1", "ERBB1")),
    ReceptorEntry("HER2",   "3PP0", "ERBB2", "rtk", ("ERBB2", "NEU")),
    ReceptorEntry("HER3",   "3KEX", "ERBB3", "rtk", ("ERBB3",)),
    ReceptorEntry("HER4",   "3BBT", "ERBB4", "rtk", ("ERBB4",)),
    ReceptorEntry("FGFR1",  "4UXQ", "FGFR1", "rtk"),
    ReceptorEntry("FGFR2",  "2PSQ", "FGFR2", "rtk"),
    ReceptorEntry("FGFR3",  "4K33", "FGFR3", "rtk"),
    ReceptorEntry("VEGFR1", "3HNG", "FLT1",  "rtk", ("FLT1",)),
    ReceptorEntry("VEGFR2", "2OH4", "KDR",   "rtk", ("KDR",)),
    # FLT4 gene_name annotation is missing on most VEGFR3 entries —
    # use full-text "FLT4" via rcsb_search fallback when re-resolving.
    ReceptorEntry("VEGFR3", "4BSK", "FLT4",  "rtk", ("FLT4",)),
    ReceptorEntry("PDGFR-alpha", "5GRN", "PDGFRA", "rtk", ("PDGFRA",)),
    ReceptorEntry("PDGFR-beta",  "3MJG", "PDGFRB", "rtk", ("PDGFRB",)),
    ReceptorEntry("MET",    "3F82", "MET",   "rtk", ("HGFR",)),
    ReceptorEntry("RET",    "5FM2", "RET",   "rtk"),
    ReceptorEntry("KIT",    "1T46", "KIT",   "rtk", ("c-KIT",)),
    ReceptorEntry("IGF1R",  "2OJ9", "IGF1R", "rtk", ("IGF-1R",)),
    ReceptorEntry("AXL",    "5U6B", "AXL",   "rtk"),
    ReceptorEntry("EPHA2",  "1MQB", "EPHA2", "rtk"),

    # ---- Hormone / Metabolic Receptors ----
    ReceptorEntry("Insulin Receptor",        "1RQQ", "INSR",  "hormone", ("INSR", "IR")),
    ReceptorEntry("GLP-1 Receptor",          "5NX2", "GLP1R", "hormone", ("GLP1R",)),
    ReceptorEntry("Glucagon Receptor",       "4ERS", "GCGR",  "hormone", ("GCGR",)),
    ReceptorEntry("GIP Receptor",            "2QKH", "GIPR",  "hormone", ("GIPR",)),
    ReceptorEntry("Leptin Receptor",         "3V6O", "LEPR",  "hormone", ("LEPR",)),
    ReceptorEntry("Melanocortin 4 Receptor", "6W25", "MC4R",  "hormone", ("MC4R",)),

    # ---- Cytokine / Interleukin / IFN / EPO ----
    ReceptorEntry("IL-2 Receptor alpha", "1Z92", "IL2RA", "cytokine", ("CD25",)),
    ReceptorEntry("IL-2 Receptor beta",  "2B5I", "IL2RB", "cytokine"),
    ReceptorEntry("IL-2 Receptor gamma", "2B5I", "IL2RG", "cytokine", ("CD132",)),
    ReceptorEntry("IL-6 Receptor",       "1N26", "IL6R",  "cytokine", ("CD126",)),
    ReceptorEntry("IL-6 gp130",          "1BQU", "IL6ST", "cytokine", ("CD130",)),
    ReceptorEntry("Erythropoietin Receptor", "1EER", "EPOR", "cytokine", ("EPOR",)),
    ReceptorEntry("IFN-gamma Receptor 1", "1FYH", "IFNGR1", "cytokine"),
    ReceptorEntry("IFN-gamma Receptor 2", "5EH1", "IFNGR2", "cytokine"),

    # ---- Immune Checkpoints ----
    ReceptorEntry("PD-1",   "3BIK", "PDCD1", "checkpoint", ("PDCD1", "CD279")),
    ReceptorEntry("CTLA-4", "1I85", "CTLA4", "checkpoint", ("CD152",)),
    ReceptorEntry("LAG-3",  "7TZG", "LAG3",  "checkpoint", ("CD223",)),
    ReceptorEntry("TIM-3",  "6DHB", "HAVCR2","checkpoint", ("HAVCR2", "CD366")),
    ReceptorEntry("TIGIT",  "3Q0H", "TIGIT", "checkpoint"),

    # ---- Toll-Like Receptors ----
    ReceptorEntry("TLR1", "2Z7X", "TLR1", "tlr"),
    ReceptorEntry("TLR2", "2Z80", "TLR2", "tlr"),
    ReceptorEntry("TLR3", "1ZIW", "TLR3", "tlr"),
    ReceptorEntry("TLR4", "3FXI", "TLR4", "tlr"),
    ReceptorEntry("TLR5", "8AR2", "TLR5", "tlr"),
    ReceptorEntry("TLR7", "5GMF", "TLR7", "tlr"),
    ReceptorEntry("TLR9", "8AR3", "TLR9", "tlr"),

    # ---- GPCRs (chemokine + monoamine + adrenergic) ----
    ReceptorEntry("CCR5",  "4MBS", "CCR5",  "gpcr_amine"),
    ReceptorEntry("CXCR4", "3ODU", "CXCR4", "gpcr_amine"),
    ReceptorEntry("CCR7",  "6QZH", "CCR7",  "gpcr_amine"),
    ReceptorEntry("Dopamine D1", "7CKW", "DRD1", "gpcr_amine", ("D1R",)),
    ReceptorEntry("Dopamine D2", "6CM4", "DRD2", "gpcr_amine", ("D2R",)),
    ReceptorEntry("5-HT1A", "7E2X", "HTR1A", "gpcr_amine"),
    ReceptorEntry("5-HT2A", "6WGT", "HTR2A", "gpcr_amine"),
    ReceptorEntry("ADRB1",  "7BU6", "ADRB1", "gpcr_amine", ("Beta-1 adrenergic",)),
    ReceptorEntry("ADRB2",  "2RH1", "ADRB2", "gpcr_amine", ("Beta-2 adrenergic",)),

    # ---- Ligand-Gated Ion Channels ----
    ReceptorEntry("NMDA GRIN1",  "1PB7", "GRIN1",  "lgic", ("NR1",)),
    ReceptorEntry("NMDA GRIN2A", "5H8H", "GRIN2A", "lgic", ("NR2A",)),
    ReceptorEntry("NMDA GRIN2B", "1OEN", "GRIN2B", "lgic", ("NR2B",)),
    ReceptorEntry("AMPA GRIA1",  "9OVV", "GRIA1",  "lgic", ("GluA1",)),
    ReceptorEntry("AMPA GRIA2",  "3RN8", "GRIA2",  "lgic", ("GluA2",)),
    ReceptorEntry("GABA-A alpha1", "6CDU", "GABRA1", "lgic"),
    ReceptorEntry("GABA-A beta2",  "6X3S", "GABRB2", "lgic"),
    ReceptorEntry("nAChR alpha7", "2MAW", "CHRNA7", "lgic"),
    ReceptorEntry("nAChR alpha4", "6CNJ", "CHRNA4", "lgic"),
    ReceptorEntry("P2X1", "9B73", "P2RX1", "lgic"),
    ReceptorEntry("P2X4", "9BQH", "P2RX4", "lgic"),
    ReceptorEntry("P2X7", "9E3M", "P2RX7", "lgic"),

    # ---- Adhesion (Integrins + Cadherins) ----
    ReceptorEntry("Integrin alphaV", "1U8C", "ITGAV", "adhesion"),
    ReceptorEntry("Integrin beta3",  "1S4X", "ITGB3", "adhesion", ("GPIIIa",)),
    ReceptorEntry("E-Cadherin", "2O72", "CDH1", "adhesion"),
    # CDH2 strict gene-name returned 0; the canonical N-cadherin
    # extracellular domain structures are tagged under different
    # source-organism rows. 1NCG is the classic EC1-EC2 crystal.
    ReceptorEntry("N-Cadherin", "1NCG", "CDH2", "adhesion", ("N-cadherin",)),

    # ---- Death Receptors ----
    ReceptorEntry("Fas",    "2NA7", "FAS",     "death", ("CD95", "TNFRSF6", "APO-1")),
    ReceptorEntry("TNFR1",  "1NCF", "TNFRSF1A","death", ("TNFR1", "CD120a")),

    # ---- Antigen Receptors (TCR + BCR subunits) ----
    ReceptorEntry("TCR alpha (TRAC)",   "3O8X", "TRAC",  "antigen"),
    ReceptorEntry("TCR beta (TRBC2)",   "3O8X", "TRBC2", "antigen"),
    ReceptorEntry("BCR Iga (CD79A)",    "7WSO", "CD79A", "antigen"),
    ReceptorEntry("BCR Igb (CD79B)",    "3KG5", "CD79B", "antigen"),
)


# Reverse lookup index built once at import — maps every alias / name /
# gene symbol (lower-cased) to a single entry. Used by the
# /api/receptors/extended/{name} endpoint to resolve user input.
_INDEX: dict[str, ReceptorEntry] = {}
for _e in EXTENDED_CATALOG:
    _INDEX[_e.name.lower()] = _e
    _INDEX[_e.gene.lower()] = _e
    for _a in _e.aliases:
        _INDEX[_a.lower()] = _e


def lookup(query: str) -> ReceptorEntry | None:
    """Resolve a free-form receptor name / gene / alias to a catalog entry.

    Case-insensitive; ignores surrounding whitespace. Returns ``None``
    when nothing matches (the caller should surface a 404 with the
    user-supplied query echoed back).
    """
    if not query:
        return None
    return _INDEX.get(query.strip().lower())


def to_dict(entry: ReceptorEntry) -> dict:
    """Catalog entry → JSON-serialisable dict for the API response."""
    return {
        "name": entry.name,
        "pdb": entry.pdb,
        "gene": entry.gene,
        "family": entry.family,
        "aliases": list(entry.aliases),
    }


def all_entries() -> list[dict]:
    return [to_dict(e) for e in EXTENDED_CATALOG]


def families_with_counts() -> list[dict]:
    """Family histogram for the UI's group-by header."""
    counts: dict[str, int] = {}
    for e in EXTENDED_CATALOG:
        counts[e.family] = counts.get(e.family, 0) + 1
    return [{"family": f, "count": counts.get(f, 0)} for f in FAMILIES]
