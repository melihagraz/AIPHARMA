# AIPHARMA reproducibility — v1.0.0 private candidate

**Status:** curated source release candidate on private GitHub branch `release/v1.0.0-candidate`. The version metadata is prepared, but there is no public GitHub Release or AIPHARMA Zenodo DOI.

## Distinct evidence levels

1. **Code-level tests:** locally provisioned source replay reported 179 passed and 3 deselected on 2026-11-14. Scientific dependencies were manually provisioned.
2. **Retrospective scientific checks:** prior local EGFR redocking and ALK5 ranking/benchmark analyses. Source scripts are present, but complete public inputs and independent reproduction are not established.
3. **Independent end-to-end reproduction:** not yet demonstrated. Requires pinned scientific engines, exact native tools, dataset acquisition rights, structures, model weights and a clean environment.

## Source layout and environment

Read [../SETUP.md](../SETUP.md) for main API and separate ADMET service invocation. Backend requirements are in `backend/requirements.txt`; scientific service requirements live in their respective directories and can conflict if forced into one environment. The checked main-API path did not establish all scientific services in a new clean environment.

The REINVENT model checkpoint is **not bundled**. Official model reference: [Zenodo 15641297](https://zenodo.org/records/15641297), official file SHA-256 `b6513ec6dbc54c87ea45cdbf9b4aaefadd7652548b74175366b27f12ec5732fe`. Never run untrusted PyTorch checkpoints without validating provenance.

## Reproduction and audit

The `scripts/` directory contains EAVE audits, ALK5 validation scripts, docking controls and release packaging utilities. Some scripts require local data and software not included in this source candidate; do not treat a missing dataset as a failed numerical replication.

Recorded scientific evidence must include fixed target/panel IDs, protein preparation, docking configurations, executable/package versions, input and output hashes, and uncertainties. The absence of all of those from a public portable archive limits what third parties can reproduce.

For the publication criteria, see [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md), [THIRD_PARTY_INVENTORY_DRAFT.md](THIRD_PARTY_INVENTORY_DRAFT.md), and [../PUBLICATION_STATUS.md](../PUBLICATION_STATUS.md). The root Apache-2.0 license covers only authorized original code; reused components retain their upstream terms.
