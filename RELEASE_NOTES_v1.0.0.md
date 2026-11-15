# AIPHARMA v1.0.0 — Source release

A provenance- and evidence-aware research software framework for AI-assisted molecular design and target-specific validation.

## Source package

- Curated FastAPI and scientific-integration code for ALK5 analogue design, REINVENT-compatible generation, ADMET prediction, docking and pocket/pharmacophore workflows.
- ALK5 reference asset derived from ChEMBL: 144 curated reference structures and 2,514 exclusion identities, accompanied by a separate data-license/attribution notice.
- Test and reproducibility utilities, third-party notices and REINVENT4 comparison reports.

## Recorded local evidence

- 179 tests passed, 3 deselected in a previously provisioned environment, not a fully independent installation.
- One ALK5 candidate was generated from the packaged reference data in a targeted smoke test.
- Prior EGFR 1M17 redocking RMSD 1.4022 Å and retrospective ALK5 Vina ROC-AUC 0.817 on 50 active/50 presumed decoy compounds.

## Important limits

- Not a turnkey full frontend/service deployment. The full React build, scientific model checkpoint files, raw DUD-E/DrugCLIP benchmark data, private outputs and all native docking/receptor files are not included.
- The exact ALK5 intermediate source and complete ChEMBL release/transformation lineage could not be recovered, and independent complete reconstruction is not shown.
- The precise historical base of the REINVENT adaptations remains uncertain; differences relative to a pinned current upstream revision are documented.
- In silico generated candidates, docking ranks and ADMET estimates are not evidence of measured potency, safety, synthesis feasibility or clinical efficacy.

## Licensing

- Original software distributed with relevant contributor authority: Apache License 2.0.
- ChEMBL-derived ALK5 dataset: CC BY-SA 3.0, with dedicated attribution; not relicensed as Apache-2.0.
- Adapted REINVENT and other third-party components retain their original license/copyright notice obligations.
- Developer reports having the contributor permissions for public release. Refer to PUBLICATION_STATUS.md for the scope and limitations of that assertion.

**Zenodo DOI (v1.0.0):** [10.5281/zenodo.24236794](https://doi.org/10.5281/zenodo.24236794). The published GitHub release is [v1.0.0](https://github.com/melihagraz/AIPHARMA/releases/tag/v1.0.0). DOI metadata was supplied after publication; this link was added to the release notes without moving the version tag.
