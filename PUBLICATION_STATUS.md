# AIPHARMA v1.0.0 — publication status

**Publication:** GitHub [v1.0.0 Release](https://github.com/melihagraz/AIPHARMA/releases/tag/v1.0.0) was published on **2026-11-15 (UTC)**, targeting commit `0d8cbe6d0bafc992cf5c12d142d18d0812a966a7`. The publisher reports Zenodo DOI [10.5281/zenodo.24236794](https://doi.org/10.5281/zenodo.24236794). The DOI was provided by the publisher; its record metadata could not be independently retrieved in this documentation update.

## Owner authorization declaration

The developer stated on 2026-11-14 that AIPHARMA was developed independently of Brown duties/resources and that publication permission has been granted for contributor/original-code rights. The intended first-party license is Apache-2.0. **These are developer assertions, not independently verified signed releases or institutional legal determinations.** No assumption of Brown ownership is made from academic affiliation.

Publishing the repository does not waive any actual contributor or institutional copyright. The publishing owner is responsible for verifying that this attestation covers *every* distributed first-party contribution and applicable employment/sponsor agreements. Where university policy applies, consult the appropriate office.

## Included and excluded

The curated source includes API and scientific integration code and a ChEMBL-derived ALK5 reference asset with separate CC BY-SA 3.0 attribution. The full frontend build, downloaded REINVENT checkpoints, DUD-E/DrugCLIP raw panels, receptor inputs and local environments are excluded. Third-party code retains its source license; see `LICENSE_SCOPE.md`, `reproducibility/THIRD_PARTY_INVENTORY_DRAFT.md`, REINVENT notices and `backend/app/data/ALK5_REFERENCE_ATTRIBUTION.md`.

## Disclosed provenance and reproducibility limitations

- The ALK5 JSON source declaration cites an intermediate `completed_panel_nodes.json` and archived ChEMBL activity records, but the original intermediate file, precise ChEMBL release and exact transformation chain were not located. This is a **documented reproducibility gap**, not evidence of an independently reconstructed dataset. The developer's general permission statement does not independently establish the rights of any unknown additional upstream source.
- The REINVENT4-adapted code follows Apache-2.0 licensing and is accompanied by attribution/change notices; the exact historical base revision remains unestablished and six file diffs from a pinned upstream comparison are documented.
- Previously reported local testing: 179 passed and 3 deselected with provisioned dependencies; ALK5 count=1 targeted design smoke test passed. This is not an independent clean-machine end-to-end replay. Scientific molecular predictions are not experimentally validated.
- Prior Gitleaks scan of the candidate reported no findings; a zero-finding scan does not prove the absence of all secret or restricted content.

## Publication provenance

- Public GitHub release: [v1.0.0](https://github.com/melihagraz/AIPHARMA/releases/tag/v1.0.0); target commit `0d8cbe6d0bafc992cf5c12d142d18d0812a966a7`.
- Zenodo DOI **as supplied by the publisher**: [10.5281/zenodo.24236794](https://doi.org/10.5281/zenodo.24236794). Whether it is the version-specific DOI versus a concept DOI should be checked against the Zenodo record metadata; the same DOI is used consistently in the current `main` citation metadata.
- The initial v1.0.0 GitHub tag and Zenodo archive were **not modified** when the DOI badge and citation files were added to `main` after publication. This preserves the archived version's source snapshot.
- This record documents the owner's release and retained limitations; it is not a legal certification or independent complete-reproduction guarantee.
