# Public reproducibility release checklist

Status: NOT COMPLETE for AIPHARMA v0.1.0-rc1. Use as a checklist, not evidence of completion.

- [x] Version file and local candidate packaging script created.
- [x] Initial citation metadata prepared (author/contributor review pending).
- [x] Unit-test command, local environment inventory and checksum generation provided.
- [x] Source-only archive excludes datasets, secrets, model weights, venvs and private reports by default.
- [x] External data/model/structure provenance documented with known SHA-256 values.
- [ ] All coauthors and code owners approve the copyright/license choice; root LICENSE file supplied.
- [ ] Upstream sources and bundled third-party code are audited for license compliance and attribution.
- [ ] Check exact installed Python scientific dependencies; resolve compatibility/conflicts; retain OS-specific reviewed lock.
- [ ] Native Vina executable version and SHA-256, receptor/grid files, weights and model versions are captured.
- [ ] Clean-environment EAVE unit tests and full scientific install tested on separate machine/container.
- [ ] Full 100-molecule panel reproduced with exact fixed panel and inputs in a clean environment.
- [ ] Publicly usable data acquisition and redistribution process established; identifiers and data hashes made public.
- [ ] Source archive manually inspected for accidental secrets/data/copyright violations.
- [ ] GitHub repo established and accessible, exact commit tagged v0.1.0-rc1 (or later version).
- [ ] GitHub release artifact SHA-256 matches verified local candidate.
- [ ] DOI from Zenodo/software archive recorded and linked from CITATION.cff.
- [ ] Manuscript availability statement updated with actual links, commit, data accessions and DOI.

**A locally packaged source candidate is not a public or independently reproduced release.**
- [ ] Root README is updated from Phase 6-era notes to reflect Phase 15 evidence and remaining limitations.
