# AIPHARMA v1.0.0 — publication status

**As of 2026-11-14: PRIVATE RELEASE CANDIDATE; NOT CLEARED FOR PUBLIC REDISTRIBUTION.**

The GitHub branch `release/v1.0.0-candidate` holds curated source code. No v1.0.0 GitHub Release, public tag, or AIPHARMA Zenodo DOI has been issued. Version and citation files intentionally indicate the intended **future** version.

## Completed or documented

- The project team states that it owns the original AIPHARMA code and has selected Apache-2.0 as the intended license for those original materials. The first-party license text and scope statement are present.
- A curated source candidate was packaged without locally downloaded model checkpoints, raw DUD-E/DrugCLIP benchmark files, local environments, or private reports.
- Prior Gitleaks checks on the selected candidate package and source-history scopes reported zero findings. A zero-finding scan is not a guarantee against all sensitive material.
- The earlier locally provisioned replay recorded 179 tests passed and 3 deselected. This is not an independent, full-stack reproduction.
- GitHub source is committed in a private release-candidate branch.

## Blocking reviews before public release

- **First-party authorization:** record approvals from every relevant code copyright holder and check any applicable contributor, employer, institutional, or sponsored-research agreement. A repository license file alone cannot establish those rights.
  - In particular, [Brown University's Copyright Ownership and Use Policy](https://policy.brown.edu/policy/copyright-ownership-and-use-policy), section 3.1.5, requires Brown software created in the course of University duties or with University funds to be disclosed to **Brown Technology Innovations (BTI)**, and requires BTI coordination for distribution other than specified academic/research distribution of University-owned software. This policy does **not** by itself prove that Brown owns AIPHARMA; the actual project circumstances must be reviewed. Record the BTI determination or other sufficient ownership evidence before public release.
- **Third-party code and assets:** reconcile every distributed copied/adapted file against upstream origin, revision, license, copyright and NOTICE requirements; audit the selected `Reinvent_App`, `Admet_App`, `Docking_App`, and `Pocket_App` files.
- **REINVENT baseline drift:** a pinned eight-file upstream comparison found six adapted modules with differences beyond attribution headers; source revision and compatibility implications must be resolved or clearly documented. See `reproducibility/REINVENT_UPSTREAM_DIFF_20261114.md`. See the draft third-party inventory.
- **Full candidate review:** scan the actual final Git tree and publishable artifacts for credentials, private data, unlicensed checkpoints, restricted structures, bundled binaries, and licensing conflicts. The existing scans have a defined scope and predate some documentation changes.
- **Documentation and packaging:** finish contributor credit, independent setup limitations, exact scientific dependencies/input acquisition, and release artifact checksum.
- **Publication:** only after the preceding gates are approved, create a public repository/tag/release and deposit a matching archive; then fill in the real DOI and publication date.

**Do not change repository visibility to public or publish the GitHub/Zenodo release solely because these documents now have version 1.0.0.** No legal approval or scientific reproducibility certification is asserted here.
