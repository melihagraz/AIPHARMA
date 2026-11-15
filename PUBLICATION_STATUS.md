# AIPHARMA v1.0.0 — publication status

**As of 2026-11-14: PRIVATE RELEASE CANDIDATE; NOT CLEARED FOR PUBLIC REDISTRIBUTION.**

The GitHub branch `release/v1.0.0-candidate` holds curated source code. No v1.0.0 GitHub Release, public tag, or AIPHARMA Zenodo DOI has been issued. Version and citation files intentionally indicate the intended **future** version.

## Completed or documented

- The developer states **AIPHARMA was developed independently of Brown duties/resources** (declaration recorded 2026-11-14); team states it owns original AIPHARMA code and has selected Apache-2.0. Root license and scope documents are present. This is not an independent legal ownership finding.
- A curated source candidate was packaged without locally downloaded model checkpoints, raw DUD-E/DrugCLIP benchmark files, local environments, or private reports.
- Prior Gitleaks checks on the selected candidate package and source-history scopes reported zero findings. A zero-finding scan is not a guarantee against all sensitive material.
- The earlier locally provisioned replay recorded 179 tests passed and 3 deselected. This is not an independent, full-stack reproduction.
- GitHub source is committed in a private release-candidate branch.

## Blocking reviews before public release

- **First-party authorization:** Independent-development declaration recorded. Confirm all actual code copyright holders/contributors and any relevant employment, grant, commissioned-work, institutional or third-party assignment obligations; record rights to license original contributions under Apache-2.0.
  - **Brown is not presumed to own this independent project.** [Brown's copyright policy](https://policy.brown.edu/policy/copyright-ownership-and-use-policy) addresses software created during duties or with University resources/funds, and contains disclosure provisions. If facts/agreements or commercial potential trigger institutional review, seek a BTI determination; do not characterize BTI approval as already obtained or automatically required solely because the developer is affiliated with Brown.
- **Third-party code and assets:** reconcile every distributed copied/adapted file against upstream origin, revision, license, copyright and NOTICE requirements; audit the selected `Reinvent_App`, `Admet_App`, `Docking_App`, and `Pocket_App` files.
- **REINVENT baseline drift:** a pinned eight-file upstream comparison found six adapted modules with differences beyond attribution headers; source revision and compatibility implications must be resolved or clearly documented. See `reproducibility/REINVENT_UPSTREAM_DIFF_20261114.md`. See the draft third-party inventory.
- **Full candidate review:** scan the actual final Git tree and publishable artifacts for credentials, private data, unlicensed checkpoints, restricted structures, bundled binaries, and licensing conflicts. The existing scans have a defined scope and predate some documentation changes.
- **Documentation and packaging:** finish contributor credit, independent setup limitations, exact scientific dependencies/input acquisition, and release artifact checksum.
- **Publication:** only after the preceding gates are approved, create a public repository/tag/release and deposit a matching archive; then fill in the real DOI and publication date.

**Do not change repository visibility to public or publish the GitHub/Zenodo release solely because these documents now have version 1.0.0.** No legal approval or scientific reproducibility certification is asserted here.
