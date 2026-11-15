# AIPHARMA — preliminary third-party license evidence (2026-11-14)

**Audit target:** private release candidate `f64547d0afafb99774d0881c2dded84d61ebd36c` (before this documentation commit). **Not a legal opinion or a completed distribution audit.** License names below were checked against the indicated project pages, not necessarily the precise historic revision copied into this project.

| Component | Upstream source | Checked license | What is present in the curated Git tree | Remaining issue |
|---|---|---|---|---|
| REINVENT4 | https://github.com/MolecularAI/REINVENT4/blob/main/LICENSE | Apache-2.0 | Eight files in `Reinvent_App/reinvent_core/`; all eight checked headers name upstream, copyright, license, source path and modification intent; `Reinvent_App/LICENSE` and `NOTICE` are present | Pin upstream revision, compare each adapted file against it, and reconcile notice wording with actual changes |
| ADMET-AI | https://github.com/swansonk14/admet_ai/blob/main/LICENSE.txt | MIT | AIPHARMA service-wrapper source, license and notice. No downloaded weights in this tree | Verify code copying versus dependency, weights/data licenses if a future bundle includes them |
| AutoDock Vina | https://github.com/forlilab/AutoDock-Vina/blob/develop/LICENSE | Apache-2.0 | Docking service source, license copy and configuration files | No Vina executable in candidate; native dependencies and exact binary SHA/version must be documented externally |
| DockStream | https://github.com/MolecularAI/DockStream | Apache-2.0 | Used for prior local comparison; modified compatibility fork is not evidenced as redistributed here | If fork is published separately, inspect modified files and notices |
| Meeko | https://github.com/forlilab/Meeko/blob/develop/LICENSE | LGPL-2.1 | Referenced as separately installed dependency, not copied into candidate | Verify installed version, precise licensing expression and binary/linked distribution obligations if applicable |
| DOCKSTRING | https://github.com/dockstring/dockstring | Apache-2.0 | A license copy and eight Vina grid configuration text files are present, but no receptor PDBQT structures | Verify provenance and whether the included grid configurations are derived from DOCKSTRING and preserve the relevant attribution |
| RDKit, PyTorch and other dependencies | Package-specific upstream projects | Various | Import/runtime requirements only | Record precise licenses for **redistributed** source/binaries or container builds; dependencies are not automatically relicensed |
| DUD-E / DrugCLIP panels | Dataset-specific original sources | Unresolved | Raw benchmark dataset mirrors not included in this selected source tree | Keep omitted until provenance and redistribution rights have been checked |
| REINVENT prior model | https://zenodo.org/records/15641297 | Independently assess model/data license | Downloaded checkpoint is **not** included | Obtain independently and verify official SHA-256 before loading |

## Findings in the repository's own notices

- `Reinvent_App/NOTICE` provides upstream attribution and a list of eight adapted modules. Header inspection verified that each corresponding file contains an attribution/change note. This was followed by a separate line-by-line comparison against a pinned upstream reference; see `REINVENT_UPSTREAM_DIFF_20261114.md`. The historical source revision used to make the adaptations remains unidentified.
- `Admet_App/NOTICE` identifies ADMET-AI as a separately installed MIT dependency; the service wrapper is in the candidate.
- The previous `Docking_App/NOTICE` referred to `bin/vina.exe` and bundled `*_target.pdbqt` receptor files. **Neither class of file was found in the private candidate Git tree**; the notice was corrected in the same documentation change.
- `Pocket_App/` contains application source, but no `Pocket_App/NOTICE` file was found in the candidate. The scope statement no longer points to a non-existent notice. Any third-party embedded code or structure requires per-file review.

## REINVENT4 license position (separate from scientific equivalence)

- The REINVENT4 upstream `LICENSE` at the pinned comparison revision `0f075b349ed22015e0406350c6f426701bb0815d` is **Apache-2.0**. It permits derivative redistribution and commercial use subject to license delivery, retained attribution/copyright notices and prominent notices of changes.
- The pinned upstream root directory has no standalone `NOTICE` file; the AIPHARMA package nevertheless includes an adapted-code NOTICE and an Apache-2.0 license copy. Each of eight adapted module headers records source/attribution and modification intent.
- This is positive evidence that basic notice requirements are **addressed**, not a certification of perfect provenance or every copyright contributor. The origin revision remains unknown and six current-source diffs are recorded. A pin/diff is a reproducibility issue, not an automatic prohibition under Apache-2.0.

## Conditions still open before public distribution

1. For eight REINVENT files, identify the exact upstream commit/version and produce a per-file original/adapted mapping with attribution and meaningful modification details.
2. Compare the whole **final** public release tree against this private candidate tree, including tracked executable/data files, not only file extensions.
3. Check attribution for docking configuration files, all other copied/modified code, and any additional distributed assets.
4. Document first-party contributor authorization and any actual contractual limits; Apache-2.0 selection for the team-owned code is already recorded.
5. Finish sensitive-file scan on final release artifacts and publish matching SHA-256 provenance before the repository is made public.

See `PUBLICATION_STATUS.md` and `reproducibility/THIRD_PARTY_INVENTORY_DRAFT.md`. Upstream license identification alone **does not constitute approval** to publish an arbitrary copied file or dataset.
