# AIPHARMA third-party software and asset inventory

**Updated:** 2026-11-14  
**State:** DRAFT — per-file redistribution review pending; NOT a public-release approval.

| Component | Known upstream / license | Candidate handling and unresolved work |
|---|---|---|
| REINVENT4 | [MolecularAI/REINVENT4](https://github.com/MolecularAI/REINVENT4), Apache-2.0 | Eight adapted core modules are present in `Reinvent_App/reinvent_core/`; retain upstream copyright, license, NOTICE and modification notices. Identify exact upstream revision and compare every adapted file before public distribution. |
| ADMET-AI | [swansonk14/admet_ai](https://github.com/swansonk14/admet_ai), MIT | `Admet_App/` is a service wrapper; its LICENSE and NOTICE are present. Verify whether any upstream code is copied, and its attribution. Downloaded model weights are not included. |
| AutoDock Vina | [ccsb-scripps/AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina), Apache-2.0 | Docking uses an external native executable; source candidate includes configuration files, not a Vina executable. Preserve notices for any adapted/bundled code and separately verify receptor/structure rights. |
| DockStream | [MolecularAI/DockStream](https://github.com/MolecularAI/DockStream), Apache-2.0 | Used for an external benchmark comparison; a modified compatibility fork was exercised locally. Do not imply that fork is distributed in this curated source package without a file-level check. |
| Meeko | [forlilab/Meeko](https://github.com/forlilab/Meeko), LGPL-2.1-or-later (verify precise installed release) | Referenced as a separately installed ligand-preparation dependency. Record version and comply with terms if distributing binaries or adapted code. |
| DOCKSTRING | [dockstring/dockstring](https://github.com/dockstring/dockstring), reported Apache-2.0 (verify exact upstream revision) | Identify any reused implementation or dataset files before redistribution; preserve applicable notices. |
| RDKit / PyTorch / scientific Python | External scientific libraries, each with independent licensing | Not relicensed by AIPHARMA. Reconcile the exact dependency/environment file with distribution method. |
| DUD-E / DrugCLIP data | Source and redistribution terms not yet cleared | Exclude raw mirrored benchmark molecules and data files from public archives until rights and provenance are documented. |
| REINVENT prior weights | [Zenodo record 15641297](https://zenodo.org/records/15641297) | Not included in source candidate. Download externally under applicable terms; verify official whole-file SHA-256 before loading. |
| Receptor / pocket structure assets | Original structural repositories to identify | Audit any structure, configuration, or adapted pocket code included in `Docking_App` or `Pocket_App`; no inferred blanket redistribution rights. |

## First-party licensing decision

The team has selected **Apache-2.0** for original AIPHARMA source code it owns or is authorized to license. The root `LICENSE` and `LICENSE_SCOPE.md` reflect that intent. This decision does **not** resolve other contributors\u0027 ownership, third-party code, or institutional agreements.

## Prior checkpoint verification

Official REINVENT prior `reinvent.prior` SHA-256: `b6513ec6dbc54c87ea45cdbf9b4aaefadd7652548b74175366b27f12ec5732fe`. A previously tested local ZIP-packaged checkpoint had a different whole-file hash but matching six ZIP-member payloads. This does not authorize redistribution or certify unrelated checkpoints.

## Required next actions

1. Map each bundled or adapted file to precise upstream revision, license, notice, and modifications; resolve unknowns.
2. Review all files from **the final Git commit/release archive**, not only the earlier 601 KB selected source tarball.
3. Record first-party contributor/copyright approvals and any institutional restrictions that apply.
4. Keep excluded datasets, weights, and sensitive files excluded unless their rights are documented.
5. Update `PUBLICATION_STATUS.md` with actual completed evidence before making the repository public or issuing a DOI.
