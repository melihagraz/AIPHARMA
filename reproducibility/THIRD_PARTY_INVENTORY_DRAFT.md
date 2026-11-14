# AIPHARMA Third-Party Software Inventory
Date: 2026-11-14
Status: Draft — package-level review pending

| Component | Upstream | License | Usage | Release handling |
|---|---|---|---|---|
| REINVENT4 | https://github.com/MolecularAI/REINVENT4 | Apache-2.0 | 8 adapted modules | Retain LICENSE, NOTICE, modification notices; identify exact upstream revision |
| ADMET-AI | https://github.com/swansonk14/admet_ai | MIT | ADMET predictions | Retain applicable copyright/license; verify bundled source |
| AutoDock Vina | https://github.com/ccsb-scripps/AutoDock-Vina | Apache-2.0 | Docking, v1.2.5 used | Record version and installation; verify if binary bundled |
| DockStream | https://github.com/MolecularAI/DockStream | Apache-2.0 | External comparison | Cite source; modified compatibility fork distributed only after separate review |
| Meeko | https://github.com/forlilab/Meeko | LGPL-2.1-or-later | Ligand preparation, v0.8.0 used | Document dependency and applicable LGPL obligations |
| DOCKSTRING | https://github.com/dockstring/dockstring | Apache-2.0 | Docking-related reference | Verify copied content and preserve required notices |
| DUD-E / DrugCLIP data | Original data source to verify | Not yet cleared | Retrospective benchmarks | Exclude raw redistributed copies pending rights review |

## Original AIPHARMA code
The project team states it owns the original AIPHARMA code.
The team must select and record its intended release license.
Any applicable contractual restrictions must be respected.

## Scientific checkpoints
REINVENT prior weights are not included in the release.
Official reference: https://zenodo.org/records/15641297
Reference SHA256:
b6513ec6dbc54c87ea45cdbf9b4aaefadd7652548b74175366b27f12ec5732fe

The local checkpoint has identical ZIP-member payloads
but a different whole-file hash.

## Remaining package checks
- Identify exact upstream revisions for copied/adapted code.
- Compare these components with the final archive manifest.
- Confirm data/model redistribution rights or exclude the assets.
- Retain applicable LICENSE and NOTICE files.
