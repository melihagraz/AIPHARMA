# Third-party components to review before distributing

Source code previously present in AIPHARMA includes:

- `Reinvent_App/LICENSE` (Apache-2.0 text; verify the code and notices that it covers).
- `Admet_App/LICENSE` (MIT notice to Kyle Swanson; preserve copyright notice).
- `Docking_App/LICENSE_VINA` (Apache-2.0), `Docking_App/LICENSE_DOCKSTRING` (review).
- RDKit, AutoDock Vina, Meeko, ADMET-AI, REINVENT4 and PyTorch depend on their own licenses and packaging terms. Bundled pretrained weights and data may have additional conditions.
- DUD-E/DrugCLIP mirrored molecular data: trace provenance and redistribution permissions rather than assuming mirror access grants redistribution.

This file is **not** a legal audit. Review bundled/reused files and record precise upstream project name, source URL, revision, applicable license and attribution in an approved notice list before public release.
