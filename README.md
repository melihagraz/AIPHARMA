# AIPHARMA

**A provenance- and evidence-aware research software framework for AI-assisted molecular design and target-specific validation.**

**Status:** v1.0.0 source release **prepared for owner-controlled public publication**. Until GitHub visibility, tag/release, and Zenodo archive are actually published, this remains a private repository without an AIPHARMA DOI. Full independent end-to-end reproduction is not certified. See [publication status](PUBLICATION_STATUS.md).

## Overview

AIPHARMA brings together target-oriented molecular design, scientific-model integrations, docking-related assessment, and provenance/reproducibility checks. The source includes:

- **Target and molecule workflows:** API orchestration and target-specific validation rather than automatically treating a generated molecule as target-active.
- **ALK5 design:** RDKit-based design around documented reference molecules, with optional paired ADMET comparison when the separate ADMET service is running.
- **Scientific engines:** interfaces to REINVENT4 generation, ADMET-AI prediction, AutoDock Vina docking, and pocket/pharmacophore workflows. These require separately configured dependencies and, where applicable, model weights and receptor inputs.
- **Evidence and audit tooling:** scripts for EGFR/ALK5 benchmarks, candidate provenance, release checks, and replay. Recorded benchmark results are not independent validation of biological activity or clinical safety.

## Repository layout

| Path | Purpose |
|---|---|
| `backend/` | FastAPI application, authentication, jobs, target/design routes, and API tests |
| `Reinvent_App/` | REINVENT-compatible generation service and adapted core modules |
| `Admet_App/` | ADMET prediction service wrapper |
| `Docking_App/` | Docking integration and configuration examples |
| `Pocket_App/` | Pocket detection and pharmacophore-related services |
| `scripts/` | Research, analysis, validation, and audit utilities |
| `reproducibility/` | Reproducibility notes, provenance, release checklist, and third-party inventory |

This is a **curated source candidate**. It does not include the full React frontend source/build, raw DUD-E/DrugCLIP panels, downloaded scientific weights, private experiment outputs, or all external receptor/structure inputs. Consequently it is **not** a self-contained turnkey installation of the entire application.

## Setup

Start with [SETUP.md](SETUP.md) and [reproducibility/README.md](reproducibility/README.md). The main API and scientific engines use different environments. For example, from an environment with the backend requirements already installed and application configuration in place:

    python -m uvicorn app.main:app --host 127.0.0.1 --port 8000 --app-dir backend

That command starts **only the API**; it does not install scientific engines, fetch model checkpoints, start ADMET, or supply the frontend. See [SETUP.md](SETUP.md) for the separately configured ADMET service and the limits of the previously checked delivery path. Do not use this research setup as an internet-facing production deployment without a security review; production mode deliberately blocks unreviewed legacy API routes.

## Recorded evidence and limits

- A locally provisioned clean-source replay on 2026-11-14 recorded **179 passing tests and 3 deselections**. Scientific packages/assets were manually supplied; this is **not** an independent clean-machine install of every service.
- EGFR 1M17 redocking in the prior local audit recorded **1.4022 Å RMSD**.
- In a prior retrospective ALK5 panel of **50 actives and 50 presumed decoys**, the AIPHARMA Vina workflow recorded ROC-AUC **0.817**. A separately adapted DockStream comparison recorded **0.815** on that panel. This small observed difference does **not** establish superiority or statistical equivalence.
- Predicted docking scores, generated molecules, and ADMET estimates do not establish experimental potency, toxicity, synthesis feasibility, or clinical utility. Reproduction of all benchmark results from scratch requires additional pinned scientific inputs and runtime details.

For what has and has not been checked, see [REPRODUCIBILITY.md](REPRODUCIBILITY.md), [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md), and the [release checklist](reproducibility/RELEASE_CHECKLIST.md).

## License, data, and citation

The intended license for **original AIPHARMA code whose copyright holders authorized it** is Apache-2.0; see [LICENSE](LICENSE) and [LICENSE_SCOPE.md](LICENSE_SCOPE.md). Third-party components retain their own terms and attribution obligations. The project developer has stated that contributor publication permission was obtained. This is an owner attestation, not independent verification of legal title. **ChEMBL-derived ALK5 data** remains subject to applicable CC BY-SA 3.0 terms and provenance disclosures, separately from Apache-2.0 original code. See [PUBLICATION_STATUS.md](PUBLICATION_STATUS.md), [ALK5 data attribution](backend/app/data/ALK5_REFERENCE_ATTRIBUTION.md) and the [third-party inventory](reproducibility/THIRD_PARTY_INVENTORY_DRAFT.md).

Raw third-party benchmark datasets and model checkpoints are not redistributed in this candidate. The official external REINVENT prior is identified in the reproducibility notes, but must be obtained under its own terms.

Citation metadata is in [CITATION.cff](CITATION.cff). A DOI and release date will be added **only after** an actual public software archive is published.
