# AI Pharma App — Setup Guide

End-to-end drug discovery platform. 6-step UI (Pocket → Pharmacophore → De Novo → ADMET → Docking → Candidate Analysis) backed by a multi-service Python stack and an LLM chat layer.

## Checked delivery path (19 October 2026)

The main API was installed in a **fresh Python 3.11 Linux/WSL venv with an empty SQLite database**. The real main application served the built React UI and called the real ADMET HTTP service. Browser sign-in → De Novo → ALK5 designs → paired ADMET → JSON download was exercised. This check covers that path; pocket detection, REINVENT, docking and LLM chat were not rerun. The scientific-model environment, cached ADMET weights and frontend dependency tree were reused. This is not a clean installation of every scientific service or validation of molecule activity/safety.

For the checked main API profile, install `requirements-delivery.txt` in a new environment. The resolved Linux/Python 3.11 package freeze is archived with `reports/alk5_delivery_20261019`. A fresh frontend install uses `npm ci`; the performed check built from the existing lockfile-backed dependency tree, rather than reinstalling it.

Use **`npm run build:delivery`** when the main API serves `dist/`. This build explicitly uses the same origin for API requests, regardless of a developer's `VITE_API_URL` in `.env`; it does not edit that file. Start the main API with the newly created environment's Python:

```bash
python -m uvicorn app.main:app --host 127.0.0.1 --port 8000 --app-dir backend
```

The ALK5 design branch needs RDKit and the bundled `backend/app/data/alk5_design_references.json`, not a REINVENT checkpoint. Paired ADMET additionally needs a running Admet_App at `ADMET_APP_URL` (default port 8030), using its scientific environment with ADMET-AI 2.0.1. A new main API environment alone does not install those model dependencies or weights. Start it separately with that environment's Python:

```bash
python -m uvicorn app.main:app --host 127.0.0.1 --port 8030 --app-dir Admet_App
```

Open the main API URL, sign in using your configured/local account, then **Analysis Results → De Novo Design → ALK5: designs from known actives**. Request designs, optionally compare their ADMET outputs with the parent, and download JSON. The separate `/results` route also exposes the design tool when the database contains no old jobs. No measured activity, safety or synthesis claim follows from generation or a model-output change.

The supervised local check used ports 18765/18766 to keep it separate from normal work. It copied source and built assets without personal `.env` files, generated a fresh database, preserved raw HTTP responses/logs and checked downloaded design/profile identities. `scripts/run_delivery_http_20261019.py` reproduces that recorded **Linux/WSL harness layout**; its scientific-runtime paths are explicit local assumptions, not a universal launcher. See the dated report for commands, scope and remaining installation limits.

## 1. Prerequisites

| Tool | Version | Why |
|---|---|---|
| **Python** | **3.11 or 3.12 for the recorded paths** | Fresh main-API delivery checked on Linux/WSL 3.11; scientific services require their own compatible dependencies/checkpoints |
| **Node.js** | 20.x or 22.x | Vite dev server + frontend build |
| **Git** | any recent | Cloning repo and managing submodules |
| **AutoDock Vina** | 1.2.7 | Pre-bundled wheels work on Windows/Linux/macOS; manual install only if pip fails |
| **OS** | Windows 10/11, macOS 12+, or Linux | Tested on Windows 11 |
| **Disk** | ~10 GB free | `.venv` ≈ 2 GB, node_modules ≈ 500 MB, ChEMBL/PDB caches grow over time |

> Check checkpoint and dependency compatibility for each scientific service. The curated ALK5 design branch does not load a REINVENT checkpoint. A fresh main-API check on Python 3.11 does not establish compatibility of every generator or newer Python release.

## 2. Clone & install

```bash
git clone <THIS_REPO_URL>
cd ai-pharma-app

# Create a clean Python 3.12 virtualenv
python3.12 -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS / Linux:
source .venv/bin/activate

# Install backend deps
pip install --upgrade pip
pip install -r requirements-delivery.txt
# Includes the main API runtime dependencies and curated-design RDKit.
# REINVENT / ADMET / pocket / docking model runtimes are separate.

# Install frontend deps
npm install
```

## 3. Environment

```bash
cp .env.example .env       # or: copy on Windows cmd, Copy-Item on PowerShell
```

Edit `.env`:

| Variable | What |
|---|---|
| `DATABASE_URL` | Default `sqlite:///./backend/enterprise.db` works. Postgres optional. |
| `SECRET_KEY` | Generate: `python -c "import secrets; print(secrets.token_urlsafe(48))"` |
| `ADMIN_PASSWORD` | Pick a strong one. Default `admin / CHANGE_ME` won't work. |
| `ANTHROPIC_API_KEY` | Optional. Get one at https://console.anthropic.com/settings/keys. Empty = chat agent disabled, rest of platform still works. |

## 4. Get the REINVENT 4 prior

The `.prior` checkpoint files are >100 MB each, so they're not in git. Two options:

**Option A — official upstream**
```bash
git clone https://github.com/MolecularAI/REINVENT4.git REINVENT4_official
# Then copy the prior:
cp REINVENT4_official/priors/reinvent.prior Reinvent_App/prior/
```

**Option B — ask whoever shared the repo** for `Reinvent_App/prior/reinvent.prior` and any `prior/finetuned/*_finetuned.prior` files via Drive/Slack.

Without a prior the De Novo step will fail at startup with "prior_loaded: False".

## 5. Start the services

You need **6 processes** running (one frontend, one backend, four microservices). All Python services use the same `.venv`. Run each in its own terminal:

| Service | Port | Command |
|---|---|---|
| **Main backend** (FastAPI) | **8000** | `.venv/Scripts/python -m uvicorn app.main:app --host 127.0.0.1 --port 8000 --app-dir backend` |
| **Yeni motor** (legacy discovery shim) | **8010** | `cd yeni_motor/backend && ../../.venv/Scripts/python -m uvicorn app.main:app --host 127.0.0.1 --port 8010` |
| **Reinvent_App** (REINVENT 4 wrapper) | **8020** | `cd Reinvent_App && ../.venv/Scripts/python -m uvicorn app.main:app --host 127.0.0.1 --port 8020` |
| **Admet_App** (ADMET-AI wrapper) | **8030** | `cd Admet_App && ../.venv/Scripts/python -m uvicorn app.main:app --host 127.0.0.1 --port 8030` |
| **Docking_App** (Vina wrapper) | **8040** | `cd Docking_App && ../.venv/Scripts/python -m uvicorn app.main:app --host 127.0.0.1 --port 8040` |
| **Pocket_App** (RDKit + Biopython pocket detection) | **8050** | `cd Pocket_App && ../.venv/Scripts/python -m uvicorn app.main:app --host 127.0.0.1 --port 8050` |
| **Frontend** (Vite) | **5173** | `npm run dev` |

> **Windows tip**: open 7 PowerShell windows, one per service. macOS/Linux: use tmux or a Procfile runner.

## 6. First-run smoke test

1. Open http://localhost:5173 — login with the admin user from your `.env`
2. Click **+ Add New Job** in the sidebar
3. **Setup card**: Target name `EGFR`, PDB ID `4ZAU`, Generate count `50` (fast smoke)
4. Click **START DE NOVO**
5. Wait ~60-90 seconds → "De Novo process summary" panel should show ~5 candidates saved
6. Go to **Analysis Results → De Novo Design** tab → click any **CHART & DETAIL** button → molecule + descriptors render
7. **Bonus**: open the **AI Pharma Assistant** chat (left sidebar) and type:
   > Generate 100 molecules against EGFR, dock top 5

   Should produce a ranked table with Vina kcal/mol values inside ~90 seconds.

If any step fails, check the per-service terminal logs first.

## 7. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Backend won't import: `ModuleNotFoundError: rdkit` | rdkit not installed in active venv | `pip install rdkit` (must be Python 3.12 venv) |
| Backend won't import: `ModuleNotFoundError: psycopg2` | Even with SQLite, SQLAlchemy dialects pre-load | `pip install psycopg2-binary` |
| `Target 'XYZ' has no recognised signal` (422) | Made-up target name | Use one of EGFR, JAK2, BACE1, KDR, MET, AKT1, LCK, KIT (bundled). Or pass `force_generic=true`. |
| `prior_loaded: False` | Missing `reinvent.prior` file | Section 4 above |
| `No Vina receptor for X` | Target not in `Docking_App/receptors/` | Either add the PDBQT or pick a bundled target |
| Step 4/5/6 dropdown empty | No completed De Novo runs yet | Run Step 3 first |
| Chat agent says `ANTHROPIC_API_KEY missing` | `.env` value empty | Get a key from Anthropic console, paste, restart backend |

## 8. Project layout

```
.
├── backend/             FastAPI app (routes, models, services)
├── src/                 React + Vite frontend
├── Reinvent_App/        REINVENT 4 HTTP wrapper
├── Admet_App/           ADMET-AI HTTP wrapper
├── Docking_App/         AutoDock Vina HTTP wrapper
├── Pocket_App/          Pocket detection + pharmacophore (RDKit + Biopython)
├── yeni_motor/          Legacy synthetic discovery engine (Node + Python)
├── reports/             Markdown research reports (committed)
├── scripts/             One-off maintenance scripts
└── STATUS.md            Project status snapshot
```

## 9. What's in the repo vs not

In git:
- All source code (backend, frontend, microservices)
- `.env.example`, `requirements.txt`, `package.json`
- Documentation (`SETUP.md`, `STATUS.md`, reports)

Not in git (you'll need to provision locally):
- `.env` — secrets, copy from `.env.example`
- `.venv/`, `node_modules/` — install fresh
- `backend/enterprise.db` — SQLite DB; created on first backend start
- `Reinvent_App/prior/reinvent.prior` — 100+ MB checkpoint; see Section 4
- `REINVENT4_official/` — upstream REINVENT 4 source; clone separately if needed

## 10. Where to read first

- [`STATUS.md`](STATUS.md) — current state, what works, what's planned
- `backend/app/api/routes_legacy.py` — all the HTTP endpoints in one file
- `backend/app/services/production_pipeline.py` — the real REINVENT 4 → ADMET → Vina pipeline
- `backend/app/services/agent_service.py` — Phase 2 LLM chat agent + tool registry
- `src/components/NewJobForm.jsx` — the 6-step iterative UI

Welcome aboard.
