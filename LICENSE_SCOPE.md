# AIPHARMA license scope

**License for first-party code:** Apache License, Version 2.0 (`LICENSE`).
**Copyright:** © 2026 Melih Agraz, for original material that he owns or
is authorized to license.

## Covered

Only original AIPHARMA source code and original project documentation for
which the named licensor actually owns the copyright or has authorization to
license. This statement does **not** establish title or override an employment,
collaboration, funding, or intellectual-property agreement.

## Excluded or separately licensed

- `Reinvent_App/`: derived/adapted REINVENT components may be subject to the
  upstream Apache-2.0 license, existing copyright and NOTICE terms.
- `Admet_App/`: preserve ADMET-AI/MIT notices as applicable; models have to be
  reviewed separately.
- `Docking_App/`: preserve relevant AutoDock Vina and Dockstring licenses and
  notices; check binaries, receptor files and any adapted code individually.
- `Pocket_App/`: check source-file provenance and redistribution terms for any
  embedded third-party code or receptor/structure assets. The curated candidate
  does not contain a `Pocket_App/NOTICE` file.
- All other third-party code, trained weights, datasets and assets retain their
  original rights and terms. A public download link does not prove that those
  materials may be redistributed.

The directories named above are **not automatically excluded from public
release**; they require a component-level license review first. This file
does not replace or change any upstream license.

**ALK5 data exception:** `backend/app/data/alk5_design_references.json` is a
curated data asset declared as derived from ChEMBL records. It is **not**
covered by the first-party Apache-2.0 software license. ChEMBL-derived
adaptations are subject to CC BY-SA 3.0 attribution/share-alike provisions;
any unverified additional source materials must be separately cleared. See
`backend/app/data/ALK5_REFERENCE_ATTRIBUTION.md` and
`reproducibility/RELEASE_RIGHTS_DECISION_20261114.md`.

## Before public release

1. Confirm rights to the original code, including any employment, sponsored
   research, university, contributor, and collaborator agreements applicable
   to AIPHARMA (Brown, Brown University Health, Giresun, or others as relevant).
2. Inspect the actual current repository, not just an old archive. Preserve
   third-party license texts, copyright notices, and attribution requirements.
3. Exclude private credentials, restricted data, unlicensed model weights,
   and files that should not be redistributed.
4. Verify install/tests; then tag a release and obtain a DOI.

Until these are completed, regard this as a **license configuration candidate**,
not proof of intellectual-property clearance or a public-ready software release.
