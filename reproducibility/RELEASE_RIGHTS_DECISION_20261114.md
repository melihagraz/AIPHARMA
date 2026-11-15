# AIPHARMA v1.0.0 — publication permissions and rights record

**As of:** 2026-11-14. **Publication decision:** the owner requested public release and stated that the necessary contributor permissions had been granted. This is a **developer/owner attestation**, not independent proof of ownership, an institutional clearance certificate, or a signed license transfer.

## Original AIPHARMA contributions

- Developer declaration: AIPHARMA was developed independently of Brown University duties/resources.
- Developer declaration: permission to publicly distribute the original/other contributor AIPHARMA code was granted.
- Chosen license for original AIPHARMA code where copyright holders have authority: **Apache License 2.0**, with its patent and attribution provisions.
- Actual contributor-by-contributor written approval, employment agreements and institutional assignment status were not independently inspected. The publishing owner remains responsible for verifying this assertion. Brown ownership or approval is **not** presumed solely from affiliation.

## REINVENT4

Upstream REINVENT4 is distributed under **Apache-2.0**, which permits redistribution and modifications subject to proper license delivery, retained notices and documentation of changes. AIPHARMA includes attribution and an adapted-code notice and records the pinned current-upstream eight-file comparison in `reproducibility/REINVENT_UPSTREAM_DIFF_20261114.md`. Six of those files differ substantively from the pinned upstream version; the exact original fork revision is unknown. There is **no assertion of strict upstream equivalence**.

## ALK5 / ChEMBL data

The ALK5 JSON references 144 ChEMBL-identified molecules and an exclusion set of 2,514 normalized parent identities. The source database is **CC BY-SA 3.0**, with license/attribution disclosed at `backend/app/data/ALK5_REFERENCE_ATTRIBUTION.md`. The derived data must not be presented as Apache-2.0-licensed code. The source `alk5_structure_mapping_20261018/completed_panel_nodes.json`, exact ChEMBL release and full input/transformation pipeline have not been recovered.

**Decision on reproducibility:** release can transparently describe an unverified historical data derivation rather than claim full reconstruction. **Decision on third-party rights:** the owner's stated permission cannot confer rights to unknown, separately restricted inputs. The publisher is responsible for ensuring no such material is included; this is not independently certified here.

## Other dependencies and distribution

The code/inventory reference Vina, Meeko, ADMET-AI, DockStream, DOCKSTRING and associated upstream distributions. Verify redistributed vendor/adapted files against their exact license text and notices; external optional dependencies and models are not automatically covered by Apache-2.0.

## Publication record

No public GitHub Release, DOI or final independent full-system certification existed at preparation time. Once the owner makes this repo public, links Zenodo first and creates `v1.0.0`, update `CITATION.cff` with the actual DOI.

References:
- https://github.com/MolecularAI/REINVENT4/blob/main/LICENSE
- https://www.ebi.ac.uk/chembl/
- https://creativecommons.org/licenses/by-sa/3.0/
- https://policy.brown.edu/policy/copyright-ownership-and-use-policy

This report is an evidence inventory and explicit record of user attestation; it is not legal advice.
