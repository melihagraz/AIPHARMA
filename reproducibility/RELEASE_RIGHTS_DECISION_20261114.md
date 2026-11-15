# AIPHARMA v1.0.0 — rights and distribution decision ledger

**Reviewed:** 2026-11-14. **Scope:** private candidate branch `release/v1.0.0-candidate`, parent commit `0af33216c4dc3ff05f36d470f27ee293091d4c70`.

**Decision:** NOT YET AUTHORIZED FOR PUBLIC DISTRIBUTION. A software license text or a successful security test cannot substitute for evidence of ownership/authorization.

| Scope | Confirmed from sources | Release status / additional action |
|---|---|---|
| Original AIPHARMA source | Team states it owns original work and has selected Apache-2.0. Root LICENSE and LICENSE_SCOPE exist. | **HOLD** until each relevant rightsholder/author approval and employment, university/funding/sponsored-work scope are documented. Do not assume repository owner owns every contribution. |
| REINVENT4 adapted eight modules | Official upstream license is Apache-2.0; eight local headers include source and change notices; local LICENSE/NOTICE exist; no root upstream NOTICE file at pinned reference revision. | **Basic license permission confirmed, remaining attribution/provenance review**: check exact source history/contributor notices and match final release tree; no algorithmic equivalence claim. |
| ChEMBL-derived ALK5 JSON | ChEMBL officially publishes source data under CC BY-SA 3.0 (commercial use included with conditions). JSON has 144 records and 2,514 exclusion identities, and identifies ChEMBL/document/activity sources. | **HOLD for entire dataset**: original intermediate and ChEMBL release/raw transformation chain unverified. Keep CC BY-SA attribution and adapted-data share-alike distinct from Apache-2.0 code. Check whether other restricted inputs or commercial calculated properties were included. |
| Other software/dependencies | Existing inventory covers ADMET-AI, Vina, Meeko, DOCKSTRING, DockStream; relevant licensing has been identified. | Check any source actually redistributed, attribution, structure/data rights, and exact final tree before asserting clearance. |

## Brown University institutional disclosure check

The [Brown University Copyright Ownership and Use Policy](https://policy.brown.edu/policy/copyright-ownership-and-use-policy) applies to relevant members of the Brown community, including academic appointees and certain affiliated-hospital faculty. Section 3.1.5 addresses computer software created in the course of University duties or using University funds, requires disclosure to [Brown Technology Innovations (BTI)](mailto:tech-innovations@brown.edu), and directs distribution of University-owned software to be coordinated with BTI except for the specified limited academic/research pathway.

**This does not establish that AIPHARMA is Brown-owned.** Whether the software was independent, developed on the contributor's personal time/resources, sponsored, assigned, or built using University resources must be factually determined. For any identified additional institutional rightsholder, obtain its approval or written statement of non-interest/authority as appropriate. Avoid uploading sensitive employment agreements to this public-intended repository.

## What would close the pending gates

1. A written **first-party rights decision** identifying actual AIPHARMA authors/copyright holders, the code they contributed, applicable employer/sponsor obligations, and permission to publish under Apache-2.0; if Brown policy may apply, record the BTI inquiry and response outside the public repo.
2. A **data-source declaration** identifying ChEMBL source release and enough evidence to show the ALK5 JSON contains no independently restricted third-party input; retain CC BY-SA attribution and modifications notice.
3. Final-file-by-file check of **the exact public candidate commit** for notices, sensitive material and dataset/model/binary rights, followed by source-archive checksum.
4. Only after signoffs, update this decision ledger to record the specific approval evidence identifiers (avoid personal/confidential details), merge/tag/release and then record the actual Zenodo DOI.

## Authority and limitations

Upstream license and documentation sources:
- https://github.com/MolecularAI/REINVENT4/blob/main/LICENSE
- https://www.ebi.ac.uk/chembl/
- https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions
- https://creativecommons.org/licenses/by-sa/3.0/
- https://policy.brown.edu/policy/copyright-ownership-and-use-policy

This is an evidence-led publication checklist and interpretation of publicly posted terms, **not a signed rightsholder release or legal determination**.
