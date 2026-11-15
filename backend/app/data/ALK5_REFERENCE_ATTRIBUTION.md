# ALK5 reference dataset — attribution, data license, and provenance

**Artifact:** `backend/app/data/alk5_design_references.json`
**SHA-256:** `6378cbe612416fef4872db936f22c3f44ecf255fcc14827d8707deb6c7b8d24b`
**Status:** Private release candidate; exact full source chain still unverified.

## Content

- 144 reference records, with ChEMBL compound identifiers, bibliographic document identifiers and activity identifiers.
- 2,514 normalized parent identities in the exclusion set.
- Current JSON review: 15 distinct document identifiers and 150 distinct activity identifiers referenced among the selected 144 records.
- The project's original selection description: 16 editable parent scaffolds plus up to 128 measured donors selected from enzyme-catalytic IC50 records; the records are not confirmed wild type. The dataset does not transfer measured parent activities to newly designed candidates.

## Source and transformation disclosure

The embedded `provenance.source` names `alk5_structure_mapping_20261018/completed_panel_nodes.json` and archived ChEMBL raw activity rows. A 2026-11-14 local search did **not** locate the original `completed_panel_nodes.json` file. The exact ChEMBL database release, archived raw source rows and full transformation script/input hashes are **not** established by the present evidence. The dataset's embedded `provenance.dataset_sha256` (`98659fa95c0968fa2d4ee3b67054ff9c96745acc90255ae6acc13c997dc7ff12`) should not be confused with this JSON file's SHA-256 above.

Declared processing includes ChEMBL-record filtering, parent selection, structure/SMILES normalization, stereochemical checks, scaffold selection and exclusion of observed parent identities. These steps represent a curated and transformed dataset, not a verbatim unmodified ChEMBL dump. The full transformation cannot currently be reconstructed from the archived inputs.

## License — separate from AIPHARMA software

**ChEMBL source database:** Creative Commons Attribution-ShareAlike 3.0 Unported (**CC BY-SA 3.0**).

- ChEMBL data: https://www.ebi.ac.uk/chembl/
- Official license explanation: https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions
- License text: https://creativecommons.org/licenses/by-sa/3.0/legalcode
- Attribution: Zdrazil B. et al. *The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods.* Nucleic Acids Research. DOI: https://doi.org/10.1093/nar/gkad1004

To the extent that this JSON is an adaptation of ChEMBL-provided database material, applicable ChEMBL-derived material and adaptation contributions must be shared under the applicable CC BY-SA terms, with appropriate attribution, source/link and identification of modifications. **The root Apache-2.0 license for original AIPHARMA code does not relicense the JSON or override other source rights.**

The ChEMBL FAQ warns that some calculated compound properties can originate from commercially licensed software. Whether any such protected feature is actually included here has **not** been established. The availability of ChEMBL IDs, SMILES and activity citations does not itself prove the origin or redistribution rights of every possible upstream input.

**Public-distribution decision:** PENDING confirmation that the entire JSON derives from licensed ChEMBL material and independently authored transformations without restricted extra inputs, plus documentation of the exact source release when recoverable. The candidate remains private pending that review.
