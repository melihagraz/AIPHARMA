# REINVENT4 adapted-source comparison — 2026-11-14

**Target repository:** `melihagraz/AIPHARMA`, private `release/v1.0.0-candidate`, base commit `6f24111a41141089c624e224c07ee59d31408296`.

**Reference repository:** [MolecularAI/REINVENT4](https://github.com/MolecularAI/REINVENT4/tree/0f075b349ed22015e0406350c6f426701bb0815d), upstream reference commit `0f075b349ed22015e0406350c6f426701bb0815d`. This is a **fixed comparison reference**, not evidence that this was the original source revision from which AIPHARMA was adapted.

**Comparison:** Source files in `Reinvent_App/reinvent_core/` were compared to the upstream paths listed in `Reinvent_App/NOTICE`. Local attribution header blocks were excluded; the remaining source lines were compared using a longest-common-subsequence line-diff. Counts indicate differences against this particular reference commit and do not imply a security or scientific defect.

| AIPHARMA file | Upstream path | Local-only lines / reference-only lines | Main difference |
|---|---|---:|---|
| `model.py` | `reinvent/models/reinvent/models/model.py` | 5 / 4 | Relative imports, metadata type reference, and changed tokenizer argument wiring |
| `rnn.py` | `reinvent/models/reinvent/models/rnn.py` | 1 / 2 | Local device argument differs from upstream `self.to(self.device)` behavior |
| `vocabulary.py` | `reinvent/models/reinvent/models/vocabulary.py` | 0 / 0 | Identical after attribution header removal |
| `dataset.py` | `reinvent/models/reinvent/models/dataset.py` | 4 / 5 | Upstream includes `isomeric` parameter, passed to SMILES randomization; local implementation differs |
| `utils.py` | `reinvent/models/reinvent/utils.py` | 1 / 1 | Upstream tensor allocation specifies `device=encoded_seqs[0].device`; local version does not |
| `conversions.py` | `reinvent/chemistry/conversions.py` | 9 / 12 | SMILES randomization implemented differently; upstream provides `isomeric`/`doRandom=True` behavior |
| `meta_data.py` | `reinvent/models/meta_data.py` | 10 / 135 | Local class is reduced for checkpoint compatibility; upstream hash/metadata utilities and default fields differ |
| `model_mode_enum.py` | `reinvent/models/model_mode_enum.py` | 0 / 0 | Identical after attribution header removal |

## Implications

- The previous notice's blanket claim that the algorithms and tokenization behavior were unchanged is **not adequately supported**. The notice has been corrected to describe the observed differences.
- Presence of an Apache-2.0 license and prominent adapted-file attribution is positive evidence for notice obligations, but does not establish that all old upstream copyright and NOTICE requirements have been satisfied.
- In particular, changes in SMILES randomization, stereochemical handling, tensor device placement and metadata loading could change behavior across environments; these deserve separate engineering scrutiny before claiming full compatibility.
- The locally observed successful two-molecule REINVENT sampling and the 179-test provisioned replay **do not show** that all upstream behaviors are equivalent.
- This comparison did **not** examine all possible historic upstream revisions, nor determine the true source commit used by the team.

## Outstanding release actions

1. Identify the original REINVENT4 revision used for each adapted module, or provide a documented best match and disclose uncertainty.
2. Preserve required upstream license, NOTICE and prominent per-file modification notices; perform a contributor provenance check.
3. If claiming compatibility, assess the scientific effects of the listed differences using targeted tests that do not silently change the release package.
4. Review the final source tree and publication approvals before public release. Keep this report with the versioned artifact.
