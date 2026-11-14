from pathlib import Path
import re,json,hashlib
p=Path.home()/'restore_full_alk5';p.mkdir(exist_ok=True)
old=(p/'original.tex').read_text(encoding='utf-8-sig');short=(p/'short.tex').read_text(encoding='utf-8-sig');s=old
changes=[]
def replace(a,b):
 global s
 if a not in s:raise RuntimeError('Missing text: '+a[:90])
 s=s.replace(a,b);changes.append(a[:100])
replace('AIPHARMA: Integrated AI-Powered Structure-Based Drug Discovery Pipeline for DE NOVO Drug Design in Pulmonary Arterial Hypertension','AIPHARMA: An Integrated Computational Drug Discovery Workflow with Endpoint-Aware Retrospective Evaluation of ALK5 Activity')
start=s.index('Pulmonary arterial hypertension (PAH) is a progressive, multifactorial cardiopulmonary disorder')
end=s.index('\\textbf{Keywords:}',start)
abstract=r'''Pulmonary arterial hypertension (PAH) motivates the investigation of signaling-associated targets, but computational candidate generation alone cannot establish therapeutic activity. Here we present AIPHARMA, a modular workflow combining pocket characterization, pharmacophore analysis, unconditional molecular generation, docking and predicted ADMET profiling. We retain the original structural and generative experiments while evaluating the limits of their interpretation through additional retrospective controls.

Internal pose controls supported pose generation for selected ALK5 complexes, whereas a separate, scaffold-excluded 36-structure ranking test did not demonstrate reliable docking discrimination: Vina ROC--AUC was 0.441 (35 scored) and DOCK6 ROC--AUC was 0.616 (34 scored). A subsequent audit showed that target-assigned IC50 records included both protein and cellular endpoints. We therefore separated assay descriptions into catalytic enzyme, biochemical binding, cellular and uncertain categories, kept explicit mutant records separate, and excluded all scaffolds from the previously revealed 36-structure set before comparing fixed fingerprint logistic regression with training-active chemical similarity.

In the catalytic stratum without an explicitly stated mutation, 345 active and 11 weak parents remained. Connected publication--scaffold five-fold evaluation yielded pooled out-of-fold ROC--AUC values of 0.918 for the model and 0.776 for similarity; model fold-average AUC was 0.816 and one fold was 0.480. Only six independent publication--scaffold components contained weak controls. The cellular signaling stratum could not support the same five-fold source-separated evaluation. Local design from 12 measured-active reference groups produced 69 analogues, none of which passed all three fixed predicted-toxicity gates. A targeted safety review did not establish adequate adjudicated DILI labels for evaluating those predictions. The study supports workflow integration and an exploratory enzyme-ranking signal; it establishes neither activity or safety of new molecules nor efficacy in PAH.

'''
s=s[:start]+abstract+s[end:]
replace('reinforcement learning; molecular docking; pharmacophore modeling; ADMET prediction; drug discovery pipeline; computational pharmacology','molecular docking; pharmacophore modeling; ADMET prediction; endpoint curation; grouped validation; computational pharmacology')
replace('The primary aim of this study is to develop an end-to-end integrated computational drug discovery pipeline capable of systematically identifying and generating novel therapeutic candidate molecules for PAH.','The primary aim of this study is to develop and critically evaluate an integrated computational workflow for generating and prioritizing molecular hypotheses around PAH-associated targets, with ALK5 as the principal retrospective evaluation target. The methodological evaluation tests whether endpoint-separated activity prediction adds ranking information beyond chemical similarity; it does not demonstrate new therapeutics for PAH.')
replace('thereby ensuring that the generated compounds remain meaningful within the biological context of the target binding pocket.','while allowing their compatibility with a target pocket to be examined as a hypothesis rather than assumed from generation alone.')
replace('ensures that candidates with unfavorable absorption, metabolism, or toxicity profiles are eliminated at early stages of the design.','applies explicit model-based filters to potential liabilities; these filters do not establish actual absorption, metabolism or toxicity and require endpoint-appropriate validation.')
replace('To reduce computational noise and improve biological interpretability, non-drug-like entities frequently present in crystallographic structures---including crystallographic water molecules, buffer components, metal ions, and solvent molecules---were systematically excluded.','For the initial pocket-residue reporting step, crystallographic water, buffer and solvent entities were excluded. This reporting convention must not be generalized to receptor preparation: functionally required ions, including the active-site zinc in the carbonic-anhydrase pose controls, were retained as documented for those experiments.')
replace('The resulting pharmacophore model serves as a structural constraint framework guiding the downstream de novo molecular generation process.','The resulting pharmacophore model provides a structure-derived interaction template for downstream inspection. The sampler used here was unconditional; the feature template was not an experimentally validated target-conditioning mechanism.')
replace('including valence checks, structural consistency verification, and elimination of chemically unstable structures prior to downstream evaluation.','including valence checks and structural consistency verification prior to downstream evaluation. These computational checks do not establish chemical stability or synthetic accessibility.')
replace('Early elimination of predicted mutagenic compounds significantly reduces the risk of downstream toxicological failure and aligns with regulatory toxicology screening practices.','The initial run used this prediction as a filtering convention; an actual reduction in toxicological failure was not measured. Later experiments additionally applied hERG and DILI model gates, whose limitations are examined below.')
replace('Independent-engine docking and biological assays remain incomplete and are not claimed.','Independent-engine calculations were subsequently completed as described in Sections 3.9--3.11; biological assays on the generated molecules remain absent.')
replace('orthogonal docking with an independent engine was not completed;','orthogonal docking was subsequently completed but did not establish reliable active--weak discrimination in the separate 36-structure test;')
replace('In particular, using Ames mutagenicity prediction as an early filtering criterion reduces potential safety risks by eliminating molecules with high genotoxic potential at an early stage.','In particular, the Ames prediction was used to flag potential genotoxic liability. Neither that prediction nor the later hERG/DILI gates establishes an actual reduction in safety risk; none of the 69 locally designed analogues passed all three fixed gates.')
replace('the candidate ranking remains provisional: exploratory DOCK6 analysis showed weak cross-engine agreement, and experimental measurements remain absent.','the candidate ranking remains provisional: agreement varied between candidate sets, the separate 36-structure ranking test did not establish reliable discrimination, and experimental measurements on the designed molecules remain absent.')
# Preserve original figure files and their contents, while separating the conceptual graphic from numbered results.
a=s.index('\\section{\\texorpdfstring{\\protect\\includegraphics')
b=s.index('\\includegraphics[width=6.4in',a)
s=s[:a]+r'''\begin{center}
\includegraphics[width=0.83\linewidth,height=0.80\textheight,keepaspectratio]{media/image1.png}
\end{center}
\noindent\emph{Graphical overview retained from the original manuscript. It depicts the conceptual pipeline and PAH motivation; it is not evidence of target-conditioned sampling, experimental safety or therapeutic efficacy. The measured implementation and later validation limits are described in the text.}

'''+s[b:]
methods=r'''\subsection{2.6 Endpoint-Aware Activity Data Curation}
The later evaluation used archived ChEMBL 37 activity, assay and document records for human ALK5/TGFBR1 (CHEMBL4439). B-type assay annotation and target-mapping confidence 9 were retained as eligibility filters, but were not treated as sufficient evidence of a biochemical endpoint. Full assay descriptions were separated into catalytic enzyme activity, biochemical binding, protein-based measurements with an incompletely specified mechanism, cellular signaling, cellular target engagement, cellular proliferation and uncertain categories. Protein production in Sf9 cells was distinguished from a readout measured in intact cells. HTRF or TR-FRET detection alone was not assumed to establish catalytic activity. Explicit T204D and other mutant descriptions were kept separate; absence of a stated mutation was not interpreted as confirmed wild-type protein. Classification decisions and their rationales are archived and remain subject to original-protocol verification.

Chemical identities underwent fragment-parent normalization, neutralization and canonical tautomer selection. Stereochemistry was removed for grouping and radius-2, 2048-bit Morgan fingerprints, limiting the ability to resolve stereospecific effects. Measurements were pooled only within the same endpoint and explicit-mutation stratum. Active parents had median exact pChEMBL at least 7; weak parents had median exact pChEMBL at most 5, or qualifying measured IC50 right-censored lower bounds of at least 10,000 nM without contradictory exact results. Censored measurements were not assigned invented point IC50 values. Intermediate and conflicting groups were excluded with an audit trail.

\subsection{2.7 Fixed Model, Similarity Reference and Grouped Evaluation}
All scaffolds in the previously revealed 36-structure docking test were excluded from the subsequent model comparison. A fixed balanced logistic regression model (C=1; radius-2, 2048-bit Morgan fingerprints) was compared with maximum Tanimoto similarity to training-fold actives. No hyperparameter search was performed. Five-fold scaffold grouping was accompanied, where feasible, by connected publication--scaffold grouping: molecules sharing a document or an exact scaffold remained in the same component. The random seed was 20261016. This is retrospective exploratory reanalysis, not a new blinded trial or a prospective test of designed molecules.

Pooled and fold-level ROC--AUC, average precision with class prevalence, and group-bootstrap uncertainty from fixed out-of-fold predictions were reported. Those intervals do not include full model-retraining or assay-classification uncertainty. Pooled decision scores may differ in scale between folds, so per-fold AUC and weak-control counts were retained. Strata with too few independent weak groups were not forced into the specified five-fold source-separated test. Model scores are not calibrated activity probabilities or quantitative IC50 predictions.

\subsection{2.8 Local Analogue Design and Safety Evidence Review}
A separate design exercise selected up to 12 measured-active parent groups with distinct scaffolds and representative experimental structures with no unresolved potential stereochemistry. Single aromatic C--H to F additions or terminal aromatic F/Cl substitutions retained the parent scaffold and at least 0.65 Morgan similarity. Previously archived experimental and generated structures were excluded after normalization. Physicochemical, PAINS, synthetic-accessibility heuristic and stereo checks preceded ADMET prediction; no synthetic route was validated. Up to six analogues per reference were selected before ADMET inspection. The previously fixed AMES, hERG and DILI gates (each below 0.5) were retained without relaxation.

For the 12 references, a targeted review of identities, original discovery sources and available safety reports separated human, animal, cell-based and ADME evidence. A negative safety report, missing data, metabolic stability or activity in a liver-related disease model was not automatically converted into a DILI-negative label. Reference-group median potency and the representative structure's individual measurement were recorded separately and were never assigned to a newly designed analogue.

'''
anchor=s.index('\\begin{center}\n\\includegraphics[width=0.83')
s=s[:anchor]+methods+s[anchor:]
results=r'''\subsection{3.11 Separate Label-Withheld Docking Recognition Test}
The original exact-IC50-only design yielded no new eligible weak controls after historical parent and scaffold exclusions. Before scoring, a documented protocol amendment admitted measured right-censored IC50 bounds of at least 10,000 nM while retaining the potency thresholds and exclusion rules. The resulting set contained 18 active and 18 weak structures on 36 distinct scaffolds; all 18 weak labels used lower-bound evidence. Randomized case identifiers and structures were supplied to scoring, with a separate hashed label key. This was code-level label withholding rather than a third-party administered blind trial. Assay heterogeneity was subsequently identified and limits interpretation of this historical benchmark.

Vina scored 35 of 36 structures and produced ROC--AUC 0.441 (fixed matched-pair-bootstrap 95\% interval 0.274--0.611), with three actives among the first ten. DOCK6 scored 34 structures and produced AUC 0.616 (0.481--0.758), with seven actives among the first ten. Secondary mean-rank consensus found five actives among the first ten. Both primary intervals included random discrimination. DOCK6 failures were preserved; Vina was attempted independently with unchanged inputs on those cases, yielding one score and one 900-second timeout. These results did not demonstrate reliable active--weak discrimination or validate new molecules. Neither score direction nor membership was changed to improve the outcome.

\subsection{3.12 Assay Endpoint Audit and Dataset Separation}
The audit classified 153 eligible assay descriptions: 50 catalytic enzyme, 23 biochemical binding, 20 protein-based but incompletely specified, 28 cellular signaling, four cellular target engagement, six proliferation and 22 uncertain. Earlier statements implying that B-type annotation and confidence 9 guaranteed biochemical-only observations are therefore corrected. Earlier numeric results remain archived as mixed-endpoint exploratory analyses. Direct binding, catalytic inhibition and cellular signal responses were not pooled in the revised comparison.

\begin{figure}[htbp]
\centering\includegraphics[width=\linewidth]{figures/figure1_workflow.png}
\caption{Figure 6. Endpoint-aware evaluation added to the original AIPHARMA workflow. Assay classification, explicit-mutant separation and publication--scaffold grouping precede interpretation of fixed-model results. This is a retrospective evaluation, not experimental validation of new candidates.}
\end{figure}
\begin{figure}[htbp]
\centering\includegraphics[width=0.95\linewidth]{figures/figure2_endpoint_counts.png}
\caption{Figure 7. Description-based endpoint counts for 153 ALK5 assays. Counts are assays, not compounds. Original-source protocol verification remains necessary.}
\end{figure}

\subsection{3.13 Endpoint-Separated Model and Similarity Comparison}
After exclusion of previously revealed holdout scaffolds, the catalytic stratum without an explicitly stated mutation contained 345 active and 11 weak normalized parents. The cellular signaling stratum contained 90 active and 19 weak parents. Separate models were fitted within each stratum; one stratum did not supply training observations to the other.

\begin{center}
\begin{tabular}{lrr}
\hline
Stratum and split & Model AUC & Similarity AUC\\
\hline
Enzyme, scaffold grouped & 0.983 & 0.883\\
Enzyme, document--scaffold grouped & 0.918 & 0.776\\
Cell signaling, scaffold grouped & 0.990 & 0.972\\
Cell signaling, document--scaffold grouped & Infeasible & Infeasible\\
\hline
\end{tabular}
\end{center}

In the enzyme document--scaffold analysis, fixed-prediction group-bootstrap 95\% intervals were 0.707--0.984 for the model and 0.400--0.936 for similarity. The exploratory AUC-difference interval was 0.014--0.425. These are conditional uncertainty summaries from fixed out-of-fold predictions, not full model-development uncertainty. Model fold AUCs were 1.000, 0.912, 0.808, 0.480 and 0.882 (mean 0.816); similarity fold-average AUC was 0.662. Only six independent document--scaffold components contained weak enzyme controls, with one to five weak observations per test fold. With this class balance, random selection would find approximately 9.69 actives among ten; ten-of-ten retrieval is not a strong success metric here.

The cellular signaling weak controls occupied only two independent source--scaffold components, so the specified source-separated five-fold evaluation was infeasible. Its high scaffold-only AUC must not be interpreted as source-independent validation. Assay-source and class-count limitations remain important even after endpoint separation. Different cohorts prevent direct comparison with earlier mixed-endpoint or docking AUC values as proof that curation alone improved accuracy.

\begin{figure}[htbp]
\centering\includegraphics[width=\linewidth]{figures/figure3_enzyme_performance.png}
\caption{Figure 8. Catalytic-enzyme comparison under document--scaffold grouping. Left: pooled out-of-fold AUC and fixed-prediction group-bootstrap intervals. Right: fold AUC with weak-control counts. Intervals omit retraining uncertainty; the ``no explicitly stated mutation'' stratum does not establish wild-type protein.}
\end{figure}

\subsection{3.14 Chemical Coverage and Local Analogue Follow-Up}
For the previous 20 generated structures, nearest-training Morgan similarities were 0.143--0.394, below the 0.534 median observed for source--scaffold-separated reference test molecules. Nine candidates lay below the fifth percentile of that reference similarity distribution. This is a chemical-coverage warning, not a calibrated prediction-error guarantee, and does not establish failure or success of any compound.

Local editing of 12 measured-active reference groups produced 69 preselected analogues. Four known or duplicate proposed structures were excluded before that screen. Of the 69 analogues, 38 passed the AMES gate, ten passed hERG and none passed DILI; no structure passed all three. None of the 12 parent references passed all three gates either. Thresholds were not relaxed and no jointly filtered series was nominated. Parent activity was not assigned to any analogue, and local atom edits plus an SA heuristic were not treated as validated synthesis routes or global novelty evidence.

\subsection{3.15 Targeted DILI Evidence Audit}
Reference identities and assay/document provenance were retained for all 12 parents. CHEMBL2364611 was identified as galunisertib and CHEMBL260015 as LY2109761. Published combination studies of galunisertib provided mixed liver-related findings: a gemcitabine-controlled study reported no apparent additional liver-toxicity difference, whereas a durvalumab combination reported liver-enzyme elevations. Combination treatment and clinical context prevent these observations from supplying a clean single-compound DILI-positive or DILI-negative label. Preclinical antifibrotic activity, metabolic stability, CYP profiling and absence of weight loss were not equated with clinical liver safety.

No adequate adjudicated binary DILI labels were established for estimating model accuracy on these 12 references. Their DILI model scores (0.841--0.999) therefore cannot establish either actual liver toxicity or model false-positive rates. Missing safety data remained unknown. The safety search was targeted rather than systematic, some evidence was abstract-only, and results from a lead compound were not transferred to other molecules solely because they appeared in the same publication. The zero-of-69 prediction-gate result remains archived without retrospective threshold adjustment.

'''
anchor=s.index('\\textbf{4. Discussion}')
s=s[:anchor]+results+s[anchor:]
newdiscussion=r'''\subsection{4.4 Endpoint-Aware Evidence and Interpretation}
The principal new methodological finding is that target assignment alone did not provide a homogeneous activity endpoint. The revised enzyme analysis suggests a ranking contribution from the fixed fingerprint model beyond a chemical-similarity reference, but this is exploratory evidence from few weak controls. One source-separated fold remained near random performance, and high pooled AUC should be read alongside fold-level variation. Explicit mutant annotation, ATP conditions, assay format and publication context may still vary within the retained stratum. Original assay protocols require further review, and additional independent weak and intermediate-activity compounds are needed.

This evidence does not retrospectively validate the unconditional generator, the previously selected 20 structures or the 69 local analogues. A model that separates strong actives from very weak controls need not provide accurate quantitative potency, activity around a decision threshold, selectivity or reliable predictions outside its chemical coverage. The DILI evidence audit also prevents interpreting model rejections as measured toxicity. A defensible contribution is a reproducible workflow with explicit failures, endpoint-aware evaluation and conditional uncertainty, rather than a validated therapeutic discovery claim.

'''
s=s.replace('\\section{5. Conclusion}',newdiscussion+'\\section{5. Conclusion}',1)
replace('Its principal contribution is a modular computational workflow accompanied by explicit limits on interpretation and a predefined route toward orthogonal and experimental validation.','Its contribution is a modular workflow whose structural and generative experiments are retained alongside negative docking evidence, an endpoint-separated model comparison and a targeted safety-evidence audit. The enzyme model showed an exploratory ranking advantage over similarity under publication--scaffold grouping, but sparse weak controls, fold variation and unverified assay details restrict the inference. No designed molecule was established as active, safe, synthesizable or therapeutically effective.')
replace('This study proposes a novel computational approach that integrates structure-based drug design and AI-based molecule generation methods for the discovery of small-molecule therapeutics targeting pulmonary arterial hypertension (PAH).','This study implements and evaluates an integrated structure-based and generative workflow motivated by PAH-associated targets, with ALK5 used for retrospective methodological evaluation. It does not establish the discovery of a therapeutic molecule.')
replace('The next validation stage will build on the completed ten-seed, ALK5 active--control and exploratory twelve-candidate DOCK6 analyses.','The next validation stage will build on the completed ten-seed sampling, independent-engine studies, separate docking recognition test, endpoint-separated model comparison and targeted safety audit. Priorities are original-source assay review, additional independent weak and intermediate enzyme controls, a time-aware split and a new test locked before further model development. The earlier revealed tests are development history, not untouched validation.')
replace('The pipeline is implemented as a single reproducible script executed in a Linux/conda environment integrating fpocket, RDKit, REINVENT, ADMET-AI and AutoDock Vina.','The pipeline and subsequent evaluations are implemented as archived scripts in Linux/conda environments integrating fpocket, RDKit, REINVENT, ADMET-AI, AutoDock Vina, DOCK6 and scikit-learn. The source package contains the original six illustrations, three endpoint-evaluation figures, generating code and machine-readable evidence packages. Analysis scripts retain documented original data-root assumptions and are not a fully installed portable software environment.')
# Preserve every original reference; add verified source pointers for newly described analyses.
refs=r'''\subsection*{Additional sources for endpoint evaluation and safety context}
Pedregosa, F., et al. (2011). Scikit-learn: Machine Learning in Python. \emph{Journal of Machine Learning Research, 12}, 2825--2830. \url{https://www.jmlr.org/papers/v12/pedregosa11a.html}.

The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Database background: \url{https://pubmed.ncbi.nlm.nih.gov/37933841/}. The present analyses use archived release 37 records, not the release described in that article.

Therapeutics Data Commons. DILI task description: \url{https://tdcommons.ai/single_pred_tasks/tox/}. Accessed 16 October 2026.

Galunisertib plus gemcitabine vs. gemcitabine for first-line treatment of patients with unresectable pancreatic cancer. Clinical context: \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC6251034/}.

Safety and activity of the TGFbeta receptor I kinase inhibitor galunisertib plus the anti-PD-L1 antibody durvalumab in metastatic pancreatic cancer. Clinical context: \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC7944986/}.

\noindent\textbf{Editorial status.} This full-length working manuscript preserves the original literature and figures. Existing references are retained, not newly certified by this revision; original-source assay verification, complete reference checking, authorship/disclosures, public data deposition and final compiled layout review remain necessary before submission.

'''
s=s.replace('\\end{document}',refs+'\\end{document}')
# Source-level preservation checks.
oldrefs=old[old.index('\\section{References}'):old.index('\\end{document}')]
assert oldrefs in s
for i in range(1,7):assert '{media/image'+str(i)+'.png}' in s
for f in ['figure1_workflow.png','figure2_endpoint_counts.png','figure3_enzyme_performance.png']:assert '{figures/'+f+'}' in s
for heading in ['1. Introduction','2.3 De Novo Molecular Design','3.10 ADMET-First','4.2 Comparison','5.2 Future Work']:
 assert heading in s
assert len(s)>len(old)
(p/'restored.tex').write_text(s,encoding='utf-8')
audit={'old_characters':len(old),'short_draft_characters':len(short),'restored_characters':len(s),'all_original_reference_text_preserved':True,'all_six_original_images_referenced':True,'three_new_figures_added':True,'new_methods':['2.6','2.7','2.8'],'new_results':['3.11','3.12','3.13','3.14','3.15'],'new_discussion':'4.4','targeted_replacements':changes,'original_sha256':hashlib.sha256(old.encode()).hexdigest(),'restored_sha256':hashlib.sha256(s.encode()).hexdigest()}
(p/'restoration_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps({k:v for k,v in audit.items() if k!='targeted_replacements'},indent=2))
