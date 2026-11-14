from pathlib import Path
import re,json
p=Path.home()/'restore_full_alk5';s=(p/'tool_input.tex').read_text(encoding='utf-8-sig');before=s
s=s.replace('AIPHARMA: An Integrated Computational Drug Discovery Workflow with Endpoint-Aware Retrospective Evaluation of ALK5 Activity','AIPHARMA: An Integrated Tool for Molecular Design and Prioritization')
a=s.index('Pulmonary arterial hypertension (PAH) motivates');b=s.index('\\textbf{Keywords:}',a)
s=s[:a]+r'''AIPHARMA is an integrated computational tool for molecular design and prioritization that combines protein input and pocket analysis, pharmacophore characterization, molecular sampling, predicted ADMET profiling and docking-based evaluation in a single workflow. The principal contribution is the implementation of this connected workflow and the traceable outputs it produces. ALK5, with PAH as the original application motivation, serves as a case study rather than evidence of therapeutic efficacy. Molecular sampling used an unconditional pretrained prior; target-specific information entered through downstream evaluation.

We evaluated tool outputs at three distinct levels: archived-output integrity, agreement with known experimental reference data, and the status of validation for newly designed molecules. Technical checks verified the unchanged locked selection, 20 unique structures, consistency with the selected prediction thresholds and the absence of source/scaffold overlap in the specified model evaluation. Repeated property prediction on an earlier 12-structure set agreed within approximately $1.2\times10^{-7}$, while a broader stereochemistry audit exposed unresolved units in four locked structures. These checks are not a same-input rerun of the entire tool. Internal ALK5 pose-recovery controls succeeded, but a separate 36-structure docking recognition test did not establish reliable discrimination (Vina AUC 0.441; DOCK6 AUC 0.616). Endpoint-separated activity modeling yielded an exploratory signal relative to chemical similarity, limited by only 11 weak enzyme controls and variable fold performance.

The case study demonstrates connected computation and auditable candidate outputs while exposing limitations of scoring, endpoint pooling and predicted-toxicity gates. None of 69 locally designed analogues passed all three fixed ADMET gates, and adequate reference DILI labels were not established for calibrating those decisions. The tool therefore supports generation and examination of molecular hypotheses; it does not establish biological activity, experimental safety or efficacy of new molecules. Output validation is reported by component and evidence level rather than as a single accuracy claim for the whole platform.

'''+s[b:]
s=s.replace('The methodological evaluation tests whether endpoint-separated activity prediction adds ranking information beyond chemical similarity; it does not demonstrate new therapeutics for PAH.','The central contribution is the AIPHARMA tool: connected computational modules, inspectable intermediate outputs and explicit candidate-selection records. Component-level validation and comparison with known experimental reference data assess the outputs without claiming new therapeutics for PAH. Endpoint-separated modeling is one evaluation of the output-prioritization problem, not a replacement for the tool contribution.')
s=s.replace('\\subsection{2.6 Endpoint-Aware Activity Data Curation}','\\subsection{2.8 Endpoint-Aware Activity Data Curation}').replace('\\subsection{2.7 Fixed Model, Similarity Reference and Grouped Evaluation}','\\subsection{2.9 Fixed Model, Similarity Reference and Grouped Evaluation}').replace('\\subsection{2.8 Local Analogue Design and Safety Evidence Review}','\\subsection{2.10 Local Analogue Design and Safety Evidence Review}')
methods=r'''\subsection{2.11 Validation of the Tool and Its Outputs}
Validation was organized by the claim being tested. Technical output checks examined the identity and immutability of the locked selection, unique molecular records, consistency between predictions and selected thresholds, candidate-count reconciliation, and group-separation records. A new audit recomputed the archived locked-file hash, uniqueness and gate checks and reconciled the enzyme evaluation counts and reported split overlaps. Earlier repeat-prediction and conformer/stereo audits were retained with their original scope. File integrity and repeated prediction do not establish full-workflow repeatability, correctness of every module, or user-interface usability.

Reference-based scientific evaluation examined crystallographic pose recovery, ranking of measured active/weak compounds, comparison between scoring engines and a simple chemical-similarity baseline, and the adequacy of assay endpoint definitions. A successful pose control tests pose recovery under its protocol, not activity ranking or affinity. The separate fingerprint-model comparison is an additional retrospective prioritization analysis; its AUC is not the accuracy of the complete AIPHARMA tool.

Validation of newly designed outputs would require independent compound measurements. No experimental activity or safety results for those outputs were available. Full same-input end-to-end reruns, controlled run-time measurements, usability evaluation and comparison with alternative integrated tools were not completed by the archived-output audit and remain distinct technical evaluation requirements. The validation summary below separates established checks, limited scientific support and missing evidence.

'''
pos=s.index('\\begin{center}\n\\includegraphics[width=0.83')
s=s[:pos]+methods+s[pos:]
# Put all detailed output evaluations beneath one validation umbrella, without deleting them.
for old in range(6,16):
 pattern=r'\\(?:textbf|subsection)\{3\.'+str(old)+r' ([^}]+)\}'
 s,n=re.subn(pattern,lambda m:'\\subsubsection{3.6.'+str(old-5)+' '+m.group(1)+'}',s)
 assert n==1,(old,n)
summary=r'''\subsection{3.6 Validation of the Tool and Its Outputs}
Table~\ref{tab:tool-validation} summarizes what the current evidence validates and what it does not. The detailed experiments follow within this section. The original structural and generation outputs are retained above as the ALK5 workflow case study. The validation evidence is intentionally component-specific: technical traceability, reference-data performance and biological confirmation are separate claims.

\begin{table}[htbp]
\centering\small
\begin{tabular}{p{0.19\linewidth}p{0.42\linewidth}p{0.25\linewidth}}
\hline
Evidence level & Completed observation & Interpretation\\
\hline
Output integrity & Locked 20-structure hash, uniqueness and selected model gates rechecked; enzyme counts and recorded source/scaffold separation reconciled & Supports consistency of archived outputs; not a full pipeline rerun\\
Repeat prediction & Earlier 12-structure predictor outputs differed by at most about $1.2\times10^{-7}$ & Supports that repeat-prediction check only\\
Stereo audit & Four locked structures contained unresolved potential stereogenic units missed by the initial screen & Exposes a limitation; complete stereo resolution was not achieved\\
Reference poses & Internal DOCK6 ALK5 pose controls: 1.392 and 0.328 A RMSD & Pose recovery for development controls; not affinity or activity validation\\
Measured activities & Separate 36-structure docking set: Vina AUC 0.441, DOCK6 0.616 & Reliable docking discrimination not established\\
Activity-model comparator & Enzyme model AUC 0.918 vs similarity 0.776; 345 active / 11 weak, variable fold results & Exploratory component analysis, not whole-tool accuracy\\
Designed molecules & No new-compound assays; 0/69 local analogues passed all three prediction gates & Activity and safety unvalidated; prediction rejection is not measured toxicity\\
\hline
\end{tabular}
\caption{Evidence map for AIPHARMA and its outputs. No row establishes experimental efficacy or safety of a designed molecule.}
\label{tab:tool-validation}
\end{table}

'''
needle='\\subsubsection{3.6.1 Redocking Pose-Recovery Results}'
assert needle in s;s=s.replace(needle,summary+needle)
s=s.replace('Sections 3.9--3.15','Section 3.6').replace('Sections 3.9--3.11','Section 3.6')
s=s.replace('The principal new methodological finding is that target assignment alone did not provide a homogeneous activity endpoint.','The principal software contribution remains the connected AIPHARMA workflow. Its output evaluation additionally showed that target assignment alone did not provide a homogeneous activity endpoint. This observation explains why component-level validation requires more than running the modules successfully.')
s=s.replace('Its contribution is a modular workflow whose structural and generative experiments are retained alongside negative docking evidence, an endpoint-separated model comparison and a targeted safety-evidence audit.','Its contribution is an integrated molecular-design and prioritization tool with preserved structural/generative case-study outputs, explicit selection records and a component-level validation evidence map. Negative docking evidence, endpoint-separated model comparisons and safety-evidence auditing characterize its current limits rather than replace the software contribution.')
s=s.replace('The work supports further method validation and hypothesis generation.','The work supports further tool evaluation and hypothesis generation.')
# Use explicit mathematical minus signs instead of Unicode text minus in the reported range.
s=s.replace('ranged from −7.8 to −10.2 kcal/mol; five candidates scored at or below −9.4 kcal/mol','ranged from $-10.2$ to $-7.8$ kcal/mol; five candidates scored at or below $-9.4$ kcal/mol')
s=s.replace('Figure 1. The AIPHARMA pipeline: target structure -\\textgreater{} fpocket pocket detection -\\textgreater{} pharmacophore modeling -\\textgreater{} REINVENT de novo design -\\textgreater{} ADMET-AI filtering -\\textgreater{} AutoDock Vina docking and ranking.','Figure 1. Conceptual AIPHARMA module overview. The measured sampler is an unconditional REINVENT prior; pocket/pharmacophore information supplies downstream evaluation context rather than conditioning generation. ADMET scores and docking ranks are computed outputs, not validated activity or safety.')
for i in range(1,7):assert '{media/image'+str(i)+'.png}' in s
assert len(re.findall(r'\\includegraphics',s))==9
assert '3.6.10 Targeted DILI Evidence Audit' in s
assert s.count('\\subsection{2.8 ')==1 and s.count('\\subsection{2.9 ')==1
assert '\\section{References}' in s
(p/'tool_updated.tex').write_text(s,encoding='utf-8');print(json.dumps({'characters_before':len(before),'characters_after':len(s),'figures_preserved':9,'validation_umbrella':'3.6','detailed_validation_subsections':10,'whole_tool_accuracy_claim':False}))
