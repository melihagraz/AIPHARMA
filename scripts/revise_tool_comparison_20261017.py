from pathlib import Path
import json,re
p=Path.home()/'restore_full_alk5';s=(p/'comparison_input.tex').read_text(encoding='utf-8-sig');before=s
# Correct three specific contradictions without removing the original background.
a=s.index('This study proposes an end-to-end, integrated, and modular computational drug discovery pipeline');b=s.index('\n\n',a)
s=s[:a]+'''This study presents AIPHARMA as an integrated tool for molecular design and prioritization, connecting structural input, pocket and pharmacophore characterization, unconditional molecular sampling, predicted properties and docking outputs. Its contribution is a concrete workflow implementation with inspectable results and component-specific evaluation, not a new docking algorithm, a demonstration of target-conditioned generation or a solution to protein flexibility. PAH provides the original application motivation and ALK5 the principal worked example. The following sections describe the tool modules, distinguish historical experiments from completed validation checks, and report both successful integration and limitations of the resulting prioritization. Reinforcement-learning approaches are background and possible extensions; they were not trained or evaluated as a comparative arm in the reported tool experiments.'''+s[b:]
a=s.index('In this context, AIPHARMA offers an end-to-end computational drug discovery pipeline');b=s.index('\\section{4.1 Pipeline Performance Evaluation}',a)
s=s[:a]+'''AIPHARMA connects structural characterization with downstream evaluation of molecules sampled from an unconditional prior. Protein-pocket information defines the setting for inspection and docking, rather than directly conditioning the generation step. The case study therefore evaluates how these modules operate together and what can be inferred from their outputs, without attributing target-specific learning to the sampler.

Pocket descriptors provide a practical interface between a supplied protein structure and subsequent prioritization (Schmidtke \\& Barril, 2010). They do not remove the rigid-receptor approximation or establish that the selected structures are biologically appropriate. The original literature discussion of protein dynamics motivates these limitations; no enhanced dynamics model, molecular-dynamics validation or general resolution of static-structure limitations is claimed here.

'''+s[b:]
# Explicitly distinguish the two historical Vina cohorts/settings.
anchor='\\subsubsection{3.6.7 Assay Endpoint Audit and Dataset Separation}'
para=r'''\paragraph{Why the two Vina benchmark estimates differ.}
The earlier AUC of 0.690 and later AUC of 0.441 were obtained on different cohorts and protocols. The first analysis began with 70 active and 70 exact-value weak controls, matched on physicochemical descriptors, and scored 139 structures at exhaustiveness 8 after one preparation failure. The separate test excluded previously recorded parents and scaffolds, retained property matching but could form only 18 pairs, and admitted measured right-censored weak bounds after the exact-only design proved infeasible. It scored 35 structures with a different documented ligand-preparation protocol and exhaustiveness 16. Assay heterogeneity was subsequently recognized in the activity labels. These estimates are neither repeated measurements of one fixed benchmark nor a controlled test of the effect of scaffold exclusion, censoring or preparation alone. The observed difference cannot be attributed to any one factor and is not evidence that the tool improved or deteriorated between runs.

'''
assert anchor in s;s=s.replace(anchor,para+anchor)
a=s.index('\\section{4.2 Comparison with Existing Approaches}');b=s.index('\\textbf{4.3 Strengths and Limitations}',a)
comparison=r'''\section{4.2 Comparison with Existing Approaches}\label{comparison-with-existing-approaches}
Integrating molecular generation with docking is established practice rather than a novelty claim for AIPHARMA. DockStream provides ligand preparation and multiple docking backends, automated benchmarking/analysis and a REINVENT scoring interface. Icolos provides configurable computational-chemistry workflows that can run independently or alongside REINVENT and include docking, molecular dynamics and quantum-chemical steps. Their developer documentation is cited below. Neither package was executed in a matched head-to-head experiment in this study.

AIPHARMA's present contribution is the particular connected workflow and its documented case study: structural input and pocket characterization, an unconditional sampler, property profiling, candidate records and downstream scoring, accompanied by an evidence map and same-host CLI repeat checks. These are implementation characteristics, not proof of a capability absent from other tools. The endpoint-separated ALK5 analysis further illustrates why execution success must be distinguished from scientific validity of a ranking.

\begin{table}[htbp]
\centering\small
\begin{tabular}{p{0.17\linewidth}p{0.43\linewidth}p{0.26\linewidth}}
\hline
Tool & Documented role & Evidence used here\\
\hline
DockStream & Ligand-preparation/docking wrapper; REINVENT docking-score integration & Developer documentation; not run here\\
Icolos & Configurable structure-based workflow manager; standalone and REINVENT-linked use & Developer documentation; not run here\\
AIPHARMA & Integrated pocket-to-candidate case-study workflow with property and docking outputs & Local source/output audits and two small CLI runs\\
\hline
\end{tabular}
\caption{Scope comparison based on the cited developer documentation and local AIPHARMA evidence. This is not a speed, predictive-accuracy, usability or feature-completeness benchmark. Lack of a feature in this table does not indicate its absence from a tool.}
\label{tab:related-tools}
\end{table}

No faster execution, better molecular activity, easier installation or superior usability relative to these systems is established. The two-run timing table measures only the specified AIPHARMA example. A future matched comparison would require aligned inputs, preparation, scoring backends, resources and failure accounting; it would still assess a defined task rather than all tool functionality. The developer repositories were consulted on 17 October 2026; this comparison is not a comprehensive survey of current software.

'''
s=s[:a]+comparison+s[b:]
s=s.replace('Activity-model comparator & Enzyme model AUC 0.918 vs similarity 0.776; 345 active / 11 weak, variable fold results','Activity-model comparator & Enzyme AUC 0.918 vs 0.776; 345 active / 11 weak; one fold AUC 0.480')
# Add traceable official sources without inventing bibliographic metadata.
refs=r'''MolecularAI. DockStream developer documentation: \url{https://github.com/MolecularAI/DockStream/blob/master/README.md}. Accessed 17 October 2026. Describes the docking wrapper, backends and REINVENT interface.

MolecularAI. Icolos developer documentation: \url{https://github.com/MolecularAI/Icolos}. Accessed 17 October 2026. Describes configurable structure-based workflows and standalone/REINVENT use.

'''
s=s.replace('\\noindent\\textbf{Editorial status.}',refs+'\\noindent\\textbf{Editorial status.}')
assert len(re.findall(r'\\includegraphics',s))==9
assert 'directly relating the molecule generation process to the structural properties' not in s
assert 'aims to overcome the classical assumptions' not in s
(p/'comparison_updated.tex').write_text(s,encoding='utf-8');print(json.dumps({'characters_before':len(before),'characters_after':len(s),'figures_retained':9,'new_comparison':'documentation-based; no head-to-head performance claim'}))
