from pathlib import Path
import re,unicodedata,json
p=Path.home()/'restore_full_alk5';s=unicodedata.normalize('NFC',(p/'editorial_input.tex').read_text(encoding='utf-8-sig'));before=s
s=s.replace('\\usepackage{graphicx}','\\usepackage{graphicx}\n\\usepackage{array,longtable}\n\\usepackage[margin=1in]{geometry}')
s=s.replace('≥',r'\ensuremath{\geq}').replace('≤',r'\ensuremath{\leq}')
s=s.replace(r'10\textsuperscript{60}',r'$10^{60}$').replace(r'10\textsuperscript{100}',r'$10^{100}$').replace(r'10\^{}60',r'$10^{60}$')
s=s.replace('Å²',r'\ensuremath{\mathrm{\AA}^{2}}').replace('Å³',r'\ensuremath{\mathrm{\AA}^{3}}').replace('Å',r'\AA{}')
s=s.replace('9729--97451','9729--9747').replace('10.1186/s13321-017-0235-}', '10.1186/s13321-017-0235-x}')
s=s.replace('(Pala et al., 2024; Zhang et al., 2024)','(Lexa \\& Carlson, 2012; Pantsar \\& Poso, 2018)')
s=s.replace('Langer, T., \\& Wolber, G. (2004). Pharmacophore definition and 3D searches. \\emph{Drug Discovery Today.}',r'Langer, T., \& Wolber, G. (2004). Pharmacophore definition and 3D searches. \emph{Drug Discovery Today: Technologies, 1}(3), 203--207. \url{https://doi.org/10.1016/j.ddtec.2004.11.015}.')
# Preserve the conceptual artwork as a supplementary item, avoiding two main workflow overviews.
pat=r'\\begin\{center\}\s*\\includegraphics\[width=0\.83\\linewidth,height=0\.80\\textheight,keepaspectratio\]\{media/image1\.png\}\s*\\end\{center\}\s*\\noindent\\emph\{Graphical overview retained.*?\n\n'
m=re.search(pat,s,re.S);assert m,'Conceptual artwork block not located';s=s[:m.start()]+s[m.end():]
extra=r'''\clearpage
\section*{Supplementary Conceptual Illustration}
\begin{center}
\includegraphics[width=0.83\linewidth,height=0.76\textheight,keepaspectratio]{media/image1.png}
\end{center}
\noindent\emph{Supplementary illustration S1. Original conceptual artwork, preserved for context. The implemented sampler is unconditional; target information enters downstream evaluation. Any graphical implication of directly conditioned generation or established activity/safety is not a validated feature. The main module overview and actual implementation are described in Figure 1 and the Methods.}
\clearpage

'''
s=s.replace('\\section{References}',extra+'\\section{References}',1)
# Convert manual captions to real floats, in source order (new figures already begin at counter5).
for idx,label,width in [(2,'pipeline',.95),(3,'pockets',.72),(4,'candidate-docking',.95),(5,'pose-overlay',.65),(6,'pose-controls',.95)]:
 pattern=r'\\includegraphics\[[^\]]+\]\{media/image'+str(idx)+r'\.png\}\s*\\emph\{Figure \d+\. (.*?)\}\s*(?=\n)'
 def repl(m):return '\\begin{figure}[htbp]\n\\centering\n\\includegraphics[width='+str(width)+'\\linewidth,keepaspectratio]{media/image'+str(idx)+'.png}\n\\caption{'+m.group(1)+'}\n\\label{fig:'+label+'}\n\\end{figure}\n'
 s,n=re.subn(pattern,repl,s,flags=re.S);assert n==1,(idx,n)
s=s.replace('(Figures 3 and 4)',r'(Figures~\ref{fig:pose-overlay} and~\ref{fig:pose-controls})')
# Long evidence map can break across pages instead of overflowing a float.
pattern=r'\\begin\{table\}\[htbp\]\s*\\centering\\small\s*\\begin\{tabular\}(.*?)\\end\{tabular\}\s*\\caption\{Evidence map(.*?)\\label\{tab:tool-validation\}\s*\\end\{table\}'
m=re.search(pattern,s,re.S);assert m,'Evidence table not found'
body=m.group(1)
body=body.replace('p{0.19\\linewidth}p{0.42\\linewidth}p{0.25\\linewidth}',r'>{\raggedright\arraybackslash}p{0.19\linewidth}>{\raggedright\arraybackslash}p{0.42\linewidth}>{\raggedright\arraybackslash}p{0.25\linewidth}')
replacement='\\begingroup\n\\small\n\\begin{longtable}'+body+'\\caption{Evidence map'+m.group(2)+'\\label{tab:tool-validation}\\\\\n\\end{longtable}\n\\endgroup'
s=s[:m.start()]+replacement+s[m.end():]
refs=r'''Berman, H. M., et al. (2000). The Protein Data Bank. \emph{Nucleic Acids Research, 28}(1), 235--242. \url{https://doi.org/10.1093/nar/28.1.235}.

Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., \& Beckstein, O. (2011). MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. \emph{Journal of Computational Chemistry, 32}, 2319--2327. \url{https://doi.org/10.1002/jcc.21787}.

Swanson, K., et al. (2024). ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical libraries. \emph{Bioinformatics, 40}(7), btae416. \url{https://doi.org/10.1093/bioinformatics/btae416}.

'''
s=s.replace('\\subsection*{Additional sources for endpoint evaluation and safety context}',refs+'\\subsection*{Additional sources for endpoint evaluation and safety context}')
assert s.count('\\includegraphics')==9
assert s.count('\\begin{figure}')==8
assert s.count('\\begin{longtable}')==s.count('\\end{longtable}')==1
assert 'Pala et al., 2024' not in s and '97451' not in s
(p/'editorial_updated.tex').write_text(s,encoding='utf-8');print(json.dumps({'figures_total':9,'main_figures':8,'supplementary_artwork':1,'evidence_map':'longtable with ragged-right columns','references_added':['Berman 2000','Michaud-Agrawal 2011','Swanson 2024'],'remaining_reference_audit':'Other original references not exhaustively verified; not submission-ready.'}))
