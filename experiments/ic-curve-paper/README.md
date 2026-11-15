# Paper build record

`verification.json` is written by `scripts/build_ic_curve_paper.py`. That script regenerates the sweep tables
(`paper/curve-sweep-data.tex`) and the sweep figure, then builds `paper/curve-lfd-bounds.pdf` with pdfLaTeX, or
with Tectonic when pdfLaTeX is not installed. The record lists the TeX engine, the LaTeX warnings, the overfull
boxes, and the SHA-256 of the inputs and outputs. It certifies that the PDF builds without unresolved citations
or references. It is not a proof or novelty certificate.
