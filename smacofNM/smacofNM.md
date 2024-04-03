---
title: |
    | Smacof at 50: A Manual
    | Part 2: Non-metric Smacof
author: 
- Jan de Leeuw - University of California Los Angeles
date: 'Started March 30 2024, Version of April 02, 2024'
output:
  bookdown::html_document2:
    keep_md: yes
    css: preamble.css
    toc: true
    toc_depth: 4
    number_sections: yes
  bookdown::pdf_document2:
    latex_engine: lualatex
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: true
    toc_depth: 4
    number_sections: yes
graphics: yes
mainfont: Times New Roman
fontsize: 12pt
bibliography: ["mypubs.bib","total.bib"]
abstract: TBD
editor_options: 
  markdown: 
    wrap: 72
---



**Note:** This is a working manuscript which will be expanded/updated
frequently. All suggestions for improvement are welcome. All Rmd, tex,
html, pdf, R, and C files are in the public domain. Attribution will be
appreciated, but is not required. The code files can be found at
<https://github.com/deleeuw/smacofCode>, the manual files at
<https://github.com/deleeuw/smacofManual>, and the example files
at <https://github.com/deleeuw/smacofExamples>.

\sectionbreak

# Introduction

pick and rank

# Loss Function
$$
\sigma(X,\delta_1,\cdots,\delta_s)=\frac{\sum_{r=1}^s\omega_r\sum_{i,j} w_{ijr}(\delta_{ijr}-d_{ij}(X))^2}{\sum_{r=1}^s\omega_r\sum_{i,j} w_{ijr}d_{ij}^2(X))}
$$
which must be minimized over $X$ and over $\delta_r\in\mathcal{K}_r$, with $\mathcal{K}_r$ pointed polyhedral convex cones, defined by a partial order 
$\leq_r$.

Minimize of $X$ for given $\delta_{ijr}$. 
$$
\sigma_R(X,\delta_1,\cdots,\delta_s)=\sum_{r=1}^s\omega_r\sum_{i,j} w_{ijr}\delta_{ijr}^2-2\sum_{r=1}^s\omega_r\sum_{i,j} w_{ijr}\delta_{ijr}d_{ij}(X)+\sum_{r=1}^s\omega_r\sum_{i,j} w_{ijr}d_{ij}^2(X))
$$
Nor use
$$
d_{ij}(X)\geq\frac{1}{d_{ij}(Y)}\text{tr}\ X'A_{ij}Y 
$$
so that
$$
\sum_{r=1}^s\omega_r\sum_{i,j} w_{ijr}\delta_{ijr}d_{ij}(X)\geq
\text{tr}\ X'B(Y)Y
$$
$$
B(Y):=\sum_{r=1}^s\omega_r\sum_{i,j} w_{ijr}\frac{\delta_{ijr}}{d_{ij}(Y)}A_{ij}
$$
Also
$$
V:=\sum_{r=1}^s\omega_r\sum_{i,j} w_{ijr}A_{ij}
$$
So that
$$
\sigma_R(X)\leq K-2\text{tr}\ X'B(Y)Y+\text{tr}\ X'VX
$$
and the smacof update over $X$ with $\text{tr}\ X'VX=1$ is the same 
as in smacofRR.

# Paired Comparisons

Positive Orthant / Absolute Value / Pairwise

@deleeuw_R_70a
@deleeuw_E_18d
@hartmann_79
@guttman_69
@johnson_73


Suppose datum $r$ says that that $(i,j)<(k,l)$. Then $w_{ijr}$ and $w_{klr}$
are non-zero and all other elements of $W_r$ are zero.
Thus
$$
w_{ijr}(\delta_{ijr}-d_{ij})^2+w_{klr}(\delta_{klr}-d_{kl})^2
$$
Must be minimized over $\delta_{ijr}\leq\delta_{klr}$. If $d_{ij}\leq d_{kl}$
then $\hat d_{ijr}=d_{ij}$ and $\hat d_{klr}=d_{kl}$, and otherwise
$$
\hat d_{ijr}=\hat d_{klr}=\frac{w_{ijr}d_{ij}+w_{klr}d_{kl}}{w_{ijr}+w_{klr}}
$$
 Thus

$$w_{ijr}(\delta_{ijr}-d_{ij})^2+w_{klr}(\delta_{klr}-d_{kl})^2$$
is zero if the order of $d_{ij}$ and $d_{kl}$ is the same as the order in the data
and 

$$
w_{ijr}(\frac{w_{ijr}d_{ij}+w_{klr}d_{kl}}{w_{ijr}+w_{klr}}-d_{ij})^2
+w_{klr}(\frac{w_{ijr}d_{ij}+w_{klr}d_{kl}}{w_{ijr}+w_{klr}}-d_{kl})^2
$$
$$
\frac{w_{ijr}w_{klr}^2}{(w_{ijr}+w_{klr})^2}(d_{ij}-d_{kl})^2+
\frac{w_{ijr}^2w_{klr}}{(w_{ijr}+w_{klr})^2}(d_{ij}-d_{kl})^2
$$
$$
\frac{w_{ijr}w_{klr}}{w_{ijr}+w_{klr}}(d_{ij}-d_{kl})^2
$$
$B$ matrix

So far we have only considered the forced-choice situation in which
the subject has to choose one of the pairs. If we allow for the alternative
that $(i,j)$ and $(k,l)$ are equally similar then we have two approaches. In the primary approach we 

$$
\begin{bmatrix}
0&0&1\\
1&0&1\\
1&0&0
\end{bmatrix}
\begin{bmatrix}
0&0&1\\
0&0&0\\
1&0&0
\end{bmatrix}
$$

# Triads

# pick 1/3 Method of triads
Pick the two most similar 
$(i,j)<(i,k)$ and $(i,j)<(j,k)$

# order 3/3 Complete method of triads
$(i,j)<(i,k)<(j,k)$

# Richardson -- 
each triple presented three times, with a different hub each time
which one of the two is maximally similar to the hub stimulus. Thus the data
is the single inequality $(i,k) < (j,k)$

# Conditional rank Orders -- Klingberg

# References
