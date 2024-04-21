---
title: |
    | Smacof at 50: A Manual
    | Part 6: Smacof with Categorical Data
author: 
- Jan de Leeuw - University of California Los Angeles
date: 'Started April 13 2024, Version of April 21, 2024'
output:
  bookdown::html_document2:
    keep_md: yes
    css: preamble.css
    toc: true
    toc_depth: 3
    number_sections: yes
  bookdown::pdf_document2:
    latex_engine: lualatex 
    includes:
      in_header: preamble.tex
    keep_tex: yes
    toc: true
    toc_depth: 3
    number_sections: yes
graphics: yes
mainfont: Times New Roman
fontsize: 12pt
bibliography: ["mypubs.bib","total.bib"]
abstract: smacofVO
---





**Note:** This is a working paper which will be expanded/updated frequently. All suggestions for improvement are welcome. 

# Introduction

## Metric Unfolding

In metric unfolding we minimize a loss function of the form
\begin{equation}
\sigma(X,Y)=\sum_{i=1}^n\sum_{j=1}^mw_{ij}(\delta_{ij}-d(x_i,y_j))^2
\end{equation}
over the $n\times p$ matrix $X$ and the $m\times p$ matrices $Y$. Both
the weights $W=\{w_{ij}\}$ and the dissimilarities $\Delta=\{\delta_{ij}\}$
are known non-negative matrices.

The multidimensional unfolding model for preference judgments is often attributed to @coombs, @kruskal_carroll, @roskam

## Non-metric Unfolding


## Simultaneous Non-Metric Unfolding

If we have data where the same individuals give preference
judgments over multiple domains, or over the same domain on multiple occasions, or over the same domain with different experimental conditions,
then we can use the loss function
$$
\sigma(X,Y_1,\cdots,Y_s)=\sum_{r=1}^s\sum_{i=1}^n\min_{\delta_i^r\in\Delta_i^r}\sum_{j=1}^{m_r}w_{ij}^r(\delta_{ij}^r-d(x_i,y_j^r))^2
$$
Note that for each occasion there are different transformations $\Delta_r$
and different matrices of column scores $Y_r$, but there is only a single matrix of row scores $X$.

$$
y_j^r=\alpha_{j}^ry_r
$$

## Homogeneity Analysis

The Gifi System (@gifi_B_90, @michailidis_deleeuw_A_98, @deleeuw_mair_A_09a) presents non-linear or non-metric versions of the classical linear multivariate analysis techniques (regression, analysis of variance, canonical analysis, discriminant analysis, principal component analysis) as special cases of Homogeneity Analysis, also known as Multiple Correspondence Analysis.

We give a somewhat non-standard introduction to homogeneity analysis here, to highlight the
similarities with unfolding and the techniques we will present later on in this paper.

The data are a number of indicator matrices $G_1,\cdots,G_s$. Indicator matrices are binary matrices, with rows that add up to one. They are
used to code categorical variables. Rows corresponds with objects
(or individuals), column with the categories (or levels) of a variable.
An element $g_{ij}$ is one in row if object $i$ is in category $j$,
and all other elements in row $i$ are zero.

Homogeneity analysis makes a joint maps in $p$ dimensions of individuals
and categories (both represented as points) in such a way that category points are close to the points for the individuals in the category. And, vice versa, individuals are close to the category points that they score in.
If there is only one variable then it is trivial to make such a 
homogeneous map. We just make sure the individual points coincide with
their category points. But there are $s>1$ indicator matrices, corresponding with $s$ categorical variables, and the solution is a compromise trying to achieve homogeneity as well as possible for all variables simultaneously.


$$
\sigma(X,Y_1,\cdots,Y_s)=\sum_{r=1}^s\sum_{i=1}^nw_i^rd^2(x_i,y_i^r)=
\sum_{r=1}^s\text{tr}\ (X-G_rY_r)'W_r(X-G_rY_r)
$$
$$
w_i^r:=\sum_{j=1}^{m^r}g_{ij}^rw_{ij}^r
$$
$$
y_{i}^r:=\sum_{j=1}^{m^r}g_{ij}^ry_r^r
$$
Constraint is that $\delta_{ij}^r$ is zero if $i$ is in category $j$ of 
variable $r$. There are no constraints on the other delta's in row $i$.
of variable $r$. Thus we want an object to coincide with all $s$ categories 
it is in. 

star plot

$$Y_r=(G_r'W_rG_r)^{-1}G_r'W_rX$$

$$
\min_Y\sigma(X,Y_1,\cdots,Y_s)=\text{tr}\ X'\left\{\sum_{r=1}^s\left\{W_r-W_rG_r(G_r'W_rG_r)^{-1}G_r'W_r\right\}\right\}X
$$

$X'W_\star X=I$


# Voronoi Loss Function

## The Unconstrained Case

Now consider the closely related problem in which we do not require,
as in homogeneity analysis, that 
$$
\delta_{i\star}^r:=\sum_{j=1}^{m_r}g_{ij}^r\delta_{ij}^r=0
$$
but we merely require that $\delta_{i\star}^r$ is less than or equal to all
$\delta_{ij}^r$ in row $i$. Formally
$$
\delta_{i\star}^r\leq\delta_{ij}^r\qquad\forall j=1,\cdots,m_r
$$
or
$$
g^r_{il}=1\Rightarrow\delta^r_{il}\leq\delta^r_{iv}\qquad\forall v\not=\ell
$$

$$
\sigma(X,Y_1,\cdots,Y_m) = 
\sum_{r=1}^s\left\{\sum_{i=1}^n\sum_{\ell=1}^{k_r}w_{ij}^r(\delta_{i\ell}^r-d(x_i,y_\ell^r))^2\right\}
$$
### ALS Stage 1

Minimizing ... over the rows $\delta_i^r$ is a monotone regression for a simple tree order. This is easily handled by using Kruskal's primary approach
to ties (@kruskal_64, @deleeuw_A_).

###ALS Stage 2

In stage 2 we do one or more metric smacof iterations for given $\Delta_r$
to decrease the loss. We more or less ignore the fact that we are dealing with a rectangular matrix and use the weights to transform the problem into a symmetric one (as in @heiser_deleeuw_A_79).

Thus for stage two purposes the loss function is
$$
\sigma(Z_1,\cdots,Z_m)=\sum_{r=1}^s\sum_{i=1}^{N_r}\sum_{j=1}^{N_r}w_{ij}^r(\delta_{ij}^r-d_{ij}(Z_r))^2,
$$
with $N_r:=n+k_r$ and 
$$
Z_r:=\kbordermatrix{
\mbox{\ }&p\\
n&X\\
k_r&Y_r}.
$$
The weights in $W_r$ are now zero for the two diagonal blocks.
$$
d^2(x_i,y_\ell^r)=\text{tr}\ Z_r'A^r_{i\ell}Z_r^{\ }
$$
$$
\sum_{i=1}^n\sum_{\ell=1}^{k_r}w_{il}^rd^2(x_i,y_\ell^r)=\text{tr}\ Z_r'V_rZ_r
$$

$$
V_r=\sum_{i=1}^n\sum_{\ell=1}^{k_r}w_{il}^rA^r_{i\ell}
$$
$$
V_r=\kbordermatrix{&n&k_r\\
n&\ddots&\Box\\
k_r&\Box&\ddots}
$$


$$
\sum_{i=1}^n\sum_{\ell=1}^{k_r}w_{il}^r\delta_{i\ell}^rd(x_i,y_\ell^r)=\text{tr}\ Z_r'B_r(Z_r)Z_r\geq\text{tr}\ Z_r'B_r(\tilde Z_r)\tilde Z_r=\text{tr}\ Z_r'V_rG(\tilde Z_r)
$$
Using standard smacof majorization theory in an iteration we must minimize

$$
\text{tr}\ (X-\tilde X_r)'V_{11}^r(X-\tilde X_r)+2\text{tr}\ (X-\tilde X_r)'V_{12}^r(Y_r-\tilde Y_r)+\text{tr}\ (Y_r-\tilde Y_r)'V_{22}^r(Y_r-\tilde Y_r)
$$
where $\tilde X_r$ and $\tilde Y_r$ are the two components of the Guttman
transform $\tilde Z_r$of the current $Z_r$.

First minimize over $Y_r$, without constraints.
$$
Y_r-\tilde Y_r=-\{V_{22}^r\}^{-1}V_{21}^r(X-\tilde X_r)
$$
and thus
$$
\text{tr}\ (X-\tilde X_r)'V_{1|2}^r(X-\tilde X_r)
$$
where
$$
V_{1|2}^r:=V_{11}^r-V_{12}^r\{V_{22}^r\}^{-1}V_{21}^r
$$
i.e. Schur complement of .. in .. Note that $V_{1|2}^r$ is doubly-centered.


Constraint either (weak)
$$
\text{tr}\ X'V_{1|2}^\star X = 1
$$
or (strong)
$$
X'V_{1|2}^\star X = I
$$
Thus (weak)
$$
X=\lambda\{V_{1|2}^\star\}^{-1}\sum_{j=1}^mV_{1|2}^r\tilde X_r
$$
with $\lambda$ chosen such that ... is satisfied. 

 Or (strong)
$$
X=\{V_{1|2}^\star\}^{-\frac12}KL'
$$
(Procrustus)

In the  unweighted case we have $V_{11}^r=k_rI$, $V_{22}^r=nI$, and
$V_{12}^r=-e_ne_{k_r}'$. Thus
$V_{1\mid2}^r=k_rJ_n$ and $V_{1\mid 2}^\star=k_\star J_n$
with $J_n=I-n^{-1}e_ne_n'$ the centering matrix of order $n.

The Guttman transform also simplifies in the unweighted case. The formulas were already given in @heiser_deleeuw_A_79.

## The Constrained Case

$$
Z_r=
\kbordermatrix{
\mbox{\ }&p\\
n&I\\
k_r&D_r^{-1}G_r'}X=H_rX
$$
$$
\sum_{j=1}^m Z_r'V_rZ_r=X'\{\sum_{j=1}^mH_r'V_rH_r\}X.
$$

## The Unweighted Case

If unweighted then $V_{11}^r=k_rI$ and $V_{22}^r=nI$. $V_{12}^r$ is $n\times k_r$ with all elements equal to $-1$. The rank of $V_r$ is $n+k_r-1$, and the only vectors in  the null-space are proportional to the vector with all elements equal to one (). 

$$V_r^+=(V_r+ee)^{-1}-\frac{1}{n+k_r}ee'$$
But 
$$V_r+ee'=(V_{11}^r+ee')\oplus(V_{22}^r+ee')=(k_rI+ee')\oplus(nI+ee')$$
$$(V_r+ee')^{-1}=(k_rI+ee')^{-1}\oplus(nI+ee')^{-1}$$
Lemma: If $A$ is a symmetric matrix of order $n$ of the form
$\alpha I+ee'$ then $(\alpha I+ee')^{-1}=$

Proof: 

@heiser_deleeuw_A_79

Sherman-Morrison

$$
(n I+ee')^{-1}=\frac{1}{n}\left\{I-\frac{1}{k_r+n}ee'\right\}$$
$$
(k_r I+ee')^{-1}=\frac{1}{k_r}\left\{I-\frac{1}{k_r+n}ee'\right\}$$

Schur complement $$V_{11}-V_{12}V_{22}^{-1}V_{21}=kI-e_ne_k'\frac{1}{n}e_ke_n'=k\{I-\frac{1}{n}e_ne_n'\}=kJ_n$$


# References
