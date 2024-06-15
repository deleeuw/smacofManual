---
title: |
    | Smacof at 50: A Manual
    | Part 8: Homogeneity Analysis with Smacof
author: 
- Jan de Leeuw - University of California Los Angeles
date: 'Started April 13 2024, Version of June 14, 2024'
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





**Note:** This is a working manuscript which will be expanded/updated
frequently. All suggestions for improvement are welcome. All Rmd, tex,
html, pdf, R, and C files are in the public domain. Attribution will be
appreciated, but is not required. The files can be found at
<https://github.com/deleeuw> in the repositories smacofCode, smacofManual,
and smacofExamples.


# Simultaneous Non-Metric Unfolding

* There are $m$ *variables*. 
* Variable $j$ has $k_j>1$ *categories*.
* There are $n$ *objects*.
* Each object defines a partial order over the categories of all variables.

In this chapter we minimize the stress loss function
\begin{equation}
\sigma(X,Y_1,\cdots,Y_m):=\sum_{j=1}^m\sum_{i=1}^n\min_{\hat d_i^j\in\Delta_i^j}\sum_{l=1}^{k_j}w_{il}^j(\hat d_{il}^j-d(x_i,y_l^j))^2
(\#eq:snmu)
\end{equation}
Note that for each object and variable there are different sets of transformations $\Delta_j$
and for each variable there different matrices of column scores $Y_j$, but there is only a single matrix of row scores $X$. Also note that index $j$, for variables,
is sometimes used as a subscript and sometimes as a superscript, depending on what
looks best.

If the $\Delta_i^j$ contain zero, then the unconstrained minimum of \@ref(eq:snmu)
is clearly zero. Collapsing all $x_i$ and all $y_l^j$ into a single point makes all distances zero, and thus makes stress zero. Some sort of normalization of either $X$ and/or the $Y_j$ is needed to prevent this trivial solution.

In fact, it is easy to see that a minimum of zero is also possible in the situation where the $\Delta_i^j$ contain the
set of all constant vectors (or all non-negative constant vectors). Collapse all
$x_i$ into a single point, and place all $y_l^j$ on a sphere around this point.
Or collapse all $y_l^j$ and put the $x_i$ on a sphere. This makes all $d(x_i,y_l^j)$
equal and thus makes stress zero.

## Homogeneity Analysis

The Gifi System (@gifi_B_90, @michailidis_deleeuw_A_98, @deleeuw_mair_A_09a) presents non-linear or non-metric versions of the classical linear multivariate analysis techniques (regression, analysis of variance, canonical analysis, discriminant analysis, principal component analysis) as special cases of Homogeneity Analysis, also known as Multiple Correspondence Analysis.

In this section we present homogeneity analysis as a special case of 
minimizing the loss function \@ref(eq:snmu). 

The data are a number of indicator matrices $G_1,\cdots,G_m$. Indicator matrices are binary matrices, with rows that add up to one. They are
used to code categorical variables. Rows corresponds with objects
(or objects), columns with the categories (or levels) of a variable.
An element $g_{il}$ is one in row if object $i$ is in category $l$,
and all other elements in row $i$ are zero.

Homogeneity analysis makes a joint maps in $p$ dimensions of objects
and categories (both represented as points) in such a way that category points are close to the points for the objects in the category. And, vice versa, objects are close to the category points that they score in.

If there is only one variable then it is trivial to make such a 
homogeneous map. We just make sure the object points coincide with
their category points. But there are $j>1$ indicator matrices, corresponding with $m$ categorical variables, and the solution is a compromise trying to achieve homogeneity as well as possible for all variables simultaneously.

In loss function \@ref(eq:snmu) applied to homogeneity analysis
the sets $\Delta_i^j$ are defined in such a way that $\hat d_{il}^j$ is zero if $i$ is in category $l$ of 
variable $j$. There are no constraints on the other $\hat d$'s in row $i$
of variable $j$. Thus for zero loss we want an object to coincide with all $m$ categories it is in. Under this definition of the $\Delta_i^j$ we have
\begin{equation}
\min_{\hat d_i^j\in\Delta_i^j}\sum_{l=1}^{k_j}w_{il}^j(\hat d_{il}^j-d(x_i,y_l^j))^2=w_{il(i,j)}^jd(x_i,y_{l(i,j)}^j)^2
(\#eq:homsnmu)
\end{equation}
where the $l(i,j)$ on the right is the index of the category of variable $j$ that object $i$ is in. 

Using indicators we can write loss function \@ref(eq:homsnmu) as
\begin{equation}
\sigma(X,Y_1,\cdots,Y_m)=
\sum_{j=1}^m\text{tr}\ (X-G_jY_j)'W_j(X-G_jY_j),
(\#eq:matsnmu)
\end{equation}
The $W_j$ are diagonal matrices with
\begin{equation}
w_i^j:=w_{il(i,j)}=\sum_{l=1}^{k_j}g^j_{il}w^j_{il}.
(\#eq:whomdef)
\end{equation}

In homogeneity analysis we minimize \@ref(eq:matsnmu) using the
explicit normalization $X'W_\star X=I$, where $W_\star$ is the
sum of the $W_j$. The solution is given by the singular value
equations
\begin{align}
X\Lambda&=W_\star^{-1}\sum_{j=1}^m W_jG_jY_j,(\#eq:homsvd1)\\
Y_j&=(G_j'W_jG_j)^{-1}G_j'W_jX,(\#eq:homsvd2)
\end{align}
where $\Lambda$ is a symmetric matrix of Lagrange multipliers. Remember that in addition we require $X'W_\star X=I$.

In homals (@gifi_B_80, @deleeuw_mair_A_09a) alternating least squares is used
to solve the equations \@ref(eq:homsvd1) and \@ref(eq:homsvd2). We start with
some initial $X$, then compute the corresponding $Y_j$ using \@ref(eq:homsvd2),
then for these new $Y_j$ we compute a new corresponding $X$ from \@ref(eq:homsvd1),
and so on. Computations are simple, because only diagonal matrices need to 
be inverted. Alternating least squares becomes reciprocal
averaging.

Alternative methods of computation (and interpretation) are suggested if we
substitute \@ref(eq:homsvd2) in \@ref(eq:homsvd1) to eliminate the $Y_j$
and obtain an equation in $X$ only. This gives
$$
W_\star X\Lambda=\sum_{j=1}^m W_jG_j(G_j'W_jG_j)^{-1}G_j'W_jX,
$$
a generalized eigenvalue equation for $X$. If we substitute \@ref(eq:homsvd1)
in \@ref(eq:homsvd2) we obtain generalized eigenvalue equations for $Y$.
$$
(G_j'W_jG_j)Y_j\Lambda=\sum_{h=1}^m G_j'W_jW_\star^{-1}W_hG_hY_h.
$$
If $k_\star$ is not too large then computing the $p + 1$ largest eigenvalues
from ..., using an efficient 


# Loss Function


Now consider the closely related problem in which we do not require,
as in homogeneity analysis, that 
$$
\sum_{l=1}^{k_j}g^j_{il}\hat d^j_{il}=0
$$
for all $i$ and $j$, but we impose the weaker condition that for all $i$ and $j$
$$
\sum_{l=1}^{k_j}g^j_{il}\hat d^j_{il}\leq \hat d^j_{i\nu}
$$
for all $\nu=1,\cdots,k_j$.

In homogeneity analysis the geometric interpretation of loss is that we
want objects to coincide with all categories they score in. The geometric interpretation of loss function ... is that we want 
objects to be closer to the categories they score in than to the categories
they do not score in. 

This can be formalized using the notion of Voronoi regions. The Voronoi region of
category $l$ of variable $j$ is the polyhedral convex set of all points of $\mathbb{R}^p$ closer to category $l$ than to any other category of variable $j$. The plot of the the $k_j$ categories of variable $j$ defines $k_j$ Voronoi regions.  Loss function
... vanishes if for each variable the $x_i$ are in the Voronoi regions of the categories they score in. This condition implies, by the way, that the interiors of the convex hulls of the $x_i$ in a given category are disjoint, and the point clouds can consequently be weakly separated by hyperplanes.
Since the category points themselves are in their Voronoi region
the convex hulls of the stars are also disjoint.

## The Guttman Transform

General theory





## The Unconstrained Case

###

Minimizing ... over the rows $\delta_i^j$ is a monotone regression for a simple tree order. This is easily handled by using Kruskal's primary approach
to ties (@kruskal_64a, @kruskal_64b, @deleeuw_A_77).

In stage 2 we do one or more metric smacof iterations for given $\Delta_j$
to decrease the loss. These smacof iterations, or Guttman transforms, more or less ignore the fact that we are dealing with a rectangular matrix and use the weights to transform the problem into a symmetric one (as in @heiser_deleeuw_A_79).

Thus for stage two purposes the loss function is
$$
\sigma(Z_1,\cdots,Z_m)=\sum_{j=1}^m\sum_{i=1}^{N_j}\sum_{j=1}^{N_j}w_{ij}^r(\delta_{ij}^r-d_{ij}(Z_j))^2,
$$
with $N_j:=n+k_j$ and

In order to not have to deal with general inverses and multiplications of gigantic symmetric matrices with largely empty diagonal blocks we solve the system
$$
\begin{bmatrix}
W_j&-W\\
-W'&W_c
\end{bmatrix}
\begin{bmatrix}
X\\Y
\end{bmatrix}
=
\begin{bmatrix}
P\\Q
\end{bmatrix}
$$
for $X$ and $Y$, with the known right hand side
$$
\begin{bmatrix}P\\Q\end{bmatrix}=
\begin{bmatrix}
B_j&-B\\
-B'&B_c
\end{bmatrix}
\begin{bmatrix}
X^{(k)}\\Y^{(k)}
\end{bmatrix}
$$

$W_jX-WY=P$ or $X=W_j^{-1}(P+WY)$. Substitute in $W_cY-W'X=Q$
to get $W_cY-W'W_j^{-1}(P+WY)=Q$ or $(W_c-W'W_j^{-1}W)Y=Q+W'W_j^{-1}P$.

Note that $W_c-W'W_j^{-1}W$ is doubly-centered and $Q+W'W_j^{-1}P$ is column-centered.

Computation of the Guttman transform requires Moore-Penrose inverses of the
matrices $V_j$, which are of order $n+k_j$ and extremely sparse. They can also be uncomfortably large because $n$ can be large. It is much more memory-friendly to solve the partioned system
$$
\begin{bmatrix}
R_W&-W\\
-W'&C_W
\end{bmatrix}
\begin{bmatrix}X\\Y\end{bmatrix}=\begin{bmatrix}
R_B&-B\\
-B'&C_B
\end{bmatrix}
$$

More explicitly we must minimize

\begin{multline}
\sum_{j=1}^m\text{tr}\ (X-\tilde X_j)'R_j(X-\tilde X_j)-2\sum_{j=1}^m\text{tr}\ (X-\tilde X_j)'W_j(Y_j-\tilde Y_j)+\\
\sum_{j=1}^m\text{tr}\ (Y_j-\tilde Y_j)'C_j(Y_j-\tilde Y_j)
\end{multline}
where $\tilde X_j$ and $\tilde Y_j$ are the two components of the Guttman
transform $\tilde Z_j$ of the current $Z_j$.

The stationary equations are
\begin{align}
Y_j&=\tilde Y_j-C_j^{-1}W_j'(X-\tilde X_j),(\#eq:stateq1)\\
X&=\{V_{11}^\star\}^{-1}\sum_{j=1}^m\left\{V_{11}^j\tilde X_j-V_{12}^j(Y_j-\tilde Y_j)\right\}.(\#eq:stateq2)
\end{align}
We could solve these equations iteratively using alternating least squares. This 
means using \@ref(eq:stateq1) to compute a new $Y$ for given $X$ and \@ref(eq:stateq2) to compute a new $X$ for given $Y$. This introduces an
infinite inner iteration loop within the main outer iteration loop. We have tried this, and it does not seem to be the way to go. The argument is the same as in 
our section ... on homogeneity analysis.

Alternatively, we can substitute \@ref(eq:stateq2) into \@ref(eq:stateq1). This
gives linear equations in the $Y_j$, which we can solve. Then use \@ref(eq:stateq2)
to compute the corresponding optimal $X$. No inner iterations are necessary. 
In the same way we can substitute \@ref(eq:stateq1) into \@ref(eq:stateq2), solve for
$X$ and compute the corresponding $Y_j$. Since in most applications the number of
objects $n$ is much larger than the total number of categories $\sum k_j$ the first
substitution of \@ref(eq:stateq2) into \@ref(eq:stateq1) seems the most 
promising. Again, the argument is the same as in 
our section ... on homogeneity analysis.

$$
C_j(Y_j-\tilde Y_j)-W_j'R_\star^{-1}\sum_{h=1}^mW_h(Y_h-\tilde Y_h)=W_j'(\tilde X_j-R_\star^{-1}\sum_{h=1}^mR_h\tilde X_h)
$$

$$
(Y_j-\tilde Y_j)-\{V_{22}^j\}^{-1}V_{21}^j\{V_{11}^\star\}^{-1}\sum_{h=1}^mV_{12}^h(Y_h-\tilde Y_h)=\{V_{22}^j\}^{-1}V_{21}^j\tilde X_j-\{V_{22}^j\}^{-1}V_{21}^j\{V_{11}^\star\}^{-1}\sum_{h=1}^mV_{11}^h\tilde X_h
$$
$$
C_j(Y_j-\tilde Y_j)-W_j'R_\star^{-1}\sum_{h=1}^mW_h(Y_h-\tilde Y_h)=-W_j\tilde X_j+W_j'R_\star^{-1}\sum_{h=1}^mR_h\tilde X_h
$$
## Missing Data


## Normalization of X

Constraint 
$$
X'V_{11}^\star X = I
$$

## The Unweighted Case

In the unweighted case we have $V_{11}^r=k_jI$, $V_{22}^r=nI$, and
$V_{12}^r=-E_{n\times k_j}$. Thus
$V_{1\mid2}^r=k_jJ_n$ and $V_{1\mid 2}^\star=k_\star J_n$
with $J_n=I_n-n^{-1}E_{n\times n}$, the centering matrix of order $n.

The Guttman transform also simplifies in the unweighted case. The formulas were already given in @heiser_deleeuw_A_79.

## Centroid Constraints on Y

$$
Z_j=
\kbordermatrix{
\mbox{\ }&p\\
n&I\\
k_j&D_j^{-1}G_j'}X=H_jX
$$
$$
\sum_{j=1}^m Z_j'V_jZ_j=X'\{\sum_{j=1}^mH_j'V_jH_j\}X.
$$

$$
X'H_j'V_jH_jX=X'\{V_{11}^j+V_{12}^jD_j^{-1}G_j+G_jD_j^{-1}V_{22}^jD_j^{-1}G_j'\}X
$$

## Rank Constraints on Y

As in homals it is possble to impose rank-one restrictions on some or all of the $Y_j$. This means we require
$$
y_l^j=\alpha_{l}^jy_l^{\ }{}
$$
Geometrically having all $y_l^j$ on a line through the origin implies that
all voronoi boundaries are hyperplanes perpendicular to that line, and consequently all voronoi regions are bounded by two parallel hyperplanes.

# Convergence and Degeneracy

# Utilities

## Object Plot Function

## Category Plots Function

## Joint Plot Function

# Examples

Presenting examples runs into some difficulties. 

## Cetacea

## Senate

## GALO




# References
