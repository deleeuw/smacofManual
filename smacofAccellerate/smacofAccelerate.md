---
title: "Accelerated Least Squares Multidimensional Scaling"
author: 
- Jan de Leeuw - University of California Los Angeles
date: 'Started July 23 2024, Version of August 20, 2024'
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
abstract: We discuss a simple accelerations of MDS smacof iterations, and compare them with recent boosted difference-of-convex algortithms.
---



**Note:** This is a working manuscript which will be expanded/updated
frequently. All suggestions for improvement are welcome. All Rmd, tex,
html, pdf, R, and C files are in the public domain. Attribution will be
appreciated, but is not required. The files can be found at
<https://github.com/deleeuw> in the repositories smacofCode, smacofManual,
and smacofExamples.

# Introduction

In this paper we study minimization of the multidimensional scaling (MDS) loss function
\begin{equation}
\sigma(X):=\frac12\mathop{\sum\sum}_{1\leq i<j\leq n} w_{ij}(\delta_{ij}-d_{ij}(X))^2
(\#eq:sdef)
\end{equation}
over all $n\times p$ *configuration* matrices $X$. Following @kruskal_64a, @kruskal_64b we call $\sigma(X)$ the *stress* of $X$. The symbol $:=$ is used for definitions.

In definition 
\@ref(eq:sdef) matrices $W=\{w_{ij}\}$ and $\Delta=\{\delta_{ij}\}$ are known non-negative, symmetric, and hollow. They contain, respectively, *weights* and *dissimilarities*. The matrix-valued function $D$, with $D(X)=\{d_{ij}(X)\}$, contains *Euclidean distances* between the rows of $X$. 

Throughout we assume, without loss of generality, that $W$ is irreducible, that $X$ is  column-centered, and that $\Delta$ is normalized by
\begin{equation}
\frac12\mathop{\sum\sum}_{1\leq i<j\leq n} w_{ij}\delta_{ij}^2=1.
(\#eq:delnorm)
\end{equation}

## Notation

It is convenient to have some matrix notation for the MDS problem.
We use the symmetric matrices $A_{ij}$, of order $n$, which have $+1$ for elements $(i,i)$ and $(j,j)$, $-1$ for elements $(i,j)$ and $(j,i)$, and zeroes everywhere else. Using unit vectors $e_i$ and $e_j$ we can write
\begin{equation}
A_{ij}:=(e_i-e_j)(e_i-e_j)'.
(\#eq:adef)
\end{equation}

Following @deleeuw_C_77 we define
\begin{equation}
\rho(X):=\mathop{\sum\sum}_{1\leq i<j\leq n} w_{ij}\delta_{ij}d_{ij}(X)=\text{tr}\ X'B(X)X,
(\#eq:rhodef)
\end{equation}
where
\begin{equation}
B(X):=\mathop{\sum\sum}_{1\leq i<j\leq n}w_{ij}\delta_{ij}r_{ij}(X)A_{ij},
(\#eq:bdef)
\end{equation}
with
\begin{equation}
r_{ij}(X)=\begin{cases}
d_{ij}^{-1}(X)&\text{ if }d_{ij}(X)>0,\\
0&\text{ if }d_{ij}(X)=0.
\end{cases}
(\#eq:rdef)
\end{equation}

Also define
\begin{equation}
\eta^2(X):=\mathop{\sum\sum}_{1\leq i<j\leq n}w_{ij}d_{ij}^2(X)=\text{tr}\ X'VX,
(\#eq:etadef)
\end{equation}
where
\begin{equation}
V:=\mathop{\sum\sum}_{1\leq i<j\leq n}w_{ij}A_{ij}.
(\#eq:vdef)
\end{equation}
Thus
\begin{equation}
\sigma(X)=1-\rho(X)+\frac12\eta^2(X)=1-\text{tr}\ X'B(X)X+\frac12\text{tr}\ X'VX.
(\#eq:sform)
\end{equation}
Both $B(X)$ and $V$ are positive semi-definite and doubly-centered. Because of the irreducibility of $W$ the matrix $V$ has rank $n-1$, with only the constant vectors in its null space. Both $\rho$ and $\eta$ are positively homogeneous convex functions, with $\eta$ being a norm on the space of column-centered configurations. 

Note that $\rho$ is continuous, but it is not differentiable at $X$ if
$d_{ij}(X)=0$ for some $(i,j)$ for which $w_{ij}\delta_{ij}>0$. Because
\begin{equation}
|d_{ij}(X)-d_{ij}(Y)|^2\leq\text{tr}\ (X-Y)'A_{ij}(X-Y)\leq 2p\|X-Y\|^2
(\#eq:lipschitz)
\end{equation}
we see that $\rho$, although not differentiable, is globally Lipschitz.

## The Guttman Transform

The *Guttman transform* of a configuration $X$, so named by @deleeuw_heiser_C_80 to honor the contribution of @guttman_68, is defined as the set-valued map
\begin{equation}
\Phi(X)=V^+\partial\rho(X),
(\#eq:phidef)
\end{equation}
with $V^+$ the Moore-Penrose inverse of $V$ and $\partial\rho(X)$ the
subdifferential of $\rho$ at $X$, i.e. the set of all $Z$ such that
$\rho(Y)\geq\rho(X)+\text{tr}\ Z'(Y-X)$
for all $Y$. Because $\rho$ is homogeneous of degree one we have that $Z\in\partial\rho(X)$
if and only if $\text{tr}\ Z'X=\rho(X)$ and
$\rho(Y)\geq\text{tr}\ Z'Y$ for all $Y$. For each $X$ 
the subdifferential $\partial\rho(X)$, and consequently the Guttman transform, is compact and convex. The map $\partial\rho$
is also positively homogeneous of degree zero, i.e. $\partial\rho(\alpha X)=\partial\rho(X)$ for all $X$ and all $\alpha\geq 0$. And consequently so is the
Guttman transform.

We start with the subdifferential of the distance function between
rows $i$ and $j$ of an $n\times p$ matrix. Straightforward calculation
gives
\begin{equation}
\partial d_{ij}(X)=\begin{cases}
\left\{d_{ij}^{-1}(e_i-e_j)(x_i-x_j)'\right\}&\text{ if }d_{ij}(X)>0,\\
\left\{Z\mid Z=(e_i-e_j)z'\text{ with }z'z\leq1\right\}&\text{ if }d_{ij}(X)=0.
\end{cases}
(\#eq:dsubsef)
\end{equation}
Thus if $d_{ij}(X)>0$, i.e. if $d_{ij}$ is differentiable at $X$,
then $\partial d_{ij}(X)$ is a singleton, containing only the gradient
at $X$.

From subdifferential calculus (@rockafellar_70, theorem 23.8 and 23.9) the subdifferential of $\rho$ is the linear combination
\begin{equation}
\partial\rho(X)=\mathop{\sum\sum}_{1\leq i<j\leq n}w_{ij}\delta_{ij}\partial d_{ij}(X)
(\#eq:subdif)
\end{equation}
Summation here is in the Minkovski sense, i.e. $\partial\rho(X)$ is the compact convex set of all linear combinations $\mathop{\sum\sum}_{1\leq i<j\leq n}w_{ij}\delta_{ij}z_{ij}$,
with $z_{ij}\in\partial d_{ij}(X)$.

It follows that
\begin{equation}
\partial\rho(X)=B(X)X+Z
(\#eq:rhosubdef)
\end{equation}
with
\begin{equation}
Z\in\mathop{\sum\sum}\{w_{ij}\delta_{ij}\partial d_{ij}(X)\mid d_{ij}(X)=0\}.
(\#eq:zsubdef)
\end{equation}
It also follows that
\begin{equation}
\partial\sigma(X)=VX-\partial\rho(X)
(\#eq:sigsubdef)
\end{equation}
Since $\sigma$ is not convex the subdifferential in \@ref(eq:sigsubdef) is
the Clarke subdifferential (@clarke_75). 

Now $X$ is a stationary point of $\sigma$ if $0\in\partial\sigma(X)$, i.e.
if and only if $X\in V^+\partial\rho(X)$. This means that stationary
points are fixed points of the Guttman transform.
A necessary condition for $\sigma$ to have a local minimum at $X$ is that $X$ is a stationary point. The condition is far from sufficient, however, since stationary points can also be saddle points or local maxima. @deleeuw_R_93c shows that stress
only has a single local maximum at the origin $X=0$, but generally there
are many saddle points.

This little excursion into convex analysis is rarely needed in practice. We 
call a configuration $X$ *friendly* if $d_{ij}(X)>0$ for all $(i,j)$ for which $w_{ij}\delta_{ij}>0$. In @deleeuw_A_88b such configurations were called *usable*, but that seems somewhat misleading. All configurations are usable. If $X$ is friendly the rank of $B(X)$ is $n-1$. More importantly,
it was shown by @deleeuw_A_84f that if $\sigma$ has local minimum at $X$
then $X$ is friendly. At friendly configurations (and thus at local minima)
$\sigma$ is differentiable, and the subdifferential 
\@ref(eq:sigsubdef) is a singleton, containing only the gradient.
Stationary points then satisfy $X=V^+B(X)X$.
If $w_{ij}\delta_{ij}=0$ for some $(i,j)$ then there can be local minima
where $\sigma$ is not differentiable. This typically happens in 
multidimensional unfolding (@mair_deleeuw_wurzer_C_15).

By the definition of the subdifferential $Z\in\partial\rho(X)$ implies $\rho(X)\geq\text{tr}\ Z'X$ and $\rho(Y)\geq\text{tr}\ Z'Y$ for all $Y$. If
$d_{ij}(X)>0$ this follows directly from the Cauchy-Schwartz inequality
\begin{equation}
d_{ij}(Y)\geq d_{ij}^{-1}(X)\text{tr}\ X'A_{ij}Y.
(\#eq:csineq)
\end{equation}
Multiplying both sides by $w_{ij}\delta_{ij}$ and summing gives
\begin{equation}
\rho(Y)\geq\text{tr}\ Y'B(X)X
(\#eq:rhoineq)
\end{equation}
for all $Y$, with equality if $Y=X$. Not "if and only if", but just "if". We
also have equality if $Y=\alpha X$ for some $\alpha\geq 0$.

Using the Guttman transform we can use \@ref(eq:rhoineq) to derive the basic smacof equality
\begin{equation}
\sigma(X)=1+\eta^2(X-\Phi(X))-\eta^2(\Phi(X))
(\#eq:smacofequality)
\end{equation}
for all $X$ and the basic smacof inequality
\begin{equation}
\sigma(X)\leq 1+\eta^2(X-\Phi(Y))-\eta^2(\Phi(Y))
(\#eq:smacofinequality)
\end{equation}
for all $X$ and $Y$.

Taken together \@ref(eq:smacofequality) and \@ref(eq:smacofinequality) imply the *sandwich inequality*
\begin{equation}
\sigma(\Phi(Y))\leq 1-\eta^2(\Phi(Y))\leq 1+\eta^2(Y-\Phi(Y))-\eta^2(\Phi(Y))=\sigma(Y).
(\#eq:sandwich)
\end{equation}
If $Y$ is not a fixed point of $\Phi$ then the second inequality in the
chain is strict and thus $\sigma(\Phi(Y))<\sigma(Y)$. As we mentioned, the
first inequality may not be strict.

It also follows
from \@ref(eq:sandwich) that $\eta^2(\Phi(Y))\leq 1$. Thus the Guttman
transform of any configuration is in a convex and compact set, in fact an ellipsoid,
containing the origin.

# One-point Iteration

## Basic Iteration

The basic smacof algorithm generates the iterative sequence
\begin{equation}
X^{(k+1)}=\Phi(X^{(k)}),
(\#eq:basic)
\end{equation}
where it is understood that we stop if $X^{(k)}$ is a fixed point. If
$X^{(k)}$ is not a fixed point it follows from \@ref(eq:sandwich) that $\sigma(X^{(k+1)})<\sigma(X^{(k)})$. Thus, without any additional
assumptions, and using basically only the Cauchy-Schwartz inequality, the algorithm either stops at a fixed point or produces
a monotonically strictly decreasing sequence of loss function values. Since stress
is bounded below by zero the sequence $\sigma(X^{k})$ converges to, say, $\sigma_\infty$.

The original derivation of the smacof algorithm (@deleeuw_C_77, @deleeuw_heiser_C_77)
used the theory of maximization a ratio of norms discussed by @robert_67. Later
derivations (@deleeuw_heiser_C_80, @deleeuw_A_88b) used the fact that \@ref(eq:smacofinequality) defines a majorization scheme for stress. Convergence
then follows from the general *majorization principle* (these days mostly known
as the *MM principle*), introduced in @deleeuw_C_94c. A recent overview of the MM approach is @lange_16.

It was also realized early on that the smacof algorithm was a special case of the 
the difference-of-convex functions algorithm (DCA), introduced by Pham Dinh Tao around 
1980. Pham Dinh also started his work in the context of ratio's of norms, using 
Robert's fundamental ideas. Around 1985 he generalized his approach to minimizing
DC functions of the form $h=f-g$, with both $f$ and $g$ convex. The basic idea
is to use the subgradient inequality $g(x)\geq g(y)+z'(x-y)$, with $z\in\partial g(x)$,
to construct the majorization $h(x):=f(x)-g(y)-z'(x-y)$. Now $h$ is obviously convex in $x$. The DC algorithm then chooses the successor of $y$ as the minimizer of this convex majorizer over $x$. In smacof the role of $f$ is played by $\eta^2$ and the role of $g$ by $\rho$. DCA is applied to MDS in @lethi_tao_01. Extensive recent surveys of the DC/DCA approach are @lethi_tao_18 and @lethi_tao_24.

Thus the smacof algorithm is both MM and DCA, which means that it inherits all
results that have been established for these more general classes of algorithms.
But additional results can be obtained by using the special properties of
the stress loss function and the smacof iterations. In the DCA context, for example, the convex subproblem that must be solved in each step is quadratic, and has the closed form solution provided by the Guttman transform. 

@deleeuw_A_88b derives some additional smacof-specific results. Using up-arrows and down-arrows for monotone convergence he shows that

* $\rho(X^{(k)})\uparrow\rho_\infty$,
* $\eta^2(X^{(k)})\uparrow\eta^2_\infty=\rho_\infty$,
* $\sigma(X^{(k)})\downarrow\sigma_\infty=1-\rho_\infty$,

and, last but not least, the sequence $\{X^{(k)}\}$ is *asymptotically regular*, i.e.
\begin{equation}
\eta^2(X^{(k+1)}-X^{(k)})\rightarrow 0.
(\#eq:etaconv)
\end{equation}

@deleeuw_A_88b argues that these convergence results are sufficient from a 
practical point of view. If we define an $\epsilon$-fixed-point as
any configuration $X$ with $\eta(X-\Phi(X))<\epsilon$ then smacof produces such an
$\epsilon$-fixed-point in a finite number of steps.

Because 

* the subdifferential is a upper semi-continuous (closed) map, 
* all iterates are in the compact set $\eta^2(X)\leq 1$, and 
* $\Phi$ is strictly monotonic (decreases stress at non-fixed points), 

it follows from theorem 3.1 of @meyer_76 that
all accumulation points are fixed points and have the same function value $\sigma_\infty$. Moreover, from theorem 26.1 of @ostrowski_73, either the sequence
converges or the accumulation points form a continuum. 

In order to prove actual convergence, additional conditions are needed.
@meyer_76 proves convergence if the number of fixed points with function value $\sigma_\infty$ is finite, or if the sequence has an accumulation point that is an isolated fixed point. Both these conditions are not met in our case, because
of rotational indeterminacy. If $X_\infty$ is a fixed point, then the
continuum of rotations of $X_\infty$ are all fixed points.

It should be emphasized that smacof can converge to stationary points that
are not local minima (and thus saddle points). Suppose all weights are equal to one,
$\delta_{12}>0$, and $\delta_{1j}=\delta_{2j}$ for all $j>2$. If $d_{12}(X)=0$ then also
$d_{12}(\Phi(X))=0$, and $d_{12}(X)$ will be zero for all iterates, and thus for all subsequential limits, which consequently cannot be local minima. Another example
is starting at $(X\mid 0)$, i.e. $X$ with some columns of zeroes added. All
updates will also have this form, and convergence again cannot be to a local minimum.
This last example can be generalized to any $X$ with rank less than $p$, because
all updates will also have rank less than $p$. As a consequence if $q>p$ 
local minima for $p-$dimensional MDS are saddle points for $q-$dimensional MDS.

In two very recent impressive papers @ram_sabach_24 and @robini_wang_zhu_24 use the powerful Kurdyka-Łojasiewicz (KL) framework (@bolte_daniilidis_lewis_07, @bolte_sabach_teboulle_14 ) to prove actual global convergence of the smacof iterates to a fixed point. We shall use the more classical local convergence analysis, based on the differentiability of the Guttman transform.

## Asymptotic Rate of Convergence

In order to study the asymptotic rate of convergence of smacof, we have to 
compute the derivative of the Guttman transform and its eigenvalues (@ortega_rheinboldt_70, chapter 10). Thus we assume we are in the neighborhood of a configuration, for example a local minmizer,  where the Guttman transform is (infinitely many times) differentiable. 

The derivative of $\Phi$ at $X$, first given in @deleeuw_A_88b, is
the linear transformation $\mathcal{D}\Phi_X$, mapping the space of column-centered
$n\times p$ matrices into itself. Its value at $H$ is equal to
\begin{equation}
\mathcal{D}\Phi_X(H)=V^+\sum w_{ij}\frac{\delta_{ij}}{d_{ij}(X)}\left\{A_{ij}H-\frac{\text{tr}\ X'A_{ij}H}{ \text{tr}\ X'A_{ij}X}A_{ij}X\right\}.
(\#eq:jacobian)
\end{equation}

It follows that $\mathcal{D}\Phi_X(X)=0$ for all $X$ and the derivative has at least one
zero eigenvalue. If we think of equation \@ref(eq:jacobian) as a linear transformation on the space of all $n\times p$ matrices, then there are an additional $p$ zero eigenvalues
corresponding with translational invariance. If we define \@ref(eq:jacobian)
on the space of column-centered matrices, then those zero eigenvalues disappear.

If $S$ is anti-symmetric and
$H=XS$ then $\text{tr}\ X'A_{ij}H=0$ and thus $\mathcal{D}\Phi_X(XS)=\Phi(X)S$.
If in addition $X$ is a fixed point then  $\mathcal{D}\Phi_X(XS)=XS$,
which means that at a fixed point $\mathcal{D}\Phi_X$ has $\frac12p(p-1)$ eigenvalues equal to one. These correspond to the rotational indeterminacy of the MDS problem and the
smacof iterations.

Since $\Phi(X)=V^+\mathcal{D}\rho(X)$ the derivative of the Guttman transform
has a simple relationship with the second derivatives of $\rho$. 
The second derivative, again from @deleeuw_A_88b, is the bilinear form defined by
\begin{equation}
\mathcal{D}^2\rho_X(G,H)=\sum w_{ij}\frac{\delta_{ij}}{d_{ij}(X)}\left\{\text{tr}\ G'A_{ij}H-\frac{\text{tr}\ G'A_{ij}X\text{tr}\ H'A_{ij}X}{d_{ij}^2(X)}\right\}=\text{tr}\ G'V\mathcal{D}\Phi_X(H).
(\#eq:hessian)
\end{equation}
Since $\rho$ is convex, all eigenvalues of $\mathcal{D}\Phi_X$ are real and
nonnegative. It also follows from \@ref(eq:hessian) that if $G$ and $H$ are eigenvectors
of the derivative $\mathcal{D}\Phi_X$ with different eigenvalues then $\text{tr}\ G'VH=0$. In addition we see from equation \@ref(eq:hessian) that
$0\lesssim\mathcal{D}^2\rho_X\lesssim B(X)$
in the Loewner sense. Because
$\mathcal{D}^2\sigma_X=V-\mathcal{D}^2\rho_X$ we have
$\mathcal{D}^2\sigma_X\gtrsim 0$ at a local minimum, and consequently
$\mathcal{D}\Phi_X\lesssim I$. Thus all eigenvalues of the derivative $\mathcal{D}\Phi_X$ at a local minimum $X$ are between zero and one.

We apply basic iterations to the two-dimensional MDS analysis of the classical color-circle example from @ekman_54, which has $n=14$ points. 
In our numerical examples we always use weights equal to one. We always start with the classical Torgerson-Gower solution and we stop if
$\sigma(X^{(k)})-\sigma(X^{(k+1)})<1e-15$. 


The fit for the Ekman example is very good and convergence is rapid. 
In iteration 56, the final iteration, stress is 2.1114112739076. The *change* $\eta(X^{(k)}-X^{(k+1)})$ is 2.26054683805646e-16. The estimated asymptotic rate of convergence or *EARC* is
the change divided by the change of the 
previous iteration. In this analysis it is 0.766978439824377.

We compute the Jacobian corresponding to the derivative $\mathcal{D}\Phi_X$ in two ways. First by using formula \@ref(eq:jacobian),
setting $H$ equal to $e_i^{\ }e_s'$ for $1,\cdots,n$ and $s=1,\cdots,p$. Second, internally as a check and precaution, by using the jacobian function from the numDeriv package (@gilbert_varadhan_19).

The eigenvalues of the derivative $\mathcal{D}\Phi_X$ are


```
##  [1]   +1.0000000000   +0.7669964993   +0.7480939418   +0.7185926294
##  [5]   +0.7007452309   +0.6920114813   +0.6859492533   +0.6593334529
##  [9]   +0.6541779410   +0.6477573343   +0.6237683213   +0.6178713316
## [13]   +0.5735285951   +0.5483330653   +0.5260355535   +0.5112510731
## [17]   +0.5064703617   +0.5059294792   +0.4919752630   +0.4827646555
## [21]   +0.4782034995   +0.4757907675   +0.4682965893   +0.4619226491
## [25]   +0.4559704884   +0.0000000000   +0.0000000000   -0.0000000000
```
Note that the largest non-trivial eigenvalue, which is another and usually better 
eastimate of the ARC, is equal to the EARC for the final iteration.

## Orthogonalized Iteration

As @deleeuw_A_88b mentions, we cannot directly apply the basic point-of-attraction theorem 10.1.3 and the equally basic linear convergence theorem 10.1.4 from @ortega_rheinboldt_70, because at a fixed point of smacof there are $\frac12 p(p-1$ eigenvalues equal to one. 

One way around this problem (@deleeuw_E_19h) is to rotate each update to orthogonality,
i.e. to principal components. Thus the update formula becomes $\Xi(X)=\Pi(\Phi(X))$, with $\Pi(X)=XL$, where $L$ are the right singular vectors of $X$.

Now clearly this modified algorithm generates the same sequence of 
stress values as basic smacof. In fact, it also generates the same
sequences of $\rho$ and $\eta$ values. Moreover $\Xi^n(X)=\Pi(\Phi^n(X))$,
which means that we can find any term of the orthogonal sequence
by orthogonalizing the corresponding term in the basic sequence.
Thus, in actual computation, there is no need to orthogonalize, we
may as well compute the basic sequence and orthogonalize after
convergence.

We compute the derivative of $\Xi$. By the chain rule
\begin{equation}
\mathcal{D}\Xi_X(H)=\mathcal{D}\Pi_{\Phi(X)}(\mathcal{D}\Phi_X(H)).
(\#eq:chain)
\end{equation}
In appendix \@ref(rotation) we show that,
if $X'XL=L\Lambda$ with $L'L=LL'=I$ and the eigenvalues in
$\Lambda$ are all different, then
\begin{equation}
\mathcal{D}\Pi_X(H)=HL+XLS,
(\#eq:pideriv)
\end{equation}
where $S$ is the anti-symmetric matrix with elements
\begin{equation}
s_{ij}=-\frac{l_i'(H'X+X'H)l_j}{\lambda_i-\lambda_j}.
(\#eq:smatdef)
\end{equation}
Thus, from equations \@ref(eq:chain) and \@ref(eq:pideriv)
\begin{equation}
\mathcal{D}\Xi_X(H)=\mathcal{D}\Phi_X(H)L+\Phi(X)LS
(\#eq:xideriv)
\end{equation}
with $L$ and $S$ computed from the singular value decomposition of $\Phi(X)$. 

At a fixed point of $\Xi$ we have both $\Phi(X)=X$ and $\Pi(X)=X$, and consequently
$L=I$ and $X'X=\Lambda$. Equation \@ref(eq:xideriv) becomes
\begin{equation}
\mathcal{D}\Xi_X(H)=\mathcal{D}\Phi_X(H)+XS,
(\#eq:xiderivfixed)
\end{equation}
where now
\begin{equation}
s_{ij}=-\frac{(H'X+X'H)_{ij}}{\lambda_i-\lambda_j}.
(\#eq:sdeffixed)
\end{equation}

Eigenvectors $H$ of $\mathcal{D}\Phi_X$ with eigenvalue one are of the form $H=XA$ with $A$ anti-symmetric.  From \@ref(eq:xiderivfixed)
\begin{equation}
\mathcal{D}\Xi_X(XA)=XA+XS,
(\#eq:xiderivasym)
\end{equation}
where
\begin{equation}
s_{ij}=-\frac{(A'\Lambda+\Lambda A)_{ij}}{\lambda_i-\lambda_j}=-a_{ij}.
(\#eq:sdefasym)\end{equation}
Thus $\mathcal{D}\Xi_X(XA)=0$.

Theoretically, however, orthogonalization gives the same 
convergence rate as the basic sequence, but the Jacobian
of $\Xi$ at a local minimum does not have the unit
eigenvalues any more. They are replaced by zeroes, reflecting
the fact that we are iterating on the nonlinear manifold
or orthogonal column-centered matrices. It is now sufficient
for linear convergence to assume that the largest
eigenvalue of the Jacobian at the solution is strictly
less than one, or alternatively assume that one of the accumulation
points is an isolated local minimum.



In iteration 55, the final iteration, stress is 2.1114112739076. The change is 3.864164805803e-16 and the EARC is 0.766992047059491. The eigenvalues of the
Jocobian are

```
##  [1]   +0.7669964950   +0.7480939419   +0.7185926296   +0.7007452320
##  [5]   +0.6920114814   +0.6859492533   +0.6593334535   +0.6541779411
##  [9]   +0.6477573344   +0.6237683215   +0.6178713316   +0.5735285953
## [13]   +0.5483330652   +0.5260355534   +0.5112510730   +0.5064703617
## [17]   +0.5059294792   +0.4919752629   +0.4827646555   +0.4782035015
## [21]   +0.4757907660   +0.4682965891   +0.4619226497   +0.4559704891
## [25]   -0.0000000000   -0.0000000000   +0.0000000000   +0.0000000000
```

## Subspace Restrictions

Instead of orthogonalizing we can also restrict $X$ to be in the subspace
of all lower triangular column-centered $n\times p$ matrices (which means $x_{ij}=0$
for all $i<j$). There are two different ways to accomplish this.

Method one uses a rotation of $X$ to lower triangular form. The R statement
we use for an nobj by ndim matrix, with nobj ≥ ndim, is


``` r
x <- x %*% qr.Q(qr(t(x[1:ndim, ])))
```
The theory
is pretty much the same as for the rotation to principal components
in the previous section.


In iteration 55, the final iteration, stress is 2.1114112739076. The "change" $\eta(X^{(k)}-X^{(k+1)})$ is 3.84673519544424e-16 and the estimate
of the asymptotic convergence ratio, the "change" divided by the "change" of the 
previous iteration, is 0.766987804354728.

The results are the same as for the basis sequence, as predicted. The eigenvalues of the Jacobian are

```
##  [1]   +0.7635162700   -0.6975184340   +0.6906648816   +0.6503675209
##  [5]   +0.6403092768   +0.6403092768   +0.6379438154   -0.6233133778
##  [9]   -0.6190294874   -0.5965727434   +0.5932846452   +0.5925642059
## [13]   -0.5809951380   +0.5658604713   +0.5431873999   -0.5421411727
## [17]   -0.5240341622   +0.5159258010   -0.4977034976   -0.4947337063
## [21]   -0.4839872579   -0.4779886296   -0.4535039166   -0.0506859167
## [25]   +0.0000000000   +0.0000000000   +0.0000000000   +0.0000000000
```

```
##  [1]   +0.7635162700   -0.6975184340   +0.6906648816   +0.6503675209
##  [5]   +0.6403092768   +0.6403092768   +0.6379438153   -0.6233133777
##  [9]   -0.6190294875   -0.5965727434   +0.5932846453   +0.5925642057
## [13]   -0.5809951380   +0.5658604713   +0.5431873999   -0.5421411727
## [17]   -0.5240341622   +0.5159258010   -0.4977034976   -0.4947337063
## [21]   -0.4839872580   -0.4779886296   -0.4535039166   -0.0506859167
## [25]   +0.0000000000   +0.0000000000   -0.0000000000   +0.0000000000
```
Again, the unit "rotation" eigenvalue from the unrotated solution has been replaced by zeroes.

Method two uses the theory of constrained smacof from @deleeuw_heiser_C_80. This
means computing the Guttman update and then projecting it on the subspace of
lower triangular matrices. We create $p$ column-centered matrices $Y_s$, of dimension $n\times(n-s)$, that satisfy $Y_s'VY_s=I$ and have their first $s-1$ rows equal to zero. Now column $s$ of $X$ is restricted to be of the form $x_s=Y_s\theta_s$. 
$$
x_s^{(k+1)}=Y_sY_s'V\Phi(X^{(k)})e_s
$$

In iteration 443, the final iteration, stress is 2.11141127390763. The "change" $\eta(X^{(k)}-X^{(k+1)})$ is 5.95185745918126e-16 and the estimate
of the asymptotic convergence ratio, the "change" divided by the "change" of the 
previous iteration, is 0.962237154391956.


```
##  [1]   +0.9622371565   +0.7669637116   +0.7480360811   +0.7105397440
##  [5]   +0.7007452013   +0.6920054108   +0.6859241126   +0.6593334167
##  [9]   +0.6541666413   +0.6476391845   +0.6234727079   +0.6172949237
## [13]   +0.5586355758   +0.5478680315   +0.5260205816   +0.5111279921
## [17]   +0.5064290866   +0.4929119272   +0.4843668597   +0.4783221142
## [21]   +0.4758041101   +0.4686440280   +0.4665512407   +0.4575763680
## [25]   -0.0000000000   -0.0000000000   -0.0000000000   +0.0000000000
```

```
##  [1]   +0.0687312255   +0.0547831223   +0.0534311487   +0.0507528389
##  [5]   +0.0500532287   +0.0494289579   +0.0489945795   +0.0470952440
##  [9]   +0.0467261887   +0.0462599417   +0.0445337648   +0.0440924945
## [13]   +0.0399025411   +0.0391334308   +0.0375728987   +0.0365091423
## [17]   +0.0361735062   +0.0352079948   +0.0345976328   +0.0341658653
## [21]   +0.0339860079   +0.0334745734   +0.0333250886   +0.0326840263
## [25]   -0.0000000000   +0.0000000000   +0.0000000000   +0.0000000000
```
# Two Point Iteration

## Relaxed

@deleeuw_heiser_C_80 suggested the "relaxed" update
\begin{equation}
\Psi(X):=2\Phi(X)-X.
(\#eq:relax)
\end{equation}
The reasoning for using \@ref(eq:relax) is two-fold. First, the smacof inequality \@ref(eq:smacofinequality) says
\begin{equation}
\sigma(X)\leq 1+\eta^2(X-\Phi(Y))-\eta^2(\Phi(Y)).
(\#eq:smaineq)
\end{equation}
If $X=\alpha\Phi(Y)+(1-\alpha)Y$ then this becomes
\begin{equation}
\sigma(\alpha\Phi(Y)+(1-\alpha)Y)\leq 1+(1-\alpha)^2\eta^2(Y-\Phi(Y))-\eta^2(\Phi(Y))
\end{equation}
If $(1-\alpha)^2\leq 1$ then 
\begin{equation}
1+(1-\alpha)^2\eta^2(Y-\Phi(Y))-\eta^2(\Phi(Y))\leq 1+\eta^2(Y-\Phi(Y))-\eta^2(\Phi(Y))=\sigma(Y)
\end{equation}
Thus updating with $X^{(k+1)}=\alpha\Phi(X^{(k)})+(1-\alpha)X^{(k)}$ is a stricly
monotone algorithm as long as $0\leq\alpha\leq 2$.

The second reason for choosing the relaxed update \@ref(eq:relax) given by @deleeuw_heiser_C_80
is that its asymptotic convergence rate is 
\begin{equation}
\max_s|2\lambda_s-1|=\max(2\lambda_{\text{max}}-1,1-2\lambda_{\text{min}}).
\end{equation}
@deleeuw_heiser_C_80 then somewhat carelessly assume that this is equal to
$2\lambda_{\text{max}}-1$ and argue that if $\lambda_{\text{max}}=1-\epsilon$
with $\epsilon$ small, as it usually is in MDS,  then 
\begin{equation}
2\lambda_{\text{max}}-1=1-2\epsilon\approx(1-\epsilon)^2=\lambda_{\text{max}}^2,
\end{equation}
so that the relaxed update requires approximately half the number of iterations of the
basic update. Despite the somewhat sloppy reasoning, the approximate halving of the number of iterations is often observed in practice. 

It turns out (@deleeuw_R_06b), however, that applying the relaxed update
has some unintended consequences, which basically imply that it should never
be used without additional precautions. Let's take a look at the
Ekman results.

```
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000
```
In iteration 25, the final iteration, stress is 2.1114112739076. The "change" $\eta(X^{(k)}-X^{(k+1)})$ is 3.76643158532133 and the estimate
of the asymptotic convergence ratio, the "change" divided by the "change" of the 
previous iteration, is 1.

The loss function value decreases and the number of iterations is reduced from 57 to 23.
But we see that $\eta(X^{(k+1)}-X^{(k)})$ does not converge to zero, and that $\sigma_k$ converges to a value which does not correspond to a local minimum of $\sigma$.

The eigenvalues of the Jacobian are 




If we check the conditions of theorem 3.1 in @meyer_76 we see that, although the 
algorithmic map is closed and the iterates are in a compact set, $\Psi$
is not strictly monotone at some non-fixed points. Suppose $X$ is a fixed point and  $\tau\not= 1$. Then
$\tau\overline{X}$ is not a fixed point of $\Psi$, because 
$\Psi(\tau\overline{X})=(2-\tau)\overline{X}$. And
\begin{equation}
\sigma(\tau\overline{X})=1-\tau\rho(\overline{X})+\frac12\tau^2\eta^2(\overline{X})=
1-\frac12\tau(2-\tau)\rho(\overline{X})=\sigma((2-\tau)\overline{X})
(\#eq:sigmatau)
\end{equation}
Thus the algorithm has convergent subsequences which may not converge to a fixed
point of $\Psi$ (and thus of $\Phi$). And indeed, the computational results show that the method produces a sequence $X^{(k)}$ with two subseqences. If $\overline{X}$ is a fixed point of $\Phi$ then there is a $\tau>0$ such that
the subsequence with $k$ even converges to $\tau\overline{X}$
while the subsequence with $k$ odd converges to $(2-\tau)\overline{X}$.

This suggests a simple fix. After convergence of the funcion values we make
a final update using $\Phi$ instead of $\Psi$. Computationally this is simple to do. If the final iteration updates $X^{(k)}$ to $X^{(k+1)}=\Psi(X^{(k)})$ then
set the final solution to the average $\frac12(X^{(k)}+X^{(k+1)})$. Making this 
adjustment at the end of the Ekman sequence gives us a final stress equal to 2.11141127390763.

## Doubling

The analysis in the previous section suggest the update function $\Psi^2$, i.e.
$$
X^{(k+1)}=\Psi(\Psi(X^{(k)}))=2\Phi(2\Phi(X^{k})-X^{{k)}})-\Psi(X^{(k)}).
$$
$$
\mathcal{D}\Psi^2_X(H)=\mathcal{D}\Psi_{\Psi(X)}(\mathcal{D}\Psi_X(H))
$$
The asymptotic convergence rate is 
$$
\max_s (2\lambda_s-1)^2=\max\{ (2\lambda_\text{max}-1)^2, (2\lambda_\text{min}-1)^2\}
$$.

With $\Psi^2$ the algorithm is everywhere strictly monotonic and
converges to a fixed point. But not all problems have disappeared.
If $X$ is a stationary point of $\sigma$, and thus a fixed 
point of $\Phi$, $\Psi$, and $\Psi^2$, then $\tau X$ is a 
fixed point of $\Psi^2$ for all $\tau>0$. Thus we cannot
exclude the possibility that the sequence converges to
a fixed point proportional to $X$, but not equal to $X$.

Here are the results for the Ekman data if we use $\Psi^2$.

```
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992
```
Again we need some adjustment. A final update using $\Phi$ will do the trick.


stress is 2.1114112739076

## Dilation

@deleeuw_R_06b discusses some other ways to fix the relaxed update problem. 
The first one, borrowed from @groenen_glunt_hayden_96, defines
Defines
$$
\Pi(X):=\frac{\rho(X)}{\eta^2(X)}X
$$
and 
$$
\Xi(X):=\Pi(\Psi(X))
$$

First, differentiate $\Pi$ of ... Using the product and quotient rules for differentiation we find
$$
\mathcal{D}\Pi_X(H)=\frac{\rho(X)}{\eta^2(X)}H+\frac{\eta^2(X)\mathcal{D}\rho_X(H)-\rho(X)\mathcal{D}\eta^2_X(H)}{\eta^4(X)}X
$$
Using
$$
\mathcal{D}\rho_X(H)=\text{tr}\ H'B(X)X
$$
$$
\mathcal{D}\eta^2_X(H)=2\text{tr}\ H'VX
$$
this becomes
$$
\mathcal{D}\Pi_X(H)=\frac{\rho(X)}{\eta^2(X)}H+\text{tr}\ H'\left\{\frac{\eta^2(X) B(X)X-2\rho(X)VX}{\eta^4(X)}\right\}X
$$


The chain rule says
$$
\mathcal{D}\Xi_X(H)=\mathcal{D}\Psi_X(H)-\frac{\text{tr}\ X'V\mathcal{D}\Psi_X(H)}{\text{tr}\ X'VX}X
$$
Since $\mathcal{D}\Psi_X(X)=-X$ we have
$$
\mathcal{D}\Xi_X(X)=-X+\frac{\text{tr}\ X'VX}{\text{tr}\ X'VX}X=0
$$
Thus the offending eigenvector $X$ of $\mathcal{D}\Psi$ is eliminated.

More generally, if $\mathcal{D}\Psi_X(H)=\lambda H$ with $H\not= X$ then
$\mathcal{D}\Xi_X(H-X)=\lambda (H-X),$
and thus $\mathcal{D}\Xi$ has the same eigenvalues as $\mathcal{D}\Psi$.


## Stabilizing

Another strategy 

\begin{equation}
\Xi(X):=\Phi(\Psi(X))
\end{equation}

$$
\mathcal{D}\Xi_X(H)=\mathcal{D}\Phi_{\Psi(X)}(\mathcal{D}\Psi_X(H))=2\mathcal{D}\Phi_{\Psi(X)}(\mathcal{D}\Phi_X(H))-\mathcal{D}\Phi_{\Psi(X)}(H)
$$

\begin{equation}
\max_s|\lambda_s(2\lambda_s-1)|
\end{equation}

# Benchmarking

We compare the eight different upgrades using the microbenchmark package (@mersmann_23).


```
## Warning in microbenchmark(smacofAccelerate(delta, ndim = 2, opt = 1, halt = 2,
## : less accurate nanosecond times to avoid potential integer overflows
```

```
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   13 sold 3.994627066568262 snew 2.111411273907599 chng 0.000000000000000 labd 0.273802752120992 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000 
## itel   25 sold 3.994627066568261 snew 2.111411273907601 chng 3.766431585321326 labd 1.000000000000000
```

```
## Unit: milliseconds
##                                                                   expr
##  smacofAccelerate(delta, ndim = 2, opt = 1, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 2, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 3, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 4, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 5, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 6, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 7, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 8, halt = 2, verbose = FALSE)
##        min        lq      mean    median        uq       max neval
##   3.233957  3.311918  3.815434  3.366777  3.471614 27.450566   100
##   4.345713  4.475683  4.754933  4.560061  4.667891  6.656883   100
##   4.018164  4.148072  4.462280  4.218900  4.319801  6.406209   100
##  26.918304 28.833947 29.255858 29.083268 29.595009 31.553272   100
##   1.658122  1.724501  1.910976  1.763000  1.819580  3.760766   100
##   1.291459  1.349535  1.457766  1.381126  1.422126  3.291521   100
##   1.687396  1.745944  1.992141  1.779400  1.841207  6.327038   100
##   1.650045  1.712324  1.966379  1.742254  1.802032  6.875741   100
```

@degruijter_67


```
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  400 sold 246.696086021759385 snew 128.883258119202168 chng 235.625655805114718 labd 1.000000000000000 
## itel  209 sold 246.696086021754525 snew 128.883258119199184 chng 0.000000000000124 labd 0.942946332991647
```

```
## Unit: milliseconds
##                                                                   expr      min
##  smacofAccelerate(delta, ndim = 2, opt = 1, halt = 2, verbose = FALSE) 44.95830
##  smacofAccelerate(delta, ndim = 2, opt = 2, halt = 2, verbose = FALSE) 62.30085
##  smacofAccelerate(delta, ndim = 2, opt = 3, halt = 2, verbose = FALSE) 58.83324
##  smacofAccelerate(delta, ndim = 2, opt = 4, halt = 2, verbose = FALSE) 57.41406
##  smacofAccelerate(delta, ndim = 2, opt = 5, halt = 2, verbose = FALSE) 20.40533
##  smacofAccelerate(delta, ndim = 2, opt = 6, halt = 2, verbose = FALSE) 14.94425
##  smacofAccelerate(delta, ndim = 2, opt = 7, halt = 2, verbose = FALSE) 23.74761
##  smacofAccelerate(delta, ndim = 2, opt = 8, halt = 2, verbose = FALSE) 21.55964
##        lq     mean   median       uq      max neval
##  46.78315 47.49755 47.22993 48.29044 50.84037   100
##  63.74924 65.67610 64.68295 66.07109 85.66126   100
##  59.44854 61.06197 60.55741 61.97714 79.67235   100
##  58.31215 60.39281 59.24900 60.91655 82.25022   100
##  21.96175 22.47479 22.34502 22.70955 43.00720   100
##  15.40913 16.39808 16.63128 17.01777 18.57173   100
##  25.20951 25.59745 25.45491 25.81770 28.72735   100
##  23.47836 23.72634 23.74736 24.06579 25.78851   100
```


```
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  103 sold 198.362049456548306 snew 78.921357319073834 chng 0.000000000000140 labd 0.875686385683123 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000 
## itel  206 sold 198.362049456547709 snew 78.921357319074076 chng 238.881384274947237 labd 1.000000000000000
```

```
## Unit: milliseconds
##                                                                   expr
##  smacofAccelerate(delta, ndim = 2, opt = 1, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 2, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 3, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 4, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 5, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 6, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 7, halt = 2, verbose = FALSE)
##  smacofAccelerate(delta, ndim = 2, opt = 8, halt = 2, verbose = FALSE)
##        min        lq      mean    median        uq       max neval
##   43.98095  47.13093  48.67453  47.61691  48.15760  72.72400   100
##   55.83741  57.50125  59.05327  58.30667  59.90717  84.83490   100
##   53.61685  54.74586  55.92418  55.44215  56.67172  61.37659   100
##  108.64250 110.47009 112.86205 111.48039 114.11499 135.28438   100
##   19.27976  19.87936  21.42617  20.32499  22.88651  46.50650   100
##   14.78153  15.42715  17.79269  15.93697  18.55365  43.31142   100
##   24.85149  25.83045  28.16856  28.41765  28.96521  51.81096   100
##   24.00677  25.20598  28.19073  27.72535  28.41659  50.75968   100
```

# (APPENDIX) Appendices {-} 

# Rotation to Principal Components {#rotation}

We transform $X$ to $\Pi(X)=XL$, where $L$ are the right singular vectors of $X$, or
equivalently the eigenvectors of $X'X$.
\begin{equation}
(X+\Delta_X)'(X+\Delta_X)(L+\Delta_L)=(L+\Delta_L)(\Lambda+\Delta_\Lambda)
\end{equation}
Only keep first order terms
$$
X'X\Delta_L+\Delta_X'XL+X'\Delta_XL=L\Delta_\Lambda+\Delta_L\Lambda
$$
Pfremultiply by $L'$. Then
$$
\Lambda L'\Delta_L+L'(\Delta_X'X+X'\Delta_X)L=\Delta_\Lambda+L'\Delta_L\Lambda
$$
Define $S=L'\Delta_L$. Then
$$
L'(\Delta_X'X+X'\Delta_X)L-\Delta_\Lambda=-(\Lambda S-S\Lambda)
$$
Also, of course,
$$
(L+\Delta_L)'(L+\Delta_L)=I
$$
and thus
$$
S+S'=0
$$
i.e. $S$ is anti-symmetric and has a zero diagonal. Taking the diagonal of both
sides of .. gives
$$
\Delta_\Lambda=\text{diag}(L'(\Delta_X'X+X'\Delta_X)L)
$$
Taking the non-diagonal gives
$$
s_{ij}=-\frac{l_i'(\Delta_X'X+X'\Delta_X)l_j}{\lambda_i-\lambda_j}
$$
Now
$$
\mathcal{D}\Pi_X(\Delta_X)=\Delta_XL+X\Delta_L=\Delta_XL+XLS
$$
To compute partial derivatives of $\Pi$ with repect to $x_{is}$
we simply take $\Delta_X=e_ie_s'$, with $e_i$ and $e_s$ unit 
vectors of length $n$ and $m$. Unit vectors have one element equal to
one in the place indicated by their subscript, all other elements
are zero.


# Rotation to Lower Triangular

$$
ut{(X+H)(K+\Delta_K)}=0
$$
$$
ut(X\Delta_K)=-ut(HK)
$$
$K'\Delta_K+\Delta_K'K=0$

$\Delta_K=KS$

Solve
$$
ut(X_1KS)=-ut(H_1K)
$$
For all $i<j$
$$
\{X_1KS\}_{ij}=\{H_1K\}_{ij}
$$
$$
S=\sum_{k<l}\alpha_{kl}(e_ke_l'-e_le_k')
$$
$U=X_1K$m $V=H_1K$

Thus for $i<j$
$$
\sum_{k<l}\alpha_{kl}\{UE_{kl}\}_{ij}=v_{ij}
$$
which is a linear system in $\frac12n(n-1)$ unknowns.

$$
\{UE_{kl}\}_{ij}=e_i'\{Ue_ke_l'-Ue_le_k'\}e_j=u_{ik}\delta^{jl}-u_{il}\delta^{jk}
$$
If $p=2$ then $(i,j)=(k,l)=(1,2)$. Thus
$$
\alpha_{12}u_{11}=v_{12}
$$
and thus
$$
S=\begin{bmatrix}0&v_{12}/u_{11}\\-v_{12}/u_{11}&0\end{bmatrix}
$$
and
$$
\mathcal{D}\Pi_X(H)=HK+XKS
$$

# Code

## smacofAccelerate.R


``` r
library(MASS)
library(microbenchmark)
library(numDeriv)

source("smacofUtils.R")
source("smacofDerivatives.R")

smacofAccelerate <- function(delta,
                             wgth = 1 - diag(nrow(delta)),
                             ndim = 2,
                             xold = smacofTorgerson(delta, ndim),
                             opt = 1,
                             halt = 0,
                             wd = 4,
                             dg = 15,
                             itmax = 1000,
                             epsx = 1e-10,
                             epsf = 1e-15,
                             verbose = 2) {
  nobj <- nrow(delta)
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  if (opt == 4) {
    bs <- smacofMakeBasis(nobj, ndim, vmat)
  }
  xold <- smacofCenter(xold)
  if ((opt == 2) || (opt == 4)) {
    xrot <- qr.Q(qr(xold[1:ndim, ]))
    xold <- xold %*% xrot
  }
  if (opt == 3) {
    xrot <- svd(xold)$v
    xold <- xold %*% xrot
  }
  dold <- as.matrix(dist(xold))
  sold <- sum(wgth * (delta - dold) ^ 2)
  cold <- Inf
  itel <- 1
  repeat {
    if (opt == 1) {
      h <- smacofOptionOne(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 2) {
      h <- smacofOptionTwo(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 3) {
      h <- smacofOptionThree(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 4) {
      h <- smacofOptionFour(xold, delta, wgth, vmat, vinv, bs)
    }
    if (opt == 5) {
      h <- smacofOptionFive(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 6) {
      h <- smacofOptionSix(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 7) {
      h <- smacofOptionSeven(xold, delta, wgth, vmat, vinv)
    }
    if (opt == 8) {
      h <- smacofOptionEight(xold, delta, wgth, vmat, vinv)
    }
    labd <- sqrt(h$cnew / cold)
    if (verbose == 2) {
      smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
    }
    if (halt == 1) {
      converge <- h$cnew < epsx
    } else {
      converge <- (sold - h$snew) < epsf
    }
    if ((itel == itmax) || converge) {
      break
    }
    itel <- itel + 1
    sold <- h$snew
    xold <- h$xnew
    cold <- h$cnew
    dold <- h$dnew
  } # end of repeat loop
  if ((verbose == 1) && (opt != 5) && (opt != 6)) {
    smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
  }
  if (opt == 5) {
    h$xnew <- (h$xnew + xold) / 2
    h$dnew <- as.matrix(dist(h$xnew))
    h$snew <- sum(wgth * (delta - h$dnew) ^ 2)
    smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
  }
  if (opt == 6) {
    bold <- -wgth * delta / (h$dnew + diag(nobj))
    diag(bold) <- -rowSums(bold)
    h$xnew <- vinv %*% bold %*% h$xnew
    h$dnew <- as.matrix(dist(h$xnew))
    h$snew <- sum(wgth * (delta - h$dnew) ^ 2)
    smacofLinePrint(itel, sold, h$snew, h$cnew, labd, wd = wd, dg = dg)
  }
  return(
    list(
      x = h$xnew,
      s = h$snew,
      d = h$dnew,
      itel = itel,
      chng = h$cnew,
      labd = labd,
      wgth = wgth,
      delta = delta
    )
  )
}

smacofOptionOne <- function(xold, delta, wgth, vmat, vinv) {
  xnew <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionTwo <- function(xold, delta, wgth, vmat, vinv) {
  ndim <- ncol(xold)
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xrot <- smacofSignEigenVectors(qr.Q(qr(t(xbar[1:ndim, ]))))
  #xrot <- qr.Q(qr(t(xbar[1:ndim, ])))
  xnew <- xbar %*% xrot
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionThree <- function(xold, delta, wgth, vmat, vinv) {
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xrot <- smacofSignEigenVectors(svd(xbar)$v)
  xnew <- xbar %*% xrot
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionFour <- function(xold, delta, wgth, vmat, vinv, bs) {
  ndim <- ncol(xold)
  nobj <- nrow(xold)
  xnew <- matrix(0, nobj, ndim)
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  for (s in 1:ndim) {
    aux <- crossprod(bs[[s]], vmat %*% xbar[, s])
    xnew[, s] <- bs[[s]] %*% aux
  }
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionFive <- function(xold, delta, wgth, vmat, vinv) {
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xnew <- 2 * xbar - xold
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionSix <- function(xold, delta, wgth, vmat, vinv) {
  nobj <- nrow(xold)
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  baux <- -wgth * delta / (daux + diag(nobj))
  diag(baux) <- -rowSums(baux)
  xbaz <- vinv %*% baux %*% xaux
  xnew <- 2 * xbaz - xaux
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionSeven <- function(xold, delta, wgth, vmat, vinv) {
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  alpa <- sum(wgth * daux * delta) / sum(wgth * daux ^ 2)
  xnew <- alpa * xaux
  dnew <- alpa * daux
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}

smacofOptionEight <- function(xold, delta, wgth, vmat, vinv) {
  nobj <- nrow(xold)
  xbar <- smacofCenter(smacofGuttman(xold, delta, wgth, vinv))
  xaux <- 2 * xbar - xold
  daux <- as.matrix(dist(xaux))
  baux <- -wgth * delta / (daux + diag(nobj))
  diag(baux) <- -rowSums(baux)
  xnew <- vinv %*% baux %*% xaux
  dnew <- as.matrix(dist(xnew))
  snew <- sum(wgth * (delta - dnew) ^ 2)
  cnew <- sum((xold - xnew) * (vmat %*% (xold - xnew)))
  return(list(
    xnew = xnew,
    dnew = dnew,
    snew = snew,
    cnew = cnew
  ))
}
```

## smacofDerivatives.R


``` r
library(numDeriv)

source("smacofUtils.R")
source("smacofPCADerivative.R")
source("smacofQRDerivative.R")

smacofRhoHessian <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  fac1 <- wgth * delta / (dmat + diag(nobj))
  fac2 <- wgth * delta / ((dmat + diag(nobj)) ^ 3)
  bmat <- -fac1
  diag(bmat) <- -rowSums(bmat)
  hess <- matrix(0, ntot, ntot)
  for (s in 1:ndim) {
    ns <- (s - 1) * nobj + 1:nobj
    hess[ns, ns] <- bmat
    for (t in 1:ndim) {
      nt <- (t - 1) * nobj + 1:nobj
      ds <- outer(x[, s], x[, s], "-")
      dt <- outer(x[, t], x[, t], "-")
      aux <- -fac2 * ds * dt
      diag(aux) <- -rowSums(aux)
      hess[ns, nt] <- hess[ns, nt] - aux
    }
  }
  return(hess)
}

smacofBasicDerivative <- function(x, h, delta, wgth, vinv, dmat) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  bmat <- wgth * delta / (dmat + diag(nobj))
  bmat <- -bmat
  diag(bmat) <- -rowSums(bmat)
  hmat <- wgth * delta / ((dmat + diag(nobj)) ^ 3)
  for (i in 1:nobj) {
    for (j in 1:nobj) {
      xhij <- sum((x[i, ] - x[j, ]) * (h[i, ] - h[j, ]))
      hmat[i, j] <- hmat[i, j] * xhij
    }
  }
  hmat <- -hmat
  diag(hmat) <- -rowSums(hmat)
  deri <- vinv %*% (bmat %*% h - hmat %*% x)
  return(deri)
}

smacofBasicJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofBasicJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    return(as.vector(xbar))
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth
  )
  return(jacob)
}

smacofPCADerivative <- function(x, h, delta, wgth, vinv, dmat) {
  xbar <- smacofGuttman(x, delta, wgth, vinv)
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- PCADerivative(xbar, dexh)
  return(deri)
}

smacofPCAJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofPCADerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofPCAJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    l <- smacofSignEigenVectors(svd(xbar)$v)
    return(xbar %*% l)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth
  )
  return(jacob)
}


smacofQRDerivative <- function(x, h, delta, wgth, vinv, dmat) {
  xbar <- smacofGuttman(x, delta, wgth, vinv)
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- QRDerivative(xbar, dexh)
  return(deri)
}

smacofQRJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofQRDerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofQRJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    l <- qr.Q(qr(t(xbar[1:ndim, ])))
    return(xbar %*% l)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth
  )
  return(jacob)
}

smacofYbasDerivative <- function(x, h, delta, wgth, vmat, vinv, dmat, bs) {
  ndim <- ncol(x)
  nobj <- nrow(x)
  ntot <- nobj * ndim
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- matrix(0, nobj, ndim)
  for (i in 1:ndim) {
    deri[, i] <- bs[[i]] %*% crossprod(bs[[i]], vmat %*% dexh[, i])
  }
  return(deri)
}

smacofYbasJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  bs <- smacofMakeBasis(nobj, ndim, vmat)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofYbasDerivative(x, h, delta, wgth, vmat, vinv, dmat, bs)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofYbasJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  bs <- smacofMakeBasis(nobj, ndim, vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth, vinv, bs) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    for (i in 1:ndim) {
      xbar[, i] <- bs[[i]] %*% crossprod(bs[[i]], xbar[, i])
    }
    return(xbar)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth,
    vinv = vinv,
    bs = bs
  )
  return(jacob)
}


smacofRelaxDerivative <- function(x, h, delta, wgth, vinv, dmat) {
  ndim <- ncol(x)
  nobj <- nrow(x)
  ntot <- nobj * ndim
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- 2 * dexh - h
  return(deri)
}

smacofRelaxJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofRelaxDerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofRelaxJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  bs <- smacofMakeBasis(nobj, ndim, vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth, vinv) {
    x <- matrix(x, nobj, ndim)
    xbar <- smacofGuttman(x, delta, wgth, vinv)
    xbaz <- 2 * xbar - x
    return(xbaz)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth,
    vinv = vinv
  )
  return(jacob)
}

smacofDoubleDerivative <- function(x, h, delta, wgth, vinv, dmat) {
  ndim <- ncol(x)
  nobj <- nrow(x)
  ntot <- nobj * ndim
  dexh <- smacofBasicDerivative(x, h, delta, wgth, vinv, dmat)
  deri <- 2 * dexh - h
  return(deri)
}

smacofDoubleJacobianFormula <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  ntot <- nobj * ndim
  dmat <- as.matrix(dist(x))
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  jacob <- matrix(0, ntot, ntot)
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- smacofDoubleDerivative(x, h, delta, wgth, vinv, dmat)
      jacob[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(jacob)
}

smacofDoubleJacobianNumerical <- function(x, delta, wgth) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  vmat <- -wgth
  diag(vmat) <- -rowSums(vmat)
  bs <- smacofMakeBasis(nobj, ndim, vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  func <- function(x, nobj, ndim, delta, wgth, vinv) {
    x <- matrix(x, nobj, ndim)
    xbar <- 2 * smacofGuttman(x, delta, wgth, vinv) - x
    xbaz <- 2 * smacofGuttman(xbar, delta, wgth, vinv) - xbar
    return(xbaz)
  }
  jacob <- jacobian(
    func,
    as.vector(x),
    nobj = nobj,
    ndim = ndim,
    delta = delta,
    wgth = wgth,
    vinv = vinv
  )
  return(jacob)
}
```

## smacofPCADerivative.R


``` r
PCADerivative <- function(x, h) {
  ndim <- ncol(x)
  eixx <- eigen(crossprod(x))
  evec <- eixx$vectors
  evec <- evec %*% diag(sign(diag(evec)))
  eval <- eixx$values
  xh <- crossprod(x, h)
  xh <- xh + t(xh)
  s <- matrix(0, ndim, ndim)
  for (i in 1:ndim) {
    for (j in 1:ndim) {
      if (i == j) {
        next
      }
      s[i, j] <- -sum(evec[, i] * (xh %*% evec[, j])) / (eval[i] - eval[j])
    }
  }
  return(h %*% evec + x %*% evec %*% s)
}

PCAJacobianFormula <- function(x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  np <- nobj * ndim
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  d <- matrix(0, np, np)
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- PCADerivative(x, h)
      d[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(d)
}

PCAJacobianNumerical <- function(x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  func <- function(x, nobj, ndim) {
    x <- matrix(x, nobj, ndim)
    l <- svd(x)$v
    l <- l %*% diag(sign(diag(l)))
    return(x %*% l)
  }
  jacob <- jacobian(func, as.vector(x), nobj = nobj, ndim = ndim)
  return(jacob)
}

PCATester <- function(x, h) {
  func <- function(x) {
    l <- svd(x)$v
    l <- l %*% diag(sign(diag(l)))
    return(x %*% l)
  }
  x0 <- func(x)
  xh <- func(x + h)
  xd <- x0 + PCADerivative(x, h)
  print(cbind(x0, xh, xd))
}
```

## smacofQRDerivative.R


``` r
QRDerivative <- function(x, h) {
  ndim <- ncol(x)
  nn <- ndim * (ndim - 1) / 2
  x1 <- x[1:ndim, ]
  h1 <- h[1:ndim, ]
  kk <- qr.Q(qr(t(x1)))
  uu <- x1 %*% kk
  vv <- -h1 %*% kk
  amat <- matrix(0, nn, nn)
  bmat <- rep(0, nn)
  indi <- 1
  for (i in 1:(ndim - 1)) {
    for (j in (i + 1):ndim) {
      bmat[indi] <- vv[i, j]
      jndi <- 1
      for (k in 1:(ndim - 1)) {
        for (l in (k + 1):ndim) {
          if (j == l) {
            amat[indi, jndi] <- amat[indi, jndi] + uu[i, k]
          }
          if (j == k) {
            amat[indi, jndi] <- amat[indi, jndi] - uu[i, l]
          }
          jndi <- jndi + 1
        }
      }
      indi <- indi + 1
    }
  }
  coef <- solve(amat, bmat)
  indi <- 1
  s <- matrix (0, ndim, ndim)
  for (i in 1:(ndim - 1)) {
    for (j in (i + 1):ndim) {
      s[i, j] <- coef[indi]
      s[j, i] <- -coef[indi]
      indi <- indi + 1
    }
  }
  deri <- h %*% kk + x %*% kk %*% s
  return(deri)
}

QRJacobianFormula <- function(x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  np <- nobj * ndim
  e <- function(i, n) {
    ifelse(i == 1:n, 1, 0)
  }
  d <- matrix(0, np, np)
  k <- 1
  for (j in 1:ndim) {
    for (i in 1:nobj) {
      h <- outer(e(i, nobj), e(j, ndim))
      r <- QRDerivative(x, h)
      d[, k] <- as.vector(r)
      k <- k + 1
    }
  }
  return(d)
}

QRJacobianNumerical <- function(x) {
  nobj <- nrow(x)
  ndim <- ncol(x)
  func <- function(x, nobj, ndim) {
    x <- matrix(x, nobj, ndim)
    l <- qr.Q(qr(t(x[1:ndim, ])))
    return(x %*% l)
  }
  jacob <- jacobian(func, as.vector(x), nobj = nobj, ndim = ndim)
  return(jacob)
}

QRTester <- function(x, h) {
  func <- function(x) {
    ndim <- ncol(x)
    l <- qr.Q(qr(t(x[1:ndim, ])))
    l <- l %*% diag(sign(diag(l)))
    return(x %*% l)
  }
  x0 <- func(x)
  xh <- func(x + h)
  xd <- x0 + QRDerivative(x, h)
  print(cbind(x0, xh, xd))
}
```

## smacofCompare.R


``` r
smacofCompare <- function(delta, ndim = 2) {
  nobj <- nrow(delta)
  wgth <- 1 - diag(nobj)
  xold <- smacofTorgerson(delta, ndim)
  return(
    microbenchmark(
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 1,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 2,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 3,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 4,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 5,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 6,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 7,
        halt = 2,
        verbose = FALSE
      ),
      smacofAccelerate(
        delta,
        ndim = 2,
        opt = 8,
        halt = 2,
        verbose = FALSE
      )
    )
  )
}
```

## smacofUtils.R


``` r
smacofMatrixPrint <- function(x,
                   digits = 10,
                   width = 15,
                   format = "f",
                   flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

smacofLinePrint <- function(itel, sold, snew, cnew, labd, wd, dg) {
  cat(
    "itel",
    formatC(itel, width = wd, format = "d"),
    "sold",
    formatC(sold, digits = dg, format = "f"),
    "snew",
    formatC(snew, digits = dg, format = "f"),
    "chng",
    formatC(cnew, digits =  dg, format = "f"),
    "labd",
    formatC(labd, digits =  dg, format = "f"),
    "\n"
  )
}

smacofMakeBasis <- function(n, ndim, vmat) {
  y <- lapply(1:ndim, function(k)
    matrix(0, n, n - k))
  for (s in 0:(ndim - 1)) {
    ns <- n - s
    aux <- qr.Q(qr(ns * diag(ns) - 1))[, -ns]
    aux <- rbind(matrix(0, s, ns - 1), aux)
    sux <- crossprod(aux, vmat %*% aux)
    y[[s + 1]] <- aux %*% smacofMatrixPower(sux, -0.5)
  }
  return(y)
}

smacofTorgerson <- function(delta, ndim) {
  n <- nrow(delta)
  dd <- delta ^ 2
  rd <- rowSums(dd) / n
  sd <- mean(dd)
  cc <- -.5 * (dd - outer(rd, rd, "+") + sd)
  ee <- eigen(cc)
  x <- ee$vectors[, 1:ndim] %*% diag(sqrt(ee$values[1:ndim]))
  return(x)
}

smacofGuttman <- function(x, delta, wgth, vinv) {
  nobj <- nrow(x)
  dmat <- as.matrix(dist(x))
  bmat <- -wgth * delta / (dmat + diag(nobj))
  diag(bmat) <- -rowSums(bmat)
  return(vinv %*% bmat %*% x)
}

smacofCenter <- function(x) {
  return(apply(x, 2, function(x)
    x - mean(x)))
}

smacofSignEigenVectors <- function(x) {
  return(x %*% diag(sign(diag(x))))
}

smacofMatrixPower <- function(s, power) {
  e <- eigen(s)
  eval <- e$values
  evec <- e$vectors
  dval <- ifelse(abs(eval) < 1e-10, 0, abs(eval) ^ power)
  return(tcrossprod(evec %*% diag(dval), evec))
}
```

# References
