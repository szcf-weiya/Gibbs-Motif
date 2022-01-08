Suppose the posterior distribution of discrete variables $(x_1,\ldots,x_N)$ is 
$$
\begin{align*}
f(x_1,\ldots,x_N) &= p_0I(x_1=a_1,\ldots, x_N=a_N)+\\
&\qquad p_1I(x_1=.,\ldots,x_N=.)+\\
&\qquad \cdots
\end{align*}
$$
where $p_0 > p_1>\cdots$, and $p_0 + p_1+\cdots = 1$, so
$$
a_1,\ldots,a_N=\arg\max f(x_1, \ldots,x_N)\,.
$$
If $p_0 \ge 0.5$, then
$$
\Pr(x_i\neq a_i) = 1-\Pr(x_i=a_i) <0.5 < p_0
$$
which implies
$$
a_i = \arg\max f_i(x_i)\,.
$$
If $p_0 < 0.5$, there might exist $b_j\neq a_j$ such that
$$
\Pr(x_j = b_j) = 1-p_0 > 0.5
$$
then
$$
\arg\max f_j(x_j) \neq \{\arg\max f(x_1,\ldots,x_N)\}_j
$$
For example,
$$
f(x_1,x_2) = 0.3I(x_1=0,x_2=0) + 0.4I(x_1=1,x_2=1)+0.3I(x_1=1,x_2=0)
$$
then
$$
f_1(x_1) = 0.7I(x_1=1),\qquad f_2(x_2) = 0.4I(x_2=1)
$$
thus,
$$
0 = \arg\max f_2(x_2)\neq \{\arg\max f(x_1,x_2)\}_2 = 1\,.
$$
