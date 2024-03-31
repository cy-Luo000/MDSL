# Dense 2-club

## Background

Find the maximum dense subgraph such that $g(E(\overline{S}))\leq f(|S|)$, and the $f(x)$ and $g(x)$ are increasing functions.

## Solution

1. $k$-defective clique: If the size of maximum $k$-defective clique is larger than $k+1$, then $diam(S)\leq 2$

2. $k$-plex: If the size of maximum $k$-plex is at least $2k-1$, then $diam(S)\leq 2$

3. degree $\gamma$-quasi clique: If the degree density of maximum degree $\gamma$-quasi clique is larger than $0.5$, then $diam(S)\leq 2$

4. density $\gamma$-quasi clique: If $S$ is a density $\gamma$-quasi clique, then the diameter of $S$ satisfied $diam(S)\leq d_{\gamma}$, where $d_{\gamma}=\lfloor |S|+\frac{1}{2}-\sqrt{\gamma|S|^2-(2+\gamma)|S|+\frac{17}{4}}\rfloor$ 

   (Here we study the $(\gamma, \lambda)$-quasi clique) with $\lambda= 0.5$ and $\gamma \geq \lambda$. Why we let $\lambda=0.5$? Because we think in a community, each member should have connection with at least half of other members(If $\lambda<0.5$, the $(\gamma,\lambda)$-quasi clique may be series of clusters, which is not we want). Besides, we hope to mine a subgraph in which each member has many connections with other members, so we should avoid adding the vertex with low connections with other nodes into the graph).

## Decomposition

$N_2(v)$：表示距离点v距离小于等于2的组成的集合

## Other Properties



## Algorithms

1. Use a heuristic search to find a initial solution $lb$
2. Cut the vertices and the edges
3. Do the graph decomposition based on the 2-hop degeneracy
4. Do step 1 and 2 for each subgraph
5. Do the branch and bound algorithm for each subgraph

## Reductions

1. If $S^*=S\cup C$ and $f(V(S^*))\geq g(E(S^*))$, then $S^*=S\cup C$ is a $ve$ dense graph.
2. If $ub(S)\leq lb$ then there is no $ve$ dense graph $S^*$ such that $S\subset S^*$ and $|S^*|> lb$. 

### *Upper Bound

We have a color based sorting upper bound, following is the steps.

1. Partition the candidate set $C$ into $\Lambda$ subsets, namely $C=C_1\cap C_2\cap...\cap C_{\Lambda}$
2. For each subset $C_i(i=1,2,...,\Lambda)$, partition the $C_i$ into independent set, namely $C_i=\pi_i^1\cap \pi_i^2\cap...\cap \pi_i^{\chi_i}$. For the independent set $\pi_i^j=\{v_1,v_2,...,v_{|\pi_i^j|}\}$, and the vertex in $\pi_i^j$ is ordered by $\bar{d}_p(\cdot)$ in increasing order. We denote the weight of  the $k$-th vertex in $\pi_i^j$ , which is $u_k$, is $w(u_k)=\bar{d}_P(u_k)+k-1$. 
3. Sort all vertices in $C$ by $w(\cdot)$ in increasing order
4. If $|C'|=l_0$, then the lower bound of the number of missing edges in $P\cup C'$ is $|r(P)+\sum_{i=1}^{l_0}w(v_i)|$. 
5. According to the problem, get an upper bound.

**Note:**  There are many ways to do the partition, following we show 2 ways.

1. We partition the $C$ into 1 subset, then the process of upper bound is Chang's sorting bound.
2. We partition the $C$ into $C=C_0\cap C_1\cap ...\cap C_{|P|}$, such that the vertex $v$ in $C_i$ satisfies $\bar{d}_p(v)=i$

**property of the upper bound**

Assume the current instance is $I=(G,P,C)$, we should find a solution $S$ such that $P\subset S\subset P\cup C$. And the vertices in candidate $C$ ordered by $w(\cdot)$

by non-decreasing way is $Ord(C)=\{v_1,v_2,...,v_{|C|}\}$. For $\forall u,v\in C$, we denote $u\prec v$ if $v$ is after $u$ in the ordered $C$. We denote $H_i=\{v_1,v_2,...,v_i\}$.

***Theorm1.***   If $ub$ satisfies that $ub$ is the largest value satisfies $g(E(\overline{P\cup H_{ub}}))\leq f(|P\cup H_{ub}|)$, then the upper bound of $S$ is $|P|+ub$.

*proof.* We use the proof by contradiction. We assume that there exists $S'$ such that $S'$ satisfies the definition and $|S'|> |P|+ub$. We denote $H'=S'\setminus P$. We denote $H'_{ij}=H'\cap \pi_i^j=\{v_{ij}^1,v_{ij}^2,...,v_{ij}^{|H'_{ij}|}\}$ and $h_{ij}=|H'_{ij}|$, then $|\bar{E}(P,H'_{ij})|+|\bar{E}(H'_{ij})|=\sum_{k=1}^{h_{ij}}\bar{d_P}(v_{ij}^k)+\sum_{k=1}^{h_{ij}}(k-1)=\sum_{k=1}^{h_{ij}}(\bar{d_P}(v_{ij}^k)+k-1)\geq \sum_{k=1}^{h_{ij}}(\bar{d_P}(u_{ij}^k)+k-1)=\sum_{k=1}^{h_{ij}} w(u_{ij}^k)$. Then we have  $|\overline{E}(S')|=\bar{E}(P)+\bar{E}(P,H)+\bar{E}(H)$

