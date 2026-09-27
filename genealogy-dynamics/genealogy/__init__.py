"""Genealogy dynamics: two-sex cluster-dynamics model of US population pedigree.

Adapts the graph-based (RAG) cluster-dynamics method to human population
genealogy. Pedigree degree is a *generation cluster* (depth d = 0..n, the
"size" analog of the defect-cluster framework); the population is split into
two populations P_pat (paternal/male) and P_mat (maternal/female), both
subsets of the cluster species class C. The cluster state (d, alpha) holds
concentration c_{d,alpha}(t) (millions), and the state vector
c = [c_{0,pat}..c_{n,pat}, c_{0,mat}..c_{n,mat}] in R^{2(n+1)} evolves as

    dc/dt = P(t) + S J(c) - D(t) c

under immigration (sources), cross-population mating/births (catalytic
binary edges), death and emigration (sex-specific sinks).
"""
