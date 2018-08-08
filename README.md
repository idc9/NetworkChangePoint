# Nonparametric change point detection in growing networks

This repository contains code to compute the change point estimator developed in [Fluctuation bounds for continuous time branching processes and nonparametric change point detection in growing networks](https://arxiv.org/abs/1808.02439) and reproduces Figure 1 for this paper. `Using the code.ipynb` shows how to run the code on your own data. The network model is described briefly below, but see the paper for more details.

## Network model
Consider the following change point detection problem where a tree is grown to size n. The parameters of the problem are 0 < gamma < 1 (the change point), and f, g which are two, unknown attachment functions. The tree starts with say one node (v_1). At each time step a new vertex is added to the tree by selecting an existing vertex to attach to (where the edge points from the selected parent vertex to the new vertex). Let out-deg(v_i) denote the out-degree of vertex v_i. For vertices v_k, 2 <= k <= n*gamma, arriving in the first phase of the process, vertex v_k connects to existing vertex v_i (1 <= i < k) with probability proportional to f(out-deg(v_i)). For vertices v_k, n gamma < k <= n, arriving in the second phase of the process, vertex v_k connects to existing vertex v_i with probability proportional to g(out-deg(v_i)). 

To summarize, the tree evolves according to (unknown) function f, for the first n gamma steps then according to (unknown) function g for the remaining steps. The goal is then to estimate the unknown change point gamma.

When there is no change point, this model referred to as a *nonuniform random recursive tree*. A famous special case of this model is [preferential attachment](https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model). 


## Figure 1
First run
```
sh run_simulation
```
to simulate the trees and compute d_n(m) (the results are saved in `data/`). Then run `Figure 1.ipynb` to create the figures. The first part is a bit slow and I suggest parallelizing the simulation across a bunch of cluster nodes. 

