This repository contains code to compute the change point estimator developed in []() and reproduces Figure 1 for this paper. `Using the code.ipynb` shows how to run the code on your own data.


# Figure 1

First grow 100 change point trees and compute d_n(m) for each of them. Warning: this is a bit slow for large values of n. I suggest parallelizing the simulation across a bunch of cluster nodes. 
```
sh run_simulation.
```
Next create the plots by running the code in `Paper figure 1.ipynb`