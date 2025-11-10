### BHAC simulation with data from the Einstein Toolkit (ET)

The part that reads in the data is done by the file 
``` mod_oneblock.t```
The subroutine `read_oneblock` is the one that reads in the data points and assignes them as either a staggered variable or non-staggered. 

#### Our set-up
We used this code to evolve a jet launched from a Hyper-massive neutron star (HMNS) remnant from a binary neutron star (BNS) merger. The original simulation was run with the ET, where the HMNS collapsed after $\sim 30ms$ (details can be found in the original paper [Mösta et al 2020](https://ui.adsabs.harvard.edu/abs/2020ApJ...901L..37M/abstract)). We used this set-up to expand the jet to further distances which was not possible with the ET. 

We created a new simulation box in BHAC, where the top layer of the ET simulation was placed at the bottom of the BHAC box, and the rest was filled by some medium. We chose a Frozen-slab approach, where in every timestep the ET cells are re-interpolated and reset to their original values. This constantly "re-injects" material in the jet.
