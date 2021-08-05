# Lab-tools

This repository includes known working versions of tools commonly used by the lab.

Known issues:

- Job submission scripts may need to be edited to accommodate the new cluster.
- Both the "default" RSFF2 script (that most people have been using) and Jovan's edited version have been included. These two files will need to be merged eventually.

Software links:

- GROMACS 2018.8: https://manual.gromacs.org/2018.8/download.html (Most recent release of the 2018 series; we are not updating to 2019+ due to removed support for NVML and group cutoff scheme, which we still use in vacuum simulations)
- PLUMED 2.5.4: https://github.com/plumed/plumed2/releases/tag/v2.5.4 (PLUMED has been updated since then, and newer versions may be worth looking into)
