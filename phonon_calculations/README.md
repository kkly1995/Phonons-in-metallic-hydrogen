This directory contains the input files to do phonon calculations. Calculations were performed using [ALAMODE](https://github.com/ttadano/alamode) v1.1.0.

For both densities there are two input files: `force_constants.in` reads forces and displacements and from them calculates the force constants with symmetrized linear regression; `band_structure.in` reads these force constants and from them produces the corresponding phonon band structure.

The forces can be found in `DFSET` and `DFSET_with_frozen_phonons`. These forces are the result of RQMC calculations (see `force_calculations`) on randomly displaced configurations, with the latter including forces from a frozen phonon calculation. Note that because they are QMC forces, they include some statistical noise.

Both input files can easily be adjusted. For example, to change which set of forces to calculate force constants from, change the `DFSET` argument in `force_constants.in`. The program also can exclude certain samples in the set of forces with the use of the `SKIP` argument. To calculate the band structure for different paths through the BZ, change the numbers and symbols under `&kpoint` in `band_structure.in`.
