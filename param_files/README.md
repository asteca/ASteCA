The files:

- `params_input.dat`
- `semi_input.dat`

located inside this folder are examples of working input parameters files
and can be copied or moved to the root folder (where `asteca.py` is located)
to be used by the code.

The `/input/CLUSTER.DAT` file contains a synthetic open cluster generated
via the [MASSCLEAN][1] package with the following parameter values:

- z =
- log(age) =
- E(B-V) =
- (m-M)o =

and serves as an example cluster to be analyzed with ASteCA.

Both the `/param_files` folder and the `CLUSTER.DAT` file can be
safely deleted.

**Important**: a properly formatted `params_input.dat` file should always be
present in the root folder for ASteCA to work.

[1]: http://www.physics.uc.edu/~bogdan/massclean/