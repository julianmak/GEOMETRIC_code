# GEOMETRIC code

Collection of code for GEOMETRIC implementations. For a brief description, there is a [ReadTheDocs](https://nemo-related.readthedocs.io/en/latest/) document. Some MITgcm versions will be put up here in due course too.

The code in this repository are unsupported and only my own attempts at backporting etc.; for the semi-officially supported ones, go to the [NEMO branch](https://forge.nemo-ocean.eu/nemo/nemo/-/merge_requests/202)

The nemo_paths and compiling option file is customised to my own laptop (where I installed libraries at slightly different places for various reasons to do with NEMO, and the compiler is fairly ancient).

* NOTE: I tend to only work with the latest version of the NEMO/MITgcm GEOMETRIC code, and it is unlikely I will actively be backporting code. The main body of the codes so far are largely the same, so those should serve as useful references if you want to backport code yourself.

[16 Oct 23]

There is a bug in the NEMO implementations of GEOMETRIC (see `ldfeke.F90`). Minimal fixes provided here.

The advection of parameterised eddy energy is missing a scale factor, so it's basically too small to do anything meaningful. Upon fixing that scale factor others show up: 

* there is an extra minus sign in the trend (which will crash the run)
* an override for the computed advective trend when linear free surface option is used (commented out for now)
* differences in the definition of the barotropic flow `[uv]n_adv` (in newer versions it seems to be depth-integrated, while older versions it is depth-averaged)

There are some minor issues that could be fixed (appropriate multiplication/divisions by TUV depth factors, using a better advection scheme such as TVD, possibly some off-by-one time-step mismatches), but the main one is that scale factor, which does make a big enough difference.

* NOTE: If you are here in relation to the Mak et al. (2022) GRL submission, with regards to full reproducibility, the folder you want is `nemo3.7-8666`, and it is the `ldfeke_full_eigen.F90` you want to use. This version is no longer used because the eigenvalue solver for computing the long Rossby phase speed is expensive, given there is a closed form approximation to the relevant eigenvalue obtained through WKB.
