# GEOMETRIC code

Collection of code for GEOMETRIC implementations. For a brief description, there is a [ReadTheDocs](https://nemo-related.readthedocs.io/en/latest/) document. Some MITgcm versions will be put up here in due course too.

The code in this repository are unsupported and only my own attempts at backporting etc.; for the semi-officially supported ones, go to the [NEMO branch](https://forge.ipsl.jussieu.fr/nemo/browser/NEMO/branches/NERC/dev_r4.0.6_GEOMETRIC)

The nemo_paths and compiling option file is customised to my own laptop (where I installed libraries at slightly different places for various reasons to do with NEMO).

* NOTE: If you are here in relation to the Mak et al. (2022) GRL submission, with regards to full reproducibility, the folder you want is `nemo3.7-8666`, and it is the `ldfeke_full_eigen.F90` you want to use. This version is no longer used because the eigenvalue solver for computing the long Rossby phase speed is expensive, given there is a closed form approximation to the relevant eigenvalue obtained through WKB.
