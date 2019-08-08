# Quant_Chem

FORTRAN 90 code for Hatree-Fock method. Only s-type integration is included. That is only s-type basis sets can be used.

The binary `scf.out` can be built by the command
```
$ make
```

Instructions for the input file `ctrfile.dat`
```
natom   = 2 ! number of atoms
charge  = 0
mult    = 1 ! spin multiplicity
kmax    = 10000 ! max iteration for SCF
eps1    = 1d-20 ! convergence criterion for density
eps2    = 1d-16 ! convergence criterion for energy
fdiis   = .TRUE. ! flag for DIIS
```

Instructions for the configuration file `config.dat`
```
! Put atoms in sequence, natom atoms in all
! For each atom, following information is necessary
            1 ! atomic number
  0.0 0.0 0.0 ! coordinate in angstrom
```

Framework written by Xinzijian Liu. Gaussian integration written by Ning Zhang. DIIS written by Yuhang Yao.
