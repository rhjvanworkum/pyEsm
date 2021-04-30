# pyEsm

pyEsm is a personal attempt at building a quantum chemistry software package in python.
The focus for this project was therefore more so about learning and not so much about building the fastest & most optimized quantum chemistry package.

#### Installation
Installation is recommended through Conda:
- Git clone the pyEsm project
- Inside pyEsm, type: 'conda env create -f environment.yml'
- next: 'conda activate pyEsmFinal'
- run 'pytest' in the terminal, to see if the environment creation was succesfull

#### Molecules
Energy calculations can be performed on Molecule objects, which can be initiated using a SMILES string. 

#### The following methods are currently implemented:
- (HF) Restricted Hartree Fock
- (MP2) Moller-Presset 2nd order perturbation theory
- (CCSD) Coupled Cluster single-double excitations
- (DFT) Density Functional theory using a LSDA functional

#### Basis-sets
The basis-set data is fetched via the basis set exchange API: https://www.basissetexchange.org/.
Most available basis-sets that can be found there will also work in pyEsm.

#### Evaluation of Molecular Gaussian Integrals
The McMurchie-Davidson scheme was used for evaluating the gaussian integrals,
largely based on Chapter 9 of the book 'Molecular Electronic-Structure Theory' by 
Helgaker, Jorgensen & Olsen.

##### For Further information see the tests folder