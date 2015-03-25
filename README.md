# cms

Measure of protein-ligand complex similarity based on contacts

## Installation

    $ cd src
    $ make cms -f makefile

## Usage

    cms --lig1 <first ligand> --prt1 <first protein> --lig2 <second ligand> --prt2 <second protein>

## Examples
    $ cd data
    $ sh cms.sh

## Dependencies

    The codes have been tested on g++ 4.2.1, but any compilers supporting c++ standard libraries should work.


## To-Do

Currently only sdf format for ligand and pdb for protein are supported.
Improvements include:

- To calculate RMSD

- To calculate the fraction of non-specific contacts