#!/usr/bin/env bash

bin=../src/cms

echo "between different conformations"
$bin -f -r -c --lig1 ./1a07C1.sdf --lig2 ./1a07C1_.sdf --prt1 ./1a07C.pdb --prt2 ./1a07C_.pdb

echo 

echo "for identical complexes"
$bin -frc --lig1 ./1a07C1.sdf --lig2 ./1a07C1.sdf --prt1 ./1a07C.pdb --prt2 ./1a07C.pdb
