#ifndef LOAD_H
#define LOAD_H

#include <string>

#include "dock.h"

using namespace std;

Ligand0* loadLigandSdf(const string sdf_path);

Protein0* loadProteinPdb(const string pdb_path);

EnePara0* loadPmf(const string para_path);

#endif /* LOAD_H */
