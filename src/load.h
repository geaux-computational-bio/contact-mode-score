#ifndef LOAD_H
#define LOAD_H

#include <string>

#include "dock.h"

using namespace std;

void loadLigandSdf(const string sdf_path, Ligand0 *lig);

void loadProteinPdb(const string pdb_path, Protein0 *prt);

void loadPmf(const string para_path, EnePara0* enepara);

#endif /* LOAD_H */
