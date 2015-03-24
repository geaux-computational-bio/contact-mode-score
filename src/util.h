#ifndef UTIL_H
#define UTIL_H

#include "dock.h"

#define INVALID_CMS 100.0

void initContactMatrix(int *ref_matrix, const Ligand0 *mylig,
                       const Protein0 *myprt, const EnePara0 *enepara);

float compareContacts(const int *ref1, const int *ref2, const int lna, const int pnp);

float calculateContactModeScore(const Ligand0 *mylig1, const Protein0 *myprt1,
                                const Ligand0 *mylig2, const Protein0 *myprt2,
                                const EnePara0 *enepara);

#endif /* UTIL_H */
