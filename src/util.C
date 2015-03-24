#include <stdlib.h>
#include <assert.h>
#include <cmath>

#include <iostream>
using namespace std;


#include "util.h"

/*
  construct the contact mode matrix based on the protein-ligand
  interface;
  return a pointer to the matrix with memory allocated
 */
void
initContactMatrix(int *ref_matrix, const Ligand0 *mylig, const Protein0 *myprt,
                  const EnePara0 *enepara)
{
  int lna = mylig->lna;
  int pnp = myprt->pnp;

  for (int l = 0; l < lna; l++) {
    const int lig_t = mylig->t[l];


    for (int p = 0; p < pnp; p++) {
      const int prt_t = myprt->t[p];

      const float dx = mylig->coord_orig.x[l] - myprt->x[p];
      const float dy = mylig->coord_orig.y[l] - myprt->y[p];
      const float dz = mylig->coord_orig.z[l] - myprt->z[p];
      const float dst = sqrtf (dx * dx + dy * dy + dz * dz);

      const float pmf0 = enepara->pmf[prt_t][lig_t][0];
      ref_matrix[l * pnp + p] = (dst <= pmf0);
    }
  }
}

/*
  compare two contact matrices and calcuate the cms value;
  will free the memory in ref1 and ref2.
 */
float
compareContacts(const int *ref1, const int *ref2, const int lna, const int pnp)
{
  int tp = 0;
  int fn = 0;
  int fp = 0;
  int tn = 0;

  for (int l = 0; l < lna; l++) {
    for (int p = 0; p < pnp; p++) {
      const int ref_val1 = ref1[l * pnp + p];
      const int ref_val2 = ref2[l * pnp + p];

      tp += (ref_val1 == 1 && ref_val2 == 1);
      fn += (ref_val1 == 1 && ref_val2 == 0);
      fp += (ref_val1 == 0 && ref_val2 == 1);
      tn += (ref_val1 == 0 && ref_val2 == 0);
    }
  }

  double d_tp = (double) tp;
  double d_fn = (double) fn ;
  double d_fp = (double) fp ;
  double d_tn = (double) tn ; 


  double cms = INVALID_CMS;

  double tmp = (d_tp + d_fp) * (d_tp + d_fn) *
    (d_tn + d_fp) * (d_tn + d_fn);
  
  if (tmp != 0.)
    cms = (d_tp * d_tn - d_fp * d_fn) / sqrtf(tmp);

  return (float) cms;
}

float
calculateContactModeScore(const Ligand0 *mylig1, const Protein0 *myprt1,
                          const Ligand0 *mylig2, const Protein0 *myprt2,
                          const EnePara0 *enepara)
{
  int lna = mylig1->lna;
  int pnp = myprt1->pnp;
  int total = lna * pnp;
  
  int *ref1 = new int[total];
  int *ref2 = new int[total];

  initContactMatrix(ref1, mylig1, myprt1, enepara);
  initContactMatrix(ref2, mylig2, myprt2, enepara);

  float cms = compareContacts(ref1, ref2, mylig1->lna, myprt1->pnp);

  delete[] ref1;
  delete[] ref2;

  return cms;
}
