#include <iostream>
#include <assert.h>
#include <stdio.h>

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

#include "dock.h"
#include "load.h"
#include "util.h"

using namespace std;

static
void printLigand(const Ligand0 *lig)
{
  assert(lig->lna != 0);
  int N = lig->lna;
  int i;
  for (i = 0; i < N; i++) {
    const LigCoord *coord = &lig->coord_orig;
    printf("%.3f\t%.3f\t%.3f\t%d\n",
           coord->x[i], coord->y[i], coord->z[i], lig->t[i]);
  }
}

static
void printProtein(const Protein0 *prt)
{
  assert(prt->pnp != 0);
  int N = prt->pnp;
  int i;
  for (i = 0; i < N; i++) {
    printf("%.3f\t%.3f\t%.3f\t%d\n",
           prt->x[i], prt->y[i], prt->z[i], prt->t[i]);
  }

}

TEST (load, Ligand)
{
  // Ligand0* lig = new Ligand0[1];

  string sdf_path = "../data/1a07C1.sdf";
  Ligand0* lig = loadLigandSdf(sdf_path);
  // printLigand(lig);

  free(lig);
}

TEST (load, Protein)
{
  Protein0 *prt = loadProteinPdb("../data/1a07C.pdb");
  EXPECT_EQ(prt->t[4], 20);
  EXPECT_EQ(prt->t[10], 5);
  EXPECT_EQ(prt->t[19], 27);
  // printProtein(prt);

  free(prt);
}

TEST (load, PMF)
{
  string para_path = "../data/paras";
  EnePara0 *enepara = loadPmf(para_path);

  free(enepara);
}

TEST (cms, usage)
{
  usage();
}
