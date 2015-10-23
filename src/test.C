#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <cmath>

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

#include "dock.h"
#include "load.h"
#include "util.h"

using namespace std;

static
void printLigand(const Ligand0 *lig) {
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
void printProtein(const Protein0 *prt) {
    assert(prt->pnp != 0);
    int N = prt->pnp;
    int i;
    for (i = 0; i < N; i++) {
        printf("%.3f\t%.3f\t%.3f\t%d\n",
               prt->x[i], prt->y[i], prt->z[i], prt->t[i]);
    }

}

TEST (load, Ligand) {
    // Ligand0* lig = new Ligand0[1];

    string sdf_path = "../data/1a07C1.sdf";
    Ligand0* lig = loadLigandSdf(sdf_path);
    // printLigand(lig);

    free(lig);
}

TEST (load, Protein) {
    Protein0 *prt = loadProteinPdb("../data/1a07C.pdb");
    EXPECT_EQ(prt->t[4], 20);
    EXPECT_EQ(prt->t[10], 5);
    EXPECT_EQ(prt->t[19], 27);
    // printProtein(prt);

    free(prt);
}

TEST (load, PMF) {
    string para_path = "../data/paras";
    EnePara0 *enepara = loadPmf(para_path);

    free(enepara);
}

TEST(calculate, MCC) {
  int tp = 0, fn = 0, fp = 0;
  int tn = 10;

  double d_tp = (double) tp;
  double d_fn = (double) fn ;
  double d_fp = (double) fp ;
  double d_tn = (double) tn ;

  double tmp = (d_tp + d_fp) * (d_tp + d_fn) *
    (d_tn + d_fp) * (d_tn + d_fn);

  double cms = 0.0;
  if (tmp != 0.)
    cms = (d_tp * d_tn - d_fp * d_fn) / sqrtf(tmp);

  EXPECT_EQ(cms, 0.0);
}


