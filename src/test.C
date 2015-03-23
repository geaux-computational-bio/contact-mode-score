#include <iostream>
#include <assert.h>
#include <stdio.h>

#include "gtest/gtest.h"
#include "gtest/internal/gtest-internal.h"

#include "dock.h"
#include "load.h"

using namespace std;

static
void printLigand(const Ligand0* lig)
{
  assert(lig->lna != 0);
  int N = lig->lna;
  int i;
  for (i = 0; i < N; i++) {
    const LigCoord * coord = &lig->coord_orig;
    printf("%.3f %.3f %.3f\n",
           coord->x[i], coord->y[i], coord->z[i]);
  }

}

TEST (load, Ligand)
{
  Ligand0 *lig = new Ligand0[1];

  string sdf_path = "../data/1a07C1.sdf";
  loadLigandSdf(sdf_path, lig);
  printLigand(lig);

  delete[] lig;
}
