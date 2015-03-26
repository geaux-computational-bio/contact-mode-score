#include <assert.h>
#include <math.h>

#include "rmsd.h"


float
calcRmsd(const Ligand0* const lig1, const Ligand0* const lig2)
{
  assert(lig1->lna == lig2->lna);
  assert(lig1->lna > 0);

  const int lna = lig1->lna;
  const LigCoord* const coord1 = &lig1->coord_orig;
  const LigCoord* const coord2 = &lig2->coord_orig;

  int i;
  float distance_squres = 0.0;
  for (i = 0; i < lna; i++) {
    const float d_x = coord1->x[i] - coord2->x[i];
    const float d_y = coord1->y[i] - coord2->y[i];
    const float d_z = coord1->z[i] - coord2->z[i];
    distance_squres += d_x * d_x + d_y * d_y + d_z * d_z;
  }

  float rmsd = sqrtf(distance_squres / (float) lna);
  return rmsd;
}
