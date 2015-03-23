#ifndef DOCK_H
#define DOCK_H

#include "size.h"

struct LigCoord
{
  float x[MAXLIG];		// Ligand x coords
  float y[MAXLIG];		// Ligand y coords
  float z[MAXLIG];		// Ligand z coords
  float center[3];		// ligand geometric center
};

struct Ligand0
{
  LigCoord coord_orig;           //                                      used

  int t[MAXLIG];		// atom type                            used
  float c[MAXLIG];		// atom charge                          used
  int n[MAXLIG];		// atom number                          used

  int lna;			// number of ligand atoms               used
  int lnb;			// number of ligand bonds               NOT USED

  float pocket_center[3];	// pocket center                        used
                                // should belong to "Protein" structure
};

#endif /* DOCK_H */
