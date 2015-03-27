#ifndef DOCK_H
#define DOCK_H

#include <string>
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
  std::string a[MAXLIG];	// atom name                            NOT USED
};

struct Protein0
{
  float x[MAXPRO];		// residue x coord                      used
  float y[MAXPRO];		// residue y coord                      used
  float z[MAXPRO];		// residue z coord                      used
  int n[MAXPRO];		// effective point number               NOT USED
  int t[MAXPRO];		// effective point type                 used
  int c[MAXPRO];		// effective point class                used
  int d[MAXPRO];		// redidue code                         used

  int pnp;			// number of protein effective points   used
  int pnr;			// number of protein residues           NOT USED

  int r[MAXPRO];		// residue number                       replaced
  int seq3[MAXPRO];		// aa sequence numbering                replaced
};

struct EnePara0
{
  float vdw[MAXTP1][MAXTP2][2];	// L-J potential                        vdw[prt_t][lig_t][]
  float ele[MAXTP3];		// electrostatic potential              ele[prt_d + 30], ele[prt_t]
  float pmf[MAXTP1][MAXTP2][2];	// contact potential                    pmf[prt_t][lig_t][]
  float hpp[MAXTP4];		// protein hydrophobicity               hpp[prt_d]
  float hpl[MAXTP2][2];		// ligand hydrophobicity                hpl[lig_t][]
  float hdb[MAXTP1][MAXTP2][2];	// ligand hydrophobicity                hdb[prt_t][lig_t][]

  float lj[3];			// L-J params
  float el[2];			// electrostatic params
  float kde;			// kde bandwidth

  float w[MAXWEI];		// weights for energy terms
  float a_para[MAXWEI];         // the a parameter in normalization
  float b_para[MAXWEI];         // the b parameter in normalization
};

struct ContactScore
{
  float cms;
  float frac;
};

#endif /* DOCK_H */
