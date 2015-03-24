#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "dock.h"
#include "load.h"
#include "util.h"


using namespace std;

static void
usage()
{
  printf("usage:\n");
  cout << "cms --lig1 <first ligand> --prt1 <first protein>";
  cout << " --lig2 <second ligand> --prt2 <second protein>" << endl;
}

int main(int argc, char *argv[])
{
  string lig1_path;
  string lig2_path;
  string prt1_path;
  string prt2_path;
  
  int c;

  while (1)
    {
      static struct option long_options[] =
        {
          {"help",    no_argument, 0, 'h'},
          {"lig1",    required_argument, 0, 'a'},
          {"lig2",    required_argument, 0, 'b'},
          {"prt1",    required_argument, 0, 'c'},
          {"prt2",    required_argument, 0, 'd'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "a:b:c:d:",
                       long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'h':
          usage();
          exit (EXIT_FAILURE);
          break;

        case 'a':
          lig1_path = optarg;
          break;

        case 'b':
          lig2_path = optarg;
          break;

        case 'c':
          prt1_path = optarg;
          break;

        case 'd':
          prt2_path = optarg;
          break;

       default:
          abort ();
        }
    }

  if (lig1_path.empty()
      || lig2_path.empty()
      || prt1_path.empty()
      || prt2_path.empty()) {
    cout << "not enough arguments!" << endl;
    usage();
    exit (EXIT_FAILURE);
  }


  EnePara0 *enepara = new EnePara0[1];

  Protein0 *prt1 = new Protein0[1];
  Ligand0 *lig1 = new Ligand0[1];
  Protein0 *prt2 = new Protein0[1];
  Ligand0 *lig2 = new Ligand0[1];


  loadProteinPdb(prt1_path, prt1);
  loadLigandSdf(lig1_path, lig1);
  loadProteinPdb(prt2_path, prt2);
  loadLigandSdf(lig2_path, lig2);

  string para_path = "../data/paras";
  loadPmf(para_path, enepara);
  
  float cms = calculateContactModeScore(lig1, prt1, lig2, prt2, enepara);
  printf("cms value:\t%f\n", cms);

  delete[] prt1;
  delete[] lig1;
  delete[] lig2;
  delete[] prt2;
  delete[] enepara;

  return 0;
}
