#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "dock.h"
#include "load.h"
#include "util.h"
#include "rmsd.h"


using namespace std;


int main(int argc, char *argv[]) {
    string lig1_path;
    string lig2_path;
    string prt1_path;
    string prt2_path;

    int cms_flag = 0;
    int rmsd_flag = 0;
    int fraction_flag = 0;


    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",    no_argument, 0, 'h'},
            {"cms", no_argument, 0, 'c'},
            {"rmsd", no_argument, 0, 'r'},
            {"fraction", no_argument, 0, 'f'},
            {"lig1",    required_argument, 0, 'a'},
            {"lig2",    required_argument, 0, 'b'},
            {"prt1",    required_argument, 0, 'd'},
            {"prt2",    required_argument, 0, 'e'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hcrf",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
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

        case 'c':
            cms_flag = 1;
            break;

        case 'r':
            rmsd_flag = 1;
            break;

        case 'f':
            fraction_flag = 1;
            break;

        case 'a':
            lig1_path = optarg;
            break;

        case 'b':
            lig2_path = optarg;
            break;

        case 'd':
            prt1_path = optarg;
            break;

        case 'e':
            prt2_path = optarg;
            break;

        default:
            abort ();
        }
    }

    if (rmsd_flag) {
        if (lig1_path.empty() || lig2_path.empty()) {
            cout << "Wrong arguments!" << endl;
            usage();
            exit (EXIT_FAILURE);
        }

        Ligand0* lig1 = loadLigandSdf(lig1_path);
        Ligand0* lig2 = loadLigandSdf(lig2_path);

        float rmsd = calcRmsd(lig1, lig2);
        printf("rmsd value:\t%f\n", rmsd);

        free(lig1);
        free(lig2);
    }

    if (cms_flag || fraction_flag) {

        if (lig1_path.empty()
                || lig2_path.empty()
                || prt1_path.empty()
                || prt2_path.empty()) {
            cout << "Wrong arguments!" << endl;
            usage();
            exit (EXIT_FAILURE);
        }

        string para_path = getenv("GEAUX_FF");
        EnePara0 *enepara = loadPmf(para_path);

        Protein0* prt1 = loadProteinPdb(prt1_path);
        Protein0* prt2 = loadProteinPdb(prt2_path);
        Ligand0* lig1 = loadLigandSdf(lig1_path);
        Ligand0* lig2 = loadLigandSdf(lig2_path);

        ContactScore cnt = calculateContactScore(lig1, prt1, lig2, prt2, enepara);

        if (cms_flag) {
            printf("cms value:\t%f\n", cnt.cms);
        }

        if (fraction_flag) {
            printf("fraction value:\t%f\n", cnt.frac);
        }

        free(prt1);
        free(prt2);
        free(lig1);
        free(lig2);
        free(enepara);
    }

    return 0;
}
