#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>

#include "load.h"
#include "data.h"

/* split the sdf file based on $$$$ */
static
vector < vector < string > >
readLigandSdfSections(const string sdf_path) {
    vector < vector < string > > sections;
    string line;
    ifstream file(sdf_path.c_str());

    if (!file.is_open()) {
        cout << "Error opening file " << sdf_path << endl;
        return sections;
    }

    vector < string > one_section;
    while(getline(file, line)) {
        one_section.push_back(line);
        if (line.compare("$$$$") == 0) {
            sections.push_back(one_section);
            one_section.clear();
        }
    }

    return sections;
}

/* split the pdb file based on ENDMDL */
static
vector < vector < string > >
readProteinPdbSections(const string pdb_path) {
    vector < vector < string > > sections;
    string line;
    ifstream file(pdb_path.c_str());

    if (!file.is_open()) {
        cout << "Error opening file " << pdb_path << endl;
        return sections;
    }

    vector < string > one_section;
    while(getline(file, line)) {
        one_section.push_back(line);
        if (line.compare("ENDMDL") == 0 ||
            line.compare("TER") == 0) {
            sections.push_back(one_section);
            one_section.clear();
        }
    }

    return sections;
}


static
void
pushProteinPoint (Protein0 * prt, int r, int n, int t, int d, int c) {
    int i = prt->pnp;

    prt->r[i] = r;
    prt->n[i] = n;
    prt->t[i] = t;
    prt->d[i] = d;
    prt->c[i] = c;

    prt->pnp += 1;
}


static
void
pushLigandPoint (Ligand0 * lig, int n, int t, float c) {
    lig->n[n] = n;
    lig->t[n] = t;
    lig->c[n] = c;
}

Ligand0*
loadLigandSdf(const string sdf_path) {
    vector < vector < string > > sects = readLigandSdfSections(sdf_path);
    assert(sects.size() > 0);


    int num_lig = 1;
    Ligand0* mylig = (Ligand0*) calloc(num_lig, sizeof(Ligand0));
    if (NULL == mylig) {
        cout << "Error: Not enough memory in loading the ligand!" << endl;
        exit(EXIT_FAILURE);
    }

    // Ligand0* mylig = new Ligand0[1];

    vector < string > sect = sects[0];
    int total_lines = sect.size();
    int lnum = 0;
    string lines[MAXSDF];

    vector<string>::iterator iter;
    for (iter = sect.begin(); iter != sect.end(); ++iter) {
        lines[lnum++] = (*iter);
    }

    mylig->lna = atoi(lines[3].substr(0, 3).c_str());
    // mylig->lnb = atoi(lines[3].substr(3, 3).c_str());

    for (lnum = 4; lnum < 4 + mylig->lna; ++lnum) {
        string line = lines[lnum];
        float x = atof(line.substr(0, 10).c_str());
        float y = atof(line.substr(10, 10).c_str());
        float z = atof(line.substr(20, 10).c_str());

        mylig->coord_orig.x[lnum - 4] = x;
        mylig->coord_orig.y[lnum - 4] = y;
        mylig->coord_orig.z[lnum - 4] = z;
    }

    float tmp2[MAXLIG];	// for OB_ATOMIC_CHARGES
    int tmp3[MAXLIG];	// for atom types

    for (lnum = 4 + mylig->lna + mylig->lnb; lnum < total_lines; lnum++) {
        if (lines[lnum].find("OB_ATOM_TYPES") != string::npos) {
            int tmp4 = 0;

            istringstream tmp5(lines[lnum + 1]);

            while (tmp5) {
                std::string tmp6;

                tmp5 >> tmp6;

                if (tmp6.length() > 0)
                    tmp3[tmp4++] = getLigCode(tmp6);
            }
        }

        else if (lines[lnum].find("OB_ATOMIC_CHARGES") != string::npos) {
            int tmp4 = 0;

            istringstream tmp5(lines[lnum + 1]);

            while (tmp5) {
                std::string tmp6;

                tmp5 >> tmp6;

                if (tmp6.length() > 0)
                    tmp2[tmp4++] = atof(tmp6.c_str());
            }
        }
    }

    /* write the properties of each atom */
    for (int i1 = 4; i1 < mylig->lna + 4; i1++) {
        pushLigandPoint(mylig, i1 - 4, tmp3[i1 - 4], tmp2[i1 - 4]);
    }

    return mylig;
}


Protein0*
loadProteinPdb(const string pdb_path) {
    vector < vector < string > > sects = readProteinPdbSections(pdb_path);
    assert(sects.size() > 0);

    int num_prt = 1;
    Protein0* myprt = (Protein0*) calloc(num_prt, sizeof(Protein0));
    if (NULL == myprt) {
        cout << "Error: Not enough memory in loading the protein!" << endl;
        exit(EXIT_FAILURE);
    }

    vector < string > sect = sects[0];

    string protein_seq1;

    int pnp = 0, pnr = 0;
    myprt->pnr = 0;

    vector < string > ::iterator iter;
    for (iter = sect.begin(); iter != sect.end(); ++iter) {
        string line = *iter;
        if (line.size() > 53) {
            if (line.substr(0, 6) == "ATOM  ") {

                std::string atom1 = line.substr(12, 4);	// 3rd column
                int residue1 = getResCode(line.substr(17, 3));	// 4th column, residue name
                int residue2 = atoi(line.substr(22, 4).c_str());	// 5th column, residue serial number

                if (atom1 == " CA ") {

                    protein_seq1.append(three2oneS(line.substr(17, 3)));

                    myprt->seq3[pnr] = residue2;	// seq3 contains residue serial number

                    pushProteinPoint(myprt, pnr, pnp++, 0, residue1, 0);

                    if (pnr > 1) {
                        pushProteinPoint(myprt, pnr - 1, pnp++, 1, residue1, 1);
                    }

                    switch (residue1) {
                    case 0:
                        pushProteinPoint(myprt, pnr, pnp++, 2, residue1, 2);
                        break;
                    case 1:
                        pushProteinPoint(myprt, pnr, pnp++, 21, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 22, residue1, 3);
                        break;
                    case 2:
                        pushProteinPoint(myprt, pnr, pnp++, 3, residue1, 2);
                        break;
                    case 3:
                        pushProteinPoint(myprt, pnr, pnp++, 18, residue1, 2);
                        break;
                    case 4:
                        pushProteinPoint(myprt, pnr, pnp++, 14, residue1, 2);
                        break;
                    case 5:
                        pushProteinPoint(myprt, pnr, pnp++, 11, residue1, 2);
                        break;
                    case 6:
                        pushProteinPoint(myprt, pnr, pnp++, 5, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 6, residue1, 3);
                        break;
                    case 7:
                        pushProteinPoint(myprt, pnr, pnp++, 9, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 10, residue1, 3);
                        break;
                    case 8:
                        pushProteinPoint(myprt, pnr, pnp++, 26, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 27, residue1, 3);
                        break;
                    case 9:
                        pushProteinPoint(myprt, pnr, pnp++, 24, residue1, 2);
                        break;
                    case 10:
                        pushProteinPoint(myprt, pnr, pnp++, 15, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 16, residue1, 3);
                        break;
                    case 11:
                        pushProteinPoint(myprt, pnr, pnp++, 4, residue1, 2);
                        break;
                    case 13:
                        pushProteinPoint(myprt, pnr, pnp++, 7, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 8, residue1, 3);
                        break;
                    case 14:
                        pushProteinPoint(myprt, pnr, pnp++, 12, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 13, residue1, 3);
                        break;
                    case 15:
                        pushProteinPoint(myprt, pnr, pnp++, 19, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 20, residue1, 3);
                        break;
                    case 16:
                        pushProteinPoint(myprt, pnr, pnp++, 23, residue1, 2);
                        break;
                    case 17:
                        pushProteinPoint(myprt, pnr, pnp++, 17, residue1, 2);
                        break;
                    case 18:
                        pushProteinPoint(myprt, pnr, pnp++, 28, residue1, 2);
                        pushProteinPoint(myprt, pnr, pnp++, 29, residue1, 3);
                        break;
                    case 19:
                        pushProteinPoint(myprt, pnr, pnp++, 25, residue1, 2);
                        break;
                    }

                    (myprt->pnr)++;

                    pnr++;
                }
            }
        }
    }

    char protein_seq2[MAXPRO];

    strcpy(protein_seq2, protein_seq1.c_str());
    // cout << myprt->pnp << endl;

    // iterating effective points
    for (int p2_i = 0; p2_i < (myprt->pnp); p2_i++) {
        int residue3 = myprt->seq3[myprt->r[p2_i]];
        int residue5 = myprt->r[p2_i];
        int point1 = myprt->t[p2_i];

        float tx1 = 0.0;	// temperary x coord
        float ty1 = 0.0;	// temperary y coord
        float tz1 = 0.0;	// temperary z coord
        float tn1 = 0.0;	// ???, what is its use?

        for (iter = sect.begin(); iter != sect.end(); ++iter) {
            string line = *iter;
            if (line.size() > 53) {
                if (line.substr(0, 6) == "ATOM  ") {
                    int residue4 = atoi(line.substr(22, 4).c_str());

                    std::string atom2 = line.substr(12, 4);

                    float tx2 = atof(line.substr(30, 8).c_str());
                    float ty2 = atof(line.substr(38, 8).c_str());
                    float tz2 = atof(line.substr(46, 8).c_str());

                    if (residue4 == residue3) {
                        if (point1 == 0 && atom2 == " CA ") {
                            tx1 += tx2;
                            ty1 += ty2;
                            tz1 += tz2;

                            tn1 += 1.0;
                        } else if (point1 == 2
                                   || point1 == 3
                                   || point1 == 4
                                   || point1 == 11
                                   || point1 == 14
                                   || point1 == 17
                                   || point1 == 18
                                   || point1 == 23 || point1 == 24 || point1 == 25) {
                            if (atom2 != " N  "
                                    && atom2 != " CA " && atom2 != " C  " && atom2 != " O  ") {
                                tx1 += tx2;
                                ty1 += ty2;
                                tz1 += tz2;

                                tn1 += 1.0;
                            }
                        } else if (point1 == 5
                                   || point1 == 7
                                   || point1 == 9
                                   || point1 == 12
                                   || point1 == 15
                                   || point1 == 19
                                   || point1 == 21 || point1 == 26 || point1 == 28) {
                            if (atom2 != " N  "
                                    && atom2 != " CA "
                                    && atom2 != " C  "
                                    && atom2 != " O  " && atom2 != " CB " && atom2 != " CG ") {
                                tx1 += tx2;
                                ty1 += ty2;
                                tz1 += tz2;

                                tn1 += 1.0;
                            }
                        } else if (point1 == 6
                                   || point1 == 8
                                   || point1 == 10
                                   || point1 == 13
                                   || point1 == 16
                                   || point1 == 20
                                   || point1 == 22 || point1 == 27 || point1 == 29) {
                            if (atom2 == " CB " || atom2 == " CG ") {
                                tx1 += tx2;
                                ty1 += ty2;
                                tz1 += tz2;

                                tn1 += 1.0;
                            }
                        }
                    }

                    if (point1 == 1 && residue5 > 0 && residue5 < pnr - 1) {
                        int residue6 = myprt->seq3[myprt->r[p2_i] - 1];

                        if ((residue4 == residue3 && atom2 == " N  ")
                                || (residue4 == residue6 && atom2 == " C  ")
                                || (residue4 == residue6 && atom2 == " O  ")) {
                            tx1 += tx2;
                            ty1 += ty2;
                            tz1 += tz2;

                            tn1 += 1.0;
                        }
                    }
                }

            } else if (line.substr(0, 6) == "ENDMDL") {
                if (tn1 > 0.0) {
                    tx1 /= tn1;
                    ty1 /= tn1;
                    tz1 /= tn1;

                    myprt->x[p2_i] = tx1;	// set the x,y,z coords
                    myprt->y[p2_i] = ty1;
                    myprt->z[p2_i] = tz1;
                }
                tx1 = 0.0;
                ty1 = 0.0;
                tz1 = 0.0;
                tn1 = 0.0;

            }
        }
    }

    return myprt;
}

EnePara0*
loadPmf(const string para_path) {

    string line1;

    ifstream d1_file(para_path.c_str());

    if (!d1_file.is_open ()) {
        cout << "cannot open energy parameter file" << endl;
        cout << "Cannot open " << para_path << endl;
    }

    EnePara0* enepara = (EnePara0*) calloc(1, sizeof(EnePara0));
    if (NULL == enepara) {
        cout << "Not enough memory for EnePara in loading " << para_path << endl;
        exit(EXIT_FAILURE);
    }

    while (getline (d1_file, line1))
        if (line1.length () > 3) {
            if (line1.substr(0, 3) == "PMF") {
                std::string dat1[5];

                int dat2 = 0;

                istringstream dat3 (line1);

                while (dat3)
                    dat3 >> dat1[dat2++];

                enepara->pmf[getPntCode (dat1[1])][getLigCode (dat1[2])][0] =
                    atof (dat1[3].c_str ());
                enepara->pmf[getPntCode (dat1[1])][getLigCode (dat1[2])][1] =
                    atof (dat1[4].c_str ());
            }
        }

    return enepara;
}
