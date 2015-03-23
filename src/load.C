#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "load.h"
#include "data.h"

static 
vector < vector < string > > 
readLigandSections(const string sdf_path) 
{
  vector < vector < string > > sections;
  string line;
  ifstream file(sdf_path);

  if (!file.is_open()) {
    cout << "Error opening file " << sdf_path << endl;
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


void
loadLigandSdf(const string sdf_path, Ligand0 *mylig)
{

  vector < vector < string > > sects = readLigandSections(sdf_path);
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
}
