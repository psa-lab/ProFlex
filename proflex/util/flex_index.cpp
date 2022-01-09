/*******************************************************************************
 *  MSU ProFlex, formerly called FIRST, is a software developed to predict and *
 *  analyze protein flexibility. * This source file is a part of MSU ProFlex. *
 *                                                                              *
 *  Copyright (C) 1997 - 2008, Michigan State University. *
 *                                                                              *
 *  This program is free software; you can redistribute to academic users only,
 ** it and/or modify it under the terms of the GNU General Public License, *
 *  version 2, as published by the Free Software Foundation. *
 *                                                                              *
 *  This program is distributed in the hope that it will be useful, * but
 *WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the * GNU General
 *Public License for more details.                                *
 *                                                                              *
 *  You should have received a copy of the GNU General Public License * along
 *with this program; if not, write to the Free Software                 *
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA, *
 *  or see http://www.gnu.org/licenses/gpl.txt *
 *******************************************************************************/

/*-----------------------------------------------------------------------------
                                    Written by S.K. Namilikonda	05/2008

   This program reads in the "<prot>_allbonds" file and for each atom
   it picks the bond with the least flexibility index over all bonds
   stemming from it and assign its value to that atom. This is a coarse
   grained description of the flexibility index over the entire protein.
   This value (which ranges between -1 (min flex) to +1 (max flex))
   is then transformed into a new number and put in the B-value column
   of a PDB file ranging between 0 to 99, so that visualization tools
   such as Insight II or PyMol can color code the protein in accordance
   to the flexibility index.

   In addition, this program also outputs PyMol-ready scripts for visualizing
   the input protein based on flex_index as well as the rigid and flexible
   clusters. By default, the hydrophobic tethers are output to the flex*.pdb
   file but will be hidden in the PyMol viewer through a command that is set
   in the PyMol scripts output.

-----------------------------------------------------------------------------*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

#define MAX_FILENAME 80

int main(int argc, char *argv[]) {
  //-------------------------------------------------- declare variables
  ifstream input;
  ifstream input1;
  ifstream input2;
  ofstream output;
  ofstream output1;
  ofstream output2;
  ofstream output3; // crude getaround for glibc error
  char intro[10];
  int count, nHP = 0;
  const char error_open[] = "\n *** main(): unable to open file: ";
  //-------------------------------------------------- preset names for now

  /*
   * Note that only 20 colors are initialized
   * Use bright colors that match the coloring scheme in hbdilute plot
   * Use pastale colors for flex clusters
   *			R	G	B
   *			-	-	-
   *	red		1	0	0
   *	cyan		0	1	1
   *	green		0	1	0
   *	yellow		1	1	0
   *	purple		0.75	0	0.75
   *	magenta		1	0	1
   *	orange		1.0	0.5	0.0
   *
        wheat		0.99	0.82	0.65
        palegreen	0.65	0.9	0.65
        lightblue	0.75	0.75	1.0
        paleyellow	1.0	1.0	0.5
        lightpink	1.00	0.75	0.87
        palecyan	0.8	1.0	1.0
        lightorange	1.0	0.8	0.5
        bluewhite	0.85	0.85	1.00
        tv_orange	1.0	0.55	0.15
   */
  string RCcolors[] = {"",        "red",      "cyan",       "green",
                       "yellow",  "purple",   "purpleblue", "brightorange",
                       "tv_red",  "tv_green", "olive",      "deeppurple",
                       "blue",    "tv_blue",  "grey70",     "limon",
                       "hotpink", "violet",   "ruby",       "yelloworange",
                       "forest"};

  string FCcolors[] = {"",           "wheat",      "palegreen",   "lightblue",
                       "paleyellow", "lightpink",  "palecyan",    "lightorange",
                       "bluewhite",  "teal",       "limegreen",   "pink",
                       "slate",      "aquamarine", "lightorange", "sand",
                       "tv_orange",  "lime",       "marine",      "splitpea",
                       "salmon"};

  //-------------------------------------------------- input/output files

  char file_weight[MAX_FILENAME]; // allbonds file
  char file_chem[MAX_FILENAME];   // proflex dataset

  char file_flexibility[MAX_FILENAME]; // flex_index and RC analysis
                                       // output in b-value and atom #
                                       // columns respectively	--- SN

  char file_flex_pml[MAX_FILENAME]; // PyMol script for visualizing
                                    // the protein colored by
                                    // flexindex values ---	SN

  char file_RC_pml[MAX_FILENAME]; // PyMol script for visualizing
                                  // the protein colored by
                                  // RC values --- SN

  char file_FC_pml[MAX_FILENAME]; // PyMol script for visualizing
                                  // the protein colored by
                                  // FC values --- SN

  char ans, answer[23];
  char sso[6], ssf[6];
  std::string fend1 = "proflexdataset";
  std::string fend2 = ".pdb";
  std::string fend3 = ".pml"; // PyMol script ----	SN
  short int num_full, num_base, num_PROFLEX;
  int flagHP = -1;
  short int noscaling = 1;
  int max_atm_num, tmp, a1, a2, clust, i;
  float wt;
  char *line = new char[90];

  //--------------------------------------------- initialize file-name arrays to
  //NULL chars

  memset((void *)file_weight, '\0', MAX_FILENAME);
  memset((void *)file_chem, '\0', MAX_FILENAME);
  memset((void *)file_flexibility, '\0', MAX_FILENAME);
  memset((void *)file_flex_pml, '\0', MAX_FILENAME); // 2008:04
  memset((void *)file_RC_pml, '\0', MAX_FILENAME);   // 2008:04
  memset((void *)file_FC_pml, '\0', MAX_FILENAME);   // 2008:05

  //-------------------------------------------------------------------
  //cosmetics

  //  if (argc >= 2 && strcmp(argv[1], "-1") == 0) {
  // Don't scale the output to 0 to 99 if "-1" is given as the command line
  // argument
  //      noscaling = 1;
  //  }

  if (argc != 2) {
    cout << "\n\tUSAGE: flex_index <allbonds filename> \n" << endl;
    exit(-1);
  }

  /* OPEN INPUT FILES */

  /*
   * Read-in all_bonds file
   */
  strcpy(file_weight, argv[1]);
  num_full = strlen(file_weight);
  num_base = -1;

  /*
   * Output tether info to the flex*.pdb file by default
   */
  flagHP = 0;

  //---------------------  End Input file setup  -----------------------

  /*  Setup OUTPUT files */

  /*
   * Figure out the prefix of "allbonds" file
   */
  for (i = 0; i < MAX_FILENAME; i++) {
    file_chem[i] = file_weight[i];
    file_flexibility[i] = file_weight[i];
    file_flex_pml[i] = file_weight[i]; // 2008:04
    file_RC_pml[i] = file_weight[i];   // 2008:04
    file_FC_pml[i] = file_weight[i];   // 2008:04

    if (file_weight[i] == '_' && file_weight[i + 1] == 'a' &&
        file_weight[i + 2] == 'l') {
      num_base = i;
      break;
    }
  }

  if (num_base == -1)
    return 0;

  strcat(file_chem, fend1.c_str());

  num_PROFLEX = strlen(file_chem);

  file_flexibility[num_base + 1] = 'f';
  file_flexibility[num_base + 2] = 'l';
  file_flexibility[num_base + 3] = 'e';
  file_flexibility[num_base + 4] = 'x';
  file_flexibility[num_base + 5] = '_';

  file_flex_pml[num_base + 1] = 'f';
  file_flex_pml[num_base + 2] = 'l';
  file_flex_pml[num_base + 3] = 'e';
  file_flex_pml[num_base + 4] = 'x';
  file_flex_pml[num_base + 5] = '_';

  file_RC_pml[num_base + 1] = 'R';
  file_RC_pml[num_base + 2] = 'C';
  file_RC_pml[num_base + 3] = '_';

  file_FC_pml[num_base + 1] = 'F';
  file_FC_pml[num_base + 2] = 'C';
  file_FC_pml[num_base + 3] = '_';

  num_full = num_full - 4; // Add number (output-id)
  for (i = 1; i <= 4; i++, num_full++) {
    file_flexibility[num_base + 5 + i] = file_weight[num_full];
    file_flex_pml[num_base + 5 + i] = file_weight[num_full];
    file_RC_pml[num_base + 3 + i] = file_weight[num_full];
    file_FC_pml[num_base + 3 + i] = file_weight[num_full];
  }

  strcat(file_flexibility, fend2.c_str()); // Add extension ".pdb"
  strcat(file_flex_pml, fend3.c_str());    // Add extension ".pml"
  strcat(file_RC_pml, fend3.c_str());      // Add extension ".pml"
  strcat(file_FC_pml, fend3.c_str());      // Add extension ".pml"

  // --- Parse the input file and store first three columns as
  // --- the two atom numbers of each bond and its flexibility weight
  // --- To do the above, figure out the max_atm_# and allocate memory
  // --- to store the final flexibility index of each atom
  count = 0;
  max_atm_num = 0;

  input.open(file_weight, ios::in);
  if (input.fail()) {
    cerr << error_open << file_weight << "\n\n";
    exit(-1);
  }

  input >> a1 >> a2 >> wt >> clust;
  while (!input.eof()) {
    count++;

    if (a1 > max_atm_num)
      max_atm_num = a1;
    if (a2 > max_atm_num)
      max_atm_num = a2;

    input >> a1 >> a2 >> wt >> clust;
  }

  input1.open(file_weight, ios::in);
  if (input1.fail()) {
    cerr << error_open << file_weight << "\n\n";
    exit(-1);
  }

  // --- Now that we have the size of the input file,
  // --- allocate memory to store all the records
  int *atm1, *atm2, *clst, *label; // label stores atom-wise clst
  int min_label = 0, max_label = 0;
  float *abs_flex, *scaled_flex;
  atm1 = new int[count];
  atm2 = new int[count];
  clst = new int[count];
  label = new int[max_atm_num];
  abs_flex = new float[count];
  scaled_flex = new float[max_atm_num];

  //-------------------------------------------------------- read data from
  //allbonds file

  for (i = 0; i < count; i++) {
    input1 >> atm1[i] >> atm2[i] >> abs_flex[i] >> clst[i];
    if (i < max_atm_num) {
      scaled_flex[i] = 1.0; // initialize to max flex wt allowed (1.0)
      label[i] = count;     // initialize to max flex
    }
  }

  //------------------------------------------------- find min flexibility index
  for (i = 0; i < count; i++) {
    if (abs_flex[i] < scaled_flex[atm1[i] - 1]) {
      scaled_flex[atm1[i] - 1] = abs_flex[i];
    }
    if (abs_flex[i] < scaled_flex[atm2[i] - 1]) {
      scaled_flex[atm2[i] - 1] = abs_flex[i];
    }

    //------------------------------------------------- find min flexibility
    //label
    if ((clst[i] < 0 && label[atm1[i] - 1] >= 0) ||
        (label[atm1[i] - 1] == count)) {
      label[atm1[i] - 1] = clst[i];
    }
    if ((clst[i] < 0 && label[atm2[i] - 1] >= 0) ||
        (label[atm2[i] - 1] == count)) {
      label[atm2[i] - 1] = clst[i];
    }
    if (clst[i] < min_label)
      min_label = clst[i];
    if (clst[i] > max_label)
      max_label = clst[i];
  }
  //----------------------------------------------- Now, scale the weight
  //factors 						  from (-1,+1) to (0,99)
  if (!noscaling) // Added as Kevin's StoneHinge fixes
  {
    for (i = 0; i < max_atm_num; i++) {
      scaled_flex[i] = (scaled_flex[i] + 1) * 49.5;
    }
  }

  //------------------------  End of flex_index calculation
  //-----------------------------

  /*    RIGID and FLEXIBLE CLUSTER IDENTIFICATION AND COLORING    */

  /*
   * clst[] contains the cluster label for the various bonds in the allbonds
   * file So, we first map the bond labels onto the atoms such that a -ve label
   * implies that an atom belongs to a rigid cluster; +ve implies flexible
   * cluster; and ZERO for uncorrelated motion. This info is in label[].
   *
   * Just as the flex_index value is assigned, each atom is assigned to the
   * most rigid cluster that it forms a part of. Now, we renumber the atoms
   * such that each cluster is assigned to a range of atom numbers. The range
   * of labels is given by min_label:max_label; The atom re-numbering proceeds
   * in the order of labels => rigid clusters followed by flexible ones.
   *
   * Then, generate two PyMol scripts that create objects for each of the rigid
   * (and flexible) clusters and colors them with a different color.
   * Note that each cluster should at least have a size of 6 atoms! Further, the
   * cluster boundaries have to be scanned carefully to determine the biggest
   * RCs and FCs
   */

  int MAX_CLUST = max_label - min_label + 10;
  int *atmIndex, *clstBoundary;
  atmIndex = new int[max_atm_num];
  clstBoundary = new int[MAX_CLUST];
  char header_written_flag, ch;
  int begin, end, size, j;

  for (i = 0; i < MAX_CLUST; i++) {
    clstBoundary[i] = 0; // Initialize clstBndry
  }

  /*
   *	Scan the label[] list and re-index the atoms! Start with rgd clusters.
   *	Then, renumber FCs. Also, make a note of all the clusters' boundaries.
   */

  size = 0;      // atm-index temp holder
  j = min_label; // clust-index
  while (size < max_atm_num && j <= max_label) {
    for (i = 0; i < max_atm_num; i++) {
      if (label[i] == j) {
        size++;
        atmIndex[i] = size;
      }
    }
    clstBoundary[j - min_label] = size;
    j++; // Next clust-index
  }

  //-----------------------  End of RC & FC Identification --------------------

  /*
   *	Write the flex.pdb with the revised atom-ids and b-value column with
   *     the re-scaled flex_index values
   */

  input2.open(file_chem, ios::in);
  if (input2.fail()) {
    cerr << error_open << file_chem << "\n\n";
    exit(-1);
  }

  output.open(file_flex_pml, ios::out);
  output1.open(file_flexibility, ios::out);
  output2.open(file_RC_pml, ios::out);
  output3.open(file_FC_pml, ios::out);

  if (output1.fail()) {
    cerr << error_open << file_flexibility << "\n\n";
    exit(-1);
  }

  count = 0;
  while (!input2.eof() && !input2.fail()) {
    input2.getline(line, 90);

    if (line != NULL &&
        (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0)) {
      if (flagHP && !strncmp(line + 13, "X ", 2) &&
          !strncmp(line + 17, "XXX", 3)) {
        nHP++;
      } else {
        for (i = 0; i < 6; i++)
          output1 << *(line + i);
        output1 << setw(5) << setiosflags(ios::fixed | ios::right)
                << atmIndex[count];
        for (i = 11; i < 54; i++)
          output1 << *(line + i);
        output1 << "  1.00 " << setw(5) << setprecision(2)
                << setiosflags(ios::showpoint | ios::fixed | ios::right)
                << *(scaled_flex + count) << endl;
      }
      count++;
    }

    if (line != NULL && strncmp(line, "REMARK", 6) == 0 && count > 0)
      break;
  }

  if (!flagHP) {
    cout << "\n\n\tThe number of ATOMS = " << count << endl;
  } else if (flagHP) {
    cout << "\n\n\tThe number of ATOMS = " << count - nHP << endl;
    cout << "\n\tThe number of Tethers omitted = " << nHP / 3 << endl;
  }

  /* #################################################################### */
  /*									*/
  /*			Output PyMol Scripts				*/
  /*									*/
  /* #################################################################### */

  /*
   * Now that the flexindex based PDB file has been written out, generate
   * a pymol script to aid visualizing it.
   */

  //  output.open(file_flex_pml,ios::out);
  output << "delete FlexIndexObj" << endl;
  output << "delete HPHOB1" << endl;
  output << "load " << file_flexibility << ", FlexIndexObj" << endl << endl;
  // Load molecule (*flex*.pdb)

  output << "cmd.color( 's170', 'FlexIndexObj and b< 40.0')" << endl;
  output << "cmd.color( 's170', 'FlexIndexObj and b= 40.0')" << endl;
  output << "cmd.color( 's200', 'FlexIndexObj and b> 40.0 and b< 42.0')"
         << endl;
  output << "cmd.color( 's200', 'FlexIndexObj and b= 42.0')" << endl;
  output << "cmd.color('s225', 'FlexIndexObj and b> 42.0 and b < 44.0')"
         << endl;
  output << "cmd.color('s225', 'FlexIndexObj and b= 44.0')" << endl;
  output << "cmd.color('s250', 'FlexIndexObj and b> 44.0 and b < 46.0')"
         << endl;
  output << "cmd.color('s250', 'FlexIndexObj and b= 46.0')" << endl;
  output << "cmd.color('s280', 'FlexIndexObj and b> 46.0 and b < 48.0')"
         << endl;
  output << "cmd.color('s280', 'FlexIndexObj and b= 48.0')" << endl;
  output << "cmd.color('s310', 'FlexIndexObj and b> 48.0 and b < 49.0')"
         << endl;
  output << "cmd.color('s310', 'FlexIndexObj and b= 49.0')" << endl;

  output << "cmd.color('grey', 'FlexIndexObj and b> 49.0 and b < 50.0')"
         << endl;
  output << "cmd.color('grey', 'FlexIndexObj and b= 50.0')" << endl;

  output << "cmd.color('s690', 'FlexIndexObj and b> 50.0 and b < 52.0')"
         << endl;
  output << "cmd.color('s690', 'FlexIndexObj and b= 52.0')" << endl;
  output << "cmd.color('s730', 'FlexIndexObj and b> 52.0 and b < 54.0')"
         << endl;
  output << "cmd.color('s730', 'FlexIndexObj and b= 54.0')" << endl;
  output << "cmd.color('s760', 'FlexIndexObj and b> 54.0 and b < 56.0')"
         << endl;
  output << "cmd.color('s760', 'FlexIndexObj and b= 56.0')" << endl;
  output << "cmd.color('s790', 'FlexIndexObj and b> 56.0 and b < 58.0')"
         << endl;
  output << "cmd.color('s790', 'FlexIndexObj and b= 58.0')" << endl;
  output << "cmd.color('s820', 'FlexIndexObj and b> 58 and b < 60.0')" << endl
         << endl
         << endl;
  output << "cmd.color('s820', 'FlexIndexObj and b= 60.0')" << endl;
  output << "cmd.color('s850', 'FlexIndexObj and b> 60.0')" << endl
         << endl
         << endl;

  // Show HPHOB as spheres
  output << "select HPHOB1, FlexIndexObj and resn \"XXX\" " << endl;
  output << "show spheres, HPHOB1" << endl;
  output << "set sphere_scale=0.4" << endl;
  output << "deselect HPHOB1" << endl;

  // Hide Tethers
  // output<<"hide everything, HPHOB1"<<endl;

  /*
   *	Write the RC visualization script
   */

  //  output2.open(file_RC_pml,ios::out);
  output2 << "delete RigidClustObj" << endl;
  output2 << "delete RC*" << endl;
  output2 << "delete smallRC*" << endl;
  output2 << "delete FC*" << endl;
  output2 << "delete HPHOB2" << endl;
  output2 << "load " << file_flexibility << ", RigidClustObj" << endl << endl;
  // Load molecule (*flex*.pdb)

  // First process the rigid clusters => -1-min_label through 0
  // Note that min_label is negative (at least -1). So, we have to
  // first start with the cluster that has label = -1 to label = min_label

  for (i = (min_label * -1) - 1, j = 1; i >= 0; i--, j++) {
    if (clstBoundary[i] - clstBoundary[i - 1] > 6) {
      output2 << "\nselect RC" << j << ", RigidClustObj and id "
              << clstBoundary[i - 1] + 1 << "-" << clstBoundary[i];
      output2 << "\ncolor " << RCcolors[(j % 20) + 1] << ", RC" << j;

      //	output2<<"\npreset.ball_and_stick('RC"<<i<<"')"<<endl;
      // 	Just show the main-chain by default as a ribbon or "cartoon" in
      // PyMol

      output2 << "\ncartoon automatic, RC" << j;
      output2 << "\nshow cartoon, RC" << j;

      output2 << endl;
    }
  }

  /*
   * Combine all the clusters of size <= 6 intoa single cluster
   */
  char flag = 'F';
  for (i = (min_label * -1) - 1, j = 1; i >= 0; i--, j++) {
    if (clstBoundary[i] - clstBoundary[i - 1] <= 6) {
      if (flag == 'F') {
        // The rest of the atoms will be colored with a default color
        output2 << "\nselect smallRCs" << j / 50 << ", RigidClustObj and id "
                << clstBoundary[i - 1] + 1 << "-" << clstBoundary[i];
        flag = 'T';
      } else {
        output2 << "+" << clstBoundary[i - 1] + 1 << "-" << clstBoundary[i];
      }
    }

    if ((j + 1) % 50 == 0) {
      flag = 'F'; // To prevent string overflow of clust indices in PyMol
      output2 << "\ncolor brown, smallRCs" << j / 50;
      output2 << "\nshow sticks, smallRCs" << j / 50 << endl;
    }
  }

  output2 << endl;

  //	All the FCs including uncorrelated motion cluster are grouped under
  //	one cluster that is colored white
  output2 << "\n\nselect FC, RigidClustObj and id "
          << clstBoundary[(min_label * -1) - 1] + 1 << "-"
          << clstBoundary[(min_label * -1) + max_label];
  output2 << "\ncolor white, FC";
  output2 << "\nshow lines, FC";

  // Show HPHOB as spheres
  output2 << "\n\nselect HPHOB2, RigidClustObj and resn \"XXX\" " << endl;
  output2 << "show spheres, HPHOB2" << endl;
  output2 << "set sphere_scale=0.4" << endl;
  output2 << "deselect HPHOB2" << endl;

  // Hide Tethers
  // output2<<"hide everything, HPHOB2"<<endl;

  /*
   *     Write the FC visualization script
   */

  //  output3.open(file_FC_pml,ios::out);
  output3 << "delete FlexClustObj" << endl;
  output3 << "delete FC*" << endl;
  output3 << "delete smallFC*" << endl;
  output3 << "delete RC*" << endl;
  output3 << "delete HPHOB3" << endl;
  output3 << "load " << file_flexibility << ", FlexClustObj" << endl << endl;
  // Load molecule (*flex*.pdb)
  // Process the flexible clusters => (min_label*-1) through
  // (max_label-min_label) Note that min_label is negative (at least -1). So, we
  // have to first start with the cluster that has label = -1 to label =
  // min_label

  for (i = (min_label * -1) + 1, j = 1; i <= max_label - min_label; i++, j++) {
    if (clstBoundary[i] - clstBoundary[i - 1] > 6) {
      output3 << "\nselect FC" << j << ", FlexClustObj and id "
              << clstBoundary[i - 1] + 1 << "-" << clstBoundary[i];
      output3 << "\ncolor " << FCcolors[(j % 20) + 1] << ", FC" << j;

      //      output3<<"\npreset.ball_and_stick('RC"<<i<<"')"<<endl;
      //      Just show the main-chain by default as a ribbon or "cartoon" in
      //      PyMol

      output3 << "\nshow sticks, FC" << j;
      output3 << "\nshow cartoon, FC" << j;
      output3 << "\ncartoon tube, FC" << j;

      output3 << endl;
    }
  }

  /*
   * Combine all the clusters of size <= 6 intoa single cluster
   */
  flag = 'F';
  for (i = (min_label * -1) + 1, j = 1; i <= max_label - min_label; i++, j++) {
    if (clstBoundary[i] - clstBoundary[i - 1] <= 6) {
      if (flag == 'F') {
        // The rest of the atoms will be colored with a default color
        output3 << "\nselect smallFCs" << j / 50 << ", FlexClustObj and id "
                << clstBoundary[i - 1] + 1 << "-" << clstBoundary[i];
        flag = 'T';
      } else {
        output3 << "+" << clstBoundary[i - 1] + 1 << "-" << clstBoundary[i];
      }
    }

    if ((j + 1) % 50 == 0) {
      flag = 'F'; // To prevent string overflow of clust indices in PyMol
      output3 << "\ncolor brown, smallFCs" << j / 50;
      output3 << "\nshow sticks, smallFCs" << j / 50;
      output3 << "\nshow cartoon, smallFCs" << j / 50;
      output3 << "\ncartoon tube, smallFCs" << j / 50 << endl;
    }
  }

  if (j % 50 != 0) // Last smallFC object
  {
    output3 << "\ncolor brown, smallFCs" << j / 50;
    output3 << "\nshow sticks, smallFCs" << j / 50;
    output3 << "\nshow cartoon, smallFCs" << j / 50;
    output3 << "\ncartoon tube, smallFCs" << j / 50 << endl;
  }

  output3 << endl;

  //      All the RCs are grouped under one cluster and colored blue
  output3 << "\n\nselect RC, FlexClustObj and id 1-"
          << clstBoundary[(min_label * -1) - 1];
  output3 << "\ncolor blue, RC";
  output3 << "\nshow lines, RC";

  output3 << endl;

  //	All uncorrelated motions are grouped as one object and colored white
  output3 << "\n\nselect dangle, FlexClustObj and id "
          << clstBoundary[(min_label * -1) - 1] + 1 << "-"
          << clstBoundary[(min_label * -1)];
  output3 << "\ncolor white, dangle";
  output3 << "\nshow lines, dangle";

  output3 << endl;

  //  output3<<"\nset stick_radius, 0.01, RC";
  //  output3<<"\nset stick_transparency, 0.8, RC";

  //	Show HPHOB as spheres
  output3 << "\n\nselect HPHOB3, FlexClustObj and resn \"XXX\" " << endl;
  output3 << "show spheres, HPHOB3" << endl;
  output3 << "set sphere_scale=0.4" << endl;
  output3 << "deselect HPHOB3" << endl;

  //	Hide Tethers
  //  output3<<"hide everything, HPHOB3"<<endl;

  delete[] label;
  delete[] clst;
  delete[] atmIndex;
  delete[] clstBoundary;
  delete[] atm1;
  delete[] atm2;
  delete[] abs_flex;
  delete[] scaled_flex;

  output.close();
  output1.close();
  output2.close();
  output3.close();

  input.close();
  input1.close();
  input2.close();

  return 0;
}

/* ------------------------- END OF CODE ----------------------------- */
