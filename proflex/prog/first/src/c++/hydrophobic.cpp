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

/********************************************************************************/
/* 05.09.02 AJR modified and improved the search methods in not_2nd_neighb, */
/* not_3rd_neighb, and find_hydrophob --> this no longer is the rate limiting */
/* subroutine. */
/* 03.20.02 AJR based upon makeHP.c of Brandon M. Hespenheide */
/*                                                                              */
/* 11.15.00 BMH */
/*                                                                              */
/* 09.06.01 AJR Modified code to not assign an Hphob contact to metal CA atoms
 */
/*          Also introduced a check for chain_ids to print or not print them as
 */
/*          dictated by the input datasetfile. */
/*                                                                              */
/* 02.14.03 BMH  Finished cleaning up code (easier to read, less logic taking */
/*               up cycles, ...). Also added some new routines for testing */
/*               alternative definitions of a hydrophobic interaction. */
/*                                                                              */
/* vdwaals_4HP_1HB.c   Written in C++ */
/*                                                                              */
/* This program modifies a *_proflexdataset file. Hydrophobic contacts    */
/* are identified using variable distance criteria. These HP-contacts are */
/* represented as hydrogen bonds, with 4 "psuedo-atoms"between the donor and */
/* the bonding hydrogen, resulting in 1 constraint for the bond. */
/*                                                                              */
/* Hydrogen bonds are identifies using the REMARK:HB lines in the input *_FIRST
 */
/* dataset file. These bonds are modifed to have a single psuedo atom, adding */
/* a degree of freedom back to the system. */
/********************************************************************************/
/*  Program modified by AJR 09.06.01 to not assign an Hphob contact to metal CA
 */
/*  atoms.  Also introduced a check for chain_ids to print or not print them as
 */
/*  dictated by the input datasetfile. */
/********************************************************************************/

/**********************************************************************/
/* Measure the distance between the atom centers of the two atoms in  */
/* question (cp and tp). The van der Waals radii are hard coded in    */
/* this subroutine. If the sum of the radii + R_factor is below a     */
/* given cutoff, accept the interaction as a hydrophobic tether.      */
/**********************************************************************/
#include <cstring>
int list::determine_vdWaals_contact(node *cp, node *tp) {

  float atom_atom_dis = 0.0, vdr1 = 0.0, vdr2 = 0.0, aad2, sum_of_vdWaals = 0.0;

  int j, so, is_contact = 1;

  char atom1[6], atom2[6];

  strcpy(atom1, cp->r1.aname);
  strcpy(atom2, tp->r1.aname);
  so = cp->r1.ri_sr_no;

  if (atom1[2] == 'S') {
    vdr1 = 1.80;
  } else if (atom1[2] == 'C') {
    vdr1 = 1.70;
    if (strcmp(atom1, "  CA ")) {
      if (!strcmp(cp->r1.field1, "HETATM")) {
        if (!strcmp(cp->r1.rname, " CA") || !strcmp(cp->r1.rname, "CA ")) {
          vdr1 = 0.0;
        }
      }
    }
  }
  if (atom2[2] == 'S') {
    vdr2 = 1.80;
  } else if (atom2[2] == 'C') {
    vdr2 = 1.70;
    if (strcmp(atom2, "  CA ")) {
      if (!strcmp(tp->r1.field1, "HETATM")) {
        if (!strcmp(tp->r1.rname, " CA") || !strcmp(tp->r1.rname, "CA ")) {
          vdr2 = 0.0;
        }
      }
    }
  }

  sum_of_vdWaals = vdr1 + vdr2;
  aad2 = 0;
  for (j = 1; j <= 3; j++) {
    aad2 += pow((cp->r1.coord[j] - tp->r1.coord[j]), 2);
  }
  atom_atom_dis = sqrt(aad2);

  /**************************************************/
  /* the value 2.2 Angstrom was selected to exclude */
  /* covalently bonded atom pairs.                  */
  /**************************************************/
  if ((atom_atom_dis < 2.20) || (atom_atom_dis > ((sum_of_vdWaals) + R_factor)))
    is_contact = 0;

  return (is_contact);
}
/**********************************************************************/

/**********************************************************************/
/* Check to see if atoms so and sf are covalently bound to eachother. */
/* Returns 1 if the atoms ARE covalently bonded, zero otherwise.      */
/* "mso" is the multiplicity of atom "so".                            */
/**********************************************************************/
int list::test_for_covalent(int so, int sf, int mso) {

  int i = 0;

  for (i = 1; i <= mso; i++) {
    if (link_noHB[so][i] == sf)
      return (sf);
  }

  return (0);
}
/**********************************************************************/

/**********************************************************************/
/* modifed 05.09.02 AJR to be a faster search. In reality you have 2  */
/* atoms, s1 & s2 -- all that needs to be searched are the set of     */
/* their neighbors, not every atom in the network.                    */
/* Returns 1 if the atoms are NOT second nearest neighbors.           */
/**********************************************************************/
int list::not_2nd_neighb(int s1, int s2) {

  int sa = 0, m = 0;

  for (m = 1; m <= mult[s1]; m++) {
    sa = link_noHB[s1][m];
    if (p[sa]->r1.aname[2] != 'H' && p[sa]->r1.aname[2] != 'D') {
      if (test_for_covalent(s2, sa, mult[s2]))
        return (0);
    }
  }

  return (1);
}
/**********************************************************************/

/**********************************************************************/
/* Check to see if atoms so and sf are third nearest neighbors (that  */
/* is, if they're connected by two intervening bonds. Returns 1 if the*/
/* atoms are NOT 3rd nearest neighbors.                               */
/**********************************************************************/
int list::not_3rd_neighb(int so, int sf) {

  int i = 0, sa = 0;

  for (i = 1; i <= mult[so]; i++) {
    sa = link_noHB[so][i];
    // exclude hydrogens, they're monovalent
    if (p[sa]->r1.aname[2] != 'H' && p[sa]->r1.aname[2] != 'D') {
      if (!not_2nd_neighb(sa, sf))
        return (0);
    }
  }

  return (1);
}
/**********************************************************************/

/**********************************************************************/
/* Alternative hydrophobic tether identifier. Only add a tether if    */
/* both hydrophobic atoms are bonded to only other carbon or sulfur   */
/* atoms.                                                             */
/**********************************************************************/
int list::only_bonded_to_Hphobes(int so, int sf) {

  int a = 0;

  for (a = 1; a <= mult[so]; a++) {
    if (p[link_noHB[so][a]]->r1.aname[2] != 'C' &&
        p[link_noHB[so][a]]->r1.aname[2] != 'S' &&

        // AJR 06.24.03 added a line to omit hydrogen atoms as counting against
        // - Sameer
        p[link_noHB[so][a]]->r1.aname[2] != 'H')

      return (0);
  }

  for (a = 1; a <= mult[sf]; a++) {
    if (p[link_noHB[sf][a]]->r1.aname[2] != 'C' &&
        p[link_noHB[sf][a]]->r1.aname[2] != 'S' &&

        // AJR 06.24.03 added a line to omit hydrogen atoms as counting against
        // - Sameer
        p[link_noHB[sf][a]]->r1.aname[2] != 'H')

      return (0);
  }

  return (1);
}
/**********************************************************************/

/**********************************************************************/
/* The following routine allows the user to select Hphohic tethers    */
/* between atoms that have a specific atomic density of Hphobic atoms */
/* within a given sphere radius.                                      */
/**********************************************************************/
int list::Hphobe_atomic_density(int so, int sf) {
  std::cerr << "Unimplemented function in hydrophobic.cpp!!!" << std::endl;
  return 0;
}
/**********************************************************************/

/**********************************************************************/
/* Identify hydrophobic tethers.                                      */
/**********************************************************************/
int list::find_hydrophob(void) {

  int i, so, sf, **temp_hpteth;
  int iix, iiy, iiz, jx, jy, jz, kx, ky, kz, nval;
  char atype1[6], atype2[6];
  ofstream xofil;

  nval = no_atoms;
  // if(no_atoms >= 20000) { nval = (int)(no_atoms/2); }
  temp_hpteth = (int **)new int *[nval + 1];
  for (i = 0; i <= nval; i++) {
    //  temp_hpteth = (int **)new int *[no_atoms+1];
    // for(i=0; i<=no_atoms;i++){
    *(temp_hpteth + i) = (int *)new int[2];
  }
  int hpcount = 0;

  /* ################################################################AJR
     03.21.02 Identify the Hydrophobic contacts and place them into the array
     hp_tether[hpcount][i]
     ######################################################################################*/
  /* Cleaned up a bit. BMH 9.25.02 */

  /************************************************************/
  /* Hphob contact found by using "Hash code data structure   */
  /************************************************************/
  for (i = 1; i <= no_atoms; i++) {
    so = p[i]->r1.ri_sr_no;
    strcpy(atype1, p[i]->r1.aname);
    if (atype1[2] == 'C' || atype1[2] == 'S') {

      iix = (int)((p[i]->r1.coord[1] - xmin) / grdlen);
      iiy = (int)((p[i]->r1.coord[2] - ymin) / grdlen);
      iiz = (int)((p[i]->r1.coord[3] - zmin) / grdlen);

      /**********************************************************************/
      /* The following code uses the "hash" routine employed in various     */
      /* spots in FIRST. The method relies on placing a 3D grid on the prot.*/
      /* where the grid spacing is grdlen (currently 4.7 Ang.). The first   */
      /* atom that is found in a given grid cube is stored in the grid[][][]*/
      /* array. Every subsequent atom in that same cube is found using the  */
      /* "linked array" named chain[] here. To find the next atom in the    */
      /* given cube, access chain[] using the first atom in that cube. To   */
      /* get the third atom, access chain[] with the second atom num, and so*/
      /* on. The three initial for loops (over jx, jy, and jz) allow you to */
      /* check all the cubes in the grid adjacent to the cube in which the  */
      /* reference atom, so, is located.                                    */
      /* The for loop searches over all 27 grid cubes adjacent to cube in   */
      /* which point "so" resides.                                          */
      /**********************************************************************/
      for (jx = -1; jx <= 1; jx++) {
        for (jy = -1; jy <= 1; jy++) {
          for (jz = -1; jz <= 1; jz++) {

            kx = (iix + jx) % 32;
            ky = (iiy + jy) % 32;
            kz = (iiz + jz) % 32;

            sf = grid[kx][ky][kz];

            while (sf != 0) {

              if (sf > so) { // prevents double counting of tethers.
                strcpy(atype2, p[sf]->r1.aname);

                if ((atype2[2] == 'C' || atype2[2] == 'S') &&
                    determine_vdWaals_contact(p[i], p[sf]) &&
                    not_2nd_neighb(so, sf) && not_3rd_neighb(so, sf) &&
                    only_bonded_to_Hphobes(so, sf)) {
                  hpcount++;
                  temp_hpteth[hpcount][0] = so;
                  temp_hpteth[hpcount][1] = sf;
#ifdef DEBUG_HPHOB
                  cout << so << " " << sf << endl;
#endif
                }
              }
              sf = chain[sf];
            }
          }
        }
      }
    }
  }

#ifdef DEBUG_HPHOB
  cout << "R_factor " << R_factor << " hphobes " << hpcount << endl;
#endif
  // exit(0);
  /************************************************************/
  /* copy the temp_hpteth table to new memory optimized array */
  /************************************************************/
  hp_tether = (int **)new int *[hpcount + 1];

  for (i = 0; i <= hpcount; i++) {
    *(hp_tether + i) = (int *)new int[2];
  }

  for (i = 1; i <= hpcount; i++) {
    hp_tether[i][0] = temp_hpteth[i][0];
    hp_tether[i][1] = temp_hpteth[i][1];
  }

  for (i = 0; i <= nval; i++) {
    delete[] * (temp_hpteth + i);
  }

  return (hpcount);
}
/**********************************************************************/

/**********************************************************************/
void list::HPaugment(node **node_p) {
  /* AJR 03.22.02 this routine not yet implemented */
  int a, i, numatoms, nHPat, *tempbond, **templink, so, sf;

  float unit_vector[3], vector_one[3], vector_two[3], vector_three[3];

  node *fresh, *endpt;

  nHPat = no_tether * 3;
  numatoms = no_atoms + no_tether * 3;

  tempbond = new int[numatoms + 1];
  templink = (int **)new int *[numatoms + 1];
  for (i = 0; i <= numatoms; i++) {
    if (i <= no_atoms) {
      tempbond[i] = mult[i];
    }
    *(templink + i) = (int *)new int[maxr];
  }

  endpt = start;

  int nps = no_atoms, rescount = 1;
  fresh = NULL;

  for (i = 1; i <= no_tether; i++) {
    so = hp_tether[i][0];
    sf = hp_tether[i][1];

    /**********************************/

    for (a = 1; a < 4; a++) {
      unit_vector[a] = node_p[so]->r1.coord[a] - node_p[sf]->r1.coord[a];
      vector_one[a] = node_p[sf]->r1.coord[a] + (0.25 * unit_vector[a]);
      vector_two[a] = node_p[sf]->r1.coord[a] + (0.5 * unit_vector[a]);
      vector_three[a] = node_p[sf]->r1.coord[a] + (0.75 * unit_vector[a]);

      fresh = new node;

      strcpy(fresh->r1.field1,
             "HETATM"); // PDB field identifier (ie. ATOM, HETATM, REMARK)
      //    strcpy(fresh->r1.sr_no,line+6,5);     // Atom number
      strcpy(fresh->r1.aname, "  X  "); // Atom name
      strcpy(fresh->r1.altloc, " ");    //
      strcpy(fresh->r1.rname, "XXX");   // Amino-acid name
      strcpy(fresh->r1.chain, " T");    // Chain ID
      // strcpy(string, line+66, 5 );          // some kind of internal FIRST
      // identifier.
      fresh->r1.mod_res_no = rescount; //
      fresh->r1.DAH_type = 'N';
      fresh->r1.coord[1] = vector_one[a];
      fresh->r1.coord[2] = vector_two[a];
      fresh->r1.coord[3] = vector_three[a];
      strcpy(fresh->r1.occupancy, "  0.00");   // Temperature factor (B_value)
      strcpy(fresh->r1.temp_string, "  0.00"); //
      fresh->r1.occ = 0.00;
      fresh->r1.temperature = 0.00;
      fresh->r1.ri_sr_no = nps++;

      // AJR 03.21.02 temporary fix:
      strcpy(fresh->r1.res_no, " 100"); // Residue number
      strcpy(fresh->r1.sr_no, " 1252"); // Atom number
      strcpy(fresh->r1.code, "    ");
      strcpy(fresh->r1.strx, "  51.000"); // X-coordinate
      strcpy(fresh->r1.stry, "  11.000"); // Y-coordinate
      strcpy(fresh->r1.strz, "  -1.000"); // Z-coordinate
      //    strcpy(fresh->r1.strx,line+30,8);     // X-coordinate
      // strcpy(fresh->r1.stry,line+38,8);     // Y-coordinate
      // strcpy(fresh->r1.strz,line+46,8);     // Z-coordinate
      // cout << fresh->r1.field1 <<"  "<<nps<<endl;

      fresh->next = fresh->prior = NULL;
    }

    i = no_tether;
    rescount++;
    //    npseudo+=3;
  }
  // cout<<no_atoms<<"  "<<nps<<" "<<rescount<<"  "<<no_tether<<endl;
  //   cout <<unit_vector[1]<<endl;
  //-----------------------------------------Calculate new (maximum)
  // multiplicity for(i=0;i<=no_atoms;i++)
  //{
  //  temp_mult[i]=mult[i];
  //}
  //  for(i=no_atoms+1;i<=
  //  mult_net = new int[no_atoms+1];
}

void list::create_newPSEUDO(node *temp1, node *temp2) {

  //  node *temp1, *temp2, *tempo;
  /*fresh = new node;
  fresh->add_newPSEUDO();
  fresh->next = fresh->prior=NULL;*/
  // int chainflg = 0;
  // float vector_four[3];

  int last_resno, a = 0;

  float unit_vector[3], vector_one[3], vector_two[3], vector_three[3];
  char str_resno[5];
  node *tempo;

  //  cout << lastno <<"\t"<<no_atoms<<endl;
  // cout <<resval<<"\t" <<largest_res<<endl;

  /*char
    cid;

  FILE
    *hp_user_file;
    */
  /*  added by AJR 09.06.01 to make things look nicer with or without chain-ids
   */
  //  if( strncmp(test_point.chain, "  ",2) == 0) {chainflg = 1;}
  /*  if( strncmp(test_point.chain, " ",1) == 0) {chainflg = 1;}
   *if(chainflg == 1){  cid = ' ';  }
   *else {    cid = 'G';  }*/

  /*int flagt = 0;

  tempo = last;
  while(flagt <2) {
    if(so == tempo->r1.ri_sr_no) {
      flagt++;
      temp1 = tempo;
    }
    if(sf == tempo->r1.ri_sr_no) {
      flagt++;
      temp2 = tempo;
    tempo=tempo->prior;
    }
    }*/

  tempo = start;
  strcpy(str_resno, tempo->r1.res_no);
  last_resno = atoi(str_resno);
  // cout << last_resno <<endl;

  for (a = 1; a < 4; a++) {
    //  p[so]->r1.coord[1]
    unit_vector[a] = temp1->r1.coord[a] - temp2->r1.coord[a];
    vector_one[a] = temp2->r1.coord[a] + (0.25 * unit_vector[a]);
    vector_two[a] = temp2->r1.coord[a] + (0.5 * unit_vector[a]);
    vector_three[a] = temp2->r1.coord[a] + (0.75 * unit_vector[a]);
    //    vector_four[a]  = temp2->r1.coord[a] + (0.8*unit_vector[a]);
    // hp_points[*hp_point_counter][a] = current_atom.coords[a] +
    // (0.5*unit_vector[a]);
    //      unit_vector[a]  = test_point.coords[a] - current_atom.coords[a];
    //      vector_one[a]   = current_atom.coords[a] + (0.2*unit_vector[a]);
    //      vector_two[a]   = current_atom.coords[a] + (0.4*unit_vector[a]);
    //      vector_three[a] = current_atom.coords[a] + (0.6*unit_vector[a]);
    //      vector_four[a]  = current_atom.coords[a] + (0.8*unit_vector[a]);
    //      hp_points[*hp_point_counter][a] = current_atom.coords[a] +
    //      (0.5*unit_vector[a]);
  }

  //  (*hp_point_counter)++;
  /*
  fprintf( output_file, "HETATM%5d  X   XXX %c%4d     %7.3f %7.3f
%7.3f  1.00  1.00      N\n", *next_atom_num, cid, *residue_number,
vector_one[0], vector_one[1], vector_one[2] ); fprintf( output_file, "HETATM%5d
X   XXX %c%4d     %7.3f %7.3f %7.3f  1.00  1.00      N\n", *next_atom_num+1,
cid, *residue_number, vector_two[0], vector_two[1], vector_two[2]
);
  fprintf( output_file, "HETATM%5d  X   XXX %c%4d     %7.3f %7.3f
%7.3f  1.00  1.00      N\n", *next_atom_num+2, cid, *residue_number,
vector_three[0], vector_three[1],vector_three[2] );
  (*residue_number)++;
*/
}

/**************************************************/
void list::write_HPCFbond(void) {

  int i, so, latom = no_atoms;
  ofstream xofil;
  xofil.open(outputfile, ios::app);

  if (!xofil) {
    cout << "Error in opening outputfile" << endl;
    exit(4);
  }

  for (i = 1; i <= no_tether; i++) {
    latom++;
    so = hp_tether[i][0];
    xofil << "REMARK:CF" << setw(6) << so << setw(6) << latom << endl;
    xofil << "REMARK:CF" << setw(6) << latom << setw(6) << latom + 1 << endl;
    xofil << "REMARK:CF" << setw(6) << latom + 1 << setw(6) << latom + 2
          << endl;
    latom = latom + 2;
  }

  xofil.close();
}
/**************************************************/
void list::write_HPHBbond(void) {

  int i, sf, patom = no_atoms;
  int ahb = nhb;
  float ephob = -9.99999;
  ofstream xofil;
  xofil.open(outputfile, ios::app);

  if (!xofil) {
    cout << "Error in opening outputfile" << endl;
    exit(4);
  }

  for (i = 1; i <= no_tether; i++) {
    patom = patom + 2;
    //    so = hp_tether[i][0];
    sf = hp_tether[i][1];
    xofil << "REMARK:HB" << setw(5) << ahb + i << setw(13) << setprecision(6)
          << ephob << setw(8) << patom << setw(8) << patom + 1 << setw(8) << sf
          << "    PH hydr phob" << endl;
    patom++;
  }

  xofil.close();
}
/**************************************************/

void list::write_HPatom(node **node_p) {

  int a, i, numatoms, nHPat, so, sf;

  float unit_vector[3], vector_one[3], vector_two[3], vector_three[3];

  node *endpt;

  nHPat = no_tether * 3;
  numatoms = no_atoms + no_tether * 3;

  ofstream xofil;

  xofil.open(outputfile, ios::app);

  if (!xofil) {
    cout << "Error in opening outputfile" << endl;
    exit(4);
  }

  endpt = start;

  int nps = no_atoms, rescount = 1, jt = 0;

  for (i = 1; i <= no_tether; i++) {
    so = hp_tether[i][0];
    sf = hp_tether[i][1];
    /**********************************/

    for (a = 0; a < 3; a++) {
      //  node_p[so]->r1.coord[1]
      unit_vector[a] =
          node_p[so]->r1.coord[a + 1] - node_p[sf]->r1.coord[a + 1];
      vector_one[a] = node_p[sf]->r1.coord[a + 1] + (0.75 * unit_vector[a]);
      vector_two[a] = node_p[sf]->r1.coord[a + 1] + (0.5 * unit_vector[a]);
      vector_three[a] = node_p[sf]->r1.coord[a + 1] + (0.25 * unit_vector[a]);
    }

    jt++;
    xofil << "HETATM" << setw(5) << nps + jt << "  X   XXX T" << setw(4)
          << rescount << setw(12) << setprecision(3)
          << setiosflags(ios::showpoint | ios::fixed | ios::right)
          << vector_one[0] << setw(8) << setprecision(3) << vector_one[1]
          << setw(8) << setprecision(3) << vector_one[2]
          << "  0.00  0.00      N" << endl;
    jt++;
    xofil << "HETATM" << setw(5) << nps + jt << "  X   XXX T" << setw(4)
          << rescount << setw(12) << setprecision(3) << vector_two[0] << setw(8)
          << setprecision(3) << vector_two[1] << setw(8) << setprecision(3)
          << vector_two[2] << "  0.00  0.00      N" << endl;
    jt++;
    xofil << "HETATM" << setw(5) << nps + jt << "  X   XXX T" << setw(4)
          << rescount << setw(12) << setprecision(3) << vector_three[0]
          << setw(8) << setprecision(3) << vector_three[1] << setw(8)
          << setprecision(3) << vector_three[2] << "  0.00  0.00      N"
          << endl;
    rescount++;
  }

  xofil.close();
}
/**************************************************/
