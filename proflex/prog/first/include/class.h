/*******************************************************************************
*  MSU ProFlex, formerly called FIRST, is a software developed to predict and  *
*  analyze protein flexibility.                                                *
*  This source file is a part of MSU ProFlex.                                  *
*                                                                              *
*  Copyright (C) 1997 - 2008, Michigan State University.                       *
*                                                                              *
*  This program is free software; you can redistribute to academic users only, *
*  it and/or modify it under the terms of the GNU General Public License,      *
*  version 2, as published by the Free Software Foundation.                    *
*                                                                              *
*  This program is distributed in the hope that it will be useful,             *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of              * 
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
*  GNU General Public License for more details.                                *
*                                                                              *
*  You should have received a copy of the GNU General Public License           *
*  along with this program; if not, write to the Free Software                 *
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA,  *
*  or see http://www.gnu.org/licenses/gpl.txt                                  *
*******************************************************************************/

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <cctype>
#include <vector>

using namespace std;

//--------------------------------------------------Define constants
#define JUNK_INT -999
#define MAXATOMS 100000

const int   lookup_atoms= 7;
const int   no_res = 20;           	// The standard 20 known residues
const int   maxnopt = 5;		//5+1: 0 - f,s; 1 - default; 
					//2,3,4 - Hbond; 5 - tethers
const int   maxr = 13;			// maximum number of ranges in cluster sizes for seting color preferences.

const float s_admaxd = 4.0;		// maximum distance allowed in H-bond between acceptor-donor pair,
					// with at least one of the atoms being Sulfur.
const float s_hamaxd = 3.0;		// maximum distance allowed in H-bond between H-acceptor pair with
			// at least one of the acceptor or donor atoms being Sulfur.

const float sb_admaxd = 4.5; 		// max distance allowed in salt bridge between acceptor-donor pair
const float sb_hamaxd = 3.5;		// max distance allowed in salt bridge between H-acceptor pair

const float admaxd = 3.5; 		// maximum distance allowed in H-bond between acceptor-donor pair.
const float hamaxd = 2.5; 		// maximum distance allowed in H-bond between H-acceptor pair.

const float raddeg = 1.745329e-2;	// conversion factor from degrees to radians
/*-------------------------	
 * 2006:01	--- Sandeep
 *
 * The default value for THETA (Donor-H-Acptr) angle is set to 110 instead of 
 * the previous value: 80, which is not a true compliant of the (D-H-A) physical
 * alignment that is ruled by underlying chemistry
 *-------------------------
 */
const float dha_ang = 110*raddeg;	// minimum allowed donor-hydrogen-acceptor angle not involving sulfur
const float s_dha_ang = 110*raddeg;	// minimum allowed donor-hydrogen-acceptor angle involving sulfur
const float sb_dha_ang = 110*raddeg;    // minimum allowed donor-hydrogen-acceptor angle in salt bridge

/*------------------------  
 * 2006:01      --- Sandeep
 * Add new constants to impose rules on H-Acptr-preAcptr angle, DELTA.
 * The threshold inthis case is set to 90 degree
 *------------------------
 */
const float hap_ang = 90*raddeg;       	// minimum allowed hydrogen-acceptor-
					// preAcptr angle not involving sulfur
const float s_hap_ang = 90*raddeg;     	// minimum allowed hydrogen-acceptor-
					// preAcptr angle involving sulfur
const float sb_hap_ang = 90*raddeg;    	// minimum allowed hydrogen-acceptor-
					// preAcptr angle in salt bridge

const float grdlen= sb_admaxd+0.1;	/* The length of one side of a cube 
					   forming a coarse grained cell
					   This distance must be bigger than 
					   the longest Donor-Acceptor distance 
					*/

const float maxdist= 2.5;		// maximum distance between pair of atoms


// Global Variables
int no_atoms = 0, nh=0, ns=0, nfile;
int *map_array;  // Renamed since it conflicts with including std::map
size_t res_count=0;
int polar_count=0;
int **link_noHB;
int **link_net;
ifstream ifil;
ifstream qfil;
ofstream ofil;
ofstream hblist;
char  *path;
char residue_lib[130],dist_lookup[130],polar_lookup[130];
char inputfile[120],summon_fortran[120];
char outputfile[130],outfile[8][130];
char answer[201];
float max_energy = -1.0;// default value in kcal/mol. maximum energy a 
			// hydorgen bond can have to be considered in FIRST analysis.
			// alternatively, this defines the energy of the weakest bond
float hpd_R_cutoff = 0.0;
int   hpd_min = 0;
			// R_factor is the distance between atom's vw-surfaces
			// to be considered
float R_factor = 0.5;	// for hydrophobic interactions.
			// Changing to 0.5 from 0.25 - Sameer, Feb 13, 2004

int   hphobs = 0;
int   HBanalysis = 10;
int   dil_scheme = 0;
int   HPs = 0;

enum Atom {C,N,O,S,H,P,Oth};

//////////////////////////////////////////////////////////////////////
class record {

 public:

// AJR 05.08.02 new variables for MODEL (NMR) record handling.
  int   modelflag;

  int   ri_sr_no;  // The atom number of the protein atoms
  int   mod_res_no;
  int   hpd;      // Number of C or S atoms within R_cutoff of current atom
  int   bonded_to_Hphobe; // state variable indicating if all the current atom's neighbors are carbon or sulfur.
  float coord[4]; // The atomic coordinates of the protein atoms. They
                  // are indexed as coord[1-3] in the program.
  float occ, temperature;

  char  field1[7],sr_no[6],aname[6],altloc[2],rname[4],chain[3];
  char  res_no[5], code[5], strx[9],stry[9],strz[9];
  char  occupancy[7], temp_string[7], DAH_type;

  void getrecord(void);

  record(void) {
    modelflag = 0;
    int a = 0;
    for( a = 0; a < 4; a++ )
      coord[a] = 0.0;
    ri_sr_no = 0;
    mod_res_no = 0;
    occ = 0.0;
    temperature = 0.0;
    strcpy(field1,"      ");
    strcpy(sr_no,"     ");
    strcpy(aname,"     ");
    strcpy(altloc," ");
    strcpy(rname,"   ");
    strcpy(chain,"  ");
    strcpy(res_no,"    ");
    strcpy(code,"    ");
    strcpy(strx,"        ");
    strcpy(stry,"        ");
    strcpy(strz,"        ");
    strcpy(occupancy,"      ");
    strcpy(temp_string,"      ");
    DAH_type=' ';
    hpd = 0;
    bonded_to_Hphobe = 0;
  }
};

// what does this do
class conect {

public:

  int ci_sr_no;
  int i_conect_sr_no[4];

  void getconect(void);

  conect(void) {
    ci_sr_no = 0;
  }

};

/**********************************************************************/
/* what does this do                                                  */
/**********************************************************************/
class residue {

public:


  int A_mult[50], A_link[50][4], natom_sum, occu[50];

  char rname[4], A_type[50][4];

  //--------------------------------------------------------------//
  // residues (res) go from 0-n;				  //
  // all the data members of res (A_type,A_mult,A_link,occu)      //
  // 						   go from 1-n;   //
  //--------------------------------------------------------------//
  void get_residue(char *temp_res);

  residue(void) {
    int a = 0;
    natom_sum = 0;
    for( a = 0; a < 50; a++ )
      occu[a] = 0;
  }

};
/**********************************************************************/

// what does this do?

class DAHspecs {
public:

  int num;

  char rname[4],aname[25][5],DAH_type[25][2];

  void get_polar_record(char *temp_res);

  DAHspecs(void) {
    num = 0;
  }

};


class Distance {
public:

  float d[4];

};


class node : public record, public residue {

public:

  record r1;
  node *next, *prior;

  void add_rec_node(void);

};


//////////////////////////////////////////////////////////////////////
class conect_node :public conect {

public:

  conect c1;
  conect_node *next, *prior;

  void add_conect_node(void);

};

/////////////////////////////////////////////////////////////////////
// Why is this public record and residue?  Placing class names after the
// class you wish to define causes the defined classes to inherit all
// of the items from the base classes (those in the list -- here it means
// that the class "list" would also include all of the items from the 
// classes record and residue
//class list : public record, public residue {
class list{

public:

  int chaincount, altloccount, no_conf[20], largest_res, chain_change_count;
  int conf_option, choice[20];
  int grid[32][32][32],*chain,*h_list,h_list_count, *mult, *mult_net;
  int prev_sel_flag, no_tether, *hb16, **hp_tether,model_flag;
  int **nrot_bond, nrot_bond_count, **hbond, *hb_type, *hb_id, nhb, *hb_pick, real_hb_no;

  float xmin, ymin, zmin, resolution, *hb_energy;

  char chain_change[7][2],*conf_chosen,prev_sel[100];
  // this is not supported by a number of compilers -- or at least not well 
  // -- vanvoor4 April 13, 2010
  // char *screen[maxnopt+1];
  std::vector<std::string> screen_msgs;

  Atom a1,a2;
  Distance distance[lookup_atoms][lookup_atoms];
  node *start, *last, *chainptr[10], *altlocptr[20],*s_ptr[1000];
  node **p, *h_ptr[10000];

  void add_record(int);
  void show_rec_list(void);
  void show_conect_list(void);
  void make_chem_rec(void);
  void cal_maxmin(void);
  void del_conf(std::string, int);
  void make_connectivity(int);
  int  check_dist(int,int,float,char *,char *,int,char *,char *,char *,char *,
		 char *,char *,char *,char *, int);
  void write_output(int,int,int,char*,char*,char*,char*,char*,char*,char *,
		    char *,char*,char*);
  void make_connection(int *,int **, int, int);
  void set_DAH_type(int, int *);
  void find_Hbonds(int);
  void check(void);
  void read_proflexdataset(void);
  void pick_hbonds(int,int); //--- 2006:03 SN
  void place_hbonds(int);
  void pick_hbonds_options(void);
// added by AJR 03.22.02 for hydrophobic tethers
  int find_hydrophob(void);
  int determine_vdWaals_contact(node *, node *);
  int test_for_covalent(int, int, int);
  int not_2nd_neighb(int, int);
  int not_3rd_neighb(int, int);
  int only_bonded_to_Hphobes(int, int);
  int Hphobe_atomic_density(int, int);
  void within_hpd_sphere( node *, node *);
  void HPaugment(node **);
  void write_HPatom(node **);
  void write_HPCFbond(void);
  void write_HPHBbond(void);
//  void create_newPSEUDO(int, int, int, int);
  void create_newPSEUDO(node *, node * );


  list(void) {
    xmin = ymin = zmin = resolution = 0.0;
    nhb = 0;
    real_hb_no = JUNK_INT; // Sameer - to track real number of h-bonds. nhb also includes hphobic tethers
    prev_sel_flag=0;
    largest_res=0;
    chaincount=0;
    altloccount=0;
    chain_change_count = 0;
    nrot_bond_count=0;
// AJR 05.08.02 new variables for MODEL (NMR) record handling.
    model_flag=0;

    for(int i=0;i<20;i++) {
      no_conf[i]=0;
    }
    p     = NULL;
    start = NULL; // beggining element of a linked list
    last  = NULL; // last element of a linked list
  }
};

//////////////////////////////////////////////////////////////////////
class conect_list : public conect {

public:

  conect_node *start, *last;

  void add_conect(void);

};

//////////////////////////////////////////////////////////////////////

// This is very annoying -- one should not be using global variables; 
// especially those global to the entire program (as opposed being global to 
// one source file) -- vanvoor4
conect_list l2;
residue *res;
DAHspecs *polar;

/*class pebble : public list
{
	public:
		int *clst, *stress **link_hinge;

		void mapmolecule(void);
}
*/
