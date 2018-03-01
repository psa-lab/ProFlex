/*******************************************************************************
*  MSU ProFlex, formerly called FIRST, is a software developed to predict and  *
*  analyze protein flexibility.                                                *
*  This source file is a part of MSU ProFlex.                                  *
*                                                                              *
*  Copyright (C) 1997 - 2007, Michigan State University.                       *
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

/*-----------------------------------------------------------------------------
                                      Re-written by S.K. Namilikonda	06/2007
			Previous versions written by D.J. Jacobs and A.J. Rader

   This program reads in the "<prot>_allbonds" file and for each atom
   it picks the bond with the least flexibility index over all bonds 
   stemming from it and assign its value to that atom. This is a coarse
   grained description of the flexibility index over the entire protein. 
   This value (which ranges between -1 (min flex) to +1 (max flex)) 
   is then transformed into a new number and put in the B-value column 
   of a PDB file ranging between 0 to 99, so that Insight II can color 
   code the protein in accordance to the flexibility index.

-----------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>

using namespace std;
#define MAXATM 10000
#define MAX_FILENAME 200

int main(void)
{
//-------------------------------------------------- declare variables
    ifstream input;
    ofstream output;
    char     intro[10];
    char     line[90];
    int      count,nHP=0;
    const char error_open[] = "\n *** main(): unable to open file: ";
//-------------------------------------------------- to store input/output data
    int      *atm1,*atm2;
    int	     max_atm_num,tmp,a1,a2,clust,i;
    float    *abs_flex,*scaled_flex,wt;
//-------------------------------------------------- preset names for now
  const char clear[] = "clear";
  const char ls_bond[] = "ls *allbond*";

//-------------------------------------------------- define arrays FORTRAN style
  char      file_weight[MAX_FILENAME];
  char      file_chem[MAX_FILENAME];
  char      file_flexibility[MAX_FILENAME];
  char      ans, answer[23];
  char      sso[6],ssf[6];
  char      *fend1 = "proflexdataset";
  char      *fend2 = ".pdb";  
  short int num_full,num_base,num_PROFLEX;
  int       flagHP = -1;

//--------------------------------------------- initialize file-name arrays to NULL chars

  memset( (void*)file_weight , '\0', MAX_FILENAME);
  memset( (void*)file_chem , '\0', MAX_FILENAME);
  memset( (void*)file_flexibility , '\0', MAX_FILENAME);
//------------------------------------------------------------------- cosmetics
  system( clear );
  cout << "\n\n     Enter one of the following \"allbonds\" files\n" << endl;
  system( ls_bond );
  cout << endl;

  cin.getline(file_weight,MAX_FILENAME);
  num_full = strlen(file_weight);
  num_base = -1;

  cout << endl;

  /*
   * Figure out the prefix of "allbonds" file
   */
  for(i=0; i<MAX_FILENAME; i++)
   {
     file_chem[i] = file_weight[i];
     file_flexibility[i] = file_weight[i];
     if( file_weight[i] == '_' && file_weight[i+1] == 'a' \
         && file_weight[i+2] == 'l')
       {
	  num_base = i;
          break;
       }
   }

   if( num_base == -1 ) return 0;

   file_chem[i+1] = '\0';

   strcat(file_chem,fend1);

   num_PROFLEX = strlen(file_chem);

   file_flexibility[num_base+1] = 'f';
   file_flexibility[num_base+2] = 'l';
   file_flexibility[num_base+3] = 'e';
   file_flexibility[num_base+4] = 'x';
   num_full = num_full - 4;
   for( i=1; i<=4; i++)
     {
       file_flexibility[num_base+4+i]= file_weight[num_full++];
     }
   strcat(file_flexibility,fend2);
   file_flexibility[num_base+13]= '\0';

   cout<< " If hydrophobic tethers are present in the "<<fend1<< " file, "<<endl;
   cout<<" do you wish to output them to the modified pdbfile: "<<file_flexibility
       <<" ? "<<endl<<" Enter (y/n) {Default is y}"<<endl;
   cin.get(answer, 3);
   cin.clear();
   ans = answer[0];
   if(ans =='n' || ans =='N') 
	flagHP = 1;
   else 
	flagHP=0;

// --- Parse the input file and store first three columns as 
// --- the two atom numbers of each bond and its flexibility weight
// --- To do the above, figure out the max_atm_# and allocate memory 
// --- to store the final flexibility index of each atom
  count = 0;
  max_atm_num = 0;

  input.open(file_weight, ios::in);
  if( input.fail() )
    {
     cerr << error_open << file_weight << "\n\n";
     exit(-1);
    }

  input >> a1 >> a2 >> wt >> clust;  
  while( !input.eof() )
   {
	count++;
 
	if(a1 > max_atm_num)
	  max_atm_num = a1;
	if(a2 > max_atm_num)
	  max_atm_num = a2;
	
	input >> a1 >> a2 >> wt >> clust;
   }

  input.close();
  input.clear();
  input.open(file_weight, ios::in);
  if( input.fail() )
    {
     cerr << error_open << file_weight << "\n\n";
     exit(-1);
    }

   // --- Now that we have the size of the input file, 
   // --- allocate memory to store all the records
   atm1 = new int [count];
   atm2 = new int [count];
   abs_flex = new float [count];
   scaled_flex = new float [max_atm_num];
   
//----------------------------------------------- read data from allbonds file
   for(i = 0; i < count; i++)
    {
       input >> atm1[i] >> atm2[i] >> abs_flex[i] >> clust;       
       if(i < max_atm_num)
	 {
	   scaled_flex[i] = 1.0; // initiailize to max flex wt allowed (1.0)
	 }
    }

    input.close();
    input.clear();
//----------------------------------------------- find min flexibility index
   for(i = 0; i < count; i++)
    {
      if(abs_flex[i] < scaled_flex[atm1[i]-1]) 	
	{
	  scaled_flex[atm1[i]-1] = abs_flex[i];
	}
      if(abs_flex[i] < scaled_flex[atm2[i]-1]) 	
	{
	  scaled_flex[atm2[i]-1] = abs_flex[i];
	}
    }
//----------------------------------------------- Now, scale the weight factors
//						  from (-1,+1) to (0,99)
   for(i = 0; i < max_atm_num; i++)
    {
      	scaled_flex[i] = (scaled_flex[i]+1)*49.5;
    }


   input.open(file_chem, ios::in);
   if( input.fail() )
     {
       cerr << error_open << file_chem << "\n\n";
       exit(-1);
     }

    output.open(file_flexibility, ios::out);
    count = 0;
    input.getline(line,90);
     while( !input.eof() )
     {
        if(strncmp(line,"ATOM",4) == 0 || strncmp(line,"HETATM",6) == 0)
         {
	  if(flagHP && !strncmp(line+13,"X ",2) && !strncmp(line+17,"XXX",3))
	    { nHP++; }
	  else
	   {
        for(i=0; i<54; i++) output << *(line+i);
        output << "  1.00 " << setw(5) << setprecision(2)
               << setiosflags(ios::showpoint | ios::fixed | ios::right)
               << *(scaled_flex+count) << endl;

	   }
           count++;
	 }
     input.getline(line,90);
     }
  input.close();
  input.clear();
  output.close();
  if(!flagHP){
     cout << "\n\n\tThe number of ATOMS = " << count << endl;
  }
  else if(flagHP) {
     cout << "\n\n\tThe number of ATOMS = " << count - nHP << endl;
     cout<<"\n\tThe number of Tethers omitted = " << nHP/3 << endl;
  }


 return 0;
}

/* ------------------------- END OF CODE ----------------------------- */


