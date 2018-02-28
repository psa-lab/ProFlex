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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#define column(n)   (linebuf+(n)-1) /* manipulate linebuf by column */

/**********************************************************************/
/**********************************************************************/
void swap_bonds( int array[400][2], int A, int B ){
  int temp = 0, a = 0;
  for( a = 0; a < 2; a++ ){
    temp = array[A][a];
    array[A][a] = array[B][a];
    array[B][a] = temp;
  }
}
/**********************************************************************/
/**********************************************************************/
void sort_bonds( int array[400][2], int start, int end ){
  
  int a=0, last=0;
  
  if( start >= end )
    return;  
  
  swap_bonds( array, start, (start+end)/2 );
  last = start;
  for( a = start+1; a <= end; a++ ){
    if( array[a][0] < array[start][0] )
      swap_bonds( array, ++start, a );
  }
  swap_bonds( array, start, last );
  sort_bonds( array, start, last-1 );
  sort_bonds( array, last+1, end   );
}
/**********************************************************************/
/**********************************************************************/
void print_connections( FILE *pdb_file, FILE *ps_file, int y_translate[10] ) {

  int    a = 0, b = 0, c = 0, d = 0, f = 0, g = 0, i = 0,
         ps_line_number, 
         total_mainchain = 0, total_backbone = 0,
         main_chain_atoms[2000], backbone_atoms[2000],
         is_atom = 0, atom_num = 0,
         atom_one_is_good = 0, atom_two_is_good = 0,
         atom_1 = 0, atom_2 = 0, color = 0, 
	 top = 0, bottom = 0,
         N_terminus = 1, Count = 0,
         bond_list[400][2], temp = 0,
         current_index = 0, increment = 0,
         short_list[400], long_list[400],
         finished = 0;

  float  bond_wt_array[4000], bond_wt = 0.0, 
         flexibilities[4000], average[4000], 
         residue_flexibility = 0.0;

  char   linebuf[150], atom[4];

  double difference = 0.0;

  FILE   *bond_wts,
         *hbond_list;

  for(a = 0; a < 4000; a++ ){
    bond_wt_array[a] = 0; 
    bond_wt_array[a] = 0;
  }
  
  bond_wts = fopen( "bond_wts", "r" );

  ps_line_number = 500;

  /******************************************************************/
  /* Read in the chem file. Store the atom number of the main-chain */
  /* N, CA, C, O atoms for referencing the bond weight file. I've   */
  /* used two arrays here. The main_chain_atoms array stores the    */
  /* atom number for the carbonyl oxygen for use in the flexibility */
  /* index calculation, however, the oxygen is not displayed in the */
  /* output. The backbone array is used for the output.             */
  /******************************************************************/
  while( fgets( linebuf, sizeof(linebuf), pdb_file) != NULL ) {
    is_atom = !strncmp( linebuf, "ATOM", 4 );
    sscanf( linebuf+13, "%3s", atom );

    if( is_atom && 
	( !strcmp( atom, "N" ) ||
	  !strcmp( atom, "CA") ||
	  !strcmp( atom, "O")  ||
	  !strcmp( atom, "C" ) ) ) {
      
      main_chain_atoms[total_mainchain] = atoi( (linebuf+5) );
      total_mainchain++;
    }

    if( is_atom && 
	( !strcmp( atom, "N" ) ||
	  !strcmp( atom, "CA") ||
	  !strcmp( atom, "C" ) ) ) {

      backbone_atoms[total_backbone] = atoi( (linebuf+5) );
      total_backbone++;
    }
  }
  /******************************************************************/
  /******************************************************************/

  /******************************************************************/
  /* Read the bond weights file. The bonds are indexed by atom num- */
  /* ber in the original .pdb file.                                 */
  /******************************************************************/
  while( fgets( linebuf, sizeof(linebuf), bond_wts ) != NULL ){

    atom_1  = atoi( column(7) );
    atom_2  = atoi( column(16));
    bond_wt = atof( column(26));

    /****************************************************************/    
    /* check to see if both atom_1 and atom_2 are main-chain atoms. */
    /****************************************************************/    
    for( i = 0; i < total_mainchain; i++ ){
      
      if( atom_1 == main_chain_atoms[i] )
	atom_one_is_good++;
      if( atom_2 == main_chain_atoms[i] )
	atom_two_is_good++;
    }
    
    if( atom_one_is_good && atom_two_is_good ){
      bond_wt_array[atom_1] += bond_wt;
      bond_wt_array[atom_2] += bond_wt;
    }
    
    atom_one_is_good = atom_two_is_good = 0;
  }  
  /******************************************************************/
  /******************************************************************/

  /******************************************************************/    
  /* compute the actual flexibility index for the main-chain atoms. */
  /* Kind of cheating here. Since the flexibility index is assigned */
  /* to the bonds of a network, the index for an atom is determined */
  /* by averaging the index of all bonds incident with a given atom */
  /* The cheat is that for main chain atoms we only average those   */
  /* bonds which originate from other main-chain atoms. Here we in- */
  /* clude the carbonyl oxygen, requiring the extra array in the    */
  /* above calculations.                                            */
  /******************************************************************/
  rewind( pdb_file );
  
  while( fgets( linebuf, sizeof( linebuf ), pdb_file ) != NULL ){
    
    is_atom = !strncmp( linebuf, "ATOM", 4 );
	
    if( is_atom ){
      sscanf( linebuf+13, "%3s", atom );
      atom_num = atoi( column(7) );
      
      average[atom_num] = 5; /* error check value */
      
      /**************************************************************/
      /* For the main-chain N and Ca atoms, there are two bonds     */
      /* incident with other main-chain atoms. Divide the summed    */
      /* bond wts by 2.                                             */
      /**************************************************************/
      if( !strcmp( atom, "N" ) ||
	  !strcmp( atom, "CA") ){
	if( N_terminus ){
	  average[atom_num] = bond_wt_array[atom_num];
	  N_terminus = 0;
	}
	else
	  average[atom_num] = bond_wt_array[atom_num] / 2;
      }
  
      /**************************************************************/
      /* Since we include the main-chain oxygen in the calculation, */
      /* divide the summed bond wts for C by 3.                     */
      /**************************************************************/
      if( !strcmp( atom, "C" ) )
	average[atom_num] = bond_wt_array[atom_num] / 3;

      bond_wt_array[atom_num] = 0;
    }
  }

  for( f = 0; f < total_backbone; f++ )
    flexibilities[f] = average[ backbone_atoms[f] ];
  /******************************************************************/
  /******************************************************************/

  /******************************************************************/
  /* Output the flexibility results in linear format.               */
  /******************************************************************/
  top    = ps_line_number + 3;
  bottom = ps_line_number - 3;
  fprintf(ps_file,"Col00\n");
  
  for( g = 0; g < total_backbone; g++ ){

    /****************************************************************/
    /* Take the average of the three main-chain atoms in an amino-  */
    /* acid, and out this single value. This routine automatically  */
    /* assumes that the first main-chain atom in the protein is a N */
    /* A safe bet, but may cause some problems.                     */
    /****************************************************************/
    residue_flexibility += flexibilities[g];
    Count++;
    
    if( Count == 3 ){
      
      /**************************************************************/
      /* Reset the values to a scale between 0 and 10               */
      /**************************************************************/
      color = ((((residue_flexibility/3) + 1 ) * 10) + 0.99);
      
      if( color < 10 )
	fprintf(ps_file,"FlexCol0%d\n", color );
      else
	fprintf(ps_file,"FlexCol%d\n",  color );
      
      fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
	      g+110, top,
	      g+113, top,
	      g+113, bottom,
	      g+110, bottom );
      
      Count = residue_flexibility = 0.0;
    }
	
  }	
  /******************************************************************/
  /******************************************************************/

  a = 0;
  /******************************************************************/
  /* This part of the code reads in a list of hydrogen bonds for the*/
  /* protein, and displays the bonds as black lines connecting the  */
  /* amino-acids in the protein.                                    */
  /******************************************************************/
  hbond_list = fopen("hbond.list", "r" );

  while( fgets( linebuf, sizeof(linebuf), hbond_list) != NULL ){
    if( atoi( column(16) ) <= 103 && atoi( column(30) ) <= 103 ) {
      bond_list[a][0] = atoi( column(16) );
      bond_list[a][1] = atoi( column(30) );
      a++;
    }
  }

  /******************************************************************/
  /* The next routine sorts the hbonds into two groups. The first   */
  /* will be 'short-range' hbonds, those which connect residues less*/
  /* than 6 amino-acids apart. These will be displayed below the    */
  /* colored flexibility diagram created above. The remaining bonds */
  /* will be displayed above the line.                              */
  /******************************************************************/

  /****************************************************************/
  /* sort the array element such that the smaller of the two num  */
  /* bers at each element is first in the list.                   */
  /****************************************************************/
  for( b = 0; b < a; b++ ){
    if( bond_list[b][0] > bond_list[b][1] ){
      temp = bond_list[b][0];
      bond_list[b][0] = bond_list[b][1];
      bond_list[b][1] = temp;
    }
  }

  /****************************************************************/
  /* Sort the list of bonds according the first array elements    */
  /****************************************************************/
  sort_bonds( bond_list, 0, a-1 );

  for( b = 0; b < a; b++ ){

    difference = bond_list[b][0] - bond_list[b][1];
    if( abs(difference) < 6 ){
      short_list[c] = b;
      c++;
    }
    else{
      long_list[d] = b;
      d++;
    }
  }
  /******************************************************************/
  /******************************************************************/

  fprintf(ps_file,"Col00\n\n");
  /******************************************************************/
  /* print the short range bonds. The first bond will range from    */
  /* n - r. The next bond may correspond to an amino acid whose num */
  /* is less than r. Displaying these with equal height lines above */
  /* the flex. diagram would be confusing. So continue looking for  */
  /* a bond that binds an amino-acid greater than r                 */
  /******************************************************************/
  current_index = 0;
  increment = 5;
  while( finished != c ){

    for( b = 0; b < c; b++ ){
      
      if( short_list[b] != -1 ) {
	if( bond_list[ short_list[b] ][0] > current_index ) {
	  fprintf(ps_file,"%4d %4d moveto\n", 110+(3*bond_list[ short_list[b] ][0]), 
		  ps_line_number-3 );
	  fprintf(ps_file,"%4d %4d lineto\n", 110+(3*bond_list[ short_list[b] ][0]), 
		  ps_line_number-3-increment );
	  fprintf(ps_file,"%4d %4d lineto\n", 110+(3*bond_list[ short_list[b] ][1]), 
		  ps_line_number-3-increment );
	  fprintf(ps_file,"%4d %4d lineto stroke\n", 110+(3*bond_list[ short_list[b] ][1]), 
		  ps_line_number-3);
	  finished++;
	  current_index = bond_list[ short_list[b] ][1];
	  short_list[b] = -1;
	}
      }
    }

    increment += 5;
    current_index = 0;

  }
  /******************************************************************/
  /******************************************************************/

  fprintf(ps_file,"FlexCol30\n");
  /******************************************************************/
  /* print the long range bonds. The first bond will range from     */
  /* n - r. The next bond may correspond to an amino acid whose num */
  /* is less than r. Displaying these with equal height lines above */
  /* the flex. diagram would be confusing. So continue looking for  */
  /* a bond that binds an amino-acid greater than r                 */
  /******************************************************************/
  current_index = finished = 0;
  increment = 5;
  while( finished != d ){

    for( b = 0; b < d; b++ ){
      
      if( long_list[b] != -1 ) {
	if( bond_list[ long_list[b] ][0] > current_index ) {
	  fprintf(ps_file,"%4d %4d moveto\n", 110+(3*bond_list[ long_list[b] ][0]), 
		  ps_line_number+3 );
	  fprintf(ps_file,"%4d %4d lineto\n", 110+(3*bond_list[ long_list[b] ][0]), 
		  ps_line_number+3+increment );
	  fprintf(ps_file,"%4d %4d lineto\n", 110+(3*bond_list[ long_list[b] ][1]), 
		  ps_line_number+3+increment );
	  fprintf(ps_file,"%4d %4d lineto stroke\n", 110+(3*bond_list[ long_list[b] ][1]), 
		  ps_line_number+3);
	  finished++;
	  current_index = bond_list[ long_list[b] ][1];
	  long_list[b] = -1;
	}
      }
    }

    if( increment/5 < 10 )
      fprintf(ps_file,"FlexCol3%d\n", increment/5 );
    else
      fprintf(ps_file,"FlexCol%d\n",  increment/5 + 30 );

    increment += 5;
    current_index = 0;

  }
  /******************************************************************/
  /******************************************************************/

  /******************************************************************/
  /* Print the flexibility scale below the connections diagram.     */
  /******************************************************************/
  /***       these need to now depend on ps_line_number!!!     ******/
  /******************************************************************/
  fprintf(ps_file, "FlexCol01\n");
  fprintf(ps_file, "110 760 120 760 120 750 110 750 Pl4\n");
  fprintf(ps_file, "FlexCol02\n");
  fprintf(ps_file, "120 760 130 760 130 750 120 750 Pl4\n");
  fprintf(ps_file, "FlexCol03\n");
  fprintf(ps_file, "130 760 140 760 140 750 130 750 Pl4\n");
  fprintf(ps_file, "FlexCol04\n");
  fprintf(ps_file, "140 760 150 760 150 750 140 750 Pl4\n");
  fprintf(ps_file, "FlexCol05\n");
  fprintf(ps_file, "150 760 160 760 160 750 150 750 Pl4\n");
  fprintf(ps_file, "FlexCol06\n");
  fprintf(ps_file, "160 760 170 760 170 750 160 750 Pl4\n");
  fprintf(ps_file, "FlexCol07\n");
  fprintf(ps_file, "170 760 180 760 180 750 170 750 Pl4\n");
  fprintf(ps_file, "FlexCol08\n");
  fprintf(ps_file, "180 760 190 760 190 750 180 750 Pl4\n");
  fprintf(ps_file, "FlexCol09\n");
  fprintf(ps_file, "190 760 200 760 200 750 190 750 Pl4\n");
  fprintf(ps_file, "FlexCol10\n");
  fprintf(ps_file, "200 760 210 760 210 750 200 750 Pl4\n");
  
  fprintf(ps_file, "100 0 translate\n");
  
  fprintf(ps_file, "FlexCol11\n");
  fprintf(ps_file, "110 760 120 760 120 750 110 750 Pl4\n");
  fprintf(ps_file, "FlexCol12\n");
  fprintf(ps_file, "120 760 130 760 130 750 120 750 Pl4\n");
  fprintf(ps_file, "FlexCol13\n");
  fprintf(ps_file, "130 760 140 760 140 750 130 750 Pl4\n");
  fprintf(ps_file, "FlexCol14\n");
  fprintf(ps_file, "140 760 150 760 150 750 140 750 Pl4\n");
  fprintf(ps_file, "FlexCol15\n");
  fprintf(ps_file, "150 760 160 760 160 750 150 750 Pl4\n");
  fprintf(ps_file, "FlexCol16\n");
  fprintf(ps_file, "160 760 170 760 170 750 160 750 Pl4\n");
  fprintf(ps_file, "FlexCol17\n");
  fprintf(ps_file, "170 760 180 760 180 750 170 750 Pl4\n");
  fprintf(ps_file, "FlexCol18\n");
  fprintf(ps_file, "180 760 190 760 190 750 180 750 Pl4\n");
  fprintf(ps_file, "FlexCol19\n");
  fprintf(ps_file, "190 760 200 760 200 750 190 750 Pl4\n");
  fprintf(ps_file, "FlexCol20\n");
  fprintf(ps_file, "200 760 210 760 210 750 200 750 Pl4\n");
  /******************************************************************/
  /******************************************************************/

  return;
}
