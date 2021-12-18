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

/************************************************************/
/* File:   text_output.h                                    */
/* Author: Sameer Arora                                     */
/*                                                          */
/* Created on November 10, 2003, 5:23 PM                    */
/* This file defines the functions used for generating      */
/* text-only output file. The output has only residue ranges*/
/* of flexible ranges along with the margin's info which    */
/* appears in the regular postscript output                 */
/************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "../include/types.h"


/**********************************************************************/
/* This routine is invoked only when user chooses text-only output.   */
/* Output is either written as Postscript or Text - not both.         */
/* detailing the flexible residue ranges only with each h-bond break. */
/* Outout is written to a text file instead of generating a postscript*/                
/**********************************************************************/

static void swap_colors(int array[MAX_CLUSTER_COUNT], int A, int B) {
  int temp;
  temp = array[A];
  array[A] = array[B];
  array[B] = temp;
}

static void quick_sort(int *array, int start, int end) {

  int a = 0, last = 0;

  if (start >= end)
    return;
  swap_colors(array, start, (start + end) / 2);
  last = start;
  for (a = start + 1; a <= end; a++) {
    if (array[a] < array[start])
      swap_colors(array, ++start, a);
  }
  swap_colors(array, start, last);
  quick_sort(array, start, last - 1);
  quick_sort(array, last + 1, end);
}

void print_text_decomp( 
                    clusters *nc_count[120], int cluster_counter, int line_number, 
                    FILE *text_file,  int y_translate[10], int number_of_residues,     
                    int hb_number, float hb_energy, float mean_coordination)
{

int 
    rigid_regions_index = 0, 
    rigid_counter = 0,
    a = 0;                                              

int * rigid_regions = NULL;

clusters  
    *nc_current = NULL;

static int 
    first_res_no,static_counter;


    /* print the Hydrogen - bond information on the left */
    //printf("\n----%d th run ------\n",++static_counter);
    fprintf( text_file, "%4d %7.3f %7.3f\t[  ", hb_number, hb_energy,mean_coordination );

    /* array rigid would store the residue numbers of residues
    which are part of any rigid cluster */

    /* thrice becuase for non-terminal residues, N, Ca, CB can appear */

    rigid_regions = malloc( 3*number_of_residues * sizeof(int));  

    /* initialize the residue numbers */
    for( 	rigid_regions_index = 0;
		rigid_regions_index <number_of_residues*3; 
		rigid_regions_index++)
    {
        rigid_regions[rigid_regions_index] = -999;
    }

    rigid_regions_index = 0;    

    if( line_number == 710 )
      first_res_no = y_translate[0];	/* First chain */
	
    for( a = 0; a < cluster_counter; a++ ){
        nc_current = nc_count[a];
        while( nc_current ) {
              rigid_regions[rigid_regions_index] = nc_current->residue_number;
              rigid_regions_index++;
              nc_current = nc_current->next_element;
         }
    }
  
  
    rigid_counter = rigid_regions_index;
    rigid_regions[rigid_counter] = first_res_no-1;

    rigid_counter++;
    rigid_regions[rigid_counter] = first_res_no + number_of_residues;
    fprintf(text_file," %d to %d , %d \n ", first_res_no , 
		first_res_no + number_of_residues-1,rigid_counter );
     
    
    quick_sort( rigid_regions, 0, rigid_counter );

    //for( a = 0;a <rigid_counter+1 ; a++)
	   //printf("%d,", rigid_regions[a]);
    
     for( a = 0;a <rigid_counter; a++){

       if( rigid_regions[a] < rigid_regions[a+1]-1)
	   //printf(" (%d--%d) ", rigid_regions[a]+1, rigid_regions[a+1]-1);
           fprintf(text_file," (%d--%d) ", rigid_regions[a]+1, rigid_regions[a+1]-1);
     }
	 
	 
    // print information onto the right margin
    
    fprintf( text_file, "  ]\n");
  free(rigid_regions);
                       
}
                       
                       
