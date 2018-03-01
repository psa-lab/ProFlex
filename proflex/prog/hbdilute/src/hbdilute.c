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
/* hbdilute.c                                        1.5.99 */
/*                                                          */
/* This program takes as input the results of FIRST_G. The  */
/* input file is named decomp_list. It contains the rigid   */
/* cluster decomposition data for the hbdilute routines of  */
/* FIRST. hbdilute.c uses this file to determine if there   */
/* was a main-chain change in the decomposition, and if so  */
/* produces a linear, 1D, description of the decomposition. */
/* Coloring of the clusters is kept consistant. If a        */
/* cluster breaks into two or more smaller clusters, the    */
/* largest of the smaller clusters is given the color of    */
/* previous cluster. The remaining new clusters are assign- */
/* ed novel colors.                                         */
/*                                                          */
/* The routines to sort the data and determine the new      */
/* colors are in hbfunctions.c. The output, which is in     */
/* postscript format, is produced by the routines in        */
/* hbpostscript.c. The data structures are defined in       */
/* types.h. The program is compiled by makefile.            */
/*                                                          */
/* revision 4.29.99                                         */
/* The program was made more modular. Moved most of the     */
/* routines that do everything to hbfunctions.c. Added lots */
/* of comment blocks to the main program. The biggest change*/
/* is that now the program is run only once. Previously,    */
/* every time the RCD changed, FIRST would run hbdilute.c   */
/* In order to keep color consistancy, there was lots of    */
/* file I/O to store the previous states of the RCD. Now,   */
/* all the data from the FIRST run is stored in decomp_list */
/* to be read by hbdilute.c at the end of the FIRST loop.   */
/* It appears to be about twice as fast, much less file I/O */
/* and a lot cleaner. The decomp_list file can potentially  */
/* get big. For a ~200 amino acid protein, it hit 1M. At the*/
/* end of the run it is deleted, but the space needs to be  */
/* there. Also, cluster.ps, the output, still overwrites any*/
/* other cluster.ps in the working directory. Will be       */
/* fixing this soon.                                        */
/*                                                          */
/* revision 5.13.99                                         */
/* Added some code, and included modified versions of the   */
/* output routines to handle multiple page landscape output */
/* The standard landscape page will hold 210 residues. If   */
/* a second or third page is needed, the output is sent to  */
/* temp_file_* files, and catted to the end of the primary  */
/* file, cluster.ps, once the pages(s) have reached the     */
/* bottom. Seems to be working. Still can't handle dimeres  */
/* or HETATMs yet though...                                 */
/*                                                          */
/* revision 6.13.99                                         */
/* The program can now handle dimers. The multipage stuff   */
/* was a hassle, but it should work for most cases. The out */
/* put is still linear, so large multisubunit proteins will */
/* require a lot of paper, as the decomposition for each    */
/* subunit will be drawn on the same line.                  */
/*                                                          */
/* addition 2.1.00                                          */
/* Added code previosly separate code to output the average */
/* flexibility profile, and to display the hbond connection */
/* s explicitly. Still working on this.                     */
/*                                                          */
/* modified output 9.17.01                                  */
/* Created a new option to output the dilution plots on a   */
/* single page, regardless of the size of the protein. The  */
/* postscript file can then be manually modified. Removed   */
/* all the misc code that had been building up, and got rid */
/* of functionality for some obsolete options.              */
/*                                                          */
/************************************************************/
/**********************************************************************/
/* addition 6.20.02 BMH                                               */
/* Switched y_translate from a single variable to an array of size 10 */
/* The program still has a max number of chains as 10, but now they   */
/* can all start with different, non 1, residue numbers.              */
/*                                                                    */
/* Added code to handle PDB files that use the "insertion code" column*/
/* of the ATOM record. Some people use insertion codes to maintain    */
/* sequence numbers according to previously published structure, for  */
/* ease in aligning the structures. hbdilute will read these in and   */
/* make appropriate spacing when printing the 1D decomp. Also, the    */
/* location of each insertion is demarked at the top of the page by an*/
/* asterisk. However, I haven't fixed the sequence numbering, and the */
/* code requires that all three backbone atoms for each inserted AA be*/
/* present in the PDB file. The multipage output does not work with   */
/* insertion codes yet.                                               */
/*                                                                    */
/* Added tick marks to the top of the output demarking the first res. */
/* and every ten residues after the first one.                        */
/*                                                                    */
/* Modified the upper bound on the number of possible clusters in the */
/* "compute_number_of_clusters" subroutine. The program runs noticably*/
/* quicker now.                                                       */
/**********************************************************************/
/**********************************************************************/
/* Additions Nov 7, 2003, Sameer Arora    			      */
/* User wants the list of flexible-region reseidue ranges in a text   */
/* file. write this info to a file here. New cmd line option has been */
/* added 't'			        			      */
/* Jan 13, 04 - Sameer Arora					      */
/* New cmd line option 'e' prints out stripe for each hbond broken    */
/* Another option 'E' to print stripe for every interval of hb_energy */
/**********************************************************************/
/**********************************************************************/
/* Additions Feb, 2008, Sandeep Namilikonda    			      */
/*
/**********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../include/types.h"
#include "../include/hbfunctions.h"
#include "../include/postscript.h"
#include "../include/connect.h"
#include "../include/graphs.h"
#include "../include/text_output.h"

#define HB_ENERGY_INTRVL 0.1

#define column(n)   (linebuf+(n)-1) /* manipulate linebuf by column */


/* #define DEBUG_SN */


/**********************************************************************/
/* The main program will read and parse the input file decomp_list,   */
/* the file produced by FIRST for the hbdilute program. The functions */
/* that do the actual parsing, and calculate the color consistancy are*/
/* in hbfunctions.c. The postscript output routines are in            */
/* in postscript.c                                                    */
/**********************************************************************/
main( int argc, char **argv ) {

    /* New additions - 2008:01	SN */
  int min_max_rsd[MAX_CHAIN_COUNT][2];
  int missing_rsd_list[MISSING_RSD_LIST_SIZE];
  int missing_rsd_count[MAX_CHAIN_COUNT];
  int insertions_prior_donor,insertions_prior_acptr,i,j,k;
  int acptr_chain_index,donor_chain_index;
  char acptr_inscode,donor_inscode;

  int
    a = 0, h = 0,
    cluster_count = 0,
    hb_number = 0,
    fit_it_all_on_one_page = 0,
    bulk_atom_label[MAXATOMS],
    ps_line_number = 0,
    ps_page_number = 0,
    atom_num = 1,
    same = 0,
    number_of_pages = 0,
    total_hbonds = 0,
    old_cluster_count = 0,
    last_line = 0,
    new_colors[MAX_CLUSTER_COUNT],
    old_colors[MAX_CLUSTER_COUNT],
    short_output = 0,
    missing_one = 0,
    long_output = 0,
    number_of_res = 0,
    y_translate[MAX_CHAIN_COUNT],
    donor_res_number = 0,
    accept_res_number = 0,
    HBtwo_start = 0, HBtwo_end = 0,
    col = 0,
    number_of_atoms = 0,
    end_of_decomp = 0,
    initial_decomp = 0,
    multiple_pages = 0,
    cluster_labels[MAX_CLUSTER_COUNT],
    cluster_sizes[MAX_CLUSTER_COUNT],
    chain_size[MAX_CHAIN_COUNT],
    number_of_chains = 0,
    side_chain_test = 0,
    single_run = 0,
    averaged_thermal = 0,
    largest_cluster_size = 0,
    largest_so_far = 0,
    largest_cluster = 0,
    total_length_of_output_line = 0,
    temp_count = 0,
    total_insertions = 0,
    insertion_res_num[100],
    number_of_insertions[100],
    dummy = 0,
    /* New Command Line options added - Sameer Arora, Jan 16 2004   	*/
    /* 't' - to write only a text file containing flexible regions only */
    /* 'e' - print every stripe for each h-bond broken			*/
    /* 'i' - to print  stripe on regular hb_energy interval		*/
    text_output_also = 0,
    output_every_stripe = 0,
    intervalled_output = 0;

  float
    aveR = 0.0,
    hb_energy = 0.0,
    hb_energy_prev = 0.0,  /* output on energy interval */
    energy_interval = 0.0,
    scale_factor = 1.0,
    mean_coordination = 0.0,
    frac_atoms_in_largest_cluster = 0.0,
    frac_flop = 0.5,
    E_array[600],
    ave_R_array[600],
    Frac_LRC_array[600],
    floppy_modes[600];

  char
    linebuf[150],
    donor_atom_type[4],
    accept_atom_type[4],
    protein_name[128],
    plots_output_file[128],
    text_output_file[128],
    *temp,
    chain_IDs[MAX_CHAIN_COUNT],
    donor_chain_ID,
    accept_chain_ID,
    insertion_chain_ID[100];

  FILE
    *ps_file,
    *ps_graphs,
    *pdb_file,
    *decomp,
    *current_file,
    *file_list[FILE_LIST_SIZE],
    *text_file,
    *test_file;

  /**********************************************************************/
  /* store the protein atom data in a linked list of type atom_list     */
  /**********************************************************************/
  atom_list
    *chem_file_head = NULL,
    *start_of_chem_file = NULL;

  /**********************************************************************/
  /* The data type cluster stores information gleaned from the FIRST    */
  /* rigid cluster decomposition data. For each rigid cluster in the    */
  /* protein (at a given Hbond conc.) this data type stores which prot. */
  /* atoms are in which cluster.                                        */
  /**********************************************************************/
  clusters
    *nc_count[MAX_CLUSTER_COUNT],
    *oc_count[MAX_CLUSTER_COUNT],
    *free_one,
    *free_next,
    *nc_current = NULL;

  /**********************************************************************/
  /* END VARIABLE DECLARATION                                           */
  /**********************************************************************/

  /**********************************************************************/
  /* INITIALIZE ARRAYS JUST IN CASE                                     */
  /**********************************************************************/
  for( h = 0; h < MAX_CLUSTER_COUNT; h++ ){
    old_colors[h] = 0;
    new_colors[h] = 0;
    oc_count[h] = NULL;
    nc_count[h] = NULL;
  }
  for( h = 0; h < FILE_LIST_SIZE; h++ ){
    file_list[h] = NULL;
  }
  for(h = 0; h < MAX_CHAIN_COUNT; h++){
    chain_size[h] = 0;
  }
  /**********************************************************************/

  /**********************************************************************/
  /* Usage check. Arguments to hbdilute are decomp_list and a single    */
  /* letter identifying the output type                                 */
  /**********************************************************************/
  if( argc <4 ){

    /* Following options either are not implemented or do not really produce a different output
       So commenting out the options' display ONLY for the meanwhile
    printf("\nUsage: hbdilute decomp_list <output type: a|b|c|e|i|m|o|s|t > <*_FIRSTdataset>\n");
    printf(" a\t- averaged_thermal\n");
    printf(" c\t- side chain test\n");
    printf(" m\t- missing one\n");
    printf(" o\t- single run\n");
    */
    printf("\nUsage: hbdilute decomp_list <output type: b|e|i|s|t > <*_proflexdataset>\n");
    printf("\noutput type options (cannot be combined) :\n");
    printf(" b\t- Fit entire plot on one page - Usually used.\n");
    printf(" e\t- Output a stripe for each h-bond broken.\n");
    printf(" i\t- Output one stripe per every energy interval change.\n");
    printf("  \t  4th (last) argument must be a +ive energy interval (e.g. 0.5).\n");
    printf(" s\t- Multi-page output.\n");
    printf(" t\t- Instead of stripes, output sequence of flexible-region residue-ranges.\n");
    printf("  \t  Output equivalent to multi-page output.\n");
    return( 0 );
  }
  /**********************************************************************/

  /**********************************************************************/
  /* The hbdilute routines included with the program FIRST allow for    */
  /* several different analysis. Depending on how FIRST was run, this   */
  /* program, hbdilute.c, will produce the appropriate output based on  */
  /* command line argument passed to this program. The following code   */
  /* sets the proper variable for producing the proper output.          */
  /**********************************************************************/
  switch( *argv[2] ){

  case 'b':
    fit_it_all_on_one_page++;
    break;
  case 's':
    short_output++;
    break;
   case 'l':
    short_output++;
    break;
  case 'm':
    missing_one++;
    break;
  case 'c':
    side_chain_test++;
    break;
  case 'o':
    single_run++;
    break;
  case 'a':
    averaged_thermal++;
    break;


    /* for text output along with short output*/
  case 't':
    text_output_also++;
    short_output++;
    break;
   /* print every stripe for each h-bond broken */
  case 'e':
    output_every_stripe++;
    break;
  case 'i':
    intervalled_output++;
    if(argc < 5) {
	printf("Enter an energy interval as the 4th argument\n");
	exit(-1);
    }
    energy_interval = (float)atof(argv[4]);
    if( energy_interval < 0 ){
	printf("Enter a positive energy interval!\n");
	exit(-1);
    }
    break;
  }
  /**********************************************************************/

  /**********************************************************************/
  /* Open the appropriate files to get and put data                     */
  /* decomp - the rigid cluster decompostion data for every H-bond      */
  /* pdb_file - the *FIRST.chem file that is output from FIRST. Contains*/
  /*            all the atomic information of the protein.              */
  /**********************************************************************/
  decomp   = fopen( argv[1], "r" );
  pdb_file = fopen( argv[3], "r" );

  /**********************************************************************/
  /* Initialize the page and line number info for the postscript        */
  /* output file. Get the number of residues and number of atoms        */
  /* from the decomp_list header. Print the header of the output        */
  /* cluster.ps. Read in the atom information of chem file.             */
  /**********************************************************************/
  fgets( linebuf, sizeof(linebuf), decomp );

  number_of_res   = atoi( column(8) );
  number_of_atoms = atoi( column(16) );
  sscanf( column(25) , "%100s", protein_name );

  temp =  strrchr( protein_name, '_' );
  if( temp )
    *temp = '\0';

  memcpy( plots_output_file, protein_name, sizeof(protein_name) );
  strcat( plots_output_file, "_plots.ps");
  strcat( protein_name,".ps");

  rewind( decomp );

    /* flexible region output - Sameer */
  if( text_output_also) {
      strcpy(text_output_file, protein_name);
      strcat(text_output_file,"_flexible.txt");
     if(open_output_file(text_output_file, &text_file) != 1)
      	return(0);
     else
      	 fprintf( text_file, "HBond# HB-Enrgy <r> \t[ \t\tFlexible Residue Ranges \t\t] \n" );

     }

  if( !fopen(protein_name,"r") ){
    ps_file   = fopen( protein_name, "w" );
  }
  else{
    printf("\n\t There is already a copy of %s in this directory.\n", protein_name );
    printf("\t The output will be saved in tmp.ps\n" );

    if( !fopen("tmp.ps", "r") ){
      ps_file   = fopen("tmp.ps" , "w" );
    }
    else{
      printf("\n\t There is also a tmp.ps in this directory. The program will\n");
      printf("\t terminate. After you rename or remove %s and tmp.ps, you\n", protein_name);
      printf("\t may run the program from the command line as follows:\n");
      printf("\t hbdilute decomp_list s *_proflexdataset\n\n");
      return(0);
    }
  }

  ps_page_number = 1;
  initial_decomp = -1;

  file_list[0] = ps_file;

  /**********************************************************************/
  /**********************************************************************/

  /**********************************************************************/
  /* reads the atom data from the *_FIRST.chem file. This file          */
  /* conforms to standard PDB file format, although the atom            */
  /* number line has been set by FIRST.                                 */
  /**********************************************************************/
  read_chem_file( &chem_file_head, pdb_file, y_translate, chain_size, chain_IDs,
		  &number_of_chains, &total_insertions, insertion_res_num,
		  number_of_insertions, insertion_chain_ID,
		  /* SN 2008:01 */ 
		  min_max_rsd,missing_rsd_list,missing_rsd_count,&total_hbonds);

  start_of_chem_file = chem_file_head;
  /**********************************************************************/

  /**********************************************************************/
  /* Calculate appropriate scale factor to fit ouput on a single page   */
  /* If the total space needed for a single line of output is < 770     */
  /* do not scale the image.                                            */
  /**********************************************************************/
  if( fit_it_all_on_one_page ) {
    for( a = 0; a <= number_of_chains; a++ )
      total_length_of_output_line += (chain_size[a]*3);
    total_length_of_output_line += (number_of_chains-1)*18 + 160;

    if( total_length_of_output_line > 740 )
      scale_factor = 740.0 / total_length_of_output_line;
  }
  /**********************************************************************/

  print_header(ps_file, scale_factor );
  print_data_headings(ps_file, long_output, number_of_res,
		      number_of_chains, chain_size );
  /**********************************************************************/

#ifdef DEBUG_SN
  printf("\nDEBUG (hbdilute.c):  number of residues: %d and chains: %d\n",\
		number_of_res, number_of_chains );
#endif

  /**********************************************************************/
  /* set the initial page and line number for the output files          */
  /**********************************************************************/
  if( number_of_chains )
    number_of_pages = ( (number_of_res+1) + ( (number_of_chains)*6 ) ) / 200 ;
  else
    number_of_pages = ( number_of_res+1 ) / 200;

  
  ps_line_number = 60;
  if( number_of_pages > 0)
    multiple_pages = 1;

  /* 2008:02 	Disable portrait format plot generation 

	  if( number_of_pages > 0 ){
	    multiple_pages = 1;
	    ps_line_number = 60;
	  }
	  else{
	    if( number_of_chains ){
	      ps_line_number = 60;
	    }
	    else{
	      if( number_of_res <= 125 )
		ps_line_number = 710;
	      else
		ps_line_number = 60;
	    }
	  }
  */

  /**********************************************************************/
  /**********************************************************************/

  /**********************************************************************/
  /* If printing all data on a single page of output, suppress routines */
  /* developed to deal with multiple pages                              */
  /**********************************************************************/
  if( fit_it_all_on_one_page )
    multiple_pages = number_of_pages = 0;
  /**********************************************************************/

  /**********************************************************************/
  /* The following code gets the cluster info for                       */
  /* the nth H-bond that caused a change in the RCD                     */
  /* The macro "column(n)" is #defined above                            */
  /**********************************************************************/
  while( fgets( linebuf, sizeof(linebuf), decomp ) != NULL ){

    /**************************************************/
    /* in the decomp_list input file, new decomp info */
    /* begins with the letter A in the first column.  */
    /* The data is terminated by an END in the row.   */
    /**************************************************/
    if( !strncmp(linebuf,"A", 1) ) {

      hb_number     = atoi( column(6)  );
      hb_energy     = atof( column(13) );
      donor_res_number = atoi( column(31) );
      sscanf( linebuf+36 ,"%3s", donor_atom_type );
      donor_chain_ID  = *column(44);
      accept_res_number = atoi( column(49) );
      sscanf( linebuf+53,"%3s", accept_atom_type );
      accept_chain_ID = *column(60);
      mean_coordination = atof( column(65) );
      largest_cluster_size = 0;
      atom_num = 1; /* set to 1 so we can index the bulk_atom_label array by the */

    }               /* chem_file atom_num which does not start from zero.        */

    if( !strncmp(linebuf,"B", 1) ) {

      aveR = atof( column(5) );
      frac_atoms_in_largest_cluster = atof( column(15) );
      frac_flop = atof( column(28) );
    }

    if( strncmp(linebuf,"A", 1) &&
	strncmp(linebuf,"B", 1) &&
	strncmp(linebuf,"HEADER",6) ) {
      
      end_of_decomp = !strncmp( linebuf, "END", 3 );
      
      /************************************************************/
      /* read in the bulk atom label data from the decomp_list    */
      /* This label assigns a unique number to each atom in the   */
      /* according to which cluster it belongs to. The clusters   */
      /* labeled from largest (1) to smallest (n).                */
      /************************************************************/
      if( !end_of_decomp ){
	for( col = 0; col < 20; col++ ) {
	  bulk_atom_label[atom_num] = atoi( linebuf+(col*5) );
	  if( atoi( linebuf+(col*5) ) == 0 )
	    break;
	  
	  if( bulk_atom_label[atom_num] == 1 )
	    largest_cluster_size++;
	  
	  atom_num++;
	}
	
      }
      
      /************************************************************/
      /* execute the following code once all the data for the     */
      /* current decomp has been read in. The data was terminated */
      /* by an END line in the decomp_list file.                  */
      /************************************************************/
      if( end_of_decomp ){
	
	/*
	  E_array[hb_number]        = hb_energy;
	  ave_R_array[hb_number]    = mean_coordination;
	  Frac_LRC_array[hb_number] = frac_atoms_in_largest_cluster;
	  floppy_modes[hb_number]   = frac_flop;
	*/
	
	/***********************************************************/
	/* The next routine computes the number of "clusters" in   */
	/* the protein. As very small clusters are unintersting, a */
	/* cutoff has been set in the routine. Only clusters with  */
	/* at least  N  mainchain atoms are counted. This also can */
	/* be very small, but it can be scaled to any desired value*/
	/***********************************************************/
	cluster_count = compute_number_of_clusters( bulk_atom_label,
						    chem_file_head,
						    cluster_labels,
						    number_of_atoms);
	
  	/************************************************************/
	/* If the decomp that is currently in memory from the decomp*/
	/* _list file is not the first one, this code gets the new  */
	/* info into linked lists of struct type: clusters. These   */
	/* lists are compared to the old_lists to determine which   */
	/* of the new_clusters are subsets of the old clusters, and */
	/* to color the 1D decomp accordingly.                      */
	/************************************************************/
	if( !initial_decomp ){
	  
	  /**********************************************************/
	  /* reset the chem_file list                               */
	  /**********************************************************/
	  chem_file_head = start_of_chem_file;
	  
	  /**********************************************************/
	  /* use the chem_file list to store                        */
	  /* the new cluster data in nc_count                       */
	  /**********************************************************/
	  set_new_decomp_info( cluster_count, bulk_atom_label, atom_num,
			       nc_count, chem_file_head, cluster_labels );
	  
	  /**********************************************************/
	  /* If there is no change in the                           */
	  /* main-chain decomposition, do                           */
	  /* not execute the remaining code                         */
	  /* compare_old_and_new_decomps returns 1 if same else 0   */
	  /**********************************************************/
	  same = compare_old_and_new_decomps( nc_count, oc_count, cluster_count,
					      cluster_sizes, old_cluster_count );
	  
	  /*
	    if( hb_number == 540 || hb_number == 539 ){
	    for( a = 0; a < cluster_count; a++ ){
	    
	    nc_current = nc_count[a];
	    
	    while( nc_current ){
	    printf("%6d %3d %c\n", nc_current->atom_number, nc_current->residue_number, nc_current->chain_ID );
	    nc_current = nc_current->next_element;
	    }
	    }
	    
	    }
	  */
	  /**********************************************************/
	  /* If the mainchain decomp did not change                 */
       	  /* reset the cluster list  nc_count for                   */
	  /* the next decomp to be read.                            */
	  /**********************************************************/
	  if( last_line )
	    {
	      continue;
	    }
	  
	  /************************************************************/
	  /* Setting same = 0 will print the decomp results for every */
	  /* hydrogen bonds, even if the result is redundant with the */
	  /* previous line.                                           */
	  /************************************************************/
	  /*same = 0;*/
	  
	  /* Print out every stripe for each bond broken.
	     This is now a command line option 'e' */
	  if( intervalled_output || output_every_stripe )
	    same = 0;
	  
	  if( same ){
	    if( cluster_count == 0 ){
	      last_line = 1;
	    }
	    else{
	      for( a = 0; a < MAX_CLUSTER_COUNT; a++ ){
		
		free_one = free_next = nc_count[a];
		
		if( free_one != NULL ){
		  
		  while( free_next != NULL ){
		    free_next = free_one->next_element;
		    free ( free_one );
		    free_one = free_next;
		  }
		  if(free_one != NULL)	/* 2006:08:02	SN  Possible "SEG FAULT" source */
  		    free ( free_one );
		}
		nc_count[a] = NULL;
		
	      }
	      same = 0;
	      
	    }
	    
	    continue;
	  }
	  
	  /*****************************************/
	  /* This checks to see which of the new   */
	  /* clusters are subsets of the old       */
	  /* clusters, and colors them accordingly */
	  /*****************************************/
	  compute_new_cluster_list(nc_count, cluster_count,
				   oc_count, old_cluster_count,
				   old_colors, new_colors );


	  /*****************************************/
	  /* if intervalled_output !=0, print  */
	  /* a stripe every HB_ENERGY_INTRVL of	   */
	  /* hb_energy                             */
	  /*****************************************/
	  /*if( intervalled_output ){
	    //printf("Prev = %6.2f  New = %6.2f\n", hb_energy_prev, hb_energy);
            if ( fabs( hb_energy_prev - hb_energy) >= HB_ENERGY_INTRVL)
      		hb_energy_prev = hb_energy;
	    else
		continue;
	  }*/

	  /*****************************************/
	  /* User wants the list of flexible-region*/
	  /* reseidue ranges in a text file.	   */
	  /* write this info to a file here	   */
	  /*****************************************/

	  if(text_output_also)
	    print_text_decomp( nc_count, cluster_count, ps_line_number,\
			       text_file,y_translate, number_of_res, \
			       hb_number, hb_energy, mean_coordination);
	  
	  
	  /* SN 2008:01	Account for inserted residues */
	  /* 		
	   * When an inserted residue is involved in an H-bond,
	   * only the insertions that precede that residue have
	   * to be accounted for!
	   */
	  insertions_prior_donor = 0;
	  insertions_prior_acptr = 0;
	  donor_inscode = ' ';
	  acptr_inscode = ' ';
	  for(i=0;i<=number_of_chains;i++)
	    {
	      if(chain_IDs[i] == accept_chain_ID)
		acptr_chain_index = i;
	      if(chain_IDs[i] == donor_chain_ID)
		donor_chain_index = i;
	    }
	  
	  /* SN 2008:01	Also, account for the missing residues! */
	  if(missing_rsd_count[donor_chain_index] > 0)
	    {
	      j = 0;
	      for ( i = 0; i < donor_chain_index; i++)
		j += missing_rsd_count[i];
#ifdef DEBUG_HBD
	      printf("\nDonor res #: %d",donor_res_number);
#endif
	      
	      for ( i = j; i < j+missing_rsd_count[donor_chain_index]; i++)
		{
		  /* Compute how many insertions are present before the first 
		   * missing rsd. This is necessary as missing_rsd_list elements
		   * do not account for insertions whereas donor_res_number
		   * and accept_res_number (the inputs from decomp_list file)
		   * are indices into the actual rsd list such that the begining
		   * of a chain and missing rsds are not accounted for. So, when
		   * we compare the missing_rsd_list elements, we have to 
		   * take into account the inserted residue count and chain start!
		   */
		  
		  insertions_prior_donor = 0;
		  
		  for( k = 0; k < total_insertions; k++ )
		    {
		      if( insertion_res_num[k] < missing_rsd_list[i] && \
			  donor_chain_ID == insertion_chain_ID[k] )
			{			
			  insertions_prior_donor += number_of_insertions[k];
			}
		    }
		  
		  if( (missing_rsd_list[i]+insertions_prior_donor) <= \
		      (donor_res_number + y_translate[donor_chain_index]-1 ))
		    {
		      donor_res_number++;
		    }
		  else
		    break;
		}
#ifdef DEBUG_HBD
	      printf("\tDonor res #: %d",donor_res_number);
#endif
	    }
	  if(missing_rsd_count[acptr_chain_index] > 0)
	    {
	      j = 0;
	      for ( i = 0; i < acptr_chain_index; i++)
		j += missing_rsd_count[i];
#ifdef DEBUG_HBD
	      printf("\nAcptr res #: %d",accept_res_number);
#endif
	      
	      for ( i = j; i < j+missing_rsd_count[acptr_chain_index]; i++)
		{
		  insertions_prior_acptr = 0;
		  for( k = 0; k < total_insertions; k++ )
		    {
		      if( insertion_res_num[k] < missing_rsd_list[i] && \
			  accept_chain_ID == insertion_chain_ID[k] )
			{			
			  insertions_prior_acptr += number_of_insertions[k];
			}
		    }
		  if( (missing_rsd_list[i]+insertions_prior_acptr) <=\
		      (accept_res_number +y_translate[acptr_chain_index]-1))
		    {
		      accept_res_number++;
		    }
		  else
		    break;		
		}
#ifdef DEBUG_HBD
	      printf("\tAcptr res #: %d\n",accept_res_number);
#endif
	    }
	  
	  
	  insertions_prior_donor = 0;
	  insertions_prior_acptr = 0;
	  for(i = 0; i < total_insertions; i++)
	    {
	      /* Check if the inserted residue precedes donor rsd */
	      /* Now, check if the donor rsd itself is an inserted rsd */
	      /* If yes then only add those insertions that precede    */
	      /* the actual inserted donor residue			   */
	      if( ((donor_res_number + y_translate[donor_chain_index] - 1\
		    - insertions_prior_donor) >= insertion_res_num[i]) \
		  && insertion_chain_ID[i] == donor_chain_ID )
		{
		  if ((donor_res_number + y_translate[donor_chain_index] - 1 \
		       - insertions_prior_donor) <= (insertion_res_num[i] \
						     + number_of_insertions[i]))
		    {
		      j = (donor_res_number + y_translate[donor_chain_index]-1) \
			- insertions_prior_donor - insertion_res_num[i];
		      
		      insertions_prior_donor += j;
		      
		      if( j > 0 )
			donor_inscode = 'A'+j-1;
		    } 
		  else
		    insertions_prior_donor += number_of_insertions[i];
		}
	      
	      if( ((accept_res_number + y_translate[acptr_chain_index] - 1\
		    - insertions_prior_acptr) >= insertion_res_num[i]) && \
		  insertion_chain_ID[i] == accept_chain_ID )
		{
		  
		  if((accept_res_number + y_translate[acptr_chain_index] - 1 \
		      - insertions_prior_acptr) <= (insertion_res_num[i] \
						    + number_of_insertions[i]))
		    {
		      j = (accept_res_number + y_translate[acptr_chain_index] - 1) - \
			insertions_prior_acptr - insertion_res_num[i];
		      
		      insertions_prior_acptr += j;
		      if( j > 0 ) 
			acptr_inscode = 'A'+j-1;
		    } 
		  else
		    insertions_prior_acptr += number_of_insertions[i];	
		}
	      
	    } /* End for */
	  
#ifdef DEBUG_HBD
	  
	  printf("%d %c",(donor_res_number + y_translate[donor_chain_index] - 1\
			  - insertions_prior_donor),donor_inscode);
	  printf("\ndonor_res_number: %d",\
		 (donor_res_number+ y_translate[acptr_chain_index] - 1));
	  printf("\tinsertions_prior_donor: %d\n", insertions_prior_donor);
	  
	  printf("%d %c",(accept_res_number + y_translate[acptr_chain_index] - 1\
			  - insertions_prior_acptr),acptr_inscode);
	  printf("\naccept_res_number: %d",\
		 (accept_res_number+ y_translate[acptr_chain_index] - 1));
	  printf("\tinsertions_prior_acptr: %d\n", insertions_prior_acptr);
	  
#endif
	  
	  


	  /*****************************************/
	  /* The following routines print the      */
	  /* decomp information to the output file */
	  /* cluster.ps in postscript format.      */
	  /* Sameer - Jan 15 2004                  */
	  /* Print only if user hasn't rqstd energy*/
	  /* intervalled-output OR if he has, the  */
	  /* interval is largen than threshold     */
	  /*****************************************/
	  
	  if( !intervalled_output || \
	      ( fabs( hb_energy_prev - hb_energy) >= energy_interval))
	    {
	      if( multiple_pages ){
		
		if(ps_line_number == 60){
		  set_temp_files( file_list, &number_of_pages, &ps_page_number );
		}
		
		print_multipage_decomp( nc_count,cluster_count,ps_line_number,\
		 ps_file,new_colors, y_translate, number_of_res, file_list,\
		 number_of_pages, number_of_chains, chain_size,\
		 chain_IDs, total_insertions, insertion_res_num,\
		 number_of_insertions, insertion_chain_ID,\
		  /* SN 2008:01 */ 
		 min_max_rsd, missing_rsd_list, missing_rsd_count );

		print_current_Hbond_info_landscape_multipage(
			ps_file,long_output,
			short_output, initial_decomp,&ps_line_number, hb_number,
			hb_energy, cluster_count,new_colors, donor_res_number,
			accept_res_number,number_of_res, scale_factor,
			&ps_page_number,missing_one, HBtwo_start, HBtwo_end,
			y_translate,donor_atom_type, accept_atom_type,
			file_list,number_of_pages, number_of_chains,
			donor_chain_ID, accept_chain_ID, chain_IDs,
			protein_name, chain_size, mean_coordination,
			/* SN 2008:01 */
			insertions_prior_donor, insertions_prior_acptr,
			donor_inscode,acptr_inscode,total_hbonds );
	      }
	      else{
		
		print_decomp( nc_count, cluster_count, ps_line_number, ps_file,
			      new_colors, y_translate, number_of_res, number_of_chains,
			      chain_size, chain_IDs, total_insertions, insertion_res_num,
			      number_of_insertions, insertion_chain_ID,
			      /* SN 2008:01 */ 
			      min_max_rsd, missing_rsd_list, missing_rsd_count );
		
		
		
		/*	2008:01:23 SN	Disable portrait output format
			
		if( !number_of_chains && (  number_of_res <= 125 ) )
		print_current_Hbond_info_portrait(ps_file, long_output, \
		short_output, initial_decomp,&ps_line_number, \
		hb_number, hb_energy, cluster_count,new_colors,\
		donor_res_number, accept_res_number,number_of_res,\
		scale_factor, &ps_page_number,missing_one, \
		HBtwo_start, HBtwo_end, y_translate,donor_atom_type,\
		accept_atom_type, number_of_chains,donor_chain_ID,\
		accept_chain_ID, chain_IDs, protein_name,\
		side_chain_test, mean_coordination, fit_it_all_on_one_page);
		else		
		*/

		print_current_Hbond_info_landscape(
			 ps_file, long_output, short_output, initial_decomp,
			 &ps_line_number, hb_number, hb_energy, cluster_count,
			 new_colors, donor_res_number, accept_res_number,
			 number_of_res, scale_factor, &ps_page_number,
			 missing_one, HBtwo_start, HBtwo_end, y_translate,
			 donor_atom_type, accept_atom_type, number_of_chains,
			 donor_chain_ID, accept_chain_ID, chain_IDs, chain_size,
			 protein_name, side_chain_test, mean_coordination,
			 fit_it_all_on_one_page,
			 /* SN 2008:01 */
			 insertions_prior_donor, insertions_prior_acptr,
			 donor_inscode,acptr_inscode,total_hbonds );	    
	      }
	      
	      //	  printf("\nPrev = %6.2f  New = %6.2f\n", hb_energy_prev, hb_energy);
	      
	      hb_energy_prev = hb_energy;
	    }
	  /*****************************************/
	  /*****************************************/
	  
	  
	  
	  /***********************************************/
	  /* reset the new cluster list. Set the current */
	  /* cluster list to the old cluster reference   */
	  /* reset variables.                            */
	  /***********************************************/
	  for( a = 0; a < MAX_CLUSTER_COUNT; a++ ){
	    oc_count[a]   = nc_count[a];
	    old_colors[a] = new_colors[a];
	    nc_count[a]   = NULL;
	    old_cluster_count = cluster_count;
	  }
	  
	}
	/************************************************************/
	/***********  End of the not initial decomp loop  ***********/
	/************************************************************/
	
	
	/************************************************************/
	/* This loop executes once, for the initial decomposition of*/
	/* the protein. All the H-bonds are present in this decomp  */
	/* and it is used as the basis for the coloring scheme.     */
	/************************************************************/
	if( initial_decomp ){

	  for( a = 0; a < cluster_count; a++ ){
	    if( cluster_sizes[a] > largest_so_far ){
	      largest_so_far = cluster_sizes[a];
	      largest_cluster = a;
	    }
	  }

	  /**********************************************************************/
	  /* I really don't know why these next four lines are here. BMH 6.20.02*/
	  /*
	  free_one = nc_count[largest_cluster];
	  while( free_one != NULL )
	  free_one = free_one->next_element;
	  free( free_one );
	  */
	  /**********************************************************************/

	  first_decomp( long_output, cluster_count, bulk_atom_label,
			chem_file_head, oc_count, old_colors, cluster_labels );

      	  /*****************************************/
      	  /* User wants the list of flexible-region*/
      	  /* residue ranges in a text file, make  */
      	  /* write to that file here        	   */
      	  /*****************************************/
      	    if(text_output_also)
	      {
		printf("\n**************************************************************");
		printf("\n* WARNING: If input protein has multiple chains then it is   *"); 
		printf("\n*          likely that the residue ranges of the chains are  *");
		printf("\n*          overlayed! e.g., chain 1: 11-100 (90 residues),   *");
		printf("\n*                           chain 2: 2-20 (19 residues), then*");
		printf("\n*          the final displayed residue range will be: 2-109  *");
		printf("\n**************************************************************\n");
      		print_text_decomp( nc_count, cluster_count, ps_line_number,
				   text_file,y_translate, number_of_res, hb_number, 
				   hb_energy, mean_coordination);
	      }
	    
	    
	    /* SN 2008:01	Account for inserted residues */
	    /* 		
	     * When an inserted residue is involved in an H-bond,
	     * only the insertions that precede that residue have
	     * to be accounted for!
	     */
	    insertions_prior_donor = 0;
	    insertions_prior_acptr = 0;
	    donor_inscode = ' ';
	    acptr_inscode = ' ';
	    
	    
	    /* if more then one page of output */
	    if( multiple_pages ){
	      
	      set_temp_files( file_list, &number_of_pages, &ps_page_number );
	      
	      print_multipage_decomp( 
			oc_count, cluster_count, ps_line_number, ps_file,
			old_colors, y_translate, number_of_res, file_list,
			number_of_pages, number_of_chains, chain_size,
			chain_IDs, total_insertions, insertion_res_num,
			number_of_insertions, insertion_chain_ID,
			/* SN 2008:01 */
			min_max_rsd, missing_rsd_list, missing_rsd_count );
	      
	      print_current_Hbond_info_landscape_multipage(
			ps_file, long_output, short_output, initial_decomp,
			&ps_line_number, hb_number, hb_energy, cluster_count,
			new_colors, donor_res_number, accept_res_number,
			number_of_res, scale_factor, &ps_page_number,
			missing_one, HBtwo_start, HBtwo_end, y_translate,
			donor_atom_type, accept_atom_type, file_list,
			number_of_pages, number_of_chains,
			donor_chain_ID, accept_chain_ID, chain_IDs,
			protein_name, chain_size, mean_coordination,
			/* SN 2008:01 */
			insertions_prior_donor, insertions_prior_acptr,
			donor_inscode,acptr_inscode,total_hbonds );
	    }
	    /********************************************************************************/
	    /* if only one page is needed to display a single linear decompostion, print it */
	    /* in portrait or landscape, depending on the size.                             */
	    /********************************************************************************/
	    else{
	      
	      print_decomp( oc_count, cluster_count, ps_line_number, ps_file,\
			    old_colors, y_translate, number_of_res, number_of_chains,\
			    chain_size, chain_IDs, total_insertions, insertion_res_num,\
			    number_of_insertions, insertion_chain_ID,\
			    /* SN 2008:01 */
			    min_max_rsd, missing_rsd_list, missing_rsd_count );
	      
	      
	      /*	2008:01:23 SN	Disable portrait output format
			
	      if( !number_of_chains && ( number_of_res <= 125 ) )
	      print_current_Hbond_info_portrait(
	      ps_file, long_output, short_output, initial_decomp,
	      &ps_line_number, hb_number, hb_energy, cluster_count,
	      new_colors, donor_res_number, accept_res_number,
	      number_of_res, scale_factor, &ps_page_number,
	      missing_one, HBtwo_start, HBtwo_end, y_translate,
	      donor_atom_type, accept_atom_type, number_of_chains,
	      donor_chain_ID, accept_chain_ID, chain_IDs, 
	      protein_name, side_chain_test, mean_coordination, 
	      fit_it_all_on_one_page );
	      else
	      {
	      if( !multiple_pages )
	      */
	      print_current_Hbond_info_landscape(
			ps_file, long_output, short_output, initial_decomp,
			&ps_line_number, hb_number, hb_energy, cluster_count,
			new_colors, donor_res_number, accept_res_number,
			number_of_res, scale_factor, &ps_page_number,
			missing_one, HBtwo_start, HBtwo_end, y_translate,
			donor_atom_type, accept_atom_type, number_of_chains,
			donor_chain_ID, accept_chain_ID, chain_IDs, chain_size,
			protein_name, side_chain_test, mean_coordination,
			fit_it_all_on_one_page,
			/* SN 2008:01 */
			insertions_prior_donor, insertions_prior_acptr, 
			donor_inscode,acptr_inscode,total_hbonds );
	      /*}*/
	    }
	    
	    /**********************************************************/
	    /* This option creates the linear flexibility diagram and */
	    /* displays the hydrogen bond connections as well.        */
	    /**********************************************************/
	    if( single_run ) {
	      rewind( pdb_file );
	      print_connections( pdb_file, ps_file, y_translate );
	    }
	    
	    old_cluster_count = cluster_count;
	    initial_decomp = 0;
	}
	/********************* END OF if(initial_decomp) ************/
	/************************************************************/
	
      } /* End of if( end_decomp ) */
      
    } /* End of 2.73333 */

  } /* End of while ( fgets ... ) */


  /**********************************************************************/
  /* the next subroutine will create 2D plots of the input data. The    */
  /* output file will have the same prefix as the dilution plot, with   */
  /* the suffix "_plots.ps".                                            */
  /**********************************************************************/
  /*
  ps_graphs = fopen( plots_output_file, "w" );
  plot_data( ps_graphs, E_array, ave_R_array, Frac_LRC_array, floppy_modes, total_hbonds );
  */
  /**********************************************************************/

  /**********************************************************************/
  /* small chuck of code to make the post script file valid and delete  */
  /* temp files    UPDATE THE COUNT TO < 15 (MAX # of files  -- Sandeep */
  /**********************************************************************/
  for( a = 1; a < FILE_LIST_SIZE; a++ ){
    current_file = file_list[a];
    if( current_file != NULL ){
      fflush( current_file );
      fclose( current_file );
    }
  }


	   /*	2008:01:23 SN	Disable portrait output format

  if( (  number_of_chains && ( (number_of_res*(number_of_chains*6) ) <= 125 ) ) ||
      ( !number_of_chains && ( number_of_res <= 125 ) ) )
    print_portrait_footer( ps_file, protein_name, ps_line_number );
  else
	  */

    print_landscape_footer( ps_file, protein_name, ps_line_number,scale_factor );

  clean_up( ps_file, number_of_pages );
  /**********************************************************************/

  fclose( ps_file );
  fclose( decomp );
  fclose( pdb_file );
  /*fclose( ps_graphs );*/

   if(text_file)
        fclose(text_file);

  return(0);
}
