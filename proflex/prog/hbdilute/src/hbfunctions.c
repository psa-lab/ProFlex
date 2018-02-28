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

/**********************************************************************/
/* These are the main functions that read in the decomp data, sort it */
/* and produce the color consistency in the decomposition output.     */
/*                                                                    */
/* BMH 6.19.02 Added code in read_chem_file to read insertion codes.  */
/* The logic isn't perfect, it requires that all three backbone atoms */
/* be present in the file (N, CA, C).                                 */
/**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "../include/types.h"
#include "../include/postscript.h"
#include "../include/text_output.h"

#define column(n)   (linebuf+(n)-1) /* manipulate linebuf by column */

#define DEBUG_SN

/***********************************************************************/
/* This routine reads the data for the main-chain atoms in the protein */
/* from the *_FIRSTdataset file and stores them in a linked list named */
/* chem_file_head. A pointer to the beginning of the list is passed    */
/* back to hbdilute.c                                                  */
/***********************************************************************/
void read_chem_file( atom_list **chem_file_head, FILE *pdb_file, int y_translate[10],
		     int chain_size[10], char chain_IDs[10], int *number_of_chains, 
		     int *total_insertions, int insertion_res_num[100], 
		     int number_of_insertions[100], char insertion_chain_ID[100],
		     /* SN 2008:01 */ 
		     	int min_max_rsd[10][2], 
			int missing_rsd_list[MISSING_RSD_LIST_SIZE],
			int missing_rsd_count[10],
			int *total_hbond_count) {

  int insertion_count[10];
  int i,prev_rsd=0,j;
			/* SN 2007:12:19 
			   - Keep track of min and max residues in each chain
			     so that no residues are missing
			   - Keep track of the total insertions in each chain
			   - Compare consecutive rsd numbers to track any 
			     missing rsds
			 */ 
  int  
    a = 0, 
    is_atom = 0,
    known_insertion = 0;

  char 
    linebuf[150], 
    atom[4], 
    chain_label;

  atom_list *temp = NULL,
            *current = NULL;

  for( a = 0; a < 100; a++)
    number_of_insertions[a] = 0;

  for( a = 0; a < MISSING_RSD_LIST_SIZE; a++ ) 
    missing_rsd_list[a] = 0; 		/* SN 2008:01 */
  

  /* SN	2007:12:19 */
  for( a = 0; a < 10; a++ ) {
    min_max_rsd[a][0] = MAX_RSD;	/* SN 2008:02 Defined in types.h */
    min_max_rsd[a][1] = 0;
    insertion_count[a] = 0;
    missing_rsd_count[a] = 0;
  }

  *total_hbond_count = 0;		/* SN 2008:05 for hbdilute plot */

  while( fgets( linebuf, sizeof(linebuf), pdb_file ) != NULL ) {
    
    /*
     * 2008:05	Count the total number of h-bonds in the file so that the
     *		hbdilution plot can have the correct index of the h-bond
     *		broken to generate the corresponding cluster analysis
     */
    if( !strncmp(linebuf,"REMARK",6) && \
	( !strncmp(linebuf+55,"HB",2) || !strncmp(linebuf+55,"SB",2) ))
    {
      (*total_hbond_count)++;
    }

    is_atom   = !strncmp( linebuf, "ATOM", 4 );
    sscanf( linebuf+13, "%3s", atom );

    if( is_atom && 
	( !strcmp( atom, "N" ) ||
	  !strcmp( atom, "CA") ||
	  !strcmp( atom, "C" ) ) ) {

      temp = ( atom_list *) malloc( sizeof( atom_list ));
      temp->atom_number = atoi( (linebuf+5)  );
      sscanf( linebuf+13, "%3s", temp->atom_type );
      sscanf( linebuf+22, "%4d", &temp->residue );
      temp->chain_ID = *column(22);
      temp->atom_or_hetatm = 1;

      /**********************************************************************/
      /* The following code was added to handle insertions that may be found*/
      /* in the PDB file. Assumes PDB file format standards.                */
      /**********************************************************************/
      if( *(linebuf+26) != ' ' ){
	
	if( *total_insertions == 0 ){
	  if( !strcmp( atom, "N" ) ){
	    insertion_res_num[*total_insertions] = temp->residue;
	    insertion_chain_ID[*total_insertions] = temp->chain_ID;
	    temp->insertion_space = 1; 
	  }
	  if( !strcmp( atom, "CA" ) )
	    temp->insertion_space = 1; 
	  if( !strcmp( atom, "C" ) ){
	    temp->insertion_space = 1; 
	    (number_of_insertions[*total_insertions])++;
	    (*total_insertions)++;
	  }
	}

      	else{
	  for( a = 0; a < *total_insertions; a++ ){
	    if( temp->residue == insertion_res_num[a] &&
		temp->chain_ID == insertion_chain_ID[a] ){
	      if( !strcmp( atom, "N" ) )
		temp->insertion_space = number_of_insertions[a] + 1; 
	      if( !strcmp( atom, "CA" ) )
		temp->insertion_space = number_of_insertions[a] + 1; 
	      if( !strcmp( atom, "C" ) ){
		temp->insertion_space = number_of_insertions[a] + 1; 
		(number_of_insertions[a])++;
	      }
	      known_insertion = -1; 
	    }
	  }
	  
	  if( !known_insertion ){
	    if( !strcmp( atom, "N" ) ){
	      insertion_res_num[*total_insertions] = temp->residue;
	      insertion_chain_ID[*total_insertions] = temp->chain_ID;
	      temp->insertion_space = 1; 
	    }
	    if( !strcmp( atom, "CA" ) )
	      temp->insertion_space = 1; 
	    if( !strcmp( atom, "C" ) ){
	      temp->insertion_space = 1; 
	      (number_of_insertions[*total_insertions])++;;
	      (*total_insertions)++;
	    }
	  }
	  known_insertion = 0; 
	}
      }
      /**********************************************************************/

      /*printf("residue: %d   total %d: number_of_insertions at residue %d: %3d\n", temp->residue, 
	     *total_insertions, insertion_res_num[(*total_insertions)-1],
	     number_of_insertions[(*total_insertions)-1] );*/

      if( *chem_file_head == NULL ){
	*chem_file_head = current = temp;
	y_translate[*number_of_chains] = current->residue;
	chain_label = temp->chain_ID;
	chain_IDs[*number_of_chains] = temp->chain_ID;
      }
      else{
	current->next_atom = temp;
	current = current->next_atom;
	current->next_atom = NULL;
      }
           
      if( !strcmp( atom, "CA") ) {
	if( current->chain_ID == chain_label ) {
	   (chain_size[*number_of_chains])++;

	   /* SN 2007:12:19 */
           if( current->residue < min_max_rsd[*number_of_chains][0] )
              min_max_rsd[*number_of_chains][0] = current->residue;
           if( current->residue > min_max_rsd[*number_of_chains][1] )
              min_max_rsd[*number_of_chains][1] = current->residue;
	   
           if( prev_rsd == 0)
	   {
	     prev_rsd = current->residue;
	     i = 0;
	   }
	   else if(current->residue == prev_rsd + 1)
	   {
	     prev_rsd++;
	   }
	   else if(current->residue > prev_rsd + 1)
	   {
             for(prev_rsd += 1;prev_rsd < current->residue; prev_rsd++,i++)
	     {
    		if(i == MISSING_RSD_LIST_SIZE)
    		{
  printf("\n**************************************************************");
  printf("\n* WARNING: Missing residue list is full!                     *");
  printf("\n*	       To prevent memory leak the program will exit now! *");
  printf("\n*	       Please change MISSING_RSD_LIST_SIZE and re-run.   *");
  printf("\n**************************************************************\n");
      		  exit(0);	
    		}
	       missing_rsd_list[i] = prev_rsd;
	     }
	   }
	   /* End of changes SN 2007:12:19 */

	}
	else{
	  chain_label = current->chain_ID;
	  (*number_of_chains)++;

	      /* SN 2007:12:19 */
              min_max_rsd[*number_of_chains][0] = current->residue;
              min_max_rsd[*number_of_chains][1] = current->residue;
	      prev_rsd = current->residue;	
	      /* End of changes SN 2007:12:19 */

	  y_translate[*number_of_chains] = current->residue;
	  (chain_size[*number_of_chains])++;
	  chain_IDs[*number_of_chains] = current->chain_ID;
	}

	if( current->residue < y_translate[*number_of_chains] )
	  y_translate[*number_of_chains] = current->residue;
	
      }
    }
  }
  
  /***************************************************************************/ 
  /*	SN 2007:12:19	 Debug missing residues => adjust chain_size!        */

#ifdef DEBUG_SN
/*  printf("\n------------- DEBUG OUTPUT - read_chem_file() - START ----------");*/
  printf("\n------- Inserted and Missing Residue Information ----------");
#endif

  for(a = 0; a < 100; a++)
  {
    if(number_of_insertions[a] > 0)
    {
#ifdef DEBUG_SN
      printf("\n\nResidue #:	Chain-ID:	# of insertions:");
      printf("\n   %d		  %c		  %d",insertion_res_num[a],\
				    insertion_chain_ID[a],\
				    number_of_insertions[a]);
#endif

      for(i = 0; i <= *number_of_chains; i++)
      {
        if(insertion_chain_ID[a] == chain_IDs[i])
         insertion_count[i] += number_of_insertions[a];
      }
    }
  }

#ifdef DEBUG_SN
  printf("\n\n");
  printf("                         Curr.                   Actual Missing\n");
  printf("Chain #: ID: Insertions: Size: Min_rsd: Max_rsd: Size:  Rsd:");

  i = j = 0; /* Used as an index into the missing rsd list hereafter*/
#endif

  for(a = 0; a <= *number_of_chains; a++)
  {
    missing_rsd_count[a] = min_max_rsd[a][1] - min_max_rsd[a][0] + 1 + \
			insertion_count[a] - chain_size[a];
#ifdef DEBUG_SN
    printf("\n  %d\t %c\t%d\t  %d\t %d\t %d\t %d\t%d",\
		a,chain_IDs[a],insertion_count[a],\
		chain_size[a],min_max_rsd[a][0],min_max_rsd[a][1],\
		min_max_rsd[a][1]-min_max_rsd[a][0]+1+insertion_count[a],\
		missing_rsd_count[a]);
#endif

    chain_size[a] = chain_size[a] + missing_rsd_count[a]; /* Actual Size */

#ifdef DEBUG_SN
    if(missing_rsd_list[i] > 0)
    {
      printf("\nMissing residue list:\n");
      
      for(i = j;i < j+missing_rsd_count[a] && i < MISSING_RSD_LIST_SIZE; i++)
      {
	if(missing_rsd_list[i] > 0)
	  printf("%d ",missing_rsd_list[i]);
      } 
      j += missing_rsd_count[a];
    }
    else
      printf("\nNo missing residues;\n");
    printf("\n");
#endif

  }

#ifdef DEBUG_SN
/*  printf("\n------------- DEBUG OUTPUT - read_chem_file() - END ----------\n");*/
  printf("\n------------- END OF INFORMATION ----------\n");
#endif

}
/**********************************************************************/
/**********************************************************************/

/**********************************************************************/
/* This rountine will determine the number of clusters in the current */
/* decomposition. Since only main-chain atoms are counted, a lower    */
/* bound on cluster size must be set. The variable is minimum_cluster */
/* _size. I have been using 3, corresponding to a single residue. The */
/* date for this calculation comes from the decomp_list file.         */
/**********************************************************************/
  int compute_number_of_clusters(int 	   atom_label[MAXATOMS], 
				 atom_list *chem_file,
				 int 	   cluster_label[MAX_CLUSTER_COUNT], 
				 int 	   number_of_atoms ){
  
  int 
    a = 0,
    tally = 0, 
    current_number = 0, 
    total_clusters = 0, 
    current_label = 1, 
    renew_list = 0, 
    check = 0,
    minimum_cluster_size = 3,
    upper_bound_on_number_of_clusters = 0;
  
  atom_list 
    *start_of_list = NULL;


  /**********************************************************************/
  /* The upper bound on the number of clusters is set to limit the      */
  /* number of times the decomp file needs to be read. Mainchain rigid  */
  /* clusters will not necessarily be labelled by consecutive integers. */
  /* side-chains and hetero groups can form rigid clusters whose label  */
  /* will (possibly) be smaller than a main-chain rigid cluster. It is  */
  /* therefore impossible to simply stop checking the decomp file the   */
  /* first time it fails to find a main-chain rigid cluster. For example*/
  /* rigid cluster 23 may consist of 10 main-chain atoms, so it will be */
  /* added to the cluster_label array. However, cluster 24 might be all */
  /* heteratoms, so we don't add it. cluster 25 could then be 9 main-   */
  /* chain atoms, which we want to pick up. Therefore, i set the upper  */
  /* bound to be 300, as a precaution to make sure we find all the main */
  /* chain clusters of size "minimum_cluster_size" or greater.          */
  /* For the dimeric protein 2cts, with 876 total residues, the largest */
  /* cluster label found was 167, so 300 should be approriate, for now. */
  /*     The size of nc_count and oc_count is related to the total      */
  /* number of clusters found ("total_clusters"), and does not depend on*/
  /* the actual value of the cluster label.                             */
  /**********************************************************************/
  upper_bound_on_number_of_clusters = 300;

  start_of_list = chem_file;

  for( a = 1; a < upper_bound_on_number_of_clusters; a++ ){

    while( chem_file ){

      if( atom_label[chem_file->atom_number] == a )
	tally++;
      
      if( tally == minimum_cluster_size ) {
	/*printf("cluster %d, label %d %d\n", total_clusters, a, chem_file->atom_number );*/
	cluster_label[total_clusters] = a;
	total_clusters++;
	tally = 0;
	chem_file = NULL;
      }
      
      if( chem_file != NULL )
	chem_file = chem_file->next_atom;
      
    }
    chem_file = start_of_list;
    tally = 0;
  }
  
  return( total_clusters );
}    

/**********************************************************************/
/* Using the chem_file_head list, this routine stores the data from   */
/* initial decompostion in the array of structures, nc_count. The ini */
/* decomposition refers to the protein with all the H-bonds present.  */
/**********************************************************************/
void first_decomp( 	int long_output, int total_clusters, 
			int bulk_atom_label[MAXATOMS],
		  	 atom_list *chem_file, 
			 clusters *new_clusters[MAX_CLUSTER_COUNT], 
			int colors[MAX_CLUSTER_COUNT],
		   	int cluster_labels[MAX_CLUSTER_COUNT] ) {

  int current_cluster=0, cluster_number=0;

  atom_list    *start_of_list = NULL;

  clusters     *nc_temp    = NULL,
               *nc_current = NULL;
  
  start_of_list = chem_file;

  if( !long_output ){

    for( current_cluster = 0; current_cluster < total_clusters; current_cluster++ ) {

      colors[current_cluster] = current_cluster+1;
      chem_file = start_of_list;

      while( chem_file ){
	
	/************************************************************/
	/* the array bulk_atom_label contains the cluster label     */
	/* assigned by first. Clusters are number consecutively, the*/
	/* largest cluster is 1, the next largest 2, ... . Which    */
	/* cluster an atom in the protein belongs to is indexed by  */
	/* bulk_atom_label. We only look at clusters of a certain   */
	/* size. The variable total_clusters indicates how many     */
	/* clusters we have in the current decomp.                  */
	/************************************************************/	  
	cluster_number = bulk_atom_label[chem_file->atom_number];

	if(cluster_number == cluster_labels[current_cluster] ){

	  nc_temp = (clusters *) malloc( sizeof( clusters ));
	 
	  if( new_clusters[current_cluster] == NULL ){
	    new_clusters[current_cluster] = nc_temp;
	    nc_current = nc_temp;
	  }
	  
	  nc_temp->atom_number     = chem_file->atom_number;
	  nc_temp->residue_number  = chem_file->residue;
	  nc_temp->insertion_space = chem_file->insertion_space;
	  sscanf( chem_file->atom_type, "%3s", nc_temp->atom_type );
	  nc_temp->chain_ID = chem_file->chain_ID;
	  nc_temp->cluster_color   = current_cluster+1;
	  
	  nc_current->next_element = nc_temp;
	  nc_current = nc_current->next_element;
	  nc_current->next_element = NULL;
	  
	}
	
	
	chem_file = chem_file->next_atom;
	
      }
    }
  }
  
}

/**********************************************************************/
/* store the atom info from the chem_file corresponding to the new    */
/* clusters in an array of linked lists called nc_count.              */
/**********************************************************************/
void set_new_decomp_info( 	int total_clusters, 
				int bulk_atom_label[MAXATOMS], 
			  	int atom_num,  
				clusters *new_decomp[MAX_CLUSTER_COUNT], 
			  	atom_list *chem_file, 
				int cluster_labels[MAX_CLUSTER_COUNT] ) {
  
  int cluster_number = 0,
      current_cluster = 0,
      found_cluster = 0,
      this_cluster = 0;
  
  atom_list *start_of_list = NULL;
  
  clusters *nc_current = NULL,
           *nc_temp    = NULL;

  start_of_list = chem_file;

  for( current_cluster = 0; current_cluster < total_clusters; current_cluster++ ) {
    
    chem_file = start_of_list;
    
    while( chem_file ){

      cluster_number = bulk_atom_label[chem_file->atom_number];

      if( cluster_number == cluster_labels[current_cluster] ){
	nc_temp = (clusters *) malloc( sizeof( clusters ));

	if( new_decomp[current_cluster] == NULL ){
	  new_decomp[current_cluster] = nc_temp;
	  nc_current = nc_temp;
	}

	nc_temp->atom_number      = chem_file->atom_number;
	nc_temp->residue_number   = chem_file->residue;
	nc_temp->insertion_space  = chem_file->insertion_space;
	sscanf( chem_file->atom_type, "%3s", nc_temp->atom_type );
	nc_temp->chain_ID = chem_file->chain_ID;
	nc_temp->next_element = NULL;
	
	nc_current->next_element = nc_temp;
	nc_current = nc_current->next_element;
	nc_current->next_element = NULL;
      }
      
      chem_file = chem_file->next_atom;
      
    }
  }
  
}
/**************************************************/

/**************************************************/
void swap_cluster_info( clusters *array[MAX_CLUSTER_COUNT], int element_A, int element_B, 
			int size[MAX_CLUSTER_COUNT] ){
  
  int temp_size[MAX_CLUSTER_COUNT];
  
  clusters *temp_struct[2];
  
  temp_size[0] = size[element_A];
  temp_struct[0] = array[element_A];

  size[element_A] = size[element_B];
  array[element_A] = array[element_B];

  size[element_B] = temp_size[0];
  array[element_B] = temp_struct[0];
  
}
/**************************************************/

/**************************************************/
void bsort( clusters *array[MAX_CLUSTER_COUNT], int size[MAX_CLUSTER_COUNT], int start, int end ){
  
  int a=0, last=0;
  
  if( start >= end )
    return;  
  swap_cluster_info( array, start, (start+end)/2, size );
  last = start;
  
  for( a = start+1; a <= end; a++ ){
    if( size[a] > size[start] ){
      swap_cluster_info( array, ++start, a, size );
    }
  }
  swap_cluster_info( array, start, last, size );
  bsort( array, size, start, last-1 );
  bsort( array, size, last+1, end );
}
/*****************************************************/

/*****************************************************/
void swap_colors( int array[MAX_CLUSTER_COUNT], int A, int B ){
  int temp;
  temp = array[A];
  array[A] = array[B];
  array[B] = temp;
}
/*****************************************************/

/*****************************************************/
void csort( int colors[MAX_CLUSTER_COUNT], int start, int end ) {

  int 
    a = 0, 
    last = 0;

  if( start >= end )
    return;
  swap_colors( colors, start, (start+end)/2 );
  last = start;
  for( a = start+1; a <= end; a++ ){
    if( colors[a] > colors[start] )
      swap_colors( colors, ++start, a );
  }
  swap_colors( colors, start, last );
  csort( colors, start, last-1 );
  csort( colors, last+1, end   );
}
/****************************************************/

/****************************************************/
/* This array is designed to compute the cluster    */
/* coloring, so as to keep it consistent through-   */
/* out the H-bond dilution plots. It an old cluster */
/* has more than one subset in the new decomposition*/
/* the largest of the new subsets is given the color*/
/* of the previous color, and a new color is assig- */
/* ned to the remaining subcluster(s)               */
/****************************************************/
void compute_new_cluster_list(clusters *nc_list[MAX_CLUSTER_COUNT], int num_new,
			      clusters *oc_list[MAX_CLUSTER_COUNT], int num_old,
			      int old_colors[MAX_CLUSTER_COUNT], int colors[MAX_CLUSTER_COUNT]){
  
  int 
    a = 0, 
    b = 0, 
    c = 0, 
    d = 0, 
    e = 0, 
    total_count = 0,
    color_count = 1,
    is_a_subset = 0, 
    subset_index = 0, 
    subset_element = 0, 
    total_atoms[MAX_CLUSTER_COUNT], 
    current_color_index = 0,
    atom_count = 0, 
    sort_color[MAX_CLUSTER_COUNT],
    t = 0, 
    new_colors[MAX_CLUSTER_COUNT];
  
  clusters 
    *new_current,
    *old_current,
    *subset_list[MAX_CLUSTER_COUNT],
    *new_list[MAX_CLUSTER_COUNT];

  /**********************************************************************/
  /* Initialize variables.                                              */
  /**********************************************************************/
  for( d = 0; d < num_old; d++ )
    sort_color[d] = old_colors[d];
  
  for( d = 0; d < MAX_CLUSTER_COUNT; d++)
    new_colors[d] = 99999;
  /**********************************************************************/

  csort( sort_color, 0, num_old-1 );
  current_color_index = sort_color[0];

  /***********************************************/
  /* First, find all new clusters which are sub- */
  /* sets of the first old cluster. If there are */
  /* more than one, sort by size and recolor. The*/
  /* new color will be 1 + the number of old     */
  /* clusters.                                   */
  /***********************************************/
  for( a = 0; a < num_old; a++ ){
    
    subset_index   = 0;
    subset_element = 0;

    /* Find all the new clusters that are subsets of old cluster "a" */
    for( b = 0; b < num_new; b++ ){
      
      is_a_subset = 0;
      new_current = nc_list[b];
      old_current = oc_list[a];
      
      while(old_current){
	
	if(new_current->residue_number == old_current->residue_number &&
	   !strcmp( new_current->atom_type, old_current->atom_type )  &&
	   ( new_current->chain_ID == old_current->chain_ID  ) &&
	   new_current->atom_number == old_current->atom_number ) {
	  is_a_subset = -1;
	  old_current = NULL;
	}
	if( old_current )
	  old_current = old_current->next_element;
      }

      if( is_a_subset ){
	subset_list[subset_index] = nc_list[b];
	new_current = nc_list[b];
	while(new_current){
	  atom_count++;
	  new_current = new_current->next_element;
	}
	total_atoms[subset_index] = atom_count;
	subset_index++;	
	atom_count = 0;
      }
    }

    /* if there is at least one subset of old cluster "a" */
    if( subset_index ){     
      
      bsort( subset_list, total_atoms, 0, subset_index-1 );

      for( t = 0; t < subset_index; t++ ){
	if( t == 0){
	  new_list[total_count] = subset_list[t];
	  new_colors[total_count] = old_colors[a];
	  total_count++;
	}
	else{    /* need to sort color_nums, and take 1+ the largest */
	  new_list[total_count] = subset_list[t];
	  new_colors[total_count] = current_color_index + color_count;
	  total_count++;
	  color_count++;
	}
      }
    }
    
  }

  /*printf("num_new %d, total_count %d\n", num_new, total_count );*/
  for( c = 0; c < num_new; c++ ) {
    nc_list[c] = new_list[c];
    colors[c]  = new_colors[c];
  }

}
/**********************************************************************/

/**********************************************************************/
/* Tallies the number of main-chain atoms in new clusters and old     */
/* clusters. If the main-chain count did not change, the output of the*/
/* current decomp is not produced.                                    */
/**********************************************************************/
int compare_old_and_new_decomps( clusters *new[MAX_CLUSTER_COUNT], clusters *old[MAX_CLUSTER_COUNT], 
				 int cluster_counter, int cluster_sizes[MAX_CLUSTER_COUNT], 
				 int old_cluster_count ) {

  int 
    a = 0, 
    size_new = 0, 
    size_old = 0, 
    this_cluster_size = 0;
  
  clusters  
    *nc_current = NULL,
    *oc_current = NULL;
 
 /*
	printf("\nCHK PT 2.73359 \n"); 
	 fflush(stdout);
 */

  for( a = 0; a < cluster_counter; a++ ){
    
    nc_current = new[a];
   
/* printf("\n DEBUG: cluster number = %d",a); */
 
    while( nc_current ){
      
      size_new++;
    /*        printf("size: %4d (%5d %4d %c)\n",size_new, nc_current->atom_number, nc_current->residue_number, nc_current->chain_ID );*/
      fflush(stdout);
      this_cluster_size++;
      nc_current = nc_current->next_element;
    }
    cluster_sizes[a] = this_cluster_size - 1;
    this_cluster_size = 0;
  }
  for( a = 0; a < old_cluster_count; a++ ){
    
    oc_current = old[a];
    
    while( oc_current ){
      size_old++;
      oc_current = oc_current->next_element;
    }
  }

/*
 printf("\nCHK PT 2.73361 \n"); 
 fflush(stdout); 
 printf("size new %d (%3d clusters) size old %d (%3d clusters)\n", size_new, cluster_counter, size_old, old_cluster_count );*/
  
  if( size_new == size_old ) {
    return(1);
  }
  return(0);
}

/**********************************************************************/
/* If there are more than 220 residues in the protein, the output is  */
/* displayed on multiple pages, in landscape format. In order to keep */
/* a clean postscript file, the output for the extra pages is sent to */
/* temporary files, one for each extra page needed. When the last     */
/* decomp of the page has finished, the tempfiles will be catted to   */
/* the output file cluster.ps, and deleted. The array file_list holds */
/* the addresses of the tempfiles.                                    */
/**********************************************************************/
void set_temp_files( FILE *file_list[FILE_LIST_SIZE], int *number_of_pages, int *ps_page_number) {

  int a = 0;

  char 
    file_name[15];

  FILE *fp;
  
  for( a = 1; a <= *number_of_pages; a++ ) {
/*    printf("setting up temp files...\n"); */
    sprintf( file_name, "%d", a );
    fp = fopen( file_name,"w+");
    (*ps_page_number)++;
    fprintf(fp, "%%%%Page:     %d   %d\n\n", *ps_page_number, *ps_page_number );
    fprintf(fp, "0 109 620 109 sequence_line\n");
    file_list[a] = fp;
  }
}

/**********************************************************************/
/* Clean up unnecessary files and stuff.                              */
/**********************************************************************/    
void clean_up( FILE *ps_file, int number_of_pages ){

  int  a=0;
  char file_name[4], linebuf[300],cmd[20];
  FILE *current_file;

  fprintf( ps_file, "showpage\n\n" );
  fflush( ps_file );

  for( a = 1; a <= number_of_pages; a++ ) {
    sprintf( file_name, "%d", a );
    current_file = fopen( file_name, "r" );
    while( fgets( linebuf, sizeof(linebuf), current_file ) != NULL ){
      fprintf(ps_file, "%s", linebuf );
    }
    fprintf(ps_file, "showpage\n");
  }
  
  fprintf(ps_file, "%%%%EOF\n\n");
  
  if( number_of_pages )    
   for( a = 1; a <= number_of_pages; a++ )
    {
      sprintf(cmd,"\\rm -f %d\n",a);
      system(cmd);
    }
}


  /*****************************************************/
  /* added to sort the rigid_regions array             */
  /*****************************************************/
  
  void quick_sort( int *array, int start, int end ) {
  
    int
      a = 0,
      last = 0;
  
    if( start >= end )
      return;
    swap_colors( array, start, (start+end)/2 );
    last = start;
    for( a = start+1; a <= end; a++ ){
      if( array[a] < array[start] )
        swap_colors( array, ++start, a );
    }
    swap_colors( array, start, last );
    quick_sort( array, start, last-1 );
    quick_sort( array, last+1, end   );
  }
 
  /**********************************************************************/
  /* This routine is invoked only when user chooses text-only output    */
  /* detailing the flexible residue ranges only with each h-bond break. */
  /* Outout is written to a text file instead of generating a postscript*/                
  /**********************************************************************/
  int open_output_file( char file_name[128], FILE** file ){
      
   
    char temp[128];
    if( !fopen(file_name,"r") ){
      *file = fopen( file_name, "w" );
      return(1);
    }
    else{
      printf("\n\t There is already a copy of %s in this directory.\n", file_name ); 
      strcpy(temp, file_name);
      strcat(temp,".tmp");
      printf("\t The output will be saved in %s\n", temp);
      
      if( !fopen(temp, "r") ){
        *file   = fopen(temp, "w" );
         return(1);
      }
      else{
        printf("\n\t There is also a tmp.ps in this directory. The program will\n");
        printf("\t terminate. After you rename or remove %s and %s, you\n", file_name,temp);
        printf("\t may run the program from the command line as follows:\n");
        printf("\t hbdilute decomp_list s *_proflexdataset\n\n");
        return(0);						       
      }
     
    }   
      return(1); 
      
  }

