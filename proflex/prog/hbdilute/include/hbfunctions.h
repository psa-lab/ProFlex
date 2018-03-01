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

#ifndef _HDFUNCTIONS_
#define _HDFUNCTIONS_

#include "types.h"

void read_chem_file( atom_list **, FILE *, int[], int[], char[], int *, 
		     int *, int[], int[], char[],
		     /* 2008:01 SN */ int [][2], int [], int [], int* );
int  compute_number_of_clusters( int [], atom_list *, int [], int );
void first_decomp( int, int, int[], atom_list *, clusters *[], int [], int [] );
void set_new_decomp_info( int, int[], int, clusters *[], atom_list *, int [] );
void compute_new_cluster_list( clusters *[], int, clusters *[], int, int [], int [] );
int  compare_old_and_new_decomps( clusters *[], clusters *[], int, int[], int );
void set_temp_files( FILE *[], int *, int * );
void clean_up( FILE *, int);
int open_output_file( char [], FILE**  );

#endif


