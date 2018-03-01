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

/* Structure declarations for hbdilute.c */

#ifndef _TYPES_H_
#define _TYPES_H_
					/* SN	2006:08	*/
#define	MAX_CLUSTER_COUNT	1000	
					/* Set based on runs for "urease" */
#define FILE_LIST_SIZE		15	
#define MAXATOMS		100000	
					/* Same as first and util programs */
#define MISSING_RSD_LIST_SIZE	100	
					/* SN 2008:01	*/
#define MAX_RSD			10000	
					/* Used in hbdilute.c */
#define MAX_CHAIN_COUNT		100

typedef struct chem_file{
  int  atom_or_hetatm;
  int  atom_number;
  char atom_type[4];
  int  residue;
  char chain_ID;
  int  insertion_space;
  struct chem_file *next_atom;
  struct chem_file *last_atom;
} atom_list;

typedef struct rcl{
  int  start_num;
  char start_type[4];
  int  end_num;
  char end_type[4];
  int  set_number;
  int  subset_number;
} rclist[30];

typedef struct n_c {
  int   atom_number;
  int   residue_number;
  int   insertion_space;
  char  atom_type[4];
  char  chain_ID;
  int   cluster_color;
  struct n_c *next_element;
} clusters;

#endif
