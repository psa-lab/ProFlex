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

#ifndef _POSTSCRIPT_
#define _POSTSCRIPT_

#include "types.h"

void print_header( FILE *, float );
void print_data_headings( FILE *, int, int, int, int[]  );
void print_current_Hbond_info_portrait( FILE *, int, int, int, int *, int, float, 
					int, int [], int, int, int, float, int *, 
					int, int, int, int[], char [], char [], int,
					char, char, char[], char[], int, float, int );

void print_current_Hbond_info_landscape( FILE *, int, int, int, int *, int, float, 
					 int, int [], int, int, int, float, int *, 
					 int, int, int, int[], char [], char [], int,
					 char, char, char[], int[], char[], int, float, int,
					 /* SN 2008:01 */ int, int, char, char, int );

void print_current_Hbond_info_landscape_multipage( FILE *, int, int, int, int *, int, float, 
						   int, int [], int, int, int, float, int *, 
						   int, int, int, int[], char [], char [],
						   FILE *[], int, int, char, char, char[], 
						   char[], int[], float,
					 /* SN 2008:01 */ int, int, char, char, int );
void print_decomp( clusters *[], int, int, FILE *, int [], int[], int, int, 
		   int[], char[], int, int[], int[], char[], 
		  /* SN 2008:01*/ int [][2], int [], int [] );
void print_multipage_decomp( clusters *[], int, int, FILE *, int [], int[], 
			     int, FILE *[], int, int, int[], char[], int, 
			     int[], int[], char[],
		  /* SN 2008:01*/ int [][2], int [], int [] );
void print_portrait_footer( FILE *, char[], int );
void print_landscape_footer( FILE *, char[], int, float ); /* SN 2008:01 */

#endif
