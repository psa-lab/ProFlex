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
/* Various bits and pieces of the output postscript file    */
/* need to be included. The lines were hard coded in to eli-*/
/* minate the need for the buch of files in a specific      */
/* directory that need to be catted to the final file in a  */
/* certain order.                                           */
/************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "../include/types.h"

#define TRUE  1
#define FALSE 0

#define OFFSET 0.5	/* Used in conjunction with start_pt for printing
			   rsd numbers and the corresponding markers */

/*
#define DEBUG_SN
*/
/*
#define DEBUG_MULTIPAGE
*/

/* 2008:01 SN	Changes
 *
 * 	Change the font to Arial from Times-Roman
 */

/************************************************************/
/*            The postscript header info                    */
/* Includes standard stuff to make the file "legal" level 2 */
/* postscript, as well as definitions for the colors used   */
/* and the geometric primitives.                            */
/************************************************************/
void print_header(FILE *ps_file, float scale ){

  time_t current_time;

  current_time = time(NULL);

  fprintf(ps_file, "%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(ps_file, "%%%%Title: HBD.ps\n");
  fprintf(ps_file, "%%%%Creator: Brandon Hespenheide\n");
  fprintf(ps_file, "%%%%CreationDate: 1999\n");
  fprintf(ps_file, "%%%%BoundingBox:0 0 620 800\n");
  fprintf(ps_file, "%%%%DocumentFonts: Arial\n");
  fprintf(ps_file, "%%%%EndComments\n");
  fprintf(ps_file, "%%%%BeginProlog\n");
  fprintf(ps_file, "/Col { sethsbcolor } bind def\n");
  fprintf(ps_file, "%% These are the colors for the flexibility scale and the\n");
  fprintf(ps_file, "%% lines that display the hydrogen bonds.\n\n");

  fprintf(ps_file, "/Blue {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/Red  {1.0000  0.0000  0.0000 setrgbcolor } def\n");

  fprintf(ps_file,"/Col00 {1.0000 0.0000 0.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col01 {0.0000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col02 {0.3000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col03 {0.6000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col04 {0.9000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col05 {0.2000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col06 {0.5000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col07 {0.8000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col08 {0.1000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col09 {0.4000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col10 {0.7000 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col11 {0.0000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col12 {0.3000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col13 {0.6000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col14 {0.9000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col15 {0.2000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col16 {0.5000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col17 {0.8000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col18 {0.1000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col19 {0.4000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col20 {0.7000 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col21 {0.0250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col22 {0.3250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col23 {0.6250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col24 {0.9250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col25 {0.2250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col26 {0.5250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col27 {0.8250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col28 {0.1250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col29 {0.4250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col30 {0.7250 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col31 {0.0250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col32 {0.3250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col33 {0.6250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col34 {0.9250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col35 {0.2250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col36 {0.5250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col37 {0.8250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col38 {0.1250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col39 {0.4250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col40 {0.7250 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col41 {0.0500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col42 {0.3500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col43 {0.6500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col44 {0.9500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col45 {0.2500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col46 {0.5500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col47 {0.8500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col48 {0.1500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col49 {0.4500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col50 {0.7500 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col51 {0.0500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col52 {0.3500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col53 {0.6500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col54 {0.9500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col55 {0.2500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col56 {0.5500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col57 {0.8500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col58 {0.1500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col59 {0.4500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col60 {0.7500 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col61 {0.0750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col62 {0.3750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col63 {0.6750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col64 {0.9750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col65 {0.2750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col66 {0.5750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col67 {0.8750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col68 {0.1750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col69 {0.4750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col70 {0.7750 1.0000 1.0000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col71 {0.0750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col72 {0.3750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col73 {0.6750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col74 {0.9750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col75 {0.2750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col76 {0.5750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col77 {0.8750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col78 {0.1750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col79 {0.4750 0.8000 0.8000 sethsbcolor } def\n");
  fprintf(ps_file,"/Col80 {0.7750 0.8000 0.8000 sethsbcolor } def\n");

  fprintf(ps_file, "/FlexCol01 {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol02 {0.1000  0.1000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol03 {0.2000  0.2000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol04 {0.3000  0.3000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol05 {0.4000  0.4000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol06 {0.4500  0.4500  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol07 {0.5000  0.5000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol08 {0.6000  0.6000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol09 {0.7000  0.7000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol10 {0.8000  0.8000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol11 {1.0000  0.8000  0.8000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol12 {1.0000  0.7000  0.7000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol13 {1.0000  0.6000  0.6000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol14 {1.0000  0.5000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol15 {1.0000  0.4500  0.4500 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol16 {1.0000  0.4000  0.4000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol17 {1.0000  0.3000  0.3000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol18 {1.0000  0.2000  0.2000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol19 {1.0000  0.1000  0.1000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol20 {1.0000  0.0000  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol21 {0.5000  0.5000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol30 {1.0000  0.0000  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol31 {0.7500  0.7500  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol32 {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol33 {1.0000  0.5000  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol34 {0.0000  1.0000  0.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol35 {0.9000  0.9000  0.9000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol36 {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol37 {0.5000  1.0000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol38 {0.5000  0.5000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol39 {0.5000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol40 {1.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol41 {1.0000  0.5000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol42 {0.2500  1.0000  0.2500 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol43 {0.0000  0.7500  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol44 {0.8000  1.0000  0.2500 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol45 {0.8000  1.0000  0.2500 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol46 {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol47 {0.5000  1.0000  0.5000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol48 {0.5000  0.5000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol49 {0.5000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_file, "/FlexCol50 {1.0000  0.0000  1.0000 setrgbcolor } def\n");


  fprintf(ps_file, "/Poly4 { moveto lineto lineto lineto fill } bind def\n");
  fprintf(ps_file, "/Pl4 { 8 copy Poly4 moveto moveto moveto moveto closepath stroke } bind def\n");
  fprintf(ps_file, "/Poly3 { moveto lineto lineto fill } bind def\n");
  fprintf(ps_file, "/Pl3 { 6 copy Poly3 moveto moveto moveto closepath stroke } bind def\n");
  fprintf(ps_file, "/Print { /Arial findfont exch scalefont setfont show } bind def  \n");
  fprintf(ps_file, "/sequence_line { moveto lineto stroke } bind def\n\n");
  fprintf(ps_file, "/CenterRot90 {\n");
  fprintf(ps_file, "  dup /Arial findfont exch scalefont setfont\n");
  fprintf(ps_file, "  exch stringwidth pop -2 div exch 3 div exch rmov\n");
  fprintf(ps_file, " } bind def\n");
  fprintf(ps_file, "/UncenterRot90 {\n");
  fprintf(ps_file, "  dup /Arial findfont exch scalefont setfont\n");
  fprintf(ps_file, "  exch stringwidth } bind def\n");
  fprintf(ps_file, "/MiniUncenterRot90 {\n");
  fprintf(ps_file, "  dup /CourierNew findfont exch scalefont setfont\n");
  fprintf(ps_file, "  exch stringwidth } bind def\n");
  fprintf(ps_file, "/Rot90 { gsave currentpoint translate 90 rotate } bind def\n");
  fprintf(ps_file, "%%%%EndProlog\n\n");
  fprintf(ps_file, "%%%%Page:    1   1\n");

  fprintf( ps_file, "470 770 moveto\n");
  fprintf( ps_file, "(%s) 10.0 Print\n", ctime(&current_time) );

#ifdef DEBUG_MULTIPAGE
	scale = 1.0;
#endif

  fprintf(ps_file, "%3.2f %3.2f scale\n", scale, scale );
  fprintf(ps_file, "save\n\n");
}
/************************************************************/


/************************************************************/
/* print a line of information at the bottom of the last    */
/* page of postscript output. for PORTRAIT style output.    */
/************************************************************/
void print_portrait_footer( FILE *ps_file, char protein_name[128], int ps_line_number ){

  ps_line_number -= 10;
  fprintf(ps_file, "Blue\n");
  fprintf(ps_file, "40 %d moveto\n", ps_line_number );
  fprintf(ps_file, "(Blue:donor) 10 Print\n");
  fprintf(ps_file, "Red\n");
  fprintf(ps_file, "90 %d moveto\n", ps_line_number );
  fprintf(ps_file, "(Red:acceptor) 10 Print\n");
  fprintf(ps_file, "Col00\n");
  fprintf(ps_file, "180 %d moveto\n", ps_line_number);
  fprintf(ps_file, "(M:main-chain   S:side-chain   W:water   H:hetero-atom) 10 Print\n");
  fprintf(ps_file, "500 %d moveto\n", ps_line_number);
  fprintf(ps_file, "(%s) 12 Print\n", protein_name );
}
/************************************************************/


/************************************************************/
/* print a line of information at the bottom of the last    */
/* page of postscript output. for LANDSCAPE style output.   */
/************************************************************/
void print_landscape_footer( FILE *ps_file, char protein_name[128], 
			     int ps_line_number, float scale_factor ){

  ps_line_number += 10;
  fprintf(ps_file, "Blue\n");	/* Col07---Incorrect color, hence change */
  fprintf(ps_file, "%d 110 moveto\n", ps_line_number);
  fprintf(ps_file, "(Blue:donor) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(Blue:donor) 10 Print\n");
  fprintf(ps_file, "grestore\n");
  fprintf(ps_file, "Red\n");		/* Col01 */
  fprintf(ps_file, "%d 170 moveto\n", ps_line_number);
  fprintf(ps_file, "(Red:acceptor) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(Red:acceptor) 10 Print\n");
  fprintf(ps_file, "grestore\n");
  fprintf(ps_file, "Col00\n");
  fprintf(ps_file, "%d 260 moveto\n", ps_line_number);
  fprintf(ps_file, "(M:main-chain   S:side-chain   W:water   H:hetero-atom) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(M:main-chain   S:side-chain   W:water   H:hetero-atom) 10 Print\n");
  fprintf(ps_file, "grestore\n");
  
  /* 2008:01:24 SN	Portrait output substituted with landscape, */
  /*			which results in filename spill out of page */
  /*			Hence, change position if scale == 1.0      */

#ifdef DEBUG_MULTIPAGE
	scale_factor = 1.0;
#endif

  if(scale_factor == 1.0)
    fprintf(ps_file, "%d 600 moveto\n", ps_line_number);
  else
    fprintf(ps_file, "%d 675 moveto\n", ps_line_number);
  fprintf(ps_file, "(%s) 12 UncenterRot90 Rot90\n", protein_name );
  fprintf(ps_file, "(%s) 12 Print\n", protein_name );
  fprintf(ps_file, "grestore\n");

  /* 2008:01:31 SN	Add text to explain the right hand side legend */
  fprintf(ps_file, "Col00\n");
  fprintf(ps_file, "%d 260 moveto\n", ps_line_number+12);
  fprintf(ps_file, "(Letter immediately before residue number at the end of line is the chain ID) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(Letter immediately before residue number at the end of line is the chain ID) 10 Print\n");
  fprintf(ps_file, "grestore\n");

  fprintf(ps_file, "Col00\n");
  fprintf(ps_file, "%d 260 moveto\n", ps_line_number+24);
  fprintf(ps_file, "(Letter immediately after residue number (if it occurs) is the insertion code) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(Letter immediately after residue number (if it occurs) is the insertion code) 10 Print\n");
  fprintf(ps_file, "grestore\n");
  
  /* 2008:05 SN		Add text to explain what H-bond index means */
  fprintf(ps_file, "Col00\n");
  fprintf(ps_file, "%d 260 moveto\n", ps_line_number+42);
  fprintf(ps_file, "(*\tIndex of the last H-bond to be broken before generating the rigidity profile of this line) 10 UncenterRot90 Rot90\n");
  fprintf(ps_file, "(*\tIndex of the last H-bond to be broken before generating the rigidity profile of this line) 10 Print\n");
  fprintf(ps_file, "grestore\n");

}
/************************************************************/


/************************************************************/
/************************************************************/
void print_data_headings( FILE *ps_file, int long_output, 
			  int number_of_residues, int number_of_chains, 
			  int chain_size[MAX_CHAIN_COUNT] ){

  int a = 0, b = 0;

 /* Disable portrait printing */

  /* Print portrait *
  if( !number_of_chains && number_of_residues <= 125 ){

    fprintf(ps_file, "20 740 moveto\n");
    fprintf(ps_file, "(Remaining) 10.0 Print\n");
    fprintf(ps_file, "20 730 moveto\n");
    fprintf(ps_file, "(H-bonds) 10.0 Print\n");
    fprintf(ps_file, "60 730 moveto\n");
    fprintf(ps_file, "(E) 10.0 Print\n");
    fprintf(ps_file, "80 730 moveto\n");
    fprintf(ps_file, "(<r>) 10.0 Print\n");
  }
  else
  */

  /* Print landscape */
  {
    fprintf(ps_file, "33 20 moveto\n");
    fprintf(ps_file, "(H-bond) 9.0 UncenterRot90 Rot90\n");
    fprintf(ps_file, "(H-bond) 9.0 Print\n");
    fprintf(ps_file, "grestore\n\n");
    fprintf(ps_file, "43 20 moveto\n");
    fprintf(ps_file, "(Index *) 9.0 UncenterRot90 Rot90\n");
    fprintf(ps_file, "(Index *) 9.0 Print\n");
    fprintf(ps_file, "grestore\n\n");
    fprintf(ps_file, "43 63 moveto\n");
    fprintf(ps_file, "(E) 9.0 UncenterRot90 Rot90\n");
    fprintf(ps_file, "(E) 9.0 Print\n");
    fprintf(ps_file, "grestore\n\n");
    fprintf(ps_file, "43 79 moveto\n");
    fprintf(ps_file, "(<r>) 9.0 UncenterRot90 Rot90\n");
    fprintf(ps_file, "(<r>) 9.0 Print\n");
    fprintf(ps_file, "grestore\n\n");
  }
}
/************************************************************/


/************************************************************/
/************************************************************/
void print_current_Hbond_info_portrait( FILE *ps_file, int long_output, int short_output,
					int first_line, int *ps_line_number, int hb_number,
					float hb_energy, int cluster_counter, int colors[30],
					int HB_start, int HB_end, int NUMBER_OF_RESIDUES,
					float scale_factor, int *ps_page_number, int missing_one,
					int HBtwo_start, int HBtwo_end, int y_translate[10],
					char donor_atom_type[4], char accept_atom_type[4],
					int number_of_chains, char donor_chain_ID,
					char accept_chain_ID, char chain_IDs[10],
					char protein_name[128], int side_chain_test, float mean_coordination,
					int fit_it_all_on_one_page ){

  int
    a = 0,
    triangle_x = 0,
    triangle_y = 0,
    this_chain = 0;

  float
    y_start = 0.0,
    y_end = 0.0;

  /* 2007:06	hb_number is the array index of the h-bond list.
		The actual hb_number as output in proflexdataset 
		is hb_number+1!!
  */
  hb_number = hb_number+1;

  fprintf(ps_file, "Col00\n");

  if( !side_chain_test ){
/*    if( !first_line ){ */
      fprintf( ps_file, "21 %3d moveto\n", *ps_line_number-4 );
      fprintf( ps_file, "(%4d) 9.0 Print\n", hb_number );
      fprintf( ps_file, "44 %3d moveto\n", *ps_line_number-4 );
      fprintf( ps_file, "(%7.3f) 9.0 Print\n\n", hb_energy );
      fprintf( ps_file, "72 %3d moveto\n", *ps_line_number-4 );
      fprintf( ps_file, "(%7.3f) 9.0 Print\n\n", mean_coordination );
/*    }
    else{

printf("\n#$@#@   \t %f",hb_energy);
        if(hb_energy == 0.0)
	 {
		fprintf(ps_file, "36 %3d moveto\n", *ps_line_number-3 );
		fprintf(ps_file, "(All %3d Hbonds) 9.0 \n", hb_number);
    		fprintf(ps_file, "(All %3d Hbonds) 9.0 Print\n", hb_number);
    		fprintf(ps_file, "grestore\n" );
	 }
	else
	 {
      		fprintf( ps_file, "21 %3d moveto\n", *ps_line_number-3 );
      		fprintf( ps_file, "(%4d) 9.0 Print\n\n", hb_number);
    		fprintf(ps_file, "%3d 44 moveto\n", *ps_line_number+4 );
		fprintf(ps_file, "(%6.3f) 9.0 Print\n", hb_energy );
      		fprintf( ps_file, "72 %3d moveto\n", *ps_line_number-3 );
      		fprintf( ps_file, "(%7.3f) 9.0 Print\n\n", mean_coordination);
	 }
    }
*/
  }

  for( a = 0; a <= number_of_chains; a++ ){

    fprintf(ps_file,"103 %3d moveto\n", *ps_line_number-((10*a)+3) );
    fprintf(ps_file,"(%c) 9.0 Print\n", chain_IDs[a] );

    /************************************************************/
    /* Print the red and blue triangles that denote where the   */
    /* broken Hbond occured, to the postscript file.            */
    /************************************************************/
    if( !first_line ) {
      if( donor_chain_ID == chain_IDs[a] ) {
	this_chain = a;
	triangle_x = ( HB_start * 3 ) + 111.5;
	triangle_y = *ps_line_number-(3+(10*a));

	fprintf(ps_file,"Blue\n");
	fprintf(ps_file,"%3d %3d %3d %3d %3d %3d Pl3\n", triangle_x, triangle_y, (triangle_x - 2),
		(triangle_y - 3), (triangle_x + 2), (triangle_y - 3) );
      }
      if( accept_chain_ID == chain_IDs[a] ){

	triangle_x = ( HB_end * 3 ) + 111.5;
	triangle_y = *ps_line_number-(3+(10*a));

	fprintf(ps_file,"Red\n");
	fprintf(ps_file,"%3d %3d %3d %3d %3d %3d Pl3\n", triangle_x, triangle_y, (triangle_x - 2),
		(triangle_y - 3), (triangle_x + 2), (triangle_y - 3) );
      }
    }
    /************************************************************/
  }

  HB_end   += (y_translate[this_chain] - 1);
  HB_start += (y_translate[this_chain] - 1);

  if( !first_line && ( short_output || missing_one || fit_it_all_on_one_page ) ){
    fprintf(ps_file,"Blue\n");

    if( donor_chain_ID == 'W' ) {
      fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+112, *ps_line_number-4 );
      fprintf(ps_file,"(W) 9.0 Print\n");
    }
    else{
      if( donor_chain_ID == 'H' ) {
	fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+112, *ps_line_number-4 );
	fprintf(ps_file,"(H) 9.0 Print\n");
      }
      else{
	if( !strcmp(donor_atom_type, "N" ) ||
	    !strcmp(donor_atom_type, "O" ) ) {
	  fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+112, *ps_line_number-4 );
	  fprintf(ps_file,"(M) 9.0 Print\n");
	}
	else{
	  fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+112, *ps_line_number-4 );
	  fprintf(ps_file,"(S) 9.0 Print\n");
	}
      }
    }

    fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+120, *ps_line_number-4 );
    fprintf(ps_file,"(%3d) 9.0 Print\n", HB_start);

    if( donor_chain_ID != 'W' && donor_chain_ID != 'H' ){
      fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+131, *ps_line_number-4 );
      fprintf(ps_file,"(%c) 9.0 Print\n", donor_chain_ID );
    }

    fprintf(ps_file,"Red\n");

    if( accept_chain_ID == 'W' ) {
      fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+140, *ps_line_number-4 );
      fprintf(ps_file,"(W) 9.0 Print\n");
    }
    else{
      if( accept_chain_ID == 'H' ) {
	fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+140, *ps_line_number-4 );
	fprintf(ps_file,"(H) 9.0 Print\n");
      }
      else{
	if( !strcmp(accept_atom_type, "N" ) ||
	    !strcmp(accept_atom_type, "O" ) ) {
	  fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+140, *ps_line_number-4 );
	  fprintf(ps_file,"(M) 9.0 Print\n");
	}
	else{
	  fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+140, *ps_line_number-4 );
	  fprintf(ps_file,"(S) 9.0 Print\n");
	}
      }
    }

    fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+148, *ps_line_number-4 );
    fprintf(ps_file,"(%3d) 9.0 Print\n", HB_end );

    if( accept_chain_ID != 'W' && accept_chain_ID != 'H' ){
      fprintf(ps_file,"%d %d moveto\n", (NUMBER_OF_RESIDUES*3)+162, *ps_line_number-4 );
      fprintf(ps_file,"(%c) 9.0 Print\n", accept_chain_ID );
    }

  }

  if( long_output )
    *ps_line_number -= 60;
  else
    *ps_line_number -= 10;

  if( short_output || missing_one ) {
    if(*ps_line_number < 60 || ( *ps_line_number - ( number_of_chains*10 ) ) < 60){
      *ps_line_number = 710;
      (*ps_page_number)++;
      fprintf(ps_file, "\nshowpage\n" );
      fprintf(ps_file,"%%%%Page:    %d   %d\n\n", *ps_page_number,
	      *ps_page_number );
    }
  }

}

/*************************************************************************/
/* If the protein has more than 125 amino acids, the output is displayed */
/* in landscape format. 250 amino-acids will fit on a page, multiple land*/
/* scape pages are used if neccesary.                                    */
/*************************************************************************/
void print_current_Hbond_info_landscape( 
			FILE *ps_file, int long_output, int short_output,
			int first_line, int *ps_line_number, int hb_number,
			float hb_energy, int cluster_counter, 
			int colors[MAX_CLUSTER_COUNT],
			int HB_start, int HB_end, int NUMBER_OF_RESIDUES,
			float scale_factor, int *ps_page_number, 
			int missing_one,
			int HBtwo_start, int HBtwo_end, 
			int y_translate[MAX_CHAIN_COUNT],
			char donor_atom_type[4], char accept_atom_type[4],
			int number_of_chains, char donor_chain_ID,
			char accept_chain_ID, char chain_IDs[MAX_CHAIN_COUNT], 
			int chain_size[MAX_CHAIN_COUNT], char protein_name[128], 
			int side_chain_test, float mean_coordination, 
			int fit_it_all_on_one_page,
			/* SN 2008:01 */
			int insertions_prior_donor,int insertions_prior_acptr,
			char donor_inscode, char acptr_inscode, 
			/* SN 2008:05 */
			int total_hbonds){

  int
    a = 0,
    b = 0,
    triangle_x = 0.0,
    summed_length = 0,
    end_pt = 0,
    start_pt = 0,
    beggining_of_current_chain = 0,
    acptr_chain = 0,donor_chain = 0;	/* SN 2008:02 
	These replace single variable this_chain, which results in error if 
	inter-chain H-bonds are present in the input, e.g., ligand to protein */

  float
    y_start=0.0,
    y_end=0.0,
    triangle_y = 0.0;

  fprintf(ps_file, "Col00\n");

  /*
   * 2007:05	SN
   * 	Print the "All 'x' Hbonds" phrase if the energy cutoff
   *    is 0 kcal/mol, i.e., the weakest Hbond would have energy
   *    equal to zero.
   */
/*printf("\n#$@#@   \t %f",hb_energy);
  if( first_line && hb_energy == 0.0){
    fprintf(ps_file, "%3d 36 moveto\n", *ps_line_number+3 );
    fprintf(ps_file, "(All %3d Hbonds) 9.0 UncenterRot90 Rot90\n", hb_number);
    fprintf(ps_file, "(All %3d Hbonds) 9.0 Print\n", hb_number);
    fprintf(ps_file, "grestore\n" );
  }
*/
  if( first_line ) {
    fprintf(ps_file, "%3d 21 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%4d) 9.0 UncenterRot90 Rot90\n", total_hbonds-hb_number);
    fprintf(ps_file, "(%4d) 9.0 Print\n", total_hbonds - hb_number );
    fprintf(ps_file, "grestore\n" );
    fprintf(ps_file, "%3d 44 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", hb_energy );
    fprintf(ps_file, "(%6.3f) 9.0 Print\n", hb_energy );
    fprintf(ps_file, "grestore\n" );
    fprintf(ps_file, "%3d 72 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", mean_coordination );
    fprintf(ps_file, "(%6.3f) 9.0 Print\n", mean_coordination );
    fprintf(ps_file, "grestore\n" );
  

    for( a = 0; a <= number_of_chains; a++ ){

      summed_length += ( (chain_size[a]*3) + 18 ); 
			/* end point for the current chain on the output file */

      for( b = 0; b < a; b++ )
	beggining_of_current_chain += (chain_size[b]*3) + 18;

      beggining_of_current_chain += 110;

      fprintf(ps_file, "Col00\n");
      fprintf(ps_file, "%3d %3d moveto\n", *ps_line_number+3, beggining_of_current_chain - 8 );
      fprintf(ps_file, "(%c) 9.0 UncenterRot90 Rot90\n", chain_IDs[a] );
      fprintf(ps_file, "(%c) 9.0 Print\n", chain_IDs[a] );
      fprintf(ps_file, "grestore\n" );

      beggining_of_current_chain = 0;
    }
  }
  else
  {
    fprintf(ps_file, "%3d 21 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%4d) 9.0 UncenterRot90 Rot90\n", total_hbonds-hb_number);
    fprintf(ps_file, "(%4d) 9.0 Print\n", total_hbonds - hb_number );
    fprintf(ps_file, "grestore\n" );
    fprintf(ps_file, "%3d 44 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", hb_energy );
    fprintf(ps_file, "(%6.3f) 9.0 Print\n", hb_energy );
    fprintf(ps_file, "grestore\n" );
    fprintf(ps_file, "%3d 72 moveto\n", *ps_line_number+4 );
    fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", mean_coordination );
    fprintf(ps_file, "(%6.3f) 9.0 Print\n", mean_coordination );
    fprintf(ps_file, "grestore\n" );
  

    for( a = 0; a <= number_of_chains; a++ ){

      summed_length += ( (chain_size[a]*3) + 18 ); 
			/* end point for the current chain on the output file */

      for( b = 0; b < a; b++ )
	beggining_of_current_chain += (chain_size[b]*3) + 18;

      beggining_of_current_chain += 110;

      fprintf(ps_file, "Col00\n");
      fprintf(ps_file, "%3d %3d moveto\n", *ps_line_number+3, beggining_of_current_chain - 8 );
      fprintf(ps_file, "(%c) 9.0 UncenterRot90 Rot90\n", chain_IDs[a] );
      fprintf(ps_file, "(%c) 9.0 Print\n", chain_IDs[a] );
      fprintf(ps_file, "grestore\n" );

      /************************************************************/
      /* Print the red and blue triangles that denote where the   */
      /* broken Hbond occured, to the postscript file.            */
      /************************************************************/
      if( donor_chain_ID == chain_IDs[a] ) {

	triangle_x = *ps_line_number + 3;
	triangle_y = ( (HB_start-1) * 3 ) + 1.5 + beggining_of_current_chain;

	fprintf(ps_file, "Blue\n");
	fprintf(ps_file, "%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n",\
		 triangle_x, triangle_y, (triangle_x + 3),
		(triangle_y - 2), (triangle_x + 3), (triangle_y + 2) );
        donor_chain = a;	/* SN 2008:02 */
      }

      if( accept_chain_ID == chain_IDs[a] ){

	if( HB_start == HB_end ){
	  triangle_x = *ps_line_number - 3;
	  triangle_y = ( (HB_end-1) * 3 ) + 1.5 + beggining_of_current_chain;

	  fprintf(ps_file,"Red\n");
	  fprintf(ps_file,"%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", triangle_x,\
		triangle_y,(triangle_x - 3),(triangle_y - 2),(triangle_x - 3),\
		(triangle_y + 2) );
	}
	else{
	  triangle_x = *ps_line_number + 3;
	  triangle_y = ( (HB_end-1) * 3 ) + 1.5 + beggining_of_current_chain;

	  fprintf(ps_file,"Red\n");
	  fprintf(ps_file,"%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", triangle_x, \
		triangle_y, (triangle_x + 3),(triangle_y - 2),(triangle_x + 3),\
		(triangle_y + 2) );
	}
	acptr_chain = a;	/* SN 2008:02 	--- replace "this_chain" */
      }
      /************************************************************/
      beggining_of_current_chain = 0;
    }

#ifdef DEBUG_INTER_CHAIN_HBOND
	printf("\nDEBUG: chain - %d, chain size: %d, chain start: %d,HB_start: %d, chain start: %d\n",\
	donor_chain,chain_size[donor_chain],beggining_of_current_chain,
	HB_start,y_translate[donor_chain]);
#endif

    HB_start += (y_translate[donor_chain] - 1);
    HB_end   += (y_translate[acptr_chain] - 1);

    end_pt = summed_length-15;

#ifdef DEBUG_SN
	printf("\nHBOND #: %d, first_line: %d",hb_number,first_line);
#endif

    if( !first_line && ( short_output || fit_it_all_on_one_page ) ){ 
    /*if( short_output || fit_it_all_on_one_page ) {*/
      fprintf(ps_file,"Blue\n");

      if( donor_chain_ID == 'W' ) {
	fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+110 );
	fprintf(ps_file,"(W) 9.0 UncenterRot90 Rot90\n" );
	fprintf(ps_file,"(W) 9.0 Print\n");
	fprintf(ps_file,"grestore\n" );
      }
      else{

      /* SN 2008:01	Erroneous! Chain_id == 'H' does not necessarily
			signify HETATM involving bond
	if( donor_chain_ID == 'H' ) {
	  fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+110 );
	  fprintf(ps_file,"(H) 9.0 UncenterRot90 Rot90\n" );
	  fprintf(ps_file,"(H) 9.0 Print\n");
	  fprintf(ps_file,"grestore\n" );
	} 
	else{
      */
	  if( !strcmp(donor_atom_type, "N" ) ||
	      !strcmp(donor_atom_type, "O" ) ) {
	    fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+110 );
	    fprintf(ps_file,"(M) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(ps_file,"(M) 9.0 Print\n");
	    fprintf(ps_file,"grestore\n" );
	  }
	  else{
	    fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+110 );
	    fprintf(ps_file,"(S) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(ps_file,"(S) 9.0 Print\n");
	    fprintf(ps_file,"grestore\n" );
	  }
/*	}*/
      }

#ifdef DEBUG_SN
	printf("\nDEBUG (postscript.c 671): %d\t%d\t%d",\
			hb_number,HB_start,insertions_prior_donor);
#endif

       /* 2008:01 SN	Change the order of HBOND info to the right side 
	* 		of the plot and account for the insertions!
	*		Old positions: 120 129 144      155 164 177
	*		New positions: 110 122 129 143  152 164 171 185
	*
	* <H|M|S|W> <rsd> <chain-id> ->	<H|M|S|W> <chain-id> <rsd> <ins code>
	*/
			/* chain-id */			

      if( donor_chain_ID != 'W' ){
	fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+122 );
	fprintf(ps_file,"(%c) 9.0 UncenterRot90 Rot90\n", donor_chain_ID );
	fprintf(ps_file,"(%c) 9.0 Print\n", donor_chain_ID );
	fprintf(ps_file,"grestore\n" );
      }

			/* RSD */		
	
      fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+129 );
      if( HB_start - insertions_prior_donor < 100 )
      {
        fprintf(ps_file,"(%4d) 9.0 UncenterRot90 Rot90\n", \
					HB_start - insertions_prior_donor );
        fprintf(ps_file,"(%4d) 9.0 Print\n", \
					HB_start - insertions_prior_donor );
      }
      else
      {
        fprintf(ps_file,"(%3d) 9.0 UncenterRot90 Rot90\n", \
					HB_start - insertions_prior_donor );
        fprintf(ps_file,"(%3d) 9.0 Print\n", \
					HB_start - insertions_prior_donor );
      }
      fprintf(ps_file,"grestore\n" );

			/* Insertion Code */

      fprintf(ps_file,"%d %d moveto\n",*ps_line_number+3, (end_pt+144 ));
      fprintf(ps_file,"(%c) 9.0 UncenterRot90 Rot90\n", donor_inscode );
      fprintf(ps_file,"(%c) 9.0 Print\n", donor_inscode );
      fprintf(ps_file,"grestore\n" );


			/* Acceptor Info. */		

      fprintf(ps_file,"Red\n");
      if( accept_chain_ID == 'W' ) {
	fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+152 );
	fprintf(ps_file,"(W) 9.0 UncenterRot90 Rot90\n" );
	fprintf(ps_file,"(W) 9.0 Print\n");
	fprintf(ps_file,"grestore\n" );
      }
      else{
/*	if( accept_chain_ID == 'H' ) {
	  fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+152 );
	  fprintf(ps_file,"(H) 9.0 UncenterRot90 Rot90\n" );
	  fprintf(ps_file,"(H) 9.0 Print\n");
	  fprintf(ps_file,"grestore\n" );
	}
	else{
*/
	  if( !strcmp(accept_atom_type, "N" ) ||
	      !strcmp(accept_atom_type, "O" ) ) {
	    fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+152 );
	    fprintf(ps_file,"(M) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(ps_file,"(M) 9.0 Print\n");
	    fprintf(ps_file,"grestore\n" );
	  }
	  else{
	    fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+152 );
	    fprintf(ps_file,"(S) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(ps_file,"(S) 9.0 Print\n");
	    fprintf(ps_file,"grestore\n" );
	  }
/*	}*/
      }

			/* Chain-id */	
	
      if( accept_chain_ID != 'W' ){
	fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+164 );
	fprintf(ps_file,"(%c) 9.0 UncenterRot90 Rot90\n", accept_chain_ID );
	fprintf(ps_file,"(%c) 9.0 Print\n", accept_chain_ID );
	fprintf(ps_file,"grestore\n" );
      }

			/* RSD */	
	
      fprintf(ps_file,"%d %d moveto\n", *ps_line_number+3, (end_pt)+171 );
      if( HB_end - insertions_prior_acptr < 100 )
      {
        fprintf(ps_file,"(%4d) 9.0 UncenterRot90 Rot90\n", \
					HB_end - insertions_prior_acptr);
        fprintf(ps_file,"(%4d) 9.0 Print\n", \
					HB_end - insertions_prior_acptr);
      }
      else
      {
        fprintf(ps_file,"(%3d) 9.0 UncenterRot90 Rot90\n", \
					HB_end - insertions_prior_acptr);
        fprintf(ps_file,"(%3d) 9.0 Print\n", \
					HB_end - insertions_prior_acptr);
      }
      fprintf(ps_file,"grestore\n" );


			/* Insertion Code */

      fprintf(ps_file,"%d %d moveto\n",*ps_line_number+3, (end_pt+186 ));
      fprintf(ps_file,"(%c) 9.0 UncenterRot90 Rot90\n", acptr_inscode );
      fprintf(ps_file,"(%c) 9.0 Print\n", acptr_inscode );
      fprintf(ps_file,"grestore\n" );


      if( side_chain_test ){
	fprintf(ps_file, "%3d 90 moveto\n", *ps_line_number+4 );
	fprintf(ps_file, "(%4d) 9.0 UncenterRot90 Rot90\n", \
						total_hbonds - hb_number );
	fprintf(ps_file, "(%4d) 9.0 Print\n",total_hbonds - hb_number );
	fprintf(ps_file, "grestore\n" );
	fprintf(ps_file, "%3d 60 moveto\n", *ps_line_number+4 );
	fprintf(ps_file, "((%4d)) 9.0 UncenterRot90 Rot90\n", HB_start );
	fprintf(ps_file, "((%4d)) 9.0 Print\n", HB_start );
	fprintf(ps_file, "grestore\n" );
	fprintf(ps_file, "%3d 30 moveto\n", *ps_line_number+4 );
	fprintf(ps_file, "(%4s) 9.0 UncenterRot90 Rot90\n", donor_atom_type );
	fprintf(ps_file, "(%4s) 9.0 Print\n", donor_atom_type );
	fprintf(ps_file, "grestore\n" );
      }

    }
  }

  if( long_output )
    *ps_line_number += 60;
  if( side_chain_test )
    *ps_line_number += 8  + ( ( 10 * number_of_chains ) );
  else
    *ps_line_number += 12;


  /* comment out the following code if you want all the output to be on a single page */
  /*
    if( short_output || missing_one || side_chain_test ) {
    if(*ps_line_number > 560 || ( ( *ps_line_number + ( number_of_chains * 10 ) ) > 560 ) ) {
    *ps_line_number = 60;
    (*ps_page_number)++;
    fprintf(ps_file, "\nshowpage\n" );
    fprintf(ps_file,"%%%%Page:    %d   %d\n\n", *ps_page_number,
    *ps_page_number );
    }
    }
  */

}
/**********************************************************************/

/**********************************************************************/
/* prints the info in landscape format, but on the right page. Allows */
/* for multipage output for proteins with more than 200 amino-acids.  */
/**********************************************************************/
void print_current_Hbond_info_landscape_multipage( 
		FILE *ps_file, int long_output, int short_output,
		int first_line, int *ps_line_number, int hb_number,
		float hb_energy, int cluster_counter, 
		int colors[MAX_CLUSTER_COUNT],
		int donor_res_num, int accept_res_num,  int NUMBER_OF_RESIDUES,
		float scale_factor, int *ps_page_number, int missing_one,
		int HBtwo_start, int HBtwo_end, 
		int y_translate[MAX_CHAIN_COUNT],
		char donor_atom_type[4], char accept_atom_type[4],
		FILE *file_list[FILE_LIST_SIZE], int number_of_pages,
		int number_of_chains, char donor_chain_ID,
		char accept_chain_ID, char chain_IDs[MAX_CHAIN_COUNT],
		char protein_name[128], int chain_size[MAX_CHAIN_COUNT],
		float mean_coordination,
		/* SN 2008:01 */
		int insertions_prior_donor,int insertions_prior_acptr,
		char donor_inscode, char acptr_inscode,
		/* SN 2008:05 */
		int total_hbonds){

  int
    a = 0,
    b = 0,
    last_page = 0,
    end_point = 0,
    output_page = 0,
    triangle_x = 0,
    end_pt=0,
    start_pt = 0,
    current_pg = 0,
    summed_length = 0,
    previous_length = 0,
    donor_pointer = 0,
    accept_pointer = 0;

  float
    y_start = 0.0,
    y_end = 0.0,
    y_position = 0.0,
    triangle_y = 0.0;

  char
    linebuf[300],
    file_name[6];

  FILE
    *current_file;


  /**********************************************************************/
  /* Print first 3 columns of output (hb number, energy, mean coord).   */

  fprintf(ps_file, "Col00\n");
  fprintf(ps_file, "%3d 21 moveto\n", *ps_line_number+4 );
  fprintf(ps_file, "(%4d) 9.0 UncenterRot90 Rot90\n", total_hbonds-hb_number );
  fprintf(ps_file, "(%4d) 9.0 Print\n", total_hbonds-hb_number );
  fprintf(ps_file, "grestore\n" );
  fprintf(ps_file, "%3d 44 moveto\n", *ps_line_number+4 );
  fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", hb_energy );
  fprintf(ps_file, "(%6.3f) 9.0 Print\n", hb_energy );
  fprintf(ps_file, "grestore\n" );
  fprintf(ps_file, "%3d 72 moveto\n", *ps_line_number+4 );
  fprintf(ps_file, "(%6.3f) 9.0 UncenterRot90 Rot90\n", mean_coordination );
  fprintf(ps_file, "(%6.3f) 9.0 Print\n", mean_coordination );
  fprintf(ps_file, "grestore\n" );
  
  /**********************************************************************/

  /**********************************************************************/
  /* Print the chain ID before each chain.                              */
  for( b = 0; b <= number_of_chains; b++ ){

    summed_length += ( chain_size[b] + ( b*6 ) ); /* Adjusted in this line+13 */
    start_pt = ( summed_length -chain_size[b] );
    current_pg = start_pt/200;
    current_file = file_list[current_pg];

    while( start_pt > 200 )
      start_pt -= 200;

    fprintf(current_file, "Col00\n");
    fprintf(current_file, "%3d %3d moveto\n", *ps_line_number+3, (start_pt*3) +110 -9 );
    fprintf(current_file, "(%c) 9.0 UncenterRot90 Rot90\n", chain_IDs[b] );
    fprintf(current_file, "(%c) 9.0 Print\n", chain_IDs[b] );
    fprintf(current_file, "grestore\n" );
    summed_length -= ( b*6 );

  }
  /**********************************************************************/

  summed_length = 0;

	/* CORRECT VALUES FOR *_res_num CAN BE INPUT FROM hbdilute.c */
	/* SIMILAR TO THE WAY THEY ARE INPUT TO print_decomp()	     */

  /**********************************************************************/
  /* ProFlex outputs the donor and acceptor residue numbers as starting */
  /* from 1, regardless of the res number the chains actually start from*/
  /* Modify these numbers to coincide with the real residue number.     */
  for( b = 0; b <= number_of_chains; b++ ){
    if( donor_chain_ID == chain_IDs[b] )
      donor_res_num  += y_translate[b] -1;
    if( accept_chain_ID == chain_IDs[b] )
      accept_res_num += y_translate[b] -1;
  }
  /**********************************************************************/

  /**********************************************************************/
  /* Print the red and blue triangles under each line showing the loc-  */
  /* ation of the hydrogen bond along the primary structure.            */
  if( !first_line ) {

    for( b = 0; b <= number_of_chains; b++ ){

/*      summed_length += (chain_size[b] + ( 6*b )); --- FAULTY!!  SN 02:08  */
	summed_length += (chain_size[b] + 6*(b>0?1:0) ); 
			/* SN 02:08	Add 6 to account for inter-chain 
			 spacing including the chain-ID (non-zero IDs only) */

      /**********************************************************************/
      /* The second condition will stop hetatm's from having triangles.     */
      /**********************************************************************/
      if( donor_chain_ID == chain_IDs[b] &&
	  donor_res_num <= chain_size[b] ) {

	output_page = ( summed_length -chain_size[b] + donor_res_num \
				- y_translate[b]) / 200.0;
	current_file = file_list[output_page];

	previous_length = summed_length - chain_size[b];

	while( previous_length > 200 ) {
	  previous_length -= 200;
	}

	donor_pointer = donor_res_num -y_translate[b];
	while( donor_pointer >= 200 ) {
	  donor_pointer -= 200;
	}

	if( ( 200 - previous_length ) > ( chain_size[b] ) )
	  y_position = ( (donor_pointer +previous_length)*3) + 111.5;
	else {
	  if( donor_pointer >= ( 200-previous_length) )
	    y_position = ( (donor_pointer -(200 -previous_length))*3 ) + 111.5;
	  else
	    y_position = ( (donor_pointer +previous_length)*3) + 111.5;
	}

	triangle_x = *ps_line_number + 3;
	triangle_y = y_position;

	fprintf(current_file, "Blue\n");
	fprintf(current_file, "%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", 
		triangle_x, triangle_y, (triangle_x + 3),
		(triangle_y - 2), (triangle_x + 3), (triangle_y + 2) );

      }
      if( accept_chain_ID == chain_IDs[b] &&
	  accept_res_num <= chain_size[b] ) {

	output_page = ((summed_length -chain_size[b] +accept_res_num - \
						y_translate[b])*1.0)/200.0;
	current_file = file_list[output_page];

    /* Do not print unnecessary debug comments - Sameer Arora, Jan 4 2004 */
#ifdef DEBUG
	printf("res num %3d output_page %d, summed %d chain_size %d y_trans %d\n", 
		accept_res_num, output_page, summed_length, chain_size[b], 
		y_translate[b] );
#endif

	previous_length = summed_length - chain_size[b];

	while( previous_length > 200 ) {
	  previous_length -= 200;
	}

	accept_pointer = accept_res_num -y_translate[b];
	while( accept_pointer >= 200 ) {
	  accept_pointer -= 200;
	}

	if( ( 200 - previous_length ) > ( chain_size[b] ) )
	  y_position = ( (accept_pointer +previous_length) *3) + 111.5;
	else {
	  if( accept_pointer >= ( 200-previous_length) )
	    y_position = ( (accept_pointer -(200-previous_length)) *3) + 111.5;
	  else
	    y_position = ( (accept_pointer +previous_length) *3) + 111.5;
	}
	triangle_x = *ps_line_number + 3;
	triangle_y = y_position;

	fprintf(current_file, "Red\n");
	fprintf(current_file, "%3d %5.2f %3d %5.2f %3d %5.2f Pl3\n", triangle_x,
		 triangle_y, (triangle_x + 3),(triangle_y - 2),(triangle_x + 3),
		 (triangle_y + 2) );

      }
    }
  }
  /**********************************************************************/

  last_page = summed_length / 200;
  current_file = file_list[last_page];
/*
  if( number_of_chains ){
*/
    while( summed_length > 200 ) {
      summed_length -= 200;
    }
    end_point = summed_length + 2;	/* Added 2 for formatting */
/*
  }
  else
    end_point = NUMBER_OF_RESIDUES - (last_page * 200);
*/

  if( !first_line && short_output ){

    fprintf(current_file, "Blue\n");

      if( donor_chain_ID == 'W' ) {
	fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+110 );
	fprintf(current_file,"(W) 9.0 UncenterRot90 Rot90\n" );
	fprintf(current_file,"(W) 9.0 Print\n");
	fprintf(current_file,"grestore\n" );
      }
      /* SN 2008:01	Erroneous! Chain_id == 'H' does not necessarily
			signify HETATM involving bond
      else{
	if( donor_chain_ID == 'H' ) {
	  fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+110 );
	  fprintf(current_file,"(H) 9.0 UncenterRot90 Rot90\n" );
	  fprintf(current_file,"(H) 9.0 Print\n");
	  fprintf(current_file,"grestore\n" );
	}
      */
      else{
	  if( !strcmp(donor_atom_type, "N" ) ||
	      !strcmp(donor_atom_type, "O" ) ) {
	    fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+110 );
	    fprintf(current_file,"(M) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(current_file,"(M) 9.0 Print\n");
	    fprintf(current_file,"grestore\n" );
	  }
	  else{
	    fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+110 );
	    fprintf(current_file,"(S) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(current_file,"(S) 9.0 Print\n");
	    fprintf(current_file,"grestore\n" );
	  }
      }

      /* 2008:01 SN	Change the order of HBOND info to the right side 
       * 		of the plot and account for the insertions!
       *		Old positions: 120 129 144      155 164 177
       *		New positions: 110 122 129 143  152 164 171 185
       *
       * <H|M|S|W> <rsd> <chain-id> ->	<H|M|S|W> <chain-id> <rsd> <ins code>
       */

      /* chain-id */		

      if( donor_chain_ID != 'W' ){
        fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+122 );
        fprintf(current_file,"(%c) 9.0 UncenterRot90 Rot90\n", donor_chain_ID );
	fprintf(current_file,"(%c) 9.0 Print\n", donor_chain_ID );
	fprintf(current_file,"grestore\n" );
      }

      /* RSD */		
	
      fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+129 );
      if( donor_res_num - insertions_prior_donor < 100 )
      {
        fprintf(current_file,"(%4d) 9.0 UncenterRot90 Rot90\n", \
			donor_res_num - insertions_prior_donor );
        fprintf(current_file,"(%4d) 9.0 Print\n", \
			donor_res_num - insertions_prior_donor );
      }
      else
      {
        fprintf(current_file,"(%3d) 9.0 UncenterRot90 Rot90\n", \
			donor_res_num - insertions_prior_donor );
        fprintf(current_file,"(%3d) 9.0 Print\n", \
			donor_res_num - insertions_prior_donor );
      }
      fprintf(current_file,"grestore\n" );

      /* Insertion Code */
	
      fprintf(current_file,"%d %d moveto\n",*ps_line_number+3, (end_point*3)+144 );
      fprintf(current_file,"(%c) 9.0 UncenterRot90 Rot90\n", donor_inscode );
      fprintf(current_file,"(%c) 9.0 Print\n", donor_inscode );
      fprintf(current_file,"grestore\n" );


      /* Acceptor Info. */	

      fprintf(current_file, "Red\n");

      if( accept_chain_ID == 'W' ) {
	fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+152 );
	fprintf(current_file,"(W) 9.0 UncenterRot90 Rot90\n" );
	fprintf(current_file,"(W) 9.0 Print\n");
	fprintf(current_file,"grestore\n" );
      }
      /* SN 2008:01	Erroneous! Chain_id == 'H' does not necessarily
			signify HETATM involving bond
      else{
	if( accept_chain_ID == 'H' ) {
	  fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, (end_point*3)+155 );
	  fprintf(current_file,"(H) 9.0 UncenterRot90 Rot90\n" );
	  fprintf(current_file,"(H) 9.0 Print\n");
	  fprintf(current_file,"grestore\n" );
	}
      */
      else{
	  if( !strcmp(accept_atom_type, "N" ) ||
	      !strcmp(accept_atom_type, "O" ) ) {
	    fprintf(current_file,"%d %d moveto\n",*ps_line_number+3, \
							(end_point*3)+152 );
	    fprintf(current_file,"(M) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(current_file,"(M) 9.0 Print\n");
	    fprintf(current_file,"grestore\n" );
	  }
	  else{
	    fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, \
							(end_point*3)+152 );
	    fprintf(current_file,"(S) 9.0 UncenterRot90 Rot90\n" );
	    fprintf(current_file,"(S) 9.0 Print\n");
	    fprintf(current_file,"grestore\n" );
	  }
      }

      /* chain-id */	
	
      if( accept_chain_ID != 'W' ){
	fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, \
							(end_point*3)+164 );
	fprintf(current_file,"(%c) 9.0 UncenterRot90 Rot90\n",accept_chain_ID);
	fprintf(current_file,"(%c) 9.0 Print\n", accept_chain_ID );
	fprintf(current_file,"grestore\n" );
      }

      /* RSD */	
	
      fprintf(current_file,"%d %d moveto\n", *ps_line_number+3, \
							(end_point*3)+171 );
      if(accept_res_num - insertions_prior_acptr < 100)
      {
        fprintf(current_file,"(%4d) 9.0 UncenterRot90 Rot90\n",\
			accept_res_num - insertions_prior_acptr );
        fprintf(current_file,"(%4d) 9.0 Print\n", \
			accept_res_num - insertions_prior_acptr );
      }
      else
      {
        fprintf(current_file,"(%3d) 9.0 UncenterRot90 Rot90\n",\
			accept_res_num - insertions_prior_acptr );
        fprintf(current_file,"(%3d) 9.0 Print\n", \
			accept_res_num - insertions_prior_acptr );
      }

      fprintf(current_file,"grestore\n" );

      /* Insertion Code */

      fprintf(current_file,"%d %d moveto\n",*ps_line_number+3, (end_point*3)+186 );
      fprintf(current_file,"(%c) 9.0 UncenterRot90 Rot90\n", acptr_inscode );
      fprintf(current_file,"(%c) 9.0 Print\n", acptr_inscode );
      fprintf(current_file,"grestore\n" );


    fprintf(current_file, "Col00\n");
  }
  /**********************************************************************/

  /**********************************************************************/
  /* Incremenet the line number. Also, if we just printed the last line */
  /* that will fit on this page, reset all variables to print at the top*/
  /* of the page, and print appropriate page markers to all postscript  */
  /* files.                                                             */
  *ps_line_number += 12;

  if(*ps_line_number > 560) {
    *ps_line_number = 60;
    (*ps_page_number)++;

    for( a = 1; a <= number_of_pages; a++ ){
      current_file = file_list[a];
      fprintf(current_file, "\nshowpage\n" );
      fflush( current_file );
      fclose( current_file );
    }
    fprintf(ps_file, "\nshowpage\n" );
    fflush( ps_file );

    for( a = 1; a <= number_of_pages; a++ ) {
      sprintf( file_name, "%d", a );
      current_file = fopen( file_name, "r" );
      while( fgets( linebuf, sizeof(linebuf), current_file ) != NULL ){
	fprintf(ps_file, "%s", linebuf );
      fflush(ps_file);
      }
    }

    fprintf(ps_file, "%%%%Page:    %d   %d\n", *ps_page_number,
	    *ps_page_number );
  }
  /**********************************************************************/

}
/**********************************************************************/

/**********************************************************************/
/* Prints the actual colorbars corresponding to the the rigid cluster */
/* decomposition. The information for each cluster in included as a   */
/* linked list. The lists for each cluster are indexed from the array */
/* nc_count[]. This routine is for single page output only (full size */
/* portrait, full size landscape, scaled "fit-it-all-on-one-page").   */
/**********************************************************************/
void print_decomp( 	clusters *nc_count[MAX_CLUSTER_COUNT], 
			int cluster_counter,
		   	int ps_line_number, FILE *ps_file, 
		   	int colors[MAX_CLUSTER_COUNT],
		   	int y_translate[MAX_CHAIN_COUNT], 
			int number_of_residues,
		   	int number_of_chains, int chain_size[MAX_CHAIN_COUNT],
			char chain_IDs[MAX_CHAIN_COUNT],
		   	int total_insertions, 
			int insertion_res_num[100],
		   	int number_of_insertions[100], 
			char insertion_chain_ID[100],
		   /* SN 2008:01 */
		   int min_max_rsd[MAX_CHAIN_COUNT][2], 
		   int missing_rsd_list[MISSING_RSD_LIST_SIZE],
		   int missing_rsd_count[MAX_CHAIN_COUNT]) {

  int
    a = 0,
    b = 0,
    c = 0,
    d = 0,
    e = 0,
    this_chain = 0,
    top = 0,
    bottom = 0,
    y_position = 0,
    start_pt = 0,
    end_pt = 0,
    previous_chains_size = 0,
    next_residue = 0,
    print_this_number = 0,
    print_number = 1,
    next_rsd_index,cumulative_insertions=0; /* 2008:01	SN */

  clusters
    *nc_current = NULL;

  char ch;

  /*************************************************************/
  /* First need to compute how the output will fit on the page */
  /* If there are <= 125 residues, use the portrait format. If */
  /* it's 125 < residues >= 250, use landscape. Greater than   */
  /* 250 residues will require multiple page printouts, there  */
  /* is an option to print it all on one page.                 */
  /* In the output, the height of the color bars is always the */
  /* same, only the width should change.                       */
  /*************************************************************/
  /*if( (  number_of_chains && ( (number_of_residues + (number_of_chains*6) ) <= 125 ) ) ||
    ( !number_of_chains && (  number_of_residues <= 125 ) ) ) {*/

  /**********************************************************************/
  /* print the decomp in single page landscape for proteins with        */
  /* 125 > amino-acids <= 175. The following code is also executed when */
  /* hbdilute.c is run with the "b" option for, fit_it_all_on_one_page. */
  /**********************************************************************/
  /* else if rsds > 125 ; previous version handled the case of 
     rsdcount < 125 differently */
  {

    start_pt = 110;

    for( a = 0; a <= number_of_chains; a++ ){

      /* 2008:01 SN	chain_size[] reflects all the residues missing 
			as well as inserted */
      end_pt = ( chain_size[a] * 3 ) + start_pt;

      /************************************************************/
      /* Draw a black line representing the primary structure.    */

      /* 2008:01 SN	Missing residues can be marked here by leaving
       * 		gaps in the black line
       */
      fprintf(ps_file, "\nCol00\n");
      fprintf(ps_file, "%3d %3d %3d %3d sequence_line\n", ps_line_number,
	      start_pt, ps_line_number, end_pt );
      fprintf(ps_file, "grestore\n");
      /************************************************************/

      /************************************************************/
      /* For landscape output, the top line of the page x=60. If  */
      /* ps_line_number == 60, print the residue numbers for each */
      /* chain.                                                   */
      if( ps_line_number == 60 ) {
        if(y_translate[a] >= 1000 && y_translate[a] < 10000)
	  fprintf(ps_file, "26 %d  moveto\n", start_pt -1 );
	else if(y_translate[a] >= 100 && y_translate[a] < 1000)
	  fprintf(ps_file, "30 %d  moveto\n", start_pt -1 );
	else
	  fprintf(ps_file, "34 %d  moveto\n", start_pt -1 );

	fprintf(ps_file, "(%d) 8.0 Print\n", y_translate[a] );
	fprintf(ps_file, "45 %7.2f 50 %7.2f sequence_line\n",\
		start_pt + OFFSET,start_pt + OFFSET); 
			/* SN 2008:01 --- correct position for chain >= 0 */

	/**********************************************************************/
	/* Display insertion code over each inserted residue in the protein.  */
	/**********************************************************************/
	cumulative_insertions = 0;
	for( d = 0; d < total_insertions; d++ ){
	  if( insertion_chain_ID[d] == chain_IDs[a] ){
            
	    for( e = 0; e < number_of_insertions[d]; e++ ){
	      /* SN --- The first residue in a chain may be an insertion! 
			and the starting residue of any chain can be > 1!
			Hence, we use the same logic for all cases and
			disable this condition check.

			if( insertion_res_num[d] != 1 ){ 
	       */
		fprintf( ps_file, "55 %5.2f  moveto\n",	start_pt + 2*OFFSET + \
			( 3*(insertion_res_num[d] + cumulative_insertions + \
						e - y_translate[a] )));
		fprintf(ps_file, "(%c) 5.0 MiniUncenterRot90 Rot90\n",(e+'A') );
		fprintf(ps_file, "(%c) 5.0 Print\n",(e+'A') );
		fprintf(ps_file, "grestore\n");
		/* fprintf(ps_file, "(*) 10.0 Print\n" ); /
	      }
	      else{
		fprintf( ps_file, "55 %5.1f  moveto\n",  start_pt + 1.5 + (3 * e) -4  );
		fprintf(ps_file, "(%c) 5.0 UncenterRot90 Rot90\n",(e+'A') );
		fprintf(ps_file, "(%c) 5.0 Print\n",(e+'A') );
		fprintf(ps_file, "grestore\n");
		/ fprintf(ps_file, "(*) 10.0 Print\n" ); /
	      } */
	    }
	    cumulative_insertions += number_of_insertions[d];	
	  }
	}
	/**********************************************************************/

	/* 2008:01 SN	MISSING residues are marked with an X */
        cumulative_insertions = 0;
	if(missing_rsd_count[a] > 0)
	{
          d = 0;
	  for(e = 0; e < a; e++)
	  {
	    d += missing_rsd_count[e];	/* Find the index of the first missing
					   residue in the current chain */
	  }

	  b = 0; 			/* insertion rsd list index  */
	  for(e = d; e < d+missing_rsd_count[a]; e++ ) 
					/* Scan the missing rsd list */
	  {
            for( ; b < total_insertions; b++ )
	    {
	      if( insertion_res_num[b] < missing_rsd_list[e] && \
		  insertion_chain_ID[b] == chain_IDs[a])
	      {
		cumulative_insertions += number_of_insertions[b];
	      }
	      if( insertion_res_num[b] > missing_rsd_list[e] )
		break;
	    }	  

	    fprintf( ps_file, "55 %5.2f  moveto\n",start_pt - 1.0 + \
			( 3*(missing_rsd_list[e] + cumulative_insertions \
						- y_translate[a] )));
	    fprintf(ps_file, "(%c) 5.0 MiniUncenterRot90 Rot90\n",('X') );
	    fprintf(ps_file, "(%c) 5.0 Print\n",('X') );
	    fprintf(ps_file, "grestore\n");	
	  }
	}

	/************************************************************/
	/* Because the first residue in the protein need not be 1 or*/
	/* a multiple of 10, find the next residue number to print  */
	/* to the top of the page such that it is at least 10       */
	/* greater than the first residue. ie. if the first residue */
	/* in the protein is 34, the next res num to print would be */
	/* 50.                                                      */
	next_residue = (y_translate[a] + 10) / 10;
	next_residue *= 10;

	/* SN 2007:12:21 */
	/* 	This is where the output residue number is fixed.   */
	/*	So, account for the insertions here!		    */
	next_rsd_index = next_residue;
	for( d = 0; d < total_insertions; d++ )
        {
	  if( insertion_chain_ID[d] == chain_IDs[a] &&
	      insertion_res_num[d] >= y_translate[a] && 
	      insertion_res_num[d] < next_residue )
          {
            next_rsd_index = next_rsd_index + number_of_insertions[d];
	  } 
	}
      
#ifdef DEBUG_SN 
	printf("\nDEBUG (postscript.c): Next residue = %d",next_residue);
#endif

	/*	
	 for( b = next_residue; b < (chain_size[a] + y_translate[a]); b+=10 ) 
	 for( b = next_residue; b < (chain_size[a] + y_translate[a]); )
	*/
	for( b = next_residue; b < (min_max_rsd[a][1] ); )
	{
	  if( print_number )
	  {
	    if( b < 100 )
	    {
	      fprintf(ps_file, "34 %d  moveto\n", \
		 start_pt -1 + (3*(next_rsd_index - y_translate[a])));

	      fprintf(ps_file, "(%d) 8.0 Print\n", b );

	      fprintf(ps_file, "45 %5.2f 50 %5.2f sequence_line\n",\
			start_pt + OFFSET +(3*(next_rsd_index - y_translate[a])),
		      	start_pt + OFFSET +(3*(next_rsd_index - y_translate[a])));
	    }
	    else{
#ifdef DEBUG_SN
		printf("\n%d\t%d",start_pt -1 +(3*(b-y_translate[a])),\
		start_pt -1 +(3*(next_rsd_index - y_translate[a])));
#endif
		if(b > 1000)
		  fprintf(ps_file, "26 %d  moveto\n", \
		     start_pt -1 + (3*(next_rsd_index - y_translate[a])));
		else
	      	  fprintf(ps_file, "30 %d  moveto\n", \
		     start_pt -1 + (3*(next_rsd_index - y_translate[a])));

	      fprintf(ps_file, "(%d) 8.0 Print\n", b );
	      fprintf(ps_file, "45 %5.2f 50 %5.2f sequence_line\n",\
			start_pt + OFFSET +(3*(next_rsd_index - y_translate[a])),
		      	start_pt + OFFSET +(3*(next_rsd_index - y_translate[a])));
	    }
	  }

	  /* SN 2007:12:21 */
	  /* 	This is where the output residue number is fixed.   */
	  /*	So, account for the insertions here!		    */
#ifdef DEBUG_SN        
	 printf("\nDEBUG (postscript.c): b = %d, index = %d",b,next_rsd_index);
#endif

	 b+=10;
	 next_rsd_index += 10;
	 for( d = 0; d < total_insertions; d++ )
         {
	   if( insertion_chain_ID[d] == chain_IDs[a] &&
	       insertion_res_num[d] >= b-10 && 
	       insertion_res_num[d] < b )
           {
              /*b = b+10+number_of_insertions[d];*/
             next_rsd_index = next_rsd_index + number_of_insertions[d];
	   } 
	 }

	}
      }
      /************************************************************/

      start_pt = end_pt + 18;	/* Update the start point of next chain */
    }

    for( a = 0; a < cluster_counter; a++ ){

      nc_current = nc_count[a];

      while( nc_current ) {
	/*
	if( nc_current != NULL )
	  printf("how did i get here if nc_current is null? %d\n", nc_current->atom_number);
	*/
	if( colors[a] < 10 )
	  fprintf(ps_file,"Col0%d\n", colors[a] );
	else
	  fprintf(ps_file,"Col%d\n", colors[a] );

	/*printf("color: %3d\n", colors[a] );*/

	if( number_of_chains ){


	  for( b = 0; b <= number_of_chains; b++ ){

	    if( nc_current->chain_ID == chain_IDs[b] ) {

	      this_chain = b;
	      top    = ps_line_number + 3;
	      bottom = ps_line_number - 3;
	      for( c = 0; c < b; c++ )
		previous_chains_size += chain_size[c];
	      y_position = (nc_current->residue_number - y_translate[b] + (6*b) + previous_chains_size) * 3;
	      previous_chains_size = 0;
	    }

	    /**********************************************************************/
	    /* Modify spacing to include insertions that share the same residue   */
	    /* numbering.                                                         */
	    /**********************************************************************/
	    for( d = 0; d < total_insertions; d++ ){
	      if( insertion_chain_ID[d] == nc_current->chain_ID &&
		  insertion_chain_ID[d] == chain_IDs[b] ){

		if( nc_current->residue_number > insertion_res_num[d] ){
		  y_position += ((number_of_insertions[d]) * 3);
		}

		if( nc_current->residue_number == insertion_res_num[d] ){

		  if( nc_current->residue_number == 1 ){
		    if( nc_current->insertion_space )
		      y_position += ((nc_current->insertion_space - 1) * 3);
		    else
		      y_position += (number_of_insertions[d]) * 3;
		  }
		  else
		    y_position += (nc_current->insertion_space) * 3;
		}
	      }

	    }
	    /**********************************************************************/

	  }
	}

	else{
	  top    = ps_line_number + 3;
	  bottom = ps_line_number - 3;
	  y_position = ((nc_current->residue_number)*3) - (y_translate[0]*3);
	}

	/*printf("y_position %3d  res_num %3d,  %d %d %d  %d\n", y_position, nc_current->residue_number,
	  total_insertions, insertion_res_num[d], number_of_insertions[d], y_translate[this_chain] );*/

	if( !strcmp( nc_current->atom_type, "N" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+110,
		  top,    (y_position)+111,
		  bottom, (y_position)+111,
		  bottom, (y_position)+110 );
	}

	if( !strcmp( nc_current->atom_type, "CA" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+111,
		  top,    (y_position)+112,
		  bottom, (y_position)+112,
		  bottom, (y_position)+111 );
	}

	if( !strcmp( nc_current->atom_type, "C" ) ){
	  fprintf(ps_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    (y_position)+112,
		  top,    (y_position)+113,
		  bottom, (y_position)+113,
		  bottom, (y_position)+112 );
	}

	nc_current = nc_current->next_element;
      }
    }
  }

}

/**********************************************************************/
/* Prints the actual colorbars corresponding to the the rigid cluster */
/* decomposition. The information for each cluster in included as a   */
/* linked list. The lists for each cluster are indexed from the array */
/* nc_count[].                                                        */
/**********************************************************************/
void print_multipage_decomp( 	clusters *nc_count[MAX_CLUSTER_COUNT], 
				int cluster_counter,
			     	int ps_line_number, FILE *ps_file, 
				int colors[MAX_CLUSTER_COUNT],
			     	int y_translate[MAX_CHAIN_COUNT], 
				int number_of_residues,
			     	FILE *file_list[FILE_LIST_SIZE], 
				int number_of_pages,
			     	int number_of_chains, 
			     	int chain_size[MAX_CHAIN_COUNT],
			     	char chain_IDs[MAX_CHAIN_COUNT], 
				int total_insertions,
			     	int insertion_res_num[100], 
			     	int number_of_insertions[100],
			     	char insertion_chain_ID[100],
		   	     /* SN 2008:01 */
			     int min_max_rsd[MAX_CHAIN_COUNT][2], 
		   	     int missing_rsd_list[MISSING_RSD_LIST_SIZE],
		   	     int missing_rsd_count[MAX_CHAIN_COUNT]) {

  int
    a = 0,
    b = 0,
    c = 0,
    d = 0,
    e = 0,
    top = 0,
    bottom = 0,
    y_position = 0,
    file_number = 0,
    end_point = 0,
    start_pt=0,
    summed_length = 0, /* length of all the chains so far */
    end_pt = 0,
    next_residue = 0,
    last_res_num_on_this_page = 0,
    current_pg = 0,
    previous_length = 0,
    seq_page = 0,
    seq_num = 0,
    space = 0,
    times = 0,
    insertion_space = 0,
    residue_pointer = 0,
    printed_first_residue_number_yet = 0,
	/* 2008:01  SN */
    next_rsd_index,cumulative_insertions=0,prev_insertions=0,
    starting_rsd,max_rsd_index,summed_insertions=0,current_rsd_position;

  clusters
    *nc_current = NULL;

  FILE
    *current_file;

  /*****************************************************/
  /* print the sequence line (black line) on all pages */
  /* and a scale at the top of the page indicating the */
  /* residue number.                                   */
  /*****************************************************/
  start_pt = 110;
  current_pg = 0;

  for( a = 0; a <= number_of_chains; a++ ){

    printed_first_residue_number_yet = FALSE;
    seq_num  = 0;
    seq_page = 0;
    current_file = file_list[current_pg];
    end_pt = ( chain_size[a] * 3 ) + start_pt;

    /**********************************************************************/
    /* If end_pt > 710, then the sequence line extends off the page. Make */
    /* accomodations to print up to y=710, then put the rest on the next  */
    /* page.                                                              */
    /**********************************************************************/
    if( end_pt > 710 ){

      /**********************************************************************/
      /* Print what you can on the current page. Keep printing the current  */
      /* chain until the end_pt will fit on the page, then code to the next */
      /* snipet of code.                                                    */
      /**********************************************************************/
      while( end_pt > 710 ){

	fprintf(current_file, "\nCol00\n");
	fprintf(current_file, "%3d %3d %3d 710 sequence_line\n", ps_line_number,
		start_pt, ps_line_number );
	fprintf(current_file, "grestore\n");

	/**********************************************************************/
	/* print the residue numbers on the top of the page.                  */
	/**********************************************************************/
	if( !printed_first_residue_number_yet ){
	  
	  if(y_translate[a] >= 1000 && y_translate[a] < 10000)
	    fprintf(current_file, "26 %d  moveto\n", start_pt -1 );
	  else
	  if(y_translate[a] >= 100 && y_translate[a] < 1000)
	    fprintf(current_file, "30 %d  moveto\n", start_pt -1 );
	  else
	    fprintf(current_file, "26 %d  moveto\n", start_pt -1 );
	  fprintf(current_file, "(%d) 8.0 Print\n", y_translate[a] );
	  fprintf(current_file, "45 %5.2f 50 %5.2f sequence_line\n",\
		  start_pt +OFFSET, start_pt +OFFSET);
	  
	  /************************************************************/
	  /* Because the first residue in the protein need not be 1 or*/
	  /* a multiple of 10, find the next residue number to print  */
	  /* to the top of the page such that it is at least 10       */
	  /* greater than the first residue. ie. if the first residue */
	  /* in the protein is 34, the next res num to print wound be */
	  /* 50.                                                      */
	  /************************************************************/
	  next_residue = (y_translate[a] + 10) / 10;
	  next_residue *= 10;
	  printed_first_residue_number_yet = TRUE;
	  
	  /************************************************************/
	  /*                       SN 2008:02                         */
	  /* Set the upper bound for residue numbering on this page.  */
	  /*							      */
	  /* Note: chain_size is the ACTUAL size taking into account  */
	  /*	     both insertions and missing residues	      */
	  /* 	      						      */
	  /* Also, last_res_num is just a count and not the actual    */
	  /* rsd_index that would be displayed!                       */
	  /************************************************************/
	  last_res_num_on_this_page = chain_size[a] - ((end_pt - 710)/3) \
	    + y_translate[a];
	  starting_rsd = y_translate[a]; /* Required for counting insertions */
	  
	  /* SN 2007:12:21 */
	  /*      This is where the output residue number is fixed.   */
	  /*      So, account for the insertions here and find what   */
	  /*      is the largest rsd that is displayed on this page   */
	  
	  next_rsd_index = next_residue;
	  cumulative_insertions = 0;
	  max_rsd_index = last_res_num_on_this_page;
	  
	  for( d = 0; d < total_insertions; d++ )
            {
              if( insertion_chain_ID[d] == chain_IDs[a] && \
		  insertion_res_num[d] >= starting_rsd )
		{
		  /* The next if-statement is to take care of 
		     the first rsd to be printed */
		  if( insertion_res_num[d] < next_residue )
		    {
		      next_rsd_index = next_rsd_index + number_of_insertions[d];
		    }
		  
		  /* Here, we track how many insertions there are before the end
		     of the current page is reached */
		  if( insertion_res_num[d] + number_of_insertions[d] \
		      <= max_rsd_index )
		    {
		      max_rsd_index -= number_of_insertions[d];
		      cumulative_insertions += number_of_insertions[d];
		    }
		  else
		    {
		      /* If all the insertions of a residue cannot be accommodated
			 then truncate this page's output to one residue less    */
		      if ( insertion_res_num[d] <= max_rsd_index )
			{
			  max_rsd_index = insertion_res_num[d] - 1;
			  last_res_num_on_this_page = max_rsd_index + \
			    cumulative_insertions;
			}
		      break;
		    }
		}
	    }
	      
	  /* Revised residue index boundary */
	  for( b = next_residue; b < max_rsd_index; )
	    {
	      if(b >= 0 && b < 100)
		fprintf(current_file, "34 %d  moveto\n", \
			start_pt -1+(3*(next_rsd_index - starting_rsd )) );
	      else
	      if(b >= 100 && b < 1000)
		fprintf(current_file, "30 %d  moveto\n", 
			start_pt -1+(3*(next_rsd_index - starting_rsd )) );
	      else
		fprintf(current_file, "26 %d  moveto\n", 
			start_pt -1+(3*(next_rsd_index - starting_rsd )) );

	      fprintf(current_file, "(%d) 8.0 Print\n", b);
	      fprintf(current_file, "45 %5.2f 50 %5.2f sequence_line\n",\
		      start_pt + OFFSET + (3*(next_rsd_index - starting_rsd)),
		      start_pt + OFFSET + (3*(next_rsd_index - starting_rsd )) );
	      
	      b+=10;
	      next_rsd_index += 10;
	      for( d = 0; d < total_insertions; d++ )
		{
		  if( 	insertion_chain_ID[d] == chain_IDs[a] &&
			insertion_res_num[d] >= b-10 &&
			insertion_res_num[d] < b )
		    {
		      /*b = b+10+number_of_insertions[d];*/
		      next_rsd_index = next_rsd_index + number_of_insertions[d];
		    }
		}
	    }
	      
	}
	else		/* After first Page! */	  
	  {
	    last_res_num_on_this_page = chain_size[a] - ((end_pt - 710)/3) \
	      + y_translate[a];
	    
	    /* SN 2007:12:21 */
	    /*      This is where the output residue number is fixed.   */
	    /*      So, account for the insertions here and find what   */
	    /*      is the largest rsd that is displayed on this page   */
	    
	    next_rsd_index = seq_num;
	    max_rsd_index = last_res_num_on_this_page - prev_insertions;
	    cumulative_insertions = 0;
	    
	    for( d = 0; d < total_insertions; d++ )
	      {
		if( insertion_chain_ID[d] == chain_IDs[a] && \
		    insertion_res_num[d] >= starting_rsd )
		  {
		    /* The next if-statement is to take care of 
		       the first rsd to be printed */
		    if( insertion_res_num[d] < seq_num )
		      {
			next_rsd_index = next_rsd_index + number_of_insertions[d];
		      }
		    
		    /* Here, we track how many insertions there are before the end
		       of the current page is reached */
		    if( insertion_res_num[d] + number_of_insertions[d] \
			<= max_rsd_index )
		      {
			max_rsd_index -= number_of_insertions[d];
			cumulative_insertions += number_of_insertions[d];
		      }
		    else
		      {
			/* If all the insertions of a residue cannot be accommodated
			   then truncate this page's output to one residue less    */
			if ( insertion_res_num[d] <= max_rsd_index )
			  {
			    max_rsd_index = insertion_res_num[d] - 1;
			    last_res_num_on_this_page = max_rsd_index + \
							cumulative_insertions;
			  }
			break;
		      }
		  }
	      }
	    
	    for( b = seq_num; b < max_rsd_index; )
	      {
		if( b >= 0 && b < 100 )
		  fprintf(current_file, "34 %d  moveto\n",\
			  start_pt -1 +(3*(next_rsd_index - starting_rsd)));
		else
		if( b >= 100 && b < 1000 )
		  fprintf(current_file, "30 %d  moveto\n",\
			  start_pt -1 +(3*(next_rsd_index - starting_rsd)));
		else
		  fprintf(current_file, "26 %d  moveto\n",\
			  start_pt -1 +(3*(next_rsd_index - starting_rsd)));

		fprintf(current_file, "(%d) 8.0 Print\n", b );
		fprintf(current_file, "45 %5.2f 50 %5.2f sequence_line\n",
			start_pt+OFFSET+(3*(next_rsd_index-starting_rsd)),
			start_pt+OFFSET+(3*(next_rsd_index-starting_rsd)));
		
		b+=10;
		next_rsd_index += 10;
		for( d = 0; d < total_insertions; d++ )
		  {
		    if(  insertion_chain_ID[d] == chain_IDs[a] &&
			 insertion_res_num[d] >= b-10 &&
			 insertion_res_num[d] < b )
		      {
			/*b = b+10+number_of_insertions[d];*/
			next_rsd_index = next_rsd_index + number_of_insertions[d];
		      }
		  }
	      }
	    
	  } /* End of if ( printed 1st rsd yet? ) */


	/**********************************************************************/
	/* Display insertion code over each inserted residue in the protein.  */
	/**********************************************************************/
	summed_insertions = 0;
	for( d = 0; d < total_insertions; d++ ){
	  
	  /* Only account for those inserted residues that are represented */
	  /* on the current page by comparing with the first and last rsds */
	  /* i.e., seq_num and last_res_num_on_this_page		     */
	  
	  if( insertion_chain_ID[d] == chain_IDs[a] && \
	      insertion_res_num[d] >= starting_rsd && \
	      insertion_res_num[d] < max_rsd_index ){
	    
	    for( e = 0; e < number_of_insertions[d]; e++ )
	      {
		fprintf( current_file, "55 %5.2f  moveto\n",start_pt+2*OFFSET \
			 + ( 3*(insertion_res_num[d] + summed_insertions + \
				e - starting_rsd )));
		fprintf(current_file, "(%c) 5.0 MiniUncenterRot90 Rot90\n",\
			(e+'A') );
		fprintf(current_file, "(%c) 5.0 Print\n",(e+'A') );
		fprintf(current_file, "grestore\n");
	      }
	    summed_insertions += number_of_insertions[d];	
	  }
	}
	/********************************************************************/

	/*************** Display Missing rsds *****************/
	/* 2008:01 SN	MISSING residues are marked with an X */
        summed_insertions = 0;
	if(missing_rsd_count[a] > 0)
	{
          d = 0;
	  for(e = 0; e < a; e++)
	  {
	    d += missing_rsd_count[e];	/* Find the index of the first missing
					   residue in the current chain */
	  }

	  c = 0; 			/* insertion rsd list index  */
	  for(e = d; e < d+missing_rsd_count[a]; e++ ) 
					/* Scan the missing rsd list */
	  {
            /* Check if 'e' is a residue that belongs to current page */
            if( missing_rsd_list[e] > starting_rsd && \
		missing_rsd_list[e] < max_rsd_index) 
            {
              for( ; c < total_insertions; c++ )
	      {
	        if( insertion_res_num[c] < missing_rsd_list[e] && \
	  	    insertion_chain_ID[c] == chain_IDs[a] && \
		    insertion_res_num[c] > starting_rsd )
	        {
	  	  summed_insertions += number_of_insertions[c];
	        }
	        if( insertion_res_num[c] > missing_rsd_list[e] )
		  break;
	      }	  

	      fprintf( current_file, "55 %5.2f  moveto\n",start_pt - 1.0 + \
			( 3*(missing_rsd_list[e] + summed_insertions \
						- starting_rsd )));
	      fprintf(current_file, "(%c) 5.0 MiniUncenterRot90 Rot90\n",('X') );
	      fprintf(current_file, "(%c) 5.0 Print\n",('X') );
	      fprintf(current_file, "grestore\n");	
	    }
          }
	}
	/**************** End MissingRSD Display **************/
	/******************************************************/
	
	current_pg++;
	start_pt = 110;
/*	end_pt   = ( end_pt - ((last_res_num_on_this_page - starting_rsd)* 3)); */
	end_pt   = end_pt - 600;
	current_file = file_list[current_pg];
	seq_page++;
	seq_num = b;
	starting_rsd = max_rsd_index+1; /* Next page starting rsd */
	prev_insertions += cumulative_insertions; /* Account for the insertions
						     from the previous page */
	
      } /* End of while ( end_pt > 710 ) */

      /**********************************************************************/
      /* Print the rest on the next page. The variable "seq_num" will set a */
      /* marker for where to begin the sequence numbering.                  */
      /**********************************************************************/
      fprintf(current_file, "\nCol00\n");
      fprintf(current_file, "%3d %3d %3d %3d sequence_line\n", ps_line_number,
	      start_pt, ps_line_number, end_pt );
      fprintf(current_file, "grestore\n");

      /****************************************************************/ 
      /* Compute the insertions before the first rsd number output and
	 also compute what the maximum rsd number is on this last page */
      
      next_rsd_index = seq_num;
      max_rsd_index = chain_size[a] + y_translate[a] - prev_insertions;
      
      for( d = 0; d < total_insertions; d++ )
        {
          if( insertion_chain_ID[d] == chain_IDs[a] && \
	      insertion_res_num[d] >= starting_rsd )
	    {
	      /* The next if-statement is to take care of 
		 the first rsd to be printed */
	      if( insertion_res_num[d] < seq_num )
		{
		  next_rsd_index = next_rsd_index + number_of_insertions[d];
		}
	      
	      max_rsd_index -= number_of_insertions[d];
	    }
        }
      /****************************************************************/ 
      

      for( c = seq_num; c < max_rsd_index; )
	{
	  if( c < 100 ){
	    fprintf(current_file, "34 %d  moveto\n", start_pt -1 +\
		    (3*(next_rsd_index - starting_rsd)) );

	    fprintf(current_file, "(%d) 8.0 Print\n", c );
	    fprintf(current_file, "45 %5.2f 50 %5.2f sequence_line\n", \
		    start_pt + OFFSET + (3*(next_rsd_index - starting_rsd)),
		    start_pt + OFFSET + (3*(next_rsd_index - starting_rsd)) );
	  }
	  else{
	    if(c >= 1000)
	      fprintf(current_file, "26 %d  moveto\n", start_pt -1 +\
		    (3*(next_rsd_index - starting_rsd)) );
	    else
	      fprintf(current_file, "30 %d  moveto\n", start_pt -1 +\
		    (3*(next_rsd_index - starting_rsd)) );

	    fprintf(current_file, "(%d) 8.0 Print\n", c );
	    fprintf(current_file, "45 %5.2f 50 %5.2f sequence_line\n", \
		    start_pt + OFFSET + (3*(next_rsd_index - starting_rsd)),
		    start_pt + OFFSET + (3*(next_rsd_index - starting_rsd)) );
	  }
	  
	  c+=10;
	  next_rsd_index += 10;
	  for( d = 0; d < total_insertions; d++ )
	    {
	      if(  insertion_chain_ID[d] == chain_IDs[a] &&
		   insertion_res_num[d] >= c-10 &&
		   insertion_res_num[d] < c )
		{
		  next_rsd_index = next_rsd_index + number_of_insertions[d];
		}
	    }
	}

	/**********************************************************************/
	/* Display insertion code over each inserted residue in the protein.  */
	/**********************************************************************/
	summed_insertions = 0;
	for( d = 0; d < total_insertions; d++ ){
	  
	  /* Only account for those inserted residues that are represented */
	  /* on the current page by comparing with the first and last rsds */
	  /* i.e., seq_num and last_res_num_on_this_page		     */
	  
	  if( insertion_chain_ID[d] == chain_IDs[a] && \
	      insertion_res_num[d] >= starting_rsd && \
	      insertion_res_num[d] < max_rsd_index )
	    {
	      for( e = 0; e < number_of_insertions[d]; e++ )
		{
		  fprintf( current_file, "55 %5.2f  moveto\n",start_pt+2*OFFSET \
			   + ( 3*(insertion_res_num[d] + summed_insertions + \
				  e - starting_rsd )));
		  fprintf(current_file, "(%c) 5.0 MiniUncenterRot90 Rot90\n",\
			  (e+'A') );
		  fprintf(current_file, "(%c) 5.0 Print\n",(e+'A') );
		  fprintf(current_file, "grestore\n");
		}
	      summed_insertions += number_of_insertions[d];	
	    }
	}
	/********************************************************************/

	/********** Display Missing rsds **********************/
	/* 2008:01 SN	MISSING residues are marked with an X */
        summed_insertions = 0;
	if(missing_rsd_count[a] > 0)
	{
          d = 0;
	  for(e = 0; e < a; e++)
	  {
	    d += missing_rsd_count[e];	/* Find the index of the first missing
					   residue in the current chain */
	  }

	  b = 0; 			/* insertion rsd list index  */
	  for(e = d; e < d+missing_rsd_count[a]; e++ ) 
					/* Scan the missing rsd list */
	  {
            /* Check if 'e' is a residue that belongs to current page */
            if( missing_rsd_list[e] > starting_rsd && \
		missing_rsd_list[e] < max_rsd_index) 

            {
              for( ; b < total_insertions; b++ )
	      {
	        if( insertion_res_num[b] < missing_rsd_list[e] && \
	  	    insertion_chain_ID[b] == chain_IDs[a] && \
		    insertion_res_num[b] > starting_rsd)
	        {
	  	  summed_insertions += number_of_insertions[b];
	        }
	        if( insertion_res_num[b] > missing_rsd_list[e] )
		  break;
	      }	  

	      fprintf( current_file, "55 %5.2f  moveto\n",start_pt - 1.0 + \
			( 3*(missing_rsd_list[e] + summed_insertions \
						- starting_rsd )));
	      fprintf(current_file, "(%c) 5.0 MiniUncenterRot90 Rot90\n",('X') );
	      fprintf(current_file, "(%c) 5.0 Print\n",('X') );
	      fprintf(current_file, "grestore\n");	
	    }
          }
	}
	/**************** End MissingRSD Display **************/

    }
    /**********************************************************************/
    /* Use this code if the sequence line will fit on the current page.   */
    /* (That is, end_pt <= 710).                                          */
    /*									  */
    /* This is necessary in the case of multiple chains wherein smaller   */
    /* chains can be totally accommodated in one page!			  */
    /**********************************************************************/
    else
      {
	fprintf(current_file, "\nCol00\n");
	fprintf(current_file, "%3d %3d %3d %3d sequence_line\n", ps_line_number,
		start_pt, ps_line_number, end_pt );
	fprintf(current_file, "grestore\n");
	
	
	if(y_translate[a] >= 0 && y_translate[a] < 100)
	  fprintf(current_file, "34 %d  moveto\n", start_pt -1 );
	else
	if(y_translate[a] >= 100 && y_translate[a] < 1000)
	  fprintf(current_file, "30 %d  moveto\n", start_pt -1 );
	else
	  fprintf(current_file, "26 %d  moveto\n", start_pt -1 );

	fprintf(current_file, "(%d) 8.0 Print\n", y_translate[a] );
	fprintf(current_file, "45 %5.2f 50 %5.2f sequence_line\n",\
		start_pt +OFFSET, start_pt +OFFSET);
	
	/************************************************************/
	/* Because the first residue in the protein need not be 1 or*/
	/* a multiple of 10, find the next residue number to print  */
	/* to the top of the page such that it is at least 10       */
	/* greater than the first residue. ie. if the first residue */
	/* in the protein is 34, the next res num to print wound be */
	/* 50.                                                      */
	/************************************************************/
	next_residue = (y_translate[a] + 10) / 10;
	next_residue *= 10;
	printed_first_residue_number_yet = TRUE;
	
	starting_rsd = y_translate[a]; /* Required for counting insertions */
	
	/* SN 2007:12:21 */
	/*      This is where the output residue number is fixed.   */
	/*      So, account for the insertions here and find what   */
	/*      is the largest rsd that is displayed on this page   */
	
	next_rsd_index = next_residue;
	cumulative_insertions = 0;
		
	for( d = 0; d < total_insertions; d++ )
	  {
	    if( insertion_chain_ID[d] == chain_IDs[a] && \
		insertion_res_num[d] >= starting_rsd && \
		insertion_res_num[d] < next_residue)
	      {
		next_rsd_index = next_rsd_index + number_of_insertions[d];
	      }
	  }
	
	/* Revised residue index boundary */
	for( b = next_residue; next_rsd_index < (chain_size[a]+y_translate[a]); )
	  {
	    if(b >= 0 && b < 100)
	      fprintf(current_file, "34 %d  moveto\n", \
		      start_pt -1+(3*(next_rsd_index - starting_rsd )) );
	    else
	    if(b >= 100 && b < 1000)
	      fprintf(current_file, "30 %d  moveto\n", \
		      start_pt -1+(3*(next_rsd_index - starting_rsd )) );
	    else
	      fprintf(current_file, "26 %d  moveto\n", 
		      start_pt -1+(3*(next_rsd_index - starting_rsd )) );

	    fprintf(current_file, "(%d) 8.0 Print\n", b);
	    fprintf(current_file, "45 %5.2f 50 %5.2f sequence_line\n",\
		    start_pt + OFFSET + (3*(next_rsd_index - starting_rsd)),
		    start_pt + OFFSET + (3*(next_rsd_index - starting_rsd )) );
	    
	    b+=10;
	    next_rsd_index += 10;
	    for( d = 0; d < total_insertions; d++ )
	      {
		if( 	insertion_chain_ID[d] == chain_IDs[a] &&
			insertion_res_num[d] >= b-10 &&
			insertion_res_num[d] < b )
		  {
		    next_rsd_index = next_rsd_index + number_of_insertions[d];
		  }
	      }
	  }
	
	/**********************************************************************/
	/* Display insertion code over each inserted residue in the protein.  */
	/**********************************************************************/
	summed_insertions = 0;
	for( d = 0; d < total_insertions; d++ )
	  {
	    /* Only account for those inserted residues of this chain */
	    
	    if( insertion_chain_ID[d] == chain_IDs[a] )
	      {
		for( e = 0; e < number_of_insertions[d]; e++ )
		  {
		    fprintf( current_file, "55 %5.2f  moveto\n",start_pt+2*OFFSET \
			     + ( 3*(insertion_res_num[d] + summed_insertions + \
				    e - y_translate[a] )));
		    fprintf(current_file, "(%c) 5.0 MiniUncenterRot90 Rot90\n",\
			    (e+'A') );
		    fprintf(current_file, "(%c) 5.0 Print\n",(e+'A') );
		    fprintf(current_file, "grestore\n");
		  }
		summed_insertions += number_of_insertions[d];	
	      }
	  }
	/********************************************************************/

	/********** Display Missing rsds **********************/
	/* 2008:01 SN	MISSING residues are marked with an X */
        summed_insertions = 0;
	if(missing_rsd_count[a] > 0)
	{
          d = 0;
	  for(e = 0; e < a; e++)
	  {
	    d += missing_rsd_count[e];	/* Find the index of the first missing
					   residue in the current chain */
	  }

	  b = 0; 			/* insertion rsd list index  */
	  for(e = d; e < d+missing_rsd_count[a]; e++ ) 
					/* Scan the missing rsd list */
	  {
            for( ; b < total_insertions; b++ )
	    {
	      if( insertion_res_num[b] < missing_rsd_list[e] && \
		  insertion_chain_ID[b] == chain_IDs[a])
	      {
		summed_insertions += number_of_insertions[b];
	      }
	      if( insertion_res_num[b] > missing_rsd_list[e] )
		break;
	    }	  

	    fprintf( current_file, "55 %5.2f  moveto\n",start_pt - 1.0 + \
			( 3*(missing_rsd_list[e] + summed_insertions \
						- y_translate[a] )));
	    fprintf(current_file, "(%c) 5.0 MiniUncenterRot90 Rot90\n",('X') );
	    fprintf(current_file, "(%c) 5.0 Print\n",('X') );
	    fprintf(current_file, "grestore\n");	
	  }
	}
	/******************************************************/

      } /* End of "else if (end_pt < 710 )" */ 

    start_pt = end_pt + 18;
  }
  /**********************************************************************/
  /* End sequence line, numbering, tick marks, etc...                   */
  /**********************************************************************/


  /**********************************************************************/
  /* The following code prints the actual colorbar text to the          */
  /* postscript file.                                                   */
  /**********************************************************************/
  for( a = 0; a < cluster_counter; a++ ){

    nc_current = nc_count[a];

    while( nc_current ) {

      summed_length = previous_length = times = 0;
      y_position = insertion_space = 0;

      /**********************************************************************/
      /* If there are 2 or more chains, find the right output file.         */
      /**********************************************************************/
      if( number_of_chains ){

	for( b = 0; b <= number_of_chains; b++ ){

/*	  summed_length += (chain_size[b] + ( 6*b )); --- FAULTY! (SN 02:08) */
	  summed_length += (chain_size[b] + 6*(b>0?1:0) ); 
			/* SN 02:08	Add 6 to account for inter-chain 
			   spacing including the chain-ID (non-zero IDs only) */

	  if( nc_current->chain_ID == chain_IDs[b] ) {

	    /**********************************************************************/
	    /* Modify spacing to include insertions that share the same residue   */
	    /* numbering.                                                         */
	    /**********************************************************************/
	    for( d = 0; d < total_insertions; d++ ){
	      if( insertion_chain_ID[d] == nc_current->chain_ID &&
		  insertion_chain_ID[d] == chain_IDs[b] ){

		if( nc_current->residue_number > insertion_res_num[d] ){
		  y_position += ((number_of_insertions[d]) * 3);
		  insertion_space += number_of_insertions[d];
		}

		if( nc_current->residue_number == insertion_res_num[d] ){

		  if( nc_current->residue_number == 1 ){
		    if( nc_current->insertion_space ){
		      y_position += ((nc_current->insertion_space - 1) * 3);
		      insertion_space += nc_current->insertion_space -1;
		    }
		    else{
		      y_position += (number_of_insertions[d]) * 3;
		      insertion_space += number_of_insertions[d];
		    }
		  }
		  else{
		    y_position += (nc_current->insertion_space) * 3;
		    insertion_space += nc_current->insertion_space;
		  }
		}
	      }

	    }
	    /**********************************************************************/

	    /**********************************************************************/
	    /* Select which temp file to print to.                                */
	    file_number = ( summed_length -chain_size[b] + \
			    nc_current->residue_number -y_translate[b] + \
			    insertion_space) / 200;

	    //printf("file number %d\n", file_number );
	    current_file = file_list[file_number];


	    /* SN 02:08 */
	    current_rsd_position = ( summed_length -chain_size[b] + \
                            nc_current->residue_number -y_translate[b] + \
                            insertion_space) % 200;
	    y_position = current_rsd_position*3;

	    /**********************************************************************/

	    /**********************************************************************/
	    /* Tell the postscript file what color to print the current atom.     */
	    if( colors[a] < 10 )
	      fprintf(current_file,"Col0%d\n", colors[a] );
	    else
	      fprintf(current_file,"Col%d\n", colors[a] );
	    /**********************************************************************/

	    /**********************************************************************/
	    /* Calculate where on the current page the rigid cluster information  */
	    /* should be printed. This will depend on how many chains have been   */
	    /* printed previously (summed_length) and how much of the current     */
	    /* chain was printed on the last page.                                */
/*   	    
	    previous_length = summed_length - chain_size[b];
	    while( previous_length > 200 ) {
	      previous_length -= 200;
	    }
*/
	    /* residue_pointer represents the residue index;
	     * 	residue_pointer -= 200 => 
	     * 		How far from the start of the current chain is the 
	     *		current residue located on the current page?	   */
/*
	    residue_pointer = (nc_current->residue_number) -y_translate[b];
	    while( residue_pointer >= 200 ) {
	      residue_pointer -= 200;
	    }
*/
	    /* 	(200 - previous_length) represents the length of the current */
	    /*	chain that was displayed on the previous page --- SN 02/08   */
/*
	    if( ( 200 - previous_length ) > ( chain_size[b] ) )
	      y_position += (residue_pointer +1 +previous_length) *3;

	    else {
	      if( residue_pointer >= ( 200-previous_length) )
		y_position += (residue_pointer -( 200 -previous_length ) ) *3;
	      / * will fit on page with the previous length * /
	      else
		y_position += (residue_pointer +previous_length) *3;
	    }
*/
	    /**********************************************************************/

	  }
	}
	top    = ps_line_number + 3;
	bottom = ps_line_number - 3;

      }
      else{
	file_number = ( nc_current->residue_number -y_translate[0] ) / 200;
	current_file = file_list[file_number];
        top    = ps_line_number + 3;
	bottom = ps_line_number - 3;
	y_position = ( nc_current->residue_number -y_translate[0] -(file_number*200) ) *3;

	if( colors[a] < 10 )
	  fprintf(current_file,"Col0%d\n", colors[a] );
	else
	  fprintf(current_file,"Col%d\n", colors[a] );
      }

      if( !strcmp( nc_current->atom_type, "N" ) ){
	fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		top,    y_position +110,
		top,    y_position +111,
		bottom, y_position +111,
		bottom, y_position +110 );
      }

	if( !strcmp( nc_current->atom_type, "CA" ) ){
	  fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    y_position +111,
		  top,    y_position +112,
		  bottom, y_position +112,
		  bottom, y_position +111 );
	}

	if( !strcmp( nc_current->atom_type, "C" ) ){
	  fprintf(current_file,"%4d %4d %4d %4d %4d %4d %4d %4d Pl4\n",
		  top,    y_position +112,
		  top,    y_position +113,
		  bottom, y_position +113,
		  bottom, y_position +112 );
	}

	nc_current = nc_current->next_element;
    }
  }
}

/**********************************************************************/
