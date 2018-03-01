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

/**********************************************************************/
void print_data_lines( float x_bottom, float x_top, float y_bottom, float y_top, 
		       float x_min, float x_max, float y_min, float y_max,
		       FILE *ps_graph, float frac_LRC[600], float mean_coord[600],
		       int total_hbonds ) {

  int
    counter = 0;

  float
    x1,x2,y1,y2,
    x_factor = 0.0,
    y_factor = 0.0;

  x_factor = ( x_bottom - x_top )    / ( x_max - x_min ); /* bottom first due to page layout */
  y_factor = ( y_top    - y_bottom ) / ( y_max - y_min );
  
/*  printf("xfactor %5.2f yfactor %5.2f  counter %d \n", x_factor, y_factor, total_hbonds  );*/
    
  for( counter = 1; counter < total_hbonds; counter++ ) {
    if(mean_coord[counter -1] != 0.00) {
      x1 = x_bottom - (frac_LRC[counter -1] * x_factor);
      x2 = x_bottom - (frac_LRC[counter] * x_factor);
      y1 = y_bottom + ( (mean_coord[counter -1] - y_min) * y_factor);
      y2 = y_bottom + ( (mean_coord[counter] - y_min) * y_factor);
      fprintf(ps_graph, "%% %5.2f %5.2f\n", frac_LRC[counter-1],  mean_coord[counter-1] );
      fprintf(ps_graph, "%5.2f %5.2f %5.2f %5.2f draw_line\n",x1,y1,x2,y2);
    }
  }
    /*    fprintf(ps_graph, "%5.2f %5.2f %5.2f %5.2f draw_line\n",
	    x_bottom - ( frac_LRC[counter-1]   * x_factor ),
	    y_bottom + ( (mean_coord[counter-1] - y_min) * y_factor ),
	    x_bottom - ( frac_LRC[counter]     * x_factor  ),
	    y_bottom + ( (mean_coord[counter] - y_min)   * y_factor ) );
	    }*/
  
}
/**********************************************************************/

/************************************************************/
/*            The postscript header info                    */
/* Includes standard stuff to make the file "legal" level 2 */
/* postscript, as well as definitions for the colors used   */
/* and the geometric primitives.                            */
/************************************************************/
void plot_data( FILE *ps_graph, float Energy[600], float Mean_Coord[600], 
		float Frac_LRC[600], float floppy_modes[600], int total_hbonds ) {

  int 
    a,
    counter = 0,
    random_flag =0;

  float
    temp_nrg,
    scale_var = 0.0,
    x_position = 0.0,
    y_position = 0.0,
    x_min = 0.0,
    x_max = 0.0,
    y_min = 0.0,
    y_max = 0.0,
    x_bottom = 0.0,
    y_bottom = 0.0,
    x_top = 0.0,
    y_top = 0.0;
  
  fprintf(ps_graph, "%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(ps_graph, "%%%%Title: HBD.ps\n");
  fprintf(ps_graph, "%%%%Creator: Brandon Hespenheide\n");
  fprintf(ps_graph, "%%%%CreationDate: 1999\n");
  fprintf(ps_graph, "%%%%BoundingBox:0 0 620 800\n");
  fprintf(ps_graph, "%%%%DocumentFonts: Times-Roman\n");
  fprintf(ps_graph, "%%%%EndComments\n");
  fprintf(ps_graph, "%%%%BeginProlog\n");
  fprintf(ps_graph, "/Col { sethsbcolor } bind def\n");
  fprintf(ps_graph, "%% These are the colors for the flexibility scale and the\n");
  fprintf(ps_graph, "%% lines that display the hydrogen bonds.\n\n");

  fprintf(ps_graph, "/Blue {0.0000  0.0000  1.0000 setrgbcolor } def\n");
  fprintf(ps_graph, "/Red  {1.0000  0.0000  0.0000 setrgbcolor } def\n");

  fprintf(ps_graph, "/Col00 {0.0000  0.0000  0.0000 sethsbcolor } def\n");
  fprintf(ps_graph, "/Col01 {1.0000  1.0000  1.0000 sethsbcolor } def\n");
  
  fprintf(ps_graph, "/Poly4 { moveto lineto lineto lineto fill } bind def\n");
  fprintf(ps_graph, "/Pl4 { 8 copy Poly4 moveto moveto moveto moveto closepath stroke } bind def\n");
  fprintf(ps_graph, "/Poly3 { moveto lineto lineto fill } bind def\n");
  fprintf(ps_graph, "/Pl3 { 6 copy Poly3 moveto moveto moveto closepath stroke } bind def\n");
  fprintf(ps_graph, "/Print { /Times-Roman findfont exch scalefont setfont show } bind def  \n");
  fprintf(ps_graph, "/draw_line { moveto lineto stroke } bind def\n\n");
  fprintf(ps_graph, "/CenterRot90 {\n");
  fprintf(ps_graph, "  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_graph, "  exch stringwidth pop -2 div exch 3 div exch rmov\n");
  fprintf(ps_graph, " } bind def\n");
  fprintf(ps_graph, "/UncenterRot90 {\n");
  fprintf(ps_graph, "  dup /Times-Roman findfont exch scalefont setfont\n");
  fprintf(ps_graph, "  exch stringwidth } bind def\n");
  fprintf(ps_graph, "/Rot90 { gsave currentpoint translate 90 rotate } bind def\n");
  fprintf(ps_graph, "/Rot135 { gsave currentpoint translate 135 rotate } bind def\n");
  fprintf(ps_graph, "%%%%EndProlog\n\n");
  fprintf(ps_graph, "%%%%Page:    1   1\n");

  
  /**********************************************************************/
  /* Setup so that quadrant 2 doesn't get printed i.e. the LRC vs. Ehb  */
  /* whenever the dilution is random (if Energy[1] > Energy[2]          */
  /*             AJR 10.31.01                                           */
  /**********************************************************************/
/*if(Energy[1] > Energy[2] || Energy[total_hbonds-1] < Energy[1]) random_flag = 1;*/
  temp_nrg = Energy[1];
  for (a=2; a<6; a++){
    if(Energy[a] < temp_nrg) random_flag = 1;
  }
  /**********************************************************************/
  /* The graphs page is diveded into 4 quadrants.       --------------  */
  /* with the following designation for the             |     ||     |  */
  /* quadrants (shown at right for letter size          |  3  ||  4  |  */
  /* paper.                                             |     ||     |  */
  /*                                                    |     ||     |  */
  /*                                                    --------------  */
  /*                                                    |     ||     |  */
  /*                                                    |  1  ||  2  |  */
  /*                                                    |     ||     |  */
  /*                                                    |     ||     |  */
  /*                                                    --------------  */
  /* The plots will be printed side-ways, for viewing                   */
  /* in landscape.                                                      */
  /**********************************************************************/
  
  /**********************************************************************/
  /* Quadrant 1: Fraction of atoms in the largest rigid cluster vs.     */
  /*             Mean coordination vs.                                  */
  /**********************************************************************/
  fprintf(ps_graph, "%%Data for plot in quadrant 1. f_LRC vs. <r>\n\n");
  fprintf(ps_graph, "40 60 290 60 draw_line\n");
  fprintf(ps_graph, "290 60 290 390 draw_line\n\n");

  fprintf(ps_graph, "162 20 moveto\n");
  fprintf(ps_graph, "(F) 14 UncenterRot90 Rot90\n");
  fprintf(ps_graph, "(F) 14 Print\n");
  fprintf(ps_graph, "grestore\n");
  fprintf(ps_graph, "168 27 moveto\n");
  fprintf(ps_graph, "(LRC) 10 UncenterRot90 Rot90\n");
  fprintf(ps_graph, "(LRC) 10 Print\n");
  fprintf(ps_graph, "grestore\n\n");

  fprintf(ps_graph, "316 220 moveto\n");
  fprintf(ps_graph, "(<r>) 14 UncenterRot90 Rot90\n");
  fprintf(ps_graph, "(<r>) 14 Print\n");
  fprintf(ps_graph, "grestore\n\n");
  
  /**********************************************************************/
  /* Print the x-axis scale and tick marks.                             */
  /**********************************************************************/
  x_position = 312;
  y_position = 53;
  
  for( counter = 230; counter < 258; counter+=4  ){
    scale_var = counter / 100.00; 
  
    fprintf(ps_graph, "%5.2f %5.2f moveto\n", x_position, y_position );
    fprintf(ps_graph, "(%3.2f) 12 UncenterRot90 Rot135\n", scale_var );
    fprintf(ps_graph, "(%3.2f) 12 Print\n", scale_var );
    fprintf(ps_graph, "grestore\n");
    
    y_position += 50;
  }

  x_position = 290; 
  for( y_position = 60; y_position <= 390; y_position+=12.5 )
    fprintf(ps_graph, "%5.2f %5.2f %5.2f %5.2f draw_line\n", 
	    x_position, y_position, x_position-5, y_position );
  for( y_position = 60; y_position <= 390; y_position+=50 )
    fprintf(ps_graph, "%5.2f %5.2f %5.2f %5.2f draw_line\n", 
	    x_position, y_position, x_position-10, y_position );

  /**********************************************************************/
  /**********************************************************************/
  /* Print the y-axis scale and tick marks.                             */
  /**********************************************************************/
  y_position = 40;
  x_position = 294;
  for( counter = 0; counter <= 10; counter++ ){
    scale_var = counter / 10.0; 

    fprintf(ps_graph, "%5.2f %5.2f moveto\n", x_position, y_position );
    fprintf(ps_graph, "(%2.1f) 12 UncenterRot90 Rot90\n", scale_var );
    fprintf(ps_graph, "(%2.1f) 12 Print\n", scale_var );
    fprintf(ps_graph, "grestore\n");
    
    x_position -= 23;
  }

  /**********************************************************************/
  /* The "x" and "y" in the next set of variable declarations refer the */
  /* coordinate placement in the postscript file, NOT the x and y values*/
  /* being plotted (ie. they're NOT the <r> and frac_LRC).              */
  /**********************************************************************/
  x_bottom = 290.0;
  x_top    = 64.0;
  y_bottom = 60.0;
  y_top    = 390.0;
  x_min    = 0.0;
  x_max    = 1.0;
  y_min    = 2.30;
  y_max    = 2.58;

  print_data_lines( x_bottom, x_top, y_bottom, y_top, 
		    x_min, x_max, y_min, y_max, 
		    ps_graph, Frac_LRC, Mean_Coord, total_hbonds );

  /**********************************************************************/
  /* END QUADRANT 1 PLOTTING                                            */
  /**********************************************************************/
  if(!random_flag) {
  /**********************************************************************/
  /* Quadrant 2:                                                        */
  /**********************************************************************/
  fprintf(ps_graph, "%%Data for plot in quadrant 2. f_LRC vs. E\n\n");
  fprintf(ps_graph, "325 60 565 60 draw_line\n");
  fprintf(ps_graph, "565 60 565 390 draw_line\n");

  fprintf(ps_graph, "457 20 moveto\n");
  fprintf(ps_graph, "(F) 14 UncenterRot90 Rot90\n");
  fprintf(ps_graph, "(F) 14 Print\n");
  fprintf(ps_graph, "grestore\n");
  fprintf(ps_graph, "464 27 moveto\n");
  fprintf(ps_graph, "(LRC) 10 UncenterRot90 Rot90\n");
  fprintf(ps_graph, "(LRC) 10 Print\n");
  fprintf(ps_graph, "grestore\n\n");

  fprintf(ps_graph, "595 198 moveto\n");
  fprintf(ps_graph, "(Energy (kcal/mol)) 12 UncenterRot90 Rot90\n");
  fprintf(ps_graph, "(Energy (kcal/mol)) 12 Print\n");
  fprintf(ps_graph, "grestore\n\n");
  
  /**********************************************************************/
  /* Print the x-axis scale and tick marks.                             */
  /**********************************************************************/
  x_position = 589;
  y_position = 48;
  
  for( counter = 0; counter <= 90; counter+=10  ){
    scale_var = counter / 10.00; 
  
    fprintf(ps_graph, "%5.2f %5.2f moveto\n", x_position, y_position );
    fprintf(ps_graph, "(-%3.2f) 12 UncenterRot90 Rot135\n", scale_var );
    fprintf(ps_graph, "(-%3.2f) 12 Print\n", scale_var );
    fprintf(ps_graph, "grestore\n");
    
    y_position += 33;
  }

  x_position = 565; 
  for( y_position = 60; y_position <= 390; y_position+=16.5 )
    fprintf(ps_graph, "%5.2f %5.2f %5.2f %5.2f draw_line\n", 
	    x_position, y_position, x_position-5, y_position );
  for( y_position = 60; y_position <= 390; y_position+=33 )
    fprintf(ps_graph, "%5.2f %5.2f %5.2f %5.2f draw_line\n", 
	    x_position, y_position, x_position-10, y_position );

  /**********************************************************************/
  /**********************************************************************/
  /* Print the y-axis scale and tick marks.                             */
  /**********************************************************************/
  y_position = 40;
  x_position = 569;
  for( counter = 0; counter <= 10; counter++ ){
    scale_var = counter / 10.0; 

    fprintf(ps_graph, "%5.2f %5.2f moveto\n", x_position, y_position );
    fprintf(ps_graph, "(%2.1f) 12 UncenterRot90 Rot90\n", scale_var );
    fprintf(ps_graph, "(%2.1f) 12 Print\n", scale_var );
    fprintf(ps_graph, "grestore\n");
    
    x_position -= 23;
  }

  x_bottom = 565.0;
  x_top    = 335.0;
  y_bottom = 60.0;
  y_top    = 390.0;
  x_min    = 0.0;
  x_max    = 1.0;
  y_min    = 0.0;
  y_max    = -10.00;

  print_data_lines( x_bottom, x_top, y_bottom, y_top, 
		    x_min, x_max, y_min, y_max, 
		    ps_graph, Frac_LRC, Energy, total_hbonds );
  }
  /**********************************************************************/
  /* END QUADRANT 2 PLOTTING                                            */
  /**********************************************************************/

  /**********************************************************************/
  /* Quadrant 3:                                                        */
  /**********************************************************************/
  fprintf(ps_graph, "30 435 290 435 draw_line\n");
  fprintf(ps_graph, "290 435 290 760 draw_line\n");

  fprintf(ps_graph, "177 388 moveto\n");
  fprintf(ps_graph, "(f) 16 UncenterRot90 Rot90\n");
  fprintf(ps_graph, "(f) 16 Print\n");
  fprintf(ps_graph, "grestore\n");

  fprintf(ps_graph, "316 590 moveto\n");
  fprintf(ps_graph, "(<r>) 14 UncenterRot90 Rot90\n");
  fprintf(ps_graph, "(<r>) 14 Print\n");
  fprintf(ps_graph, "grestore\n\n");
  
  /**********************************************************************/
  /* Print the x-axis scale and tick marks.                             */
  /**********************************************************************/
  x_position = 312;
  y_position = 428;
  
  for( counter = 230; counter < 258; counter+=4  ){
    scale_var = counter / 100.00; 
  
    fprintf(ps_graph, "%5.2f %5.2f moveto\n", x_position, y_position );
    fprintf(ps_graph, "(%3.2f) 12 UncenterRot90 Rot135\n", scale_var );
    fprintf(ps_graph, "(%3.2f) 12 Print\n", scale_var );
    fprintf(ps_graph, "grestore\n");
    
    y_position += 50;
  }

  x_position = 290; 
  for( y_position = 435; y_position <= 760; y_position+=12.5 )
    fprintf(ps_graph, "%5.2f %5.2f %5.2f %5.2f draw_line\n", 
	    x_position, y_position, x_position-5, y_position );
  for( y_position = 435; y_position <= 760; y_position+=50 )
    fprintf(ps_graph, "%5.2f %5.2f %5.2f %5.2f draw_line\n", 
	    x_position, y_position, x_position-10, y_position );

  /**********************************************************************/
  /**********************************************************************/
  /* Print the y-axis scale and tick marks.                             */
  /**********************************************************************/
  y_position = 410;
  x_position = 294;
  for( counter = 0; counter <= 8; counter++ ){
    scale_var = counter / 100.0; 

    fprintf(ps_graph, "%5.2f %5.2f moveto\n", x_position, y_position );
    fprintf(ps_graph, "(%3.2f) 12 UncenterRot90 Rot90\n", scale_var );
    fprintf(ps_graph, "(%3.2f) 12 Print\n", scale_var );
    fprintf(ps_graph, "grestore\n");
    
    x_position -= 25;
  }

  /**********************************************************************/
  /* The "x" and "y" in the next set of variable declarations refer the */
  /* coordinate placement in the postscript file, NOT the x and y values*/
  /* being plotted (ie. they're NOT the <r> and frac_LRC).              */
  /**********************************************************************/
  x_bottom = 290.0;
  x_top    = 90.0;
  y_bottom = 435.0;
  y_top    = 760.0;
  x_min    = 0.00;
  x_max    = 0.08;
  y_min    = 2.30;
  y_max    = 2.58;

  print_data_lines( x_bottom, x_top, y_bottom, y_top, 
		    x_min, x_max, y_min, y_max, 
		    ps_graph, floppy_modes, Mean_Coord, total_hbonds );

  
  fprintf(ps_graph, "showpage\n");
  fprintf(ps_graph, "%%%%EOF\n");
}
/************************************************************/
