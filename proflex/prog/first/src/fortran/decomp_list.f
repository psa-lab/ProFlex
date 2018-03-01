cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  MSU ProFlex, formerly called FIRST, is a software developed to predict and  c
c  analyze protein flexibility.                                                c
c  This source file is a part of MSU ProFlex.                                  c
c                                                                              c
c  Copyright (C) 1997 - 2008, Michigan State University.                       c
c                                                                              c
c  This program is free software; you can redistribute to academic users only, c
c  it and/or modify it under the terms of the GNU General Public License,      c
c  version 2, as published by the Free Software Foundation.                    c
c                                                                              c
c  This program is distributed in the hope that it will be useful,             c
c  but WITHOUT ANY WARRANTY; without even the implied warranty of              c 
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               c
c  GNU General Public License for more details.                                c
c                                                                              c
c  You should have received a copy of the GNU General Public License           c
c  along with this program; if not, write to the Free Software                 c
c  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA,  c
c  or see http://www.gnu.org/licenses/gpl.txt                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
      subroutine decomp_list(oldlabel,natom,clst,nclst)
c ------------------------------------------------------------------------------
c                                             LAST UPDATED:       April 21, 1999
c                                             PROGRAM WRITTEN:    April 21, 1999
c                                       program written by:  Brandon Hespenheide
c                                                                     Don Jacobs
c ------------------------------------------------------------------------------
c                               Description
c     This subroutine generates a file with all the rigidity information 
c compressed into a single array of numbers corresponding to the 
c <xxxx_FIRST.chem>  PDB file.
c ------------------------------------------------------------------------------
c INPUT: 
c    1) clst and oldlabel
c ------------------------------------------------------------------------------
c OUTPUT:  
c    1) creates a compressed version of the decomp* file for use with the
c       program HBdilute.c. The bulk atom label for each atom is delimited
c       by commas. Maximum line length is 100 characters. 
c ------------------------------------------------------------------------------
      integer      s_old,s_new
      integer      oldlabel(natom),clst(nclst)
      character*80 decomp_file,bondwt_file
      common/hbanalys/ HBanalysis,decomp_file,bondwt_file

c ------------------------------------------------------------ open fdecomp file
 1000 format(i4,':',i4,':',i4,':',i4,':',i4,':',i4,':',i4,':',i4,':',
     &       i4,':',i4,':',i4,':',i4,':',i4,':',i4,':',i4,':',i4,':',
     &       i4,':',i4,':',i4,':',i4)
 2000 format('END')

      open(8,file=decomp_file,ACCESS='APPEND')

c --------------------------------------------------- write rigidity information
      do s_old=1,natom,20
         write(8,1000) clst(s_old),clst(s_old+1),clst(s_old+2),
     &        clst(s_old+3),clst(s_old+4),clst(s_old+5),
     &        clst(s_old+6),clst(s_old+7),clst(s_old+8),
     &        clst(s_old+9),clst(s_old+10),clst(s_old+11),
     &        clst(s_old+12),clst(s_old+13),clst(s_old+14),
     &        clst(s_old+15),clst(s_old+16),clst(s_old+17),
     &        clst(s_old+18),clst(s_old+19)
      enddo
      write(8,2000) 
      
      close( 8 )

      return
      end
