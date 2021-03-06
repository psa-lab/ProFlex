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
      subroutine out_decomp(fdecomp,oldlabel,natom,clst,nclst)
c ------------------------------------------------------------------------------
c                                             LAST UPDATED:        July 17, 1997
c                                             PROGRAM WRITTEN:     July 17, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
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
c    1) creates  fdecomp  file. 
c ------------------------------------------------------------------------------
      integer      s_old
      integer      oldlabel(natom),clst(nclst)
      character*80 fdecomp
c ------------------------------------------------------------ open fdecomp file
      open(2,file=fdecomp,status='new')
      rewind(2)
c --------------------------------------------------- write rigidity information
         do s_old=1,natom
         write(2,*) clst(s_old) 
         enddo
      close(2)
      return
      end
