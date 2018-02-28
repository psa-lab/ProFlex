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
      subroutine decomp_L3dbb(natom,tag,mark,laman,nlaman)
c ------------------------------------------------------------------------------
c PROGRAM WRITTEN BY: Donald J. Jacobs                              May 20, 1998
c ------------------------------------------------------------------------------
c                               Description
c        This program uniquely defines each Laman subgraph by the lowest bulk
c  site label within the region. 
c ------------------------------------------------------------------------------
c INPUT:
c    1) natom
c    2) laman()
c ------------------------------------------------------------------------------
c OUTPUT:
c    1) nlaman
c    2) laman()
c ------------------------------------------------------------------------------
      integer laman(-1:natom),mark(-1:natom),s,tag
      do s=1,natom
      laman(s) = laman( laman(s) )
      enddo
      nlaman = 0
      tag = tag + 1
         do s=1,natom
         label = laman(s)
            if( label .ne. s ) then
               if( mark(label) .lt. tag ) then
               mark(label) = tag
               nlaman = nlaman + 1
               endif
            endif
         enddo
      return
      end
