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
 
      subroutine decomp_C3dbb(natom,nih,ncmode,iflop)
c ------------------------------------------------------------------------------
c                                       PROGRAM WRITTEN:            June 4, 1998
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c      This subroutine takes a generic three dimensional bond-bending network
c and identifies all collective motions. The pebble index is modified as
c additional torsion constraints are placed at each hinge joint in the molecule.
c ------------------------------------------------------------------------------
c INPUT: 
c    1) The network topology  --> described by linkcf, point, multcf  
c ------------------------------------------------------------------------------
c OUTPUT: 
c    1) nih
c    2) ncmode
c    3) label_hinge() =>  labels subgraphs defining regions of collective motion
c ------------------------------------------------------------------------------
      include 'set_parameter'
      integer so,sf,tag,tc,tcmol,th_min,th_max
      integer pebble(6,0:maxatm),torsion(2,nbmax)
      integer block(-1:maxatm),shell(0:maxatm),btrack(0:maxatm)
      integer   link_hinge(nb2max),label_hinge(nbmax)               
      common/search/   btrack
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/cmotions/ label_hinge,link_hinge
         do tc=th_min,th_max
            if( label_hinge(tc) .lt. 0 ) then
            so = torsion(1,tc)
            sf = torsion(2,tc)
            kflop = iflop
            call lock_hinge(so,sf,iflop)
            if( iflop .lt. kflop ) torsion(2,tc) = -sf
            endif
         enddo
         do tc=th_min,th_max
         label = label_hinge(tc)
         if( label .gt. 0 ) label_hinge(tc) = label_hinge(label)
         enddo
      nih = 0
      ncmode = 0
      tag = tag + 1
         do tc=th_min,th_max
         label = label_hinge(tc)
            if( label .lt. 0 ) then
            nih = nih + 1
            else
               if( block(label) .ne. tag ) then
               ncmode = ncmode + 1
               block(label) = tag
               endif
            endif
         enddo
      return
      end

