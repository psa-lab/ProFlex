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
      subroutine decomp_H3dbb(nhinge)
c ------------------------------------------------------------------------------
c PROGRAM WRITTEN BY: Donald J. Jacobs                             June  3, 1998
c ------------------------------------------------------------------------------
c                               Description
c      This program takes the rigid cluster decomposition of a three dimensional
c  bond-bending network and counts and stores the location of all hinge joints.
c ------------------------------------------------------------------------------
c INPUT:
c    1) natom,linkcf(),point(),multcf() 
c    2) clst() 
c ------------------------------------------------------------------------------
c OUTPUT:
c    1) nhinge
c    2) th_min,th_max,torsion()   for (th= th_min  to  th= th_max)
c ------------------------------------------------------------------------------
      include   'set_parameter'
      integer   so,sf,tc,tcmol,th_max,th_min
      integer   linkcf(nb2max),point(maxatm)
      integer   multcf(maxatm)
      integer   clst(maxatm),torsion(2,nbmax)
      integer   link_hinge(nb2max),label_hinge(nbmax)                  
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/topology/ natom,linkcf,point,multcf
      common/clusters/ clst
      common/cmotions/ label_hinge,link_hinge
      indexf = point(natom) + multcf(natom)
         do index=1,indexf
         link_hinge(index) = -1
         enddo
      th_min = tcmol + 1
      tc = tcmol
         do so=1,natom-1
            if( multcf(so) .gt. 1 ) then
            iclst = clst(so)
            index = point(so)
               do j=1,multcf(so)
               index = index + 1
               sf = linkcf(index)
                  if( sf .gt. so ) then
                     if( iclst .ne. clst(sf) ) then
                        if( multcf(sf) .ne. 1 ) then
                        tc = tc + 1
                        torsion(1,tc) = so
                        torsion(2,tc) = sf
                        label_hinge(tc) = -1
                        link_hinge(index) = tc
                        index2 = point(sf)
                           do jf=1,multcf(sf)
                           index2 = index2 + 1
                              if( linkcf(index2) .eq. so ) then
                              link_hinge(index2) = tc
                              endif
                           enddo
                        endif
                     endif
                  endif
               enddo
            endif
         enddo 
      th_max = tc
      nhinge = th_max - th_min + 1
      return
      end

