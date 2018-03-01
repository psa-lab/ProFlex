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

      subroutine decomp_R3dbb(n_one,n_two,nclst)
c ------------------------------------------------------------------------------
c PROGRAM WRITTEN BY: Donald J. Jacobs                              Feb  5, 1997
c                                                   Last Modified:  May 15, 1998
c ------------------------------------------------------------------------------
c                               Description
c       This program uses the current connectivity of the bond-bending network,
c the current state of the pebble index (characterizing the network rigidity) 
c and the prior rigid cluster decomposition (defining known rigid substructures)
c to obtain the current rigid cluster decomposition.
c ------------------------------------------------------------------------------
c INPUT:
c    1) natom                          => number of atoms
c    2) pebble()                       => defines the state of rigidity.
c    3) linkcf(),point(),multcf()      => network connectivity
c    4) clst()                         => prior rigid cluster decomposition
c ------------------------------------------------------------------------------
c OUTPUT:
c    1) nclst                          => number of distinct rigid clusters
c    1) clst()                         => new rigid cluster decomposition
c ------------------------------------------------------------------------------
      include   'set_parameter'
      integer   s,sf,sj,smin,snew,so,stest,tag,floppy
      integer   linkcf(nb2max),point(maxatm),clst(maxatm)
      integer   pebble(6,0:maxatm),shell(0:maxatm),block(-1:maxatm)
      integer   btrack(0:maxatm) 
      integer   multcf(maxatm)

      common/rigidity/ pebble,block,tag,shell,kmax
      common/topology/ natom,linkcf,point,multcf
      common/search/   btrack
      common/clusters/ clst
c ==============================================================================
      n_one = 0
      n_two = 0
      nclst = 0
      do so=1,natom
         if( multcf(so) .gt. 1 ) then
            if( clst(so) .eq. so ) then
            nclst = nclst + 1
            call group_max_dof(so)
            block(so) = irigid
            floppy = - tag
            shell(1) = so
            kmax = 1
            k = 0
  100          k = k + 1
               if( k .le. kmax ) then
               s = shell(k) 
               indexs = point(s)
                  do js=1,multcf(s)
                  indexs = indexs + 1
                  stest = linkcf(indexs) 
                     if( multcf(stest) .ne. 1 ) then
                        if( block(stest) .ne. irigid ) then
                        smin = clst(stest) 
                           if( block(smin) .eq. irigid ) then
                           kmax = kmax + 1
                           clst(stest) = so
                           shell(kmax) = stest
                           block(stest) = irigid
                           elseif( (      smin  .ge. so    ) 
     &                      .AND.  (block(smin) .ne. floppy) ) then
      tag = tag + 1
      jmax = kmax + 1
      block(stest) = tag
      btrack(stest) = -1
      shell(jmax) = stest
      j = kmax
  200    j = j + 1
         if( j .le. jmax ) then
         sj = shell(j)
            do jj=1,6
            snew = pebble(jj,sj)
               if( block(snew) .lt. tag ) then
               sf = clst(snew)
                  if(  (      sf  .lt. so    )
     &            .OR. (block(sf) .eq. floppy) ) then
  300                block( clst(sj) ) = floppy
                     sj = btrack(sj)
                     if( sj .gt. 0 ) goto 300
                  goto 1000                                           
                  else
                  jmax = jmax + 1
                  shell(jmax) = snew
                  block(snew) = tag
                  btrack(snew) = sj
                  endif
               elseif( snew .lt. 0 ) then
  400          block( clst(sj) ) = floppy
               sj = btrack(sj)
               if( sj .gt. 0 ) goto 400
               goto 1000                                             
               endif
            enddo
         goto 200
         else
            do j=kmax+1,jmax
            clst( shell(j) ) = so
            block( shell(j) ) = irigid
            enddo
         kmax = jmax
         endif
 1000 continue
                           endif
                        endif
                     endif
                  enddo
               goto 100
               endif
               do k=1,kmax
               block( shell(k) ) = 0
               enddo
            endif
         elseif( multcf(so) .eq. 1 ) then
         sf = linkcf(point(so) + 1)
            if( multcf(sf) .eq. 1 ) then
            if( sf .gt. so ) n_two = n_two + 1
            endif
         else
         n_one = n_one + 1
         endif
      enddo
      nclst = n_one + n_two + nclst
      return
      end
