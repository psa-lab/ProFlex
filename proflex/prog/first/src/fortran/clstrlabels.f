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
      subroutine clstrlabels(size,chain,nclst,bin,maxsize,clst,
     &                 newlabel,natom,mult,pointr,maxr,nc2,nc1)
c ------------------------------------------------------------------------------
c                                             LAST UPDATED:        July 15, 1997
c                                             PROGRAM WRITTEN:     July 15, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
  
c       This subroutine sorts out the clusters in descending order in size (# of
c atoms in cluster) and then defines a cluster label map so that the 1st cluster
c is of biggest size and the last cluster is the smallest size. It also counts
c the number of clusters within each size-range and sets up the range pointer.
c ------------------------------------------------------------------------------
c INPUT: 
c    1)  natom   --->  number of atoms in the macromolecule.
c    2)  nclst   --->  number of rigid clusters in the macromolecule. 
c    3)  maxr    --->  maximum number of ranges in size used. 
c    4)  various working arrays  
c ------------------------------------------------------------------------------
c OUTPUT: 
c    1) newlabel()     defines new cluster labeling map. 
c    2) nc1,nc2        defines two special cluster labels in the stack.
c ------------------------------------------------------------------------------
      integer   r,pointr(maxr),mult(maxr)
      integer   newlabel(natom),clst(natom)
      integer   size(nclst),chain(nclst),bin(maxsize)
      do isize=1,maxsize
      bin(isize) = -1
      enddo
         do nc=1,nclst
         chain(nc) = -1
         enddo
c ------------------------------------------------------- group by cluster sizes
      do nc=1,nclst
      isize = size(nc) 
      if( bin(isize) .gt. 0 ) chain(nc) = bin(isize)
      bin(isize) = nc
      enddo
         do r=1,maxr
         mult(r) = 0
         enddo
c ------------------------------------------ set up map for new cluster labeling
      kc = 0
      do isize=maxsize,3,-1
      nc = bin(isize)
         if( nc .gt. 0 ) then
c ------------------------------------------------------ determine range of size
            if( isize .lt. 6 ) then
            r = 1
            elseif( isize .lt. 11 ) then
            r = 2
            elseif( isize .lt. 16 ) then
            r = 3
            elseif( isize .lt. 21 ) then
            r = 4
            else
            r = 5
            endif
c -------------------------------------- relabel all clusters of degenerate size
   50    continue
         kc = kc + 1 
         newlabel(nc) = kc
         mult(r) = mult(r) + 1
         nc = chain(nc)
         if( nc .gt. 0 ) goto 50
         endif
      enddo
c --------------------------------- consider isolated bonds and sites separately
      r = 1
      nc2 = kc + 1
c ---------------------------------- relabel all clusters of degenerate size = 2
      nc = bin(2)
         if( nc .gt. 0 ) then
  100    continue
         kc = kc + 1 
         newlabel(nc) = kc
         mult(r) = mult(r) + 1
         nc = chain(nc)
         if( nc .gt. 0 ) goto 100
         endif
      nc1 = kc + 1
c ---------------------------------- relabel all clusters of degenerate size = 1
      nc = bin(1)
         if( nc .gt. 0 ) then
  150    continue
         kc = kc + 1 
         newlabel(nc) = kc
         mult(r) = mult(r) + 1
         nc = chain(nc)
         if( nc .gt. 0 ) goto 150
         endif
      pointr(maxr) = 0
         do r=maxr,2,-1
         pointr(r-1) = pointr(r) + mult(r) 
         enddo
      return
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine shellsrt(n,a,b,c)
c     subroutine shellsrt(n,a,b)
      integer n
      integer a(n),b(n),c(n)
c     sorts an array a(1:n) into descending numerical order by Shell's method (diminishing
c     increment sort. n is input a is replaced on output by its sorted rearrangement.
c     array b(1:n) is rearranged in the process
      
      integer i,j,inc
      integer v,w,x
      
      inc=1
 1    inc=3*inc+1
      if(inc.le.n) goto 1
 2    continue
      inc=inc/3
      do i=inc+1,n
         v=a(i)
         w=b(i)
         x=c(i)
         j=i
 3       if(a(j-inc).lt.v) then
c     3           if(a(j-inc).gt.v) then           
            a(j) = a(j-inc)
            b(j) = b(j-inc)
            c(j) = c(j-inc)
            j=j-inc
            if(j.le.inc) goto 4
            goto 3
         endif
 4       a(j)=v
         b(j)=w
         c(j)=x
      enddo
      
      if(inc.gt.1) goto 2
      
      return
      end
