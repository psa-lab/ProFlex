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
      subroutine colorclustr(clst,chain,natom,bin,color,nclst,
     &                       mult,pointr,cname)
c ------------------------------------------------------------------------------
c                                             LAST UPDATED:        July 15, 1997
c                                             PROGRAM WRITTEN:     July 15, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     This subroutine assigns a color to each rigid cluster such that no 
c neighboring cluster (shared by a hinge joint) has the same color. 
  
c ------------------------------------------------------------------------------
c INPUT: 
c    1) clst        ---> cluster labels 
c    2) natom       ---> number of atoms in the macromolecule. 
c    3) nclst       ---> number of distinct rigid clusters.
c    4) pointr()    ---> index for range of cluster size
c    5) mult()      ---> number of clusters within given range. 
c ------------------------------------------------------------------------------
c OUTPUT:  
c    1) color(nc)   ---> a color code defining the color of the nc-th cluster.
c ------------------------------------------------------------------------------
c                              Variable List 
 
c chain() = linked array to group all sites belonging to a given rigid cluster. 
c choice(icolor,ichoice) = color code with the icolor-th level of preference for
c     the ichoice-th arrangement of colors.  
c cname(i) = the name of the i-th color
c color(nc) = a color code for the nc-th rigid cluster. 
c ctag = a running dummy index to tag which colors cannot be used. 
c icolor = a dummy coloring index.
c maxch = maximum number of preference sets to choose from. 
c maxcolor = maximum number of different colors used. 
c maxr =  the maximum number of different ranges. 
c mult(r) = the number of rigid clusters within the r-th range of sizes. 
c nc = generic cluster label from 1 to nclst.
c nogood(icolor)= flag for determining which colors have already been considered
c pointr(r) = a pointer index to determine all clusters within range r.
c r = range index. Various cluster sizes are placed into a pre-specified range.
c bin() = used to bin all BULK sites within a given cluster.  

c ==============================================================================
      include      'set_parameter'
      integer      ctag,r,s,sf,so
      integer      linkref(nb2max),point(nsmax)
      integer      clst(natom),chain(natom),bin(nclst)
      integer      mult(maxr),pointr(maxr),nogood(0:maxcolor)
      integer      multref(nsmax),color(nclst) 
      integer      choice(maxcolor,maxch)
      character*13 cname(maxcolor)
      common/topology/ ns,linkref,point,multref
c ==============================================================================
c ------------------------------------------------------------ format statements
 6000 format(5x,'ERROR: Coloring scheme failed',
     &          '  --->  colorclustr.f')
      do s=1,natom
      chain(s) = -1
      enddo
         do nc=1,nclst
         bin(nc) = -1
         color(nc) = 0 
         enddo
c ----------------------------- group BULK sites belonging to same rigid cluster
      do s=1,natom
      nc = clst(s)
      if( bin(nc) .gt. 0 ) chain(s) = bin(nc) 
      bin(nc) = s
      enddo
c -------------------------------------------------- color scheme initialization
c                                                                        color type
      cname(1) = 'yellow       '
      cname(2) = 'cyan         '
      cname(3) = 'magenta      '
      cname(4) = 'green        '                                        
      cname(5) = 'red          '
      cname(6) = 'blue         '
      cname(7) = 'grey         '
      cname(8) = 'white        '
      cname(9) = 'black        '
      cname(10)= 'orange       ' 
c     cname(11)= '[180,90,180] '                                         pink
c              etc ....
c     maxcolor = 10                                      change in set_parameter
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c ------------------------------------------------------- initialize constraints
         do icolor=0,maxcolor
         nogood(icolor) = 0
         enddo
      ctag = 0
c ---------------------------------------------------- define preference choices
c     range -> 1
      choice(1,1) = 5                                                    red
      choice(2,1) = 1                                                    yellow
      choice(3,1) = 4                                                    green
      choice(4,1) = 10                                                   orange
      choice(5,1) = 8                                                    white
      choice(6,1) = 3                                                    magenta
      choice(7,1) = 7                                                    grey
      choice(8,1) = 2                                                    cyan
      choice(9,1) = 6                                                    blue
c ------------------------------------------------------------------------------
c     range -> 2
      choice(1,2) = 5                                                    red
      choice(2,2) = 1                                                    yellow
      choice(3,2) = 4                                                    green
      choice(4,2) = 10                                                   orange
      choice(5,2) = 3                                                    magenta
      choice(6,2) = 8                                                    white
      choice(7,2) = 7                                                    grey
      choice(8,2) = 2                                                    cyan
      choice(9,2) = 6                                                    blue
c ------------------------------------------------------------------------------
c     range -> 3
      choice(1,3) = 1                                                    yellow
      choice(2,3) = 4                                                    green
      choice(3,3) = 3                                                    magenta
      choice(4,3) = 10                                                   orange
      choice(5,3) = 7                                                    grey
      choice(6,3) = 8                                                    white
      choice(7,3) = 2                                                    cyan
      choice(8,3) = 6                                                    blue
      choice(9,3) = 5                                                    red
c ------------------------------------------------------------------------------
c     range -> 4
      choice(1,4) = 3                                                    magenta
      choice(2,4) = 7                                                    grey
      choice(3,4) = 4                                                    green
      choice(4,4) = 10                                                   orange
      choice(5,4) = 4                                                    cyan
      choice(6,4) = 5                                                    white
      choice(7,4) = 6                                                    yellow
      choice(8,4) = 9                                                    blue
      choice(9,4) = 8                                                    red
c ------------------------------------------------------------------------------
c     range -> 5       1st choice when r=5 
      choice(1,5) = 9                                                    black
      choice(2,5) = 6                                                    blue
      choice(3,5) = 2                                                    cyan
      choice(4,5) = 3                                                    magenta
      choice(5,5) = 7                                                    grey
      choice(6,5) = 10                                                   orange
      choice(7,5) = 4                                                    green
      choice(8,5) = 1                                                    yellow
      choice(9,5) = 5                                                    red
c ------------------------------------------------------------------------------
c     range -> 5       2nd choice when r=5 
      choice(1,6) = 6                                                    blue
      choice(2,6) = 2                                                    cyan
      choice(3,6) = 9                                                    black
      choice(4,6) = 3                                                    magenta
      choice(5,6) = 7                                                    grey
      choice(6,6) = 10                                                   orange
      choice(7,6) = 4                                                    green
      choice(8,6) = 1                                                    yellow
      choice(9,6) = 5                                                    red
c ------------------------------------------------------------------------------
c     range -> 5       3rd choice when r=5 
      choice(1,7) = 2                                                    cyan
      choice(2,7) = 6                                                    blue
      choice(3,7) = 9                                                    black
      choice(4,7) = 3                                                    magenta
      choice(5,7) = 7                                                    grey
      choice(6,7) = 10                                                   orange
      choice(7,7) = 4                                                    green
      choice(8,7) = 1                                                    yellow
      choice(9,7) = 5                                                    red
c ------------------------------------------------------------------------------

c ----------------------------------------------------------- color the clusters 
      r = 5
      ipr = pointr(r)
         do ko=1,mult(r)
         nco = ipr + ko
         so = bin(nco)
         ctag = ctag + 1
  100       continue
c ---------------------------------- check over all nearest neighbors to site so
            ipo = point(so)
               do jo=1,multref(so)
               sf = linkref(ipo + jo)
c ------------------------------------------- only real atoms need consideration
                  if( sf .le. natom ) then
                  ncf = clst(sf)
c ----------------------------------------------------------- tag color conflict
                     if( nco .ne. ncf ) then
                     nogood( color(ncf) ) = ctag
                     endif
                  endif
               enddo
c ------------------------------------------ consider BULK site in rigid cluster
            so = chain(so)
            if( so .gt. 0 ) goto 100
c ----------------------------------------------- determine choice of preference
         j = r + mod(ko,3) 
c -------------------------------- try a series of colors in order of preference
            do i=1,9
            icolor = choice(i,j) 
            if( nogood(icolor) .ne. ctag ) goto 150
            enddo
c ------------------------------------------- ERROR: NO GOOD COLOR CAN BE FOUND!
                              write(6,*)
                              write(6,*)
                              write(6,6000)
                              stop
c -------------------------------------------------- assign color to cluster nco
  150    color(nco) = icolor
         enddo
c ----------------------------------------------- continue with ranges r=4,3,2,1
      do r=4,1,-1
      ipr = pointr(r)
         do ko=1,mult(r)
         nco = ipr + ko
         so = bin(nco)
         ctag = ctag + 1
  200       continue
c ---------------------------------- check over all nearest neighbors to site so
            ipo = point(so)
               do jo=1,multref(so)
               sf = linkref(ipo + jo)
c ------------------------------------------- only real atoms need consideration
                  if( sf .le. natom ) then
                  ncf = clst(sf)
c ----------------------------------------------------------- tag color conflict
                     if( nco .ne. ncf ) then
                     nogood( color(ncf) ) = ctag
                     endif
                  endif
               enddo
c ------------------------------------- consider next BULK site in rigid cluster
            so = chain(so)
            if( so .gt. 0 ) goto 200
c -------------------------------- try a series of colors in order of preference
            do i=1,9
c ----------------------------------------------------- choice of preference = r
            icolor = choice(i,r) 
            if( nogood(icolor) .ne. ctag ) goto 250
            enddo
c -------------------------------------------- ERROR: NO GOOD COLOR CAN BE FOUND
                              write(6,*)
                              write(6,*)
                              write(6,6000)
                              stop
c -------------------------------------------------- assign color to cluster nco
  250    color(nco) = icolor
         enddo
      enddo
      return
      end

