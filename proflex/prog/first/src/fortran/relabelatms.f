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
      subroutine relabelatms(nhinge,clst,chain,oldlabel,newlabel,
     &   natom,bin,color,nclst,paint,ncolor,calow,cabig,cslow,csbig)
c------------------------------------------------------------------------------
c AJR 03.27.02 
c This version no longer uses CFB "pseudo"-atoms called "BOND" sites below, 
c that portion has been commented out, just much simpler this way.
c
c ------------------------------------------------------------------------------
c                                             LAST UPDATED:       March 27, 2002
c                                             PROGRAM WRITTEN:     July 15, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     This subroutine relabels the atom serial numbers and the BOND-site labels 
c in order to obtain the following labeling structure: 
c
c      consecutive sequence of numbers  ----> 
c 1 2 3                ...              natom ... natom+nhinge  ... natom+nbonds
c            physical atoms               |              BOND sites            |
c ---------|-----------|---------|--------|-----------|------|-----|-----|-----|
c   C1     |   C2      |   ...   |  Cn    | colorless |  C1  |  C2 | ... | Cn  | 
c  | | | |   | | | | |            |  | | |              | | | | | |       | | |  
c
c   Ci = the i-the color
c   Each range of colors is subdivided such that all sites that are acting as 
c   BULK atoms or as NON-hinge bonds within a given rigid cluster are grouped. 
c   The colorless BOND sites represent hinge joints and their associated labels 
c   have the range:   (natom+1 .LE. s .LE. natom+nhinge)
c   NOTE:  natom+nhinge = pointbrc(1)
c ------------------------------------------------------------------------------
c INPUT: 
c    1) nclst         ---> number of distinct rigid clusters.
c    2) natom         ---> number of atoms in macromolecule. 
c    3) nhinge        ---> number of hinge joints in macromolecule. 
c    4) bin(),chain() ---> groups of BULK sites within rigid clusters. 
c ------------------------------------------------------------------------------
c OUTPUT:  
c    1) paint(),ncolor
c    2) newlabel()
c    3) oldlabel()
c    4) pointarc() and pointbrc()
c    5) calow(),cabig(),cslow(),scbig()
c ------------------------------------------------------------------------------
c                              Variable List 
c bin(nc) =    input: used to bin all BULK sites within rigid cluster nc.
c          as output: used to define pointarc which gives the beginning new 
c          atom label for the BULK atoms in rigid cluster nc. 
c calow(icolor) = the smallest "new atom label" assigned to the icolor-th color.
c cabig(icolor) = the biggest "new atom label" assigned to the icolor-th color.
c cslow(icolor) = the smallest "new site label" assigned to the icolor-th color.
c csbig(icolor) = the biggest "new site label" assigned to the icolor-th color.
c          Note that all so called site labels > natom  which can run up to ns. 
c chain(nc) =  input: linked array to group all BULK sites of rigid cluster nc. 
c          as output: pointbrc(nc) defines a beginning label for a BOND site 
c          that is a NON-hinge within rigid cluster nc. 
c color(nc) = a color code for the nc-th rigid cluster. 
c ctag = a runing dummy index to tag which colors cannot be used. 
c icolor = a dummy coloring index.
c iorder = a dummy index used to order color groups.
c nc = generic cluster label from 1 to nclst.
c ncolor = number of distinct colors used in the FIRST analysis. 
c nogood(icolor) = flag for determining which colors have already been considered.
c paint(iorder) = defines which color to paint. 
      include   'set_parameter'
      integer   ctag,s,sf,snew,so
      integer   linkref(nb2max),point(nsmax)
      integer   chain(natom),clst(natom)
      integer   oldlabel(natom),newlabel(natom),bin(nclst)
      integer   calow(maxcolor),cabig(maxcolor),nogood(0:maxcolor)
      integer   cslow(0:maxcolor),csbig(0:maxcolor)
      integer   multref(nsmax),color(nclst),paint(maxcolor)
      common/topology/ ns,linkref,point,multref
c ------------------------------------------------------------ format statements
 6000 format(5x,'ERROR: Relabeling atoms by color failed  ',
     &          '--->  relabelatms.f')
      ctag = 0
         do icolor=0,maxcolor
         nogood(icolor) = 0
         enddo
c -------------------------------------- initialize order to paint macromolecule
      ncolor = 0
      ctag = ctag + 1
         do nc=1,nclst
         icolor = color(nc)
            if( nogood(icolor) .lt. ctag ) then
            nogood(icolor) = ctag
            ncolor = ncolor + 1
            paint(ncolor) = icolor
            endif
         enddo
c ------------------------- relabel atoms in groups of color while consecutively 
c                                      grouping atoms in the same rigid clusters
      sf = 0
      do iorder=1,ncolor
c --------------------------------------------------- paint each color in groups
      icolor = paint(iorder)
      calow(icolor) = sf + 1
c ----------------------------------------------- search over all rigid clusters
         do nc=1,nclst
c ---------------------------------------------- select only the icolor-th paint
            if( color(nc) .eq. icolor ) then    
            so = bin(nc)
c ------------------------------------------- define pointarc() where bin() will
c                                                   be passed back as pointarc()
            bin(nc) = sf
c ------------------------------------------ relabel all BULK sites in the nc-th
c                                             cluster having the icolor-th paint
   50       continue
c ----------------------------------------------------- increment new atom label
            sf = sf + 1
c ------------------------------------------------- define newlabel and oldlabel
            newlabel(so) = sf
            oldlabel(sf) = so
c ------------------------------------- consider next BULK site in rigid cluster
            so = chain(so)
            if( so .gt. 0 ) goto 50
            endif
         enddo
      cabig(icolor) = sf
      enddo
c ------------------------------------------------------------------ ERROR CHECK
                  if( sf .ne. natom ) then
                  write(6,*)
                  write(6,*)
                  write(6,6000)
                  stop
                  endif
c ----------------------------------------------------------- initialize chain()
c                                                   regard chain() as pointbrc()
      do nc=1,nclst
      chain(nc) = -1
      enddo
c ----------------------------------------------------------- initialize cslow()
         do icolor=1,maxcolor
         cslow(icolor) = -1
         enddo
c$$$c ---------------------------------------- relabel BOND-sites in groups of color 
c$$$      index = natom + nhinge
c$$$c ------------------------------------- initialize range of colorless BOND sites
c$$$      cslow(0) = natom + 1
c$$$      csbig(0) = index
c$$$c ----------------------- use new labeling scheme to generate the "proper order"
c$$$         do snew=1,natom
c$$$         s = oldlabel(snew) 
c$$$         nc = clst(s)
c$$$         icolor = color(nc)
c$$$         ip = point(s)
c$$$c -------------------------------------- record index on first cluster encounter
c$$$c                                           and index+1 on first color encounter
c$$$         if( chain(nc) .lt. 0 ) then
c$$$         chain(nc) = index
c$$$         if( cslow(icolor) .lt. 0 ) cslow(icolor) = index + 1
c$$$         endif
c$$$c --------------------------------- search over all BOND sites within cluster nc
c$$$            do jn=1,multref(s)
c$$$            so = linkref(ip + jn)
c$$$c ------------------------------------------------ prevent double counting bonds
c$$$               if( s .gt. so ) then
c$$$c -------------------------------------------------------- count NON-hinge bonds
c$$$               if( nc .eq. clst(so) ) index = index + 1
c$$$               endif
c$$$            enddo
c$$$         csbig(icolor) = index 
c$$$         enddo
      return
      end
