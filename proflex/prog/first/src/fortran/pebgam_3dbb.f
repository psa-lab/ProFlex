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
      subroutine pebgam_3dbb(iflop)
c ------------------------------------------------------------------------------
c                                       PROGRAM MODIFIED:          Sept 28, 2000
c                                       PROGRAM WRITTEN:            Feb  4, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c      This subroutine takes a generic three dimensional bond-bending network
c and identifies all the rigidity information and transmits this information 
c through the pebble index. 
c ------------------------------------------------------------------------------
c INPUT: 
c    1) The network topology  --> described by linkcf, point, multcf  
c ------------------------------------------------------------------------------
c OUTPUT: 
c    1) iflop,pebble() 
c ------------------------------------------------------------------------------
c                              Variable List 
c block(s) = a marker array that locates all previously tagged sites that were 
c     searched when looking for a free pebble.
c btrack( ) = a BACKWARD TRACK that is followed when rearranging pebbles also 
c     used as a marker array like block().
c id = dimension of network. Note that here id=3 and cannot be changed!
c iflop = current number of floppy modes in the network
c linkcf() = central-force nearest neighbor table. The j-th nearest neighbor to
c      atom, so, is given by: linkcf(point(so) + j) 
c iskiptag = a large value of tag such that tag < skiptag always. When the value
c     of pebble() is -1, indicating a free pebble, block() will not prevent the
c     detection of a free pebble as all tags are less than skiptag.  
c stress(s) = label for BULK site s within a given Laman subgraph.
c maxatm = maximum number of atoms that can be handled. 
c multcf(s) = final number of nearest central-force neighbors to atom s.
c natom =  number of atoms in the macromolecule
c nb1 = 1 = # of bonds to be placed between two sites (acting as rigid bodies)
c nb5 = 5 = # of bonds to be placed between two sites (acting as rigid bodies)
c pebble(,) = defines the current directed graph of the network
c point() = used to give appropriate index in linkcf()
c s,so,sf = site labels
c shell(j) = stores the list of sites checked in a breath-first search
c tag = a reference mark to which site was checked during a pebble search, and 
c     it is a running dummy index incremented by one before each new search.
c tc = a dummy index used for torsion constraints that lock dihedral motions.
c tcmol = number of torsion constraints to be applied within the molecule,
c      usually due to a peptide or resonant bond.
c th_min = (tcmol + 1) = the minimum value for index th
c th_max = (tcmol + nh) = the maximum value for index th
c torsion(j,tc) = <so-sf>  for j=1,2 respectively. Defines a pair of sites where
c      a torsion constraint is to be placed.
      include   'set_parameter'
      integer   s,so,sf,tag,tc,tcmol,th_min,th_max
      integer   linkcf(nb2max),point(maxatm),torsion(2,nbmax)
      integer   multcf(maxatm)
      integer   pebble(6,0:maxatm),stress(-1:maxatm)
      integer   block(-1:maxatm),shell(0:maxatm),btrack(0:maxatm)
      common/search/   btrack
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/topology/ natom,linkcf,point,multcf
      common/subgraph/ stress

c ------------------------------------------------------------------- initialize 
      block(-1) = iskiptag
      block(0) = iskiptag
      iflop = id*natom
      tag = 0
c ------------------------------------------------ pebble() = 0 => non-existence
      pebble(1,0) = 0
      pebble(2,0) = 0
      pebble(3,0) = 0
      pebble(4,0) = 0
      pebble(5,0) = 0
      pebble(6,0) = 0
c ------------------------------------------------- pebble() = -1 => free pebble
      do s=1,natom
         pebble(1,s) = -1
         pebble(2,s) = -1
         pebble(3,s) = -1
         pebble(4,s) =  0
         pebble(5,s) =  0
         pebble(6,s) =  0
         multcf(s)   =  0 
         block(s)    =  0
      enddo
c ------------------------------------------------------------- build up network
      index_max = 0
      do so=1,natom-1
         index_min = index_max + 1
         index_max = point(so+1)
c ----------------------------------------------- place CF neighbors of atom, so
         do indexo=index_min,index_max
            sf = linkcf(indexo)
c ------------------------------------------------ prevent double counting bonds
            if( sf .gt. so ) then
c -------------------------------------------------------- place CF-bond <so,sf>  
               if(     multcf(so) .eq. 0 ) then                                
                  pebble(1,so) = sf                                               
                  pebble(2,so) = sf                                               
                  pebble(3,so) = sf                                               
                  pebble(4,so) = sf                                               
                  pebble(5,so) = sf                                               
                  if(  multcf(sf) .gt. 1 ) then                                
                     iflop = iflop - 3                                            
                  elseif( multcf(sf) .eq. 1 ) then                             
                     pebble(6,sf) = -1                                 
                     iflop = iflop - 2                                            
                  else                                                         
                     pebble(4,sf) = -1                                            
                     pebble(5,sf) = -1                                            
                     iflop = iflop - 1                                            
                  endif                                                        
               elseif( multcf(sf) .eq. 0 ) then                                
                  pebble(1,sf) = so                                               
                  pebble(2,sf) = so                                               
                  pebble(3,sf) = so                                               
                  pebble(4,sf) = so                                               
                  pebble(5,sf) = so                                               
                  if( multcf(so) .eq. 1 ) then                                 
                     pebble(6,so) = -1                              
                     iflop = iflop - 2                                            
                  else                                                         
                     iflop = iflop - 3                                            
                  endif                                                        
               else                                                            
                  if( stress_label(so).ne.stress_label(sf) ) then              
                     if( multcf(so) .eq. 1 ) then                              
                        pebble(6,so) = -1                            
                        iflop = iflop + 1             
                     endif                                                     
                     if( multcf(sf) .eq. 1 ) then                              
                        pebble(6,sf) = -1                           
                        iflop = iflop + 1                 
                     endif                                                     
                     call cover_bonds(so,sf,nb5,iflop)                            
                  endif                                                        
               endif                                                           
               multcf(so) = multcf(so) + 1
               multcf(sf) = multcf(sf) + 1
            endif
         enddo
      enddo
c --------------------------------------------- place dihedral angle constraints
      do tc=1,tcmol
         so = torsion(1,tc)
         sf = torsion(2,tc)
         if((multcf(so) .gt. 0) .AND.( multcf(sf) .gt.0)) then
            call cover_bonds(so,sf,nb1,iflop)
         endif
      enddo

      return
      end
