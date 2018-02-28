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

      subroutine mapmolecule()
c ------------------------------------------------------------------------------
c                                       PROGRAM WRITTEN:           July  8, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     Given the topology of a macromolecule, forming a three dimensional generic
c bond-bending network, this subroutine obtains the Floppy Inclusion and Rigid
c Substructure Topography (FIRST) of the molecule. 
c ------------------------------------------------------------------------------
c INPUT: 
c   1) The network topology  --> described by linkcf, point, multcf  
c ------------------------------------------------------------------------------
c OUTPUT: 
c   1) nflop     --> total number of independent degrees of freedom
c   2) iflop     --> # of floppy modes (internal independent degrees of freedom)
c   3) n_one     --> number of isolated atoms
c   4) n_two     --> number of isolated dimers
c   5) nbody     --> number of disconnected bodies (three or more atoms)
c   6) nclst     --> number of distinct rigid clusters 
c   7) nstress   --> number of distinct Laman subgraphs each defining a disjoint
c                    stressed region, not accounting for molecular chemistry.
c   8) nhinge    --> number of hinge joints within the macromolecule

c   9) nih       --> number of independent hinge joints (simple rotormers)
c  10) ncmode    --> number of distinct (uncoupled) regions each defining a 
c                    collective motion of at least two hinge joints.
c  11) clst()    --> assignment of unique rigid cluster labels to BULK sites. 
c  12) torsion() --> list of all hinge joints between rigid clusters
c  13) stress()  --> assignment of unique Laman subgraph labels to BULK sites. 
c  14) label_hinge()  --> subgraph assignment of unique collective motion labels  
c ------------------------------------------------------------------------------
      include   'set_parameter'
      integer   s,tag,tcmol,th_max,th_min
      integer   linkcf(nb2max),point(maxatm),torsion(2,nbmax)
      integer   multcf(maxatm)
      integer   pebble(6,0:maxatm)
      integer   block(-1:maxatm),shell(0:maxatm),btrack(0:maxatm)
      integer   clst(maxatm),stress(-1:maxatm)
      integer   link_hinge(nb2max),label_hinge(nbmax)               
      common/search/   btrack
      common/numbers1/ nflop,n_one,n_two,nbody,nclst,nstress
      common/numbers2/ iflop,nhinge,nih,ncmode
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/topology/ natom,linkcf,point,multcf
      common/subgraph/ stress
      common/cmotions/ label_hinge,link_hinge
      common/clusters/ clst
      stress(-1) = -1
      stress(0) = 0
         do s=1,natom
         clst(s) = s
         stress(s) = s
         enddo
      iflop = 0
      kflop = 0
      nflop =0
      call pebgam_3dbb(nflop)
      call decomp_R3dbb(n_one,n_two,nclst)
      call decomp_H3dbb(nhinge)
      call decomp_L3dbb(natom,tag,block,stress,nstress)
      kflop = nflop
      call decomp_C3dbb(natom,nih,ncmode,kflop)
      kflop = kflop - 3*n_one - 5*n_two
      nbody = kflop/6
      iflop = nflop - 3*n_one - 5*n_two - 6*nbody
      return
      end
