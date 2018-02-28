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
      subroutine place_Hbond()
c ------------------------------------------------------------------------------
c                                             LAST UPDATED:        June 25, 1997
c                                             PROGRAM WRITTEN:     June  8, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c The linknoh(), pointer() and multnoh() arrays are used in addition to the list 
c of Hydrogen bonds to determine the complete topology of the macromolecule. 
c ------------------------------------------------------------------------------
c INPUT: 
c    1) khb 
c    2) nsnoh,linknoh,pointer,multnoh 
c ------------------------------------------------------------------------------
c OUTPUT: 
c    1) linkref,point,multref 
c ------------------------------------------------------------------------------
c                              Variable List 
c hbond(j,nhb) = <so-s-sf>  for j=1,2,3 respectively. so => donor atom label, 
c     s => Hydrogen label and sf => acceptor atom label.
c khb = number of Hydrogen bonds placed in the molecule.
c linknoh() = network connectivity information (nearest neighbor table) of the 
c     macromolecule when NO Hydorogen bonds are present.
c linkref() = network connectivity information (nearest neighbor table) of the 
c     macromolecule with (khb) Hydorogen bonds are present.
c maxatm = maximum number of atoms allowed in the macromolecule.
c maxh = maximum number of H-bonds that can be considered.
c multnoh(s) = number of CF bonds incident to site s when no H-bonds are present
c multref(s) = number of CF bonds incident to site s when H-bonds are present.
c nb2max = twice the maximum number of central-force bonds in the BB-network.
c nhb = maximum number of [potential] Hydrogen bonds.
c nsnoh,ns = number of atoms (sites) in the macromolecule (bond-bending network)
c pointer() = used to give appropriate index in linknoh()
c point() = used to give appropriate index in linkref()
c ==============================================================================
      include      'set_parameter'
      integer      hb,s,so,sf
      integer      linknoh(nb2max),linkref(nb2max)
      integer      pointer(maxatm),point(maxatm) 
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      integer      multnoh(maxatm),multref(maxatm)
      dimension    hb_energy(maxh)
      character*1  hb_type(maxh)
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/network0/ nsnoh,linknoh,pointer,multnoh
      common/topology/ ns,linkref,point,multref

 3333 format(i8,2x,i8)
 3334 format('Atom: ',i8,' has incorrect multiplicity of ',i8)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     The topology of the system is stored in a linear array,
c     either "linkref" for the entire system, or "linknoh"
c     for the system without hydrogen bond constraints added.
c     The array "point(so)" points to the position in the linear
c     array linkref where the atoms connected to so are listed.
c     For example, point(33) may resolve to 421, which can be
c     used to index the array linkref. linkref(421) lists the
c     atom number of an atom connected to atom 33. If atom 33
c     has 4 bonds, then multref(33) will equal 4. Therefore, 
c     to find the atom numbers of the atoms connected to atom 33
c     we access linkref() from point(33) to point(33) + multref(33); 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     "multref(s)" tells us how many bonds atoms "s" has.
c     "multnoh(s)" is the array before we add hydrogens, so
c     we need to copy it to the complete array, "multref(s)".
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do s = 1, nsnoh
         multref(s) = multnoh(s)
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     For each hydrogen bond, increment the multiplicity (the 
c     number of bonds each atom makes).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do hb = 1,khb
         s  = hbond(2,hb_pick(hb))
         sf = hbond(3,hb_pick(hb))
         multref(s)  = multref(s)  + 1
         multref(sf) = multref(sf) + 1
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Augment the "point()" array to account for the new 
c     bonds added by the hydrogen bonds. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      point(1) = 0
      do s = 1,nsnoh-1
         point(s+1) = point(s) + multref(s)
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copy linknoh into linkref.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do so = 1, nsnoh
         ipo = pointer(so) 
         ipoint = point(so)
         do j = 1, multnoh(so)
            linkref(ipoint +j) = linknoh(ipo +j)
         enddo
         multref(so) = multnoh(so)
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Add H-Acceptor bonds
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do hb=1,khb
         s  = hbond(2,hb_pick(hb))
         sf = hbond(3,hb_pick(hb))
         multref(s)  = multref(s)  + 1
         multref(sf) = multref(sf) + 1
         linkref( point(s)  + multref(s) )  = sf
         linkref( point(sf) + multref(sf) ) = s 
      enddo
      ns = nsnoh

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     AJR 02.28.02 error check max multiplicity. A given atom 
c     can't have more than "izcf" bonds. The constant "izcf" 
c     is set in the "set_parameter" file. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do s = 1, ns
         if(multref(s) .gt. izcf) then
            write(6,3334) s,multref(s)
            write(36,3333) s,multref(s)
c            stop
         endif
      enddo

      return
      end
