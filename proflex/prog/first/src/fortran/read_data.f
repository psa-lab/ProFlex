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
      subroutine read_data(firstdata)
c ------------------------------------------------------------------------------
c                                             PROGRAM WRITTEN:    April 11, 1999
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c                                Last modified on:   March 14, 2001 by A.J. Rader
c ------------------------------------------------------------------------------
c                               Description
c This subroutine reads the <file>_FIRSTdataset file to collect all information
c relevant to the chemistry and connectivity of the macromolecule.
c ------------------------------------------------------------------------------
c INPUT: 
c    1) <file>_FIRSTdata
c ------------------------------------------------------------------------------
c OUTPUT: 
c    1) linknoh(), pointer(), multnoh(), tcmol, torsion()
c    2) linkref(), point(),   multref()
c    3) hb_bond(,),  hb_energy(), 
c ------------------------------------------------------------------------------
c                              Variable List 
c aname(s) = A four character name for atom s. 
c bval(s) = the "B-value" or the "temperature factor" found in the PDB file.
c chainid(g) = defines the chain that group g belongs to. 
c rec(s) = Single character description of PDB record, (eg. A => ATOM  or 
c          H => HETATM)
c freq(s) = gives the occupancy value found in the PDB records.
c g,g_old = a distinct group number for a known or unknown group of atoms. 
c gname(g) = three character name for group g. 
c gnum(g) = The group number for group g. (In general group numbers
c           are repeated when there are more than one chain.)
c idres(g) = identification of residue or group.
c linknoh(pointer(so)+j) = sf such that <so,sf> is a central-force bond in the
c     macromolecule and atom, sf, is the j-th neighbor to site so. 
c linkref() is used as a temporary holding array.
c locgroup(s) = distinct group number, g, for which atom s belongs.
c multnoh(so) = number of nearest central-force neighbors to atom so.
c natom = number of atoms in the macromolecule. These atoms are listed in the
c     <file>_FIRSTdataset  data file, although not necessarily in the original 
c     .pdb file because Hydrogen atoms may have been added. 
c nbonds = number of covalent bonds in the network.
c ngrp = total number of groups in the macromolecule. 
c nhet = total number of HETATM-atoms in the macromolecule. 
c nwater = total number of water molecules floating around.
c maxatm = maximum number of atoms that can be handled.
c maxgrp = maximum number of distinct groups of atoms allowed.
c nbmax = maximum number of central-force bonds that can be handled.
c nb2max = 2*nbmax = twice the maximum # of CF bonds that can be handled.
c nsnoh = natom = total number of atoms (sites) in macromolecule (network).
c pointer(so) = an index for array linknoh() to get connectivity information.
c     No hydrogen bonds (noh) are included at this stage. 
c s,so,sf = atom labels. 
c slow(g),sbig(g) = the range of site labels (slow .LE. s .LE. sbig) that
c     belong to group g.
c tcmol = number of torsion constraints to be applied within the molecule,
c     usually due to a peptide or resonant bond.
c th_max = treated here as a dummy variable: part of "dihedral" common group.
c th_min = treated here as a dummy variable: part of "dihedral" common group.
c torsion(j,tc) = <so-sf>  for j=1,2 respectively. Defines a pair of atoms where
c      a torsion constraint is to be placed.
c s,smax = atom label. 
c slow(g),sbig(g) = the range of site labels (slow .LE. s .LE. sbig) that
c                   belong to group g. 
c ------------------------------------------------------------------------------
      include      'set_parameter'
      integer      g,g_old,hb,s,sf,smax,so,tcmol,th_max,th_min
      integer      gnum(maxgrp),slow(maxgrp),sbig(maxgrp)
      integer      locgroup(maxatm),pointer(maxatm),point(maxatm)
      integer      linknoh(nb2max),linkref(nb2max),torsion(2,nbmax)
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      integer      multnoh(maxatm),multref(maxatm),nhp
      dimension    sxyz(3,maxatm),freq(maxatm),bval(maxatm)
      dimension    hb_energy(maxh)
      character*80 firstdata
      character*77 record
      character*20 ctemp
      character*5  name3,idres(maxgrp)
      character*4  name1,aname(maxatm)
      character*3  name2,gname(maxgrp)
      character*1  rec(maxatm),chainid(maxgrp),hb_type(maxh)
      character*1  cid
      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval
      common/groups/   ngrp,nres,nhet,nwater,chainid,
     &                 idres,gname,gnum,slow,sbig
      common/network0/ nsnoh,linknoh,pointer,multnoh
      common/topology/ ns,linkref,point,multref
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
c ------------------------------------------------------------ format statements
 6000 format(a77)
 6005 format(5x,'ERROR in subroutine read_data.f')
 6010 format(5x,'maximum number of groups exceeded:   maxgrp = ',i7)
 6015 format(5x,'increase  maxgrp  in the  set_parameter  file')
 6020 format(5x,'maximum number of atoms exceeded:   maxatm = ',i8)
 6025 format(5x,'increase  maxatm  in the  set_parameter  file')
 6030 format(5x,'STOP: Please check original pdb file:  ',
     &          'No groups found!')
 6035 format(5x,'ERROR:  valence of atom ',i5,' is ',
     &          'greater than 4')
 6040 format(5x,' Atom  ','Group',' Atom# ',' Group#',' Chain')
 6045 format(6x,    a4,   3x,a3,  2x,  i5,2x,    a5,4x,   a1)
 6050 format(5x,'maximum number of central force constraints ',
     &          'exceeded:  nbmax = ',i8)
 6055 format(5x,'increase  nbmax  in the  set_parameter  file')
 6060 format(5x,'Self connection ignored!')
 6065 format(5x,'Duplicate connection ignored!')
 6070 format(5x,'HYDROGEN BOND ---> Self connection ignored!')
 6075 format(5x,'HYDROGEN BOND ---> Duplicate connection ignored!')
 6080 format(5x,'Hydrogen and covalent BOND overlap detected!')
 6085 format(5x,'WARNING:  Low number of hydrogen bonds found! ',
     &                                         'H-bond # =',i5)
 6090 format(5x,'Type  c  to continue:')
 6095 format(5x,'WARNING: unidentified type of hydrogen bond ')
 6100 format(5x,'Donor    =',i7)
 6105 format(5x,'Hydrogen =',i7)
 6110 format(5x,'Acceptor =',i7)
 6115 format(5x,'ERROR: Current hydrogen bond ID =',i6,
     &                                         ' (set ID > 0)')
 6120 format(5x,'ERROR: Current hydrogen bond ID =',i6)
 6121 format(5x,'       ID should be consecutive between [1,maxh]')
 6122 format(5x,'       maxh =',i6,' Can changed in parameter file')
 6125 format(5x,'ERROR: hydrogen bond ID is not unique.  ID =',i7)
      nhet = 0
      ngrp = 0
      tcmol = 0
      nwater = 0
      nbonds = 0
      nhp = 0
      g_old = -100000
      smax = -100000
         do hb=1,maxh
         hb_pick(hb) = -1
         enddo
c --------------------------------------------------------------- open PDB files
      open(1,file=firstdata,status='old')
      rewind(1)
c ------------------------------------------------------ initialize multiplicity
         do s=1,maxatm
         multnoh(s) = 0
         enddo
c -------------------------------------------------------------------- read data
      iflag = 1
   10 read(1,6000,end=11) record 
c ----------------------------------------------------------------- ATOM records
cTOM      1  N   MET A   1      -6.017 -15.110  26.296  1.00 41.45 
cccccc12345_name_RES_C1234x___4321.1234312.1234321.123321.12321.12_1234
c ------------------------------------------------------------------------------
         if( record(1:6) .eq. 'ATOM  ' ) then
         read(record,
     &   '(6x,i5,1x,a4,1x,a3,1x,a1,a5,3x,3f8.3,2f6.2,1x,i4)')
     &   s,name1,name2,cid,name3,x,y,z,occupny,bvalue,g
            if( s .gt. smax ) then
            smax = s
               if( smax .gt. maxatm ) then
               write(6,6005)
               write(6,6020) maxatm
               write(6,6025)
               stop
               endif
            endif
c -------------------------------------------------- record information on group
               if( g .ne. g_old )  then
               g_old = g
               nres = nres + 1
               ngrp = ngrp + 1
               slow(ngrp) = s
               gnum(ngrp) = g
               gname(ngrp) = name2
               chainid(ngrp) = cid
               idres(ngrp) = name3
               endif
            if( ngrp .gt. maxgrp ) then
            write(6,6005)
            write(6,6010) maxgrp
            write(6,6015)
            stop
            endif
         rec(s) = 'A'
         aname(s) = name1
         sxyz(1,s) = x
         sxyz(2,s) = y
         sxyz(3,s) = z
         freq(s) = occupny
         bval(s) = bvalue
         locgroup(s) = ngrp
c --------------------------------------------------------------- HETATM records
         elseif( record(1:6) .eq. 'HETATM' ) then
         nhet = nhet + 1
c ------------------------------------------------- convert to standard notation
            if(  record(18:20) .eq. 'HOH'
     &      .OR. record(18:20) .eq. 'H2O'
     &      .OR. record(18:20) .eq. 'OH2'
     &      .OR. record(18:20) .eq. 'WAT'
     &      .OR. record(18:20) .eq. 'DOD'
     &      .OR. record(18:20) .eq. 'D2O'
     &      .OR. record(18:20) .eq. 'h2o'
     &      .OR. record(18:20) .eq. 'hoh'
     &      .OR. record(18:20) .eq. 'wat'
     &      .OR. record(18:20) .eq. 'dod'
     &      .OR. record(18:20) .eq. 'd2o'  ) then
            nwater = nwater + 1
            record(18:20) = 'HOH'
c Added the following line for use with hbdilute programs. BMH 3.30.00
            record(22:22) = 'W'
            endif
         read(record,
     &   '(6x,i5,1x,a4,1x,a3,1x,a1,a5,3x,3f8.3,2f6.2,1x,i4)')
     &   s,name1,name2,cid,name3,x,y,z,occupny,bvalue,g
            if( s .gt. smax ) then
            smax = s
               if( smax .gt. maxatm ) then
               write(6,6005)
               write(6,6020) maxatm
               write(6,6025)
               stop
               endif
            endif
c -------------------------------------------------- record information on group
               if( g .ne. g_old )  then
               g_old = g
               ngrp = ngrp + 1
               slow(ngrp) = s
               gnum(ngrp) = g
               gname(ngrp) = name2
               chainid(ngrp) = cid
               idres(ngrp) = name3
               endif
            if( ngrp .gt. maxgrp ) then
            write(6,6005)
            write(6,6010) maxgrp
            write(6,6015)
            stop
            endif
         rec(s) = 'H'
         aname(s) = name1
         sxyz(1,s) = x
         sxyz(2,s) = y
         sxyz(3,s) = z
         freq(s) = occupny
         bval(s) = bvalue
         locgroup(s) = ngrp
c ----------------------------------------------------- read CF bond connections
         elseif( record(1:9) .eq. 'REMARK:CF' ) then
         read(record,'(9x,2i6)') so,sf
         multnoh(so) = multnoh(so) + 1
         multnoh(sf) = multnoh(sf) + 1
         nbonds = nbonds + 1
c ------------------------------------------------ ERROR check if too many bonds
            if( nbonds .gt. nbmax ) then
            write(6,*)
            write(6,6005) 
            write(6,6050) nbmax
            write(6,6055) 
            stop
            endif
         k2 = 2*nbonds
         k1 = k2 - 1
         linkref(k1) = so
         linkref(k2) = sf
c ---------------------------------------------------------- read TF constraints
         elseif( record(1:9) .eq. 'REMARK:TF' ) then
         read(record,'(9x,2i6)') so,sf
         tcmol = tcmol + 1
         torsion(1,tcmol) = so
         torsion(2,tcmol) = sf
c ---------------------------------------------------------- read HB constraints
         elseif( record(1:9) .eq. 'REMARK:HB' ) then
         nhb = nhb + 1
c000000001111111111222222222233333333334444444444555555555566666666667
c234567890123456789012345678901234567890123456789012345678901234567890
cxxxxxxxxiiiiibfffff.fffffbiiiiiiibiiiiiiibiiiiiiibbbbcccccccccccccccccccc
cEMARK:HB  446     -9.00000    1606    1607    1521    SB no energy
cEMARK:HB  447    -9.00000     159     161    1927    SB no energy
cEMARK:HB  448     -9.00000    1606    1608    1521    SB no energy
         read(record,'(9x,i5,1x,f12.5,3(1x,i7),4x,a20)') 
     &                 id_hb,   ftemp, so,s,sf,   ctemp
         hbond(1,nhb) = so
         hbond(2,nhb) = s 
         hbond(3,nhb) = sf
            if( ctemp(1:2) .eq. 'HB' ) then
            hb_type(nhb) = 'H'
            elseif( ctemp(1:2) .eq. 'SB' ) then
            hb_type(nhb) = 'I'
c ------------ added next two lines to accomodate hydrophobic tethers AJR 03.14.01
            elseif( ctemp(1:2) .eq. 'PH' ) then
            hb_type(nhb) = 'B'
            nhp = nhp +1
            else
            iflag = -1
            write(6,6095) 
            write(6,6100) so
            write(6,6105) s
            write(6,6110) sf
            hb_type(nhb) = 'U'
            endif
         hb_energy(nhb) = ftemp
            if( id_hb .lt. 1 ) then
            write(6,*)
            write(6,6115) id_hb
            stop
            elseif( id_hb .gt. maxh ) then
            write(6,*)
            write(6,6120) id_hb
            write(6,6121) 
            write(6,6122) maxh
            stop
            elseif( hb_pick(id_hb) .gt. 0 ) then
            write(6,*)
            write(6,6125) id_hb
            stop
            else
            hb_pick(id_hb) = 1
            endif
         hb_id(nhb) = id_hb
         endif
c ------------------------------------------------------- skip all other records
      goto 10
   11 natom = smax
      close(1)
c -------------------------------------------------- No groups found in PDB file
         if( ngrp .eq. 0 ) then
         write(6,*)
         write(6,6030)
         write(6,*)
         stop
         endif
c ----------------------------------------------- record last atom in each group
      sbig(ngrp) = natom
         do g=1,ngrp-1
         sbig(g) = slow(g+1) - 1
         enddo
c ----------------------------------------------- define pointer() for linknoh()
      pointer(1) = 0
         do s=1,natom-1
         pointer(s+1) = pointer(s) + multnoh(s)
         multnoh(s) = 0
         enddo
      multnoh(natom) = 0
c ----------------------------------------------------- construct neighbor table
      k2max = 2*nbonds
      do 50 k2=2,k2max,2
      k1 = k2 - 1
      so = linkref(k1)
      sf = linkref(k2)
c ----------------------------------------------------- self-connect ERROR check
         if( so .eq. sf ) then
         write(6,*) 
         write(6,6060) 
         write(6,6040) 
         go = locgroup(so)
         write(6,6045) aname(so),gname(go),so,idres(go),chainid(go)
         goto 50
         endif
c ------------------------------------------------------- DEGENERACY ERROR check
      index = pointer(so)
         do j=1,multnoh(so)
         index = index + 1
         s = linknoh(index)
            if( s .eq. sf ) then
            write(6,*) 
            write(6,6065) 
            write(6,6040) 
            go = locgroup(so)
            gf = locgroup(sf)
            write(6,6045)aname(so),gname(go),so,idres(go),chainid(go)
            write(6,6045)aname(sf),gname(gf),sf,idres(gf),chainid(gf)
            goto 50
            endif
         enddo
c --------------------------------------------------- augment connectivity table
      multnoh(so) = multnoh(so) + 1
      multnoh(sf) = multnoh(sf) + 1
      linknoh( pointer(so) + multnoh(so) ) = sf
      linknoh( pointer(sf) + multnoh(sf) ) = so
   50 continue
c ---------------------------------------------- test hydrogen bond connectivity
      nedges = 0
         do s=1,maxatm
         multref(s) = multnoh(s)
         nedges = nedges + multnoh(s)
         enddo
      nedges = nedges/2 + nhb
c ------------------------------------------------ ERROR check if too many edges
         if( nedges .gt. nbmax ) then
         write(6,*)
         write(6,6005) 
         write(6,6050) nbmax
         write(6,6055) 
         stop
         endif
c ---------------------------------------------------------- expand multiplicity
      do hb=1,nhb
      s = hbond(2,hb)
      sf = hbond(3,hb)
      multref(s) = multref(s) + 1
      multref(sf) = multref(sf) + 1
      enddo
c ------------------------------------------------- define point() for linkref()
      point(1) = 0
         do s=1,natom-1
         point(s+1) = point(s) + multref(s)
         multref(s) = 0
         enddo
      multref(natom) = 0
c ----------------------------- construct neighbor table for hydrogen bonds only
      do 100 hb=1,nhb
      so = hbond(2,hb)
      sf = hbond(3,hb)
c ----------------------------------------------------- self-connect ERROR check
         if( so .eq. sf ) then
         write(6,*) 
         write(6,6070) 
         write(6,6040) 
         go = locgroup(so)
         write(6,6045) aname(so),gname(go),so,idres(go),chainid(go)
         goto 100
         endif
c ------------------------------------------------------- DEGENERACY ERROR check
      index = point(so)
         do j=1,multref(so)
         index = index + 1
         s = linkref(index)
            if( s .eq. sf ) then
            write(6,*) 
            write(6,6075) 
            write(6,6040) 
            go = locgroup(so)
            gf = locgroup(sf)
            write(6,6045)aname(so),gname(go),so,idres(go),chainid(go)
            write(6,6045)aname(sf),gname(gf),sf,idres(gf),chainid(gf)
            goto 100
            endif
         enddo
c --------------------------------------------------- augment connectivity table
      multref(so) = multref(so) + 1
      multref(sf) = multref(sf) + 1
      linkref( point(so) + multref(so) ) = sf
      linkref( point(sf) + multref(sf) ) = so
  100 continue
c --------------------------------------------------- augment the covalent bonds 
      do so=1,natom
      ipo = point(so)
      index = pointer(so)
         do jo=1,multnoh(so)
         index = index + 1
         sf = linknoh(index)
            if( sf .gt. so ) then
c ------------------------------------------------------- DEGENERACY ERROR check
               do j=1,multref(so)
               ip = ipo + j
               s = linkref(ip)
                  if( s .eq. sf ) then
                  write(6,*) 
                  write(6,6080) 
                  write(6,6040) 
                  go = locgroup(so)
                  gf = locgroup(sf)
                  write(6,6045) aname(so),gname(go),so,
     &                            idres(go),chainid(go)
                  write(6,6045) aname(sf),gname(gf),sf,
     &                            idres(gf),chainid(gf)
                  goto 150
                  endif
               enddo
c --------------------------------------------------- augment connectivity table
            multref(so) = multref(so) + 1
            multref(sf) = multref(sf) + 1
            linkref( point(so) + multref(so) ) = sf
            linkref( point(sf) + multref(sf) ) = so
  150       continue
            endif
         enddo
      enddo
c ---------------------------------------------------- fix up hydrogen bond type
         do hb=1,nhb
            if( hb_type(hb) .eq. 'H' ) then
            so = hbond(2,hb)
            sf = hbond(3,hb)
               if( (aname(so)(2:2) .eq. 'S') 
     &        .OR. (aname(sf)(2:2) .eq. 'S') ) then
               hb_type(hb) = 'S'
               endif
            endif
         enddo
c ---------------------------- sort hydrogen bonds from highest to lowest energy
         if( nhb .gt. 1 ) then
         call sort_HB(nhb,1)
         endif
c ----------------------------------------------- define miscellaneous constants
      khb = nhb
c ------------------ DEBUG h-bond info output  - SN 2006:03 ------------
c      write (6,*) khb
c      write (6,*) nhp
c --------------------------- END DEBUG --------------------------------
      ns = natom
      nsnoh = natom
      return
      end

      subroutine sort_HB(low,high)
c ------------------------------------------------------------------------
      include         'set_parameter'
      integer         high,high2
      integer         hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      dimension       hb_energy(maxh)
      character*1     hb_type(maxh),ctemp
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
	itest = high - low
	if(itest .eq. 0) return
	istep = itest/iabs(itest)
	high2 = high - istep
        flag = hb_energy(low)
           do i=low,high,istep
           if( flag .lt. hb_energy(i) ) flag = hb_energy(i)
           enddo
        flag = flag + 100.0e0 
        do i=low,high,istep
        temp = flag
           do j=low,high2,istep
           k = j + istep
              if( hb_energy(j) .gt. hb_energy(k) ) then
              temp = hb_energy(j)
              hb_energy(j) = hb_energy(k)
              hb_energy(k) = temp
              itemp = hb_id(j)
              hb_id(j) = hb_id(k)
              hb_id(k) = itemp
c              itemp = hb_pick(j)
c              hb_pick(j) = hb_pick(k)
c              hb_pick(k) = itemp
              ctemp = hb_type(j)
              hb_type(j) = hb_type(k)
              hb_type(k) = ctemp
                 do l=1,3
                 itemp = hbond(l,j)
                 hbond(l,j) = hbond(l,k)
                 hbond(l,k) = itemp
                 enddo
              endif
           enddo
        if(temp .eq. flag) return
        enddo
	return
	end
