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
      subroutine out_PDBfile(fgrphcs,smap,matom,clst,nclst,wt)
c      subroutine out_PDBfile(fgrphcs,smap,matom,
c     &                       clst,pointbrc,multbrc,nclst)
c ------------------------------------------------------------------------------
c AJR 03.27.02 eliminated the use of CFB (BOND) atoms and CONECT records
c also added the flexibility index into the b-value column.
c ------------------------------------------------------------------------------
c                                             LAST UPDATED:       April 12, 2002
c                                              revisions by:      A.J. Rader
c                                                                rader@pa.msu.edu
c                                             PROGRAM WRITTEN:     JUNE 23, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     This subroutine creates a PDB formated file in which RasMol can read. It
c places an atom between every pair of atoms connected by a central-force bond.
c In addition, it connects all donor-acceptor pairs that form H-bonds that have
c been placed. The group number for the fictitious atoms, bisecting each CF
c bond, is set as 9999, the group name is called CFB for central force bond and 
c the sites are named  BOND. The coordinate of each BOND site is in the middle 
c of the two atoms forming the CF bond. The labeling of the BOND atoms is done
c with the help of pointbrc() and multbrc(), both of which are previously 
c determined in the relabelatms() subroutine. The exact labels is not 
c prespecified, however, the grouping of labels is. 

  
c ------------------------------------------------------------------------------
c INPUT: 
c    1) output file name
c ------------------------------------------------------------------------------
c OUTPUT:  
c    1) creates  fgrphcs
c ------------------------------------------------------------------------------
c                              Variable List 
c pointbrc(),multbrc() used to define the BOND site indices. 
c smap(s) = new set of site labels defined from the original set to optimize the
c        coloring scheme by constructing a canonical coloring sequence.
c ==============================================================================
      include      'set_parameter'
      integer      g,s,sf,so
      integer      gnum(maxgrp),slow(maxgrp),sbig(maxgrp)
      integer      linkref(nb2max),point(nsmax),locgroup(nsmax)
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      integer      smap(matom),testwt
      integer      clst(nclst),pointbrc(nclst),multbrc(nclst)
      integer      multref(nsmax),torsion(nb2max) 
      dimension    sxyz(3,nsmax),freq(maxatm),bval(maxatm)
      dimension    w(3),hb_energy(maxh)
      real         zooc,wt(matom)
      character*80 fgrphcs,fwebpdb
      character*5  idres(maxgrp)
      character*4  aname(nsmax),anam
      character*3  gname(maxgrp),gnam
      character*1  rec(maxatm),chainid(maxgrp),idchain, hb_type(maxh)

      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval
      common/groups/   ngrp,nres,nhet,nwater,chainid,
     &                 idres,gname,gnum,slow,sbig
      common/topology/ ns,linkref,point,multref
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag

c ------------------------------------------------------------ format statements 
 5000 format(5x,'ERROR: Torsion constraint inconsistency ',
     &          'detected in out_PDBfile.f')


 6000 format('ATOM  ',i5,1x,a4,1x,a3,1x,a1,a5,3x,3f8.3)
 6005 format('HETATM',i5,1x,a4,1x,a3,1x,a1,a5,3x,3f8.3)
 6001 format('ATOM  ',i5,1x,a4,1x,a3,1x,a1,a5,3x,3f8.3,f6.2,i6)
 6006 format('HETATM',i5,1x,a4,1x,a3,1x,a1,a5,3x,3f8.3,f6.2,i6)

c 6000 format('ATOM  ',i5,1x,a4,1x,a3,1x,a1,a5,3x,3f8.3,2f6.2) 
c 6006 format('HETATM',i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
 6010 format('CONECT',i5,i5)
      
      anam = 'BOND'
      gnam = 'CFB'
      idchain = ' '
      
      nwts  = 0
      wtmin = 0.0
      wtmax = 0.0
      zocc  = 1.00
      
      open(3, file=fgrphcs,status='unknown')
      rewind(3)
c ---------------------------------------------------- output all atoms
      do g=1,ngrp
         do s=slow(g),sbig(g)
            testwt = int(1000*wt(s))
            if( rec(s) .eq. 'H') then
               write(3,6006) smap(s),aname(s),gname(g),chainid(g),
     &              idres(g),(sxyz(j,s), j=1,3),zocc,testwt
            else
               write(3,6001) smap(s),aname(s),gname(g),chainid(g),
     &              idres(g),(sxyz(j,s), j=1,3),zocc,testwt
            endif
         enddo
      enddo

      close(3)
      return
      end







