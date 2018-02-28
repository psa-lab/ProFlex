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


      subroutine prunesetup(avgR,xclst,xxf,mode)
c ==============================================================================
c                                             PROGRAM WRITTEN:  Aug. 26, 2000
c                                              Last Modified :  Oct.  9, 2001
c                                       program written by:     A. J. Rader
c                                                               rader@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     This function, prunesetup(), and the one that follows, hbprune(ihb), use
c the same set of subroutines I wrote for the iterative stripping process that 
c I have been conducting on proteins.  These two functions return a real value 
c for  the average coordination after all dangling ends and non-hydrogen bonded 
c sidechain rings have been removed.  This leaves a network of atoms with 
c coordination of 2 or higher.  By removing the non-hydrogen bonded ring 
c sidechains, extra degrees of flexibility not associated with the bulk 
c protein are also removed.  
c ==============================================================================

c INPUT: 
c        common/numbers1/ nflop,n_one,n_two,nbody,nclst,nstress
c        common/numbers2/ iflop,nhinge,nih,ncmode

c  1) natom   --> number of atoms in the macromolecule
c  2) khb     --> number of hydrogen bonds realized in the macromolecule
c  3) tcmol   --> number of torsional constraints
c  4) nflop   --> number of independent degrees of freedom
c  5) nclst   --> number of rigid clusters 
c  6) n_one   --> number of isolated atoms  (no dimension;  a single point)
c  7) n_two   --> number of isolated dimers (one dimensional)
c  8) nbody   --> number of three dimensional objects
c  9) nstress --> number of network induced stressed regions
c 10) iflop   --> number of floppy modes (internal independent deg. of freedom)
c 11) nhinge  --> number of hinge joints in the macromolecule
c 12) nih     --> number of independent hinge joints
c 13) ncmode  --> number of collective motions having two or more hinges
c ------------------------------------------------------------------------------
c OUTPUT: 
c    1) writes basic INPUT numbers (see above) to screen. 
c    2) archive = Y or N  to save data files or not.
c       IF Y:
c    3) writes all output files (see above):  outfile(j)  for j=1 to j=11
c ------------------------------------------------------------------------------
c                              Variable List 
c aname(s) = a four character name of atom s.
c archive,answer0 = Y or N answers to questions from screen in interactive mode.
c calow(icolor) = the smallest "new atom label" assigned to the icolor-th color.
c cabig(icolor) = the biggest "new atom label" assigned to the icolor-th color.
c chainid(g) = defines the chain that group g belongs to.
c choice(icolor,ichoice) = color code with the icolor-th level of preference for
c        the ichoice-th arrangement of colors.
c clst(s) = cluster label for BULK site s. Note all sites act as a BULK site to
c        one and only one unique rigid cluster. 
c cname(i) = the name of the i-th color
c color(nc) = a color code for the nc-th rigid cluster.
c cslow(icolor) = the smallest "new site label" assigned to the icolor-th color.
c csbig(icolor) = the biggest "new site label" assigned to the icolor-th color.
c        Note that all so called site labels > natom, corresponding to bonds.  
c fdecomp,fgrphcs,fscript = output file names.
c g = a generic dummy index for a group.
c gname(g) = three character name of group g.
c gnum(g) = The group number for group g. (In general group numbers
c        are repeated when there are more than one chain.)
c h = an index for HET-atoms
c hb = an index for Hydrogen bond.  (1 .LE. hb .LE. khb)
c hbond(j,hb) = <so-s-sf>  for j=1,2,3 respectively. so => donor atom label,
c        s => Hydrogen label and sf => acceptor atom label.
c ichoice = a dummy index defining the type of color preference table.
c icolor = a dummy index defining a color of interest. (icolor=0 => colorless.) 
c khb = number of Hydrogen bonds that have been placed in molecule.
c linkref() = Central-force nearest neighbor table of the bond-bending network.
c        The j-th nearest neighbor to site so is given by: linkref(point(so)+j)
c locgroup(s) = distinct group number, g, for which atom s belongs.
c maxch = maximum number of preference sets to choose from.
c maxcolor = maximum number of different colors used.
c maxr =  the maximum number of different ranges.
c mult(r) = the number of rigid clusters within the r-th range of sizes.
c multarc(nc) = number of BULK sites contained in the nc-th rigid cluster. 
c multbrc(nc) = # of all NON-hinge bonds contained in the nc-th rigid cluster. 
c multref(s) = # of nearest neighbors that site s has.
c natom = number of atoms in the macromolecule. Not ghost sites. 
c nc = cluster label running from 1 to nclst.
c nc1 = the "new" cluster label such that  size(nc) = 1 for (nc .GE. nc1) 
c nc2 = the "new" cluster label such that  size(nc) .LE. 2 for (nc .GE. nc2) 
c nclst = maximum number of rigid clusters found in macromolecule. 
c        Note that (nclst .LE. maxatm).
c newlabel(s) = gives the new atom label (1 .le. s .le. natom) for simplifying 
c        graphics. Note the fgrphcs  PDB file consist of the new labels. 
c nhinge = number of hinge joints found in macromolecule. 
c nogood(icolor) =flag for determining which colors have already been considered
c nsc = number of lines to define selection criteria
c oldlabel() = inverse of newlabel() such that s = oldlabel( newlabel(s) ). 
c paint(iorder) = defines which color to paint.
c point() = used to give an appropriate index in linkref()
c pointarc(nc) = is used to give an appropriate index in determining BULK atoms 
c        such that: (pointarc(nc)+1 .le. s_new .le. pointarc(nc)+multarc(nc))
c pointbrc(nc) = is used to give an appropriate index in determining NON-hinge 
c        BOND sites: (pointbrc(nc)+1 .le. s_new .le. pointbrc(nc)+multbrc(nc))
c pointr(r) = a pointer index to determine all clusters within range r.
c r = an index for either a residue or for the range in rigid cluster size. 
c s,so,sf = site labels
c size(nc) = size of nc-th rigid cluster.
c slow(g),sbig(g) = the range of site labels (slow .LE. s .LE. sbig) that
c        belong to group g.
c sxyz(jj,s) = the jj-th component of the s-th site's coordinate.
c mode  = the mode this program is run in. If mode = 0, then  don't print to
c     pruned.dat the output.  If mode =1 then do print pruned.dat (ajr 11/5/01)
c ==============================================================================
      include      'set_parameter'
      integer      si,sf0,totcoord,numatom,xmsf,msf,lmsf,totatm
      integer      hbA,hbH,xA,xH,site,gD,gA
      integer      halfcoord,uF,ncoord(izcf),mode
      integer      g,h,hb,r,s,sf,so,tc,tcmol,th_max,th_min
      integer      gnum(maxgrp),slow(maxgrp),sbig(maxgrp)
      integer      linkref(nb2max),point(maxatm),locgroup(maxatm)
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      integer      multref(maxatm),pmult(maxatm),maxmult(maxatm)
      integer      plink(nb2max),pointer(maxatm),mbod(maxatm)
      integer      gflag(maxgrp),gp,ps1,ps2,co1,pnf(natom)
      integer      torsion(2,nbmax),clst(maxatm)
      integer      orgsites,orgpoint(maxatm),orglink(nb2max)
      integer      orgmult(maxatm),size(maxatm)
      integer      multarc(maxatm),multbrc(maxatm)
      real         avgR,aveR,xf,xxf,xclst
      dimension    sxyz(3,maxatm),freq(maxatm),bval(maxatm)
      dimension    hb_energy(maxh)
      character*5  idres(maxgrp)
      character*4  aname(maxatm),file_id,atm,atA,atD
      character*3  gname(maxgrp)
      character*1  rec(maxatm),chainid(maxgrp),hb_type(maxh)
      character*10 fname
      common/numbers1/ nflop,n_one,n_two,nbody,nclst,nstress
      common/numbers2/ iflop,nhinge,nih,ncmode
      common/clusters/ clst
      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval
      common/groups/   ngrp,nres,nhet,nwater,chainid,idres,
     &                 gname,gnum,slow,sbig
      common/topology/ ns,linkref,point,multref
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/ptopo/    pmult,plink,pointer,maxmult,gflag

 7200 format(9i8)
 7588 format(f9.6,2(2x,f11.9),2x,f9.5)
      fname ='pruned.dat'
      ihb = hb_pick(khb)
c ------------------------------------ make neighbor table with current hbonds
c     and store original topology into org* arrays
      orgsites = ns
      do so = 1,natom
         pmult(so) = multref(so)
         maxmult(so) = multref(so)
         pointer(so) = point(so)
         mbod(so) = 0
         orgmult(so) = multref(so)
         orgpoint(so) = point(so)
         do k = 1,multref(so)
            plink(pointer(so) +k) = linkref(point(so)+k)
            orglink(orgpoint(so) +k) = linkref(point(so)+k)
         enddo
      enddo
      do i = 1,izcf
         ncoord(i) = 0
      enddo
c -------------------------------------------------- initial pruning
c                                   sidechain rings: setup and removal
      do ig = 1,ngrp
         gflag(ig) = 0
      enddo
      do jhb = 1,khb
         lhb = hb_pick(jhb)
         do i=1,3,2
            g = locgroup(hbond(i,lhb)) 
            if(( gname(g) .eq. 'PHE').or. (gname(g) .eq. 'TRP') .or. 
     &           ( gname(g) .eq. 'HIS').or. (gname(g) .eq. 'TYR')) then
               atm = aname(hbond(i,lhb))
               if((atm .ne. ' N  ') .and.(atm .ne. ' O  ') .and. 
     &              (atm .ne. ' C  ') .and.(atm .ne. ' CA ') .and. 
     &              (atm .ne. ' OXT')) gflag(g) = gflag(g)+1
c -------------------------(added 05.02.01 ajr)----------------------------------
c                             to cope with hydrophobic tether pseudoatoms
c                            NOTE: must change code if # of psuedoatoms increases
            elseif( gname(g) .eq. 'XXX') then
               ps2 = hbond(1,lhb)
               ps1 = plink(pointer(ps2)+1)
               if(ps1 .eq. hbond(2,lhb)) ps1 = plink(pointer(ps2)+2)
               co1 = plink(pointer(ps1)+1)
               if(co1 .eq. ps2) co1 = plink(pointer(ps1)+2)
               g = locgroup(co1)
               if(( gname(g) .eq. 'PHE').or. (gname(g) .eq. 'TRP') .or.
     &            ( gname(g) .eq. 'HIS').or. (gname(g) .eq. 'TYR')) then
                  atm = aname(hbond(i,lhb))
                  if((atm .ne. ' N  ') .and.(atm .ne. ' O  ') .and. 
     &                 (atm .ne. ' C  ') .and.(atm .ne. ' CA ') .and.
     &                 (atm .ne. ' OXT')) gflag(g) = gflag(g)+1
               endif
c----------------------(added 05.02.01 ajr)-------------------------------------
            endif
         enddo
      enddo

      do ig=1,ngrp
         if(((gname(ig) .eq. 'PHE').or.(gname(ig) .eq. 'TRP') .or. 
     &        ( gname(ig) .eq. 'HIS').or.(gname(ig) .eq. 'TYR'))
     &        .and. (gflag(ig).eq.0) ) gflag(ig) = -1
      enddo
      
      do si=1,natom
         ncoord(pmult(si)) = ncoord(pmult(si))+1
      enddo
      do j = 1,izcf
         totatm= totatm +ncoord(j)
         totcoord = totcoord + j*ncoord(j)
      enddo
      call hbprunesc(aname,locgroup)
c                                    pruning dangling ends and cleanup
      call hbpruneit(natom)
      call hbcondenseMult(natom)
c ------------------------------------- identify connected sites
      do si = 1,natom
         pnf(si) =0
      enddo
      call pclazz(pnf,natom)
c -------------------------------------- calculate number of bodies
      ibod = 0
      maxbod =0
      do si = 1,natom
         mbod(si) = 0
      enddo
      do si =1,natom
         mbod(pnf(si))=mbod(pnf(si))+1
      enddo

      do si = 1,natom
         if(mbod(si) .gt. 0) then
            if(mbod(si).gt. maxbod) maxbod = mbod(si)
            ibod=ibod+1
         endif
      enddo
c -------------------------------- removing smaller isolated bodies
      if(ibod .gt.1) then
         do so =1,natom
            if(mbod(pnf(so)).eq. 1) then
               do i=1,pmult(so)
                  plink(pointer(so)+i) = -1
               enddo
               pmult(so) = 0
               ibod=ibod-1
            endif
         enddo
      endif
c -----------------------------------------error check
      if(ibod .gt. 1) write(18,2000) ibod
 2000 format(5x,'Abandon all hope, the number of bodies ',i4,' is gt 1')
c ------------------------------------------------- initialize
      do i=1,izcf
         ncoord(i) = 0
      enddo
      totcoord =0
      totatm =0
c --------------------- check on number of torsional constraints
      jtcmol = 0
      do tc = 1,tcmol
         so = torsion(1,tc)
         sf = torsion(2,tc)
         if((pmult(so).ne.0).and.(pmult(sf).ne.0)) then
            do lt =1,pmult(so)
               si = plink(pointer(so)+lt)
               if(si.eq.sf) jtcmol = jtcmol +1
            enddo
         endif
      enddo
c ------------------------------------------------ calculate <r>
      do si=1,natom
         ncoord(pmult(si)) = ncoord(pmult(si))+1
      enddo
      do j = 1,izcf
         totatm= totatm +ncoord(j)
         totcoord = totcoord + j*ncoord(j)
      enddo
      aveR = float(totcoord)/(float(totatm))
      xpep = float(jtcmol)/float(totatm)
      avgR = (aveR+10.0*xpep)/(1.0+4.0*xpep)
c ======================================================
c AJR 10.05.01 changes made to run the pebblegame twice 
c ----------------------------------- run pebblegame on skelton network
      do so = 1,natom
         multref(so) = pmult(so)
         point(so) = pointer(so)
         do k =1,orgmult(so)
c  must run the entire range to make sure that empty spots get mapped too       
c         do k = 1,pmult(so)
            linkref(point(so) +k) = plink(pointer(so)+k)
         enddo
      enddo
      ns = natom
      call mapmolecule()   
c --------------- calculate lgst_clst   
c ----------------------------- count # of BULK & SURFACE sites for each cluster
      do nc=1,nclst
         size(nc) = 0
         multarc(nc) = 0
         multbrc(nc) = 0
      enddo
      do so=1,natom
         ipo = point(so)
         nco = clst(so)
c ---------------------------------------------- count BULK sites in cluster nco
         multarc(nco) = multarc(nco) + 1
         do jo=1,multref(so)
            indexo = ipo + jo
            sf = linkref(indexo)
c ------------------------------------------------- prevent double counting bond
            if( sf .gt. so ) then
               ncf = clst(sf)
               if( nco .ne. ncf ) then
c ----------------------------------------------------------- count SURFACE site
                  size(nco) = size(nco) + 1
                  size(ncf) = size(ncf) + 1
               else
c --------------------------------------------------------- count BULK BOND-site
                  multbrc(nco) = multbrc(nco) + 1
               endif
            endif
         enddo
      enddo
c --------------------------------------------- calculate number of hinge joints
c                  2*nhinge= sum_{nc=1}^{nclst} [# of SURFACE sites in clst(nc)]
      nhinge = 0
      maxsize = -1
      do nc=1,nclst
         nhinge = nhinge + size(nc)
c -------------------------------------------------- size = SURFACE + BULK sites
         size(nc) = size(nc) + multarc(nc)
c ----------------------------------------------------------- find maximum size
         if( maxsize .lt. size(nc) ) maxsize = size(nc)
      enddo
      nhinge = nhinge/2
      lg_clst = maxsize
      num_clst = nclst
      xf = float(iflop)/(3.0*float(totatm))
      xxf = xf/(1.0+4.0*xpep)
      xclst = float(lg_clst)/float(totatm)

c      open(24,file=fname,access='append',status='unknown')
      if(mode .eq. 1) then
      open(24,file='pruned.dat',access='append',status='unknown')
      write(24,7588) avgR,xxf,xclst,hb_energy(hb_pick(khb))
      close(24)
      endif
c ---------------------------------- restore original topology
      ns = orgsites
      do so = 1,ns
         multref(so) = orgmult(so)
         point(so) = orgpoint(so)
         do k = 1,orgmult(so)
            linkref(point(so) +k) = orglink(orgpoint(so)+k)
         enddo
      enddo
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hbprune(ihb,avgR,xclst,xxf,mode)

      include      'set_parameter'
      integer      si,sf0,totcoord,numatom,xmsf,msf,lmsf,totatm
      integer      hbA,hbH,xA,xH,site,gD,gA,hbD
      integer      halfcoord,uF,ncoord(izcf)
      integer      g,h,hb,r,s,sf,so,tc,tcmol,th_max,th_min
      integer      gnum(maxgrp),slow(maxgrp),sbig(maxgrp)
      integer      linkref(nb2max),point(maxatm),locgroup(maxatm)
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      integer      multref(maxatm),pmult(maxatm),maxmult(maxatm)
      integer      plink(nb2max),ppnt(maxatm),mbod(maxatm)
      integer      gflag(maxgrp),gp,ps1,ps2,co1,pnf(natom)          
      integer      torsion(2,nbmax),newlabel(maxatm)
      integer      orgsites,orgpoint(maxatm),orglink(nb2max)
      integer      orgmult(maxatm),size(maxatm)
      integer      multarc(maxatm),multbrc(maxatm),clst(maxatm)
      real         aveR,xpep,avgR,xf,xxf,xclst
      dimension    sxyz(3,maxatm),freq(maxatm),bval(maxatm)
      dimension    hb_energy(maxh)
      character*5  idres(maxgrp)
      character*4  aname(maxatm),file_id,atm,atA,atD
      character*3  gname(maxgrp)
      character*1  rec(maxatm),chainid(maxgrp),hb_type(maxh)
      character*10 fname
*
      LOGICAL      TMPOP

      common/numbers1/ nflop,n_one,n_two,nbody,nclst,nstress
      common/numbers2/ iflop,nhinge,nih,ncmode
      common/clusters/ clst
      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval
      common/groups/   ngrp,nres,nhet,nwater,chainid,idres,
     &                 gname,gnum,slow,sbig
      common/topology/ ns,linkref,point,multref
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/ptopo/    pmult,plink,ppnt,maxmult,gflag

 7588 format(f9.6,2(2x,f11.9),2x,f9.5)
      fname = 'pruned.dat'
c ----------------------------------------- prune for given Hbond, ihb
      lhb = hb_pick(ihb)
      hbH = hbond(2,lhb)
      hbA = hbond(3,lhb)
      hbD = hbond(1,lhb)

      INQUIRE(FILE='tmpout',OPENED=TMPOP)
      if(.NOT.TMPOP) then
          open (77,FILE='tmpout',STATUS='UNKNOWN')
      endif
 1214 format (I6,2X,I6,2X,I6)
 1215 format ('gD: ',I4,', atD:',A4,', gflag: ',I4)
      write (77,*)
      write (77,*) 'H-bond atoms'
      write (77,1214) hbD,hbH,hbA


      if((pmult(hbH) .gt. 0).AND.(pmult(hbA).gt. 0)) then
         do iH = 1,pmult(hbH)
            if(plink(ppnt(hbH)+iH) .eq. hbA) xH = iH
         enddo
         plink(ppnt(hbH)+xH) = -1
         pmult(hbH) = pmult(hbH) -1

         do iA = 1,pmult(hbA)
            if( plink(ppnt(hbA)+iA) .eq. hbH) xA = iA
         enddo
         plink(ppnt(hbA)+xA) = -1
         pmult(hbA) = pmult(hbA) -1
         call hbcondenseMult(natom)

         if( pmult(hbH) .eq. 1) call hbsiteprune(hbH)
         if( pmult(hbA) .eq. 1) call hbsiteprune(hbA)
         call hbcondenseMult(natom)

         gA = locgroup(hbA)
         atA = aname(hbA)
         if( gflag(gA) .gt. 0) then
            if((atA.eq.' ND1').or.(atA.eq.' NE1').or.(atA.eq.' NE2')
     &           .or.(atA.eq.' OH ')) gflag(gA) = gflag(gA) -1
            if(gflag(gA) .eq. 0) call hbringbreak(gA,aname)
         endif


* Below, the gflag array was resulting in some out-of-bounds memory access
* for reasons unclear to me. A crude fix for this problem is the two write
* statements, which enable modification of gflag without any problem. 
* A possible reason for why the fix works is that the write statements
* allocate a temp buffer which helps avoid this problem.

         gD = locgroup(hbD)
         atD = aname(hbD)
         if( gflag(gD) .gt. 0) then
            if((atD.eq.' ND1').or.(atD.eq.' NE1').or.(atD.eq.' NE2')
     &           .or.(atD.eq.' OH ')) then
      write (77,1215) gD,atD,gflag(gD)
               gflag(gD) = gflag(gD) -1
      write (77,1215) gD,atD,gflag(gD)
               if(gflag(gD) .eq. 0) call hbringbreak(gD,aname)
            endif
         endif
         call hbcondenseMult(natom)
      endif

*      write (77,*) 'Made it this far-2'
c ------------------------------------- identify connected sites
      do si = 1,natom
         pnf(si) =0
      enddo
      call pclazz(pnf,natom)
c -------------------------------------- calculate number of bodies
      ibod = 0
      maxbod =0
      jmbod = 0
c ibod = number of bodies, jmbod = number of nontrivial bodies 
c                                      (i.e more than 1 atom)
      do si = 1,natom
         mbod(si) = 0
      enddo
c                                initialize body label,mbod(s)
      do si = 1, natom
         mbod(pnf(si))=mbod(pnf(si))+1
      enddo
      do si = 1,natom
         if(mbod(si) .gt. 0) then
            if(mbod(si).gt. maxbod) maxbod = mbod(si)
            ibod=ibod+1
c            if((mbod(si) .gt. 1).AND.( mbod(si) .lt. maxbod) ) then 
c               jmbod = jmbod +1
c            endif
         endif
      enddo
c -------------------------------- removing smaller isolated bodies (single sites)
      if(ibod .gt. 1) then
         do so =1,natom
            if(mbod(pnf(so)) .eq. 1) then
               do i=1,pmult(so)
                  plink(ppnt(so)+i) = -1
               enddo
               pmult(so) = 0
               ibod = ibod -1
            elseif(mbod(pnf(so)) .lt. maxbod) then
               jmbod = jmbod+1
c               write(23,*) lhb,so,pnf(so),mbod(pnf(so))
c               stop
            endif
         enddo
      endif
c ----------- AJR 03.25.02: remove small isolated rings (bodies)
      if(ibod .gt. 1) then
c         write(65,*) ihb,ibod,maxbod,jmbod
         do so = 1,natom 
            lbody = mbod(pnf(so))
            if((lbody .gt. 1) .AND.(lbody .ne. maxbod)) then
               do i= 1,pmult(so)
                  plink(ppnt(so) +i) = -1
               enddo
               pmult(so) = 0
               jmbod = jmbod -1
            endif
         enddo
      endif
      if(jmbod .eq. 0) ibod = 1

c        if( lhb .lt. nhb) pause 'Made it this far'
c -----------------------------------------error check
      if(ibod .gt. 1) write(17,2000) ibod
 2000 format(5x,'Abandon all hope, the number of bodies ',i4,' is gt 1')
      if(ibod .gt. 2) then
         write(65,*) ihb,ibod,maxbod,jmbod
c         write(14,*) natom, (plink(ppnt(natom)+i),i=1,pmult(natom))
         do s = 1,natom
            write(15,*) s,mbod(s),pnf(s),mbod(pnf(s)),pmult(s)
         enddo
         stop
      endif

c ------------------------------------------------- initialize
      do i=1,izcf
         ncoord(i) = 0
      enddo
      totcoord =0
      totatm = 0
c --------------------- check on number of torsional constraints
      jtcmol = 0
      do tc = 1,tcmol
         so = torsion(1,tc)
         sf = torsion(2,tc)
         if((pmult(so).ne.0).and.(pmult(sf).ne.0)) then
            do lt =1,pmult(so)
               si = plink(ppnt(so)+lt)
               if(si.eq.sf) jtcmol = jtcmol +1
            enddo
         endif
      enddo

c ---------------------------------------------------- calculate <r>
      do si=1,natom
         ncoord(pmult(si)) = ncoord(pmult(si))+1
      enddo

      do j = 1,izcf
         totatm = totatm + ncoord(j)
         totcoord = totcoord + j*ncoord(j)
      enddo

      aveR = (1.0*totcoord)/(1.0*totatm)
      xpep = (1.0*jtcmol)/(1.0*totatm)
      avgR = (aveR+10.0*xpep)/(1.0+4.0*xpep)
c
c ------------------------------------------------ 2007:06:27	SN
c 7000 format (i5,1x,i5,1x,i5,1x,i5,1x,f9.5,1x,f9.5)
c      write (*,7000) ihb,totcoord,jtcmol,totatm,
c     *               hb_energy(hb_pick(ihb)),avgR
c
c	proflex5 is including a few more atoms and a few more bonds, 
c	at least at the later steps of HBdilution, than proflex4 did, 
c	but this does not change the fundamental results of either 
c	HBdilute or flexibility index values, so we did not track 
c	them down. This behavior was noticed when proflex 4 and 5 were 
c	run on 1bck (no HOH) with an E-cutoff of -2.3 kcal/mol and with
c	bond filter choices being: 5 for proflex5 and options 2,3,4,6,7,
c	and 15 for proflex4.
c
c ======================================================
c AJR 10.05.01 changes made to run the pebblegame twice 
c ----------------------------------- run pebblegame on skelton network
c     and store original topology into org* arrays
      orgsites = ns
      do so = 1,ns
         orgmult(so) = multref(so)
         orgpoint(so) = point(so)
         do k = 1,multref(so)
            orglink(orgpoint(so) +k) = linkref(point(so)+k)
         enddo
      enddo
      do so = 1,natom
         multref(so) = pmult(so)
         point(so) = ppnt(so)
         do k =1,orgmult(so)
c  must run the entire range to make sure that empty spots get mapped too
c         do k = 1,pmult(so)
            linkref(point(so) +k) = plink(ppnt(so)+k)
         enddo
      enddo
      ns = natom
      nclst = 0
      call mapmolecule()    
c      if(ihb .eq. 15) then
c         do si=1,natom
c            write(45,*) si,(linkref(point(si)+i),i=1,multref(si))
c         enddo
c         stop
c      endif
c         write(71,*) ihb,n_one,n_two,nbody
c         write(72,*) totatm,(ncoord(j),j=1,izcf)
c --------------- calculate lgst_clst   
c ---------------------------------------------------- temporarily change clst()
      k = 0
      do so=1,natom
         if( pmult(so) .eq. 1 ) then
            sf = plink( ppnt(so) + 1 )
            clst(so) = clst(sf)
         endif
         if( clst(so) .eq. so ) then
            k = k + 1
            newlabel( clst(so) ) = k
         endif
      enddo

      
      if( k .ne. nclst ) then
         write(16,*) k,nclst
c         stop
      endif

      do so=1,natom
         clst(so) = newlabel( clst(so) )
      enddo
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ----------------------------- count # of BULK & SURFACE sites for each cluster
c      write(6,*) nclst,iflop
      do nc=1,nclst
         size(nc) = 0
         multarc(nc) = 0
         multbrc(nc) = 0
      enddo
      do so=1,natom
         ipo = ppnt(so)
         nco = clst(so)
c ---------------------------------------------- count BULK sites in cluster nco
         multarc(nco) = multarc(nco) + 1
cc                           must cover all of the space (gaps created by pruning
cc               do jo=1,maxmult(so)
         if(pmult(so) .gt. 0) then
         do jo =1,pmult(so)
c         do jo=1,multref(so)
            indexo = ipo + jo
            sf = plink(indexo)
c            sf = linkref(indexo)
c ------------------------------------------------- prevent double counting bond
            if( sf .gt. so ) then
               ncf = clst(sf)
               if( nco .ne. ncf ) then
c ----------------------------------------------------------- count SURFACE site
                  size(nco) = size(nco) + 1
                  size(ncf) = size(ncf) + 1
               else
c --------------------------------------------------------- count BULK BOND-site
                  multbrc(nco) = multbrc(nco) + 1
               endif
            endif
         enddo
      endif
      enddo
c --------------------------------------------- calculate number of hinge joints
c                  2*nhinge= sum_{nc=1}^{nclst} [# of SURFACE sites in clst(nc)]
      nhinge = 0
      maxsize = -1
      do nc=1,nclst
         nhinge = nhinge + size(nc)
c -------------------------------------------------- size = SURFACE + BULK sites
         size(nc) = size(nc) + multarc(nc)
c ----------------------------------------------------------- find maximum size
         if( maxsize .lt. size(nc) ) maxsize = size(nc)
      enddo
      nhinge = nhinge/2
      lg_clst = maxsize
      num_clst = nclst
      xf = float(iflop)/(3.0*float(totatm))
      xxf = xf/(1.0+4.0*xpep)
      xclst = float(lg_clst)/float(totatm)
c      write(24,*) avgR,xxf,xclst,hb_energy(hb_pick(ihb))
c      open(24,file=fname,access='append',status='unknown')
      if( mode .eq. 1 ) then 
         open(24,file='pruned.dat',access='append',status='unknown')
         write(24,7588) avgR,xxf,xclst,hb_energy(hb_pick(ihb))
         close(24)
      endif
c      write(83,*) avgR,jtcmol,totatm,ypep
c      write(14,*) avgR,iflop,xf,xxf
c      write(15,*) avgR,lg_clst,xclst,num_clst
c ---------------------------------- restore original topology
      ns = orgsites
      do so = 1,ns
         multref(so) = orgmult(so)
         point(so) = orgpoint(so)
         do k = 1,orgmult(so)
            linkref(point(so) +k) = orglink(orgpoint(so)+k)
         enddo
      enddo
      return 
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine pclazz(nf,n)
      external neighbor
      integer n,nf(n)
      
      nf(1) = 1
      do jj =2,n
         nf(jj) = jj
         do kk=1,jj-1
            nf(kk) = nf(nf(kk))            
            if( neighbor(jj,kk).eq. 1) nf(nf(nf(kk)))=jj
c            if( neighb(jj,kk).eq. 1) nf(nf(nf(kk)))=jj
         enddo
      enddo
      do jj =1,n
         nf(jj)=nf(nf(jj))
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function neighbor(s1,s2)
      include      'set_parameter'
      integer      s1,s2
      integer      plink(nb2max),point(maxatm)
      integer      pmult(maxatm),maxmult(maxatm),gflag(maxgrp)
      common/ptopo/    pmult,plink,point,maxmult,gflag
      
      neighbor = 0
      do i = 1,pmult(s1)
         if(s2 .eq. plink(point(s1) +i)) neighbor = 1
      enddo
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hbpruneit(ns)

      include      'set_parameter'
      integer      s,so,msf,xsf,lxsf,sf,sf0,xso
      integer      plink(nb2max),point(maxatm)
      integer      pmult(maxatm),maxmult(maxatm),gflag(maxgrp)
      real         avgR
      common/ptopo/    pmult,plink,point,maxmult,gflag
c -------------------------------------------------- prune off deadends
      nATM = ns
      do s = 1,nATM
         so = s
   10    nbond = 0
         do icount =1,maxmult(so)
            if(plink(point(so)+icount) .gt. 0) nbond = nbond+1
         enddo
         if( nbond .eq. 1) then
            do mso = 1,maxmult(so)
               sf = plink(point(so) +mso)
               if(sf .ne. -1) xso = mso
            enddo
            sf = plink(point(so)+xso)
            plink(point(so)+xso) = -1
            do msf=1,maxmult(sf)
               if(plink(point(sf)+msf) .eq. so) xsf = msf
            enddo
            plink(point(sf)+xsf) = -1
            pmult(so) = pmult(so) -1
            pmult(sf) = pmult(sf) -1
            
            mbond = 0
            do jcount=1,maxmult(sf)
               if(plink(point(sf)+jcount) .gt. 0) mbond = mbond+1
            enddo
            if(mbond .eq. 1) then
               so = sf
               goto 10
            endif
         endif
      enddo
      ns = nATM
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hbcondenseMult(ns)
      include      'set_parameter'
      integer      plink(nb2max),point(maxatm),so,sf,sfnew
      integer      maxmult(maxatm),gflag(maxgrp),pmult(maxatm)
      integer      s, smax, sp, sfhold, holding(izcf)
      common/ptopo/ pmult,plink,point,maxmult,gflag

      do s=1,ns
         if(pmult(s) .gt. 0) then
            sp = point(s)
            do k =1,maxmult(s)
               holding(k) = 0
            enddo
            kcount =0
            do i =1,maxmult(s)
               sf = plink(sp+i)
               if(sf .gt. 0) then
                  kcount = kcount+1
                  holding(kcount) = sf
               endif
            enddo
            jcount =0
            do i=1,pmult(s)
               jcount = jcount+1
               plink(sp+i) = holding(jcount)
            enddo
            do i = pmult(s)+1,maxmult(s)
               plink(sp+i) = -1
            enddo
c  error check for hbcondenseMult
            nill = 0
            do j =1,pmult(s)
               so = plink(point(s)+j)
               if(so .eq. -1) nill = nill +1
            enddo
            if(nill .gt. 0) write(33,*) s,nill,maxmult(s),pmult(s)
         endif
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hbsiteprune(s1)

      include      'set_parameter'
      integer      s1,s2,x2,sf,x1
      integer      plink(nb2max),point(maxatm)
      integer      pmult(maxatm),maxmult(maxatm),gflag(maxgrp)
      real         avgR
      common/ptopo/    pmult,plink,point,maxmult,gflag


   10 do im = 1,maxmult(s1)
         if(plink(point(s1)+im) .ne. -1) x1 = im
      enddo
      s2 = plink(point(s1)+x1)
      do m2=1,maxmult(s2)
         if(plink(point(s2)+m2) .eq. s1) x2 = m2
      enddo
      plink(point(s2)+x2) = -1
      plink(point(s1)+x1) = -1
      pmult(s1) = pmult(s1) -1
      pmult(s2) = pmult(s2) -1
      if(pmult(s2) .eq. 1) then 
         s1 = s2 
         goto 10
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hbprunesc(aname,group)
      include      'set_parameter'
      integer      s,msf,xsf,lxsf,sf,g,gf
      integer      plink(nb2max),point(maxatm),group(maxatm)
      integer      slow(maxgrp),sbig(maxgrp),gflag(maxgrp),gnum(maxgrp)
      integer      pmult(maxatm),maxmult(maxatm)
      real         avgR
      character*4  aname(maxatm),atm,ats
      character*3  gname(maxgrp)
      character*5  idres(maxgrp)
      character*1  chainid(maxgrp)
      common/ptopo/    pmult,plink,point,maxmult,gflag
      common/groups/   ngrp,nres,nhet,nwater,chainid,idres,
     &                 gname,gnum,slow,sbig
c----------------------------------------- prune off non-Hbonded side chains
      do g = 1,ngrp
         if(gflag(g) .eq. -1) then
            do s=slow(g),sbig(g)
               ats = aname(s)
               if((ats .ne. ' N  ') .and.(ats .ne. ' O  ') .and. 
     &            (ats .ne. ' C  ') .and.(ats .ne. ' CA ') .and.
     &            (ats .ne. 'OXT')) then
                  do ix = 1,maxmult(s)
                     sf = plink(point(s)+ix)
                     if(sf.ne. -1) then
                        atm = aname(sf)
                        gf = group(sf)
                        if(gf .eq. g) then
                           if((atm .ne. ' N  ') .and.(atm .ne. ' O  ') 
     &                  .and. (atm .ne. ' C  ') .and.(atm .ne. ' CA ') 
     &                  .and. (atm .ne. 'OXT')) then
                              pmult(s) = pmult(s) -1
                              plink(point(s)+ix) = -1
                              do ixf=1,maxmult(sf)
                                 if(plink(point(sf)+ixf) .eq. s) then
                                    pmult(sf) = pmult(sf) -1
                                    plink(point(sf)+ixf) = -1
                                 endif
                              enddo
                           endif
                        endif
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hbringbreak(g,aname)
      include      'set_parameter'
      integer      g,gfval,sg0,sgX,s,sbk,srm,xsd,site
      integer      plink(nb2max),point(maxatm)
      integer      slow(maxgrp),sbig(maxgrp),gflag(maxgrp),gnum(maxgrp)
      integer      pmult(maxatm),mult(maxatm)
      real         avgR
      character*4  aname(maxatm),atm
      character*3  gname(maxgrp)
      character*5  idres(maxgrp)
      character*1  chainid(maxgrp)
      common/ptopo/    pmult,plink,point,mult,gflag
      common/groups/   ngrp,nres,nhet,nwater,chainid,idres,
     &                 gname,gnum,slow,sbig

      jcount = 0
   10 do s = slow(g),sbig(g)
         if(pmult(s).ne.0) sbk = s
      enddo
      mbk = pmult(sbk)
      do isk = 1,mbk
         srm = plink(point(sbk)+isk)
         do msd = 1,pmult(srm)
            if(plink(point(srm)+msd) .eq. sbk) xsd = msd
         enddo
         pmult(sbk) = pmult(sbk) -1
         pmult(srm) = pmult(srm) -1
         plink(point(sbk)+isk) = -1
         plink(point(srm)+xsd) = -1
         if(pmult(srm) .eq. 1) call hbsiteprune(srm)
      enddo
      jcount=jcount+1
c ------------------------------------------------ breaks 2nd ring in TRP      
      if((jcount .lt. 2).and.(gname(g) .eq. 'TRP')) goto 10
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine breakthis(ihb)
      include      'set_parameter'
      integer      hbH,hbA,xH,xA,ns,hb_id(maxh),hb_pick(maxh)
      integer      hbond(3,maxh)
      integer      mult(maxatm),pointer(maxatm),link(nb2max)
      dimension    hb_energy(maxh)
      character*1  hb_type(maxh)
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/topology/ ns,link,pointer,mult
      
      jhb = hb_pick(ihb)
      hbH = hbond(2,jhb)
      hbA = hbond(3,jhb)
      
C ERROR check:
      if( (mult(hbH).lt. 1). OR.(mult(hbA).lt. 1)) then
         write(56,*) 'Error: low multiplicity in hbond to break.'
         stop
      endif
      do iH = 1,mult(hbH)
         if( link(pointer(hbH)+iH) .eq. hbA) xH = iH
      enddo
      link(pointer(hbH)+xH) = -1
      if( xH .ne. mult(hbH)) call shrinkM(hbH)
      mult(hbH) = mult(hbH) -1
         
      do iA = 1,mult(hbA)
         if( link(pointer(hbA)+iA) .eq. hbH) xA = iA
      enddo
      link(pointer(hbA)+xA) = -1
      if( xA .ne. mult(hbA)) call shrinkM(hbA)
      mult(hbA) = mult(hbA) -1

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine shrinkM(s)

      include 'set_parameter'
      integer s,mult(maxatm),sp,s2,hold(izcf)
      integer point(maxatm),link(nb2max),ns
      common/topology/ ns,link,point,mult

      sp = point(s)
      kcount =0
      do k =1,mult(s)
         s2 = link(sp+k)
         if(s2 .gt. 0) then
            kcount = kcount+1
            hold(kcount) = s2
         endif
      enddo
      do j =1,kcount
         link(sp+j) = hold(j)
      enddo
      do i=kcount+1,mult(s)
         link(sp+i) = -1
      enddo

      return
      end


