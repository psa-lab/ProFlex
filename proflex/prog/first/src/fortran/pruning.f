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

      subroutine pruning()
c ------------------------------------------------------------------------------
c                                             PROGRAM WRITTEN:  October 5, 1999
c                                              Last Modified :  Sept.  28, 2001
c                                       program written by:     A. J. Rader
c                                                               rader@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     This subroutine is a modification of the standard output file, outputfirst
c which was written by Don Jacobs.  The goal here is to have a place where during
c the running of FIRST the user can have an option of running the program a 
c multiple number of times where the current process is somehow related to the 
c previous run.  This is a more interactive or iterative version of FIRST that 
c due to its nature would be overly cumbersome to generate all of the standard 
c FIRST output.  However, like the standard outputfirst file this program 
c calculates the cluster statistics for rigid clusters, Laman subgraphs, and 
c collective motions.  It does not give the typical FIRST summary or RasMol 
c script files.  Specifically this program prunes the protein down to a skeleton
c of atoms that are all 2 or higher-fold coordinated.  Then this skeleton, or
c pruned network is passed into the pebble game for analysis.  This process is 
c completed iteratively by removing hydrogen bonds one at a time, from weakest
c to strongest.
c ------------------------------------------------------------------------------
c INPUT: 
c        common/select/   nsc,criteria
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
c ==============================================================================
      include      'set_parameter'
      integer      si,sf0,totcoord,numatom,xmsf,msf,lmsf,totatm
      integer      hbA,hbH,xA,xH,site,gD,gA
      integer      halfcoord,uF,ncoord(izcf)
      integer      link_hinge(nb2max),label_hinge(nbmax)                  
      integer      g,h,hb,r,s,sf,so,tc,tcmol,th_max,th_min
      integer      gnum(maxgrp),slow(maxgrp),sbig(maxgrp)
      integer      linkref(nb2max),point(maxatm),locgroup(maxatm)
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      integer      torsion(2,nbmax)
      integer      newlabel(maxatm),oldlabel(maxatm)
      integer      redundant(maxatm)
      integer      size(maxatm),multarc(maxatm)
      integer      clst(maxatm),stress(-1:maxatm)
c      integer      multhld(maxatm),pointhld(maxatm),linkhld(nb2max)
      integer      pointarc(maxatm),pointbrc(maxatm),multbrc(maxatm)
      integer      mult(maxr),pointr(maxr)
      integer      cslow(0:maxcolor),csbig(0:maxcolor)
      integer      calow(maxcolor),cabig(maxcolor)
      integer      multref(maxatm),color(maxatm),paint(maxcolor)
      integer      maxmult(maxatm),lg_clst,num_clst
      integer      gflag(maxgrp),gp,ps1,ps2,co1
      integer      mbod(maxatm),tbod(-1:maxatm)
      real         avgR,xclst
      integer      nf(natom)
c      logical      neighbor
c      integer      neighbor,nf(natom)
c      dimension    nf(natom)
      dimension    sxyz(3,maxatm),freq(maxatm),bval(maxatm)
      dimension    hb_energy(maxh)
      character*80 input_data,outfile(12),criteria(99)
      character*80 fscript,fgrphcs,fdecomp,line
c      character*80 fmaxw,favgr,fisop,fmeane
      character*2  fkase

      character*80 frcrd
      character*11 fcoo
c      character*9  fmax,favg,fmean
c      character*8  fiso
      character*13 cname(maxcolor)
      character*5  idres(maxgrp)
      character*4  aname(maxatm),file_id,atm,atA,atD
      character*3  gname(maxgrp)
      character*1  rec(maxatm),chainid(maxgrp),hb_type(maxh)
      character*1  digit(0:9)
c      character*1  archive,answer0

      common/fnames/   n_root,input_data,outfile,file_id,digit


      common/select/   nsc,criteria
      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval
      common/groups/   ngrp,nres,nhet,nwater,chainid,idres,
     &                 gname,gnum,slow,sbig
      common/numbers1/ nflop,n_one,n_two,nbody,nclst,nstress
      common/numbers2/ iflop,nhinge,nih,ncmode
      common/topology/ ns,linkref,point,multref
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/subgraph/ stress
      common/cmotions/ label_hinge,link_hinge
      common/clusters/ clst
      common/prune/    maxmult

c ------------------------------------------------------------ format statements

 5000 format(5x,'ERROR!  number of floppy modes calculated ',
     &          'by sum rule = ',i7)
 5005 format(5x,'ERROR!  number of stressed regions calculated ',
     &          'by summing cluster numbers = ',i7)
 5010 format(5x,'ERROR!  number of independent hinges calculated ',
     &          'by summing cluster numbers = ',i7)
 5015 format(5x,'ERROR!  number of collective modes calculated ',
     &          'by summing cluster numbers = ',i7)
 5020 format(5x,'ERROR!  number of independent DOF calculated ',
     &          'by summing cluster numbers = ',i7)
c 6200 format(a4,1x,a3,1x,a1,i6,15x,a4,1x,a3,1x,a1,i6,3x,e14.7)
 6255 format(14x,i9,10x,i9)
c 6355 format(10x,i9,10x,i9,7x,i7)
c 6360 format(3i10)
 6400 format(2i10,5x,f9.5)
c 6450 format(5x,'Please wait:  Writing grphcs, decomp & ',
c     &          'RasMol script files')
 5500 format(5x,'Please select from the following options: ')
 5501 format(15x,'(1)   No pruning')
 5502 format(15x,'(2)   Standard Pruning (keep side chains)')
 5503 format(15x,'(3)   Full Pruning (remove dead-end side chains).')
 5505 format(i4)
 5600 format(6x,'Sorry that option is currently not working. ',
     &           'option (2) has been chosen instead.')
c 6455 format(5x,'These results are not being saved!')
c 6500 format(10x,' Archive decomposition?   Enter Y/N  ',
c     &             '--> Y is default') 
 6505 format(a1)
 7000 format(2x,11(i7))
c 7001 format(2x,10(i7))
 7002 format(f8.5,2(2x,i6,2x,i8),3x,f8.5)
 7003 format(f8.5,2(2x,i8,2x,i6))
 7100 format(f8.5,3i8,2f9.5)
 7183 format(f8.5,2i8,2f9.5)
 7200 format(9i8)
 7383 format(f8.5,8i8)
 7483 format(f8.5,4(2x,f9.5),3x,2i8)
 7583 format(f8.5,4(2x,f9.5))
 7585 format(f8.5,2x,f9.5,f11.5)
 7588 format(f9.6,2(2x,f11.9),2x,f9.5)
 7589 format(9x,2x,11x,2x,11x,2x,9x)
      call system( 'clear' )
      write(6,*)
      write(6,5501)
      write(6,5502)
      write(6,5503)
      read(5,5505) kopt
      if(kopt .eq. 1) fkase ='c1'
      if(kopt .eq. 2) fkase ='c2'
      if(kopt .eq. 3) fkase ='c3'
c---------------------------------------------------------------------------c
c                           Iterative removal of H-bonds.                   c
c                                                                           c
c This is the standard H-bond dilution scheme where hydrogen bonds are      c
c removed one at a time starting with the lowest energy.                    c
c                                                                           c
c===========================================================================c

      khb0 = khb
      call place_Hbond()
c ----------------------------------- setup to not remove hydrophobics 
        nphob=0
        khbond=0
        nsbdg=0

        do i = 1,khb0
	if(hb_type(hb_pick(i)) .eq. 'B') nphob= nphob+1
        if(hb_type(hb_pick(i)) .eq. 'I') nsbdg=nsbdg+1
        if(hb_type(hb_pick(i)) .eq. 'H') khbond=khbond+1
	enddo
c        write(99,*) khb0,nphob,khbond,nsbdg

c      n_end = n_root+11
c      fmax =  '_maxwell.'
c      favg =  '_average.'
c      fiso =  '_isopep.'
c      fmean = '_meanNRG.'
      fcoo = '_meancoord.'
c      fmaxw = input_data(1:n_root)//fmax//fkase
c      favgr = input_data(1:n_root)//favg//fkase
c      fisop = input_data(1:n_root)//fiso//fkase
c      fmeane = input_data(1:n_root)//fmean//fkase
      frcrd = input_data(1:n_root)//fcoo//fkase

c      open(80,file=fmaxw,status='UNKNOWN',access='APPEND')
cc      open(80,file="maxwell.dat",status='UNKNOWN', access='APPEND')
c      open(81,file=favgr,status='UNKNOWN',access='APPEND')
c      open(83,file=fisop,status='UNKNOWN',access='APPEND')
c      open(85,file=fmeane,status='UNKNOWN',access='APPEND')
      open(88,file=frcrd,status='UNKNOWN',access='APPEND')
c	  write(88,7589) '<r>',' ','F_Floppy_Modes',' ','F_Atms_Lgst_Clst',
	  write(88,*) '<r>',' ','F_Floppy_Modes',' ','F_Atms_Lgst_Clst',
     &' ','HB_Enrgy'


c--------------------------------------- make neighbor table with current hbonds
         do so = 1,natom
         maxmult(so) = multref(so)
         mbod(so) = 0
         enddo

c ----------------------------------------- break 1 HB at a time,weakest first
         kindex = khb0-nphob
         do ihb = 0,kindex
c         write(92,*) hb_pick(ihb)
         write(6,*) ihb,kindex,khb0,nphob
c ------------------------------------------------- initialize
            do i=1,izcf
            ncoord(i) = 0
            enddo
         totcoord =0
         totatm =0

         if( ihb .eq. 0) then             
c -------------------------------------------------- initial pruning
            if( kopt .eq. 2) then
               call pruneit()
               call condenseMult()
            elseif(kopt .eq. 3) then
               do ig = 1,ngrp
                  gflag(ig) = 0
               enddo
               do jhb = 1,khb0
                  lhb = hb_pick(jhb)
                   do i=1,3,2
                     g = locgroup(hbond(i,lhb)) 
         if(( gname(g) .eq. 'PHE').or. (gname(g) .eq. 'TRP') .or. 
     &      ( gname(g) .eq. 'HIS').or. (gname(g) .eq. 'TYR')) then
            atm = aname(hbond(i,lhb))
            if((atm .ne. ' N  ') .and.(atm .ne. ' O  ') .and. 
     &         (atm .ne. ' C  ') .and.(atm .ne. ' CA ') .and. 
     &         (atm .ne. ' OXT')) gflag(g) = gflag(g)+1
c --------------------(added 05.02.01 ajr)-------------------------------------
c             to cope with hydrophobic tether pseudoatoms
c         NOTE: must change code if # of psuedoatoms increases
         elseif( gname(g) .eq. 'XXX') then
            ps2 = hbond(1,lhb)
            ps1 = linkref(point(ps2)+1)
            if(ps1 .eq. hbond(2,lhb)) ps1 = linkref(point(ps2)+2)
            co1 = linkref(point(ps1)+1)
            if(co1 .eq. ps2) co1 = linkref(point(ps1)+1)
            g = locgroup(co1)
            if(( gname(g) .eq. 'PHE').or. (gname(g) .eq. 'TRP') .or. 
     &        ( gname(g) .eq. 'HIS').or. (gname(g) .eq. 'TYR')) then
               atm = aname(hbond(i,lhb))
            if((atm .ne. ' N  ') .and.(atm .ne. ' O  ') .and. 
     &         (atm .ne. ' C  ') .and.(atm .ne. ' CA ') .and. 
     &         (atm .ne. ' OXT')) gflag(g) = gflag(g)+1
            endif
c -------------------------(added 05.02.01 ajr)-------------------------------------
         endif
                  enddo
               enddo
               do ig=1,ngrp
               if(((gname(ig) .eq. 'PHE').or.(gname(ig) .eq. 'TRP') .or. 
     &            ( gname(ig) .eq. 'HIS').or.(gname(ig) .eq. 'TYR'))
     &            .and. (gflag(ig).eq.0) ) gflag(ig) = -1
               enddo
               call prunesc(aname,gflag,locgroup)
               call pruneit()
               call condenseMult()
            endif

c --------------------- check on number of torsional constraints
            jtcmol = 0
            do tc = 1,tcmol
               so = torsion(1,tc)
               sf = torsion(2,tc)
               if((multref(so).ne.0).and.(multref(sf).ne.0)) then
                  do lt =1,multref(so)
                     si = linkref(point(so)+lt)
                     if(si.eq.sf) jtcmol = jtcmol +1
                  enddo
               endif
            enddo
c ------------------------------------------------ calculate <r>
         do si=1,natom
         ncoord(multref(si)) = ncoord(multref(si))+1
         enddo
         do j = 1,izcf
         totatm= totatm +ncoord(j)
         totcoord = totcoord + j*ncoord(j)
         enddo

         avgR = float(totcoord)/(float(totatm))

c -------------------------------------- macromolecule decomposition 
c         do ilip =1,natom
c            write(98,*) ilip,point(ilip),multref(ilip)
c         enddo
         call mapmolecule()
c         do ilip =1,natom
c            write(99,*) ilip,point(ilip),multref(ilip)
c        enddo
         else
c ----------------------------------------- break 1 HB at a time,weakest first
         lhb = hb_pick(ihb)
         hbH = hbond(2,lhb)
         hbA = hbond(3,lhb)
         if((multref(hbH) .gt. 0).AND.(multref(hbA).gt. 0)) then
            do iH = 1,multref(hbH)
            if(linkref(point(hbH)+iH) .eq. hbA) xH = iH
            enddo
         linkref(point(hbH)+xH) = -1
         multref(hbH) = multref(hbH) -1
         
         do iA = 1,multref(hbA)
            if( linkref(point(hbA)+iA) .eq. hbH) xA = iA
         enddo
         linkref(point(hbA)+xA) = -1
         multref(hbA) = multref(hbA) -1
         call condenseMult()
         if( kopt .ne. 1) then 
         if( multref(hbH) .eq. 1) call siteprune(hbH)
         if( multref(hbA) .eq. 1) call siteprune(hbA)
         call condenseMult()
         if(kopt .eq. 3) then
            hbD = hbond(1,lhb)
            gA = locgroup(hbA)
            atA = aname(hbA)
            if( gflag(gA) .gt. 0) then
               if((atA.eq.' ND1').or.(atA.eq.' NE1').or.(atA.eq.' NE2')
     &            .or.(atA.eq.' OH ')) gflag(gA) = gflag(gA) -1
               if(gflag(gA) .eq. 0) call ringbreak(gA,aname)
            endif
            gD = locgroup(hbD)
            atD = aname(hbD)
            if( gflag(gD) .gt. 0) then
               if((atD.eq.' ND1').or.(atD.eq.' NE1').or.(atD.eq.' NE2')
     &            .or.(atD.eq.' OH ')) then
               gflag(gD) = gflag(gD) -1
               if(gflag(gD) .eq. 0) call ringbreak(gD,aname)
               endif
            endif
         call condenseMult()
         endif
         endif
         endif

c --------------------------------- identify the connected backbone
         do si = 1,natom
            nf(si) = 0
         enddo
         call eclazz(nf,natom)
c -------------------------------------- calculate number of bodies
         do si = 1,natom
            mbod(si) = 0
         enddo
c                        that was to initialize the body label, mbod(s)
         do si = 1,natom
            mbod(nf(si)) = mbod(nf(si)) +1
         enddo
c         if(ihb.eq. 1) then 
c            do i = 1,natom 
c               write(10,*) i,nf(i),mbod(i)
c            enddo
c         endif
         ibod=0
         jmbod = 0
c ibod = number of bodies, jmbod = number of nontrivial bodies (i.e more than 1 atom)
         maxbod = 0
         do si = 1,natom
            if(mbod(si) .gt. 0) then
               if(mbod(si) .gt. maxbod) maxbod = mbod(si)
c               jbod(mbod(si)) = jbod(mbod(si))+1
               ibod = ibod+1
               if((mbod(si) .gt. 1).AND.( mbod(si) .lt. maxbod) ) then 
                  jmbod = jmbod +1
               endif
cc               write(11,*) si,mbod(si)
ccc               write(11,*) si,mbod(si),jbod(mbod(si))
            endif
         enddo
c         write(65,*) ihb,ibod,maxbod
cc         stop

c ------ removing smaller isolated bodies
c -- first identify the bodies and their sizes. mbod(nf(s)) = size(that body).
         if(kopt .ne.1) then
            if (ibod .gt. 1) then
               do so = 1,natom 
                  if(mbod(nf(so)) .eq. 1) then
c     if(mbod(so) .ne. maxbod) then
c     if(mbod(so) .gt. 1) then
                     do i= 1,multref(so)
                        linkref(point(so) +i) = -1
                     enddo
                     multref(so) = 0
                     ibod = ibod -1
c     AJR 03.25.02   now remove isolated rings if not part of the largest body (maxbod)
                  elseif(mbod(nf(so)).ne. maxbod) then     
                     jmbod = jmbod +1
                  endif
cc     write(88,*) so, mbod(so),nf(so),mbod(nf(so)),multref(so)
cc             c  remove rings
cc             ieliminate = mbod(nf(so))
cc                      do i=1,multref(so)
ccc                   sf = linkref(point(so)+i)
ccc                   if(multref(sf) .ge. 2) then                      
cc                         linkref(point(so)+i) = -1
cc                      enddo
c                      do i=1,multref(so)
c                         linkref(point(so)+i) = -1
c                      enddo
c                      multref(so) =0
c                      ibod=ibod-1
c                  endif
               enddo               

c     AJR 03.25.02   now remove isolated rings if not part of the largest body (maxbod)
            endif
            if(ibod .gt. 1) then
c              write(65,*) ihb,ibod,maxbod,jmbod
               do so = 1,natom 
                  lbody = mbod(nf(so))
                  if((lbody .gt. 1) .AND.(lbody .ne. maxbod)) then
                     do i= 1,multref(so)
                        linkref(point(so) +i) = -1
                     enddo
                     multref(so) = 0
                     jmbod = jmbod -1
                  endif
               enddo
            endif
c            write(65,*) ihb,ibod,maxbod,jmbod
            if(jmbod .eq. 0) ibod = 1

c -----------------------------------------error check
c            if(ibod .gt. 1) write(17,2100) ibod
            if(ibod .gt. 1) write(17,2100) ibod,lhb
 2100 format(5x,'Abandon all hope, the number of bodies ', i4,   
     &       ' is gt 1. ', i5) 
c     &       ' is gt 1.') 
            if(ibod .gt. 2) then
               do s =1,natom
                  write(66,*) s,mbod(s),nf(s),mbod(nf(s)),multref(s)
               enddo
            endif
         endif
c         write(56,*) ihb,ibod,maxbod
c         if(ihb.eq. 103) then 
c            do i = 1,natom 
c               write(20,*) i,multref(i),nf(i)
c            enddo
c         endif
c -------------------------------------- macromolecule decomposition 
c         write(70,*) ihb,n_one,n_two,nbody
c         write(6,*) natom,ns,totatm
c      if(ihb .eq. 15) then
c$$$      if(lhb .eq. 23) then
c$$$         do si=1,natom
c$$$cc            write(99,*) si,point(si),multref(si)
c$$$            write(95,*) si,(linkref(point(si)+i),i=1,multref(si))
c$$$         enddo
c$$$      elseif(lhb.eq. 24) then
c$$$         do si=1,natom
c$$$            write(96,*) si,(linkref(point(si)+i),i=1,multref(si))
c$$$         enddo
c$$$         stop
c$$$      endif         

         call mapmolecule()
c         write(71,*) ihb,n_one,n_two,nbody
         
c         if(ihb.eq. 103) then 
c            do i = 1,natom 
c               write(21,*) i,multref(i),nf(i)
c            enddo
c         endif
c ------------------------------------------------- initialize
         do i=1,izcf
            ncoord(i) = 0
         enddo
         totcoord =0
         totatm =0
         sum1 = 0.0
c --------------------- check on number of torsional constraints
            jtcmol = 0
            do tc = 1,tcmol
               so = torsion(1,tc)
               sf = torsion(2,tc)
               if((multref(so).ne.0).and.(multref(sf).ne.0)) then
                  do lt =1,multref(so)
                     si = linkref(point(so)+lt)
                     if(si.eq.sf) jtcmol = jtcmol +1
                  enddo
               endif
            enddo
c --------------------------------------------------- calculate <r>

         do si=1,natom
            ncoord(multref(si)) = ncoord(multref(si))+1
         enddo
         if(kopt.ne.1) then
            if(ncoord(1) .gt. 0) write(98,*) ncoord(1)
         endif
         do j = 1,izcf
            totatm= totatm +ncoord(j)
            totcoord = totcoord + j*ncoord(j)
            sum1 = sum1+2.5*j*ncoord(j)
         enddo
         if(kopt .eq. 1) then
            if(totatm .ne. natom - n_one) then
            write(15,*) ihb,' ERROR: natom mixup. ',totatm,natom - n_one
            do j=1,izcf
	       write(15,*) j, ncoord(j)
            enddo
            write(15,*)
            do k =1,natom
               if( multref(k) .gt. 8) write(15,*) k, multref(k)
	    enddo
            stop
            endif
         endif
         avgR = float(totcoord)/(float(totatm))
      endif
      n1 = ncoord(1)
      isum1 = nint(sum1)
c      write(75,*) aveR,nclst,iflop,nflop
c      write(72,*) totatm,(ncoord(j),j=1,izcf)
c -------------------------------------------------------- write to output files
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++ temporary code: REMOVE
c ---------------------------------------------------- temporarily change clst()
      k = 0
      ierr =0
         do so=1,natom
            if( multref(so) .eq. 1 ) then
               if(kopt .ne. 1) then 
               write(12,*) ihb,' Error: mult=1',so
               ierr = ierr +1
c               stop
               else
               sf = linkref( point(so) + 1 )
               clst(so) = clst(sf)
               endif
            endif
            if( clst(so) .eq. so ) then
            k = k + 1
            newlabel( clst(so) ) = k
            endif
         enddo      
            if(ierr .ne. 0) stop

            if( k .ne. nclst ) then
            write(16,*) k,nclst
            stop
            endif
      
            do so=1,natom
            clst(so) = newlabel( clst(so) )
            enddo
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
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
c                           must cover all of the space (gaps created by pruning
c               do jo=1,maxmult(so)
            if(multref(so) .gt. 0) then
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
c            write(55,*) ihb,lg_clst,num_clst
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
cc --------------------------------------------------- initialize histogram bins
cc                                           use oldlabel() for binning the data
c            do s=1,maxsize
c            oldlabel(s) = 0
c            enddo
cc ------------------------------------------- calculate cluster size histogram 
c         do nc=1,nclst
c         oldlabel( size(nc) ) = oldlabel( size(nc) ) + 1
c         enddo
c
cc -------------------------- write to analysis file the cluster size statistics
c      write(60,*) 
c      write(60,6250) 
c         do s=1,maxsize
c         if( oldlabel(s) .ne. 0 ) write(60,6255) s,oldlabel(s)
c         enddo
cc ---------------------------------------------------------- GLOBAL ERROR CHECK
         kflop = 3*natom + nhinge
            do nc=1,nclst
            isize = size(nc)
               if( isize .gt. 2 ) then
               kflop = kflop - (3*isize - 6)
               elseif( isize .eq. 2 ) then
               kflop = kflop - 1
               endif
            enddo
c --------------------------------------------------------------- ERROR FOUND!!!
            if( kflop .ne. nflop ) then
               write(56,*) 
               write(56,*) 
               write(56,*) 
               write(56,5000) kflop
               stop
            endif
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ FINISHED ERROR CHECK

c --------------------------- relabel rigid clusters in descending order on size
c                                 use pointarc and pointbrc as dummy work arrays
c                                      set up pointers for range in cluster size
c                                                                 define nc2,nc1
      call clstrlabels(size,pointarc,nclst,pointbrc,maxsize,clst,
     &                 newlabel,natom,mult,pointr,maxr,nc2,nc1)
c     &                 newlabel,natom,mult,pointr,maxr,nc2,nc1)
c ---------------------------------------------------------------------- relabel
      do s=1,natom
      clst(s) = newlabel( clst(s) )
      enddo
c --------------------------- use pointarc,pointbrc,oldlabel as temporary arrays
         do nc=1,nclst
         kc = newlabel(nc)
         pointarc(kc) = multarc(nc)
         pointbrc(kc) = multbrc(nc)
         oldlabel(kc) = size(nc)
         enddo
c ------------------------------------ transfer back into the appropriate arrays
      do nc=1,nclst
      multarc(nc) = pointarc(nc)
      multbrc(nc) = pointbrc(nc)
      size(nc) = oldlabel(nc) 
      enddo
c --------------------------------------------------------- color rigid clusters
c                    use pointBrc as a chain for linking BULK sites in a cluster
c                                use pointArc to bin all BULK sites in a cluster
c                                                                 define cname()
c      call colorclustr(clst,pointbrc,natom,pointarc,color,nclst,
c     &                       mult,pointr,cname)
c ---------------------------------------------------------------- relabel sites 
c                                                     for quick graphics display
c            INPUT:  use pointBrc as a chain for linking BULK sites in a cluster
c            INPUT:              use pointArc to bin all BULK sites in a cluster
c             NOTE:   pointArc,pointBrc are set up from subroutine colorclustr()
c           OUTPUT: pointArc,pointBrc,oldlabel,newlabel,paint,ncolor,calow,cabig 
c                                                                and cslow,csbig

c      call relabelatms(nhinge,clst,pointbrc,oldlabel,newlabel,natom,
c     &    pointarc,color,nclst,paint,ncolor,calow,cabig,cslow,csbig)

cc --------------------------------------------- write corresponding output files
c         fdecomp = outfile(3)
c         fgrphcs = outfile(6)
c         fscript = outfile(7)
c         itemp = natom
c         write(6,6450)
c         call out_script(fscript,fgrphcs,pointarc,multarc,
c     &     pointbrc,multbrc,nclst,calow,cabig,cslow,csbig,
c     &     paint,ncolor,cname,mult,pointr,nc1,nc2)
c         call out_PDBfile(fgrphcs,newlabel,
c     &                    itemp,clst,pointbrc,multbrc,nclst)
c         call out_decomp(fdecomp,oldlabel,natom,clst,nclst)
c ---------------------------------------- calculate Nr (# redundant constraints)
c     commented out to see if I can get an accurate calc of Nr for all cases
c         if(kopt .eq. 1) then 
c         nf1 = (29*ncoord(7) + 19*ncoord(5) + 9*ncoord(3) +ncoord(1))/2
c         nf2 = 17*ncoord(8) +12*ncoord(6) + 7*ncoord(4) + 2*ncoord(2)
c         uF = 3*totatm - nf1 - nf2 - tcmol -6*nbody
c         nrsum = iflop - uF
cc            do so = 1,natom
cc            ncoord(multref(so)) = ncoord(multref(so))+1

c         else
c ------------------------------------------------ calculate bond weight factors
         do s=1,natom
         multarc(s) = 0
         multbrc(s) = 0
         size(s) = 0
         redundant(s) = 0
         enddo
         do so=1,natom
         mc = 0
         index = point(so)
         istress = stress(so)
         multarc(istress) = multarc(istress) + 1                  
c                           must cover all of the space (gaps created by pruning
c            do jo=1,maxmult(so)
            do jo=1,multref(so)
            index = index + 1
            sf = linkref(index)
            if(sf .ne. -1) then
c -------------------------------------------------------- allow double counting
               if( istress .eq. stress(sf) ) then 
               size(istress) = size(istress) + 1    
               mc = mc + 1
               endif
            endif
            enddo
c         multbrc(istress) = multbrc(istress) + mc*(mc-1)/2   
c ------------- change per DonJacobs
         multbrc(istress) = multbrc(istress) + 2*mc-3   
         enddo
c ---------------------------------------------- finalize CF and BB constraints
         do s=1,natom
         size(s) = size(s)/2
         multbrc(s) = multbrc(s) + size(s) 
         enddo
c ------------------------------------------- augment dihedral angle constraints
            do tc=1,tcmol
            so = torsion(1,tc)
            sf = torsion(2,tc)
            istress = stress(so)
               if( istress .eq. stress(sf) ) then
               multbrc(istress) = multbrc(istress) + 1
               endif
            enddo
         nr=0
         do so=1,natom-1
         index = point(so)
         istress = stress(so)
            do jo=1,multref(so)
            index = index + 1
            sf = linkref(index)
c ------------------------------------------------------ prevent double counting
               if( sf .gt. so ) then
                  if( istress .eq. stress(sf) ) then 
                  nr = multbrc(istress) - (3*multarc(istress) - 6) 
               if( redundant(istress) .lt. nr) redundant(istress) = nr
                  endif
               endif
            enddo
         enddo
         nrsum = 0
         do jstress = 1,natom
            nrsum = redundant(jstress) + nrsum
         enddo
c         endif
         xnr = float(nrsum)/(3.0*float(totatm))
         xf = float(iflop)/(3.0*float(totatm))
c         write(83,*) avgR,jtcmol,totatm,xpep
c         write(81,7100) avgR,iflop,nrsum,totatm,xf,xnr

c         write(86,7100) avgR,nrsum,iflop,nstress,xnr,xf
c ######## uncommenting the following lines will give 
c                 maxsize = size (# of sites) for largest stressed cluster size  
cc -------------------------------------------------- temporarily change stress()
c      k = 0
c         do so=1,natom
c            if( stress(so) .eq. so ) then
c            k = k + 1
c            newlabel(so) = k
c            endif 
c         enddo
cc ---------------------------------------------- use nclst and clst() variables
c      nclst = k
c         do so=1,natom
c         clst(so) = newlabel( stress(so) )
c         enddo
cc ---------------------------------------- count # of BULK sites for each region
c            do nc=1,nclst
c            multarc(nc) = 0
c            multbrc(nc) = 0
c            size(nc) = 0
c            enddo
c         nhinge = 0
c         do so=1,natom
c         ipo = point(so)
c         nco = clst(so)
cc ---------------------------------------------- count BULK sites in cluster nco
c         multarc(nco) = multarc(nco) + 1
c            do jo=1,multref(so)
c            indexo = ipo + jo
c            sf = linkref(indexo)
cc ------------------------------------------------- prevent double counting bond
c               if( sf .gt. so ) then
c               ncf = clst(sf) 
c                  if( nco .eq. ncf ) then
cc --------------------------------------------------------- count BULK BOND-site
c                  multbrc(nco) = multbrc(nco) + 1
c                  else
cc ----------------------------------------------------------- bond is unstressed
c                  nhinge = nhinge + 1
cc ----------------------------------------------------------- count SURFACE site 
c                  size(nco) = size(nco) + 1
c                  size(ncf) = size(ncf) + 1
c                  endif
c               endif
c            enddo
c         enddo
cc <<<<<<<<<<<<<<<<<< LAMAN SUBGRAPH analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cc ----------------------------------------- calculate number of unstressed bonds
c         maxsize = -1
cc ------------------------------------------------------------ find maximum size 
c            do nc=1,nclst
c            s = multarc(nc)
c            if( maxsize .lt. s ) maxsize = s
c            enddo
cc ---------------------------------------------------- initialize histogram bins
cc                                            use oldlabel() for binning the data
c            do s=1,maxsize
c            oldlabel(s) = 0
c            enddo
cc --------------------------------------------- calculate cluster size histogram 
c         do nc=1,nclst
c         s = multarc(nc)
c         size(nc) = size(nc) + s                       
c         oldlabel(s) = oldlabel(s) + 1
c         enddo
cc --------------------------- write to analysis file the cluster size statistics
c      write(60,*) 
c      write(60,6300) 
c      kstress = 0
c         do s=2,maxsize
c         kstress = kstress + oldlabel(s)
c         if( oldlabel(s) .ne. 0 ) write(60,6255) s,oldlabel(s)
c         enddo
c ^^^^^^^^^^^^^^^^^^^^^^^^ STOP uncommenting to ouput Lgst Stressed Clst size
cc ----------------------------------------------------------- GLOBAL ERROR CHECK
c             if( kstress .ne. nstress ) then
cc --------------------------------------------------------------- ERROR FOUND!!!
c             write(6,*) 
c             write(6,*) 
c             write(6,5005) kstress
c             stop
c             endif
cc -------------------------------------------------------- find NEW maximum size 
c         maxsize = -1
c            do nc=1,nclst
c            if( maxsize .lt. size(nc) ) maxsize = size(nc)
c            enddo
cc --------------------------- relabel rigid clusters in descending order on size
cc                                 use pointarc and pointbrc as dummy work arrays
cc                                      set up pointers for range in cluster size
cc                                                                 define nc2,nc1
c      call clstrlabels(size,pointarc,nclst,pointbrc,maxsize,clst,
c     &                 newlabel,natom,mult,pointr,maxr,nc2,nc1)
cc ---------------------------------------------------------------------- relabel
c      do s=1,natom
c      clst(s) = newlabel( clst(s) )
c      enddo
cc --------------------------- use pointarc,pointbrc,oldlabel as temporary arrays
c         do nc=1,nclst
c         kc = newlabel(nc)
c         pointarc(kc) = multarc(nc)
c         pointbrc(kc) = multbrc(nc)
c         oldlabel(kc) = size(nc)
c         enddo
cc ------------------------------------ transfer back into the appropriate arrays
c      do nc=1,nclst
c      multarc(nc) = pointarc(nc)
c      multbrc(nc) = pointbrc(nc)
c      size(nc) = oldlabel(nc) 
c      enddo
cc --------------------------------------------------------- color rigid clusters
c                    use pointBrc as a chain for linking BULK sites in a cluster
c                                use pointArc to bin all BULK sites in a cluster
c                                                                 define cname()
c      call colorclustr(clst,pointbrc,natom,pointarc,color,nclst,
c     &                       mult,pointr,cname)
cc ---------------------------------------------------------------- relabel sites 
c                                                     for quick graphics display
c            INPUT:  use pointBrc as a chain for linking BULK sites in a cluster
c            INPUT:              use pointArc to bin all BULK sites in a cluster
c             NOTE:   pointArc,pointBrc are set up from subroutine colorclustr()
c           OUTPUT: pointArc,pointBrc,oldlabel,newlabel,paint,ncolor,calow,cabig 
c                                                                and cslow,csbig

c      call relabelatms(nhinge,clst,pointbrc,oldlabel,newlabel,natom,
c     &    pointarc,color,nclst,paint,ncolor,calow,
c     &    cabig,cslow,csbig)
cc --------------------------------------------- write corresponding output files
c         fdecomp = outfile(4)
c         call out_decomp(fdecomp,oldlabel,natom,clst,nclst)
cc --------------------------- write to analysis file the cluster size statistics
c         do s=th_min,th_max
c         multbrc(s) = 0
c         size(s) = 0
c         enddo
c         do tc=th_min,th_max
c         s = label_hinge(tc)
c            if( s .gt. 0 ) then
c            size(s) = size(s) + 1
c            if( torsion(2,tc) .lt. 0 ) multbrc(s) = multbrc(s) + 1
c            endif
c         enddo
c         do nh=1,nhinge0
c         oldlabel(nh) = 0
c         enddo
c            do so=th_min,th_max
c            s = size(so) 
c            if( s .gt. 0 ) oldlabel(s) = oldlabel(s) + 1
c            enddo
cc ----------------------------------------------------------- GLOBAL ERROR CHECK
c                  k = 0
c                  do nh=1,nhinge0
c                  k = k + nh*oldlabel(nh)
c                  enddo
c            k = nhinge0 - k
c            if( k .ne. nih ) then
c            write(6,*) 
c            write(6,*) 
c            write(6,5010) k
c            stop
c            endif
c      write(60,*) 
c      write(60,6350) 
c      kcmode = 0
c      idof = nih
c         do s=th_min,th_max
c            if( size(s) .gt. 0 ) then
c            kcmode = kcmode + 1
c            idof = idof + multbrc(s)
c            write(60,6355) kcmode,size(s),multbrc(s)
c            endif
c         enddo
cc ----------------------------------------------------------- GLOBAL ERROR CHECK
c             if( kcmode .ne. ncmode ) then
cc --------------------------------------------------------------- ERROR FOUND!!!
c             write(6,*) 
c             write(6,*) 
c             write(6,5015) kcmode
c             stop
c             endif
c                if( idof .ne. iflop ) then
cc --------------------------------------------------------------- ERROR FOUND!!!
c                write(6,*) 
c                write(6,*) 
c                write(6,5020) idof
c                stop
c                endif
c -------------------------------------------------------- write to bond_wt file
c <<<<<<<<<<<<<<<<<<<<            WRITES Flexibility index to bond_wt file 
c <<<<<<<<<<<<<<<<<<<<                        for underconstrained bonds
c         do tc=th_min,th_max
c         so = torsion(1,tc)
c         sf = iabs( torsion(2,tc) )
c         s = label_hinge(tc)
c            if( s .gt. 0 ) then
c            x = float( multbrc(s) )/float( size(s) )
c            else
c            x = 1.0e0
c            endif
c         write(59,6400) so,sf,x
c         enddo
cc ---------------------------------------- temporarily make a fictitious clst()
c            do s=1,natom
c            clst(s) = 1
c            enddo
c         nclst = 1
c         do tc=th_min,th_max
c         newlabel(tc) = -1
c         enddo
c            do tc=th_min,th_max
c            s = label_hinge(tc)
c               if( s .gt. 0 ) then
c                  if( newlabel(s) .lt. 0 ) then
c                  nclst = nclst + 1
c                  newlabel(s) = nclst
c                  endif
c               so = torsion(1,tc)
c               sf = iabs( torsion(2,tc) )
c               clst(so) = newlabel(s)
c               clst(sf) = newlabel(s)
c               endif
c            enddo
c      fdecomp = outfile(2)
c      open(99,file=fdecomp,status='new')
c      rewind(99)
c         do tc=th_min,th_max
c         so = torsion(1,tc)
c         sf = iabs( torsion(2,tc) )
c         s = label_hinge(tc)
c            if( s .gt. 0 ) then
c            label = newlabel(s) - 1
c            else
c            label = -1
c            endif
c         write(99,6360) so,sf,label
c         enddo
c      close(99)
cc ---------------------------------------- count # of BULK sites for each region
c            do nc=1,nclst
c            multarc(nc) = 0
c            multbrc(nc) = 0
c            size(nc) = 0
c            enddo
c         nhinge = 0
c         do so=1,natom
c         ipo = point(so)
c         nco = clst(so)
cc ---------------------------------------------- count BULK sites in cluster nco
c         multarc(nco) = multarc(nco) + 1
c            do jo=1,multref(so)
c            indexo = ipo + jo
c            sf = linkref(indexo)
cc ------------------------------------------------- prevent double counting bond
c               if( sf .gt. so ) then
c               ncf = clst(sf) 
c                  if( nco .eq. ncf ) then
cc --------------------------------------------------------- count BULK BOND-site
c                  multbrc(nco) = multbrc(nco) + 1
c                  else
cc ----------------------------------------------------------- bond is unstressed
c                  nhinge = nhinge + 1
cc ----------------------------------------------------------- count SURFACE site 
c                  size(nco) = size(nco) + 1
c                  size(ncf) = size(ncf) + 1
c                  endif
c               endif
c            enddo
c         enddo
cc ----------------------------------------- calculate number of unstressed bonds
c         maxsize = -1
cc ------------------------------------------------------------ find maximum size 
c            do nc=1,nclst
c            s = multarc(nc)
c            if( maxsize .lt. s ) maxsize = s
c            enddo
c ---------------------------------------------------- initialize histogram bins
c                                            use oldlabel() for binning the data
c            do s=1,maxsize
c            oldlabel(s) = 0
c            enddo
c --------------------------------------------- calculate cluster size histogram 
c                                                                BULK sites only
c         do nc=1,nclst
c         s = multarc(nc)
c         size(nc) = size(nc) + s                           
c         oldlabel(s) = oldlabel(s) + 1
c         enddo
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
cc         write(82,7200) ihb,(ncoord(j),j=1,izcf)
cc         write(85,7002) avgR,iflop,totatm,ncmode,maxsize,xf
cc   above option is original prune15, below is option with number of bonds
c         halfcoord = totcoord/2
cc         jtest = iflop+6*nbody+jtcmol-nrsum+n1+isum1-6*totatm
         xpep = jtcmol/(3.0*totatm)
         xbod = 2.0*nbody/totatm
         xn1 = n1/(3.0*totatm)
c      write(80,7483) avgR,xf,xf+xpep+xbod+xn1,xpep,xbod,totatm,n1

      ypep = float(jtcmol)/float(totatm)
      aveR = (avgR+10.0*ypep)/(1.0+4.0*ypep)
      xxf = xf/(1.0+4.0*ypep)
c      xclst = float(lg_clst)/float(natom)
      xclst = float(lg_clst)/float(totatm)
      if(kopt .eq. 1) then
         xn1 = 0.4*float(n1)/float(totatm)
         aveR = aveR+xn1
      endif
c      write(83,*) aveR,jtcmol,totatm,ypep
c      write(84,*) aveR,iflop,xf,xxf
c      write(85,*) aveR,lg_clst,xclst,num_clst
c      write(22,*) aveR,xf/(1.0+4.0*ypep),lg_clst,xclst,num_clst
      if(ihb .eq. 0) then
         write(88,7588) aveR,xxf,xclst,hb_energy(hb_pick(ihb+1))
      else
         write(88,7588) aveR,xxf,xclst,hb_energy(hb_pick(ihb))
      endif

c ajr 9.10.01 following 2 files along with 81 & 80 were normal output until 
c today.  Now outfile (88) meancoord contains isopep corrected <r> and f,
c size of largest cluster/natom and energy of bond removed.
c      write(83,7583) aveR,xf/(1.0+4.0*ypep),avgR,xf,ypep         
c      write(85,7583) aveR,xf/(1.0+4.0*ypep),hb_energy(hb_pick(ihb))         

cc      write(83,7483) avgR,xf,xpep,xbod,xf+xpep+xbod+xn1,totatm,n1
cc      write(83,7383) avgR,jtest,iflop,nrsum,jtcmol,totatm,n1,isum1,nbody
cc  ^^---- jtest should be exactly zero in above counting check.
cc         write(83,7183) avgR,totatm,jtcmol,xnr,xf
c         write(85,7002) avgR,iflop,totatm,ncmode,halfcoord,xf
cc         write(8,*) ihb,natom,ns,totatm,avgR
      enddo 
c      close(80)
c      close(81)
c      close(83)
c      close(85)
      close(88)
      return
      end

      subroutine eclazz(nf,n)
      external neighb
      integer n,nf(n)
      
      nf(1) = 1
      do jj =2,n
         nf(jj) = jj
         do kk=1,jj-1
            nf(kk) = nf(nf(kk))            
            if( neighb(jj,kk).eq. 1) nf(nf(nf(kk)))=jj
         enddo
      enddo
      do jj =1,n
         nf(jj)=nf(nf(jj))
c         write(8,*) jj,nf(jj)
      enddo
      return
      end

c      subroutine ebody(nf,n,mult,maxn)
c      external neighb
c      integer n,nf(n),mult(maxn)
      
c      do ii=1,n
c         if(mult(ii) .gt. 0) nf(ii) = ii
c         nf(ii) = -1
cc         if(mult(ii) .gt. 0) nf(ii) = 0
c      enddo

c      do jj = 1,n
c         if(nf(jj) .ne. -1) then
c            nf(jj) = jj
c            if(jj .ne. 1) then
c               do kk=1,jj-1
c                  nf(kk) = nf(nf(kk))            
c                  if( neighb(jj,kk).eq. 1) nf(nf(nf(kk)))=jj
c               enddo
c            endif
c         endif
c      enddo
c      do jj =1,n
c         nf(jj)=nf(nf(jj))
cc         write(8,*) jj,nf(jj)
c      enddo
c      return
c      end

      integer function neighb(s1,s2)
c      logical function neighbor(s1,s2)
      include      'set_parameter'
      integer      s1,s2
      integer      link(nb2max),point(maxatm)
      integer      mult(maxatm)

      common/topology/ ns,link,point,mult      
      neighb = 0
c      neighbor = .FALSE.
      do i = 1,mult(s1)
         if(s2 .eq. link(point(s1) +i)) neighb = 1
c         if(s2 .eq. link(point(s1) +i)) neighbor = .TRUE.
      enddo
      
      end


      subroutine pruneit()
      include      'set_parameter'
      integer      s,so,msf,xsf,lxsf,sf,sf0,xso
      integer      link(nb2max),point(maxatm)
      integer      mult(maxatm),maxmult(maxatm)
      common/topology/ ns,link,point,mult
      common/prune/    maxmult
c -------------------------------------------------- prune off deadends
      nATM = ns
      do 20 s = 1,nATM
         so = s
   10    if(mult(so) .eq. 1) then                
            do mso = 1,maxmult(so)
               sf = link(point(so) +mso)
               if(sf .ne. -1) xso = mso
            enddo
            sf = link(point(so)+xso)
            link(point(so)+xso) = -1
            do msf=1,maxmult(sf)
               if(link(point(sf)+msf) .eq. so) xsf = msf
            enddo
            link(point(sf)+xsf) = -1
            mult(so) = mult(so) -1
            mult(sf) = mult(sf) -1
            if(mult(sf) .eq. 1) then
               so = sf
               goto 10
            endif
         endif
   20 continue
      ns = nATM
      return
      end

      subroutine condenseMult()
      include      'set_parameter'
      integer      link(nb2max),point(maxatm),so,sf,sfnew
      integer      sp,nill
      integer      mult(maxatm),maxmult(maxatm),holding(izcf)
      common/topology/ ns,link,point,mult
      common/prune/    maxmult

c     changed condenseMult to match changes in hbcondenseMult
c     AJR 09.28.01, what follows is the good version and error check.
      do s =1,ns
         if(mult(s) .gt. 0) then
            sp = point(s)
            do k =1,mult(s)
               holding(k) = 0
            enddo
            kcount =0
            do i =1,maxmult(s)
               sf = link(sp+i)
               if(sf .gt. 0) then
                  kcount = kcount+1
                  holding(kcount) = sf
               endif
            enddo
            jcount =0
            do i=1,mult(s)
               jcount = jcount+1
               link(sp+i) = holding(jcount)
            enddo
            do i = mult(s)+1,maxmult(s)
               link(sp+i) = -1
            enddo
c     error check for hbcondenseMult
            nill = 0
            do j =1,mult(s)
               so = link(point(s)+j)
               if(so .eq. -1) nill = nill+1
            enddo
            if(nill .gt. 0) write(33,*) s,nill,maxmult(s),mult(s)
         endif
      enddo     
      return
      end

c      do so=1,ns
c         if(maxmult(so) .gt. 1) then
c         do mso = 1,maxmult(so)-1
c            sf = link(point(so)+mso)
c            if(sf .eq. -1) then
c               sfnew = link(point(so)+mso+1)
c               if(sfnew .ne. -1) then
c                  link(point(so)+mso+1) = sf
c                  link(point(so)+mso) = sfnew
c               endif
c            endif
c         enddo
c         endif
c      enddo
c      return
c      end


      subroutine siteprune(s1)
      include      'set_parameter'
      integer      s1,s2,x2,sf,x1
      integer      link(nb2max),point(maxatm)
      integer      mult(maxatm),maxmult(maxatm)
      common/topology/ ns,link,point,mult
      common/prune/    maxmult

   10 do im = 1,maxmult(s1)
         if(link(point(s1)+im) .ne. -1) x1 = im
      enddo
      s2 = link(point(s1)+x1)
      do m2=1,maxmult(s2)
         if(link(point(s2)+m2) .eq. s1) x2 = m2
      enddo
      link(point(s2)+x2) = -1
      link(point(s1)+x1) = -1
      mult(s1) = mult(s1) -1
      mult(s2) = mult(s2) -1
      if(mult(s2) .eq. 1) then 
         s1 = s2 
         goto 10
      endif
      return
      end

      subroutine prunesc(aname,gflag,group)
      include      'set_parameter'
      integer      s,msf,xsf,lxsf,sf,g,gf
      integer      link(nb2max),point(maxatm),group(maxatm)
      integer      slow(maxgrp),sbig(maxgrp),gflag(maxgrp),gnum(maxgrp)
      integer      mult(maxatm),maxmult(maxatm)
      character*4  aname(maxatm),atm,ats
      character*3  gname(maxgrp)
      character*5  idres(maxgrp)
      character*1  chainid(maxgrp)
      common/topology/ ns,link,point,mult
      common/prune/    maxmult
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
                  sf = link(point(s)+ix)
                  if(sf.ne. -1) then
                  atm = aname(sf)
                  gf = group(sf)
               if(gf .eq. g) then
                  if((atm .ne. ' N  ') .and.(atm .ne. ' O  ') .and. 
     &               (atm .ne. ' C  ') .and.(atm .ne. ' CA ') .and.
     &               (atm .ne. 'OXT')) then
                     mult(s) = mult(s) -1
                     link(point(s)+ix) = -1
                     do ixf=1,maxmult(sf)
                        if(link(point(sf)+ixf) .eq. s) then
                           mult(sf) = mult(sf) -1
                           link(point(sf)+ixf) = -1
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

      subroutine ringbreak(g,aname)
      include      'set_parameter'
      integer      g,gfval,sg0,sgX,s,sbk,srm,xsd,site
      integer      link(nb2max),point(maxatm)
      integer      slow(maxgrp),sbig(maxgrp),gflag(maxgrp),gnum(maxgrp)
      integer      mult(maxatm),maxmult(maxatm)
      character*4  aname(maxatm),atm
      character*3  gname(maxgrp)
      character*5  idres(maxgrp)
      character*1  chainid(maxgrp)
      common/topology/ ns,link,point,mult
      common/prune/    maxmult
      common/groups/   ngrp,nres,nhet,nwater,chainid,idres,
     &                 gname,gnum,slow,sbig

      jcount = 0
   10 do s = slow(g),sbig(g)
         if(mult(s).ne.0) sbk = s
      enddo
      mbk = mult(sbk)
      do isk = 1,mbk
         srm = link(point(sbk)+isk)
         do msd = 1,mult(srm)
            if(link(point(srm)+msd) .eq. sbk) xsd = msd
         enddo
         mult(sbk) = mult(sbk) -1
         mult(srm) = mult(srm) -1
         link(point(sbk)+isk) = -1
         link(point(srm)+xsd) = -1
         if(mult(srm) .eq. 1) call siteprune(srm)
      enddo
      jcount=jcount+1
c ------------------------------------------------ breaks 2nd ring in TRP      
      if((jcount .lt. 2).and.(gname(g) .eq. 'TRP')) goto 10
      return
      end
