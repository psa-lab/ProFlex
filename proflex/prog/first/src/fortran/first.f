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
      program  first
c ------------------------------------------------------------------------------
c  Address inquiries to:  proflex@sol.bch.msu.edu
c
c                                                  PROGRAM WRITTEN:  7.24.97
c                                                  LAST MODIFIED:    3.22.02 AJR
c                                                  LAST MODIFIED:    5.  .08 SkN
c ------------------------------------------------------------------------------
c                              DESCRIPTION

c     The FIRST program is a utility tool to determine the Floppy Inclusion and
c Rigid Substructure Topography of macromolecules such as proteins.
c ------------------------------------------------------------------------------
c                               FIRST1_2 
c INPUT:   
c    1) <file>_FIRSTdataset

c OUTPUT:
c    outfile(1) = "wxyz_h-bonds.ijkl"
c    outfile(2) = "wxyz_fdecomp.ijkl"
c    outfile(3) = "wxyz_rdecomp.ijkl"
c    outfile(4) = "wxyz_sdecomp.ijkl"
c    outfile(5) = "wxyz_bond_wt.ijkl"
c    outfile(6) = "wxyz_graphic.ijkl"
c    outfile(7) = "wxyz_Rscript.ijkl"
c    outfile(8) = "wxyz_analysis.log"
c    outfile(9) = "wxyz_fig_ijkl.pdb"
c    outfile(10)= "wxyz_fig_ijkl.htm"
c    outfile(11)= "wxyz_txt_ijkl.htm"
c ==============================================================================
c                          List of Definitions
c ------------------------------------------------------------------------------
c                             Variable List 
c aname(s) = A four character name describing atom s. 
c bval(s) = the "B-value" or the "temperature factor" found in the PDB file.
c chainid(g) = defines the chain that group g belongs to.
c rec(s) = Single character discription of PDB record, (eg. A => ATOM,
c             H => HETATM)
c freq(s) = gives the occupancy value found in the PDB records.
c g = distinct group number for a set of atoms, for example may be a residue.
c gname(g) = A three character name describing group g. 
c gnum(g) = a "group number" as specified in the original pdb file for group g. 
c      In general group numbers are repeated when there is more than one chain.
c hb = a dummy index used for Hydrogen bonds
c hbond(j,hb) = <so-s-sf>  for j=1,2,3 respectively. so => donor atom label,
c      s => Hydrogen atom label and sf => acceptor atom label.
c khb = number of Hydrogen bonds placed in macromolecule. 
c linkref() = nearest neighbor table of the bond-bending network.
c        The j-th nearest neighbor to site so is given by: linkref(point(so)+j)
c linknoh() = j-th central-force nearest neighbor to atom s in the molecule 
c        without Hydrogen bonds present.  
c locgroup(s) = distinct group number, g, for which atom s belongs.
c maxgrp = maximum number of distinct groups of atoms allowed.
c maxh = maximum number of Hydrogen bonds that can be considered.
c multnoh(s) = # of [covalent] nearest neighbors that atom s has. Note that no
c        Hydrogen bonds are included.  
c multref(s) = # of nearest neighbors that atom s has, including H-bonds.
c natom = # of atoms in the given structure of the macromolecule. 
c nclst = number of rigid clusters, including isolated sites and dimers
c nfile = current file number in the series of data files
c ngrp = total number of distinct groups found in the macromolecule.
c nh = number of hinge joints within the molecule.
c nhb = number of [potential] Hydrogen bonds.
c n_one = number of isolated sites
c nsc = number of lines describing selection criteria
c n_two = number of isolated dimers
c het1 = first index for a non-standard group corresponding to a
c        HETATM-field in the pdb formated  XXXX_FIRSTdataset input file. 
c nhet = number of HETATM type atoms.
c nwater = number of water molecules floating around.
c point() = used to give appropriate index in linkref()
c pointer() = used to give appropriate index in linknoh()
c s,so,sf = site labels within the bond bending network.
c slow(g),sbig(g) = the range of site labels (slow .LE. s .LE. sbig) that
c      belong to group g. 
c tc = a dummy index used for torsion constraints that lock dihedral motions.
c tcmol = number of torsion constraints to be applied within the molecule, 
c      usually due to a peptide or resonant bond.
c th = a dummy index used for torsion constraints that lock dihedral angle 
c      motions about hinge joints, to characterize collectve motions.
c th_min = (tcmol + 1) = the minimum value for index th
c th_max = (tcmol + nh) = the maximum value for index th
c torsion(j,tc) = <so-sf>  for j=1,2 respectively. Defines a pair of sites where
c      a torsion constraint is to be placed. Also used to define hinge joints
c      (between th_min  and  th_max) before the torsion constraint is added.
c ------------------------------------------------------------------------------
c ==============================================================================
c                         DEFINE DATA STRUCTURE
      include      'set_parameter'
      integer      tcmol,th_min,th_max
      integer      locgroup(maxatm),slow(maxgrp),sbig(maxgrp)
      integer      gnum(maxgrp)
      integer      linknoh(nb2max),pointer(maxatm)
      integer      linkref(nb2max),point(maxatm)
      integer      multnoh(maxatm),multref(maxatm)
      integer      torsion(2,nbmax), runtype
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      integer      clst(maxatm),stress(-1:maxatm)
      integer      link_hinge(nb2max),label_hinge(nbmax)               
      dimension    sxyz(3,maxatm),freq(maxatm),bval(maxatm)
      dimension    hb_energy(maxh)
      character*80 firstdata,outfile(12),criteria(99)
      character*80 decomp_file, bondwt_file
      character*5  idres(maxgrp)
      character*4  aname(maxatm),file_id
      character*3  gname(maxgrp)
      character*1  rec(maxatm),chainid(maxgrp),hb_type(maxh)
      character*1  archive,digit(0:9)

      common/select/   nsc,criteria
      common/fnames/   n_root,firstdata,outfile,file_id,digit
      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval
      common/groups/   ngrp,nres,nhet,nwater,chainid,
     &                 idres,gname,gnum,slow,sbig
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/numbers1/ nflop,n_one,n_two,nbody,nclst,nstress
      common/numbers2/ iflop,nhinge,nih,ncmode
      common/hbanalys/ HBanalysis,decomp_file,bondwt_file
      common/network0/ nsnoh,linknoh,pointer,multnoh
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/topology/ ns,linkref,point,multref
      common/subgraph/ stress
      common/cmotions/ label_hinge,link_hinge
      common/clusters/ clst

c ------------------------------------------------------------ format statements
 4000 format(4x,'ANALYSIS MENU')
 4004 format(2x,'What would you like to perform:')
 4001 format(4x,'(1) Flexibility and rigidity analysis')
 4002 format(4x,'(2) Hydrogen bond dilution')
 4003 format(4x,'(3) Bond stripping')
 4005 format(i4)
 4103 format(4x,'Sorry that option is not currently available.')
 6005 format(5x,'Enter    c    to continue the FIRST analysis')
 6100 format(a1)
      
      iuseopt = 0

c -------------------------------------------------------- automatic file naming
c added 'iuseopt' AJR 03.22.02 to allow noninteractive use of FIRST
      call get_files(nfile,iuseopt)
c ------------------------------------------------------------------- enter data
      call read_data(firstdata)
c --------------------------------------------------------------- select H-bonds
      call read_hb_list(outfile(1))
c -------------------------------------- inserted option to allow for hbdilution
c =================AJR 03.22.02 
c 'iuseopt' bypasses runtype options to allow noninteractive use of FIRST
c
      if(iuseopt .eq. 1) then
         runtype = 2
      elseif(iuseopt .eq. 2) then
         runtype = 1
      else
         call system('clear')
 100     write(6,*)
         write(6,4000)
         write(6,*)
         write(6,4004)
         write(6,*)
         write(6,4001)
         write(6,4002)
c         write(6,4003)		--- 2006:03	Sandeep (disable bnd striping)
         read(5,4005) runtype
         if((runtype .lt. 1) .or. (runtype .gt. 2)) goto 100
      endif
c  the code below can be removed once option 3 is working and distributed.
c      if(runtype .eq. 3) then
c         write(6,4103)
c         goto 100
c      endif
c runtype = 1. Run FIRST according to the original design.
      if(runtype .eq. 1) then 
c ------------------------------------ make neighbor table with selected H-bonds
         call place_Hbond()
c -------------------------------------------------- macromolecule decomposition
         call mapmolecule()
c ------------------------------------------------------------ summarize results
         call outputfirst(archive,iuseopt)
c ----------------------------------------------------- update output file names

c runtype = 2.  The user selected the HBdilution set of routines. All of these 
c subroutines are located in a file named hbdilute.f . There are several different 
c options now, all of which rely on removing hydrogen bonds from the flexibility 
c calculations in a continuous manner, and observing the resultant effects on 
c protein flexibility. In turn, the output for most of the hbdilute routines is 
c created by the anciliary program hbdilute.c written in C. 
      elseif(runtype .eq. 2) then
c =================Sameer Arora 01.04.04 
c 'iuseopt' bypasses runtype options to allow noninteractive use of FIRST
c	 call hbdilute()
         call hbdilute(iuseopt)

c runtype =3.  Again protein structures have hydrogen bonds removed one at a 
c time and various schemes of stripping away dangling ends can be envoked by 
c choosing this option.
      elseif(runtype .eq. 3) then 
         call pruning()
      endif
      stop
      end
      

