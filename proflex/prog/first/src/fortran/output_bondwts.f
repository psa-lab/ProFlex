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

c Revision 1.2  1999/04/02 21:12:47  don
c removed reference to motion() as this array was never used.
c
c Revision 1.1  1999/03/30 22:22:40  don
c Initial revision
c
c Revision 1.2  1999/03/17  16:36:03  jacobs
c Fixed up some format statements with missing commas.
c Modified the  hbonds  common block to incorporate energy information.
c Fixed up the calling of subroutine out_PDBfile() to eliminate direct
c reference to hbond() list by passing hbonds common block instead.
c Also, modified the way the list of H-bonds is recorded per realization.
c
c Revision 1.1  1998/12/20  23:06:12  jacobs
c Initial revision
c
c
      subroutine output_bondwts(archive)
c ------------------------------------------------------------------------------
c                                             PROGRAM WRITTEN:     July  8, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     This subroutine writes to standard output (the screen) a short summary of 
c the FIRST results. It determines the cluster statistics for the rigid clusters,
c Laman subgraphs, and collective motions. It generates a coloring assignment to 
c the rigid cluster decomposition, Laman subgraph decompostion and collective 
c motion decomposition. At this time, it constructs an "optimal" site-relabeling 
c scheme so that the RasMol script file can be constructed in a very few number 
c of lines. These results are outputed to the following files:

c outfile(1)  = "h-bonds.ijkl" = List of all selected H-bonds
c outfile(3)  = "decompR.ijkl" = single column file defining rigid cluster 
c               decomposition in terms of the atom labels from _FIRST.chem
c outfile(4)  = "decompL.ijkl" = single column file defining the Laman subgraph
c               decomposition in terms of the atom labels from _FIRST.chem
c outfile(2)  = "decompC.ijkl" = single column file defining the collective 
c               motions in terms of the atom labels from _FIRST.chem
c outfile(5)  = "bond_wt.ijkl" = a four column file: so  sf  f_i  r_i which 
c               defines for each CF-bond, <so,sf>, a local degree of floppiness
c               and a local degree of redundancy.
c outfile(6)  = "grphcsR.ijkl" = new PDB file to plot rigid cluster decomposition
c               with the atom labels scrambled to simplify script file. 
c outfile(10) = "grphcsL.ijkl" = new PDB file to plot Laman subgraph decomposition
c               with the atom labels scrambled to simplify script file. 
c outfile(9)  = "grphcsC.ijkl" = new PDB file to plot the collective motions
c               with the atom labels scrambled to simplify script file. 
c outfile(7)  = "scriptR.ijkl" = corresponding script file for grphcsR
c outfile(11) = "scriptL.ijkl" = corresponding script file for grphcsL
c outfile(12) = "scriptC.ijkl" = corresponding script file for grphcsC
c outfile(8) = "analysis.log" = an ongoing summary file
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
      integer      link_hinge(nb2max),label_hinge(nbmax)                  
      integer      g,h,hb,r,s,sf,so,tc,tcmol,th_max,th_min
      integer      gnum(maxgrp),slow(maxgrp),sbig(maxgrp)
      integer      linkref(nb2max),point(maxatm),locgroup(maxatm)
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      integer      torsion(2,nbmax)
      integer      newlabel(maxatm),oldlabel(maxatm)
      integer      size(maxatm),multarc(maxatm)
      integer      clst(maxatm),stress(-1:maxatm)
      integer      pointarc(maxatm),pointbrc(maxatm),multbrc(maxatm)
      integer      mult(maxr),pointr(maxr)
      integer      cslow(0:maxcolor),csbig(0:maxcolor)
      integer      calow(maxcolor),cabig(maxcolor)
      integer      multref(maxatm),color(maxatm),paint(maxcolor)
      dimension    sxyz(3,maxatm),freq(maxatm),bval(maxatm)
      dimension    hb_energy(maxh)
      character*80 input_data,outfile(12),criteria(99)
      character*80 fscript,fgrphcs,fdecomp,line
      character*13 cname(maxcolor)
      character*5  idres(maxgrp)
      character*4  aname(maxatm),file_id
      character*3  gname(maxgrp)
      character*1  rec(maxatm),chainid(maxgrp),hb_type(maxh)
      character*1  digit(0:9),archive,answer0
      common/fnames/   n_root,input_data,outfile,file_id,digit
      common/select/   nsc,criteria
      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval
      common/groups/   ngrp,nres,nhet,nwater,chainid,
     &                 idres,gname,gnum,slow,sbig
      common/numbers1/ nflop,n_one,n_two,nbody,nclst,nstress
      common/numbers2/ iflop,nhinge,nih,ncmode
      common/topology/ ns,linkref,point,multref
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/subgraph/ stress
      common/cmotions/ label_hinge,link_hinge
      common/clusters/ clst
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
 6000 format('---------------------< SELECTION CRITERIA >',
     &       '-------------------------------- ',a4)
 6005 format(26x,'SELECTION CRITERIA')
 6010 format(a80)
 6095 format('------------------------------------------',
     &       '---------------------------')
 6100 format(26x,'BASIC SUMMARY')
 6105 format(2x,i8,' Atoms')
 6110 format(2x,i8,' Hydrogen bonds')
 6115 format(2x,i8,' Torsional constraints')
 6120 format(2x,i8,' Independent degrees of freedom (DOF)')
 6125 format(2x,i8,' Rigid clusters')
 6130 format(2x,i8,' Isolated atoms')
 6135 format(2x,i8,' Isolated Dimers')
 6140 format(2x,i8,' Three-dimensional objects')
 6145 format(2x,i8,' Network induced stress regions')
 6150 format(2x,i8,' Floppy modes (internal independent DOF)')
 6155 format(2x,i8,' Hinge joints')
 6160 format(2x,i8,' Independent hinge joints')
 6165 format(2x,i8,' Collective motions (with two or more hinges)')
 6200 format(a4,1x,a3,1x,a1,i6,15x,a4,1x,a3,1x,a1,i6,3x,e14.7)
 6250 format('---------- Rigid cluster Size ----- ',
     &       'Cluster # ------------')
 6255 format(14x,i9,10x,i9)
 6300 format('---------- Stress Region Size ----- ',
     &       'Cluster # ------------')
 6350 format('---- Collective Motion ---- # of hinges ',
     &       '---- # of DOF ----')
 6355 format(10x,i9,10x,i9,7x,i7)
 6360 format(3i10)
 6400 format(2i10,5x,f9.5)
 6450 format(5x,'Please wait:  Writing grphcs, decomp & ',
     &          'RasMol script files')
 6500 format(10x,' Archive decomposition?   Enter Y/N  ',
     &             '--> Y is default') 
 6505 format(a1)

      archive = 'y'
      if( archive .eq. 'n' ) archive = 'N'
      if( archive .ne. 'N' ) then
         nhinge0 = nhinge
c     -------------------------------------------------------- write to output files
         open(59,file="bond_wts",access='append')
c     ----------------------------------------------------- write selection criteria
         
c     ---------------------------------------------------- temporarily change clst()
         k = 0
         do so=1,natom
            if( multref(so) .eq. 1 ) then
               sf = linkref( point(so) + 1 )
               clst(so) = clst(sf)
            endif
            if( clst(so) .eq. so ) then
               k = k + 1
               newlabel( clst(so) ) = k
            endif
         enddo
         
         if( k .ne. nclst ) then
            write(6,*) k,nclst
            stop
         endif
         
         do so=1,natom
            clst(so) = newlabel( clst(so) )
         enddo
         
c     ----------------------------- count # of BULK & SURFACE sites for each cluster
         do nc=1,nclst
            size(nc) = 0
            multarc(nc) = 0
            multbrc(nc) = 0
         enddo
         do so=1,natom
            ipo = point(so)
            nco = clst(so)
c     ---------------------------------------------- count BULK sites in cluster nco
            multarc(nco) = multarc(nco) + 1
            do jo=1,multref(so)
               indexo = ipo + jo
               sf = linkref(indexo)
c     ------------------------------------------------- prevent double counting bond
               if( sf .gt. so ) then
                  ncf = clst(sf) 
                  if( nco .ne. ncf ) then
c     ----------------------------------------------------------- count SURFACE site 
                     size(nco) = size(nco) + 1
                     size(ncf) = size(ncf) + 1
                  else
c     --------------------------------------------------------- count BULK BOND-site
                     multbrc(nco) = multbrc(nco) + 1
                  endif
               endif
            enddo
         enddo
c     --------------------------------------------- calculate number of hinge joints
c     2*nhinge= sum_{nc=1}^{nclst} [# of SURFACE sites in clst(nc)]
         nhinge = 0
         maxsize = -1
         do nc=1,nclst
            nhinge = nhinge + size(nc) 
c     -------------------------------------------------- size = SURFACE + BULK sites
            size(nc) = size(nc) + multarc(nc)
c     ------------------------------------------------------------ find maximum size 
            if( maxsize .lt. size(nc) ) maxsize = size(nc) 
         enddo
         nhinge = nhinge/2
c     ---------------------------------------------------- initialize histogram bins
c     use oldlabel() for binning the data
         do s=1,maxsize
            oldlabel(s) = 0
         enddo
c     --------------------------------------------- calculate cluster size histogram 
         do nc=1,nclst
            oldlabel( size(nc) ) = oldlabel( size(nc) ) + 1
         enddo
         
c     ----------------------------------------------------------- GLOBAL ERROR CHECK
         kflop = 3*natom + nhinge
         do nc=1,nclst
            isize = size(nc)
            if( isize .gt. 2 ) then
               kflop = kflop - (3*isize - 6)
            elseif( isize .eq. 2 ) then
               kflop = kflop - 1
            endif
         enddo
c     --------------------------------------------------------------- ERROR FOUND!!!
         if( kflop .ne. nflop ) then
            write(6,*) 
            write(6,*) 
            write(6,*) 
            write(6,5000) kflop
            stop
         endif
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ FINISHED ERROR CHECK

c     --------------------------- relabel rigid clusters in descending order on size
c     use pointarc and pointbrc as dummy work arrays
c     set up pointers for range in cluster size
c     define nc2,nc1
         call clstrlabels(size,pointarc,nclst,pointbrc,maxsize,clst,
     &        newlabel,natom,mult,pointr,maxr,nc2,nc1)
c     ---------------------------------------------------------------------- relabel
         do s=1,natom
            clst(s) = newlabel( clst(s) )
         enddo
c     --------------------------- use pointarc,pointbrc,oldlabel as temporary arrays
         do nc=1,nclst
            kc = newlabel(nc)
            pointarc(kc) = multarc(nc)
            pointbrc(kc) = multbrc(nc)
            oldlabel(kc) = size(nc)
         enddo
c     ------------------------------------ transfer back into the appropriate arrays
         do nc=1,nclst
            multarc(nc) = pointarc(nc)
            multbrc(nc) = pointbrc(nc)
            size(nc) = oldlabel(nc) 
         enddo
c     --------------------------------------------------------- color rigid clusters
c     use pointBrc as a chain for linking BULK sites in a cluster
c     use pointArc to bin all BULK sites in a cluster
c     define cname()
         call colorclustr(clst,pointbrc,natom,pointarc,color,nclst,
     &        mult,pointr,cname)
c     ---------------------------------------------------------------- relabel sites 
c     for quick graphics display
c     INPUT:  use pointBrc as a chain for linking BULK sites in a cluster
c     INPUT:              use pointArc to bin all BULK sites in a cluster
c     NOTE:   pointArc,pointBrc are set up from subroutine colorclustr()
c     OUTPUT: pointArc,pointBrc,oldlabel,newlabel,paint,ncolor,calow,cabig 
c     and cslow,csbig
         
         call relabelatms(nhinge,clst,pointbrc,oldlabel,newlabel,natom,
     &        pointarc,color,nclst,paint,ncolor,calow,cabig,cslow,csbig)
         
c     ------------------------------------------------ calculate bond weight factors
         do s=1,natom
            multarc(s) = 0
            multbrc(s) = 0
            size(s) = 0
         enddo
         do so=1,natom
            mc = 0
            index = point(so)
            istress = stress(so)
            multarc(istress) = multarc(istress) + 1                  
            do jo=1,multref(so)
               index = index + 1
               sf = linkref(index)
c     -------------------------------------------------------- allow double counting
               if( istress .eq. stress(sf) ) then 
                  size(istress) = size(istress) + 1    
                  mc = mc + 1
               endif
            enddo
c     multbrc(istress) = multbrc(istress) + mc*(mc-1)/2  
            multbrc(istress) = multbrc(istress) + 2*mc-3   
         enddo
c     ---------------------------------------------- finialize CF and BB constraints
         do s=1,natom
            size(s) = size(s)/2
            multbrc(s) = multbrc(s) + size(s) 
         enddo
c     ------------------------------------------- augment dihedral angle constraints
         do tc=1,tcmol
            so = torsion(1,tc)
            sf = torsion(2,tc)
            istress = stress(so)
            if( istress .eq. stress(sf) ) then
               multbrc(istress) = multbrc(istress) + 1
            endif
         enddo
         do so=1,natom-1
            index = point(so)
            istress = stress(so)
            do jo=1,multref(so)
               index = index + 1
               sf = linkref(index)
c     ------------------------------------------------------ prevent double counting
               if( sf .gt. so ) then
                  if( istress .eq. stress(sf) ) then 
                     nr = multbrc(istress) - (3*multarc(istress) - 6) 
                     x = - float(nr)/float( size(istress) )
                     if( x .gt. -1.0e-5 ) write(89,*) so,sf,istress,
     &                   size(istress),multbrc(istress),multarc(istress)
                     write(59,6400) so,sf,x
                  endif
               endif
            enddo
         enddo
c     -------------------------------------------------- temporarily change stress()
         k = 0
         do so=1,natom
            if( stress(so) .eq. so ) then
               k = k + 1
               newlabel(so) = k
            endif 
         enddo
c     ---------------------------------------------- use nclst and clst() variables
         nclst = k
         do so=1,natom
            clst(so) = newlabel( stress(so) )
         enddo
c     ---------------------------------------- count # of BULK sites for each region
         do nc=1,nclst
            multarc(nc) = 0
            multbrc(nc) = 0
            size(nc) = 0
         enddo
         nhinge = 0
         do so=1,natom
            ipo = point(so)
            nco = clst(so)
c     ---------------------------------------------- count BULK sites in cluster nco
            multarc(nco) = multarc(nco) + 1
            do jo=1,multref(so)
               indexo = ipo + jo
               sf = linkref(indexo)
c     ------------------------------------------------- prevent double counting bond
               if( sf .gt. so ) then
                  ncf = clst(sf) 
                  if( nco .eq. ncf ) then
c     --------------------------------------------------------- count BULK BOND-site
                     multbrc(nco) = multbrc(nco) + 1
                  else
c     ----------------------------------------------------------- bond is unstressed
                     nhinge = nhinge + 1
c     ----------------------------------------------------------- count SURFACE site 
                     size(nco) = size(nco) + 1
                     size(ncf) = size(ncf) + 1
                  endif
               endif
            enddo
         enddo
c     ----------------------------------------- calculate number of unstressed bonds
         maxsize = -1
c     ------------------------------------------------------------ find maximum size 
         do nc=1,nclst
            s = multarc(nc)
            if( maxsize .lt. s ) maxsize = s
         enddo
c     ---------------------------------------------------- initialize histogram bins
c     use oldlabel() for binning the data
         do s=1,maxsize
            oldlabel(s) = 0
         enddo
c     --------------------------------------------- calculate cluster size histogram 
         do nc=1,nclst
            s = multarc(nc)
            size(nc) = size(nc) + s                       
            oldlabel(s) = oldlabel(s) + 1
         enddo
c     --------------------------- write to analysis file the cluster size statistics
         kstress = 0
         do s=2,maxsize
            kstress = kstress + oldlabel(s)
         enddo
c     ----------------------------------------------------------- GLOBAL ERROR CHECK
         if( kstress .ne. nstress ) then
c     --------------------------------------------------------------- ERROR FOUND!!!
            write(6,*) 
            write(6,*) 
            write(6,5005) kstress
            stop
         endif
c     -------------------------------------------------------- find NEW maximum size 
         maxsize = -1
         do nc=1,nclst
            if( maxsize .lt. size(nc) ) maxsize = size(nc)
         enddo
c     --------------------------- relabel rigid clusters in descending order on size
c     use pointarc and pointbrc as dummy work arrays
c     set up pointers for range in cluster size
c     define nc2,nc1
         call clstrlabels(size,pointarc,nclst,pointbrc,maxsize,clst,
     &        newlabel,natom,mult,pointr,maxr,nc2,nc1)
c     ---------------------------------------------------------------------- relabel
         do s=1,natom
            clst(s) = newlabel( clst(s) )
         enddo
c     --------------------------- use pointarc,pointbrc,oldlabel as temporary arrays
         do nc=1,nclst
            kc = newlabel(nc)
            pointarc(kc) = multarc(nc)
            pointbrc(kc) = multbrc(nc)
            oldlabel(kc) = size(nc)
         enddo
c     ------------------------------------ transfer back into the appropriate arrays
         do nc=1,nclst
            multarc(nc) = pointarc(nc)
            multbrc(nc) = pointbrc(nc)
            size(nc) = oldlabel(nc) 
         enddo
c     --------------------------------------------------------- color rigid clusters
c     use pointBrc as a chain for linking BULK sites in a cluster
c     use pointArc to bin all BULK sites in a cluster
c     define cname()
         call colorclustr(clst,pointbrc,natom,pointarc,color,nclst,
     &        mult,pointr,cname)
c     ---------------------------------------------------------------- relabel sites 
c     for quick graphics display
c     INPUT:  use pointBrc as a chain for linking BULK sites in a cluster
c     INPUT:              use pointArc to bin all BULK sites in a cluster
c     NOTE:   pointArc,pointBrc are set up from subroutine colorclustr()
c     OUTPUT: pointArc,pointBrc,oldlabel,newlabel,paint,ncolor,calow,cabig 
c     and cslow,csbig
         
         call relabelatms(nhinge,clst,pointbrc,oldlabel,newlabel,natom,
     &        pointarc,color,nclst,paint,ncolor,calow,
     &        cabig,cslow,csbig)

c     --------------------------- write to analysis file the cluster size statistics
         do s=th_min,th_max
            multbrc(s) = 0
            size(s) = 0
         enddo
         do tc=th_min,th_max
            s = label_hinge(tc)
            if( s .gt. 0 ) then
               size(s) = size(s) + 1
               if( torsion(2,tc) .lt. 0 ) multbrc(s) = multbrc(s) + 1
            endif
         enddo
         do nh=1,nhinge0
            oldlabel(nh) = 0
         enddo
         do so=th_min,th_max
            s = size(so) 
            if( s .gt. 0 ) oldlabel(s) = oldlabel(s) + 1
         enddo
c     ----------------------------------------------------------- GLOBAL ERROR CHECK
         k = 0
         do nh=1,nhinge0
            k = k + nh*oldlabel(nh)
         enddo
         k = nhinge0 - k
         if( k .ne. nih ) then
            write(6,*) 
            write(6,*) 
            write(6,5010) k
            stop
         endif
         kcmode = 0
         idof = nih
         do s=th_min,th_max
            if( size(s) .gt. 0 ) then
               kcmode = kcmode + 1
               idof = idof + multbrc(s)
            endif
         enddo
c     ----------------------------------------------------------- GLOBAL ERROR CHECK
         if( kcmode .ne. ncmode ) then
c     --------------------------------------------------------------- ERROR FOUND!!!
            write(6,*) 
            write(6,*) 
            write(6,5015) kcmode
            stop
         endif
         if( idof .ne. iflop ) then
c     --------------------------------------------------------------- ERROR FOUND!!!
            write(6,*) 
            write(6,*) 
            write(6,5020) idof
            stop
         endif
c     -------------------------------------------------------- write to bond_wt file
         do tc=th_min,th_max
            so = torsion(1,tc)
            sf = iabs( torsion(2,tc) )
            s = label_hinge(tc)
            if( s .gt. 0 ) then
               x = float( multbrc(s) )/float( size(s) )
            else
               x = 1.0e0
            endif
            write(59,6400) so,sf,x
         enddo
c     ---------------------------------------- temporarily make a fictitious clst()
         do s=1,natom
c  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> edit 02.25.02 AJR_BMH: output dangling ends to bond_wt file
            if(multref(s) .eq. 1) then
               x = 1.0e0
               write(59,6400) s,linkref(point(s)+1),x
            endif
c  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end output edit 02.25.02 AJR_BMH
            clst(s) = 1
         enddo
         nclst = 1
         do tc=th_min,th_max
            newlabel(tc) = -1
         enddo
         do tc=th_min,th_max
            s = label_hinge(tc)
            if( s .gt. 0 ) then
               if( newlabel(s) .lt. 0 ) then
                  nclst = nclst + 1
                  newlabel(s) = nclst
               endif
               so = torsion(1,tc)
               sf = iabs( torsion(2,tc) )
               clst(so) = newlabel(s)
               clst(sf) = newlabel(s)
            endif
         enddo
         fdecomp = outfile(2)
         do tc=th_min,th_max
            so = torsion(1,tc)
            sf = iabs( torsion(2,tc) )
            s = label_hinge(tc)
            if( s .gt. 0 ) then
               label = newlabel(s) - 1
            else
               label = -1
            endif
         enddo
c     ---------------------------------------- count # of BULK sites for each region
         do nc=1,nclst
            multarc(nc) = 0
            multbrc(nc) = 0
            size(nc) = 0
         enddo
         nhinge = 0
         do so=1,natom
            ipo = point(so)
            nco = clst(so)
c     ---------------------------------------------- count BULK sites in cluster nco
            multarc(nco) = multarc(nco) + 1
            do jo=1,multref(so)
               indexo = ipo + jo
               sf = linkref(indexo)
c     ------------------------------------------------- prevent double counting bond
               if( sf .gt. so ) then
                  ncf = clst(sf) 
                  if( nco .eq. ncf ) then
c     --------------------------------------------------------- count BULK BOND-site
                     multbrc(nco) = multbrc(nco) + 1
                  else
c     ----------------------------------------------------------- bond is unstressed
                     nhinge = nhinge + 1
c     ----------------------------------------------------------- count SURFACE site 
                     size(nco) = size(nco) + 1
                     size(ncf) = size(ncf) + 1
                  endif
               endif
            enddo
         enddo
c     ----------------------------------------- calculate number of unstressed bonds
         maxsize = -1
c     ------------------------------------------------------------ find maximum size 
         do nc=1,nclst
            s = multarc(nc)
            if( maxsize .lt. s ) maxsize = s
         enddo
c     ---------------------------------------------------- initialize histogram bins
c     use oldlabel() for binning the data
         do s=1,maxsize
            oldlabel(s) = 0
         enddo
c     --------------------------------------------- calculate cluster size histogram 
c     BULK sites only
         do nc=1,nclst
            s = multarc(nc)
            size(nc) = size(nc) + s                           
            oldlabel(s) = oldlabel(s) + 1
         enddo
c     ----------------------------------------------------------- GLOBAL ERROR CHECK
         kcmode = -1                                       
         do s=2,maxsize
            kcmode = kcmode + oldlabel(s)
         enddo
         if( kcmode .ne. ncmode ) then
c     --------------------------------------------------------------- ERROR FOUND!!!
            write(6,*) 
            write(6,*) 
            write(6,5015) kcmode
            stop
         endif
c     -------------------------------------------------------- find NEW maximum size 
         maxsize = -1
         do nc=1,nclst
            if( maxsize .lt. size(nc) ) maxsize = size(nc)
         enddo
c     --------------------------- relabel rigid clusters in descending order on size
c     use pointarc and pointbrc as dummy work arrays
c                                      set up pointers for range in cluster size
c                                                                 define nc2,nc1
c      call clstrlabels(size,pointarc,nclst,pointbrc,maxsize,clst,
c     &                 newlabel,natom,mult,pointr,maxr,nc2,nc1)
c ---------------------------------------------------------------------- relabel
c      do s=1,natom
c      clst(s) = newlabel( clst(s) )
c      enddo
c --------------------------- use pointarc,pointbrc,oldlabel as temporary arrays
c         do nc=1,nclst
c         kc = newlabel(nc)
c         pointarc(kc) = multarc(nc)
c         pointbrc(kc) = multbrc(nc)
c         oldlabel(kc) = size(nc)
c         enddo
c ------------------------------------ transfer back into the appropriate arrays
c      do nc=1,nclst
c      multarc(nc) = pointarc(nc)
c      multbrc(nc) = pointbrc(nc)
c      size(nc) = oldlabel(nc) 
c      enddo
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
c
c      call relabelatms(nhinge,clst,pointbrc,oldlabel,newlabel,natom,
c     &    pointarc,color,nclst,paint,ncolor,calow,
c     &    cabig,cslow,csbig)
c --------------------------------------------- write corresponding output files
c         fdecomp = outfile(4)
c         fgrphcs = outfile(8)
c         fscript = outfile(11)
c         itemp = natom
c         call out_script(fscript,fgrphcs,pointarc,multarc,
c     &     pointbrc,multbrc,nclst,calow,cabig,cslow,csbig,
c     &     paint,ncolor,cname,mult,pointr,nc1,nc2)
c         call out_PDBfile(fgrphcs,newlabel,
c     &                    itemp,clst,pointbrc,multbrc,nclst)
         write(59,*)'END'
         close(59)
      else
c     ------------------------------------------------------ do not writeout results
      endif
      
      return
      end
      
