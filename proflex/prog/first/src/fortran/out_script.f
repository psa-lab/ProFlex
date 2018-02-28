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
      subroutine out_script(fscript,fgrphcs,pointarc,multarc,
     &        pointbrc,multbrc,nclst,calow,cabig,cslow,csbig,
     &        paint,ncolor,cname,multr,pointr,nc1,nc2)
c ------------------------------------------------------------------------------
c AJR 03.27.02 eliminated the use of "CFB" or "BOND" pseudoatoms
c AJR 04.16.02 replaced references to "CFB" and "HCFB" since pseudoatoms are no longer here.
c ------------------------------------------------------------------------------
c                                             LAST UPDATED:       April 16, 2002
c                                             PROGRAM WRITTEN:     July  9, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c This subroutine generates a RasMol script file to color code the macromolecule 
c ------------------------------------------------------------------------------
c INPUT: 
c    1) output and graphics file names
c    2) nclst       ---> number of distinct rigid clusters.
c    3) pointarc(),multarc() 
c    4) pointbrc(),multbrc() 
c ------------------------------------------------------------------------------
c OUTPUT:  
c    1) creates  fscript  file. 
c ------------------------------------------------------------------------------
c                              Variable List 
 
c calow(icolor) = the smallest "new atom label" assigned to the icolor-th color.
c cabig(icolor) = the biggest "new atom label" assigned to the icolor-th color.
c cslow(icolor) = the smallest "new site label" assigned to the icolor-th color.
c csbig(icolor) = the biggest "new site label" assigned to the icolor-th color.
c     Note that all so called site labels > natom  which can run up to ns. 
c cname(i) = the name of the i-th color
c icolor = a dummy coloring index.
c maxcolor = maximum number of different colors used. 
c ==============================================================================
      include      'set_parameter'
      integer      so,sf
      integer      pointarc(nclst),multarc(nclst)
      integer      pointbrc(nclst),multbrc(nclst)
      integer      calow(maxcolor),cslow(0:maxcolor)
      integer      cabig(maxcolor),csbig(0:maxcolor)
      integer      multr(maxr),pointr(maxr)
      integer      paint(maxcolor)
      character*80 fscript,fgrphcs
      character*13 cname(maxcolor)
c ------------------------------------------------------------ format statements
 6000 format(2x,'#  RasMol 2.6 beta-2 Script file: ',
     &          '---- by Donald J. Jacobs')
 6001 format(2x,'#                                 ',
     &          '        Michigan State University')
 6002 format(2x,'#                                 ',
     &          '        jacobs@@pa.msu.edu')
 6005 format(2x,'#  Floppy Inclusion and Rigid Substructure Topography')
 6010 format(2x,'#  RasMol script file for FIRST analysis')
 6015 format(2x,'#  PDB graphics file required as input')
c ------------------------------------------------------------------- initialize
 6020 format(2x,'zap')
 6025 format(2x,'reset')
 6030 format(2x,'echo  Loading ',a80)
 6035 format(2x,'load ',a80)
 6040 format(2x,'set background white')
 6045 format(2x,'wireframe ',f4.2)
 6050 format(2x,'select CFB')
 6055 format(2x,'wireframe off')
c
c AJR 03.27.02 ----- all you need are comments for each atom in each color group
c
 5000 format(2x,'select all')
 5100 format(2x,'spacefill ',f4.2)
 5060 format(2x,'select atomno>=',i7,' and atomno<=',i7)
 5065 format(2x,'color ',a13)
 5115 format(2x,'select backbone')
 5118 format(2x,'wireframe ',f4.2)
 5120 format(2x,'color bond black')
 5125 format(2x,'wireframe 120')
 5400 format(2x,'restrict backbone')
c --------------------------------------------------------------- generate cloud
 6059 format(2x,'echo "Coloring the FIRST analysis"')
 6060 format(2x,'#select atomno>=',i7,'  and  atomno<=',i7)
 6061 format(2x,'#select atomno>=',i7,'  and  atomno<=',i7)
 6065 format(2x,'#color ',a13)
 6066 format(2x,'#color ',a13)
 6070 format(2x,'define ATM  not CFB')
 6075 format(2x,'select ATM')
 6080 format(2x,'#set radius 0.5')
 6085 format(2x,'#dots 1000')
c ------------------------------------------------------------------ color atoms
 6090 format(2x,'color cpk')
 6100 format(2x,'spacefill ',f4.2)
c ------------------------------------------------------------------ color bonds
 6105 format(2x,'select CFB, (atomno>=',i7,' and  atomno<=',i7,')')
c 6105 format(2x,'select 9999, (atomno>=',i7,' and  atomno<=',i7,')')
 6110 format(2x,'color bond ',a13)
c ------------------------------------------------------------ backbone emphasis
 6115 format(2x,'select backbone,CFB')
c -------------------------------------------------------------------- H-bonding
 6120 format(2x,'color bond black')
 6125 format(2x,'dash')
c ------------------------------------------------------- echo group definitions
 6495 format(2x,'echo "---------------------------------',
     &                   '----------------------------"')
 6500 format(2x,'echo "Setting up user defined groups:"')
 6505 format(2x,'echo "Group Name     Definition"')
 6510 format(2x,'echo "   CFB    ---> sites representing ',
     &                                    'CENTRAL-FORCE BOND"')
 6515 format(2x,'echo "   ATM    ---> all ATOMS within the ',
     &                                         'macromolecule"')
 6520 format(2x,'echo "   HJ     ---> all HINGE JOINTS between ',
     &                                    'the rigid clusters"')
c 6525 format(2x,'echo "   ARC',i1,'   ---> bulk ATOMS within RIGID',
 6525 format(2x,'echo "    RC',i1,'   ---> all ATOMS within RIGID',
     &                                       ' CLUSTER ',i1,'"')
c 6530 format(2x,'echo "   ARC',i2,'  ---> bulk ATOMS within RIGID ',
 6530 format(2x,'echo "    RC',i2,'   ---> all ATOMS within RIGID',
     &                                        'CLUSTER ',i2,'"')
 6535 format(2x,'echo "   BRC',i1,'   ---> bulk BONDS within RIGID',
     &                                       ' CLUSTER ',i1,'"')
 6540 format(2x,'echo "   BRC',i2,'  ---> bulk BONDS within RIGID ',
     &                                        'CLUSTER ',i2,'"')
 6545 format(2x,'echo "   RC',i1,'    ---> all atoms and bonds ',
     &                           'within RIGID CLUSTER ',i1,'"')
 6550 format(2x,'echo "   RC',i2,'   ---> all atoms and bonds ',
     &                           'within RIGID CLUSTER ',i2,'"')
 6555 format(2x,'echo "   CLEAN  ---> everything except isolated ',
     &                                    'sites and bonds"')
 6560 format(2x,'echo "   SITE1  ---> only the isolated sites"')
 6565 format(2x,'echo "   BOND1  ---> only the isolated bonds"')
 6570 format(2x,'echo "   TETHER ---> all hydrophobic pseudoatoms"')
      
c ------------------------------------------------------------------ definitions
 7000 format(2x,'define HCFB (atomno>=',i7,' and  atomno<=',i7,')')
c 7005 format(2x,'define ARC',i1,' (atomno>=',i7,
 7005 format(2x,'define RC',i1,' (atomno>=',i7,
     &          ' and  atomno<=',i7,')')
c 7010 format(2x,'define ARC',i2,' (atomno>=',i7,
 7010 format(2x,'define RC',i2,' (atomno>=',i7,
     &          ' and  atomno<=',i7,')')
 7015 format(2x,'define BRC',i1,' (atomno>=',i7,
     &          ' and  atomno<=',i7,')')
 7020 format(2x,'define BRC',i2,' (atomno>=',i7,
     &          ' and  atomno<=',i7,')')
 7025 format(2x,'define RC',i1,' ARC',i1,',BRC',i1,', within(0.9, ',
     &          '(within(0.9,ARC',i1,') and HCFB) )')
 7030 format(2x,'define RC',i2,' ARC',i2,',BRC',i2,', within(0.9, ',
     &          '(within(0.9,ARC',i2,') and HCFB) )')
c ------------------------------------------------------------- define utilities
 7035 format(2x,'define SITE1 (atomno>=',i7,
     &                   ' and atomno<=',i7,')')
 7040 format(2x,'define aa (atomno>=',i7,
     &                   ' and atomno<=',i7,')')
 7045 format(2x,'define BOND1 aa, (atomno>=',i7,
     &                   ' and atomno<=',i7,')')
 7050 format(2x,'define bb SITE1, BOND1')
 7055 format(2x,'define CLEAN not bb')
c 7060 format(2x,'define HJ HCFB, (within(0.9,HCFB) and ATM)')
 7100 format(2x,'select atm')
 7105 format(2x,'spacefill off')
c 7110 format(2x,'restrict cfb,backbone')
 7110 format(2x,'restrict 9999,backbone')
 7200 format(2x,'select hetero and not *.X')

c ------------------------------------------------------------ open fscript file
      open(2,file=fscript,status='unknown') 
      rewind(2)
c ----------------------------------------------------- write Header information
      write(2,*)
      write(2,6000)
      write(2,6001)
      write(2,6002)
      write(2,*) 
      write(2,*) 
      write(2,6005)
      write(2,6010) 
      write(2,6015) 
c --------------------------------------------------------- basic setup commands
      write(2,*)
      write(2,6020)
      write(2,6025)
      write(2,6030) fgrphcs
      write(2,6035) fgrphcs
c -------------------------------------------------- custom design macromolecule
      write(2,*)
      write(2,6040) 
      write(2,6045) 0.07 
c -------------------------------------------------- declare working on coloring
      write(2,*)
      write(2,6059)
c ------------------------------------------------------ Generate cloud of color
      write(2,*)
c ------------------------------------------------------------ color all regions
         do iorder=1,ncolor
         icolor = paint(iorder) 
         write(2,5060) calow(icolor),cabig(icolor)
         write(2,5065) cname(icolor)
c         write(2,6060) calow(icolor),cabig(icolor)
c         write(2,6065) cname(icolor)
         enddo
c ---------------------------------------------- do not put a cloud on CFB sites
c             do iorder=1,ncolor
c             icolor = paint(iorder) 
c             write(2,6061) cslow(icolor),csbig(icolor)
c             write(2,6066) cname(icolor)
c             enddo
c ------------------------------------------------------------------- apply DOTS
c      write(2,6070)
c      write(2,6075)
c      write(2,6080) 
c      write(2,6085)
c      write(2,*) 
c ------------------------------------------------------- color atoms & set size
c      write(2,*) 
c      write(2,6090)
c      write(2,6100) 0.3 
c ------------------------------------------------------------------ color bonds
c      write(2,*) 
c         do iorder=1,ncolor
c         icolor = paint(iorder) 
c         write(2,6105) calow(icolor),cabig(icolor)
c         write(2,6110) cname(icolor)
c         enddo
c ----------------------------------------------------------- emphasize backbone
      write(2,*) 
      write(2,5115)
c      write(2,6115)
      write(2,6045) 0.12
c ---------------------------------------------------------------- place H-bonds
      write(2,*) 
      write(2,5000)
      write(2,5100) 0.12 
c      write(2,5120)
c      write(2,5125)
c      write(2,6075)
c      write(2,6120)
c      write(2,6125)
c ------------------------------------------------------ eliminate CFB-CFB bonds
c      write(2,*) 
c      write(2,6050) 
c      write(2,6055) 
c --------------------------------------------------- define user defined groups
      write(2,*) 
      write(2,6500) 
      write(2,6505) 
c ---------------------------------------- first define rigid cluster properties
c      write(2,7000) cslow(0),csbig(0)
      write(2,6495) 
      nc_cut = pointr(4) + multr(4)
         do nc=1,nc_cut
c ------------------------------------------------------------------- define ARC
         so = pointarc(nc) + 1
         sf = pointarc(nc) + multarc(nc)
            if( nc .lt. 10 ) then
            write(2,6525) nc,nc
            write(2,7005) nc,so,sf
            elseif( nc .lt. 100 ) then
            write(2,6530) nc,nc
            write(2,7010) nc,so,sf
            endif
c ------------------------------------------------------------------- define BRC
c         so = pointbrc(nc) + 1
c         sf = pointbrc(nc) + multbrc(nc)
c            if( nc .lt. 10 ) then
c            write(2,6535) nc,nc
c            write(2,7015) nc,so,sf
c            elseif( nc .lt. 100 ) then
c            write(2,6540) nc,nc
c            write(2,7020) nc,so,sf
c            endif
c -------------------------------------------------------------------- define RC
c            if( nc .lt. 10 ) then
c            write(2,6545) nc,nc
c            write(2,7025) nc,nc,nc,nc 
c            elseif( nc .lt. 100 ) then
c            write(2,6550) nc,nc
c            write(2,7030) nc,nc,nc,nc
c            endif
         enddo
c ---------------------------------------------------- echo standard definitions
      write(2,*) 
      write(2,6495) 
c      write(2,6500) 
c      write(2,6505) 
c      write(2,6495) 
c      write(2,6510) 
c      write(2,6515) 
c      write(2,6520) 
c      write(2,6555) 
c      write(2,6560) 
c      write(2,6565)
c      write(2,6570)  
c         if( nc1 .le. nclst ) then
cc ----------------------------------------------------------------- define SITE1
c         so = pointarc(nc1) + 1 
c         sf = pointarc(nclst) + 1                                  
c         write(2,7035) so,sf
cc ----------------------------------------------------------------- define BOND1  
c         so = pointarc(nc2) + 1 
c         sf = pointarc(nc1)                                      
c         write(2,7040) so,sf
c         so = pointbrc(nc2) + 1 
c         sf = pointbrc(nc1)                               
c         write(2,7045) so,sf
c         elseif( nc2 .le. nclst ) then
cc ----------------------------------------------------------------- define SITE1
c         so = cslow(0)
c         sf = cslow(0) - 1 
c         write(2,7035) so,sf
cc ----------------------------------------------------------------- define BOND1  
c         so = pointarc(nc2) + 1 
c         sf = pointarc(nclst) + 2                    
c         write(2,7040) so,sf
c         so = pointbrc(nc2) + 1 
c         sf = pointbrc(nclst) + 1                         
c         write(2,7045) so,sf
c         else
c         so = cslow(0)
c         sf = cslow(0) - 1 
c         write(2,7035) so,sf
c         write(2,7040) so,sf
c         write(2,7045) so,sf
c         endif
c ----------------------------------------------------------------- define CLEAN
c      write(2,7050)
c      write(2,7055)
cc -------------------------------------------------------------------- define HJ
c      write(2,7060)
c ---------------------------------------------- add in usual selection criteria
      write(2,*)
c      write(2,7100)
c      write(2,7105)
c      write(2,*)
c      write(2,7110)
      write(2,5400)
      write(2,7200)
      write(2,5118) 0.09
      write(2,5100) 0.3
c ----------------------------------------------------------- close fscript file
      close(2)
      return
      end
