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
      subroutine out_chime(pointarc,multarc,pointbrc,multbrc,nclst,
     &              calow,cabig,paint,ncolor,cname,multr,pointr,nc1,nc2)

c ------------------------------------------------------------------------------
c                                             LAST UPDATED:       May    8, 2002
c                                             PROGRAM WRITTEN:    March 27, 2002
c                                       program written by:          A. J. Rader
c                                                               rader@pa.msu.edu
c     based upon the out_script file 
c ------------------------------------------------------------------------------
c                               Description
c This subroutine generates Chime script files to color code the macromolecule 
c via a webrowswer.
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
      integer      calow(maxcolor),cabig(maxcolor)
      integer      multr(maxr),pointr(maxr)
      integer      paint(maxcolor)
      character*80 fscript,fgrphcs
      character*80 input_data,outfile(12)
      character*13 cname(maxcolor),ajc(8)
      character*20 otitle,oname
      character*4 file_id
      character*1 digit(0:9)
      
      common/fnames/  n_root,input_data,outfile,file_id,digit
      mleng = n_root+13
 1111 format(a20)

c-----------------------------------------------newcolors
      ajc(1) = '[Xae00fe]    '
      ajc(2) = '[X808090]    '
      ajc(4) = 'orange       '                                        
      ajc(3) = 'greenblue    '
      ajc(5) = '[XAABB00]    '
      ajc(6) = 'brown        '
      ajc(7) = 'purple       '
      ajc(8) = 'cyan         ' 

      do iorder=1,ncolor
         ic = paint(iorder) 
      if(cname(ic) .eq. 'black        ') cname(ic) = '[X787878]    '
      enddo     

c ------------------------------------------------------------ format statements
c$$$ 6000 format(2x,'#  RasMol 2.6 beta-2 Script file: ',
c$$$     &          '---- by Donald J. Jacobs')
c$$$ 6001 format(2x,'#                                 ',
c$$$     &          '        Michigan State University')
c$$$ 6002 format(2x,'#                                 ',
c$$$     &          '        jacobs@@pa.msu.edu')
c$$$ 6005 format(2x,'#  Floppy Inclusion and Rigid Substructure Topography')
c$$$ 6010 format(2x,'#  RasMol script file for FIRST analysis')
c$$$ 6015 format(2x,'#  PDB graphics file required as input')
c ------------------------------------------------------------------- initialize
c
c AJR 03.27.02 ----- all you need are comments for each atom in each color group
c
 1000 format('<html>')
 1010 format('<title>',a20,'</title>')
 1015 format('<body>')
 1020 format('<body bgcolor="#FFFFFF">')
 1029 format('<embed src="')
 1030 format('" align=abscenter width=100% height=97% spiny=20 ',
     &   'startspin = true color3d=cpk display=spacefill name="rasmol"')
 1034 format(3x,'script="set background white;')
 1035 format('load ',a80,';')
 1036 format(6x,'select all;')
 1136 format(6x,'cartoon off;')
 1037 format(6x,'spacefill 0.5; wireframe 0.2;')
c 1037 format(6x,'spacefill;')
 1038 format(6x,'spacefill off;')
 1500 format(6x,'select atomno>=',i7,' and atomno<=',i7,';')
 1550 format(6x,'color ',a13,';')
 2005 format(6x,'define RC',i1,' (atomno>= ',i6,
     &          ' and  atomno<= ',i6,');')
 2010 format(6x,'define RC',i2,' (atomno>= ',i6,
     &          ' and  atomno<= ',i6,');')
 2015 format(6x,'select RC',i1,';')
 2020 format(6x,'select RC',i2,';')
 3000 format(6x,'restrict protein;">')
 3100 format('</body>')
 3110 format('</html>')
 5125 format(6x,'wireframe 100;')
 5400 format(6x,'restrict backbone;')
 3150 format(6x,'spin off;')
c ------------------------ ajr 04.01.01 new for 2-frame version
 4000 format('<embed type="application/x-spt" width=15 height=15 ',
     &       'align=bottom button="push" target="molecule"')
 4010 format('</embed>')
 4100 format('Default: color by Rigid Cluster Decomposition (RCD)')
 4150 format('Stop spinning')
 4200 format('RCD of backbone only')
 4210 format('Show backbone only')
 4300 format('Flexibility Index coloring')
 5999 format('">')

 4810 format('<p>')
 4811 format('</p>')

c ------------------------------------------------------------ open fscript file
c              actually both the fig_XXXX.htm and txt_XXXX.htm chime script files

      fscript =outfile(10)
      open(2,file=fscript,status='unknown') 
      rewind(2)
      open(3,file=outfile(11),status='unknown',access='append')
c ----------------------------------------------------- write Header information

c|--------------------------------------------------------------------------------|
c|  Create one file with both embedded rasmol controls and the structure - Sameer |
c|--------------------------------------------------------------------------------|
      write(3,1000)
      write(3,1020)
      write(3,*)
      write(3,*)
      write(3,1029)      
      if(mleng .le. 20) then
         write(3,1111) outfile(9)         
      else
         write(3,*) outfile(9)         
      endif
      write(3,1030)
      write(3,1034)
      write(3,1036)
      write(3,1037)
c ------------------------------------------------------------ color all regions
      do iorder=1,ncolor
         icolor = paint(iorder) 
         write(3,1500) calow(icolor),cabig(icolor)
         write(3,1550) cname(icolor)
      enddo
c ----------------------------------------------------------- emphasize backbone
      write(3,*)       
c ---------------------------------------- first define rigid cluster properties
      nc_cut = pointr(4) + multr(4)
      do nc=1,nc_cut
c ------------------------------------------------------------------- define RC
         so = pointarc(nc) + 1
         sf = pointarc(nc) + multarc(nc)
         if( nc .lt. 10 ) then
            write(3,2005) nc,so,sf
         elseif( nc .lt. 100 ) then
            write(3,2010) nc,so,sf
         endif
      enddo
c ---------------------------------------------------- echo standard definitions
      write(3,*)
c  ajr 03.29.02 only recolor rigid clusters if more than 3 RCs
c                        this really requires rewriting colorclstr.f to
      if(nc_cut .gt. 3) then
         do nc = 4,nc_cut
            so = pointarc(nc) + 1
            sf = pointarc(nc) + multarc(nc)
            if(nc .lt. 10) then
               write(3,2015),nc
               jco = nc-3
               write(3,1550),ajc(jco)
            elseif(nc .lt. 100) then
               write(3,2020),nc
c  quick fix:  ajr 05.08.02 
               if(nc .lt. 12) then
                  jco = nc -3
               elseif(nc .lt. 20) then
                  jco = nc -11
               else
                  jco = int(nc/100)
               endif
               write(3,1550),ajc(jco)
            endif
         enddo
      endif
      write(3,*)
      write(3,3000)
c -------------------------------------------------- finish fscript file

      write(3,4010)
c     write(3,4100)
c     write(3,4811)
c ----------------------------------------------------------- close fscript file
c|--------------------------------------------------------------------------------|
c| 		Single file fix ends - Sameer					  |
c|--------------------------------------------------------------------------------|


c ------------------------------------------ write left (molecule frame)
      write(2,1000)
      write(2,1020)
      write(2,*)
      write(2,*)
      write(2,1029)
      if(mleng .le. 20) then
         write(2,1111) outfile(9)
      else
         write(2,*) outfile(9)
      endif
      write(2,1030)
      write(2,1034)
      write(2,1036)
      write(2,1037)
c ---------------------------------  write right (color scheme frame)
c     write(3,1000)
c     write(3,1015)

      write(3,*)
      write(3,4810)
      write(3,4000)
      write(3,1034)
      write(3,1036)
      write(3,1136)
      write(3,1037)

c ------------------------------------------------------------ color all regions
      do iorder=1,ncolor
         icolor = paint(iorder) 
         write(2,1500) calow(icolor),cabig(icolor)
         write(2,1550) cname(icolor)
         write(3,1500) calow(icolor),cabig(icolor)
         write(3,1550) cname(icolor)
      enddo
c ----------------------------------------------------------- emphasize backbone
      write(2,*) 
c ---------------------------------------- first define rigid cluster properties
      nc_cut = pointr(4) + multr(4)
      do nc=1,nc_cut
c ------------------------------------------------------------------- define RC
         so = pointarc(nc) + 1
         sf = pointarc(nc) + multarc(nc)
         if( nc .lt. 10 ) then
            write(2,2005) nc,so,sf
            write(3,2005) nc,so,sf
         elseif( nc .lt. 100 ) then
            write(2,2010) nc,so,sf
            write(3,2010) nc,so,sf
         endif
      enddo
c ---------------------------------------------------- echo standard definitions
      write(2,*) 
      write(3,*)
c  ajr 03.29.02 only recolor rigid clusters if more than 3 RCs
c                        this really requires rewriting colorclstr.f to
      if(nc_cut .gt. 3) then
         do nc = 4,nc_cut
            so = pointarc(nc) + 1
            sf = pointarc(nc) + multarc(nc)
            if(nc .lt. 10) then
               write(2,2015),nc
               write(3,2015),nc
               jco = nc-3
               write(2,1550),ajc(jco)
               write(3,1550),ajc(jco)
            elseif(nc .lt. 100) then
               write(2,2020),nc
               write(3,2020),nc
c  quick fix:  ajr 05.08.02 
               if(nc .lt. 12) then
                  jco = nc -3
               elseif(nc .lt. 20) then
                  jco = nc -11
               else
                  jco = int(nc/100)
               endif
               write(2,1550),ajc(jco)
               write(3,1550),ajc(jco)
            endif
         enddo
      endif
      write(2,*)
      write(2,3000)
      write(3,*)
      write(3,3000)
c -------------------------------------------------- finish fscript file

      write(3,4010)
      write(3,4100)
      write(3,4811)
      write(2,*)
      write(2,3100)
      write(2,3110)
c ----------------------------------------------------------- close fscript file
      close(2)
      close(3)
      return
      end

c ====================================================== AJR 04.04.02
c program to create the chime script for flexibility index mapping
c

      subroutine chime_flex(matom,wt)

      include      'set_parameter'
      integer      so,sf,red,gre,blu,topc
      real         wt(matom)
      integer      wtbin(2001),wtgrp(2001),test,tempwt
c      integer      wtbin(-1000:1000),wtgrp(0:2000)
      character*80 input_data,outfile(12)

      character*4 file_id
      character*1 digit(0:9)
      
      common/fnames/  n_root,input_data,outfile,file_id,digit


 2000 format('<p>')
 2100 format('<embed type="application/x-spt" width=15 height=15 ',
     &       'align=bottom  button="push" target="molecule"')
 2134 format(3x,'script="set background white;')
 2135 format('load ',a80,';')
 2136 format(6x,'select all;')     
 2137 format(6x,'spacefill;')
 2138 format(6x,'spacefill off; wireframe off;')
 2200 format(6x,'cartoon;')
c 2200 format(6x,'ribbon;')
 2210 format(6x,'restrict backbone;')
 2220 format(6x,'cartoon off;')
c 2220 format(6x,'ribbon off;')
c 4100 format('default view: color by rigid decomposition')
 3000 format(6x,'select temperature=',i5,';')
 3010 format(6x,'color [',i3,',',i3,',',i3,'];')
 3100 format(6x,'color ',a13,';')
 3150 format(6x,'spin off;')
 4150 format('Stop spinning')
 4210 format('Show backbone only')
 4300 format('Color by Flexibility Index')
 4400 format('Turn off ribbon')
 5125 format(6x,'wireframe 100;')

 3800 format('">')
 3900 format('</embed>')
 4000 format('</p>')
 5000 format('</body>')
 5001 format('</html>')

      mleng = n_root+13
      topc = 255
      nrigid = 0
      nflex = 0

      do ip = 1,2001
         wtgrp(ip) = 9999
         wtbin(ip) = 0
c      do ip = -1000,1000
c         wtgrp(ip+1001) = 9999
c         wtbin(ip+1001) = 0
      enddo
      do s=1,matom
         tempwt = int(1000*wt(s))
         wtbin(tempwt+1001) = wtbin(tempwt+1001) + 1
      enddo
      nwts = 0
c 2121 format(i6,f8.4,i6)
c      do ii = -1000,1000
      do ii = 1,2001
         ip = ii - 1001
         if(wtbin(ii) .gt. 0) then
            nwts = nwts+1
            test = nwts+1001
            wtgrp(test) = ip
         endif
      enddo
c 2222 format(6x,'made it into chime_flex')
c      write(6,2222) 
c ----------------------------------------------------- open chime command file htm
      open(3,file=outfile(11),status='unknown',access='append')
c ----------------------------- add flex index options to this file.      
      write(3,2000)
      write(3,2100)
      write(3,2134)
      write(3,2136)
      write(3,2138)
      write(3,2200)
      write(3,2210)
      do ii = 1,nwts
         iwt = wtgrp(ii+1001)
c         write(6,*) iwt
         write(3,3000),iwt
         if(iwt .eq. 0) then 
            red = 164
            blu = 164
            gre = 164
         elseif(iwt .lt. -200) then
            red = 0
            gre = 0
            blu = topc + (iwt * 0.155)
c            write(6,*) red,gre,blu
         elseif(iwt. lt. 0) then
            red = 0
            gre = topc + iwt
c            gre = gre -iwt/10
            blu = topc           
         elseif(iwt .gt. 570) then
            red = topc - (iwt *0.015)
c         elseif(iwt .gt. 255) then
c            red = topc - (iwt * 0.155)
            gre = 0
            blu = 0
         elseif(iwt .gt. 340) then
            red = topc
            gre = 0
            blu = (iwt*0.75) - topc
         elseif(iwt .gt. 0) then 
            red = topc
            gre = topc - (iwt* 0.75)
            blu = 0
         endif
         write(3,3010) red,gre,blu
      enddo

      write(3,3800)
      write(3,3900)
      write(3,4300)
      write(3,4000)

      write(3,*)
      write(3,2000)
      write(3,2100)

      write(3,2134)
      write(3,2136)
      write(3,2138)
      write(3,2220)
      write(3,2210)
      write(3,5125)
      write(3,3800)
      write(3,3900)
      write(3,4210)
      write(3,4000)

 8000 format('  Right click on the molecule for other display options.')
      write(3,*)
      write(3,2000)
      write(3,8000)
      write(3,4000)
      write(3,5000)
      write(3,5001)

      close(3)
      return
      end




