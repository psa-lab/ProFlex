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
c      subroutine get_files(nfile)
      subroutine get_files(nfile,iuseopt)
c ------------------------------------------------------------------------------
c Revision: 08.28.02 AJ Rader --
c     added output for the "allbond.ijkw" data see notes in outputfirst.f
c Revision: 04.12.02 AJ Rader --
c     incorporated iuseopt to eliminte user input questions
c     added chime output files: outfile(9) through outfile(11)
c
c                                             PROGRAM WRITTEN:      May 29, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     This subroutine is called first by 'first'. Its purpose is to automate
c the selection of output file names following a name format convention such 
c that a data set can be generated. Also, this gives a mechanism for passing
c data from the C++ program into the FORTRAN program. Eventually, everything
c will be written in C++. For the moment, an inefficient but robust method
c to determine the file names for the data series is implemented. The input
c is obtained from a special file that links FORTRAN with C++.
c ------------------------------------------------------------------------------
c INPUT: 
c    1) name of the original proflexdataset file as read from a temporary file. 
c    2) nfile
c ------------------------------------------------------------------------------
c INTERNAL: 
c    1) scans the current directory and selects appropriate output file names.
c    2) determines the next sequence number ranging from 1 to 9999. 
c    3) removes temporary file 
c -------------------------------------------- file names stemming from wxyz.pdb
c   input_data(1:n_end) = "wxyz_proflexdataset"
c   outfile(1)(1:n_end) = "wxyz_h-bonds.ijkl"
c   outfile(2)(1:n_end) = "wxyz_fdecomp.ijkl"
c   outfile(3)(1:n_end) = "wxyz_rdecomp.ijkl"
c   outfile(4)(1:n_end) = "wxyz_sdecomp.ijkl"
c   outfile(5)(1:n_end) = "wxyz_bond_wt.ijkl"
c   outfile(6)(1:n_end) = "wxyz_graphic.ijkl"
c   outfile(7)(1:n_end) = "wxyz_Rscript.ijkl"
c   outfile(8)(1:n_end) = "wxyz_analysis.log"
c                                                    AJR 04.01.02
c   outfile(9)(1:n_end) = "wxyz_fig_ijkl.pdb"
c   outfile(10)(1:n_end)= "wxyz_fig_ijkl.htm"
c   outfile(11)(1:n_end)= "wxyz_txt_ijkl.htm"
c                                                     AJR 08.28.02
c   outfile(12)(1:n_end)= "wxyz_allbond.ijkl"
c ------------------------------------------------------------------------------
c                              Variable List 
c file_id = a string variable representing the ijkl identification schem.
c nfile = current file number in the series of output data files.
c n_base = length of string up to the period just before the ijkl file id.
c n_end = full length of standard length of string output/input file name.
c n_root = length of string defining the root name of protein. (wxyz).
c ==============================================================================
      character*80 input_data,outfile(12)
      character*13 analysis
      character*21 dataset
      character*16 hbond,harlem
      character*9  fdecomp,rdecomp,sdecomp
      character*9  bondwt,graphic,Rscript
      character*4  file_id,htme,pdbe
      character*5  figu,text,main
      character*1  blank,digit(0:9)
      integer      i,j
      common/fnames/ n_root,input_data,outfile,file_id,digit
c ------------------------------------------------------------ format statements
 6000 format(a80)
 6005 format(5x,'ERROR: --> getfiles.f  Cannot find'
     *          ,' _proflexdataset')
 6010 format(5x,'Cannot have any characters after' 
     *          ,' "_proflexdataset"')
 6015 format(5x,'File name should be less than 80 characters long')
 6020 format(76x,i4)
c --------------------------------------------- initialize some string variables
      hbond   = '_h-bonds_SEfilt.'
      bondwt  = '_bond_wt.'
      rdecomp = '_rdecomp.'
      sdecomp = '_sdecomp.'
      fdecomp = '_fdecomp.'
      graphic = '_graphic.'
      Rscript = '_Rscript.'
      dataset = '_proflexdataset'
      analysis = '_analysis.log'
      harlem  = '_allbond.'

      digit(0) = '0'
      digit(1) = '1'
      digit(2) = '2'
      digit(3) = '3'
      digit(4) = '4'
      digit(5) = '5'
      digit(6) = '6'
      digit(7) = '7'
      digit(8) = '8'
      digit(9) = '9'
      blank = ' '
c AJR 04.01.02
      pdbe = '.pdb'
      htme = '.htm'
      text = '_txt_'
      figu = '_fig_'
      main = '_main'
c ------------------------------- open temporary file generated from C++ program
      open(1,file='qXyZaB.proflexdataset',status='old')
      rewind(1)
c ------------------------------------------------------ read name of input file
      read(1,6000) input_data
      read(1,*) nfile
c  AJR 03.22.02  the value of iuseopt is setup to allow for non interactive runs
c  of FIRST in the standard mode. (later kopt = 1)
      read(1,*) iuseopt
      close(1)
         do i=1,68
         j = i + 14
c -------------------------------------------------------- calculate root length
            if( input_data(i:j) .eq. dataset ) then
            n_root = i-1
            n_end = n_root + 15
            k = n_end + 1
               if( input_data(k:k) .ne. blank ) then
               write(6,6005)
               write(6,6010)
               call system('\\rm qXyZaB.proflexdataset')
               stop
               endif
            goto 50
            endif
         enddo
c ------------------------------------------------------ ERROR: cannot find file
      write(6,6005)
      write(6,6015)
      stop
   50 continue
      call system('\\rm qXyZaB.proflexdataset')
c -------------------------------------- construct the set of input/output names
      i4 = nfile/1000
      i3 = (nfile - i4*1000)/100
      i2 = (nfile - i4*1000 - i3*100)/10
      i1 =  nfile - i4*1000 - i3*100 - i2*10
         do iout=1,12
c         do iout=1,8
         outfile(iout)(1:n_root) = input_data(1:n_root)
         enddo
      i = n_root + 1

c ---------------------- 2006:03:02 SN 
c 	Due to renaming of files, the lengths of various files have changed
c	Hence, the following changes.
      j = i + 15
      outfile(1)(i:j)  = hbond
      outfile(1)(j+1:j+1)  = digit(i4)
      outfile(1)(j+2:j+2)  = digit(i3)
      outfile(1)(j+3:j+3)  = digit(i2)
      outfile(1)(j+4:j+4)  = digit(i1)


      j = i + 8
      outfile(2)(i:j) = fdecomp
      outfile(3)(i:j) = rdecomp
      outfile(4)(i:j) = sdecomp
      outfile(5)(i:j) = bondwt
      outfile(6)(i:j) = graphic
      outfile(7)(i:j) = Rscript
      outfile(12)(i:j) = harlem
      j = j + 1
      outfile(2)(j:j) = digit(i4)
      outfile(3)(j:j) = digit(i4)
      outfile(4)(j:j) = digit(i4)
      outfile(5)(j:j) = digit(i4)
      outfile(6)(j:j) = digit(i4)
      outfile(7)(j:j) = digit(i4)
      outfile(12)(j:j) = digit(i4)
      j = j + 1
      outfile(2)(j:j) = digit(i3)
      outfile(3)(j:j) = digit(i3)
      outfile(4)(j:j) = digit(i3)
      outfile(5)(j:j) = digit(i3)
      outfile(6)(j:j) = digit(i3)
      outfile(7)(j:j) = digit(i3)
      outfile(12)(j:j) = digit(i3)
      j = j + 1
      outfile(2)(j:j) = digit(i2)
      outfile(3)(j:j) = digit(i2)
      outfile(4)(j:j) = digit(i2)
      outfile(5)(j:j) = digit(i2)
      outfile(6)(j:j) = digit(i2)
      outfile(7)(j:j) = digit(i2)
      outfile(12)(j:j) = digit(i2)
      j = j + 1
      outfile(2)(j:j) = digit(i1)
      outfile(3)(j:j) = digit(i1)
      outfile(4)(j:j) = digit(i1)
      outfile(5)(j:j) = digit(i1)
      outfile(6)(j:j) = digit(i1)
      outfile(7)(j:j) = digit(i1)
      outfile(12)(j:j) = digit(i1)

      j = i + 12
      outfile(8)(i:j) = analysis

c ------------------- 2006:03:02  SN's changes end ---------------------

c AJR 08.28.02 additions for harlem output
c      outfile(13)(1:n_root) = input_data(1:n_root)
c      outfile(12)(i:n_end) = harlem

c AJR 04.01.02 additions for chime output
      k = i+4 
c      outfile(9)(i:k)  = main
      outfile(9)(i:k) = figu
      outfile(10)(i:k) = figu
      outfile(11)(i:k) = text
c      outfile(12)(i:k) = text
      k = k+1
      outfile(9)(k:k) = digit(i4)
      outfile(10)(k:k) = digit(i4)
      outfile(11)(k:k) = digit(i4)
c      outfile(12)(k:k) = digit(i4)
      k = k+1
      outfile(9)(k:k) = digit(i3)
      outfile(10)(k:k) = digit(i3)
      outfile(11)(k:k) = digit(i3)
c      outfile(12)(k:k) = digit(i3)
      k = k+1
      outfile(9)(k:k) = digit(i2)
      outfile(10)(k:k) = digit(i2)
      outfile(11)(k:k) = digit(i2)
c      outfile(12)(k:k) = digit(i2)
      k = k+1
      outfile(9)(k:k) = digit(i1)
      outfile(10)(k:k) = digit(i1)
      outfile(11)(k:k) = digit(i1)
c      outfile(12)(k:k) = digit(i1)
      k = k+1
      n_end = n_root + 13
      outfile(9)(k:n_end) = pdbe
      outfile(10)(k:n_end) = htme
      outfile(11)(k:n_end) = htme
c      outfile(12)(k:n_end) = htme

      file_id(1:1) = digit(i4)
      file_id(2:2) = digit(i3)
      file_id(3:3) = digit(i2)
      file_id(4:4) = digit(i1)


c     The following code is used for FIRSTweb cgi scripts
c      open(10,file='output_files.html',status='unknown')
c      do a = 1, 12
c         write(10,6000) outfile(a)
c      enddo
c      
c      write(10,6000) input_data
c      close(10)

      return
      end




