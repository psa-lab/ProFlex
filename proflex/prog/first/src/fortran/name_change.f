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
      subroutine name_change(nfile)
c ------------------------------------------------------------------------------
c                                             PROGRAM WRITTEN:     JUNE 23, 1997
c                                       program written by:     Donald J. Jacobs
c                                                              jacobs@@pa.msu.edu
c ------------------------------------------------------------------------------
c                               Description
c     This subroutine is called by FIRST1_0 whenever the user wants to continue
c to analize the macromolecule. Its purpose is to automate the selection of 
c output file names following a name format convention such that a data set can 
c be generated. In particular, the names of three file types must be updated by
c incrementing the file index,  nfile. 
c ------------------------------------------------------------------------------
c INPUT: 
c    1) nfile ---> the old file index
c    2) six old file names. 
c ------------------------------------------------------------------------------
c OUTPUT:  
c    1) new file index
c    2) six new file names. 
c ------------------------------------------------------------------------------
c                              Variable List 
 
c nfile = current file number in the series of data files
c outfile(1)  = "wxyz_h-bonds.ijkl"
c outfile(2)  = "wxyz_fdecomp.ijkl"
c outfile(3)  = "wxyz_rdecomp.ijkl"
c outfile(4)  = "wxyz_sdecomp.ijkl"
c outfile(5)  = "wxyz_bond_wt.ijkl"
c outfile(6)  = "wxyz_graphic.ijkl"
c outfile(7)  = "wxyz_Rscript.ijkl"
c ==============================================================================
      character*80 input_data,outfile(12)
      character*4  file_id
      character*1  digit(0:9)
      common/fnames/ n_root,input_data,outfile,file_id,digit
c ------------------------------------------------------------ format statements
 6000 format(5x,'CONGRATULATIONS you have reached the 9999-th ',
     &          'combination!!')
 6005 format(5x,'SORRY, but more combinations cannot be added ',
     &          'to the existing set of files')
 6010 format(5x,'To continue; copy the reference files with a ',
     &          'new root name')
c ------------------------------------------------ add a new run to the data set
      nfile = nfile + 1
         if( nfile .gt. 9999 ) then
         call system( 'clear' )
         write(6,*)
         write(6,*)
         write(6,6000)
         write(6,6005)
         write(6,6010)
         stop
         endif
c ------------------------------------------ construct a set of new output names
      i4 = nfile/1000
      i3 = (nfile - i4*1000)/100
      i2 = (nfile - i4*1000 - i3*100)/10
      i1 =  nfile - i4*1000 - i3*100 - i2*10
      j = n_root + 10 
      outfile(1)(j+7:j+7) = digit(i4)
      outfile(2)(j:j) = digit(i4)
      outfile(3)(j:j) = digit(i4)
      outfile(4)(j:j) = digit(i4)
      outfile(5)(j:j) = digit(i4)
      outfile(6)(j:j) = digit(i4)
      outfile(7)(j:j) = digit(i4)
      j = j + 1
      outfile(1)(j+7:j+7) = digit(i3)
      outfile(2)(j:j) = digit(i3)
      outfile(3)(j:j) = digit(i3)
      outfile(4)(j:j) = digit(i3)
      outfile(5)(j:j) = digit(i3)
      outfile(6)(j:j) = digit(i3)
      outfile(7)(j:j) = digit(i3)
      j = j + 1
      outfile(1)(j+7:j+7) = digit(i2)
      outfile(2)(j:j) = digit(i2)
      outfile(3)(j:j) = digit(i2)
      outfile(4)(j:j) = digit(i2)
      outfile(5)(j:j) = digit(i2)
      outfile(6)(j:j) = digit(i2)
      outfile(7)(j:j) = digit(i2)
      j = j + 1
      outfile(1)(j+7:j+7) = digit(i1)
      outfile(2)(j:j) = digit(i1)
      outfile(3)(j:j) = digit(i1)
      outfile(4)(j:j) = digit(i1)
      outfile(5)(j:j) = digit(i1)
      outfile(6)(j:j) = digit(i1)
      outfile(7)(j:j) = digit(i1)

      file_id(1:1) = digit(i4)
      file_id(2:2) = digit(i3)
      file_id(3:3) = digit(i2)
      file_id(4:4) = digit(i1)
      return
      end
