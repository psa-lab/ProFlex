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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Print the header for the decomp_list file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine print_decomp_list_header(nres,outfile,decomp_file)

      include      'set_parameter'
      integer      nres, natom, locgroup(maxatm)
      character*1  rec(maxatm)
      character*4  aname(maxatm)
      character*80 decomp_file, outfile
      dimension    sxyz(3,maxatm),freq(maxatm),bval(maxatm)

      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval

 100  format('HEADER',1x,i7,1x,i7,3x,a100)

      open(91,file=decomp_file,status='new',access='sequential')
      write(91,100) nres, natom, outfile
      close( 91 )
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Print the A: blah blah stuff to the decomp file describing
c the particular Hbond that was just diluted from the network
c --------------------------------------------------------  
c     modified 10.10.01 by AJR
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine print_Hbond_info_to_decomp_file(jhb,hps,locgroup,
     &     nres,gnum,aname,chainid,decomp_file,mean_coord,
     &     frac_in_lrgest_clstr, frac_floppy)

      include      'set_parameter'
      integer      khb, hbond(3,maxh), hb_pick(maxh), locgroup(maxatm)
      integer      nres, gnum(maxgrp), f_so, f_sf, f_go, f_gf, hps
      integer      hb_number,hb_id(maxh),jhb
      real         mean_coord, frac_in_lrgest_clstr
      dimension    hb_energy(maxh)
      character*1  chainid(maxgrp),hb_type(maxh)
      character*4  aname(maxatm)
      character*80 decomp_file

      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag

      
 100  format('A:',1x,i7,1x,f9.5,1x,i7,1x,i7,a4,a4,1x,i7,a4,a4,4x,f5.3)
 101  format('B:',1x,f9.5)
 102  format('B:',1x,f9.6,2(2x,f11.9),2x,f9.5)
      open(91,file=decomp_file,access='append')

      f_so = hbond(1,hb_pick(jhb) )
      f_sf = hbond(3,hb_pick(jhb) )
      f_go = locgroup(f_so)
      f_gf = locgroup(f_sf)
c-------------------------------------------------------------------------------------
c AJR 10.19.01 now hb_number is the number of hbonds remaining after the jhb-th hbond
c is removed.  The E listed is the energy of the removed bond, jbh, and the <r> is for the
c skeleton network with that hydrogen bond removed.
c-----------------------------------------------------------------------------------
  
      hb_number = jhb -hps -1

      write(91,100) hb_number, hb_energy(hb_pick(jhb) ),
     &     nres, gnum(f_go), aname(f_so), chainid(f_go), 
     &     gnum(f_gf), aname(f_sf), chainid(f_gf), mean_coord
      write(91,102) mean_coord,frac_in_lrgest_clstr,frac_floppy,
     &     hb_energy(hb_pick(jhb))

      close(91)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Print a list of the hydrogen bonds according to the array
c     hb_pick.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine print_Hbond_list(khb, hbond, hb_pick, locgroup,
     &     aname, gnum )

      include      'set_parameter'      
      integer      i,khb,hbond(3,maxh), hb_pick(maxh), locgroup(maxatm)
      integer      gnum(maxgrp)
      character*4  aname(maxatm)

 100  format(i6,2x,a4,2x,i6,2x,a4,2x,i6)
      
      open(20, file='hbond.list',status='unknown' )
      do i = 1, khb
         f_so = hbond( 1, hb_pick(i) )
         f_sf = hbond( 3, hb_pick(i) )
         f_go = locgroup(f_so)
         f_gf = locgroup(f_sf)
         write(20,100) i, aname( f_so ), gnum( f_go ),
     &        aname( f_sf ), gnum( f_gf )
      enddo
      close(20)

      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
cccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine addthisHB(ihb,hbond,hb_pick)
      include      'set_parameter'
      integer      hbH,hbA,ns,hb_pick(maxh)
      integer      hpt,apt,hbond(3,maxh)
      integer      mult(maxatm),pointer(maxatm),link(nb2max)
      common/topology/ ns,link,pointer,mult
      
      jhb = hb_pick(ihb)
      hbH = hbond(2,jhb)
      hbA = hbond(3,jhb)
      
      hpt = pointer(hbH)
      apt = pointer(hbA)
      mult(hbH) = mult(hbH)+1
      mult(hbA) = mult(hbA)+1
      link(hpt+mult(hbH)) = hbA
      link(apt+mult(hbA)) = hbH

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     random function from Don Jacobs, it was used in all glass network
c     calculations.               AJR 10.09.01

      function ranx(i)
      i = i*1566083941 +1
      ranx = float(i)*2.328306e-10+0.5
      return
      end   
