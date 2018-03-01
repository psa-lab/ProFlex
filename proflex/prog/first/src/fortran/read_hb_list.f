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
      subroutine read_hb_list(fname)
      include      'set_parameter'
      integer      hbond(3,maxh),hb_id(maxh),hb_pick(maxh)
      dimension    hb_energy(maxh)
      character*1  hb_type(maxh)
      character*80 fname,criteria(99),line
c ===========================================================================
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/select/   nsc,criteria

      open(25,file=fname,status='unknown')
      rewind(25)
      khb = 0
      nhp_flag = 0
   10 read(25,*,end=11) itemp
c ---- 2006:03 SN --- 	 Here, make a note of any tethers being processed
      if (hb_type(itemp)(1:1).eq.'B') then
         nhp_flag = 1
      endif
      khb = khb + 1
      hb_pick(khb) = itemp
      goto 10
   11 continue
      close(25)

      open(26,file="qXyZaB.criteria",status='unknown')
      rewind(26)
      nsc = 0
   20 read(26,6000,end=21) line
 6000 format(a80)
      nsc = nsc + 1
      criteria(nsc) = line
      goto 20
   21 continue
      close(26)
      call system( '\\rm qXyZaB.criteria' )

      return
      end
