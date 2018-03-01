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
c ------------------------------------------------------------------------------
c     PROGRAM WRITTEN BY:                        Donald J. Jacobs   May 12, 1998
c ------------------------------------------------------------------------------
      subroutine group_max_dof(s0)
      include 'set_parameter'
      integer s0,tag
      integer pebble(6,0:maxatm),shell(0:maxatm),block(-1:maxatm)
      integer btrack(0:maxatm) 
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack
      nget = 0
      if( pebble(1,s0) .gt. 0 ) nget = nget + 1
      if( pebble(2,s0) .gt. 0 ) nget = nget + 1
      if( pebble(3,s0) .gt. 0 ) nget = nget + 1
      if( pebble(4,s0) .gt. 0 ) nget = nget + 1
      if( pebble(5,s0) .gt. 0 ) nget = nget + 1
      if( pebble(6,s0) .gt. 0 ) nget = nget + 1
      btrack(s0) = -1
         do iget=1,nget
         tag = tag + 1
         block(s0) = tag
         call collect1(s0)
         enddo
      return
      end
      subroutine free_max_dof(s0,s1,np)
      include 'set_parameter'
      integer s0,s1,tag
      integer pebble(6,0:maxatm),shell(0:maxatm),block(-1:maxatm)
      integer btrack(0:maxatm) 
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack
      call group_max_dof(s0)
      np = 0
      if( pebble(1,s1) .lt. 0 ) np = np + 1
      if( pebble(2,s1) .lt. 0 ) np = np + 1
      if( pebble(3,s1) .lt. 0 ) np = np + 1
      if( pebble(4,s1) .lt. 0 ) np = np + 1
      if( pebble(5,s1) .lt. 0 ) np = np + 1
      if( pebble(6,s1) .lt. 0 ) np = np + 1
      nget = 6 - np
      btrack(s1) = -1
         do iget=1,nget
         tag = tag + 1
         block(s0) = tag
         block(s1) = tag
         call collect1(s1)
            if( nsfail .lt. 0 ) then
            np = np + 1
            else
            return
            endif
         enddo
      return
      end
      subroutine cover_bonds(s0,s1,nb,iflop)
      include 'set_parameter'
      integer s,s0,s1,smin,tag
      integer pebble(6,0:maxatm),shell(0:maxatm),block(-1:maxatm)
      integer btrack(0:maxatm),laman(-1:maxatm) 
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack
      common/subgraph/ laman
      call group_max_dof(s0)
      nget = nb
      if( pebble(1,s1) .lt. 0 ) nget = nget - 1
      if( pebble(2,s1) .lt. 0 ) nget = nget - 1
      if( pebble(3,s1) .lt. 0 ) nget = nget - 1
      if( pebble(4,s1) .lt. 0 ) nget = nget - 1
      if( pebble(5,s1) .lt. 0 ) nget = nget - 1
      if( pebble(6,s1) .lt. 0 ) nget = nget - 1
      btrack(s1) = -1
         do iget=1,nget
         tag = tag + 1
         block(s0) = tag
         block(s1) = tag
         call collect1(s1)
            if( nsfail .gt. 0 ) then
            shell(0) = s0                      
            shell(1) = s1                             
            smin = s0
            iflop = iflop -        (nb-(nget-iget+1))
               if( nsfail .lt. nsfail_min ) then  
                  do k=0,nsfail
                  s = shell(k)
   50                if( laman(s) .lt. s ) then
                     s = laman(s)
                     goto 50
                     endif
                  if( s .lt. smin ) smin = s
                  enddo
               else
               k = -1
  100          k = k + 1
                  if( k .le. nsfail ) then
                  s = shell(k)
  150                if( laman(s) .lt. s ) then
                     s = laman(s)
                        if( block(s) .lt. tag ) then
                        call expand_laman(s)                  
                        endif
                     goto 150
                     endif
                  if( s .lt. smin ) smin = s
                  goto 100
                  endif
               endif
            is1 = laman(s0)                                                      
            laman(s0) = smin                                                     
               if( is1 .lt. s0 ) then                                            
  101          is0 = is1                                                         
               is1 = laman(is0)                                                  
               laman(is0) = smin                                                 
               if( is1 .lt. is0 ) goto 101                                       
               endif                                                             
               do k=1,nsfail
               s = shell(k)
               is1 = laman(s)                                                    
               laman(s) = smin                                                   
                  if( is1 .lt. s ) then                                          
  102             is0 = is1                                                      
                  is1 = laman(is0)                                               
                  laman(is0) = smin                                              
                  if( is1 .lt. is0 ) goto 102                                    
                  endif                                                          
               pebble(1,s) = s0
               pebble(2,s) = s0
               pebble(3,s) = s0
               pebble(4,s) = s0
               pebble(5,s) = s0
               pebble(6,s) = s0
               enddo
            return
            endif
         enddo
      iflop = iflop - nb
      np = nb - 1
         if( pebble(1,s1) .lt. 0 ) then
         np = np - 1
         pebble(1,s1) = s0                                               
         endif
      if( np .lt. 0 ) return
         if( pebble(2,s1) .lt. 0 ) then
         np = np - 1
         pebble(2,s1) = s0                                              
         endif
      if( np .lt. 0 ) return
         if( pebble(3,s1) .lt. 0 ) then
         np = np - 1
         pebble(3,s1) = s0                                               
         endif
      if( np .lt. 0 ) return
         if( pebble(4,s1) .lt. 0 ) then
         np = np - 1
         pebble(4,s1) = s0                                               
         endif
      if( np .lt. 0 ) return
         if( pebble(5,s1) .lt. 0 ) then
         np = np - 1
         pebble(5,s1) = s0                                               
         endif
      if( np .lt. 0 ) return
         if( pebble(6,s1) .lt. 0 ) then
         np = np - 1
         pebble(6,s1) = s0                                               
         endif
      return
      end
      subroutine lock_hinge(s0,s1,iflop)
      include   'set_parameter'
      integer   s,s0,s1,sa,tag,tc,tc_min
      integer   linknet(nb2max),ptr_net(maxatm)
      integer   multnet(maxatm)
      integer   pebble(6,0:maxatm),shell(0:maxatm),block(-1:maxatm)
      integer   btrack(0:maxatm)
      integer   link_hinge(nb2max),label_hinge(nbmax) 

      common/topology/ nsnet,linknet,ptr_net,multnet
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack
      common/cmotions/ label_hinge,link_hinge
      
      call group_max_dof(s0)
      nget = 1
      if( pebble(1,s1) .lt. 0 ) nget = nget - 1
      if( pebble(2,s1) .lt. 0 ) nget = nget - 1
      if( pebble(3,s1) .lt. 0 ) nget = nget - 1
      if( pebble(4,s1) .lt. 0 ) nget = nget - 1
      if( pebble(5,s1) .lt. 0 ) nget = nget - 1
      if( pebble(6,s1) .lt. 0 ) nget = nget - 1
      btrack(s1) = -1
         do iget=1,nget
         tag = tag + 1
         block(s0) = tag
         block(s1) = tag
         call collect1(s1)
            if( nsfail .gt. 0 ) then
            shell(0) = s0                          
            shell(1) = s1                  
            kmax = -1
            tc_min = iskiptag
               do k=0,nsfail
               so = shell(k)
               index = ptr_net(so)
                  do jo=1,multnet(so) 
                  index = index + 1
                  tc = link_hinge(index)
                     if( tc .gt. 0 ) then
                     sf = linknet(index)
                        if( block(sf) .eq. tag ) then
                        label = label_hinge(tc)
                           if( label .lt. 0 ) then
                          if( tc_min .gt. tc ) tc_min = tc
                           else
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c                                    patch modified by D. J. Jacobs Feb 08, 2002
c This patch makes sure that the minium collective motion label is found for 
c any other collective motion that intersects with the current failed pebble
c search.  A new function was created called  min_hinge_label()
c ------------------------------------------------------------------------------
                           label = min_hinge_label(label)
                           if( tc_min .gt. label ) tc_min = label
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           endif
                        endif
                     endif
                  enddo
               enddo
            label_hinge(tc_min) = tc_min
               do k=0,nsfail
               so = shell(k)
               index = ptr_net(so)
                  do jo=1,multnet(so) 
                  index = index + 1
                  tc = link_hinge(index)
                     if( tc .gt. 0 ) then
                     sf = linknet(index)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c                                    patch modified by D. J. Jacobs Feb 08, 2002
c This subroutine updates all other collective motions that intersect with the
c current failed pebble search. The previous version correctly identified the
c collective motion labels but did not update them properly. 
c ------------------------------------------------------------------------------
                        if( block(sf) .eq. tag ) then
                        label = label_hinge(tc)
                           if( label .lt. 0 ) then
                           label_hinge(tc) = tc_min
                           else
                           label_hinge(min_hinge_label(tc)) = tc_min
                           endif
                        endif
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     endif
                  enddo
               enddo
            return
            endif
         enddo
      iflop = iflop - 1
      np = 0
         if( pebble(1,s1) .lt. 0 ) then
         np = np - 1
         pebble(1,s1) = s0                                          
         endif
      if( np .lt. 0 ) return
c --------------------------
         if( pebble(2,s1) .lt. 0 ) then
         np = np - 1
         pebble(2,s1) = s0                                   
         endif
      if( np .lt. 0 ) return
c --------------------------
         if( pebble(3,s1) .lt. 0 ) then
         np = np - 1
         pebble(3,s1) = s0                                  
         endif
      if( np .lt. 0 ) return
c --------------------------
         if( pebble(4,s1) .lt. 0 ) then
         np = np - 1
         pebble(4,s1) = s0                              
         endif
      if( np .lt. 0 ) return
c --------------------------
         if( pebble(5,s1) .lt. 0 ) then
         np = np - 1
         pebble(5,s1) = s0                     
         endif
      if( np .lt. 0 ) return
c --------------------------
         if( pebble(6,s1) .lt. 0 ) then
         np = np - 1
         pebble(6,s1) = s0             
         endif
      return
      end

      function min_hinge_label(label_o)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c                                    patch modified by D. J. Jacobs Feb 08, 2002
c This subroutine was added in to make sure that the minimum collective motion
c labels are found.
c ------------------------------------------------------------------------------
c       Given an initial hinge label, the minimum hinge label will be produced.
c ==============================================================================
      include 'set_parameter'
      integer link_hinge(nb2max),label_hinge(nbmax)

      common/cmotions/ label_hinge,link_hinge
c ==============================================================================
c ------------------------------------------------ determine minimum hinge label
      min_hinge_label = label_o
   50    if( label_hinge(min_hinge_label) .eq. min_hinge_label ) return
         min_hinge_label = label_hinge(min_hinge_label)
         goto 50
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       subroutine collect1(s0)
      include 'set_parameter'
      integer s,s0,sa,sb,stest,tag
      integer pebble(6,0:maxatm),shell(0:maxatm),block(-1:maxatm)
      integer btrack(0:maxatm) 
      common/rigidity/ pebble,block,tag,shell,kmax
      common/search/   btrack
      kmax = 1      
      stest = pebble(1,s0)
         if( block(stest) .lt. tag ) then
         kmax = kmax + 1
         shell(kmax) = stest
         block(stest) = tag
         btrack(stest) = s0
         endif
      stest = pebble(2,s0)
         if( block(stest) .lt. tag ) then
         kmax = kmax + 1
         shell(kmax) = stest
         block(stest) = tag
         btrack(stest) = s0
         endif
      stest = pebble(3,s0)
         if( block(stest) .lt. tag ) then
         kmax = kmax + 1
         shell(kmax) = stest
         block(stest) = tag
         btrack(stest) = s0
         endif
      stest = pebble(4,s0)
         if( block(stest) .lt. tag ) then
         kmax = kmax + 1
         shell(kmax) = stest
         block(stest) = tag
         btrack(stest) = s0
         endif
      stest = pebble(5,s0)
         if( block(stest) .lt. tag ) then
         kmax = kmax + 1
         shell(kmax) = stest
         block(stest) = tag
         btrack(stest) = s0
         endif
      stest = pebble(6,s0)
         if( block(stest) .lt. tag ) then
         kmax = kmax + 1
         shell(kmax) = stest
         block(stest) = tag
         btrack(stest) = s0
         endif
      k = 1                                           
  100    k = k + 1
         if( k .gt. kmax ) return
         s = shell(k)
            do j=1,6
            stest = pebble(j,s)
               if( block(stest) .lt. tag ) then
               kmax = kmax + 1
               shell(kmax) = stest
               block(stest) = tag
               btrack(stest) = s
               else
                  if( stest .lt. 0 ) then
                  kmax = -1
                  sa = btrack(s)
                  pebble(j,s) = sa
  200                sb = btrack(sa)
                        if(     pebble(1,sa) .eq. s ) then
                        pebble(1,sa) = sb
                        elseif( pebble(2,sa) .eq. s ) then
                        pebble(2,sa) = sb
                        elseif( pebble(3,sa) .eq. s ) then
                        pebble(3,sa) = sb
                        elseif( pebble(4,sa) .eq. s ) then
                        pebble(4,sa) = sb
                        elseif( pebble(5,sa) .eq. s ) then
                        pebble(5,sa) = sb
                        else
                        pebble(6,sa) = sb
                        endif
                     if( sb .lt. 0 ) return
                     s  = btrack(sb)
                        if(     pebble(1,sb) .eq. sa ) then
                        pebble(1,sb) = s
                        elseif( pebble(2,sb) .eq. sa ) then
                        pebble(2,sb) = s
                        elseif( pebble(3,sb) .eq. sa ) then
                        pebble(3,sb) = s
                        elseif( pebble(4,sb) .eq. sa ) then
                        pebble(4,sb) = s
                        elseif( pebble(5,sb) .eq. sa ) then
                        pebble(5,sb) = s
                        else
                        pebble(6,sb) = s
                        endif
                     if( s  .lt. 0 ) return
                     sa = btrack(s)
                        if(     pebble(1,s) .eq. sb) then
                        pebble(1,s) = sa
                        elseif( pebble(2,s) .eq. sb) then
                        pebble(2,s) = sa
                        elseif( pebble(3,s) .eq. sb) then
                        pebble(3,s) = sa
                        elseif( pebble(4,s) .eq. sb) then
                        pebble(4,s) = sa
                        elseif( pebble(5,s) .eq. sb) then
                        pebble(5,s) = sa
                        else
                        pebble(6,s) = sa
                        endif
                     if( sa .lt. 0 ) return
                     goto 200
                  endif
               endif
            enddo
         goto 100
      end
      subroutine expand_laman(sf)
      include 'set_parameter'
      integer s,sf,stest,tag
      integer pebble(6,0:maxatm),shell(0:maxatm),block(-1:maxatm)
      common/rigidity/ pebble,block,tag,shell,kmax
      k = kmax
      kmax = kmax + 1
      shell(kmax) = sf
      block(sf) = tag
  100    k = k + 1
         if( k .gt. kmax ) return
         s = shell(k)
         stest = pebble(1,s)
            if( block(stest) .lt. tag ) then
            kmax = kmax + 1
            shell(kmax) = stest
            block(stest) = tag
            endif
         goto 100
      end
      function stress_label(so)
      include 'set_parameter'
      integer so,laman(-1:maxatm)
      common/subgraph/ laman
      stress_label = laman(so)
   50    if( laman(stress_label) .eq. stress_label ) return
         stress_label = laman(stress_label)
         goto 50
      end
      subroutine shorten_stress(so,smin)
      include 'set_parameter'
      integer s0,s1,so,smin,stress(-1:maxatm)
      common/subgraph/ stress 
      s1 = stress(so)
      stress(so) = smin
         if( s1 .lt. so ) then
  100       s0 = s1
            s1 = stress(s0)
            stress(s0) = smin
            if( s1 .lt. s0 ) goto 100
         endif
      return
      end
      subroutine update_laman(natom,laman)
      integer s,laman(-1:natom)
         do s=1,natom
         laman(s) = laman( laman(s) )
         enddo
      return
      end
      subroutine check_5dof(s0,s1,iflop)
      include 'set_parameter'
      integer s0,s1,tag
      integer pebble(6,0:maxatm),shell(0:maxatm),block(-1:maxatm)
      integer btrack(0:maxatm)
      common/rigidity/ pebble,block,tag,shell,nsfail
      common/search/   btrack
      call group_max_dof(s0)
      np = 0
      if( pebble(1,s1) .lt. 0 ) np = np + 1
      if( pebble(2,s1) .lt. 0 ) np = np + 1
      if( pebble(3,s1) .lt. 0 ) np = np + 1
      if( pebble(4,s1) .lt. 0 ) np = np + 1
      if( pebble(5,s1) .lt. 0 ) np = np + 1
      if( pebble(6,s1) .lt. 0 ) np = np + 1
      btrack(s1) = -1
      nget = 5 - np
         do iget=1,nget
         tag = tag + 1
         block(s0) = tag
         block(s1) = tag
         call collect1(s1)
            if( nsfail .lt. 0 ) then
            np = np + 1
            else
            iflop = iflop - np
            return
            endif
         enddo
      iflop = iflop - 5
      return
      end
