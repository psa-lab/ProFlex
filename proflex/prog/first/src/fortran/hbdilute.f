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


      subroutine hbdilute(iuseopt)
c                                             PROGRAM WRITTEN     March 18, 1999
c                                                    author: Brandon Hespenheide
c                                             Last modification:  Oct.  16, 2001 
c                                                         by AJ Rader
c                                             Last modification:  Apr.  21, 2008 
c                                                         by SK Namilikonda
c ------------------------------------------------------------------------------
c                              DESCRIPTION
c
c Options front end for the various routines using the hbdilute.c program
c ------------------------------------------------------------------------------
      include      'set_parameter'
      integer      tcmol,th_min,th_max,l, hps
      integer      locgroup(maxatm),slow(maxgrp),sbig(maxgrp)
      integer      gnum(maxgrp),hb_pick(maxh), temp_pick(maxh)
      integer      HBanalysis, itemp, fix_list, start
      integer      linknoh(nb2max),pointer(maxatm)
      integer      multnoh(maxatm), all_hbonds, this_bond
      integer      pmult(maxatm),plink(nb2max),ppt(maxatm)
      integer      maxmult(maxatm),gflag(maxgrp)
      integer      linkref(nb2max),point(maxatm), tally
      integer      multref(maxatm), break_it 
      integer      hbond(3,maxh),torsion(2,nbmax)
      integer      clst(maxatm),stress(-1:maxatm)
      integer      link_hinge(nb2max),label_hinge(nbmax)
      integer      output_format, increment_tally
      integer      no_new_info, HB_three, start_range, end_range
      integer      one, two, f_so, f_sf, f_go, f_gf, next_bond
      integer      hb_temp, switch, cluster_count,thb
      integer      iseed, winsize, window_prob
      integer      counter, g, a
      integer      number_of_runs, stored_values(maxh), on_off(maxh)
      integer      hb_id(maxh)
      integer      khb0, total_increments, hb
      integer      file_name_start, file_name_end
      real*4       initial_temp, final_temp, stepsize, temp_var
      real         mean_coord, factor, stop_point
      real         current_value, current_low_energy
      real         frac_lgst_clst,frac_floppy
      real         sum_of_f, sum_of_mean_coord
      real*4       cutoff, E_b, prefactor
      real*4       random, temperature
      double precision constant, P_open, test, randnum, partition_func
      double precision this_T_part_func
      dimension    sxyz(3,maxatm),freq(maxatm),bval(maxatm)
      dimension    hb_energy(maxh), hbond_exps(maxh)
      character*100 runline
      character*80 outfile(12),decomp_file,input_data,bondwt_file
      character*80 hphob_file
      character*5  idres(maxgrp)
      character*4  aname(maxatm),file_id
      character*4  donor, acceptor
      character*3  gname(maxgrp), name
      character*1  archive, hb_type(maxh), more_options, single_page
      character*1  rec(maxatm),chainid(maxgrp),remove_residues
      character*1  keep_mc,removeh,chain1,nowats,chain2,digit(0:9)
*
      LOGICAL      TMPOP

      common/fnames/   n_root,input_data,outfile,file_id,digit
      common/atomic/   natom,locgroup,rec,aname,sxyz,freq,bval
      common/groups/   ngrp,nres,nhet,nwater,chainid,
     &                 idres,gname,gnum,slow,sbig
      common/hbonds/   nhb,hbond,khb,hb_id,hb_pick,hb_energy,hb_type,
     &                 nhp,nhp_flag
      common/hbanalys/ HBanalysis,decomp_file,bondwt_file
      common/numbers2/ iflop,nhinge,nih,ncmode
      common/numbers1/ nflop,n_one,n_two,nbody,nclst,nstress
      common/network0/ nsnoh,linknoh,pointer,multnoh
      common/dihedral/ tcmol,th_min,th_max,torsion
      common/topology/ ns,linkref,point,multref
      common/ptopo/    pmult,plink,ppt,maxmult,gflag
      common/subgraph/ stress
      common/cmotions/ label_hinge,link_hinge
      common/clusters/ clst
      
c ==============================================================================
c      fraction_in_lrgest_clstr = 0.0

c ------------------------------------------------------------ format statements
 100  format('A:',1x,i7,1x,f9.5,1x,i7,1x,i7,a4,a4,1x,i7,a4,a4,4x,f5.3)
 2100 format('123456789-- Enter a 9-digit random seed.')
 2105 format(i9)
 6203 format(i5,1x,f9.5,1x,i5,1x,a4,1x,a1,1x,i5,1x,a4,1x,a1)
 7000 format(2x,'Which hydrogen bond dilution analysis would you like?')
 7001 format(3x,'(1) Standard hydrogen bond dilution, removing weakest',
     &          'H-bonds, one at a time')
c 7005 format(3x,'(2) Random dilution using a user defined window size.')
 7005 format(3x,'(2) Random dilution over all H-bonds.')
 7006 format(7x,'NOTE: This option is *not* recommended;')
 7007 format(7x,'      It may be used to probe the influence of H-bond')
 7008 format(7x,'      density, as opposed to strength, on rigidity')
 7018 format(3x,'(3) Flexible region and rigid cluster analysis at a')
 7019 format(3x,'    specific energy in the H-bond dilution.')
c 7020 format(3x,'4.  Thermal Dilution.')
 7010 format(i7)
 7200 format(5x,'Would like to see extra preprocessing options?')
 7201 format(a1)
 8100 format(2x,'Please enter an integer (0 - 32768) to seed')
 8101 format(2x,'the random number generator.')
 8105 format(i10)
 8110 format(2x,'Enter the size of the window from which to select')
 8111 format(2x,'the next Hbond(RETURN to select over all bonds)')
 8115 format(i4)
 8300 format(8x,'Would you like to remove H-bonds from the network')
 8301 format(a1)
 8303 format(8x,'based on H-bond number <array index>?')
 8302 format(8x,'Enter <y>es <n>o')
 8310 format(2x,'You must enter the number of the hydrogen bond as it')
 8311 format(2x,'occurs in the program. The hydrogens are sorted from')
 8312 format(2x,'weakest, which will be given a number corresponding')
 8313 format(2x,'to the total number of H-bonds in the protein, to the')
 8314 format(2x,'strongest, which will have a Hbond number of 1. You')
 8315 format(2x,'should probably run one of the other options to')
 8316 format(2x,'determine which H-bonds to break.')
 8317 format(2x,'Enter the number of the H-bond to remove. 0 to quit')
 8318 format(i7)
 8350 format(8x,'Would you like to keep all the main-chain bonds?')
 8351 format(8x,'If you answer yes, the main-chain to main-chain')
 8352 format(8x,'H-bonds will not be removed from the network.')
 8353 format(8x,'Enter <y>es <n>o')
 8355 format(a1)
 8356 format(2x,'here',1x,a4)
 8390 format(8x,'Would you like to remove H-bonds from a given')
 8391 format(8x,'residue range? <y>es or <n>o')
 8392 format(a1)
 8395 format(2x,'Would you like the output on a single page <y>or<n>?')
 8396 format(a1)
 8400 format(2x,'Enter the start of the range or zero to quit')
 8401 format(i7)
 8402 format(2x,'Enter the end of the range')
 8403 format(i7)
 8600 format(2x,'By what step size would you like to increment the',
     &     ' the temperature at each step of simulation')
 8601 format(2x,'during the simulation.')
 8602 format(f12.8)
 8603 format(2x,'How many runs would you like to average over at each', 
     &     ' temperature')
 8604 format(i6)
 8605 format(2x,'At what temperature would you like to begin the',
     &     ' simulation?')
 8606 format(f12.8)
 8607 format(2x,'At what temperature would you like to end the',
     &     ' simulation?')
 8608 format(f12.8)
 8609 format(2x,'Enter a g value (Boltzmann scaling factor):')
 8610 format(i6)
 8620 format(f8.6,2x,f8.6,2x,i5)
 8621 format(2x,f8.6,2x,f8.6,2x,f12.6)
c 9010 format(2x,'At which H-bond would you like to stop?')
c9010 format(2x,'Enter the H-bond number at which you would like to stop    
c    &(All h-bonds with energies lower than this will be kept and the an
c    &alysis run on this network) ')
 9010 format(2x,'Enter the H-bond energy at which you would like to', 
     &       ' perform')
 9020 format(2x,'analysis of flexible regions and rigid clusters' 
     &         ,' (in -Kcal/mol, e.g. -1.2): ')
c9010 format(2x,'Enter the number of hydrogen bond (in FIRSTdataset) at 
c    &which you would like to stop (All h-bonds stronger than this bond 
c    &will be included in the h-bond dilution) : ')

c --- 2006:03 Sandeep - Change the type of stop_point to real so that
c ---			user can input energy instead of H-bond number
 9011 format(f8.6)
 9500 format(a1)
c -------------------------------------------------------- automatic file naming

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc                   Print the options screen                 ccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c =================Sameer Arora 01.04.04 
c 'iuseopt' bypasses runtype options to allow noninteractive use of FIRST
c  noninteractive equivalent is "Standard Hydrogen Bond Dilution" of the 
c  interactive process. iuseopt should be 1 for non-interactive run

 150  write (*,*)
      if(iuseopt .eq. 1) then
	      HBanalysis = 1
      else
              call system('clear')
   	
	      write(6,*)
	      write(6,7000)
	      write(6,*)
	      write(6,7001)
	      write(6,*)
	      write(6,7005)
	      write(6,7006)
	      write(6,7007)
	      write(6,7008)
	      write(6,*)
*	      write(6,7018)
*	      write(6,7019)
*	      write(6,*)
c	      write(6,7020)
c	      write(6,*)

	      read(5,7010) HBanalysis
      endif      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c If no valid options are selected, bring up the options screen again.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      if( (HBanalysis .le. 0) .or. ( HBanalysis .gt. 2 ) ) then 
         goto 150
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     output files
      decomp_file(1:11) = "decomp_list"
      bondwt_file(1:8)  = "bond_wts"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Switch the Hbond arrays so that they are from highest energy to 
c lowest energy.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do switch = 1, nhb, 1
         temp_pick(switch) = hb_pick(switch)
      enddo

      do switch = 1, khb, 1
         hb_pick( switch ) = temp_pick( (khb-switch)+1 )         
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initialize Some Variables. Ingrained C programming
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mean_coord = 0.00
      frac_floppy = 0.00
      frac_lgst_clst = 0.00
      hps = 0
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     count how many hp-interactions there are.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do counter = 1, khb
         if( hb_energy(hb_pick(counter)) .le. -9.99999 ) then
c --- 2006:04 Sandeep	- Changed the cutoff from -9.9998 as some 
c ---			  salt bridges (SB) have energy 'e' such that
c ---			  -9.9998 <= e < -9.99999
            hps = hps + 1
         endif
      enddo

c      write(6,*) khb, ' Hydrogen bonds'
c      write(6,*) hps, ' Hydrophobic interactions.'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ask if user wants the output on a single page
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c =================Sameer Arora 01.04.04 
c 'iuseopt' bypasses runtype options to allow noninteractive use of FIRST
c  noninteractive mode should choose multi-page output      
c This question is appropriate only for HBanalysis=1, not for other
c options - so shifting it below

c      if(iuseopt .eq. 1) then
c	      single_page = 'N'
c      else
c	     write(6, 8395)
c	     read(5, 8396) single_page
c      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Option 1: incremental dilution of the H-bond network. FIRST is executed each 
c     time a bond is removed. Only changes in the decomposition are output.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if( HBanalysis .eq. 1 ) then


      if(iuseopt .eq. 1) then
	      single_page = 'N'
      else
	     write(6, 8395)
	     read(5, 8396) single_page
      endif
         khb0 = khb

         call print_decomp_list_header(nres,outfile,decomp_file,natom)
         call place_Hbond()
         call prunesetup(mean_coord,frac_lgst_clst,frac_floppy,0)
         call mapmolecule()
c call print_Hbond_info_to_decomp_file(khb0+1 -> 2007:06 SN --- out-of-bound array access!
         call print_Hbond_info_to_decomp_file(khb0, hps, locgroup,
     &        nres, gnum, aname, chainid, decomp_file, mean_coord,
     &        frac_lgst_clst, frac_floppy)
         call output_HBdilute(archive)


* SN 2008:05    Display a message to the user if non-interactive mode chosen!
         if (iuseopt.eq.1) then
             write(6,*)
             write(6,*) 'H-bond dilution in progress. Please wait... '
         end if


         do ihb = khb0-1, hps+1, -1
            call breakthis(ihb)

c1214 format ('H-bond: ',I6)
c     write(*,*) 
c     write(*,*) 'CK 1'
c     write (*,1214) ihb

            call hbprune(ihb,mean_coord,frac_lgst_clst,frac_floppy,0)

c     write(*,*) 'CK 2'

* SN 2008:02 	Disable this if non-interactive mode chosen!
            if (iuseopt.ne.1) then            
             write(6,*) '# hbonds present: ',ihb-hps
            end if
            call mapmolecule()
            call print_Hbond_info_to_decomp_file(ihb, hps, locgroup,
     &           nres, gnum, aname, chainid, decomp_file, mean_coord,
     &           frac_lgst_clst, frac_floppy)
            call output_HBdilute(archive)
         enddo

      INQUIRE(FILE='tmpout',OPENED=TMPOP)
      if(TMPOP) then
          call system('\\rm tmpout')
      endif      

c --- 2006:03 Sandeep	$FIRSTPTB is replaced with $PROFLEX_HOME
c         runline(1:44) = "$FIRSTPTB/hbdilute/bin/hbdilute decomp_list "
      runline(1:48) = "$PROFLEX_HOME/hbdilute/bin/hbdilute decomp_list "

c         file_name_start = 47
         file_name_start = 51

c --- 2006:03 Sandeep 	'_FIRSTdataset' replaced by '_proflexdataset'
c ---			Hence, change 13 to 15 below:
c         file_name_end   = 47 +n_root +15
         file_name_end   = 51 +n_root +15
c 
         runline(file_name_start:file_name_end) = input_data
         runline(file_name_end+1:100) = " "
         
         if((single_page .eq. 'Y').or.(single_page .eq. 'y')) then 
c            runline(45:46) = "b "
            runline(49:50) = "b "
            write(6,*)
            write(6,*) runline
            call system(runline)
         else
            runline(49:50) = "s "
            write(6,*)
            write(6,*) runline
            call system(runline)
         endif
ccc Remove only bond_wts file, not decomp_list- Sameer , Mar 11, 2004
c        call system('rm -f bond_wts decomp_list')
         call system('rm -f bond_wts')
c -----------------------------------------------------------------
c Sandeep's changes:    1) Open new h-phob file
c                       2) Write h-bond info into 58 and h-phob info into 57
c
         hphob_file = outfile(1)
         hphob_file(n_root+4:n_root+7) = 'phob'
         open(57,file=hphob_file,status='UNKNOWN')

         open(58,file=outfile(1),status='UNKNOWN')
         rewind(58)
c ----------------------------------------------------------------- write H-bond
c --- At this point, hblist is in reverse order of ids
c --- Write hphobs
         if(nhp_flag .eq. 0) then
c --- 2006:04	- Sandeep	If (filter tethers option was chosen) THEN
            do hb=khb0,1,-1
               write(58,*) hb_id( hb_pick(hb) )
            enddo
         else
            do hb=nhp,1,-1
               write(57,*) hb_id( hb_pick(hb) )
            enddo
c --- Now write h-bonds
            do hb=khb0,nhp+1,-1
               write(58,*) hb_id( hb_pick(hb) )
            enddo
         endif
c -----------------------------------------------------------------
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Option 2: Allows for random selection of the next hydrogen bond to remove
c     The random seed, and window size are input by the user. The random number
c     code is architecture dependent, so be sure to check the code when compiling
c     on a new maching. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
      if( HBanalysis .eq. 2 ) then

         khb0 = khb
         write(6,2100)
         read(5,2105) iseed
cc Sameer 12.Mar.04 - the only time this option is random is when window includes
cc all bonds. So why bother asking for window size
c        write(6,8110)
c        write(6,8111)
c        read(5,8115) win_siz
c        window_size = win_siz
         window_size = 0
         
c         khb  = khb0	- 2006:03 Sandeep	- REDUNDANT STATEMENT
         
         call print_decomp_list_header(nres,outfile,decomp_file,natom) 
         call place_Hbond()
         call prunesetup(mean_coord,frac_lgst_clst,frac_floppy,1)
         call mapmolecule()
         call print_Hbond_info_to_decomp_file(khb+1, hps, locgroup,
     &        nres, gnum, aname, chainid, decomp_file, mean_coord,
     &        frac_lgst_clst, frac_floppy)
         
         call output_HBdilute(archive)
         
         if( window_size .eq. 0 ) then
            window_size = khb
            write(6,*) 'Random selection over all Hbonds'
         else
            write(6,*) 'Random selection over a window of',
     &           window_size, ' bonds.'
         endif
         
         do khb = khb0, hps+1, -1
            
            random = ranx(iseed)
            if(random .eq. 1.0 ) random = 0.999999
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     When the number of hbonds remaining is less than the window size, set
c     the window size to be the number of hbonds remaining. 
            if( (khb - hps) .lt. window_size ) then
               window_size = khb-hps
            endif
            
            ihb = int( window_size * random )
            next_bond = khb - ihb
            
            hb_temp = hb_pick( khb )
            hb_pick( khb ) = hb_pick( next_bond )
            hb_pick( next_bond ) = hb_temp
            
            write(6,*) 'index ', khb-hps, ' chosen ', next_bond-hps
            
            call breakthis(khb)
            call hbprune(khb,mean_coord,frac_lgst_clst,frac_floppy,0)
            call mapmolecule()
            call print_Hbond_info_to_decomp_file(khb, hps, locgroup,
     &           nres, gnum, aname, chainid, decomp_file, mean_coord,
     &           frac_lgst_clst, frac_floppy)
            call output_HBdilute(archive)
         enddo

c --- 2006:03 Sandeep   $FIRSTPTB is replaced with $PROFLEX_HOME
      runline(1:48) = "$PROFLEX_HOME/hbdilute/bin/hbdilute decomp_list "
         file_name_start = 51
         file_name_end   = 51 +n_root +15
         runline(file_name_start:file_name_end) = input_data
         runline(file_name_end+1:100) = " "

         if((single_page .eq. 'Y').or.(single_page .eq. 'y')) then
            runline(49:50) = "b "
            write(6,*) runline
            call system(runline)
         else
            runline(49:50) = "s "
            write(6,*) runline
            call system(runline)
         endif
         
ccc Remove only bond_wts file, not decomp_list- Sameer , Mar 11, 2004
c        call system('rm -f bond_wts decomp_list')
         call system('rm -f bond_wts')
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Option 3: This option will output the regular FIRST file for
c     creating the pretty graphics files. Only works with the standard
c     dilution data. You can give it a hydrogen bond which will act
c     as a cutoff value. All h-bonds With energies lower than this
c     one will be kept, and the analysis run on this network. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if( HBanalysis .eq. 3 ) then

         write(6,9010)
         write(6,9020)
         read(5,9011) stop_point

c --- 2006:03 SN	- By now, the h-bond list is stored in the reverse order
c ---			  of their ids. So, hb_pick(1) ----> last tether atom #
c ---			  hb_pick(nhp+1)-> last h-bond in _proflexdataset & so on
c ---			- khb currently will be equal to (nhp + # of filtered hbnd)
         do counter = nhp+1, khb
            if( hb_energy(hb_pick(counter)) >=  stop_point ) then
               goto 1000
            endif
         enddo

1000     if (counter < khb) then
            khb = counter-1
c --- 2006:03 SN	- Above, khb has to be carefully adjusted so as to
c ---			  produce the right output in outputfirst.f
c ---			- Note that counter corresponds to the first H-bond
c ---			  encountered that violates the energy cutoff and hence,
c ---			  has to be discarded too
         endif

         call place_Hbond()
         call mapmolecule()

c --- 2006:03 Sandeep	- The hb_id() is in decreasing order of h-bond ids
c ---			  Hence, set a flag's value accordingly so that
c ---			  outputfirst() can write out info in the right order
         nhp_flag = 2

         call outputfirst(archive,iuseopt)
         
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Option 4: Thermal Dilution Routine. 
c
c     This routine uses a function to determine the probability that a
c     hydrogen bond exists based on its energy and the current temp. 
c     At a given temperature, the routine will cycle over "number_of_runs"
c     many FIRST analyses. The results are averaged. The users sets 
c     the initial and final temperatures, and the temperature 
c     increment. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
      if( HBanalysis .eq. 4 ) then

         constant = 503.467961

         write(6,2100)
         read(5,2105) iseed
         write(6,8600)
c        write(6,8601)
         read(5,8602) stepsize
         write(6,8603)
         read(5,8604) number_of_runs
         write(6,8605)
         read(5,8606) initial_temp
         write(6,8607)
         read(5,8608) final_temp
         write(6,8609)
         read(5,8610) g

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Store the original values of the hb_pick array so that we can mess 
c  them up stochastically, and have a reference to the original array
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         all_hbonds = khb
         call srand( iseed )               

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     store the original bond network so it can be 
c     reset later. There are 'hps' hydrophobic 
c     interactions, so start accessing the hbonds
c     at hps+1. 
cccccccccccccccccccccccccccccccccccccccccccccccccc
         do counter = hps+1, all_hbonds
            stored_values(counter) =  hb_pick(counter)
         enddo

         open( 40, file="f_VS_mean_coord",status='new')
         write(40,*) initial_temp, final_temp, stepsize,
     &        number_of_runs, g

         temp_var = (final_temp-initial_temp)/stepsize
         total_increments = int(temp_var)

         write(40,*) final_temp, initial_temp, stepsize, temp_var,
     &        total_increments
         write(6,*) total_increments
         total_increment = total_increments + 1
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Loop over the temperature range in increments
c     of "stepsize". 
cccccccccccccccccccccccccccccccccccccccccccccccccc
         do increment_tally = 0, total_increments

            temperature = initial_temp + (increment_tally*stepsize)

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Prepare output files. c
cccccccccccccccccccccccccccccccccccccccccccccccccc
            write(6,*)'Current Temp', temperature
            
            open(23,file=bondwt_file,status='new',
     &           access='sequential')
            write(23,*)'INITIAL', initial_temp
            write(23,*)'FINAL', final_temp 
            write(23,*)'GVALUE', g
            write(23,*)'DELTATEMP',stepsize
            write(23,*)'TEMP', temperature
            close(23)

ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Do "number_of_runs" many FIRST analyses. 
cccccccccccccccccccccccccccccccccccccccccccccccccc
            do increment = 1, number_of_runs

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Reset the array of which bonds are selected 
c     to be there, to all 'on'.
cccccccccccccccccccccccccccccccccccccccccccccccccc
               do counter = hps+1, all_hbonds
                  on_off(counter) = 1
               enddo

cccccccccccccccccccccccccccccccccccccccccccccccccc               
c     Stochastically select which bonds 'exist' 
c     according to temp and energy. 
cccccccccccccccccccccccccccccccccccccccccccccccccc
               do a = hps+1, all_hbonds
                  
                  E_b = 0 - hb_energy(stored_values(a))
                  
                  factor = E_b / temperature
                  if( factor .ge. 0.1 ) then
                     P_open = 1.0
                     goto 134
                  endif
                  exponential = exp( factor * constant )
                  P_open = exponential / ( g + exponential )
                  

c     This function returns a R.N. between 0 and 1.0
 134              random = ranx(iseed)

c                  write(6,*) P_open, random, factor, 
c     &                 hb_energy(stored_values(a))

                  if( P_open .le. random ) then
                     on_off(a) = 0
                  elseif( P_open .gt. random ) then
                     on_off(a) = 1
                  endif
                  
               enddo

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Modify the hb_pick array to include only 
c     those bonds that 'exist'.
cccccccccccccccccccccccccccccccccccccccccccccccccc
               tally = hps+1
               khb = all_hbonds
               do counter = hps+1, all_hbonds
                  if(on_off(counter) .eq. 1 ) then
                     hb_pick(tally) = stored_values(counter)
                     tally = tally + 1
                  endif
                  if(on_off(counter) .eq. 0 ) then
                     khb = khb - 1
                  endif
               enddo
            
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Error check. 
cccccccccccccccccccccccccccccccccccccccccccccccccc
               if( (tally-1) .ne. khb ) then
                  write(6,*) 'there is an error in the program.'
               endif

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Run the pebble game on the network. 
cccccccccccccccccccccccccccccccccccccccccccccccccc
               call place_Hbond()
               call prunesetup(mean_coord,frac_lgst_clst,frac_floppy,0)
               call mapmolecule()
               call output_bondwts(archive)

c               write(6,*) 'TEMP ', temperature, ' <r> ', mean_coord,
c     &              ' f ', frac_floppy
               sum_of_f          = sum_of_f + frac_floppy
               sum_of_mean_coord = sum_of_mean_coord + mean_coord
c               write(40,8620) frac_floppy,mean_coord,tally-hps-1
               
            enddo   
            
            open(23,file=bondwt_file,access='append')
            write(23,*)'STOP'
            close(23)
            
c            write(40,*)
            write(40,8621) sum_of_f/number_of_runs, 
     &           sum_of_mean_coord/number_of_runs, temperature
c            write(40,*)
c            write(40,*)

            sum_of_f          = 0.0
            sum_of_mean_coord = 0.0
            
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     Create the 'stripy plot' line for the 
c     current temperature and delete temp files. 
cccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc flex_var2 was a temporary program ( Brandon informed in email)
ccccccccc Commenting this call - Sameer Arora, Mar 6, 2004
c           call system('$FIRSTPTB/hbdilute/bin/flex_var2 bond_wts
c    &           STdataset')
            call system('\\rm bond_wts')
            
         enddo   
         
         close(40)
         call system('\\rm old_colors black_line line ')
         
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
 3000 return
      end


