/*******************************************************************************
*  MSU ProFlex, formerly called FIRST, is a software developed to predict and  *
*  analyze protein flexibility.                                                *
*  This source file is a part of MSU ProFlex.                                  *
*                                                                              *
*  Copyright (C) 1997 - 2008, Michigan State University.                       *
*                                                                              *
*  This program is free software; you can redistribute to academic users only, *
*  it and/or modify it under the terms of the GNU General Public License,      *
*  version 2, as published by the Free Software Foundation.                    *
*                                                                              *
*  This program is distributed in the hope that it will be useful,             *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of              *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
*  GNU General Public License for more details.                                *
*                                                                              *
*  You should have received a copy of the GNU General Public License           *
*  along with this program; if not, write to the Free Software                 *
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA,  *
*  or see http://www.gnu.org/licenses/gpl.txt                                  *
*******************************************************************************/
/* THIS IS THE C++ CODE FOR pick_Hbonds() subroutine                  */
/*  2006:01:26  --- Sandeep (SN)   (WRITTEN AFRESH RETAINING RELEVANT PREVIOUS CODE)

    usage >= 1 --> Keep all H-bonds and filter by H-bond cutoff energy

    Changes to be implemented:
    --------------------------
	- First, filter the H-bond list through tight geometric rules for
	  H-A, D-A distances, D-H-A and H-A-P angles. Further, user can 
	  specify more stringent ruleset.
        - Filtering options are made more hierarchical i.e., first start with
          all H-bonds, ask the user if Water molecules have to removed, then
          the side chain H-bonds, and finally, HETATM H-bonds (aside from water)
        - Next, using the filtered H-bonds list, ask the user for an energy cutoff, 
	  which is then used as an additional energy-based filter
    STEPS:
    ------
        1) Disable all unused options by removing the correponding code
        2) Make changes to some global variables decalred in "class.h"
        3) Change the initial screening options
	4) Impose stronger stereochemical selection criteria to select
	   the initial "DEFAULT" list of H-bonds
	5) Provide options for further stronger geometric filter
	6) Provide filters based on the type of bond (HET, tether, etc.)
	7) Output proflexdataset
	8) Energy-based filtering with the user-input or default cutoff
	9) Output allbonds_SEfilt
	10) Output summary file that documents all the above changes incorporated

	2006:03:15	--- added extra p_flag to facilitate stringent stereo
			    chemical filtering by user when '-p' flag is selected

	2006:04:10	--- Changes to non-H2O HETatm filter
			--- sidechain-non-water HETATM filter added 
			    This matters when the user only selects to remove 
			    nonwater-HETATM H-bonds as this set of bonds were
			    being left out unscreened.
			    Further,	
			    : the atom name field format is as follows:
				'[][][atom-tye][position][label]' ex: '  NH1','  CD '
			--- So, the filter should just check for the atom-type and not
			    the whole atom name
			
*/
/*  --- 	modifications on 04.17.01 by AJR ---
    I remedied the usagse of screening options 11,16,and 17.  Also
    attempting to add an option:  input file with desired h-bond list.

    AJR 03.25.02
    usage == 1 --> -non -h runoptions 
    usage == 2 --> -non -p runoptions
    usage == 0 --> otherwise

	-non => non-interactive H-bond dilution => default energy filter
*/


/**********************************************************************/
#include <cstring>

void list::pick_hbonds(int usage, int p_flag ) {
//---------------------------------------------------------Declare variables
  float xda,xha,adha,ahap; //DONOR-H-ACPTR-PREAC distance and angles
  float s_xda,sb_xda,s_xha,sb_xha;
  float min_theta,hb_ecut,energy_max,hd2,ha2,ad2,costheta,ha2max,ad2max;
  float hp2,ap2,cosdelta,min_delta,min_cosdelta,min_costheta; // pre-acptr related
  int i,j,khb,khb_min,khb_max,khb2;
  int flag=0;
  int so,s,sf,nsc,sp; // sp -- preacptr
  char my_answer,ans;
  float energy_input;
  ofstream fptr;
  ifstream pacc_ifil; //--- 2006:03 SN

  /************************************************************/
  /* Allocate memory for hb_pick[]                            */
  /************************************************************/
 if(p_flag == 1)       	//--- 2006:03 SN --- *COMPLICATION*  When '-p' option is
                        // chosen, real_hb_no gives the actual HB+SB count and
                        // nhb is the cumulative of tether and real_hb_no
    {
        khb = real_hb_no;
        no_tether = nhb - khb;
 	nhb = khb;
    }
 else			//--- 2006:03 SN --- '-h' or options other than '-p'
  khb = nhb;

  hb_pick = new int[nhb+1];

// --- SN 2006:02       --- Initialize hb_pick[]
  for(i=0;i<=nhb;i++)
     hb_pick[i]=1;

 //------------- DEFAULT OPTION:-Keep all hydrogen bonds => khb = nhb
 // This variable is now decremented once per bond whose hb_pick[] is made 0
 // ------- SN 2006:03:01


//-----------------------------------------Initialize assumed geometric criteria
//----------------------------------------- read from the ../../include/class.h
  xda=admaxd;
  xha=hamaxd;
  adha=dha_ang;
  //2006:01     --- added by SN
  ahap=hap_ang;
 
  //Delta and Theta anglular constraints are assumed to be the same for SB and
  //Sulfur-based H-bonds. So, I only implement the distance-based constraints
  s_xda=s_admaxd;  
  sb_xda=sb_admaxd;
  s_xha=s_hamaxd;
  sb_xha=sb_hamaxd;

  // If user chooses interactive filtering then
  if(usage == 0)
  {
    //----------------------------------------------------------
    //--- Ask user if stronger geometric constraints are to be imposed
    //--- If so, read the constraints and impose the constraints
   do{
    system("clear");
    cout<<"\t\tH-BOND/SALT BRIDGE BOND LENGTH & ANGLE CRITERIA"<<endl<<endl;
    cout<<"Default H-bond & salt bridge list has been compiled and written to"<<endl;
    cout<<"'_proflexdataset' file. It includes intra- and inter-molecular"<<endl;
    cout<<"H-bonds meeting the following default stereochemical criteria:"<<endl;
    cout<<"\n\n\tHydrogen-Acceptor (H-A) distance <= 2.5 Angstroms";
      cout<<"\n\t                with Sulfur atom <= 3.0 Angstroms";
      cout<<"\n\t                for Salt bridge  <= 3.5 Angstroms"<<endl;
      cout<<"\n\tDonor-Acceptor (D-A) distance    <= 3.5 Angstroms";
      cout<<"\n\t             with Sulfur atom    <= 4.0 Angstroms";
      cout<<"\n\t              for Salt bridge    <= 4.5 Angstroms"<<endl;
      cout<<"\n\tTheta (D-H-A angle)              >= 110 Degrees";
      cout<<"\n\tDelta (H-A-Pre-acceptor angle)   >= 90  Degrees";
      cout<<"\n\nIndices of H-bonds filtered more stringently (see below) are written to"; 
      cout<<"\n'_h-bonds_SEfilt' file."; 
      cout<<"\n\nDo you wish to impose more stringent rules? (y/n): ";
    cin>>my_answer;

   }while(my_answer != 'n' && my_answer != 'N' && my_answer != 'y' &&
          my_answer != 'Y');

  } // End interactive mode
  else
  {
    my_answer = 'n';
  }
 		 
  if(my_answer == 'y' || my_answer == 'Y')
    {
      system("clear");
      
      //------------------------------------------------------------Read xha
      cout<<"\n\tEnter H-A Distance (<= 2.5 A): ";
      cin>>xha;
      
      while(xha > hamaxd || xha < 0)
	{
	  cout<<endl<<"\t\tInvalid value; must be <= 2.5 A";
          cout<<endl<<"\t\tPlease enter another value: ";
	  cin>>xha;
	}

      //------------------------------------------------------------Read s_xha
      cout<<"\n\tEnter H-A Distance for H-bonds with Sulfur (<= 3.0 A): ";
      cin>>s_xha;

      while(s_xha > s_hamaxd || s_xha < 0)
        {
          cout<<endl<<"\t\tInvalid value; must be <= 3.0 A";
          cout<<endl<<"\t\tPlease enter another value: ";
          cin>>s_xha;
        }

      //------------------------------------------------------------Read sb_xha
      cout<<"\n\tEnter H-A Distance for salt bridges (<= 3.5 A): ";
      cin>>sb_xha;

      while(sb_xha > sb_hamaxd || sb_xha < 0)
        {
          cout<<endl<<"\t\tInvalid value; Must be <= 3.5 A";
          cout<<endl<<"\t\tPlease enter another value: ";
          cin>>sb_xha;
        }

      //------------------------------------------------------------Read xda
      cout<<"\n\tEnter D-A Distance (<= 3.5): ";
      cin>>xda;
      
      while(xda > admaxd || xda < 0)
	{
	  cout<<endl<<"\t\tInvalid value; <= 3.5 A";
          cout<<endl<<"\t\tPlease enter another value: ";
	  cin>>xda;
	}

      //------------------------------------------------------------Read s_xda
      cout<<"\n\tEnter D-A Distance for H-bonds with Sulfur (<= 4.0): ";
      cin>>s_xda;

      while(s_xda > s_admaxd || s_xda < 0)
        {
          cout<<endl<<"\t\tInvalid value; must be <= 4.0 A";
          cout<<endl<<"\t\tPlease enter another value: ";
          cin>>s_xda;
        }

      //------------------------------------------------------------Read sb_xda
      cout<<"\n\tEnter D-A Distance for salt bridges (<= 4.5): ";
      cin>>sb_xda;

      while(sb_xda > sb_admaxd || sb_xda < 0)
        {
          cout<<endl<<"\t\tInvalid value; must be <= 4.5 A";
          cout<<endl<<"\t\tPlease enter another value: ";
          cin>>sb_xda;
        }


      //------------------------------------------------------------Read theta	
      cout<<"\n\tEnter Theta Angle (>= 110): ";
      cin>>adha;
      adha = adha*raddeg;
      /* 
       * 2006:01	-- Sandeep
       * Ensure that input angle < 180 and > THETA_limit and theta input by user
       */ 
      while( (adha > 3.14159) || (adha < dha_ang) )
	{
	  cout<<endl<<"\tInvalid value; must be >= 110 "<<endl;
          cout<<endl<<"\t\tPlease enter another value: ";
	  cin>>adha;
	  adha=adha*raddeg;
	}        
      //------------------------------------------------------------Read delta
      cout<<"\n\tEnter Delta Angle (>= 90): ";
      cin>>ahap;
      ahap = ahap*raddeg;
      
      while( (ahap > 3.14159) || (ahap < hap_ang) )
	{
	  cout<<endl<<"\tInvalid value; must be >= 90"<<endl;
          cout<<endl<<"\t\tPlease enter another value: ";
	  cin>>ahap;
	  ahap=ahap*raddeg;
	}
      
      //------------------------------------------------------------Impose new filters
      
      // --- SN 2006:02	--- Initial Geometric filtering is performed in
      //			    find_hbonds.cpp 
      // 			--- Perform new filtering only if more stringent 
      //			--- geometric constraints have been set
      
      // --- SN 2006:03	--- IF p_flag is set THEN
      //			    READ preacptr info from "preacptr_info" file

      if(p_flag == 1)
	pacc_ifil.open("preacptr_info",ios::in);
      
      for(i=1;i<=nhb;i++)
	{
	if(hb_pick[i] == 1)
	 {
          /* --- SN 2007:02	  --- hbond[i][0] is loaded in
				      find_hbonds.cpp. So, we need not
				      read the preacptr_info file

	  // --- SN 2006:03       --- IF p_flag is set THEN
          //                          READ preacptr info from "preacptr_info" file
	  if(p_flag == 1 && !pacc_ifil.eof())
	    {
	      pacc_ifil.getline(line,10);
	      strncpy(pacc_val,line,5);
	      hbond[i][0]=atoi(pacc_val);
	    }
          */

	  /************************************************************/
	  /* Perform stereochemical filtering to obtain the initial   */
	  /* list of H-bonds                                          */
	  /************************************************************/
	  sp=hbond[i][0];
	  so=hbond[i][1];
	  s=hbond[i][2];
	  sf=hbond[i][3];

#ifdef CHKPT
	  cout<<"\n"<<i<<"sp = "<<sp<<"\tsf = "<<sf<<"\ts = "<<s<<"\tso = "<<so;
#endif

	  //----------------------------------------------------Determine geometrical case
	  min_theta = adha;
	  min_delta = ahap;
	  
	  hd2=0;
	  ha2=0;
	  ad2=0;
	  hp2=0;
	  ap2=0;
	  //-----------------------Calculate H-donor, H-acceptor, acceptor-donor distances
	  for(j=1;j<=3;j++)
	    {
	      hd2+=pow((p[so]->r1.coord[j]-p[s]->r1.coord[j]),2);
	      ha2+=pow((p[sf]->r1.coord[j]-p[s]->r1.coord[j]),2);
	      ad2+=pow((p[sf]->r1.coord[j]-p[so]->r1.coord[j]),2);
              /* 2007:02 SN */
	      if(sp != 0) {
	        hp2+=pow((p[s]->r1.coord[j]-p[sp]->r1.coord[j]),2);
	        ap2+=pow((p[sf]->r1.coord[j]-p[sp]->r1.coord[j]),2);
	      }	
	    }
	  //----------------------------------------------------------Calculate cos(angle)
	  costheta = (hd2+ha2-ad2)/sqrt(4*hd2*ha2);
	  min_costheta = cos(min_theta);
          /* 2007:02 SN */
          if(sp != 0) {
	    cosdelta = (ha2+ap2-hp2)/sqrt(4*ha2*ap2);
	    min_cosdelta = cos(min_delta);
          }
	  //----------------------------------------------------------Geometric filters
                  if(hb_type[i]<4)
                    {
                      ha2max = xha*xha;
                      ad2max = xda*xda;
                    }
                  else if(hb_type[i]==4)
                    {
                      ha2max = sb_xha*sb_xha;
                      ad2max = sb_xda*sb_xda;
                    }
                  else if(hb_type[i]==5)
                    {
                      ha2max = sb_xha*sb_xha;
                      ad2max = sb_xda*sb_xda;
                    }
                  else
                    {
                      ha2max = s_xha*s_xha;
                      ad2max = s_xda*s_xda;
                    }

/*
 * 2007:06:07   SN      --- Salt bridges do not have to be subjected
 *                          to angular check and can be accepted if H-A or
 *                          D-A distance is acceptable
 */
      if((hb_type[i] == 4 && ha2 > ha2max && ad2 > ad2max) || \
	 (hb_type[i] != 4 && (ad2 > ad2max || ha2 > ha2max || \
	   costheta > min_costheta || (sp != 0 && cosdelta > min_cosdelta))))
					/* 2007:02 SN --- Ensure existence of 
					   preacptr before testing for delta */
	    {
	      khb--;
	      hb_pick[i] = 0;
	    }
	 } // --- End if(hb_pick(i) == 1) 
	} // --- End for(i) 
      /* ------------------- Geometric filter implementation ends ---------- */

#ifdef CHKPT
      cout<<"\n---------------- DEBUG END ---------------\n";
#endif

      if(p_flag == 1)
	pacc_ifil.close();

    } // End if(my_answer == 'y' ||)

  /*
   * END OF STEREO-CHEMICAL FILTERING: So, display the HB and HP counts
   */
    if(usage == 0)	// interactive mode
    {
      cout<<"\n\tNumber of H-bonds remaining: "<<khb;
      cout<<"\n\tNumber of hydrophobic tethers remaining: "<<no_tether<<endl;
      cout<<"\n\n\nPlease enter 'c' and press enter to continue:";
      cin>>my_answer;
    }


  /************************************************************/
  /* Define and accept screen options for identifying Hbonds  */
  /************************************************************/
  while( flag == 0 ) {
    if(usage == 0) {
      pick_hbonds_options();
      
      cout << endl << endl << endl << "\t\tSELECTED SCREENING OPTIONS" << endl << endl;
      
      //--- 2006:03 SN - If default (option '5' or choice[1]) is selected then make
      //---		 other choices zero and proceed as if 'f' were selected 
      if(choice[1] == 1)
	{
	  cout<<"\t"<<setw(2) <<maxnopt<< " "<<screen_msgs[1]<<endl;
	  for( i = 2; i <= maxnopt; i++ )
	     choice[i] = 0;
	}
      else
	{
          for( i = 2; i <= maxnopt; i++ ) {
   	    if(choice[i]==1) {
	      cout<<"\t"<<setw(2) <<i-1<< " "<<screen_msgs[i]<<endl;
	    }
          }
	}

      cin.ignore(80,'\n');
      flag = 1;
      ans = '@';
      while( ans != 'Y' &&
	     ans != 'y' &&
	     ans != 'N' &&
	     ans != 'n') {

	cout << endl << "\tDo you accept these options? (y)es or (no): ";
	ans = cin.get();

	if( ans == '\n' )
	  cout << "\t\tPlease enter \"y\" or \"n\"." << endl;
	else{
	  if( ans == 'n' ||
	      ans == 'N' )
	    flag = 0;
	  cin.ignore( 80, '\n' );
	}
      }
    }
    else if(usage >= 1)	{
      for( i = 2; i <= maxnopt; i++ ) {
	choice[i]=0;
      }
      choice[1] = 1;
      flag =1;
    }
  }


  for(i = 1; i <= nhb; i++)
    {
     if(hb_pick[i] == 1)
     {
      so = hbond[i][1];
      sf = hbond[i][3];

      /************************************************************/
      /* Implement the filters chosen above:                      */
      /************************************************************/

      /* -------- SN 2006:03 ------------
	 Removing *ALL* water-bound H-bonds is default. So, the below filter
	 has to be implemented unless the user chooses not to by selecting choice[2]
      */
      if(choice[2] == 0) 
	//----------------------------------------- Remove *ALL* H-bonds involving Water
	//----------------------------------------- Old choices: 5,8, and 10
	{
	  //--------------------------------------- 5)   mainchain-WATERatom H-bonds
	  //--------------------------------------- 8)   sidechain-WATERatom H-bonds
	  //--------------------------------------- 10)  WATERatom-WATERatom H-bonds

          if(( p[so]->r1.field1[0]!='H' && strcmp(p[sf]->r1.rname,"HOH")==0 &&
               (strcmp(p[so]->r1.aname,"  N  ")==0
                || strcmp(p[so]->r1.aname,"  O  ")==0)
              )
             ||(p[sf]->r1.field1[0]!='H' && strcmp(p[so]->r1.rname,"HOH")==0 &&
                (strcmp(p[sf]->r1.aname,"  N  ")==0
                 || strcmp(p[sf]->r1.aname,"  O  ")==0)
                )
             )
	    {
	      khb--;
	      hb_pick[i] = 0;
	    }
	  else
            if(((strcmp(p[sf]->r1.rname,"HOH")==0 && p[so]->r1.field1[0]!='H') &&
                (strcmp(p[so]->r1.aname,"  N  ")!=0 && strcmp(p[so]->r1.aname,"  O  ")!=0))
               ||((strcmp(p[so]->r1.rname,"HOH")==0 && p[sf]->r1.field1[0]!='H')&&
                  (strcmp(p[sf]->r1.aname,"  N  ")!=0 && strcmp(p[sf]->r1.aname,"  O  ")!=0))
               )
	      {
		khb--;
		hb_pick[i]=0;
	      }
	    else
	      if(strcmp(p[sf]->r1.rname,"HOH")==0 && strcmp(p[so]->r1.rname,"HOH")==0)
		{
		  hb_pick[i]=0;
		  khb--;
		}
	}
      
      //------------ Remove *ALL* H-bonds involving side-chains
      //------------ Old choices: 3,6, and 7
      if(choice[3] == 1) 
	{
	  //------------------------------------ 3) mainchain-sidechain H-bonds
	  //------------------------------------ 6) sidechain-sidechain H-bonds
	  //------------------------------------ 7) sidechain-HET_NWatoms H-bonds
	  
	  if((p[so]->r1.field1[0]!='H' && p[sf]->r1.field1[0]!='H') 
	     && (((strcmp(p[so]->r1.aname,"  N  ") == 0 
		   || strcmp(p[so]->r1.aname,"  O  ")==0) 
		  && (strcmp(p[sf]->r1.aname,"  N  ") !=0 
		      && strcmp(p[sf]->r1.aname,"  O  ")!=0)
		  )
		 || ((strcmp(p[sf]->r1.aname,"  N  ")==0
		      || strcmp(p[sf]->r1.aname,"  O  ")==0)
		     && (strcmp(p[so]->r1.aname,"  N  ")!=0 
			 && strcmp(p[so]->r1.aname,"  O  ")!=0)
		     )
		 )
	     )
	    {
	      khb--;
	      hb_pick[i] = 0;
	    }
	  else
	    if((p[so]->r1.field1[0]!='H' && p[sf]->r1.field1[0]!='H')
	       && ((strcmp(p[so]->r1.aname,"  N  ")!=0 
		    && strcmp(p[so]->r1.aname,"  O  ")!=0) 
		   && (strcmp(p[sf]->r1.aname,"  N  ")!=0 
		       && strcmp(p[sf]->r1.aname,"  O  ")!=0))
	       )
	      {
		khb--;
		hb_pick[i] = 0;
	      }
	    else
	      if(((p[so]->r1.field1[0]!='H' && p[sf]->r1.field1[0]=='H') &&
		  (strcmp(p[so]->r1.aname,"  N  ")!=0 && strcmp(p[so]->r1.aname,"  O  ")!=0))
		 ||((p[sf]->r1.field1[0]!='H' && p[so]->r1.field1[0]=='H') &&
		    (strcmp(p[sf]->r1.aname,"  N  ")!=0 && strcmp(p[sf]->r1.aname,"  O  ")!=0))
		 )
		{
		  if( !(strcmp(p[sf]->r1.rname,"HOH")==0 ||
			strcmp(p[so]->r1.rname,"HOH")==0)   )
		    {
		      hb_pick[i]=0;
		      khb--;
		    }
		}

	} // --- End choice[3] ---

      //---------------- Remove *ALL* H-bonds with non-Water HETATOMs
      //---------------- Old choices: 4 and 9
      if(choice[4] == 1) 
	{
	  //---------------------------------- 4) mainchain-HET_NWatom H-bonds
	  //---------------------------------- 10) sidechain-HET_NWatom H-bonds 
	  //---------------------------------- 9) HET_atom-HET_NWatom H-bonds

	  //--- 2006:04 Sandeep		- The atom-type is what has to be compared
	  //---				  Note that PDB file has atom-type+positional
	  //---				  label as the atom name ex: NH1 and not just N
	  //---				- Hence, changed the conditions below

	  if( ((p[so]->r1.field1[0]!='H' && p[sf]->r1.field1[0]=='H') &&
	       (strcmp(p[so]->r1.aname,"  N  ")==0||strcmp(p[so]->r1.aname,"  O  ")==0))||
	      ((p[sf]->r1.field1[0]!='H' && p[so]->r1.field1[0]=='H') &&
	       (strcmp(p[sf]->r1.aname,"  N  ")==0||strcmp(p[sf]->r1.aname,"  O  ")==0))
	      )
	    {
	      if( !(strcmp(p[sf]->r1.rname,"HOH")==0 ||
		    strcmp(p[so]->r1.rname,"HOH")==0)   )
		{
		  hb_pick[i]=0;
		  khb--;
		}
	    }
	  else if(choice[3] == 0){
              if(((p[so]->r1.field1[0]!='H' && p[sf]->r1.field1[0]=='H') &&
                  (strcmp(p[so]->r1.aname,"  N  ")!=0 && strcmp(p[so]->r1.aname,"  O  ")!=0))
                 ||((p[sf]->r1.field1[0]!='H' && p[so]->r1.field1[0]=='H') &&
                    (strcmp(p[sf]->r1.aname,"  N  ")!=0 && strcmp(p[sf]->r1.aname,"  O  ")!=0))
                 )
                {
                  if( !(strcmp(p[sf]->r1.rname,"HOH")==0 ||
                        strcmp(p[so]->r1.rname,"HOH")==0)   )
                    {
                      hb_pick[i]=0;
                      khb--;
                    }
                }
	      } //End of if(choice[3] == 0)	
	  else
	    if(p[so]->r1.field1[0]=='H' && p[sf]->r1.field1[0]=='H')
	      {
		if( !(strcmp(p[sf]->r1.rname,"HOH")==0 &&
		      strcmp(p[so]->r1.rname,"HOH")==0)   )
		  {
		    hb_pick[i]=0;
		    khb--;
		  }
	      }
	} // --- End of choice[4]
     } // End of if(hb_pick[i] == 1)
    } // End of for(i = ... )


    //----------- END OF H-BOND FILTERING -------------
    // Display the number of H-bonds remaining
//    if(usage >= 1)
//    	cout<<"\nNumber of H-bonds remaining after removal of H-bonds with waters (Default option): "<<khb<<endl;
//    else
//    	cout<<"\nNumber of H-bonds remaining: "<<khb<<endl;
	

      // -------------------- Remove hydrophobic tethers ---------------      
//      if(choice[5] == 1) 
//	{
	  // Set a flag here or just use choice[5] as a check
	  // while building the hbonds file
	  // ********* SEE BELOW ******************
//          cout<<"\nNumber of H-phobic tethers remaining: 0"<<endl;
//	}
//      else
//          cout<<"\nNumber of H-phobic tethers remaining: "<<no_tether<<endl;


//------------------ Prompt user for energy filter ------------------
//------------------ Implement energy filter ------------------------
      khb_max=nhb;
      for(i=1;i<=nhb;i++)
        {
          if(hb_pick[i]>0)
            {
              khb_max=i;
              break;
            }
        }
      khb_min=1;
      for(i=nhb;i>=1;i--)
        {
          if(hb_pick[i]>0)
            {
              khb_min=i;
              break;
            } 
        }   
      energy_max=hb_energy[khb_max];

	if(usage >= 1) {
            hb_ecut=100;
            energy_max = max_energy;
          }
          else{
            
		system("clear");

	   // Display the number of H-bonds remaining
	    cout<<"\nNumber of H-bonds remaining: "<<khb<<endl;
	      if(choice[5] == 1)
        	{
          	  cout<<"\nNumber of H-phobic tethers remaining: 0"<<endl;
	        }
	      else
        	  cout<<"\nNumber of H-phobic tethers remaining: "<<no_tether<<endl;

            cout<<endl<<"\t\t Filter on Hydrogen Bond Energy:";
	    cout<<endl<<"\t\t ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

            cout<<endl<<"\tAll current hydrogen bonds have energies between:"<<endl;
            cout<<"\t"<<setiosflags(ios::showpoint|ios::fixed)<<setprecision(4)
                  //            <<setw(13)<<hb_energy[khb_min]<<" <= H-Bond energies <= "
                <<setw(11)<<hb_energy[khb_min]<<" Kcal/Mol  to  "
                <<setiosflags(ios::showpoint|ios::fixed)<<setprecision(4)
                << hb_energy[khb_max]<< " Kcal/Mol" << endl;

//            hb_ecut = 100;
            cout<<"\n\tEnter a maximum acceptable hydrogen bond energy (in Kcal/Mol)";
	    cout<<"\n\t(-1.0 is a reasonable cutoff in general): ";
            cin >> energy_input;
            energy_max = energy_input;
          }


          if(energy_max > hb_energy[khb_max]) 
		energy_max=hb_energy[khb_max];

//          khb2 = int( (khb*hb_ecut)/100);
//          khb=0;

          for(i=nhb;i>=1;i--)
            {
              if(hb_pick[i] == 1)
                {
                 if(hb_energy[i] > energy_max)
		   {
                     hb_pick[i]=0;
		     khb--;
		   }
                }
            }

// ------------------------------- END OF ENERGY-BASED FILTERING --------------------

	  if(usage >= 1)
	    cout<<"\nNumber of H-bonds remaining after filtering with default energy cutoff of "<<energy_max<<" Kcal/Mol: "<<khb<<endl;





//Output allbonds_SEfilt --- Generated by main.cpp
//Actually, output by FORTRAN code




//------------------ Output "hbonds_SEfilt" file --------------------
      khb2=0;
      for(i=1;i<=nhb;i++)
        {
          if(hb_pick[i]>0)
            {
              khb2++;
              hb_pick[khb2]=i;
            }
        }
      if(khb2!=khb)
        {
          cout<<"Error: Inconsistency in number of selected hydrogen bonds\n";
          cout<<"khb2="<<khb2<<" khb="<<khb<<endl;
          exit(33);
        }

  //**********************************************************************
  // The following chunk of code generates a single column list of integers
  // that correspond to specific hydrogen bonds and hydrophobic tethers in
  // a proflexdataset file. The FORTRAN routines read in this list, and use
  // it to select which hydrogen bonds to include in the analysis. All the
  // selection of hbonds has been moved from the FORTRAN to the C++.
  // Fixed 9.4.02 BMH.
  //**********************************************************************
  //    usage == 1 --> -non -h runoptions
  //    usage == 2 --> -non -p runoptions
  //    usage == 0 --> otherwise
  //**********************************************************************

  hblist.open(outfile[0],ios::out);

    for( i = 1; i <= khb; i++){
      j = hb_pick[i];
      hblist << hb_id[j] << endl;
    }

    if(choice[5] != 1) //--- INCLUDE TETHERS
      {
        for( i = 1; i <= no_tether; i++)    {
           hblist << nhb+i << endl;
        }
      }

  hblist.close();
  hblist.clear();



//--- Output the criteria file with summary of all the options chosen ---

//aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
//		place_hbonds(khb);
//-----------------------------------------------------Record selection criteria
        fptr.open("qXyZaB.criteria",ios::out);
	nsc=0;
	setiosflags(ios::showpoint|ios::fixed);
	setprecision(2);
	setw(6);
//	fptr<<"---------------------------------------------------------------- "<<nfile<<endl;
	fptr<<"        Hydrogen bond geometric selection criteria:"<<endl;
	fptr<<"        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
        fptr<<"        Donor-Acceptor distance <= "<<xda<<" Ang. without sulfur"<<endl;
        fptr<<"        Donor-Acceptor distance <= "<<s_xda<<" Ang. with sulfur"<<endl;
        fptr<<"        Donor-Acceptor distance <= "<<sb_xda<<" Ang. salt bridge"<<endl;
        fptr<<"        Hydrogen-Acceptor distance <= "<<xha<<" Ang. without sulfur"<<endl;
	fptr<<"        Hydrogen-Acceptor distance <= "<<s_xha<<" Ang. with sulfur"<<endl;
        fptr<<"        Hydrogen-Acceptor distance <= "<<sb_xha<<" Ang. salt bridge"<<endl;
 	fptr<<"        Donor-Hydrogen-Acceptor angle >= "<<adha/raddeg<<" Deg."<<endl;
        fptr<<"        Hydrogen-Acceptor-Preacceptor angle >= "<<ahap/raddeg<<" Deg."<<endl<<endl;


	fptr<<"        Hydrogen bond energetic selection criteria:"<<endl;
	fptr<<"        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
  	fptr<<"        Maximum allowed hydrogen bonding energy = "<< setprecision(4) << energy_max<<" Kcal/Mol"<<endl;
	fptr<<endl;

	fptr<<"        Non-covalent bond selection criteria:"<<endl;
	fptr<<"        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
	nsc+=17;

	if(!usage){
	  if(choice[2] == 0) 
	    {
		fptr<<"        Default: Remove *ALL* H-bonds involving waters"<<endl;
		nsc++;
	    }
	  else
	    {
                fptr<<"        "<<screen_msgs[2]<<endl;
  	  	nsc+=3;
	    }

	  for(i=3;i<=maxnopt;i++)
	    {
	      if(choice[i]==1)
		{
		  fptr<<"        "<<screen_msgs[i]<<endl;
		  nsc++;
		}
	    }
	}
	else {
	  fptr<<"        Keep all hydrogen bonds and hydrophobic tethers except"<<endl;
          fptr<<"        those including waters ----> Default assumption"<<endl;
	  fptr<<"        Filter on Donor-Hydrogen-Acceptor angle"<<endl;
	  nsc+=3;
	}

	fptr.close();

} // --- End pick_hbonds()












/************************************************************/
/* There are several options available to provide the user  */
/* with control over how hydrogen bonds are identified in a */
/* protein. This routine will list the options to the       */
/* screen, and store the selected options for later use.    */
/************************************************************/
void list::pick_hbonds_options(void) {

  int i,ch,isumrule;
  int flag=0,choice_count;
  char ans;

  /************************************************************/
  /* Initialize choices                                       */
  /************************************************************/
  for( i = 0; i <= maxnopt; i++) {
    choice[i] = 0;
  }

  isumrule=11;

  /************************************************************/
  /* Display screening option menu                            */
  /************************************************************/
  // screen[1] = "\tDefaults: (not selecting 1,2,3, or 4)\n\t\t Will keep all hydrophobic tethers \n\t\t and H-bonds except those including water";
  screen_msgs.push_back("Finished selecting options; Continue to next step");
  screen_msgs.push_back("Use defaults (keep hydrophobic tethers and H-bonds except those including water)");
  screen_msgs.push_back("Keep *ALL* H-bonds involving WATER \n\t\t(Select this ONLY if you have included only\n\t\t buried water molecules in the input PDB file)");
  screen_msgs.push_back("Remove *ALL* H-bonds involving sidechain atoms");
  screen_msgs.push_back("Remove *ALL* H-bonds involving non-water HETATOMs");
  screen_msgs.push_back("Remove hydrophobic tethers");

  /************************************************************/
  /* This while loop will keep printing the options until the */
  /* user enters "f" for finished selecting options.          */
  /************************************************************/
  flag = 0;
  choice_count = 0;
  while( flag == 0 ) {

    system("clear");

    /************************************************************/
    /* Print the hydrogen bond screening options.               */
    /************************************************************/

    //--- 2006:03 SN - change the menu name
    //cout<<endl<<"\t\t\tH-BOND SCREENING OPTIONS MENU"<<endl<<endl;
    cout<<endl<<"\t\t\tNON-COVALENT BOND SCREENING OPTIONS MENU"<<endl<<endl;
    // SN 2006:02	--- Default option has to be displayed irrespective of what
    //			    options will be selected by the user
    for( i = 2; i <= maxnopt; i++ ) { 
      if(choice[i]==0) {
        cout<<"\t"<<setw(2)<<i-1<<" : "<<screen_msgs[i]<<endl;
      }
    }
    cout<<"\t"<<setw(2)<<maxnopt<<" : "<<screen_msgs[1]<<endl;
//    cout<<endl<<screen_msgs[1]<<endl;

    cout<<endl<<"\t f : "<<screen_msgs[0];
    cout<<endl<<"\t s : Stop program";
    cout << endl << endl;

    do{
	    cout<<"\tEnter Option: ";
	    cin >> ans;
    
	    //--- 2006:03 SN	Ensure that user selects at least one option
	    if(choice_count == 0 && (ans == 'f' || ans == 'F'))
		{
	  cout<<"\nPlease select at least one valid (1-5) option before finish!\n\n";
		}
	    else
		break;
    }while(1);

    //--- 2006:03 SN	- If user selects '5' (defaults) prog should behave like 'f'
    //---		- If other options are selected along with '5', make others
    //---		  zero and initiate 'f'
    if(ans == '5'){
      flag = 1;
      choice[1] = 1;
      choice_count++;
    }
    else if( ans == 'f' || ans == 'F'){
      flag = 1;
    }  
    else if( ans == 's' || ans == 'S' ) {
      exit(-1);
    }
    else {
      ch = ans - 48; // 48 is the ASCII code for '0'

      /************************************************************/
      /* Check for invalid options.                               */
      /************************************************************/
      if( ch < 1 || ch > maxnopt) {
        cout<<endl<<"\tInvalid option!"<<endl;
      }
      else
        if(choice[ch]==1 || choice[ch]==-1)    // OR !=0 ???? aaaaaaaaaaaaaaa
          cout<<endl<<"\tYou already chose this option!"<<endl;
        else {
         //------------------------------------------------Record option
          choice[ch+1]=1;
	  choice_count++;
        }

    } // else ASCII

  } // End while(!flag)
    
} // End pick_hbond_options()

