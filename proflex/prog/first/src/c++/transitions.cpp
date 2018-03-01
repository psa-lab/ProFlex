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

/* THIS SUBROUTINE IS THE C++ VERSION OF getfiles.f */


/**********************************************************************/
/* Create the outputfiles used by FIRST.                              */
/* Modus Operandi: parsing the global variable outputfile (defined in */
/* class.h), this function looks for the file *_analysis.log. From    */
/* this file, it looks for line containing "SELECTION CRITERIA" phrase*/
/* At the end of this line is the 'run number' which can be maximum 9999*/
/**********************************************************************/
void getfiles(void) {

  int len, i, i1, i2, i3, i4, itemp, n_root;

  char analfile[100],line[90],c1,c2,c3,c4;

  ifstream fptr;

  len = strlen(outputfile);
  n_root = len-15;
  // open _analysis.log file
  strncpy( analfile, outputfile, n_root );
  analfile[n_root]='\0';

  strcat(analfile,"_analysis.log");
  //cout<<endl<<"analfile = "<<analfile<<endl<<endl;
  fptr.open(analfile,ios::in);

  // Determine nfile
  nfile = 1;
  //cout<<endl<<"nfile="<<nfile<<endl;


  if( !fptr.fail() ) {

    fptr.getline(line,90);
    while( !fptr.eof() ) {
      if( !strncmp( line+23,"SELECTION CRITERIA",18) ) {
	nfile++;
      }
      fptr.getline( line, 90 );
      fptr.getline( line, 90 );
    }

  }

  // Construct the set of input/output names
  i4 = nfile/1000;
  itemp = nfile%1000;
  i3 = itemp/100;
  itemp = itemp%100;
  i2 = itemp/10;
  i1 = itemp%10;

  c4='0'+i4;
  c3='0'+i3;
  c2='0'+i2;
  c1='0'+i1;

  for( i = 0; i < 8; i++ ) {
    strncpy(outfile[i],outputfile,n_root);
    outfile[i][n_root]='\0';
  }

// -- SN 02:2006 - Changing the names of some of the output files
  strcat(outfile[0],"_h-bonds_SEfilt.");
  strcat(outfile[1],"_bond_wt.");
  strcat(outfile[2],"_rdecomp.");
  strcat(outfile[3],"_sdecomp.");
  strcat(outfile[4],"_fdecomp.");
  strcat(outfile[5],"_graphic.");
  strcat(outfile[6],"_Rscript.");
  strcat(outfile[7],"_analysis.log");

  outfile[0][n_root+16] =c4;
  outfile[0][n_root+17] =c3;
  outfile[0][n_root+18] =c2;
  outfile[0][n_root+19] =c1;
  outfile[0][n_root+20] ='\0';

  for( i = 1; i < 7; i++ ) {
    outfile[i][n_root+9]  = c4;
    outfile[i][n_root+10] = c3;
    outfile[i][n_root+11] = c2;
    outfile[i][n_root+12] = c1;
    outfile[i][n_root+13] = '\0';
  }

}
/**********************************************************************/


/**********************************************************************/
/* FIRST uses it's own version of the pdb file, which includes all of */
/* the structural information, as well as the information on which    */
/* atoms are involved in central-force (CF) bond, torsional-force (TF)*/
/* bonds, and hydrogen bonds (HB). These data are stored in the file  */
/* with an extension _proflexdataset.                                   */
/**********************************************************************/
void list::read_proflexdataset(void) {

  int **temp_link, i, j, k, mlt, iflag;
  int flag1=0, so, s, sf, flag2=0;
  int **temp_hbond, *temp_hb_type, *temp_hb_id, l, id_hb;

  float *temp_hb_energy;

  char line[90],string[20];

  ifstream fptr;

  p = (node **) new node*[100000];

  // Open the *_proflexdataset file
  fptr.open( outputfile, ios::in );
  if( fptr.fail() ) {
    cout<<endl<<"\tproflexdataset file not found!"<<endl;
    exit(8);
  }

  int unknown=0;
  node *fresh, *temp;

  //------since read_proflexdataset has been invoked, we know user is using previous dataset - Sameer 12.Mar.04
  real_hb_no = 0;

  while( !fptr.eof() ) {

    //------------------------------------Read ATOM or HETATM or TER
    fptr.getline(line,90);
    if( !strncmp( line, "ATOM",   4 ) ||
	!strncmp( line, "HETATM", 6 ) ||
	!strncmp( line, "TER",    3 ) ) {

      fresh = new node;
      if( !strncmp(line,"TER",3) ) {
	l=strlen(line);
	for( i = 0; i < (85-l); i++ ) {
	  line[l+i]=' ';
	}
	line[84]='\0';
      }

      strncpy(fresh->r1.field1,line,6);      // PDB field identifier (ie. ATOM, HETATM, REMARK)
      strncpy(fresh->r1.sr_no,line+6,5);     // Atom number
      strncpy(fresh->r1.aname,line+11,5);    // Atom name
      strncpy(fresh->r1.altloc,line+16,1);   //
      strncpy(fresh->r1.rname,line+17,3);    // Amino-acid name
      strncpy(fresh->r1.chain,line+20,2);    // Chain ID
      strncpy(fresh->r1.res_no,line+22,4);   // Residue number
      strncpy(fresh->r1.code,line+26,4);
      strncpy(fresh->r1.strx,line+30,8);     // X-coordinate
      strncpy(fresh->r1.stry,line+38,8);     // Y-coordinate
      strncpy(fresh->r1.strz,line+46,8);     // Z-coordinate
      strncpy(fresh->r1.occupancy,line+54,6);// Temperature factor (B_value)
      strncpy(fresh->r1.temp_string,line+60,6);     //
      strncpy(string, line+66, 5 );          // some kind of internal FIRST identifier.
      fresh->r1.mod_res_no = atoi(string);   //
      fresh->r1.DAH_type = line[72];

      fresh->r1.coord[1]    = atof(fresh->r1.strx);
      fresh->r1.coord[2]    = atof(fresh->r1.stry);
      fresh->r1.coord[3]    = atof(fresh->r1.strz);
      fresh->r1.occ         = atof(fresh->r1.occupancy);
      fresh->r1.temperature = atof(fresh->r1.temp_string);
      fresh->r1.ri_sr_no    = atoi(fresh->r1.sr_no);

      if( !strncmp(line,"ATOM",4) || !strncmp(line,"HETATM",6) ) {
	no_atoms += 1;
	p[no_atoms] = fresh;
      }
      else
	fresh->r1.ri_sr_no = 0;

      if( start == NULL ) {
	start = fresh;
	last  = fresh;
      }
      else {
	fresh->next  = start;
	start->prior = fresh;
	start = fresh;
      }
    }

    // Read central-force bonds(CF) from proflexdataset file
    else if( !strncmp(line,"REMARK:CF",9) ) {
      if( flag1 == 0 ) {

	mult = new int[no_atoms+1];
	nrot_bond = (int**) new int*[no_atoms+1];
	temp_link = (int**) new int*[no_atoms+1];

	for( i = 0; i <= no_atoms; i++ ) {
	  nrot_bond[i] = (int*) new int[3];
	  temp_link[i] = (int*) new int[maxr];
	  mult[i]=0;
	}
	flag1 = 1;
      }

      strncpy( string, line+10, 5 );
      so = atoi(string);
      strncpy(string,line+16,5);
      sf = atoi(string);
      if( so == sf ) {
	cout<<endl<<"\tCovalent bond self connection ignored!"<<endl;
	cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl;
	write_output(1,sf,0,p[sf]->r1.aname,p[sf]->r1.rname,p[sf]->r1.res_no,
		     p[sf]->r1.code,p[sf]->r1.chain,p[sf]->r1.aname,p[sf]->r1.aname,
		     p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname);
      }
      else {
	mult[so]++;
	temp_link[so][mult[so]] = sf;
	mult[sf]++;
	temp_link[sf][mult[sf]] = so;
      }
    }

    // Read Torsional Force constraints (TF) from proflexdataset file
    else if( !strncmp(line,"REMARK:TF",9) ) {
      strncpy(string,line+10,5);
      so = atoi(string);
      strncpy(string,line+16,5);
      sf = atoi(string);

      nrot_bond_count++;
      nrot_bond[nrot_bond_count][1]=so;
      nrot_bond[nrot_bond_count][2]=sf;
    }

    // Read Hydrogen Bond constraints from proflexdataset file
    else if( !strncmp( line, "REMARK:HB", 9 ) ) {

      nhb++; // by incrementing here, all subsequent arrays are indexed from 1 !

      /* of these HB records, count only the hydrogen bonds and salt bridges i.e. every HB record
       * except the hydrophobic tethers - Sameer 12.Mar.04
       */
      if( strncmp( line+55, "PH hydr", 7 ) )
			real_hb_no++;



      if( flag2 == 0 ) {
	temp_hb_energy = new float[no_atoms+1];
	temp_hb_type   = new int[no_atoms+1];
	temp_hb_id     = new int[no_atoms+1];
	temp_hbond     = (int**) new int*[no_atoms+1];
	for( i = 0; i <= no_atoms; i++ )
	  temp_hbond[i] = (int*) new int[4];

	flag2 = 1;
      }

      // id's are listed 1 - number of H-bonds in *proflexdataset

      // the explicit atoi's and atof's are OK for now, the dataset file
      // has distinct whitespace delimiters between data. BMH 3.27.00

      id_hb = atoi( line+9 );

      temp_hb_energy[nhb] = atof( line+14 );

      temp_hbond[nhb][1]  = atoi( line+30 );

      temp_hbond[nhb][2]  = atoi( line+38 );

      temp_hbond[nhb][3]  = atoi( line+46 );

      if( !strncmp(line+55,"HB Dsp2 Asp2",12) )
	temp_hb_type[nhb] = 0;
      else if( !strncmp(line+55,"HB Dsp2 Asp3",12) )
	temp_hb_type[nhb] = 1;
      else if( !strncmp(line+55,"HB Dsp3 Asp2",12) )
	temp_hb_type[nhb] = 2;
      else if( !strncmp(line+55,"HB Dsp3 Asp3",12) )
	temp_hb_type[nhb] = 3;
      else if( !strncmp(line+55,"SB no energy",12) )
	temp_hb_type[nhb] = 4;
      else if( !strncmp(line+55,"PH hydr phob",12) ) {
	temp_hb_type[nhb] = 5;
	//BMH removed next line. Count tethers in pick_hbonds line 959.
	// no_tether++;
      }
      else {
	temp_hb_type[nhb] = 6;	// unknown HB type!
	unknown++;
      }

      // Error checks
      if( id_hb < 1 ) {
	cout<<endl<<"\a\tERROR: Current hydrogen bond ID = "<<id_hb<<endl;
	cout<< "\t       FIRST must terminate until ID label is changed."<<endl;
	cout<< "\t       (Set ID > 0)"
	    <<endl<<endl;
	exit(25);
      }
      if( id_hb >= no_atoms ) {
	cout<<endl<<endl
	    <<"\a\tERROR: Current hydrogen bond ID = "<<id_hb<<endl;
	cout<< "\t       FIRST must terminate until ID label is changed."<<endl;
	cout<< "\t       (Set ID < " << no_atoms << ")"
	    <<endl<<endl;
	exit(26);
      }
      temp_hb_id[nhb] = id_hb;
    }
  }

  // Close _proflexdataset file
  fptr.close();

  // Check if the file was empty
  if( no_atoms == 0 ) {
    cout<<"\a\tThere were no atom records in the file!"<<endl<<endl;
    exit(11);  }

  //--------------------------------------------Check if any unknown H-bonds detected
  if( unknown > 0 ) {
    cout << endl
	 << "\a\tNumber of unknown hydrogen bond types detected = "
	 << unknown << endl;
    cout << "\tGeometrical screening will proceed as if a salt bridge."
	 << endl;
    cout << "\tEnter  c  to continue: ";
    char ans;
    cin >> answer;

    if( strlen(answer) == 1 )
      ans = answer[0];
    else ans = '@';

    if( ans == 'c' )
      ans = 'C';
    if( ans != 'C' )
      exit(5);
  }

  // Marking the begining of chains
  temp = last;
  chainptr[0] = temp;

  while( temp != start ) {
    temp = temp->prior;
    if( strcmp(temp->r1.chain, chainptr[chaincount]->r1.chain) ) {
      chainptr[chaincount+1] = temp;
      ++chaincount;
    }
  }

  // Allocate memory for permanent arrays
  link_noHB = (int **) new int*[no_atoms+1];

  for( i = 0; i <= no_atoms; i++ )
    link_noHB[i] = (int *) new int[mult[i]+1];

  mult_net  = new int[no_atoms+1];
  hb_energy = new float[nhb+1];
  hb_type   = new int[nhb+1];
  hb_id     = new int[nhb+1];
  hbond     = (int **) new int *[nhb+1];
  for( i = 0; i <= nhb; i++ )
    hbond[i] = (int *) new int[4];

  // Copy temp arrays to permanent arrays
  for( i = 1; i <= no_atoms; i++ ) {
    for( j = 1; j <= mult[i]; j++ )
      link_noHB[i][j] = temp_link[i][j];
  }

  for( i = 1; i <= no_atoms; i++ ) {
    mlt = 0;
    for( j = 1; j <= mult[i]; j++ ) {

      // Covalent degeneracy error check
      iflag = -1;
      for( k = j+1; k <= mult[i]; k++ )	{
	if(link_noHB[i][j] == link_noHB[i][k])
	  iflag = 1;
      }
      if( iflag < 0 ) {
	mlt++;
	link_noHB[i][mlt] = link_noHB[i][j];
      }
      else {
	so = i;
	sf = link_noHB[i][j];
	if( so < sf ) {
	  cout<<endl<<"\tDuplicate covalent bond connection ignored!"<<endl;
	  cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl;
	  write_output(1,so,sf,p[so]->r1.aname,p[so]->r1.rname,p[so]->r1.res_no,
		       p[so]->r1.code,p[so]->r1.chain,p[sf]->r1.aname,
		       p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
		       p[sf]->r1.chain);
	}
      }
    }
    mult[i] = mlt;
  }

  for( i = 1; i <= nhb; i++ ) {
    hb_energy[i] = temp_hb_energy[i];
    hb_type[i]   = temp_hb_type[i];
    hb_id[i]     = temp_hb_id[i];
    hbond[i][1]  = temp_hbond[i][1];
    hbond[i][2]  = temp_hbond[i][2];
    hbond[i][3]  = temp_hbond[i][3];
  }

  // Check for HB ID degeneracy
  for( i = 1; i <= no_atoms; i++ ) {
    temp_hb_id[i] = 0;
  }

  // ERROR check.
  for( i = 1; i <= nhb; i++ ) {
    temp_hb_id[ hb_id[i] ]++;
    if( temp_hb_id[ hb_id[i] ] > 1 ) {
      cout << endl << endl
	   << "\a\tERROR: Degeneracy found in hydgrogen bond ID = "
	   << hb_id[i] << " " << i <<  " " << temp_hb_id[ hb_id[i] ]
	   << endl << "\t       FIRST must terminate until ID labeling is corrected!"
	   << endl << endl;
      exit(27);
    }
  }

  // Expand multiplicity
  for( i = 0; i <= no_atoms; i++ ) {
    mult_net[i]=mult[i]; //--Copy original multiplicities into new array
  }
  for( i = 1; i <= nhb; i++ ) {
    s  = hbond[i][2];	 //--include multiplicites due to hbonds
    sf = hbond[i][3];
    mult_net[s]++;
    mult_net[sf]++;
  }

  // Allocate memory for link_net[][]
  link_net =(int **) new int*[no_atoms+1];
  for( i = 0; i <= no_atoms; i++ )
    *(link_net+i)=(int *) new int[mult_net[i]+1];

  // Expand connectivity table
  for( i = 1; i <= no_atoms; i++ ) {
    for( j = 1; j <= mult[i]; j++ ) {
      link_net[i][j]=link_noHB[i][j];
    }
    mult_net[i] = mult[i];
    delete [] *(temp_link+i);
  }

  for( i = 1; i <= nhb; i++ ) {
    s  = hbond[i][2];	//--include multiplicites due to hbonds
    sf = hbond[i][3];
    mult_net[s]++;
    mult_net[sf]++;
    link_net[s][mult_net[s]]   = sf;
    link_net[sf][mult_net[sf]] = s;
  }

  delete [] temp_link;
  delete [] temp_hb_type;
  delete [] temp_hb_id;
  delete [] temp_hb_energy;
}
/**********************************************************************/


/**********************************************************************/
/* what does this do?                                                 */
/**********************************************************************/
void list::check(void)
{
  int i,j,k,m,iflag,mlt,so,sf;
  int tmp1,tmp2,tmp3,tmp4,tmp5;
  int link_counter = 0; // To count links before flagging valency warning - Sameer
  float tmp;
  //--------------------------------------------No groups in PDB file
  if(no_atoms==0)
    {
      cout<<endl<<"STOP: Please check original PDB file!  No atoms found!"<<endl;
      exit(11);
    }
  for(i=1;i<=no_atoms;i++)
    {
      for(j=1;j<=mult_net[i];j++)
	{
	  //-------------------------------------------Self connect error check
	  if(link_net[i][j]==i)
	    {
	      so = i;
	      sf = link_net[i][j];
	      cout<<endl<<"Self connection ignored!"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl;
	    write_output(1,sf,0,p[sf]->r1.aname,p[sf]->r1.rname,p[sf]->r1.res_no,
		  p[sf]->r1.code,p[sf]->r1.chain,p[sf]->r1.aname,p[sf]->r1.aname,
		  p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname);

	    //----------------------------------Remove that bond
	    for(k=j;k<=mult_net[so]-1;k++)
	      {
		link_net[so][k]=link_net[so][k+1];
	      }
	    link_net[so][mult_net[so]]=-1;
	    mult_net[so]--;
	    }
	}
    }
  for(i=1;i<=no_atoms;i++)
    {
      mlt = 0;
      for(j=1;j<=mult_net[i];j++)
	{
	  //-------------------------------------------Degeneracy error check
	  iflag = -1;
	  for(k=j+1;k<=mult_net[i];k++)
	    {
	      if( link_net[i][j] == link_net[i][k] ) iflag = 1;
	    }
	  if( iflag < 0 )
	    {
	      mlt++;
	      link_net[i][mlt]=link_net[i][j];
	    }
	  else
	    {
	      so=i;
	      sf=link_net[i][j];
	      if( so < sf )
		{
		 cout<<endl<<"\tHydrogen-bond overlap connection ignored!"<<endl;
		 cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl;
	   write_output(1,so,sf,p[so]->r1.aname,p[so]->r1.rname,p[so]->r1.res_no,
			p[so]->r1.code,p[so]->r1.chain,p[sf]->r1.aname,
			p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
			p[sf]->r1.chain);
//--------------------------------------------------- Remove H-bond from the list
	   m = 0;
	   iflag = -1;
	   for(k=1;k<=nhb;k++)
	     {
	       if( (so == hbond[k][2] && sf == hbond[k][3]) ||
		   (sf == hbond[k][2] && so == hbond[k][3]) )
		 {
		   if( iflag < 0 )
		     {
		       iflag = 1;
		       m++;
		       hbond[m][1] = hbond[k][1];
		       hbond[m][2] = hbond[k][2];
		       hbond[m][3] = hbond[k][3];
		       hb_energy[m] = hb_energy[k];
		       hb_type[m] = hb_type[k];
		       hb_id[m] = hb_id[k];
		     }
		 }
	       else
		 {
		   m++;
		   hbond[m][1] = hbond[k][1];
		   hbond[m][2] = hbond[k][2];
		   hbond[m][3] = hbond[k][3];
		   hb_energy[m] = hb_energy[k];
		   hb_type[m] = hb_type[k];
		   hb_id[m] = hb_id[k];
		 }
	     }
	   nhb = m;
		}
	    }
	}
      mult_net[i] = mlt;
    }
  int flag7 = 0;
  //------------------------------------------Valence check Begins
  for(i=1;i<=no_atoms;i++)
    {
	/*--------------------------- Real Valency check  ------------------*/
    /* The valency is taken as the no of bonds of an atom. This is fine, as long as we don't add
   	 * up the bonds to pseudo atoms as well. This becoms a problem when re-using proflexdataset in a run
  	 * as then bonds to psuedo atoms too leads FIRST to decide valency as incorrectly high, and
   	 * give these warnings - Sameer
   	 */

		link_counter = 0;
      	for(j=1;j<=mult[i];j++)
		{
	  		sf=link_noHB[i][j];
			if( strcmp(p[sf]->r1.rname,"XXX") != 0 )
							link_counter++;
			//cout<< "["<<p[sf]->r1.rname<<" : " << link_counter<<"]\t";
		}
	/*--------------------------- Real Valency check ends - Sameer  ------------------*/
      if(p[i]->r1.aname[2]=='H' || p[i]->r1.aname[2]=='D')
	{
		if(link_counter > 1)
	  //if(mult[i]>1)
	    {
	      so=i;
	      if(flag7==0)
		{
		  //system("clear");
		  // cout<<endl<<endl<<"\t\t    Checking valencies"<<endl;
		  cout<<endl<<endl<<"\t\t     Valency warnings"<<endl;
		  cout<<"---------------------------------------------------------------"<<endl;
		  flag7=1;
		}
	      cout<<endl<<"\tValency of atom "<<i<<" is "<<mult[i]<<", greater than its maximum valency (1)"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl;
	      write_output(1,so,0,p[so]->r1.aname,p[so]->r1.rname,p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname);
	      for(j=1;j<=mult[i];j++)
		{
		  sf=link_noHB[i][j];
		  write_output(1,sf,0,p[sf]->r1.aname,p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,p[sf]->r1.chain,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname);
		}
                                //cout << endl;
	    }
	}
      else if(p[i]->r1.aname[2]=='O')
	{
		if(link_counter > 2)
	  //if(mult[i]>2)
	    {
	      so=i;
	      if(flag7==0)
		{
		  //system("clear");
		  //cout<<endl<<endl<<"\t\t    Checking valencies"<<endl;
		  cout<<endl<<endl<<"\t\t     Valency warnings"<<endl;
		  cout<<"---------------------------------------------------------------"<<endl;
		  flag7=1;
		}
	      cout<<endl<<"\tValency of atom "<<i<<" is "<<mult[i]<<", greater than its maximum valency (2)"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl;
	      write_output(1,so,0,p[so]->r1.aname,p[so]->r1.rname,p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname);
	      for(j=1;j<=mult[i];j++)
		{
		  sf=link_noHB[i][j];
		  write_output(1,sf,0,p[sf]->r1.aname,p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,p[sf]->r1.chain,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname);
		}
                                //cout << endl;
	    }
	}
      else if( p[i]->r1.aname[2]=='N' &&
	       ( p[i]->r1.mod_res_no != 1 || ( p[i]->r1.mod_res_no == 1 &&
		 strncmp(p[i]->r1.aname+1," N  ",4)!=0 ) )  &&
	       !(strcmp(p[i]->r1.rname,"LYS")==0 &&
		 strncmp(p[i]->r1.aname+1," NZ ",4)==0)   )
	{
		if(link_counter > 3)
	  //if(mult[i]>3)
	    {
	      so=i;
	      if(flag7==0)
		{
		  //system("clear");
		  //cout<<endl<<endl<<"\t\t    Checking valencies"<<endl;
		  cout<<endl<<endl<<"\t\t     Valency warnings"<<endl;
		  cout<<"---------------------------------------------------------------"<<endl;
		  flag7=1;
		}
	      cout<<endl<<"\tValency of atom "<<i<<" is "<<mult[i]<<", greater than its maximum valency (3)"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl;
	      write_output(1,so,0,p[so]->r1.aname,p[so]->r1.rname,p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname);
	      for(j=1;j<=mult[i];j++)
		{
		  sf=link_noHB[i][j];
		  write_output(1,sf,0,p[sf]->r1.aname,p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,p[sf]->r1.chain,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname);
		}
                                //cout << endl;
	    }
	}
      else
	{
		if(link_counter > 4)
	  	//if(mult[i]>4)
	    {

	      	so=i;
	      	if(flag7==0)
			{
		  		//system("clear");
			  	//cout<<endl<<endl<<"\t\t    Checking valencies"<<endl;
			  	cout<<endl<<endl<<"\t\t     Valency warnings"<<endl;
			  cout<<"---------------------------------------------------------------"<<endl;
			  flag7=1;
			}
	      cout<<endl<<"\tValency of atom "<<i<<" is "<<mult[i]<<", greater than its maximum valency (4)"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl;
	      write_output(1,so,0,p[so]->r1.aname,p[so]->r1.rname,p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname,p[so]->r1.aname);
	      for(j=1;j<=mult[i];j++)
		{
		  sf=link_noHB[i][j];
		  write_output(1,sf,0,p[sf]->r1.aname,p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,p[sf]->r1.chain,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname,p[sf]->r1.aname);
		}
                                //cout << endl;
	    }
	}
    } // valence checks end
  if(flag7==1)
    {
      cout<<"---------------------------------------------------------------"<<endl;
      /*
	cout<<"Type  s  (Default --> continue) to stop: ";
	cin.seekg(0,ios::end);
	cin.clear();
	cin.get(answer,4);
	cin.seekg(0,ios::end);
	cin.clear();
	if( strlen(answer) == 1 ) ans = answer[0];
	else ans = '@@';
	if( ans == 's' ) ans = 'S';
	if(ans=='S') exit(-1);
      */
      cout<<endl<<endl;
    }
  //--------------------------------------------Fix up hydrogen bond type
  for(i=1;i<=nhb;i++)
    {
      so = hbond[i][1];
      sf = hbond[i][3];
      if(p[so]->r1.aname[2]=='S'||p[sf]->r1.aname[2]=='S')
	{
	  hb_type[i]+=6;
	}
    }
//------------------------------------Sort H-bonds from highest to lowest energy
//                                                    Use a Shell sort algorithm
  hb_energy[0]=1e25;
  nhb++;
  for(int gap = nhb/2; gap > 0; gap = gap == 2 ? 1 : (int)(gap/2.2))
    {
      for(i=gap;i<nhb;i++)
	{
	  tmp=hb_energy[i];
	  tmp1=hbond[i][1];
	  tmp2=hbond[i][2];
	  tmp3=hbond[i][3];
	  tmp4=hb_type[i];
	  tmp5=hb_id[i];
	  j=i;

	  for(;j>=gap && tmp > hb_energy[j-gap]; j-=gap)
	    {
	      hb_energy[j]=hb_energy[j-gap];
	      hbond[j][1]=hbond[j-gap][1];
	      hbond[j][2]=hbond[j-gap][2];
	      hbond[j][3]=hbond[j-gap][3];
	      hb_type[j]=hb_type[j-gap];
	      hb_id[j]=hb_id[j-gap];
	    }
	  hb_energy[j]=tmp;
	  hbond[j][1]=tmp1;
	  hbond[j][2]=tmp2;
	  hbond[j][3]=tmp3;
	  hb_type[j]=tmp4;
	  hb_id[j]=tmp5;
	}
    }
  nhb--;
//--------------------------------------------------------------Warning message
  if(nhb < (no_atoms)/65)
    {
      cout<<endl<<endl
	  <<"WARNING:  Low number of hydrogen bonds found!  H-bond # = "
	  <<nhb<<endl;
      if(nhb == 0) { exit(55);}
    }
}
/**********************************************************************/


/**********************************************************************/
/* Include the hydrogen bonds in the list of network connections.      */
/**********************************************************************/
void list::place_hbonds(int khb) {

  int i=0, j = 0, s = 0, sf = 0;

  // Copy link_noHB to link_net
  for( i = 1; i <= no_atoms; i++ ) {
    for( j = 1; j <= mult[i]; j++)
      link_net[i][j] = link_noHB[i][j];
    mult_net[i]=mult[i];	//------------------------Also copy multiplicity
  }

  // Update multiplicity and link table with H-bonds present
  for( i = 1; i <= khb; i++ ) {
    s  = hbond[hb_pick[i]][2];
    sf = hbond[hb_pick[i]][3];
    mult_net[s]++;
    mult_net[sf]++;
    link_net[s][mult_net[s]]   = sf;
    link_net[sf][mult_net[sf]] = s;
  }
}
/**********************************************************************/


