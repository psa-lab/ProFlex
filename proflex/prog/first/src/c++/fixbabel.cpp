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

/**********************************************************************/
/* THIS IS A C++ CODE FOR FIX_BABEL                                   */
/*                                                                    */
/* The following routines create a connectivity map (a linked list of */
/*                                                                    */
/*    Dec 31 2003, Sameer : A temporary fix to prevent debug msgs     */
/*    on stdout using DEBUG define                                    */
/* arrays) of a protein, based on the chemistry of the covalent, salt */
/* and hydrogen bonds present.                                        */ 
/*    Last revision 6.27.00 by AJR                                    */ 
/*    PorousCity aka FIRSTweb version, revised AJR 4.10.01            */
/*    modified 03.20.02 by AJR --added search for hydrophobics       */
/**********************************************************************/
/**********************************************************************/
void residue::get_residue(char *temp_res) {

  int i,j;

  strcpy( rname, temp_res );
  ifil >> natom_sum;
  occu[0] =- 1;

  for( i = 1; i <= natom_sum; i++ ) {
    occu[i] =- 1;
    ifil >> A_type[i];
    if( strlen(A_type[i]) == 1 ) {
      A_type[i][1] =  ' ';
      A_type[i][2] = ' ';
      A_type[i][3] = '\0';
    }
    else {
      if( strlen(A_type[i]) == 2 ) {
	A_type[i][2] = ' ';
	A_type[i][3] = '\0';
      }
    }		
    ifil >> A_mult[i];
    for( j = 1; j <= A_mult[i]; j++ ) {
      ifil >> A_link[i][j];
    }
  }

}
/**********************************************************************/


/**********************************************************************/
void list::make_connectivity(int usage) {

  int natom = 0, i, j, k, flag = 0, xh_count = 0, atom_found = 0; 
  int ma_count = 0, rnumber, nsulfur = 0, ns_bond, nhydro = 0; 
  int flag8 = 0, **temp_link,loop_flag=1,flag17=0, oxtflag=0;
  int so,sf,prev_occu=-1,bond_oxt[500][2],n_oxt=0,oflag[chaincount];
  int flag7=0,iix,iiy,iiz,ix,iy,iz,ua_count=0,count=0;
  int label[1000], m,a,jx,jy,jz,flag9=0,flag10,kx,ky,kz,flag12=0;
  int **u_nrot_bond,u_nrot_bond_count=0,ret_val=0,q,**chk_bond;
  int chk_bond_count=0,**poor_bond,poor_bond_count=0;
  int *skip_nhydro,skip_nhydro_count=0,*skip_hydro,skip_hydro_count=0;
  int *missing,all_iso=0,O_iso=0,n;
  float temp,dist[1000],**sxyz,no_bonds=0,*chk_bond_dist,*poor_bond_dist,factor;
  char name[4],type[3],**atype,ans, line[80], **missing_atoms;
  std::string dash1 ="------------------------------------------------------------";
  std::string dash2 ="-----------";
  std::string dash3 ="--------------------";

  node *tempo,*res_ptr,*A_ptr[1000],*unknown_atom[100000],*x_hydro[10000];
  
  ofstream xofil,yofil,zofil,vofil,wofil,uofil,test_file;
  ifstream sofil;
  
  ofil<<"*********** "<<inputfile<<" ***********"<<"\n"; 

  /************************************************************/
  /* Initialize variables.                                    */
  p              = (node **) new node *[no_atoms+10];
  mult           = new int[no_atoms+1];
  h_list         = new int[nh];
  missing        = new int[no_atoms+1];
  chain          = new int[no_atoms+10];
  skip_nhydro    = new int[no_atoms+1];
  skip_hydro     = new int[no_atoms+1];
  chk_bond       = (int **) new int *[no_atoms+1];
  poor_bond      = (int **) new int *[no_atoms+1];
  chk_bond_dist  = new float[no_atoms+1];
  poor_bond_dist = new float[no_atoms+1];
  atype          =(char**)  new  char *[no_atoms+1];
  missing_atoms  =(char**)  new char *[no_atoms+1];
  nrot_bond      =(int**)   new int *[no_atoms+1];
  u_nrot_bond    =(int**)   new int *[1000];
  sxyz           =(float**) new float *[no_atoms+1];
  temp_link      =(int**)   new int *[no_atoms+1];
  
  for( i = 0; i <= no_atoms; i++)	{
    *(chk_bond+i)  = (int*)   new int[3]; 
    *(poor_bond+i) = (int*)   new int[3]; 
    *(atype+i)     = (char*) new char[6];
    *(missing_atoms+i) = (char*) new char[6];
    *(nrot_bond+i) = (int*)   new int[3];
    *(sxyz+i)      = (float*) new float[4];
    *(temp_link+i) = (int*)   new int[maxr];
    
    mult[i] = 0;
    missing[i] = -1;
  }
  
  for( i = 0; i < 1000; i++ ) {
    *(u_nrot_bond+i) = (int*) new int[3];
    for( j = 0; j < 3; j++ )
      u_nrot_bond[i][j] = 0;
  }
  
  for( i = 0; i < nh; i++ ) 
    h_list[i] = 0;
  
  
  for( i = 0; i < 32; i++ ) {
    for( j = 0; j < 32; j++ ) {
      for( k = 0; k < 32; k++ )
	grid[i][j][k] = 0;
    }
  }
  
  xmin = xmin-grdlen;
  ymin = ymin-grdlen;
  zmin = zmin-grdlen;
  xmin = xmin-0.2;
  ymin = ymin-0.2;
  zmin = zmin-0.2;
  
  for( i = 0; i <= no_atoms; i++ ) {
    for( j = 0; j < maxr; j++ ) {
      temp_link[i][j] = -1;
    }
  }
  
  for( i = 0; i <= no_atoms; i++ ) {
    for( j = 0; j < 3; j++ ) {
      chk_bond[i][j] = -1;
      poor_bond[i][j] = -1;
    }
  }
  
  for ( i = 0; i < chaincount; i++) 
    oflag[i] = 0;
  /************************************************************/

  /************************************************************/
  // Reading in the Lookup table for distances between atom pairs	
  dist_lookup[0] = '\0';
  strcat(dist_lookup,path);
  strcat(dist_lookup,"/first/lib/dist_lookup.lib");
  sofil.open(dist_lookup,ios::in);
  if( sofil == NULL ) {
    cout<<endl<<"\aDistance Look-up file not found!"<<endl;
    cout<<"Expected location:  Path = " << path <<endl;
    cout<<"Define new Path in class.h"<<endl<<endl<<endl;
    exit(7);
  }
  sofil.getline(line,80);
  sofil.getline(line,80);

  for( k = 0; k < lookup_atoms; k++) {
    for( i = 0; i < 4; i++) {
      sofil >> line;
      for( j = 0; j < lookup_atoms; j++) {
	sofil >> distance[k][j].d[i];
      }
    }
  }

  /************************************************************/
  /* Setting the resolution for the distances.                */
  for( k = 0; k < lookup_atoms; k++) {
    for( j = 0; j < lookup_atoms; j++) {
      factor = (resolution/100)*(distance[k][j].d[2]-distance[k][j].d[1]);
      distance[k][j].d[0] -= 1.1*factor;
      distance[k][j].d[1] -= factor;
      distance[k][j].d[2] += factor;
      distance[k][j].d[3] += 1.1*factor;
    }
  }

  /************************************************************/
  /* Start processing each group.                             */
  /************************************************************/
  tempo = last; //reset the linked list of atom records to the beginning

  while(1) {

    if( strncmp(tempo->r1.field1,"TER",3) == 0 ) {
      if( tempo == start ) {
	loop_flag = 0;
	break;
      }
      tempo = tempo->prior;
      if( tempo == start ) {
	loop_flag = 0;
	break;
      }
    }
    if( loop_flag == 0 )
      break;

    flag17 = 0;
    flag = 0;
    flag8 = 0;
    natom = 0;
    nhydro = 0;
    res_ptr = tempo;

    /************************************************************/
    /* Compare the current atom's residue type to see if it is  */
    /* one of the "standard" residue types stored in            */
    /* /lib/residue.lib                                         */
    for(size_t III = 0; III < res_count; III++) {
      if( strcmp(res_ptr->r1.rname,res[III].rname) == 0 ) {
	flag = 1;
	rnumber = III;
	break;
      }
    }	
    /************************************************************/
    
    /************************************************************/
    /* If the atom belongs to a residue type listed in /lib/    */
    /* residue.lib .                                            */
    /************************************************************/
    if( flag == 1 ) {
      if( prev_occu == 1 ) {
	res[rnumber].occu[0] = 1;
      }
      while(tempo->r1.mod_res_no == res_ptr->r1.mod_res_no && 
	    strcmp(tempo->r1.chain,res_ptr->r1.chain) == 0 ) {

	//------------------------------------------------proces record of KNOWN group
	//--------------------------------------------------for non-hyro
	strcpy( atype[tempo->r1.ri_sr_no], tempo->r1.aname );
	p[tempo->r1.ri_sr_no] = tempo;
	for( i = 1; i < 4; i++) {
	  sxyz[tempo->r1.ri_sr_no][i] = tempo->r1.coord[i];
	}

	/************************************************************/
	/* Initialize the linked-list hash table for searching for  */
	/* nearby atoms in 3D.                                      */
	/************************************************************/
	iix = (int)((tempo->r1.coord[1]-xmin)/grdlen);
	iiy = (int)((tempo->r1.coord[2]-ymin)/grdlen);
	iiz = (int)((tempo->r1.coord[3]-zmin)/grdlen);
	ix = iix%32;
	iy = iiy%32;
	iz = iiz%32;
	chain[tempo->r1.ri_sr_no] = grid[ix][iy][iz];
	grid[ix][iy][iz] = tempo->r1.ri_sr_no;

	if( tempo->r1.aname[2] != 'H' && 
	    tempo->r1.aname[2] != 'D') {
	  // AJR 05.13.02 deal with labeling of known amino acid residue as O' instead of O: relabel
	  if( strcmp(tempo->r1.aname,"  O' ") == 0 && 
	      tempo->r1.field1[0] == 'A') { 
	    strcpy(p[tempo->r1.ri_sr_no]->r1.aname,"  O  ");
	    strcpy(atype[tempo->r1.ri_sr_no],tempo->r1.aname);
	  }
	  name[0] = tempo->r1.aname[2];
	  name[1] = tempo->r1.aname[3];
	  name[2] = tempo->r1.aname[4];
	  name[3] = '\0';
		  
	  for( i = 1; i <= res[rnumber].natom_sum; i++) {
	    if(strcmp( name, res[rnumber].A_type[i]) == 0 && 
	       res[rnumber].occu[i] != 1) {
	      flag17 = 1;
	      A_ptr[i] = tempo;
	      atom_found = 1;   //---------known non-hydrogen atom
	      res[rnumber].occu[i] = 1;
	      natom++;
	      if(tempo->r1.aname[2]=='S') {   //identify sulfur atoms
		nsulfur++;
		if( nsulfur > 1000 ) {
		  cout<<"Number of S atoms greater than 1000"<<endl;
		  exit(-1);
		}
		s_ptr[nsulfur] = tempo;
	      }
	      break;
	    }
	  }
	  
	  if( flag17 != 1 ) {
	    for( i = 1; i <= res[rnumber].natom_sum; i++) {
	      type[0] = res[rnumber].A_type[i][0];
	      type[1] = res[rnumber].A_type[i][1];
	      type[2] = '\0';
	      if( strncmp(name,type,2) == 0 && 
		  res[rnumber].occu[i] != 1) {
		A_ptr[i] = tempo;
		atom_found = 1;   //---------known non-hydrogen atom
		res[rnumber].occu[i] = 1;
		natom++;
		if( tempo->r1.aname[2] == 'S') { //identify sulfur atoms
		  nsulfur++;
		  if( nsulfur > 1000 ) {
		    cout<<"Number of S atoms greater than 1000"<<endl;
		    exit(-1);
		  }
		  s_ptr[nsulfur]=tempo;
		}
		break;
	      }
	    }	
	  }
	  flag17 = 0;

	  if( atom_found != 1 ) { //-------------------extra atoms
	    if( strcmp(tempo->r1.aname,"  OT1") == 0 ||
		strcmp(tempo->r1.aname,"  OT2") == 0 || 
		strcmp(tempo->r1.aname,"  O''") == 0) { // AJR 05.13.02 just rewrite as OXT
	      strcpy(p[tempo->r1.ri_sr_no]->r1.aname,"  OXT");
	    }
	    if(strcmp(tempo->r1.aname,"  OXT")!=0 ) {
	      /*if(strcmp(tempo->r1.aname,"  OXT")!=0 && 
		strcmp(tempo->r1.aname,"  OT1")!=0 && 
		strcmp(tempo->r1.aname,"  OT2")!=0)*/
	      ua_count++;
	      unknown_atom[ua_count]=tempo;
	    }
	    else {
	      n_oxt++;
	      so = tempo->r1.ri_sr_no;
	      for( j = 1; j <= res[rnumber].natom_sum; j++) {
		if( strcmp(res[rnumber].A_type[j],"C  ") == 0 ) {
		  sf = A_ptr[j]->r1.ri_sr_no;
		}
	      }	
	      if( n_oxt > 500 ) {
		cout<<"No. of OXT atoms greater than 500"<<endl;
		exit(12);
	      }
	      bond_oxt[n_oxt][0] = so;
	      bond_oxt[n_oxt][1] = sf;
	    }
	  }
	  atom_found = 0;
	}
	else {
	  nhydro++;
	  h_ptr[nhydro]=tempo;
	}
	if( tempo != start ) {
	  tempo = tempo->prior;
	  if( strcmp(tempo->r1.field1,"TER   ") == 0 ) {
	    if( tempo == start ) {
	      loop_flag = 0;
	      break;
	    }
	    tempo = tempo->prior;
	    if( tempo == start ) {
	      loop_flag = 0;
	      break;
	    }
	  }
	}
	else {
	  loop_flag = 0;
	  break;
	}
      }
//------------------End of group - Now work with known group------------------//
			
//--------------------------------------connectivity table for NON-Hydrogen atoms
      for( i = 1; i <= res[rnumber].natom_sum; i++) {
	if( res[rnumber].occu[i] == 1 ) {
	  for( j = 1; j <= res[rnumber].A_mult[i]; j++) {
	    if( res[rnumber].occu[res[rnumber].A_link[i][j]] == 1 ) {
	      so = A_ptr[i]->r1.ri_sr_no;
	      sf = A_ptr[res[rnumber].A_link[i][j]]->r1.ri_sr_no;
	      temp = 0;
	      for( k = 1; k < 4; k++) {
		temp += pow((sxyz[so][k]-sxyz[sf][k]),2);
	      }
	      temp = sqrt(temp);

	      // check distance criteria before connecting 2 residues
	      if( i == 1 ) {
		if( strcmp(p[so]->r1.chain, p[sf]->r1.chain) == 0 ) {
		  ret_val=check_dist(so,sf,temp,atype[so],atype[sf],1,p[so]->r1.rname,
				     p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,
				     p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
				     p[sf]->r1.chain,usage);
		  if( ret_val == 1 ) {
		    make_connection(mult,temp_link,so,sf);	
		    // record N=C Non-rotating bond between consecutive residues
		    no_bonds += 1;
		    nrot_bond_count++;
		    nrot_bond[nrot_bond_count][1]=so;
		    nrot_bond[nrot_bond_count][2]=sf;
		  }
		  else {
		    poor_bond_count++;
		    poor_bond[poor_bond_count][1]=so;
		    poor_bond[poor_bond_count][2]=sf;
		    poor_bond_dist[poor_bond_count]=temp;
		  }
		}
		else {
		  ret_val=check_dist(so,sf,temp,atype[so],atype[sf],3,p[so]->r1.rname,
				     p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,
				     p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
				     p[sf]->r1.chain,usage);
		  if( ret_val == 1 ) {
		    make_connection(mult,temp_link,so,sf);	
		    // record N=C Non-rotating bond between consecutive residues
		    no_bonds += 1;
		    nrot_bond_count++;
		    nrot_bond[nrot_bond_count][1]=so;
		    nrot_bond[nrot_bond_count][2]=sf;
		  }
		  else {
		    poor_bond_count++;
		    poor_bond[poor_bond_count][1]=so;
		    poor_bond[poor_bond_count][2]=sf;
		    poor_bond_dist[poor_bond_count]=temp;
		  }
		}
		
	      }	
	      else {
		ret_val=check_dist(so,sf,temp,atype[so],atype[sf],1,p[so]->r1.rname,
				   p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,
				   p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
				   p[sf]->r1.chain,usage);
		
		if( ret_val == 1 ) {
		  // Special case: ARG for non-rotating bonds
		  if( strcmp( res[rnumber].rname, "ARG" ) == 0 ) {
		    if( strcmp(atype[so],"  CZ ") == 0 ) {
		      if(strcmp(atype[sf],"  NH1") == 0 || strcmp(atype[sf],"  NH2") == 0 || 
			 strcmp(atype[sf],"  NE ") == 0) {
			nrot_bond_count++;
			nrot_bond[nrot_bond_count][1] = so;
			nrot_bond[nrot_bond_count][2] = sf;
		      }
		    }
		    if( strcmp(atype[sf],"  CZ ") == 0 ) {
		      if(strcmp(atype[so],"  NH1") == 0 || 
			 strcmp(atype[so],"  NH2") == 0 || 
			 strcmp(atype[so],"  NE ") == 0) {
			nrot_bond_count++;
			nrot_bond[nrot_bond_count][1]=so;
			nrot_bond[nrot_bond_count][2]=sf;
		      }
		    }
		  }
		  // End Special case: ARG for non-rotating bonds
			
		  make_connection( mult, temp_link, so, sf );	
		  no_bonds += 1;
		}
		else {
		  poor_bond_count++;
		  poor_bond[poor_bond_count][1] = so;
		  poor_bond[poor_bond_count][2]=sf;
		  poor_bond_dist[poor_bond_count]=temp;
		}
	      }	
	    }
	  }
	}
	else {
	  ma_count++;	//---------------missing atoms
	  missing[ma_count] = res_ptr->r1.ri_sr_no;
	  strcpy( missing_atoms[ma_count], res[rnumber].A_type[i] );
	}
      }
      
      /************************************************************/
      /* Place Hydrogen atoms.                                    */
      /************************************************************/
      for( i = 1; i <= nhydro; i++ ) {
	mult[h_ptr[i]->r1.ri_sr_no] = 0;
      }
      for( i = 1; i <= nhydro; i++) {
	flag9 = 0;
	so = h_ptr[i]->r1.ri_sr_no;
	count = 0;
	for( j = 1; j <= res[rnumber].natom_sum; j++) {
	  if( res[rnumber].occu[j] == 1 ) {
	    sf = A_ptr[j]->r1.ri_sr_no;
	    temp = 0;
	    for( k = 1; k < 4; k++) {
	      temp += pow((h_ptr[i]->r1.coord[k] - A_ptr[j]->r1.coord[k]),2);
	    }
	    temp = sqrt(temp);

// 2008:10 SN	Add a diagnostic message to the user if bond length < 1.0
	    if( temp >= 0.85 && temp < 0.995  && \
		strncmp(A_ptr[j]->r1.aname,"  C",3)==0) {
		cout<< "\nWARNING: Distance between "<<so<<" and "<<sf<<" is";
		cout<< temp<<" A (between 0.85 and 1.0 A)!";
		cout<< "\n         ProFlex will still include the hydrogen";
		cout<< " atom, "<<so<<" in its analysis."<<endl;
	    }


	    ret_val=check_dist(so,sf,temp,atype[so],atype[sf],2,p[so]->r1.rname,
			       p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,
			       p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
			       p[sf]->r1.chain,usage);

	    if( ret_val == 1 ) {
	      count++;
	      dist[count]=temp;
	      label[count]=sf;
	    }
	  }
	}
	for( j = 1; j <= 2; j++) {
	  for( k = j+1; k <= count; k++) {
	    if( dist[j] > dist[k] ) {
	      temp = dist[j];
	      dist[j] = dist[k];
	      dist[k] = temp;
	      a = label[j];
	      label[j] = label[k];
	      label[k] = a;
	    }
	  }
	}
	sf = label[1];
	if( count > 1 ) {
	  flag9 = 1;			
	}
	else {
	  if( count == 1 ) {
	    make_connection( mult, temp_link, so, sf );	
	    no_bonds++;
	  }
	}
	if( flag9 == 1 || 
	    count == 0) {
	  xh_count++;
	  x_hydro[xh_count] = h_ptr[i];
	}
      }
      if( (res[rnumber].occu[3]) == 1 ) {		
	A_ptr[0] = A_ptr[3];
	prev_occu=1;
      }
      for( i = 0; i <= res[rnumber].natom_sum; i++) {
	res[rnumber].occu[i]=-1;
      }
    }
    else { // Process record of UNKNOWN group
      while( tempo->r1.mod_res_no == res_ptr->r1.mod_res_no && 
	    strcmp( tempo->r1.chain,res_ptr->r1.chain) == 0 ) {
	p[tempo->r1.ri_sr_no] = tempo;
	strcpy(atype[tempo->r1.ri_sr_no],tempo->r1.aname);
	for( i = 1; i < 4; i++ ) {
	  sxyz[tempo->r1.ri_sr_no][i]=tempo->r1.coord[i];
	}
	iix=(int)((tempo->r1.coord[1]-xmin)/grdlen);
	iiy=(int)((tempo->r1.coord[2]-ymin)/grdlen);
	iiz=(int)((tempo->r1.coord[3]-zmin)/grdlen);
	ix=iix%32;
	iy=iiy%32;
	iz=iiz%32;
	chain[tempo->r1.ri_sr_no]=grid[ix][iy][iz];
	grid[ix][iy][iz]=tempo->r1.ri_sr_no;
	if(tempo->r1.aname[2]=='H' || tempo->r1.aname[2]=='D') {
	  nhydro++;
	  h_ptr[nhydro]=tempo;
	  xh_count++;
	  x_hydro[xh_count]=tempo;
	}
	else {
	  natom++;
	  A_ptr[natom]=tempo;
	  ua_count++;
	  unknown_atom[ua_count]=tempo;
	}
	if( tempo != start ) {
	  tempo=tempo->prior;
	}
	else {
	  loop_flag = 0;
	  break;
	}
      }
      for( i = 1; i <= natom; i++ ) {
	if(A_ptr[i]->r1.aname[2] == 'S' ) {
	  nsulfur++;
	  if( nsulfur > 1000 ) {
	    cout<<"Number of S atoms greater than 1000"<<endl;
	    exit(-1);
	  }
	  s_ptr[nsulfur]=A_ptr[i];
	}
      }
      prev_occu=-1;
    }
  }
  
  // Unknown Non-hydrogen atom connectivity using "Hash code data structure"
  for( i = 1; i <= ua_count; i++) {
    so=unknown_atom[i]->r1.ri_sr_no;
    iix=(int)((unknown_atom[i]->r1.coord[1]-xmin)/grdlen);
    iiy=(int)((unknown_atom[i]->r1.coord[2]-ymin)/grdlen);
    iiz=(int)((unknown_atom[i]->r1.coord[3]-zmin)/grdlen);
    
    //---------------------------------------go to the beginning of the residue
    flag9 = 0;
    for( jx = -1; jx <= 1; jx++) {
      for( jy = -1; jy <= 1; jy++) {
	for( jz = -1; jz <= 1; jz++) {
	  kx = (iix+jx)%32;
	  ky = (iiy+jy)%32;
	  kz = (iiz+jz)%32;
		  
	  sf = grid[kx][ky][kz];
		  
	  while( sf != 0 ) {
	    if( atype[sf][2] != 'H' && 
		atype[sf][2] != 'D') {
	      for( j = 1; j <= mult[so]; j++) {
		if( temp_link[so][j] == sf ) {
		  flag12 = 1;
		  flag9  = 1;
		  break;
		}
	      }
	      if( flag12 != 1 ) {
		temp = 0;
		for( k = 1; k < 4; k++) {
		  temp += pow((sxyz[so][k]-sxyz[sf][k]),2);
		}
		temp = sqrt(temp);
		ret_val=check_dist(so,sf,temp,atype[so],atype[sf],0,p[so]->r1.rname,
				   p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,
				   p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
				   p[sf]->r1.chain,usage);
		
		if( ret_val == 1 ) {
		  /************************************************************/
		  /* Check for non rotating  Note: Bonds are non rotating if  */
		  /* the 2 atoms belong to diff. groups                       */
		  if(p[so]->r1.mod_res_no != p[sf]->r1.mod_res_no) {
		    u_nrot_bond_count++;
		    u_nrot_bond[u_nrot_bond_count][1] = so;
		    u_nrot_bond[u_nrot_bond_count][2] = sf;
		  }
		  make_connection(mult,temp_link,so,sf);	
		  no_bonds+=1;
		  flag12=0;
		  chk_bond_count++;
#ifdef DEBUG
printf("1. chk_bond_count =%d so %d sf %d dist %f\n",chk_bond_count,so,sf,temp);
#endif
		  chk_bond[chk_bond_count][1]=so;
		  chk_bond[chk_bond_count][2]=sf;
		  chk_bond_dist[chk_bond_count]=temp;
		  
		  flag9=1;
		  sf=chain[sf];
		}
		else {
		  sf=chain[sf];
		}
	      }
	      else {
		sf=chain[sf];
		flag12=0;	
	      }
	    }
	    else {
	      sf=chain[sf];
	    }
	  }
	}
      }
    }
      
    if( flag9 != 1 ) {
      skip_nhydro_count++;
      skip_nhydro[skip_nhydro_count]=so;
    }
  }			
  
  // Sulfur-Sulfur atom connectivity using "Hash code data structure"
  ns_bond = 0;
  for( i = 1; i <= nsulfur; i++ ) {
    so = s_ptr[i]->r1.ri_sr_no;
    iix=(int)((s_ptr[i]->r1.coord[1]-xmin)/grdlen);
    iiy=(int)((s_ptr[i]->r1.coord[2]-ymin)/grdlen);
    iiz=(int)((s_ptr[i]->r1.coord[3]-zmin)/grdlen);
    
    //---------------------------------------go to the beginning of the residue
    for( jx = -1; jx <= 1; jx++) {
      for( jy = -1; jy <= 1; jy++) {
	for( jz = -1; jz <= 1; jz++) {
	  kx = (iix+jx)%32;
	  ky = (iiy+jy)%32;
	  kz = (iiz+jz)%32;
		  
	  sf = grid[kx][ky][kz];
	  while( sf != 0 ) {
	    if( atype[sf][2] == 'S' && 
	        so != sf) {
	      for( j = 1; j <= mult[so]; j++) {	
		if( temp_link[so][j] == sf ) {
		  flag12 = 1;
		  break;
		}
	      }
	      if( flag12 != 1 ) {
		temp = 0;
		for( k = 1; k < 4; k++) {
		  temp+=pow((sxyz[so][k]-sxyz[sf][k]),2);
		}
		temp=sqrt(temp);
		ret_val=check_dist(so,sf,temp,atype[so],atype[sf],0,p[so]->r1.rname,
				   p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,
				   p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
				   p[sf]->r1.chain,usage);
		
		if( ret_val == 1 ) {
		  /************************************************************/
		  /* Check for non rotating. Note: Bonds are non rotating if  */
		  /* the 2 atoms belong to diff. groups                       */
		  /************************************************************/
		  if(p[so]->r1.mod_res_no != p[sf]->r1.mod_res_no) {
		    u_nrot_bond_count++;
		    u_nrot_bond[u_nrot_bond_count][1]=so;
		    u_nrot_bond[u_nrot_bond_count][2]=sf;
		  }
		  make_connection(mult,temp_link,so,sf);	
		  no_bonds+=1;
		  ns_bond++;
		  flag12=0;
		  chk_bond_count++;
#ifdef DEBUG
printf("1. chk_bond_count =%d so %d sf %d dist %f\n",chk_bond_count,so,sf,temp);
#endif
		  chk_bond[chk_bond_count][1]=so;
		  chk_bond[chk_bond_count][2]=sf;
		  chk_bond_dist[chk_bond_count]=temp;
		  sf=chain[sf];
		}
		else {
		  sf = chain[sf];
		}
	      }
	      else {
		sf = chain[sf];
		flag12 = 0;
	      }		
	    }
	    else {
	      sf = chain[sf];
	    }
	  }
	}
      }
    }
  }
  
  // Extra & unknown hydrogen atom connectivity using "Hash code data structure"
  for(i=1;i<=xh_count;i++)
    {
      so=x_hydro[i]->r1.ri_sr_no;

      iix=(int)((x_hydro[i]->r1.coord[1]-xmin)/grdlen);
      iiy=(int)((x_hydro[i]->r1.coord[2]-ymin)/grdlen);
      iiz=(int)((x_hydro[i]->r1.coord[3]-zmin)/grdlen);
      count=0;
      flag9=0;
      for(jx=-1;jx<=1;jx++)
	{
	  for(jy=-1;jy<=1;jy++)
	    {
	      for(jz=-1;jz<=1;jz++)
		{
		  kx=(iix+jx)%32;
		  ky=(iiy+jy)%32;
		  kz=(iiz+jz)%32;
		  
		  sf=grid[kx][ky][kz];
		  
		  while(sf!=0)
		    {
		      if(atype[sf][2]!='H' && atype[sf][2]!='D')
			{
			  temp=0;
			  for(k=1;k<4;k++)
			    {
			      temp+=pow((sxyz[so][k]-sxyz[sf][k]),2);
			    }
			  temp=sqrt(temp);
       ret_val=check_dist(so,sf,temp,atype[so],atype[sf],0,p[so]->r1.rname,
			  p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,
			  p[sf]->r1.rname,p[sf]->r1.res_no,p[sf]->r1.code,
			  p[sf]->r1.chain,usage);

                          if(ret_val==1) {
			    count++;
			    dist[count]=temp;
			    label[count]=sf;
			  }
			  sf=chain[sf];
			}
		      else
			{
			  sf=chain[sf];
			}
		    }
		  
		}
	    }
	}
      for(j=1;j<=2;j++)
	{
	  for(k=j+1;k<=count;k++)
	    {
	      if(dist[j]>dist[k])
		{
		  temp=dist[j];
		  dist[j]=dist[k];
		  dist[k]=temp;
		  a=label[j];
		  label[j]=label[k];
		  label[k]=a;
		}
	    }
	}
      sf=label[1];
      if(count>1)
	{
	  {
	    if(flag7==0)
	      {
                cout<<endl<<dash1<<dash3<<endl;
		cout<<endl<<"\tHydrogen Atom Connectivity"<<endl;	
		cout<<endl<<"Warning: The following hydrogen atoms can be" 
                    << " connected to more than 1 atom"<<endl;
		cout<<endl<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
		flag7=1;
	      }	
	write_output(1,so,0,atype[so],p[so]->r1.rname,p[so]->r1.res_no,
		     p[so]->r1.code,p[so]->r1.chain,atype[so],atype[so],
		     atype[so],atype[so],atype[so]);
	    flag9=1;
	  }
	}
      else 
	{
	  if(count==1)
	    {
       /*---------------------------------------------check for non rotating 
	 Note: Bonds are non rotating if the 2 atoms belong to diff. groups  */

	      if(p[so]->r1.mod_res_no != p[sf]->r1.mod_res_no)
		{
		  u_nrot_bond_count++;
		  u_nrot_bond[u_nrot_bond_count][1]=so;
		  u_nrot_bond[u_nrot_bond_count][2]=sf;
		}
	      make_connection(mult,temp_link,so,sf);	
	      no_bonds+=1;
	    }
	}
      if(flag9==1 || count==0)
	{
	  skip_hydro_count++;
	  skip_hydro[skip_hydro_count]=so;
	}
    }			  
//--------------------------------------------------------------Add in OXT atoms
	
  for( i = 1; i <= n_oxt; i++) {
    so = bond_oxt[i][0];
    sf = bond_oxt[i][1];
    make_connection(mult,temp_link,so,sf);	
    no_bonds++;
  }
  
  if(usage == 0) {
//----------------------------------------------------------------Check valencies
  flag7=0;
  for(i=1;i<=no_atoms;i++)
    {
      flag10=0;
      flag9=0;
      if(atype[i][2]=='H' || atype[i][2]=='D')
	{
	  if(mult[i]>1)
	    {
	      flag10 = -1;
	      so=i;
	      if(flag7==0)
		{
		  //system("clear");
		  cout<<endl<<endl<<"\t\t    Checking valencies"<<endl;	
		  cout<<dash1<<dash3<<endl;
		  flag7=1;
		}	
	      cout<<endl<<"Valency of atom "<<i<<" is "<<mult[i]
		  <<", greater than its maximum valency (1)"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	      write_output(1,so,0,atype[so],p[so]->r1.rname,p[so]->r1.res_no,
			   p[so]->r1.code,p[so]->r1.chain,atype[so],atype[so],
			   atype[so],atype[so],atype[so]);
	      for(j=1;j<=mult[i];j++)
		{
		  if(strcmp(p[i]->r1.res_no,p[temp_link[i][j]]->r1.res_no)!=0)
		    {
		      if(flag9==0)
			{
		    cout<<endl<<"Possible incorrect bond(s) with atom(s):"<<endl;
			  flag9=1;
			}
		      sf=temp_link[i][j];
		      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
		      write_output(1,sf,0,atype[sf],p[sf]->r1.rname,p[sf]->r1.res_no,
				   p[sf]->r1.code,p[sf]->r1.chain,atype[sf],atype[sf],
				   atype[sf],atype[sf],atype[sf]);
		      ans = '@';
		      while( ans != 'Y' && 
			     ans != 'y' && 
			     ans != 'N' && 
			     ans != 'n') {
			cout<<"Do you want to keep the bond? (y/n): ";
			ans = cin.get();
			if( ans == '\n' )
			  cout << "Please enter \"y\" or \"n\"." << endl;
			else
			  cin.ignore(80, '\n');
		      }	

		      if(ans=='N' || ans=='n')
			{
			  for(k=j;k<=mult[so]-1;k++)
			    {
			      temp_link[so][k]=temp_link[so][k+1];
			    }
			  temp_link[so][mult[so]]=-1;
			  mult[so]--;
			  for(m=1;m<=mult[sf];m++)
			    {
			      if(temp_link[sf][m]==so)
				{
				  for(q=m;q<=mult[sf]-1;q++)
				    {
				      temp_link[sf][q]=temp_link[sf][q+1];
				    }
				  temp_link[sf][mult[sf]]=-1;
				}	
			    }
			  mult[sf]--;
			  no_bonds--;
			  for(k=1;k<=chk_bond_count;k++)
			    {
			      if(chk_bond[k][1]==so || chk_bond[k][1]==sf)
				{
				  chk_bond[k][1]=-1;
				}
			    }
			  for(k=1;k<=u_nrot_bond_count;k++)
			    {
			      if(u_nrot_bond[k][1]==so || u_nrot_bond[k][1]==sf)
				{
				  u_nrot_bond[k][1]=0;
				}
			    }
			  flag10=1;
			  j--;
			}
		    }
		}
	      if( mult[i] > 1 ) 
		{
		  for(j=1;j<=mult[i];j++)
		    {
		     if(strcmp(p[i]->r1.res_no,p[temp_link[i][j]]->r1.res_no)==0)
		       {
			 if(flag9==0)
			   {
		    cout<<endl<<"Possible incorrect bond(s) with atom(s):"<<endl;
			     flag9=1;
			   }
			 sf=temp_link[i][j];
			 cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
		write_output(1,sf,0,atype[sf],p[sf]->r1.rname,p[sf]->r1.res_no,
			     p[sf]->r1.code,p[sf]->r1.chain,atype[sf],atype[sf],
			     atype[sf],atype[sf],atype[sf]);
			 ans = '@';
			 while( ans != 'Y' && 
				ans != 'y' && 
				ans != 'N' && 
				ans != 'n') {
			   cout<<"Do you want to keep the bond? (y/n): ";
			   ans = cin.get();
			   if( ans == '\n' )
			     cout << "Please enter \"y\" or \"n\"." << endl;
			   else
			     cin.ignore(80, '\n');
			 }	
			 if(ans=='N' || ans=='n')
			   {
			     for(k=j;k<=mult[so]-1;k++)
			       {
				 temp_link[so][k]=temp_link[so][k+1];
			       }
			     temp_link[so][mult[so]]=-1;
			     mult[so]--;
			     for(m=1;m<=mult[sf];m++)
			       {
				 if(temp_link[sf][m]==so)
				   {
				     for(q=m;q<=mult[sf]-1;q++)
				       {
					 temp_link[sf][q]=temp_link[sf][q+1];
				       }
				     temp_link[sf][mult[sf]]=-1;
				   }	
			       }
			     mult[sf]--;
			     no_bonds--;
			     for(k=1;k<=chk_bond_count;k++)
			       {
				 if(chk_bond[k][1]==so || chk_bond[k][1]==sf)
				   {
				     chk_bond[k][1]=-1;
				   }
			       }
			     for(k=1;k<=u_nrot_bond_count;k++)
			       {
			      if(u_nrot_bond[k][1]==so || u_nrot_bond[k][1]==sf)
				   {
				     u_nrot_bond[k][1]=0;
				   }
			       }
			     flag10=1;
			     j--;
			   }
		       }
		    }
		}
	    }
	}
//      else if(atype[i][2]=='O' && atype[i][3] == ' ')		
//	--- 2007:06 SN	---	Discrepancy in valency checks between v4 and v5
      else if(atype[i][2]=='O')
	{		
	  if(mult[i]>2)
	    {
	      flag10 = -1;
	      so=i;
	      if(flag7==0)
		{
		  //system("clear");
		  cout<<endl<<endl<<"\t\t    Checking valencies"<<endl;	
		  cout<<dash1<<dash3<<endl;
		  flag7=1;
		}	
	      cout<<endl<<"Valency of atom "<<i<<" is "<<mult[i]
		  <<", greater than its maximum valency (2)"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	      write_output(1,so,0,atype[so],p[so]->r1.rname,p[so]->r1.res_no,
			   p[so]->r1.code,p[so]->r1.chain,atype[so],atype[so],
			   atype[so],atype[so],atype[so]);
	      for(j=1;j<=mult[i];j++)
		{
		  if(strcmp(p[i]->r1.res_no,p[temp_link[i][j]]->r1.res_no)!=0)
		    {
		      if(flag9==0)
			{
		    cout<<endl<<"Possible incorrect bond(s) with atom(s):"<<endl;
			  flag9=1;
			}
		      sf=temp_link[i][j];
		      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
		write_output(1,sf,0,atype[sf],p[sf]->r1.rname,p[sf]->r1.res_no,
			     p[sf]->r1.code,p[sf]->r1.chain,atype[sf],atype[sf],
			     atype[sf],atype[sf],atype[sf]);
		      ans = '@';
		      while( ans != 'Y' && 
			     ans != 'y' && 
			     ans != 'N' && 
			     ans != 'n') {
			cout<<"Do you want to keep the bond? (y/n): ";
			ans = cin.get();
			if( ans == '\n' )
			  cout << "Please enter \"y\" or \"n\"." << endl;
			else
			  cin.ignore(80, '\n');
		      }	
		      if(ans=='N' || ans=='n')
			{
			  for(k=j;k<=mult[so]-1;k++)
			    {
			      temp_link[i][k]=temp_link[i][k+1];
			    }
			  temp_link[i][mult[i]]=-1;
			  mult[i]--;
			  for(m=1;m<=mult[sf];m++)
			    {
			      if(temp_link[sf][m]==so)
				{
				  for(q=m;q<=mult[sf]-1;q++)
				    {
				      temp_link[sf][q]=temp_link[sf][q+1];
				    }
				  temp_link[sf][mult[sf]]=-1;
				}	
			    }
			  mult[sf]--;
			  no_bonds--;
			  for(k=1;k<=chk_bond_count;k++)
			    {
			      if(chk_bond[k][1]==so || chk_bond[k][1]==sf)
				{
				  chk_bond[k][1]=-1;
				}
			    }
			  for(k=1;k<=u_nrot_bond_count;k++)
			    {
			      if(u_nrot_bond[k][1]==so || u_nrot_bond[k][1]==sf)
				{
				  u_nrot_bond[k][1]=0;
				}
			    }
			  flag10=1;
			  j--;
			}
		    }
		}
	      if( mult[i] > 2 ) 
		{
		  for(j=1;j<=mult[i];j++)
		    {
		     if(strcmp(p[i]->r1.res_no,p[temp_link[i][j]]->r1.res_no)==0)
			{
			  if(flag9==0)
			    {
		    cout<<endl<<"Possible incorrect bond(s) with atom(s):"<<endl;
			      flag9=1;
			    }
			  sf=temp_link[i][j];
			  cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
		write_output(1,sf,0,atype[sf],p[sf]->r1.rname,p[sf]->r1.res_no,
			     p[sf]->r1.code,p[sf]->r1.chain,atype[sf],atype[sf],
			     atype[sf],atype[sf],atype[sf]);
			  ans = '@';
			  while( ans != 'Y' && 
				 ans != 'y' && 
				 ans != 'N' && 
				 ans != 'n') {
			    cout<<"Do you want to keep the bond? (y/n): ";
			    ans = cin.get();
			    if( ans == '\n' )
			      cout << "Please enter \"y\" or \"n\"." << endl;
			    else
			      cin.ignore(80, '\n');
			  }	
			  if(ans=='N' || ans=='n')
			    {
			      for(k=j;k<=mult[so]-1;k++)
				{
				  temp_link[i][k]=temp_link[i][k+1];
				}
			      temp_link[i][mult[i]]=-1;
			      mult[i]--;
			      for(m=1;m<=mult[sf];m++)
				{
				  if(temp_link[sf][m]==so)
				    {
				      for(q=m;q<=mult[sf]-1;q++)
					{
					  temp_link[sf][q]=temp_link[sf][q+1];
					}
				      temp_link[sf][mult[sf]]=-1;
				    }	
				}
			      mult[sf]--;
			      no_bonds--;
			      for(k=1;k<=chk_bond_count;k++)
				{
				  if(chk_bond[k][1]==so || chk_bond[k][1]==sf)
				    {
				      chk_bond[k][1]=-1;
				    }
				}
			      for(k=1;k<=u_nrot_bond_count;k++)
				{
			       if(u_nrot_bond[k][1]==so || u_nrot_bond[k][1]==sf)
				 {
				   u_nrot_bond[k][1]=0;
				 }
				}
			      flag10=1;
			      j--;
			    }
			}
		    }
		}
	    }
	}
//      else if( (atype[i][2]=='N' && atype[i][3] == ' ') && 
//	--- 2007:06 SN	---	Discrepancy in valency checks between v4 and v5
      else if( atype[i][2]=='N' && 
	       ( p[i]->r1.mod_res_no != 1 || ( p[i]->r1.mod_res_no == 1 && 
		 strncmp(p[i]->r1.aname+1," N  ",4)!=0 ) )  &&
	         !(strcmp(p[i]->r1.rname,"LYS")==0 && 
		 strncmp(p[i]->r1.aname+1," NZ ",4)==0)   )
	{		
	  if(mult[i]>3)
	    {
	      flag10 = -1;
	      so=i;
	      if(flag7==0)
		{
		  //system("clear");
		  cout<<endl<<endl<<"\t\t    Checking valencies"<<endl;	
		  cout<<dash1<<dash3<<endl;
       		  flag7=1;
		}	
	      cout<<endl<<"Valency of atom "<<i<<" is "<<mult[i]
		  <<", greater than its maximum valency (3)"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	      write_output(1,so,0,atype[so],p[so]->r1.rname,p[so]->r1.res_no,
			   p[so]->r1.code,p[so]->r1.chain,atype[so],atype[so],
			   atype[so],atype[so],atype[so]);
	      for(j=1;j<=mult[i];j++)
		{
		  if(strcmp(p[i]->r1.res_no,p[temp_link[i][j]]->r1.res_no)!=0)
		    {
		      if(flag9==0)
			{
		    cout<<endl<<"Possible incorrect bond(s) with atom(s):"<<endl;
		    flag9=1;
			}
		      sf=temp_link[i][j];
		      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
		write_output(1,sf,0,atype[sf],p[sf]->r1.rname,p[sf]->r1.res_no,
			     p[sf]->r1.code,p[sf]->r1.chain,atype[sf],atype[sf],
			     atype[sf],atype[sf],atype[sf]);
		      ans = '@';
		      while( ans != 'Y' && 
			     ans != 'y' && 
			     ans != 'N' && 
			     ans != 'n') {
			cout<<"Do you want to keep the bond? (y/n): ";
			ans = cin.get();
			if( ans == '\n' )
			  cout << "Please enter \"y\" or \"n\"." << endl;
			else
			  cin.ignore(80, '\n');
		      }	
		      if(ans=='N' || ans=='n')
			{
			  for(k=j;k<=mult[so]-1;k++)
			    {
			      temp_link[i][k]=temp_link[i][k+1];
			    }
			  temp_link[i][mult[i]]=-1;
			  mult[i]--;
			  for(m=1;m<=mult[sf];m++)
			    {
			      if(temp_link[sf][m]==so)
				{
				  for(q=m;q<=mult[sf]-1;q++)
				    {
				      temp_link[sf][q]=temp_link[sf][q+1];
				    }
				  temp_link[sf][mult[sf]]=-1;
				}	
			    }
			  mult[sf]--;
			  no_bonds--;
			  for(k=1;k<=chk_bond_count;k++)
			    {
			      if(chk_bond[k][1]==so || chk_bond[k][1]==sf)
				{
				  chk_bond[k][1]=-1;
				}
			    }
			  for(k=1;k<=u_nrot_bond_count;k++)
			    {
			      if(u_nrot_bond[k][1]==so || u_nrot_bond[k][1]==sf)
				{
				  u_nrot_bond[k][1]=0;
				}
			    }
			  flag10=1;
			  j--;
			}
		    }
		}
	      if( mult[i] > 3 ) 
		{
		  for(j=1;j<=mult[i];j++)
		    {
		     if(strcmp(p[i]->r1.res_no,p[temp_link[i][j]]->r1.res_no)==0)
		       {
			 if(flag9==0)
			   {
		    cout<<endl<<"Possible incorrect bond(s) with atom(s):"<<endl;
		    flag9=1;
			   }
			 sf=temp_link[i][j];
			 cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
		write_output(1,sf,0,atype[sf],p[sf]->r1.rname,p[sf]->r1.res_no,
			     p[sf]->r1.code,p[sf]->r1.chain,atype[sf],atype[sf],
			     atype[sf],atype[sf],atype[sf]);
			 ans = '@';
			 while( ans != 'Y' && 
				ans != 'y' && 
				ans != 'N' && 
				ans != 'n') {
			   cout<<"Do you want to keep the bond? (y/n): ";
			   ans = cin.get();
			   if( ans == '\n' )
			     cout << "Please enter \"y\" or \"n\"." << endl;
			   else
			     cin.ignore(80, '\n');
			 }	
			 if(ans=='N' || ans=='n')
			   {
			     for(k=j;k<=mult[so]-1;k++)
			       {
				 temp_link[i][k]=temp_link[i][k+1];
			       }
			     temp_link[i][mult[i]]=-1;
			     mult[i]--;
			     for(m=1;m<=mult[sf];m++)
			       {
				 if(temp_link[sf][m]==so)
				   {
				     for(q=m;q<=mult[sf]-1;q++)
				       {
					 temp_link[sf][q]=temp_link[sf][q+1];
				       }
				     temp_link[sf][mult[sf]]=-1;
				   }	
			       }
			     mult[sf]--;
			     no_bonds--;
			     for(k=1;k<=chk_bond_count;k++)
			       {
				 if(chk_bond[k][1]==so || chk_bond[k][1]==sf)
				   {
				     chk_bond[k][1]=-1;
				   }
			       }
			     for(k=1;k<=u_nrot_bond_count;k++)
			       {
			       if(u_nrot_bond[k][1]==so || u_nrot_bond[k][1]==sf)
				 {
				   u_nrot_bond[k][1]=0;
				 }
			       }
			     flag10=1;
			     j--;
			   }
		       }
		    }
		}
	    }
	}
      else /* If Not 'C' 'N' 'O' 'H' Then ( metals! ligand! )*/
	{		
	  if(mult[i]>4)
	    {
	      flag10 = -1;
	      so=i;
	      if(flag7==0)
		{
		  //system("clear");
		  cout<<endl<<endl<<"\t\t    Checking valencies"<<endl;	
		  cout<<dash1<<dash3<<endl;
		  flag7=1;
		}	
	      cout<<endl<<"Valency of atom "<<i<<" is "<<mult[i]
		  <<", greater than its maximum valency (4)"<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	      write_output(1,so,0,atype[so],p[so]->r1.rname,p[so]->r1.res_no,
			   p[so]->r1.code,p[so]->r1.chain,atype[so],atype[so],
			   atype[so],atype[so],atype[so]);
	      for(j=1;j<=mult[i];j++)
		{
		  if(strcmp(p[i]->r1.res_no,p[temp_link[i][j]]->r1.res_no)!=0)
		    {
		      if(flag9==0) {
			cout<<endl<<"Possible incorrect bond(s) with atom(s):"<<endl;
			flag9=1;
		      }
		      sf=temp_link[i][j];
		      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
		      write_output(1,sf,0,atype[sf],p[sf]->r1.rname,p[sf]->r1.res_no,
				   p[sf]->r1.code,p[sf]->r1.chain,atype[sf],atype[sf],
				   atype[sf],atype[sf],atype[sf]);
		      ans = '@';
		      while( ans != 'Y' && 
			     ans != 'y' && 
			     ans != 'N' && 
			     ans != 'n') {
			cout<<"Do you want to keep the bond? (y/n): ";
			ans = cin.get();
			if( ans == '\n' )
			  cout << "Please enter \"y\" or \"n\"." << endl;
			else
			  cin.ignore(80, '\n');
		      }	
		      if( ans == 'N' || ans == 'n' )
			{
			  for(k=j;k<=mult[i]-1;k++)
			    {
			      temp_link[i][k]=temp_link[i][k+1];
			    }
			  temp_link[i][mult[i]]=-1;
			  mult[i]--;
			  for(m=1;m<=mult[sf];m++)
			    {
			      if(temp_link[sf][m]==so)
				{
				  for(q=m;q<=mult[sf]-1;q++)
				    {
				      temp_link[sf][q]=temp_link[sf][q+1];
				    }
				  temp_link[sf][mult[sf]]=-1;
				}	
			    }
			  mult[sf]--;
			  no_bonds--;
			  for(k=1;k<=chk_bond_count;k++)
			    {
			      if(chk_bond[k][1]==so || chk_bond[k][1]==sf)
				{
				  chk_bond[k][1]=-1;
				}
			    }
			  for(k=1;k<=u_nrot_bond_count;k++)
			    {
			      if(u_nrot_bond[k][1]==so || u_nrot_bond[k][1]==sf)
				{
				  u_nrot_bond[k][1]=0;
				}
			    }
			  flag10=1;
			  j--;
			}
		    }
		}
	      if( mult[i] > 4 )  
		{
		  for(j=1;j<=mult[i];j++)
		    {
		     if(strcmp(p[i]->r1.res_no,p[temp_link[i][j]]->r1.res_no)==0)
		       {
			 if( flag9 == 0) {
			   cout<<endl<<"Possible incorrect bond(s) with atom(s):"<<endl;
			   flag9=1;
			 }
			 sf=temp_link[i][j];
			 cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
			 write_output(1,sf,0,atype[sf],p[sf]->r1.rname,p[sf]->r1.res_no,
				      p[sf]->r1.code,p[sf]->r1.chain,atype[sf],atype[sf],
				      atype[sf],atype[sf],atype[sf]);
			 ans = '@';
			 while( ans != 'Y' && 
				ans != 'y' && 
				ans != 'N' && 
				ans != 'n') {
			   cout<<"Do you want to keep the bond? (y/n): ";
			   ans = cin.get();
			   if( ans == '\n' )
			     cout << "Please enter \"y\" or \"n\"." << endl;
			   else
			     cin.ignore(80, '\n');
			 }	
			 if( ans == 'N' || ans == 'n' )
			   {
			     
			     for(k=j;k<=mult[i]-1;k++)
			       {
				 temp_link[i][k]=temp_link[i][k+1];
			       }
			     temp_link[i][mult[i]]=-1;
			     mult[i]--;
			     for(m=1;m<=mult[sf];m++)
			       {
				 if(temp_link[sf][m]==so)
				   {
				     for(q=m;q<=mult[sf]-1;q++)
				       {
					 temp_link[sf][q]=temp_link[sf][q+1];
				       }
				     temp_link[sf][mult[sf]]=-1;
				   }	
			       }
			     mult[sf]--;
			     no_bonds--;
			     for(k=1;k<=chk_bond_count;k++)
			       {
				 if(chk_bond[k][1]==so || chk_bond[k][1]==sf)
				   {
				     chk_bond[k][1]=-1;
				   }
			       }
			     for(k=1;k<=u_nrot_bond_count;k++)
			       {
			       if(u_nrot_bond[k][1]==so || u_nrot_bond[k][1]==sf)
				 {
				   u_nrot_bond[k][1]=0;
				 }
			       }
			     flag10=1;
			     j--;
			   }
		       }
		    }
		}
	    }
	}	
      if(flag10==1)
	{
	  cout<<endl<<"New valency = "<<mult[i]<<endl;
	}
      else if(flag10==-1)
	{
	  cout<<endl<<"Unchanged valency = "<<mult[i]<<endl;
	}
    }	
  if( flag7 == 1 ) 
    {
      cout<<endl<<endl<<"\t\tValency check is finished! "<<endl;	
      cout<<"---------------------------------------------------------------"<<endl;
      cout<<"Type \"s\"  to stop or \"c\" to continue: ";

      ans = cin.get();
      if( ans == 's' || ans == 'S' ) {
	cout << "proflex process killed by user.\n" << endl; 
	exit(-1); 
      }
      if( ans != '\n' )
	cin.ignore(80,'\n');

      cout << endl << endl;
    }
  } // ends the if(usage == 0) condition
//-------------copy the temp. connectivety table to new memory optimized table.
//------------------------------allocating memory for the new connectivity table
  link_noHB=(int**)new int*[no_atoms+1];
  for(i=0;i<=no_atoms;i++)
    {
      *(link_noHB+i)=(int*) new int[mult[i]+1];
    }
//---------------------------------copying the temp. table to the new one	
  for(i=1;i<=no_atoms;i++)
    {
      for(j=1;j<=mult[i];j++)
	{
	  link_noHB[i][j]=temp_link[i][j];
	}
    }
  
//--------------------------------------------------Unknown Non Rotating Bonds
  flag9 = 0;
  if( usage == 0 ) {
    if( u_nrot_bond_count > 0 ) {
      for( i = 1; i <= u_nrot_bond_count; i++)	{
	if( u_nrot_bond[i][1] != 0) {
	  // changed 11.02.01 AJR to only prompt for this double bond lock if a peptide bond.
	  if( (atype[u_nrot_bond[i][1]][2] == 'C' && atype[u_nrot_bond[i][2]][2] != 'N') ||
	      (atype[u_nrot_bond[i][1]][2] == 'N' && atype[u_nrot_bond[i][2]][2] != 'C') ) {
	    so = u_nrot_bond[i][1];
	    sf = u_nrot_bond[i][2];
	    
	    /************************************************************/
	    /* Print the warning message title only once.               */
	    if( flag9 == 0 ) {
	      cout<<endl<<"\tDihedral Angle Constraints:"<<endl;
	      cout<<endl<<"Should the following bond be considered as "
		  <<"Non-Rotating?"<<endl;	
	      flag9 = 1;
	    }
	    
	    cout<<endl<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	    write_output(1,so,sf,atype[so],p[so]->r1.rname,p[so]->r1.res_no,
			 p[so]->r1.code,p[so]->r1.chain,atype[sf],p[sf]->r1.rname,
			 p[sf]->r1.res_no,p[sf]->r1.code,p[sf]->r1.chain);
	    cout<<"\tType y to LOCK the dihedral angle or any key to keep unlocked: ";
	    
	    ans = cin.get();
	    if(ans == 'y' || ans == 'Y' ) {
	      nrot_bond_count++;
	      nrot_bond[nrot_bond_count][1]=u_nrot_bond[i][1];
	      nrot_bond[nrot_bond_count][2]=u_nrot_bond[i][2];
	    }
	    if( ans != '\n' )
	      cin.ignore( 80, '\n' );

	  }
	}
      }
    }	
  } // ends the if(usage == 0) condition
  
  //---count no. of isolated water atoms, isolated Hydrogens and all isolated atoms
  for(i=1;i<=no_atoms;i++)
    {
      if(mult[i]==0)
	{
	  all_iso++;
	  if( (strcmp(p[i]->r1.rname,"HOH")==0 || 
	       strcmp(p[i]->r1.rname,"DOD")==0) && p[i]->r1.aname[2]=='O')
	    {
	      O_iso++;
	    }
	}
    }
//------------------------------------Set DAH_type-------------------------------
  
  h_list_count = 0;
  set_DAH_type(usage, mult);
  flag8 = 0;
  #ifdef DEBUG
  cout << "do i get here " << endl;
  #endif
  for( i = 0; i < chaincount; i++) {
    tempo = chainptr[i];
    while( strcmp( tempo->r1.aname,"  N  ") != 0 ) { // catches end of ATOM records
      if( tempo == start ) {
	flag8 = 1;
	break;
      }
      tempo = tempo->prior;
    }
    
    //  modify C-terminus
    if( flag8 != 1 ) {
#ifdef DEBUG
      cout << "in flag8 != 1" << endl;
#endif
      tempo->r1.DAH_type = 'C';
      tempo = chainptr[i];
      if( i != 0 ) {
	tempo = tempo->next;
	oxtflag = 0;
	// Check whether the residue is known or unknown
	for( j = 0; j < 20; j++ ) {
	  if(strcmp(tempo->r1.rname,res[j].rname)==0) {
	    n = tempo->r1.mod_res_no;
	    while(tempo->r1.mod_res_no==n) {
	      if(strcmp(tempo->r1.aname,"  OXT")==0 || 
		 strcmp(tempo->r1.aname,"  OT ")==0 || 
		 strcmp(tempo->r1.aname,"  OT1")==0 || 
		 strcmp(tempo->r1.aname,"  OT2")==0) {
		tempo->r1.DAH_type='E';
		oxtflag=n;
#ifdef DEBUG
		cout << "found OXT? ?" << oxtflag << endl;
#endif
	      }
	      tempo=tempo->next;
	    }
	    tempo=tempo->prior;
	    while(strcmp(tempo->r1.aname,"  O  ")!=0)	tempo=tempo->prior;
	    if( oxtflag == 0 ) { 
	      if( usage == 0 ) {
		cout << endl << "\tThe final residue, # " << tempo->r1.res_no << " of chain ID "
		     << tempo->r1.chain <<endl << "\twas not identified as a terminal oxygen "
		     << "(OXT)."<<endl<<"\t" 
		     << "Enter (y)es if the main chain oxygen in this residue"<<endl 
		     << "\tshould be regarded as a charged acceptor at the C-terminus: "<< endl;
		
		ans = cin.get();
		
		if( ans == 'y'  ||
		    ans == 'Y' ) {
		  tempo->r1.DAH_type='E';			      
		}
		if( ans != '\n' )
		  cin.ignore( 80, '\n' );
		
	      }
	      else {
		tempo->r1.DAH_type='A';
	      }
	    }
	    else {
	      if(strcmp(tempo->r1.aname,"  O  ") == 0 )
		tempo->r1.DAH_type='E';
	    }
	  }
	}
      }
    }
  }
  flag8 = 0;
  tempo = start;
  while( strcmp( tempo->r1.field1,"ATOM  ") != 0 )
    {
      if(tempo==last)
	{
	  flag8=1;
	  break;
	}
      tempo=tempo->next;
    }
  if(flag8!=1)
	{
	  for(j=0;j<20;j++)
	    {
	      if(strcmp(tempo->r1.rname,res[j].rname)==0)
		{
		  n=tempo->r1.mod_res_no;
		  while(tempo->r1.mod_res_no==n)
		    {
		      if(strcmp(tempo->r1.aname,"  OXT")==0 || 
			 strcmp(tempo->r1.aname,"  OT ")==0 || 
			 strcmp(tempo->r1.aname,"  OT1")==0 || 
			 strcmp(tempo->r1.aname,"  OT2")==0)
			{
			  tempo->r1.DAH_type='E';
			}
		      if(strcmp(tempo->r1.aname,"  O  ")==0 )
			{
			  /*if(usage == 0){
      cout << endl << "\tEnter (n) if the oxygen in residue " << tempo->r1.rname 
           << " at the end" << endl << "\tof chain ID [" << tempo->r1.chain 
           << "] and residue # " << tempo->r1.res_no << " should NOT be"  
           << endl << "\tregarded as a charged acceptor at a C-terminus: " ;
       
       cin >> answer;
       if( strlen(answer) == 1 ) ans = answer[0];
       else ans = 'Y';
       if( ans == 'n' ) ans = 'N';
       if(ans!='N') 
	 tempo->r1.DAH_type='E';
			}
			
			else tempo->r1.DAH_type='E';*/
			  tempo->r1.DAH_type='E';   
			}      
		      if(tempo==last)
			{
			  break;
			}
		      tempo = tempo->next;
		    }
		  break;
		}
	    }
	}

  for( i = 1; i <= no_atoms; i++ ) {
    
    if( p[i]->r1.aname[2] == 'H' || p[i]->r1.aname[2] == 'D' ) {
      
      so = p[i]->r1.ri_sr_no;
      if( mult[i]==1 ) {
	sf = link_noHB[so][1];
	if( p[sf]->r1.DAH_type=='B' || p[sf]->r1.DAH_type=='C' || 
	    p[sf]->r1.DAH_type=='D' ) 
	  {
	    p[so]->r1.DAH_type='V';
	    h_list[h_list_count] = so;
	    h_list_count++;
	  }
	else {
	  p[so]->r1.DAH_type='N';
	}
      }
      else {
	p[so]->r1.DAH_type='N';
      }
    }
  }

  #ifdef DEBUG
  cout << "one" << endl;
  #endif
// search for hydrophobic tethers AJR 03.20.02 using new function
  no_tether = find_hydrophob();
		
  #ifdef DEBUG
  cout << "three" << endl;
  #endif
// Write all record info into _proflexdataset.
  make_chem_rec();
  #ifdef DEBUG
  cout << "four" << endl;
  #endif
  
// Write three new pseudoatoms for each Hydrophobic Tethers 
  write_HPatom(p);

#ifdef DEBUG
  cout << "five" << endl;
#endif

  l2.add_conect();

#ifdef DEBUG
  cout << "six" << endl;
#endif
// Append CF connectivity to _proflexdataset file

  yofil.open(outputfile,ios::app);
  yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
  yofil<<"REMARK:cf    so    sf     (atom-label pairs listing central-force "
       <<"bonds)"<<endl;
  for(i=1;i<=no_atoms;i++)
    {
      for(j=1;j<=mult[i];j++)
	{
	  if(link_noHB[i][j]>i)
	    {
	    yofil<<"REMARK:CF "<<setw(5)<<i<<" "<<setw(5)<<link_noHB[i][j]<<endl;
	    }
	}
    }
//  write CF bonds for Hydrophobic tethers.
  write_HPCFbond();
  yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
//-------------------------------Append TF connectivity to the _proflexdataset file
  flag9=0;
  for(i=1;i<=nrot_bond_count;i++)
    {
      if(flag9==0)
	{
	  yofil<<"REMARK:tf    so    sf     (atom-label pairs listing locked "
               <<"dihedral angles)"<<endl;
	  flag9=1;
	}
      yofil<<"REMARK:TF "<<setw(5)<<nrot_bond[i][1]<<" "<<setw(5)
	   <<nrot_bond[i][2]<<endl;
    }
  yofil<<"REMARK:L:"<<dash1<<dash2<<endl;

//----------------------------------Writing CHECK_BONDING into _proflexdataset file
  if(chk_bond_count>0)
    {
    yofil<<"REMARK:w:check_bonding:  Dist.  A1#   A1  Res1 Res1# C1 A2#   A2  "
         <<"Res2 Res2# C2"<<endl; 
      for(i=1;i<=chk_bond_count;i++)
	{
	  if(chk_bond[i][1]!=-1)
	    {
	      so=chk_bond[i][1];
	      sf=chk_bond[i][2];
	      yofil<<"REMARK:W:check_bonding:"
		   <<setiosflags(ios::showpoint|ios::fixed)
		   <<setprecision(5)<<setw(5)<<chk_bond_dist[i]<<setw(5)<<so
		   <<atype[so]<<"  "<<p[so]->r1.rname<<" "<<p[so]->r1.res_no
		   <<p[so]->r1.code[0]<<" "<<p[so]->r1.chain; 
	      yofil<<setw(5)<<sf<<atype[sf]<<"  "<<p[sf]->r1.rname<<" "
		   <<p[sf]->r1.res_no<<p[sf]->r1.code[0]<<" "<<p[sf]->r1.chain
		   <<endl;
	    }
	}
      yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
	}
//--------------------------------Write number of isolated water in _proflexdataset
  if(all_iso>0)
    {
    yofil<<"REMARK:w:isolated_atoms: Number of isolated atoms = "<<all_iso<<endl;
     if(O_iso>0)
       {
       yofil<<"REMARK:w:isolated_water: Number of isolated water atoms = "
	    <<O_iso<<endl;
       }
     if(skip_hydro_count>0)
       {
       yofil<<"REMARK:w:isolated_hydrogen: Number of isolated hydrogen atoms = "
	    <<skip_hydro_count<<endl;
       }
     yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
    }
  
//------------------------------Write number of disulfide bonds in _proflexdataset
  if(ns_bond>0)
    {
      yofil<<"REMARK:w:disulfide_bonds: Number of disulfide bonds = "<<ns_bond
	   <<endl;
      yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
    }
//--------------------------------Write SKIP_NONHYDROGEN atoms into _proflexdataset
  if(skip_nhydro_count>0)
    {
      yofil<<"REMARK:w:skip_nonhydrogen: Atom# Atom  Res  Res# Chain"<<endl;
      for(i=1;i<=skip_nhydro_count;i++)
	{
	  so=skip_nhydro[i];
	  yofil<<"REMARK:W:skip_nonhydrogen:"<<setw(5)<<so<<" "<<atype[so]<<"  "
	       <<p[so]->r1.rname<<"  "<<p[so]->r1.res_no<<p[so]->r1.code[0]
	       <<p[so]->r1.chain<<endl;
	}
      yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
    }		

  /************************************************************/
  /* Writing POOR_BOND into _proflexdataset file.               */
  /************************************************************/
  if( poor_bond_count > 0 ) {
    yofil<<"REMARK:w:poor_bond:     Dist.  A1#  A1   Res1 Res1# C1  A2#  A2  "
	 <<"Res2 Res2# C2"<<endl; 
    for( i = 1; i <= poor_bond_count; i++ ) {
      if( poor_bond[i][1] != -1 ) {
	so = poor_bond[i][1];
	sf = poor_bond[i][2];
	yofil<<"REMARK:W:poor_bond:"<<setiosflags(ios::showpoint|ios::fixed)
	     <<setprecision(5)<<setw(10)<<poor_bond_dist[i]<<setw(5)<<so
	     <<atype[so]<<"  "<<p[so]->r1.rname<<"  "<<p[so]->r1.res_no
	     <<p[so]->r1.code[0]<<" "<<p[so]->r1.chain; 
	yofil<<setw(5)<<sf<<atype[sf]<<"  "<<p[sf]->r1.rname<<"  "
	     <<p[sf]->r1.res_no<<p[sf]->r1.code[0]<<" "<<p[sf]->r1.chain
	     <<endl;
      }
    }
    yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
  }
  /************************************************************/
  
  /************************************************************/
  /* Writing SKIP_HYDROGEN into _proflexdataset.                */
  /************************************************************/
  if( skip_hydro_count > 0 ) {
    yofil<<"REMARK:w:skip_hydrogen: Atom# Atom  Res  Res# Chain"<<endl;
    for( i = 1; i <= skip_hydro_count; i++ ) {
      so = skip_hydro[i];
      yofil<<"REMARK:W:skip_hydrogen:"<<setw(5)<<so<<" "<<atype[so]<<"  "
	   <<p[so]->r1.rname<<"  "<<p[so]->r1.res_no<<p[so]->r1.code[0]
	   <<p[so]->r1.chain<<endl;
    }
    yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
  }
  /************************************************************/

  /************************************************************/
  /* Print messages to screen: Missing atoms                  */
  /************************************************************/
  if( ma_count > 0 ) {
    yofil<<"REMARK:w:missing_atoms:  Atom  Res  Res# Chain"<<endl;
    for( i = 1; i <= ma_count; i++) {
      so = missing[i];
      yofil<<"REMARK:W:missing_atoms:"<<"   "<<missing_atoms[i]<<"  "
	   <<p[so]->r1.rname<<" "<<p[so]->r1.res_no<<p[so]->r1.code<<"  "
	   <<p[so]->r1.chain<<endl;
    }
    yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
  }	
  /************************************************************/

  /************************************************************/
  /* Writing chain-change info. to _proflexdataset              */
  /************************************************************/
  if( chain_change_count > 0 ) {
    cout<<endl<<dash1<<dash3<<endl;
    cout<<"The following chains were relabelled:"<<endl;
    yofil<<"REMARK:w:chain_id: Old Cid  New Cid"<<endl;
    cout<<"chain_id: Old Cid  New Cid"<<endl;
    for(i=1;i<=chain_change_count;i++) {
      yofil<<"REMARK:W:chain_id:     "<<chain_change[i][0]<<"        "
	   <<chain_change[i][1]<<endl;
      cout<<"              "<<chain_change[i][0]<<"        "
	  <<chain_change[i][1]<<endl;
    }
    yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
  }
  /************************************************************/

  /************************************************************/
  /* In the case where alternative side-chain positions are   */
  /* reported in the input PDB file, only one conformation is */
  /* used. Here, that conformation is being written to the    */
  /* proflexdataset file.                                       */
  /************************************************************/
  if(altloccount>0) {
    if(conf_option==1) {
      yofil<<"REMARK:w:conformation_auto: Residue Residue no. Chain "
	   <<"Conformation"<<endl;
      for(i=0;i<altloccount;i++) {
	yofil<<"REMARK:W:conformation_auto:"<<"   "<<altlocptr[i]->r1.rname
	     <<"    "<<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code[0]
	     <<"        "<<altlocptr[i]->r1.chain<<"       "<<conf_chosen[i]
	     <<endl;
      }
    }
    else {
      yofil<<"REMARK:w:conformation_user: Residue Residue no. Chain "
	   <<"Conformation"<<endl;
      for(i=0;i<altloccount;i++) {
	yofil<<"REMARK:W:conformation_user:"<<"   "<<altlocptr[i]->r1.rname
	     <<"    "<<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code[0]
	     <<"        "<<altlocptr[i]->r1.chain<<"       "<<conf_chosen[i]
	     <<endl;
	
      }
    }
    yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
  }
  /************************************************************/
  
  /************************************************************/
  /* Writing total no. of atoms to _proflexdataset file.        */
  /************************************************************/
  yofil<<"REMARK:s:total_atoms: "<<no_atoms<<endl;		
  yofil<<"REMARK:L:"<<dash1<<dash2<<endl;
  yofil.close();
  /************************************************************/

  /************************************************************/
  /* Free up some memory from variables no longer needed.     */
  /************************************************************/
  for( i = 0; i <= no_atoms; i++) {
    delete [] *(chk_bond+i); 
    delete [] *(poor_bond+i); 
    delete [] *(atype+i);
    delete [] *(missing_atoms+i);
    delete [] *(nrot_bond+i);
    delete [] *(sxyz+i);
    delete [] *(temp_link+i);
  }
  
  for( i = 0; i < 1000; i++) {
    delete [] *(u_nrot_bond+i);
  }
  delete [] chk_bond;
  delete [] poor_bond;
  delete [] sxyz;
  delete [] atype;
  delete [] temp_link;
  delete [] nrot_bond;
  delete [] u_nrot_bond;
  delete [] missing_atoms;
  delete [] missing;
  delete [] skip_nhydro;
  delete [] skip_hydro;
  delete [] chk_bond_dist;
  delete [] poor_bond_dist;
  /************************************************************/
  
  xofil.close();
  zofil.close();
}

/**********************************************************************/
int list::check_dist(int so,int sf,float d, char *atom1, char *atom2, 
		     int type, char *rname1,char *res_no1,char *code1,
		     char *chain1,char *rname2, char *res_no2,char *code2,
		     char *chain2, int usage ) {

  char 
    ans;

  switch( atom1[2] ) {
  case 'C':
    {a1 = C;
    break;}
  case 'N':
    {a1=N;
    break;}
  case 'O':
    {a1 = O;
    break;}
  case 'S':
    {a1 = S;
    break;}
  case 'H':
    {a1 = H;
    break;}
  case 'D':
    {a1 = H;
    break;}
  case 'P':
    {a1 = P;
    break;}
  default:
    {a1 = Oth;
    break;}
  }
 
/*
 * 2006:07:31	SN
 *
 * To handle metal atoms, we include the following specific atom-type check 
 * List of metals handled: Co, Cu, Fe, K, Mn, Mg, Na, Ni, Zn, and calcium CA 
 * Note the position switch for calcium
 */
  if((atom1[2] == 'C' && (atom1[3] == 'O' || atom1[3] == 'U')) 
     || (atom1[2] == 'N' && (atom1[3] == 'I' || atom1[3] == 'A'))
     || (atom1[2] == 'M' && (atom1[3] == 'G' || atom1[3] == 'N'))
     || (atom1[2] == 'F' && atom1[3] == 'E') 
     || (atom1[2] == 'Z' && atom1[3] == 'N')
     || atom1[2] == 'K' 
     || (atom1[3] == 'C' && atom1[4] == 'A'))
     a1 = Oth; 
 
  switch( atom2[2] ) {
  case 'C':
    {a2 = C;
    break;}
  case 'N':
    {a2=N;
    break;}
  case 'O':
    {a2 = O;
    break;}
  case 'S':
    {a2 = S;
    break;}
  case 'H':
    {a2 = H;
    break;}
  case 'D':
    {a1 = H;
    break;}
  case 'P':
    {a2 = P;
    break;}
  default:
    {a2 = Oth;
    break;}
  }

/*
 * 2006:07:31	SN
 *
 * To handle metal atoms, we include the following specific atom-type check 
 * List of metals handled: Co, Cu, Fe, K, Mn, Mg, Na, Ni, Zn
 */
  if((atom2[2] == 'C' && (atom2[3] == 'O' || atom2[3] == 'U')) 
     || (atom2[2] == 'N' && (atom2[3] == 'I' || atom2[3] == 'A'))
     || (atom2[2] == 'M' && (atom2[3] == 'G' || atom2[3] == 'N'))
     || (atom2[2] == 'F' && atom2[3] == 'E') 
     || (atom2[2] == 'Z' && atom2[3] == 'N')
     || atom2[2] == 'K'
     || (atom2[3] == 'C' && atom2[4] == 'A'))
     a2 = Oth; 
 
  if( type == 1 ) {
    if(d<distance[a1][a2].d[0])	{ //------dist<d_min
      if(usage) { return(0);}//------------don't make connection
      else {
	cout<<endl<<"Caution: The distance of [d ="<<setiosflags(ios::showpoint | ios::fixed)
	    <<setprecision(2)<<setw(5)<<d<<" Ang.] seems low!"<<endl; 
	cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	write_output(1,so,sf,atom1,rname1,res_no1,code1,chain1,atom2,rname2,res_no2,code2,chain2);
	cout<<"Type  n  NOT to connect: ";
	cin >> answer;

	if( strlen(answer) == 1 ) ans = answer[0];
	else ans = 'Y';
	if( ans == 'n' ) ans = 'N';
	if(ans!='N') {
	  return (1);//-------------make connection
	}
	else {
	  return (0);//------------don't make connection
	}
      }
    }
    else {
      if(d>distance[a1][a2].d[3]) {	//----------dist>d_max
	if(d < maxdist) {
	  if(usage) {return(0);}//------------don't make connection
	  else{
	    cout<<endl<<"Caution: The distance of [d ="<<setiosflags(ios::showpoint|ios::fixed)
		<<setprecision(3)<<setw(5)<<d<<" Ang.] seems high!"<<endl; 
	    cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	    write_output(1,so,sf,atom1,rname1,res_no1,code1,chain1,atom2,rname2,res_no2,code2,chain2);
	    cout<<"Type  n  NOT to connect: ";
	    cin >> answer;
	    if( strlen(answer) == 1 ) ans = answer[0];
	    else ans = 'Y';
	    if( ans == 'n' ) ans = 'N';
	    if(ans!='N') {
	      return (1);//----------make connection
	    }
	    else {
	      return (0);//---------don't make conection
	    }
	  }
	}
	else {
	  return (0);//------------don't make connection
	}	
      }
      else {
	return (1);
      }
    }
  }
  else if( type == 2 ) {
    if(d<distance[a1][a2].d[0] || d>distance[a1][a2].d[3]) {//------(dist < d_min) or (dist > d_max)
      return (0);	//--------------don't make connection
    }
    else {
      return (1);	//--------------make connection
    }	
  }
  else if( type == 3 ) {
    if(d<distance[a1][a2].d[0] || d>distance[a1][a2].d[3]) { //------(dist < d_min) or (dist > d_max)
      return (0);	//--------------don't make connection
    }
    else {
      cout<<endl<<"Atoms belong to different chains!"<<endl;
      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
      write_output(1,so,sf,atom1,rname1,res_no1,code1,chain1,atom2,rname2,res_no2,code2,chain2);
      cout<<"Type  n  NOT to connect: ";
      cin >> answer;
      if( strlen(answer) == 1 ) ans = answer[0];
      else ans = 'Y';
      if(ans!='N' && ans!='n') {
	return (1);//----------make connection
      }
      else {
	return (0);//---------don't make conection
      }
    }	
  }
  else {	
    if(d<distance[a1][a2].d[0] || d>distance[a1][a2].d[3]) {//------(dist < d_min) or (dist > d_max)
      return (0);	//--------------don't make connection
    }
    else {
      if(d<distance[a1][a2].d[1]) { //------d_min < dist < d_low
	cout<<endl<<"Caution: The distance of [d ="<<setiosflags(ios::showpoint|ios::fixed)
	    <<setprecision(2)<<setw(5)<<d<<" Ang.] seems low!"<<endl; 
	/*cout<<"\t"<<setw(5)<<so<<atom1<<rname1<<res_no1<<endl; 
	  cout<<"\t"<<setw(5)<<sf<<atom2<<rname2<<res_no2<<endl;*/ 
	cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	write_output(1,so,sf,atom1,rname1,res_no1,code1,chain1,atom2,rname2,res_no2,code2,chain2);
	return (1);
      }
      else {
	if(d>distance[a1][a2].d[2]) { //----------d_big < dist < d_max
	  cout<<endl<<"Caution: The distance of [d ="<<setiosflags(ios::showpoint|ios::fixed)
	      <<setprecision(2)<<setw(5)<<d<<" Ang.] seems high!"<<endl; 
	  cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	  write_output(1,so,sf,atom1,rname1,res_no1,code1,chain1,atom2,rname2,res_no2,code2,chain2);
	  return (1); //----make connection
	}
	else {
	  return (1);	//----------make connection
	}
      }
    }
  }
}
/**********************************************************************/	
			

/**********************************************************************/
/* This function had the first variable renamed to my_mult since it is   */
/* that a local variable is desired -- anyhow that would have been the   */
/* way the function treated mult in the past since it shadowed the class */
/* variable with the same name -- vanvoor4 April 13, 2010                */
void list::make_connection(int *my_mult, int **temp_link, int so, int sf) {
  if( my_mult[so] >= maxr -1 ) {
    // 2007:03 SN --- "Multiplicity" is a misnomer. Change it to "Valency"
    cout << "Valency of atom # " << so << " is (" << my_mult[so] << ") greater than ";
    cout << maxr-1 << "; therefore EXIT!" << endl;
    exit(20);
  }
  if( my_mult[sf] >= maxr -1 ) {
    cout << "Valency of atom # " << sf << " is (" << my_mult[sf] << ")greater than ";
    cout << maxr-1 << "; therefore EXIT!" << endl;
    exit(20);
  }
  my_mult[so]++;
  my_mult[sf]++;
  temp_link[so][my_mult[so]] = sf;
  temp_link[sf][my_mult[sf]] = so;
}
/**********************************************************************/

/**********************************************************************/
/* Subroutine to format output to stdout. Based on commonly outputting*/
/* these variables at various points during FIRST.                    */
/**********************************************************************/
void list::write_output(int type, int so, int sf, char *atom1, char *rname1, 
			char *res_no1,char *code1, char *cid1, char *atom2, 
			char *rname2, char *res_no2,char *code2,char *cid2) {

  if( type == 1 ) { //---------write to screen(one below the other)
    cout<<"\t"<<setw(5)<<so<<"  "<<atom1<<"  "<<rname1<<"    "<<res_no1<<code1<<"  "<<cid1<<endl; 
    if( sf != 0 ) {
      cout<<"\t"<<setw(5)<<sf<<"  "<<atom2<<"  "<<rname2<<"    "<<res_no2<<code2<<"  "<<cid2<<endl; 
    }	
  }
	
  if( type == 2) { //---------write to screen(side by side)
    cout<<"\t"<<setw(5)<<so<<"  "<<atom1<<"  "<<rname1<<"    "<<res_no1<<code1<<"  "<<cid1; 
    if( sf != 0 ) {
      cout<<"\t"<<setw(5)<<sf<<"  "<<atom2<<"  "<<rname2<<"    "<<res_no2<<code2<<"  "<<cid2<<" : "; 
    }	
  }
}
/**********************************************************************/


