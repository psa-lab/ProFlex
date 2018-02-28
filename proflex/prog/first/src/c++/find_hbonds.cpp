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

/* THIS IS THE C++ CODE FOR FIND_HBONDS */
/*       Last revised:  20.03.06 by Sandeep Namilikonda  */
/*       Last revised:  10.01.01 by AJR  */
/*  uses  super-new and sb1 energy fcns */
/*  involves the -non option and is setup for porous linux g++  */
/*  3.19.01 added the webfirst DAHtype loop tags */

/**********************************************************************/
/* For each line in the polar_lookup file, read the hydrogen bonding  */
/* properties of the atoms (ie. B = "both donor and acceptor", N = "no*/
/* hydrogen bonding"). This routine uses operator overloading, which  */
/* I don't like, but hey it's a free world, and I'm not that familiar */
/* with C++ yet. Each time the ifstream "ifil" is accessed with an    */
/* overload operator (it's first accessed in main.cpp, then here, but */
/* it's a global var.) it set the next set of characters to the       */
/* variable. For instance, "num" below is declared to be a global var */
/* of type int in "class.h". ifil then sets the next white-space-de-  */
/* limited characters to num. This requires that the format of the    */
/* infile be strictly adhered to.                                     */
/**********************************************************************/
void DAHspecs::get_polar_record(char *temp_res) {

  int 
    i, j;

  strcpy(rname,temp_res);
  ifil >> num;

  for( i = 0; i < num; i++) {
    ifil >> aname[i];
 
    for(j = 0; j < 5; j++) {
      if( aname[i][j] == '_' ) 
	aname[i][j] = ' ';
    }
    ifil >> DAH_type[i];
  }

}
/**********************************************************************/

/**********************************************************************/
/* What does this do?                                                 */
/* This routine examines every atom in the input file to see if it can*/
/* be a hydrogen bond donor or acceptor. Default atom types are set   */
/* automatically. Unknown types are flagged (using the "warn1[] array */
/* The user is prompted to identify the hydrogen bonding status of    */
/* unknown atom types.                                                */
/**********************************************************************/
/* This function had the last variable renamed to mult_array since it is   */
/* that a local variable is desired -- anyhow that would have been the   */
/* way the function treated mult in the past since it shadowed the class */
/* variable with the same name -- vanvoor4 April 13, 2010                */
void list::set_DAH_type(int usage, int *mult_array) {

  int 
    i, j, 
    flag = 0, 
    flag1 = 0, 
    loop_flag = 1, 
    flag2 = 0, 
    so, 
    sf;
  
  char 
    ans = 'p';

  node 
    *tempo, 
    *temp;

  residue 
    *ptr[no_res],
    *t;

  /**********************************************************************/
  /* Declare and initialize the warn1[] array. Values of -1 indicate    */
  /* the indicated atom is of an unknown (to FIRST) type.               */
  /**********************************************************************/
  int warn1[no_atoms+1];
  for( i = 0; i < no_atoms; i++) {
    warn1[i] = 0;
  }

  /**********************************************************************/
  /* begin stepping through each atom to asign a DAH_type. last is      */
  /* really 1st atom of pdb. tempo is a linked list, of type "node",    */
  /* with all the atoms from the pdb file.                              */
  /**********************************************************************/
  tempo = last;  
  for( i = 0; i < no_res; i++) {
    ptr[i] = &res[i];
  }

  while( loop_flag ) {
    
    if( strncmp(tempo->r1.field1,"TER",3) == 0 ) {
      if(tempo==start) {
	break;
      }
      tempo=tempo->prior;
    }
    flag = 0;

    /**********************************************************************/
    /* Check to see if the residue name of the current atom is the same as*/
    /* that in ptr[]?                                                     */
    /**********************************************************************/
    if(strcmp(ptr[0]->rname,tempo->r1.rname) == 0 ) {
      flag = 1;
    }
    else {  // flag Known residue types
      for( i = 0; i < no_res; i++) {
	if(strcmp(ptr[i]->rname,tempo->r1.rname) == 0 ) {
	  flag = 1;
	  t = ptr[i];
	  ptr[i] = ptr[0];
	  ptr[0] = t;
	  break;
	}
      }
    }

    if( flag == 1 ) {
      /**********************************************************************/
      /* If flag == 1, then the current residue is a standard amino acid.   */
      /* The main-chain and side-chain hydrogen bonding atoms are set here. */
      /**********************************************************************/

      /* Set DAH atom types for SERINE */
      if(strcmp(tempo->r1.rname,"SER")==0 || 
	 strcmp(tempo->r1.rname,"THR")==0){
	if(strncmp(tempo->r1.aname+2,"OG",2) == 0 ) {
	  tempo->r1.DAH_type='B';
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ") == 0 ) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}
      }

      /* Set DAH atom types for CYSTIENE */
      else if(strcmp(tempo->r1.rname,"CYS")==0) {
	if(strncmp(tempo->r1.aname+2,"SG",2)==0) {
	/* Brandon pointed out that 'free thiol' should be treated as sp3 donor and acceptor
	 * Probable fix for Future, after discussions and test cases -Sameer
	 * POssible fix: tempo->r1.DAH_type='N';
	 */
	  tempo->r1.DAH_type='D';
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}
      }

      /* METHIONINE */
      else if(strcmp(tempo->r1.rname,"MET")==0) {
	if(strncmp(tempo->r1.aname+2,"SD",2)==0) {
	  tempo->r1.DAH_type='A';
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}
      }	
      
      /* ASPARTATE */
      else if(strcmp(tempo->r1.rname,"ASP")==0) {	
	if(strncmp(tempo->r1.aname+2,"OD",2)==0) {
	  tempo->r1.DAH_type='E';
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}
      }	

      /* ASPARAGINE */
      else if(strcmp(tempo->r1.rname,"ASN")==0) {	
	if(strncmp(tempo->r1.aname+2,"OD",2)==0) {
	  tempo->r1.DAH_type='A';
	}
	else {
	  if(strncmp(tempo->r1.aname+2,"ND",2)==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else {
	    if(strcmp(tempo->r1.aname,"  N  ")==0) {
	      tempo->r1.DAH_type='D';
	    }
	    else if(tempo->r1.aname[2]=='O') {
	      tempo->r1.DAH_type='A';
	    }
	    else {
	      tempo->r1.DAH_type='N';
	    }
	  }
	}
      }	

      /* GLUTAMATE */
      else if(strcmp(tempo->r1.rname,"GLU")==0) {	
	if(strncmp(tempo->r1.aname+2,"OE",2)==0) {
	  tempo->r1.DAH_type='E';
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}
      }	

      /* GLUTAMINE */
      else if(strcmp(tempo->r1.rname,"GLN")==0) {	
	if(strncmp(tempo->r1.aname+2,"NE",2)==0) {
	  tempo->r1.DAH_type='D';
	}
	else {
	  if(strncmp(tempo->r1.aname+2,"OE",2)==0) {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    if(strcmp(tempo->r1.aname,"  N  ")==0) {
	      tempo->r1.DAH_type='D';
	    }
	    else if(tempo->r1.aname[2]=='O') {
	      tempo->r1.DAH_type='A';
	    }
	    else {
	      tempo->r1.DAH_type='N';
	    }
	  }
	}
      }	

      /* LYSINE */
      else if(strcmp(tempo->r1.rname,"LYS")==0) {	
	if(strncmp(tempo->r1.aname+2,"NZ",2)==0) {
	  tempo->r1.DAH_type='C';
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}	
      }	

      /* ARGININE */
      else if(strcmp(	tempo->r1.rname,"ARG")==0) {	
	if(strncmp(tempo->r1.aname+2,"NE",2)==0 
	   ||strncmp(tempo->r1.aname+2,"NH",2)==0) {	
	  tempo->r1.DAH_type='C';
	}	
	else {	
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else {
	    if(tempo->r1.aname[2]=='O') {
	      tempo->r1.DAH_type='A';
	    }
	    else {
	      tempo->r1.DAH_type='N';
	    }
	  }
	}
      }

      /* HISTIDINE */
      else if(strcmp(tempo->r1.rname,"HIS")==0) {
	if(strncmp(tempo->r1.aname+2,"ND",2)==0 || 
	   strncmp(tempo->r1.aname+2,"NE",2)==0) {
	  so=tempo->r1.ri_sr_no;
	  for(i=1;i<=mult_array[so];i++) {
	    sf=link_noHB[so][i];
	    for(temp=last;temp->r1.ri_sr_no!=sf;temp=temp->prior);
	    if(temp->r1.aname[2]=='H' ||temp->r1.aname[2]=='D') {
	      tempo->r1.DAH_type='C';
	      flag1=1;
	      break;
	    }
	  }
	  if(flag1!=1) {
	    tempo->r1.DAH_type='A';
	  }
	  flag1=0;
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}
      }

      /* TYROSINE */
      else if(strcmp(tempo->r1.rname,"TYR")==0) {
	if(strncmp(tempo->r1.aname+2,"OH",2)==0) {
	  tempo->r1.DAH_type='B';
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}
      }

      /* TRYPTOPHAN */
      else if(strcmp(tempo->r1.rname,"TRP")==0) {
	if(strncmp(tempo->r1.aname+2,"NE",2)==0) {
	  tempo->r1.DAH_type='D';
	}
	else {
	  if(strcmp(tempo->r1.aname,"  N  ")==0) {
	    tempo->r1.DAH_type='D';
	  }
	  else if(tempo->r1.aname[2]=='O') {
	    tempo->r1.DAH_type='A';
	  }
	  else {
	    tempo->r1.DAH_type='N';
	  }
	}
      }

      /* PROLINE */
      else if(strcmp(tempo->r1.rname,"PRO")==0) {
	if(tempo->r1.aname[2]=='O') {
	  tempo->r1.DAH_type='A';
	}
	else {
	  tempo->r1.DAH_type='N';
	}
      }

      /* for all other known groups */
      else {
	if(strcmp(tempo->r1.aname,"  N  ")==0) {
	  tempo->r1.DAH_type='D';
	}
	else if(tempo->r1.aname[2]=='O') {
	  tempo->r1.DAH_type='A';
	}
	else {
	  tempo->r1.DAH_type='N';
	}
      }
    
      /* Point to the next atom */
      if( tempo != start ) {
	tempo = tempo->prior;
      }
      else {
	loop_flag = 0;
	break; 
      }
    }
    
    /**********************************************************************/
    /* Evaluate the DAH type for unknown groups.                          */
    else {
      if(tempo->r1.aname[2] == 'C') {
	tempo->r1.DAH_type = 'N';
      }
      else { // if not a carbon 
	for( i = 0; i < polar_count; i++) {
	  if( strcmp( tempo->r1.rname, polar[i].rname ) == 0 ) {
	    for( j = 0;  j < polar[i].num; j++) { 
	      if( strcmp( polar[i].aname[j], tempo->r1.aname+1 ) == 0 ) {
		tempo->r1.DAH_type = polar[i].DAH_type[j][0];
		goto ID_successfull;
	      }
	    }
	  }
	}
	if( tempo->r1.aname[2] != 'H' && 
	    tempo->r1.aname[2] != 'D' ) {
	  so = tempo->r1.ri_sr_no;
	  warn1[so] = -1;

	  if( usage == 0 ) {
	    if( flag2 == 0 ) {
	      cout<<endl<<"---------------------------------------------------------"; 
	      cout<<endl<<"\tIdentification of H-bond Donors, Acceptors etc."<<endl<<endl;
	      cout<<"\t(d) H-bond Donor \n\t(a) H-bond Acceptor \n\t(b) Both H-bond Donor and H-bond Acceptor \n\t(c) ";
	      cout<<"Charged Donor \n\t(e) Charged Acceptor \n\t(n) None "<<endl;
	      cout<<"---------------------------------------------------------"<<endl;
	      cout<<"\tEnter (a,b,c,d,e or n) for the following atoms."<<endl<<endl;
	      cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
	      flag2 = 1;
	    }	
	    
	    ans = '@';
	    while( ans!='A' && ans!='B' && ans!='C' && ans!='D' &&
		   ans!='E' && ans!='N' && ans!='a' && ans!='b' && 
		   ans!='c' && ans!='d' && ans!='e' && ans!='n') {

	      write_output(2, so, 0, tempo->r1.aname, tempo->r1.rname, tempo->r1.res_no,
			   tempo->r1.code, tempo->r1.chain,tempo->r1.aname, tempo->r1.aname,
			   tempo->r1.aname, tempo->r1.aname, tempo->r1.aname);
	      cout<<"  Enter (a,b,c,d,e or n): ";
	      ans = cin.get();
	      
	      if( ans == '\n' )
		cout << "\tInvalid Option." << endl;
	      else
		cin.ignore( 80, '\n' );
	    }
	    tempo->r1.DAH_type=toupper(ans);

	  }
	  else {
	    /************************************************************/
	    /* Execute this loop if the user chose the -non command line*/
	    /* option. Here, unknown DAH types will be set to hard-coded*/
	    /* defaults.                                                */
	    /************************************************************/
	    // changed default DAH_type assignment of unknown ligands %%% AJR 05.07.02

	    int  hmult = 0;
	    if( tempo->r1.aname[2] == 'N' ) { 
	      for( i = 1; i <= mult_array[so]; i++) {
		sf = link_noHB[so][i];
		for( temp = last; temp->r1.ri_sr_no != sf; temp = temp->prior );
		if( temp->r1.aname[2] == 'H' || 
		    temp->r1.aname[2] == 'D') {
		  hmult++;
		}
	      }
	      if(hmult >= 3){
		tempo->r1.DAH_type='C';
	      }
	      else if(hmult == 0){
		tempo->r1.DAH_type='A';
	      }
	      else {
		tempo->r1.DAH_type='D';
	      }
	    }
	    else if( tempo->r1.aname[2] == 'O' || 
		     tempo->r1.aname[2] == 'S' ) {
	      for( i = 1; i <= mult_array[so]; i++ ) {
		sf = link_noHB[so][i];
		for( temp = last; temp->r1.ri_sr_no != sf; temp = temp->prior );
		if( temp->r1.aname[2] == 'H' ||
		    temp->r1.aname[2] == 'D' ) {
		  hmult++;
		}
	      }
	      if( hmult >= 1 ){
		tempo->r1.DAH_type='B';
	      }
	      else {
		tempo->r1.DAH_type='A'; 
	      }
	    }
	    else {
	      tempo ->r1.DAH_type='N';
	    }
	  }
	  // %%%%%%%%%%%%%%%%%%%%%%%%% AJR 05.07.02
	}
      }
    
    ID_successfull: 
      if( tempo != start ) {	
	tempo = tempo->prior;
	if( strcmp(tempo->r1.field1,"TER   ") == 0 ) {
	  if( tempo == start ) {
	    loop_flag = 0;
	    break;
	  }
	  tempo = tempo->prior;
	}
      }
      else {
	loop_flag=0;
	break;
      }

    }
  }
  /* END LOOP checking DAH type for all atoms                           */
  /**********************************************************************/
}					
/**********************************************************************/

/**********************************************************************/
/* Based on a given set of hard-coded geometric criteria, identify all*/
/* possible hydrogen bonds.                                           */
/**********************************************************************/
/* #define DEBUG */
void list::find_Hbonds(int usage) {

   int i,j,s,so,sf,iix,iiy,iiz,jx,jy,jz,kx,ky,kz,flag=0,tmp1,tmp2,tmp3,tmp4;
   int tsy,tsz,tmp0; // SN 02:2006
   int ibad,ibado,hb=0,flag2=0,nd,na,icase,sy,sx,sz,ijk;
   int **temp_hbond,*temp_hb_type;
   float *temp_hb_energy;
   float xda,xha,adha,eflag,hd2,ha2,ad2,hddist,hadist,addist;
   float angle,temp,costheta,angular,ay2,aydist,hy2,cosphi,factor;
   float uid2,uid10,tmp,radial,adsb1;
   float bx,by,bz,ax,ay,az,dx,dy,dz,cx,cy,cz,amag,dmag,temp0,ttemp0,tcosphi,t2;
   float theta,thetapi,factx=6.0;
   char ans;
   ofstream yofil,pacc_ofil;

#ifdef CHKPT
   ofstream hbnd_filt_log_ofil;
#endif

   mult_net = new int[no_atoms+1];
   temp_hb_energy = new float[no_atoms];
   temp_hb_type = new int[no_atoms];
   temp_hbond = (int **) new int *[no_atoms];
   for(i=0;i<no_atoms;i++)
      {
      *(temp_hbond+i) = (int *) new int[4];
      }
   for( i = 0; i < h_list_count; i++) {
     s = h_list[i];
     so=link_noHB[s][1];
     iix=(int)((p[so]->r1.coord[1]-xmin)/grdlen);
     iiy=(int)((p[so]->r1.coord[2]-ymin)/grdlen);
     iiz=(int)((p[so]->r1.coord[3]-zmin)/grdlen);
      
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
	       flag=0;
	       if( p[sf]->r1.DAH_type=='A')
	       	  {
		  if(p[so]->r1.aname[2]=='S' || p[sf]->r1.aname[2]=='S')
		     {                         	
 //--------------------------------------------------- H-bond with Sulfur
       		     xda=s_admaxd;
		     xha=s_hamaxd;
		     adha=s_dha_ang;
		     eflag=-1;
		     }
		  else	
//---------------------------------------------- H-bond without Sulfur
		     {
		     xda=admaxd;
		     xha=hamaxd;
		     adha=dha_ang;
		     eflag=1.0e-20;
		     }
		  }
	       else if(p[sf]->r1.DAH_type=='B' && so!=sf)
	       	  {
		  if(p[so]->r1.aname[2]=='S' || p[sf]->r1.aname[2]=='S')
		     {
//--------------------------------------------------- H-bond with Sulfur
     	       	     xda=s_admaxd;
		     xha=s_hamaxd;
		     adha=s_dha_ang;
		     eflag=-1;
		     }
		  else	//---------------------------- H-bond without Sulfur
		     {
		     xda=admaxd;
		     xha=hamaxd;
		     adha=dha_ang;
		     eflag=1.0e-20;
		     }
		  }
	       else if(p[sf]->r1.DAH_type=='E')
	       	  {
		  if(p[so]->r1.DAH_type=='C')	//-----------------Salt Bridge
		     {
		     xda=sb_admaxd;
		     xha=sb_hamaxd;
		     adha=sb_dha_ang;
		     eflag=-6;
		     }
		  else	//------------------------------- H-bond without sulfur
		     {
		     xda=admaxd;
		     xha=hamaxd;
		     adha=dha_ang;
		     eflag=1.0e-20;
		     }
		  }
	       else
		  {
		  sf=chain[sf];
		  flag=1;
		  }
	       if(flag==0)
		  {
//---------------------------------------------------Calculate H-donor distance
		  hd2=0;
		  for(j=1;j<=3;j++)
		     {
		     hd2+=pow((p[so]->r1.coord[j]-p[s]->r1.coord[j]),2);
		     }
		  hddist=sqrt(hd2);
//------------------------------------------------Calculate H-acceptor distance
		  ha2=0;
		  for(j=1;j<=3;j++)
		     {
		     ha2+=pow((p[sf]->r1.coord[j]-p[s]->r1.coord[j]),2);
		     }
		  hadist=sqrt(ha2);
//--------------------------------------------Calculate Acceptor-Donor distance
		  ad2=0;
		  for(j=1;j<=3;j++)
		     {
		     ad2+=pow((p[sf]->r1.coord[j]-p[so]->r1.coord[j]),2);
		     }
		  addist=0;
//------------------------------------------------------------Check Angle
       		  temp=0;
		  angle=0;
		  temp=(hd2+ha2-ad2)/(2*hddist*hadist);
		  angle=acos(temp);
		  if(angle < adha) 
		     {
		     sf=chain[sf];
		     }
		  else
		     {
//----------------------------------------------------Check H-Acceptor distance
	       	     if(hadist < xha)	//---------------Record H-bond
		       	{
		       	nhb++;
		       	temp_hbond[nhb][1]=so;
		       	temp_hbond[nhb][2]=s;
		       	temp_hbond[nhb][3]=sf;
		       	temp_hb_energy[nhb]=eflag;
		       	}
		     else
		       	{
		       	addist=sqrt(ad2);
		        if(addist < xda)	//------------Record H-bond
		       	   {
			   nhb++;
			   temp_hbond[nhb][1]=so;
			   temp_hbond[nhb][2]=s;
			   temp_hbond[nhb][3]=sf;
			   temp_hb_energy[nhb]=eflag;
			   }
			}
		     sf=chain[sf];
		     }
		  }
	       }
	    }
	 }
      }
   }
//--------------------------- remove H-bonds between residues i and (i-1,i,i+1)

#ifdef CHKPT
   cout<<"\nH-bond #: "<<nhb;
printf("\n\t\t\tCHECK POINT 1"); 
for(int count=1;count<=nhb;count++)
  {
    printf("\n%d --- %d\t%d\t%d\t%f",count,temp_hbond[count][1],temp_hbond[count][2],temp_hbond[count][3]);
  }
#endif

   int flag2b=-1;
   j = 0;
   for(i=1;i<=nhb;i++)
      {
      flag2=0;
      so = temp_hbond[i][1];
      sf = temp_hbond[i][3];
//      itemp = abs( p[so]->r1.mod_res_no - p[sf]->r1.mod_res_no );
/* -----------------------------------------------------------------AJR 5/12/00
	REMOVED explicit i to i, i+/-1 HBond taboo.  The modified Mayo energy
	function we now use insures that these non-physical Hbonds will have
	an energy of 0.0 kcal/mol and can be screened against by setting Ecut.

      if( itemp < 2 )          // prescreening
         {
         if( strcmp(p[so]->r1.chain,p[sf]->r1.chain)==0 ) 
            {
            if( strcmp(p[so]->r1.aname,"  N  ")==0 && ( 
                strcmp(p[sf]->r1.aname,"  O  ")==0 || 
                strcmp(p[sf]->r1.aname,"  OXT")==0 || 
                strcmp(p[sf]->r1.aname,"  OT1")==0 || 
                strcmp(p[sf]->r1.aname,"  OT2")==0    ) ) 
                {
                   tmp1 = atoi( p[so]->r1.res_no );              
                   tmp2 = atoi( p[sf]->r1.res_no );              
                   if( abs( tmp2 - tmp1 ) < 2 )
                      {
                      flag2 = 1;
		      for(k=0;k<no_res;k++)
			 {
                         if(strcmp(res[k].rname,p[so]->r1.rname)==0)
                            {
                            flag2++;
                            }
                         if(strcmp(res[k].rname,p[sf]->r1.rname)==0)
                            {
                            flag2++;
                            }
                         }
                      }
                }
            }
         }*/
// -----------------------------------------------------------------AJR 5/12/00
      if( flag2 == 0 ) 
         {
         j++;
         temp_hbond[j][1] = temp_hbond[i][1]; 
         temp_hbond[j][2] = temp_hbond[i][2];
         temp_hbond[j][3] = temp_hbond[i][3]; 
         temp_hb_energy[j]= temp_hb_energy[i];
         }
      else if( flag2 < 3 ) 
         { 
         if( flag2b < 0 )
            {                 
 cout<<endl<<"---------------------------------------------------------------";
 cout<<endl<<"\tVerification required on questionable hydrogen bonds."<< endl;
 cout<<   "\tEach hydrogen bond satisfies all geometrical criteria." << endl;
 cout<<"----------------------------------------------------------------";
 cout << endl;
            flag2b = 1;
            }
         if(usage == 0) {
 cout << endl<<"\t   Donor: "<<p[so]->r1.aname<<" "<<p[so]->r1.rname<<" "
      << p[so]->r1.res_no<<p[so]->r1.code<<endl;
 cout << "\tAcceptor: "<<p[sf]->r1.aname<<" "<<p[sf]->r1.rname<<" "
      << p[sf]->r1.res_no<<p[sf]->r1.code<<endl;
 cout << "\tType  n  NOT to consider as a hydrogen bond: ";
 //cin >> answer;
            cin.seekg(0,ios::end);
            cin.clear();
            cin.get(answer,4);
            cin.seekg(0,ios::end);
            cin.clear();
            if( strlen(answer) == 1 ) ans = answer[0];
            else ans = '@';
         }
         else ans = 'n';
            if( ans == 'n' ) ans = 'N';
            if( ans != 'N' )
               {
               j++;
               temp_hbond[j][1] = temp_hbond[i][1];
               temp_hbond[j][2] = temp_hbond[i][2];
               temp_hbond[j][3] = temp_hbond[i][3];
               temp_hb_energy[j]= temp_hb_energy[i];
	       }
            }
         }
         nhb = j;
//-----------------------------------------Calculate new (maximum) multiplicity
	for(i=0;i<=no_atoms;i++)
	   {
	   mult_net[i]=mult[i];
	   }
	for(i=1;i<=nhb;i++)
	   {
	   s=temp_hbond[i][2];
	   sf=temp_hbond[i][3];
	   mult_net[s]++;
	   mult_net[sf]++;
	   }
//-------------------------------------------- Allocate memory for link_net[][]
	
	link_net = (int **) new int *[no_atoms+1];
	for(i=0;i<=no_atoms;i++)
	{
		*(link_net+i)=(int *) new int[mult_net[i]+1];
	}

//-------------------------------------------------Copy link_noHB to link_net
			
 	for(i=1;i<=no_atoms;i++)
	   {
	   for(j=1;j<=mult[i];j++)
	      {
	      link_net[i][j]=link_noHB[i][j];
	      }
	   mult_net[i]=mult[i];
	   }
//------------------------------------------------------- Add H-acceptor bonds
	flag=0;
	ibad=nhb+1;
	ibado=ibad;
	hb++;
	while(hb<ibad)
	{
	    so=temp_hbond[hb][1];
	    s=temp_hbond[hb][2];
	    sf=temp_hbond[hb][3];
	    for(j=1;j<=mult_net[s];j++)
	       {
	       if(link_net[s][j]==sf)
	          {
		  ibad--;
		  temp_hbond[hb][1]=temp_hbond[ibad][1];
		  temp_hbond[hb][2]=temp_hbond[ibad][2];
		  temp_hbond[hb][3]=temp_hbond[ibad][3];
		  temp_hbond[ibad][1]=so;
	       	  temp_hbond[ibad][2]=s;
	       	  temp_hbond[ibad][3]=sf;
		  flag=1;
		  }
	       }
	    if(flag!=1)
	       {
	       mult_net[s]++;
	       mult_net[sf]++;
	       link_net[s][mult_net[s]]=sf;
	       link_net[sf][mult_net[sf]]=s;
	       hb++;
	       }
	    flag=0;
	}

//-------------------------Record no. of valid H-bonds and report the bad ones
        flag2=0;
	if(ibad<ibado)
	{
	   for(i=ibad;i<=ibado-1;i++)
	      {
	      if(flag2==0)
		 {

        cout<<endl<<"Warning: A covalent or Hydrogen CF-bond "<<endl;
       	cout<<"overlaps with a hydrogen-acceptor CF-bond"<<endl;
       	cout<<"\t Atom#  Atom  Res    Res#    Chain"<<endl; 
       	       	 flag2=1;
		 }	
	      so=temp_hbond[i][2];  // Hydrogen-Acceptor overlap (not donor)
	      sf=temp_hbond[i][3];
	      write_output(1,so,sf,p[so]->r1.aname,p[so]->r1.rname,
		    p[so]->r1.res_no,p[so]->r1.code,p[so]->r1.chain,
		    p[sf]->r1.aname,p[sf]->r1.rname,p[sf]->r1.res_no,
		    p[sf]->r1.code,p[sf]->r1.chain);
	      }	
	   cout<<endl<<endl;
	   cout<<"\tError found: Please check the original PDB file"<<endl;
	   cout<<"\t             and the corresponding proflexdataset." << endl;
           cout<<"\t             Error must be corrected to continue."<<endl;
	   exit(31);
	}
//------------------------STATUS CHECK --- 2006:03:01 SN
// At this point, the energies are initialized to either 0.0 or -6.0 for SB,
// and -1.0 for some of them (WHY -1.0 and for what?)


//------------------------Calculate the energy for each H-bond or Salt bridge
	for(i=1;i<=nhb;i++)
	   {
	   so=temp_hbond[i][1];
	   s=temp_hbond[i][2];
	   sf=temp_hbond[i][3];
//--------------------------------------------Calculate Acceptor-Donor distance
       	   ad2=0;	
	   for(j=1;j<=3;j++)
	      {
	      ad2+=pow((p[sf]->r1.coord[j]-p[so]->r1.coord[j]),2);
	      }
//------------------------------------------Calculate theta:angle between D-H-A
	   hd2=0;
	   for(j=1;j<=3;j++)
	      {
	      hd2+=pow((p[so]->r1.coord[j]-p[s]->r1.coord[j]),2);
	      }
	   hddist=sqrt(hd2);

	   ha2=0;
	   for(j=1;j<=3;j++)
	      {
	      ha2+=pow((p[sf]->r1.coord[j]-p[s]->r1.coord[j]),2);
	      }
	   hadist=sqrt(ha2);
	       
	   costheta=(hd2+ha2-ad2)/(2*hddist*hadist);
/* ----------- inserted and altered lines that change the theta dependence 
   from cos^2 to cos^2*exp(-(pi-theta)^6)		     AJR 6.8.00	       */
	   angular=costheta*costheta;
           theta = acos(costheta);
	   thetapi = 180*raddeg - theta;
           factor = pow(thetapi,factx);
           angular *= exp(-factor);
//	   angular=costheta*costheta;
/* ----------- inserted and altered lines that change the theta dependence
		from cos^2 to 1-(pi-theta)^2					*/
/*           theta = acos(costheta);
	   thetapi = 180*raddeg - theta;
	   if( thetapi > 1.0) 
	      { 
              angular = 0;
	      }
	   else	
              {
	      angular = 1 - thetapi*thetapi;	
              }
*/
//--------------------------------------------------take each case separately
	   if(temp_hb_energy[i]<-2)
	      {			//------------------Salt bridge E-function
//----- modified 6.9.00 by AJR to give sb1 --> a simple SB E-fcn
		//	      temp_hb_energy[i]=-9;

	      adsb1 = sqrt(ad2)+3.750e-1;
	      uid2=10.240e0/pow(adsb1,2);
	      uid10=pow(uid2,5);
	      temp_hb_energy[i]=10.0e0*uid10*(5.0e0*uid2-6.0e0);
	      		    // temp_hb_energy[i]=radial*angular;
	      temp_hb_type[i]=4;

//------------ STATUS CHECK 1
//printf("\n\t\t CHECKPOINT 1");
//printf("\nEnergy < -2 KCal/Mol : h-bond #: %d",i);

                        //----- Added on 2006:03:01 by SN to fix Preacptr
                        if(p[sf]->r1.aname[4]=='1')
                           sy = sf - 1; //-- This is *VERY* SPECIFIC to PDB file
                                        //   numbering scheme ex: CD OE1 OE2
                                        //   ex: CG OD1 OD2
                        else if(p[sf]->r1.aname[4]=='2')
                           sy = sf - 2; //-- ex: HIS  bonds: CG-ND1 CD2-NE2
                                        //       CB CG ND1 CD2 CE1 NE2
                        else
                           sy = sf - 1; //-- ex: Carbonyl Carbon C O
                                        //--     CYS: CB SG     SER: CB OG
                                        //--     MET: CG SD

                        temp_hbond[i][0]=sy;


	      }
	   else
	      {	 //-------H-bond E-function without Sulfur Hybridization table 
       				//-----------Donor hybridization table	      
       	      if(mult[so]<2)
	       	 {
		 nd=1;
		 }
	      else if(mult[so]<3)
	       	 {
		 nd=1;
		 if(p[so]->r1.aname[2]=='N')
		    {
		    nd=0;
		    }
		 if(p[so]->r1.aname[2]=='O' &&strcmp(p[so]->r1.rname,"TYR")==0)
		    {
		    nd=0;
		    }
		 }
	      else if(mult[so]<4)
	       	 {
		 nd=1;
		 if(p[so]->r1.aname[2]=='N')
		    {
		    nd=0;
		    }
		 }
	      else
	       	 {
		 nd=1;
		 }
	    //----------------------------check acceptor
	      if(mult[sf]<1)
	       	 {
		   na = 1;
		   /* This is probably wrong if(nd=0) */
		   if( nd == 0 ) {
		     angular *= angular;
		   }
		 }
	      else if(mult[sf]<2)
	       	 {		//------------------Determine hybridization
		 na=1;
		 if(p[sf]->r1.aname[2]=='S' || p[sf]->r1.aname[2]=='O')
		    {
		    na=0;
		    }
		    icase=na+(2*nd);
		    if(icase==1)
		       {
		       angular*=angular;

                        //----- Added on 2006:03:01 by SN to fix Preacptr
                        if(p[sf]->r1.aname[4]=='1')
                           sy = sf - 1; //-- This is *VERY* SPECIFIC to PDB file
                                        //   numbering scheme ex: CD OE1 OE2
                                        //   ex: CG OD1 OD2
                        else if(p[sf]->r1.aname[4]=='2')
                           sy = sf - 2; //-- ex: HIS  bonds: CG-ND1 CD2-NE2
					//	 CB CG ND1 CD2 CE1 NE2
                        else
                           sy = sf - 1; //-- ex: Carbonyl Carbon C O
					//--	 CYS: CB SG	SER: CB OG
					//--	 MET: CG SD

                        temp_hbond[i][0]=sy;

		       }
		    else
		       {
	       //-------------------------Calculate phi: angle between H-A-Y
	       	       sy=link_noHB[sf][1];
                // --- 2006:02:13 SN    - Fix for DELTA filter
		// ---------------------- Store the preacceptor info.
                       temp_hbond[i][0] = sy;

		       ay2=0;
		       for(j=1;j<=3;j++)
		       	  {
			  ay2+=pow((p[sf]->r1.coord[j]-p[sy]->r1.coord[j]),2);
			  }
		       aydist=sqrt(ay2);
		       hy2=0;
		       for(j=1;j<=3;j++)
		       	  {
			  hy2+=pow((p[s]->r1.coord[j]-p[sy]->r1.coord[j]),2);
			  }
//		         hydist=sqrt(hy2);
		       cosphi=(ha2+ay2-hy2)/(2*hadist*aydist);
		       if(icase==2)
		       	  {
			  angular*=pow(cosphi,2);
			  }
		       else if(icase==3)
		       	  {
			  temp=acos(cosphi)-109.5e0*raddeg;
			  angular*=pow(cos(temp),2);
			  }
		       else   //---------------------SP2 donor and SP2 acceptor
			  {
			  factor=cosphi*cosphi;
			  if(mult[so]>1 && mult[sy]>1)
			     {	//---------------------typical expected case
			     sz=link_noHB[sy][1];
			     if(sz==sf)
			       	{
			       	sz=link_noHB[sy][2];
			        }
			     bx=p[sf]->r1.coord[1]-p[sy]->r1.coord[1];
			     by=p[sf]->r1.coord[2]-p[sy]->r1.coord[2];
			     bz=p[sf]->r1.coord[3]-p[sy]->r1.coord[3];
			     cx=p[sz]->r1.coord[1]-p[sy]->r1.coord[1];
			     cy=p[sz]->r1.coord[2]-p[sy]->r1.coord[2];
			     cz=p[sz]->r1.coord[3]-p[sy]->r1.coord[3];
			     ax=by*cz-bz*cy;
			     ay=bz*cx-bx*cz;
	              	     az=bx*cy-by*cx;
			     amag=sqrt(ax*ax + ay*ay + az*az);
			     ax=ax/amag;
			     ay=ay/amag;
			     az=az/amag;
		//----------------------------Calculate plane normal for X-D-H
#ifdef DEBUG
        printf("\n\n\nphi = %6.3f \t phi term = %f",acos(cosphi)/raddeg,factor);
        printf("\ndonor: %d %s %s \t acptr: %d %s %s",
                so,p[so]->r1.aname,p[so]->r1.rname,
                sf,p[sf]->r1.aname,p[sf]->r1.rname);
#endif
	       		     for(j=1;j<=mult[so];j++)
			       	{
			       	sx=link_noHB[so][j];
			       	if(sx!=s && p[sx]->r1.aname[2]!='H')
			     	   {
				   bx=p[s]->r1.coord[1]-p[so]->r1.coord[1];
				   by=p[s]->r1.coord[2]-p[so]->r1.coord[2];
				   bz=p[s]->r1.coord[3]-p[so]->r1.coord[3];
				   cx=p[sx]->r1.coord[1]-p[so]->r1.coord[1];
				   cy=p[sx]->r1.coord[2]-p[so]->r1.coord[2];
				   cz=p[sx]->r1.coord[3]-p[so]->r1.coord[3];
				   dx=by*cz-bz*cy;
				   dy=bz*cx-bx*cz;
				   dz=bx*cy-by*cx;
				   dmag=sqrt(dx*dx + dy*dy + dz*dz);
				   dx=dx/dmag;
				   dy=dy/dmag;
				   dz=dz/dmag;
				   temp=pow((ax*dx + ay*dy + az*dz),2);

#ifdef DEBUG
        printf("\npre-donor: %d %s %s",sx,p[sx]->r1.aname,p[sx]->r1.rname);
        printf("\ngamma = %6.3f \t gamma term = %f",
		acos(ax*dx + ay*dy + az*dz)/raddeg,temp);
#endif
				   if(temp>factor)
				      {
				      factor=temp;
#ifdef DEBUG
        printf("\tGAMMA CHOSEN!");
#endif
				      }
				   }
			    	}
			     }
			  angular*=factor;
			  }
		       }
		    }
		 else if(mult[sf]<3)
		    {
		       //------------------------------Determine hybridization
		    na=1;
		    if(p[sf]->r1.aname[2]=='N')
		       {
		       na=0;
		       }
		    icase=na+2*nd;
		    if(icase==1)
		       {
		       angular*=angular;

//------------ STATUS CHECK 1.6
//printf("\n\t\t CHECKPOINT 1.6");
//printf("\nAcptr: 'O',  h-bond #: %d",i);

                        //----- Added on 2006:03:01 by SN to fix Preacptr
                        if(p[sf]->r1.aname[4]=='1')
                           sy = sf - 1; //-- This is *VERY* SPECIFIC to PDB file
                                        //   numbering scheme ex: CD OE1 OE2
                                        //   ex: CG OD1 OD2
                        else if(p[sf]->r1.aname[4]=='2')
                           sy = sf - 2; //-- ex: HIS  bonds: CG-ND1 CD2-NE2
                                        //       CB CG ND1 CD2 CE1 NE2
                        else
                           sy = sf - 1; //-- ex: Carbonyl Carbon C O
                                        //--     CYS: CB SG     SER: CB OG
                                        //--     MET: CG SD

                        temp_hbond[i][0]=sy;
		       }
		    else
		       {
		       	//-----------------------------Base atom1 
		       sy=link_noHB[sf][1];
		       ay2=0;
		       for(j=1;j<=3;j++)
		       	  {
			  ay2+=pow((p[sf]->r1.coord[j]-p[sy]->r1.coord[j]),2);
			  }
		       aydist=sqrt(ay2);
		       hy2=0;
		       for(j=1;j<=3;j++)
		       	  {
			  hy2+=pow((p[s]->r1.coord[j]-p[sy]->r1.coord[j]),2);
			  }
		       cosphi=(ha2+ay2-hy2)/(2*hadist*aydist);
		              //-----------------------------Base atom2     
       		       sz=link_noHB[sf][2];	
      		       ay2=0;
		       for(j=1;j<=3;j++)
		       	  {
			  ay2+=pow((p[sf]->r1.coord[j]-p[sz]->r1.coord[j]),2);
			  }
		       aydist=sqrt(ay2);
		       hy2=0;
		       for(j=1;j<=3;j++)
		       	  {
			  hy2+=pow((p[s]->r1.coord[j]-p[sz]->r1.coord[j]),2);
			  }
		       temp0=(ha2+ay2-hy2)/(2*hadist*aydist);
		       tcosphi=cosphi;
		       ttemp0=temp0;
// -- 2006:02:13 SN	-- Use similar mechanism as one used for picking 
//			   'DELTA' for picking the preacceptor atom
			tsy=sy;
			tsz=sz;

	  	     	if(p[sy]->r1.aname[2]=='H')
		       	  {
			  tcosphi=temp0;
			  tsy=sz; // SN 02:2006
			  }
			if(p[sz]->r1.aname[2]=='H')
			  {
			  ttemp0=cosphi;
			  tsz=sy; // SN 02:2006
			  }
			  cosphi=tcosphi;
			  temp0=ttemp0;
			  sy = tsy;
			  sz = tsz;

			  if(icase==2)
			     {
			     if(fabs(temp0)<fabs(cosphi))
			       	{
			       	cosphi=temp0;
				sy = sz; // SN 02:2006
			       	}
			     angular=angular*(pow(cosphi,2));
                                temp_hbond[i][0] = sy; // SN 02:2006
			     }
			  else if(icase==3)
			     {
			     temp=acos(cosphi)-109.5e0*raddeg;
			     t2=acos(temp0)-109.5e0*raddeg;
			     if(fabs(temp)<fabs(t2))
			      	{
			      	temp=t2;
				sy = sz; // SN 02:2006
			       	}
			     angular=angular*(pow(cos(temp),2));
				temp_hbond[i][0] = sy; // SN 02:2006
			     }
			  else
			     {
//------------ STATUS CHECK 1.7
//printf("\n\t\t CHECKPOINT 1.7");
//printf("\nAcptr: 'O',  h-bond #: %d",i);

 //--------------------------------------------------SP2 Donor and SP2 acceptor
       			     if(fabs(temp0)<fabs(cosphi))
			       	{
			       	cosphi=temp0;
			       	}
			     factor=pow(cosphi,2);
			     if(mult[so]>1)
			       	{
 //---------------------------------------------Calculate plane angle for Z-A-Y
	       			bx=p[sz]->r1.coord[1]-p[sf]->r1.coord[1];
	       			by=p[sz]->r1.coord[2]-p[sf]->r1.coord[2];
	       			bz=p[sz]->r1.coord[3]-p[sf]->r1.coord[3];
	       			cx=p[sy]->r1.coord[1]-p[sf]->r1.coord[1];
	       			cy=p[sy]->r1.coord[2]-p[sf]->r1.coord[2];
	       			cz=p[sy]->r1.coord[3]-p[sf]->r1.coord[3];
	       			ax=by*cz-bz*cy;
	       			ay=bz*cx-bx*cz;
	       			az=bx*cy-by*cx;
	       			amag=sqrt(ax*ax + ay*ay + az*az);
	       			ax=ax/amag;
	       			ay=ay/amag;
	       			az=az/amag;
//---------------------------------------------Calculate plane angle for X-D-H
		       		for(j=1;j<=mult[so];j++)
		       	       	   {
				   sx=link_noHB[so][j];
				   if(sx!=s && p[sx]->r1.aname[2]!='H')
				      {
				      bx=p[s]->r1.coord[1]-p[so]->r1.coord[1];
				      by=p[s]->r1.coord[2]-p[so]->r1.coord[2];
				      bz=p[s]->r1.coord[3]-p[so]->r1.coord[3];
				      cx=p[sx]->r1.coord[1]-p[so]->r1.coord[1];
				      cy=p[sx]->r1.coord[2]-p[so]->r1.coord[2];
				      cz=p[sx]->r1.coord[3]-p[so]->r1.coord[3];
				      dx=by*cz-bz*cy;
				      dy=bz*cx-bx*cz;
			              dz=bx*cy-by*cx;
				      dmag=sqrt(dx*dx + dy*dy + dz*dz);
				      dx=dx/dmag;
				      dy=dy/dmag;
				      dz=dz/dmag;
				      temp=pow((ax*dx + ay*dy + az*dz),2);
#ifdef DEBUG
        printf("\nphi = %6.3f   gamma = %6.3f",
                acos(cosphi)/raddeg,acos(ax*dx + ay*dy + az*dz)/raddeg);
        printf("\ndonor: %d %s %s \t acptr: %d %s %s\n",
                so,p[so]->r1.aname,p[so]->r1.rname,
                sf,p[sf]->r1.aname,p[sf]->r1.rname);
#endif

				      if(temp>factor)
				       	 {
					 factor=temp;
					 }
				      }
				   }
				}	
			     angular*=factor;

                             temp_hbond[i][0]=sy; //-- Added by SN Preacptr FIX


			     } //--- End of if(mult[so]>1)

			  } //--- End of else if(icase != 3 or 2)
//-------- STATUS CHECK 1.7 end
		        }
		     else
		        {
	       //-------------------------------------Determine Hybridization
	       		na=1;
		       	icase=na+2*nd;
		       	if(icase==1)
		       	   {
                         angular*=angular;
//  AJR 10.01.01  identified the below as bad code, replace with above line.
//                       angular=angular*2;

//------------ STATUS CHECK 1.8
//printf("\n\t\t CHECKPOINT 1.8");
//printf("\nh-bond #: %d",i);

                        //----- Added on 2006:03:01 by SN to fix Preacptr
                        if(p[sf]->r1.aname[4]=='1')
                           sy = sf - 1; //-- This is *VERY* SPECIFIC to PDB file
                                        //   numbering scheme ex: CD OE1 OE2
                                        //   ex: CG OD1 OD2
                        else if(p[sf]->r1.aname[4]=='2')
                           sy = sf - 2; //-- ex: HIS  bonds: CG-ND1 CD2-NE2
                                        //       CB CG ND1 CD2 CE1 NE2
                        else
                           sy = sf - 1; //-- ex: Carbonyl Carbon C O
                                        //--     CYS: CB SG     SER: CB OG
                                        //--     MET: CG SD

                        temp_hbond[i][0]=sy;


			   }
			else
		       	   {
	     //-------------------------------------Determine smallest cosphi
			   temp=1;
			   tsy = 0; // SN 02:2006
			   for(ijk=1;ijk<=mult[sf];ijk++)
			       {
			       sy=link_noHB[sf][ijk];
				if(tsy == 0) // SN 02:2006
				{
				  tsy = sy;
				}
			       if(p[sy]->r1.aname[2] !='H')
			       	  {
				  ay2=0;
				  for(j=1;j<=3;j++)
				     {
			ay2+=pow((p[sf]->r1.coord[j]-p[sy]->r1.coord[j]),2);
				     }
				  aydist=sqrt(ay2);
				  hy2=0;
				  for(j=1;j<=3;j++)
				      {
			hy2+=pow((p[s]->r1.coord[j]-p[sy]->r1.coord[j]),2);
		       		      }
			          cosphi=(ha2+ay2-hy2)/(2*hadist*aydist);
		       	          if(fabs(temp)<fabs(cosphi))
				     {
				     cosphi=temp;
				     tsy=sy; // SN 02:2006
				     }
				  temp=cosphi;
				  }
			       }
			    cosphi=temp;
			    sy = tsy; // SN 02:2006
			    temp_hbond[i][0] = sy; // SN 02:2006
			    temp=acos(cosphi)-109.5e0*raddeg;
			    angular=angular*pow(cos(temp),2);
			    }
			}
	 //----------------------------------------Evaluate energy function
	       	
		     uid2=7.840e0/ad2;
		     uid10=pow(uid2,5);
		     radial=8.0e0*uid10*(5.0e0*uid2-6.0e0);
		     temp_hb_energy[i]=radial*angular;
	      //---------------------------------------Define type of HB
		     temp_hb_type[i]=na+2*nd;
		}
	}
//---------------------------------Sort H-bonds from highest to lowest energy 
//                                                 Use a Shell sort algorithm

#ifdef CHKPT
//---------------------------------STATUS CHECK 2
printf("\n\t\t\tCHECK POINT 2"); 
for(int count=1;count<=nhb;count++)
  {
    printf("\n%d --- %d\t%d\t%d\t%f",count,temp_hbond[count][1],temp_hbond[count][2],temp_hbond[count][3],temp_hb_energy[count]);
  }
#endif

     temp_hb_energy[0]=1e25;
     nhb++;
     for(int gap = nhb/2; gap > 0; gap = gap == 2 ? 1 : (int)(gap/2.2))
	{
       	for(i=gap;i<nhb;i++)
       	   {
	   tmp=temp_hb_energy[i];
	   tmp0=temp_hbond[i][0]; // SN 02:2006
	   tmp1=temp_hbond[i][1];
	   tmp2=temp_hbond[i][2];
	   tmp3=temp_hbond[i][3];
	   tmp4=temp_hb_type[i];
	   j=i;
	       	
	   for(;j>=gap && tmp > temp_hb_energy[j-gap]; j-=gap)
	      {
	       	temp_hb_energy[j]=temp_hb_energy[j-gap];
		temp_hbond[j][0]=temp_hbond[j-gap][0]; // SN 02:2006
	       	temp_hbond[j][1]=temp_hbond[j-gap][1];
	       	temp_hbond[j][2]=temp_hbond[j-gap][2];
	       	temp_hbond[j][3]=temp_hbond[j-gap][3];
	       	temp_hb_type[j]=temp_hb_type[j-gap];
	      }
	   temp_hb_energy[j]=tmp;
	   temp_hbond[j][0]=tmp0; // SN 02:2006
	   temp_hbond[j][1]=tmp1;
	   temp_hbond[j][2]=tmp2;
	   temp_hbond[j][3]=tmp3;
	   temp_hb_type[j]=tmp4;
	   }
	}
     nhb--;

#ifdef CHKPT
//---------------------------------STATUS CHECK 3
printf("\n\t\t\tCHECK POINT 3");
for(int count=1;count<=nhb;count++)
  {
    printf("\n%d --- %d\t%d\t%d\t%d\t%f",count,temp_hbond[count][1],\
	temp_hbond[count][2],temp_hbond[count][3],temp_hbond[count][0],\
	temp_hb_energy[count]);   }
#endif


/**********************************************************/
// --- SN 02:2006/06:2007
// ---  Implement initial geometric filters here ( PARTLY REDUNDANT 
//						   see lines 687-722 )
// ---  Update nhb
// ---  and writeout the filtered output into the proflexdataset
// --- SN 05:2007
// --- 	Output all the H-bonds filtered warning messages
//	into a file named: filtered_Hbond_log	--- FOR DEBUGGING PURPOSE ONLY!
/**********************************************************/

#ifdef CHKPT
  hbnd_filt_log_ofil.open("filtered_Hbonds.log",ios::out);
#endif

  xda = admaxd;
  xha = hamaxd;

  float min_theta = dha_ang,cosdelta,alt_cosdelta;
  float min_delta = hap_ang,min_costheta,min_cosdelta,tmp_delta,alt_delta;
  float ap2,hp2,ad2max,ha2max;
  int   *tmp_hb_pick,sp,khb,tmp_sp;
  char  preacptr_flag;

  /************************************************************/
  /* Allocate memory for hb_pick[]                            */
  /************************************************************/
  tmp_hb_pick = new int[nhb+1];

#ifdef CHKPT
  cout<<"\n\n****************************************************";
    cout<<"\n*          H-bond Length and Angle Check!          *";
    cout<<"\n* (H-bonds are indexed as per the _proflexdataset) *";
  
  hbnd_filt_log_ofil<<"\n\n****************************************************";
    hbnd_filt_log_ofil<<"\n*          H-bond Length and Angle Check!          *";
    hbnd_filt_log_ofil<<"\n* (H-bonds are indexed as per the _proflexdataset) *";
#endif

  khb=nhb;
  for(i=1;i<=nhb;i++)
    {
      /************************************************************/
      /* Perform stereochemical filtering to obtain the initial   */
      /* list of H-bonds                                          */
      /************************************************************/
      so=temp_hbond[i][1];
      s=temp_hbond[i][2];
      sf=temp_hbond[i][3];
      sp=temp_hbond[i][0];

      /************************************************************/
      /* 2007:02:28	SN					  */
      /* 							  */
      /* Preacptr atom can be missing if 'H' atoms are not present*/
      /* on water molecules. So, perform preacptr-based filtering */
      /* IF AND ONLY IF preactpr atom number is within the valid  */
      /* limit of 1 to no_atoms					  */
      /************************************************************/
      if (sp > 0 && sp < MAXATOMS) {
         preacptr_flag = 'T'; }
      else {
         preacptr_flag = 'F';
	 temp_hbond[i][0] = 0; // initialize to zero
         sp = 0;
      }

      tmp_hb_pick[i] = 1;	

      //----------------------------------------------------Determine geometrical case

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
          /* 2007:02 SN*/
          if(preacptr_flag == 'T') {
            hp2+=pow((p[s]->r1.coord[j]-p[sp]->r1.coord[j]),2);
            ap2+=pow((p[sf]->r1.coord[j]-p[sp]->r1.coord[j]),2);
          }
        }
      //----------------------------------------------------------Calculate cos(angle)
      costheta = (hd2+ha2-ad2)/sqrt(4*hd2*ha2);
      min_costheta = cos(min_theta);

      /* 2007:02 SN*/
      if(preacptr_flag == 'T') {
        cosdelta = (ha2+ap2-hp2)/sqrt(4*ha2*ap2);
        min_cosdelta = cos(min_delta);

	if(strcmp(p[sf]->r1.rname,"HOH") == 0) {
          tmp_delta = acos((double)cosdelta)/raddeg;
          if(tmp_delta < 90) {
	     if(p[sp-1]->r1.aname[2] == 'H')
		  tmp_sp = sp-1;
             else
          	  tmp_sp = sp+1;
#ifdef CHKPT
	     cout<<"\n\n#"<<i<<" preacptr: "<<sp<<", and alt_preacptr: "<<tmp_sp<<endl;
	     cout<<"Delta: "<<tmp_delta<<endl;
//<<" and cosdelta: "<<cosdelta<<endl;
//
#endif
             hp2=0;
	     ap2=0;
             for(j=0;j<3;j++) {
	    	hp2+=pow((p[s]->r1.coord[j]-p[tmp_sp]->r1.coord[j]),2);
            	ap2+=pow((p[sf]->r1.coord[j]-p[tmp_sp]->r1.coord[j]),2);
	     }
	     alt_cosdelta = (ha2+ap2-hp2)/sqrt(4*ha2*ap2);
             alt_delta = acos((double)alt_cosdelta)/raddeg;
#ifdef CHKPT
	     cout<<"alt_delta: "<<alt_delta<<endl;
//" and alt_cosdelta: "<<alt_cosdelta<<endl;
#endif

// If the other H atom in water gives a more obtuse angle 
// then choose that atom as preacptr
	     if(alt_delta > tmp_delta) {
                 sp = tmp_sp;
		 temp_hbond[i][0] = sp;
		 cosdelta = alt_cosdelta;
	     }
          }		
	} // End of 'HOH' handling
      } // End of preacptr flag set

#ifdef CHKPT
// ------ SN 2006:03 ------ DEBUG INFO -------
printf("\n%d   H-donor dist: %f  D-A dist: %f  theta: %f  delta: %f",\
	i,sqrt(hd2),sqrt(ad2),acos((double)costheta)/raddeg,\
	acos((double)cosdelta)/raddeg);
#endif

      //----------------------------------------------------------Geometric filters
                  if(temp_hb_type[i]<4)
                    {
                      ha2max = hamaxd*hamaxd;
                      ad2max = admaxd*admaxd;
                    }
                  else if(temp_hb_type[i]==4)
                    {
                      ha2max = sb_hamaxd*sb_hamaxd;
                      ad2max = sb_admaxd*sb_admaxd;
                    }
                  else if(temp_hb_type[i]==5)
                    {
                      ha2max = sb_hamaxd*sb_hamaxd;
                      ad2max = sb_admaxd*sb_admaxd;
                    }
                  else
                    {
                      ha2max = s_hamaxd*s_hamaxd;
                      ad2max = s_admaxd*s_admaxd;
                    }

/*
 * 2007:07:07	SN	--- Salt bridges do not have to be subjected 
 *			    to angular check and can be accepted if H-A or
 *			    D-A distance is acceptable
 */
      if((temp_hb_type[i] != 4 && (ad2 > ad2max || ha2 > ha2max)) || \
	 (temp_hb_type[i] == 4 && ha2 > ha2max && ad2 > ad2max)) 
	{
          khb--;
          tmp_hb_pick[i]=0;

#ifdef CHKPT

          cout<<"\n\n*   WARNING: Preferred distance range violated.  ";
          cout<<"\n*            Filtered out the following bond:   "<<i;
	  cout<<"\n*            Donor: "<<so<<" H: "<<s<<" Acptr: "<<sf;
          cout<<"\n"<<i<<"  H-acceptor dist: "<<sqrt(ha2)<<", Maximum allowed: "<<sqrt(ha2max);
	  cout<<"\n"<<i<<"  donr-acptr dist: "<<sqrt(ad2)<<", Maximum allowed: "<<sqrt(ad2max);	

          hbnd_filt_log_ofil<<"\n\n*   WARNING: Preferred distance range violated.  ";
          hbnd_filt_log_ofil<<"\n*            Filtered out the following bond:   "<<i;
	  hbnd_filt_log_ofil<<"\n*            Donor: "<<so<<" H: "<<s<<" Acptr: "<<sf;
          hbnd_filt_log_ofil<<"\n"<<i<<"  H-acceptor dist: "<<sqrt(ha2)<<", Maximum allowed: "<<sqrt(ha2max);
	  hbnd_filt_log_ofil<<"\n"<<i<<"  donr-acptr dist: "<<sqrt(ad2)<<", Maximum allowed: "<<sqrt(ad2max);

#endif
        }
      else if(temp_hb_type[i] != 4 && \
	      (costheta > min_costheta || \
         	/* 2007:02 SN --- Ensure that preacptr 
				  exists before checking delta*/
                (preacptr_flag == 'T' && cosdelta > min_cosdelta)))
	{	
          khb--;
          tmp_hb_pick[i]=0;

#ifdef CHKPT

          cout<<"\n\n*   WARNING: Preferred angular range violated.   ";
          cout<<"\n*            Filtered out the following bond:   "<<i;
	  cout<<"\n*            Donor: "<<so<<" H: "<<s<<" Acptr: "<<sf;
          cout<<"\n"<<i<<"  donr-H-acptr angle: "<<acos((double)costheta)/raddeg;
          cout<<", Minimum theta required: 110";
          cout<<"\n"<<i<<"  H-acptr-preacptr angle: "<<acos((double)cosdelta)/raddeg;
          cout<<", Minimum delta required: 90";

          hbnd_filt_log_ofil<<"\n\n*   WARNING: Preferred angular range violated.   ";
          hbnd_filt_log_ofil<<"\n*            Filtered out the following bond:   "<<i;
	  hbnd_filt_log_ofil<<"\n*            Donor: "<<so<<" H: "<<s<<" Acptr: "<<sf;
          hbnd_filt_log_ofil<<"\n"<<i<<"  donr-H-acptr angle: "<<acos((double)costheta)/raddeg;
          hbnd_filt_log_ofil<<", Minimum theta required: 110";
          hbnd_filt_log_ofil<<"\n"<<i<<"  H-acptr-preacptr angle: "<<acos((double)cosdelta)/raddeg;
          hbnd_filt_log_ofil<<", Minimum delta required: 90";

#endif

	}
    } // --- End for(i)
  /* ----- Geometric filter implementation ends ----- */

#ifdef CHKPT
  cout<<"\n\n*   End of default H-bond length and angle check!    *";
  cout<<"\n******************************************************\n";

  hbnd_filt_log_ofil<<"\n\n*   End of default H-bond length and angle check!    *";
  hbnd_filt_log_ofil<<"\n******************************************************\n";

  /*
   * 2007:07:07 SN	--- Make a note in proflexdataset that there were
   *			    WARNING messages and that they can be read from
   *			    filtered_Hbonds.log file
   */
/*
  if(khb < nhb)
   {
	yofil.open(outputfile,ios::app);
	yofil<<"REMARK:w:Hbnd_filtered: See filtered_Hbonds.log";
	yofil.close();
   }
*/

#endif
	

 /************************************************/
 /* Create memory for permanent arrays		 */
 /* --- 2006:03:07 SN	--- 			 */
 /* khb is the # of h-bonds picked at this point */
 /* Hence, allocate arrays of size khb+1 	 */
 /* and not nhb+1				 */
 /************************************************/

     hb_energy = new float[khb+1];
     hb_type = new int[khb+1];
     hb_id = new int[khb+1];
     hbond = (int **) new int *[khb+1];
     for(i=0;i<=khb;i++)
	{
		*(hbond+i) = (int *) new int[4];
	}
//--------------------------------------Copy temp arrays into permanent arrays
     j = 0;
     for(i=1;i<=nhb;i++)
	{
	  if(tmp_hb_pick[i] == 1)
	    {
		j++;
	        hb_energy[j]=temp_hb_energy[i];
		hb_type[j]=temp_hb_type[i];
		hbond[j][0]=temp_hbond[i][0]; // SN 02:2006
		hbond[j][1]=temp_hbond[i][1];
		hbond[j][2]=temp_hbond[i][2];
		hbond[j][3]=temp_hbond[i][3];
	    }
	}

     if(j != khb)
	{
	  cout<<"\nERROR: Initial stereochemical filter failure";
	  exit(1);
	}

// --- Free the temp array used		- SN
     delete [] tmp_hb_pick;

//---------------------- Write *FILTERED* data in _proflexdataset
//---- 2006:03:15 SN --- Write out preacptr info to preacptr_info file
     pacc_ofil.open("preacptr_info",ios::out);

     yofil.open(outputfile,ios::app); //--- 2006:03 SN
     if(khb>0)
     {
	yofil<<"REMARK:hb   ID   E(Kcal/Mol)  Donor  Hydrogen Acceptor  "
	     <<"Description"<<endl;
	for(i=1;i<=khb;i++)
       	   {

           pacc_ofil<<setw(5)<<hbond[i][0]<<endl; //--- 2006:03 SN

	   hb_id[i]=i;
   yofil<<"REMARK:HB"<<setw(5)<<i<<setiosflags(ios::showpoint|ios::fixed)
	<<setprecision(5)<<setw(13)<<hb_energy[i]<<"   "<<setw(5)
	<<hbond[i][1]<<"   "<<setw(5)<<hbond[i][2]<<"   "
       	<<setw(5)<<hbond[i][3];
       	   switch(hb_type[i])
	       	{
	       	case 0: yofil<<"    HB Dsp2 Asp2"<<endl;
		        break;
	       	case 1: yofil<<"    HB Dsp2 Asp3"<<endl;
		        break;
		case 2: yofil<<"    HB Dsp3 Asp2"<<endl;
	       		break;
	       	case 3: yofil<<"    HB Dsp3 Asp3"<<endl;
	       		break;
	       	case 4: yofil<<"    SB no energy"<<endl;
	       		break;
	       	default: cout<<"\tError: Unknown hydrogen bond type at i="<<i
                             <<endl<<endl;
	        exit(32);
	                break;
		}
	    }
/************************** AJR 03.22.02 *********************/
 write_HPHBbond();

 yofil<<"REMARK:L:-----------------------------------------------------------"
      <<"------------"<<endl;
	}
     yofil<<"END"<<endl;
     yofil.close();

     pacc_ofil.close();

#ifdef CHKPT
     hbnd_filt_log_ofil.close();
#endif
//----------------------------------------------Delete temp. arrays
     delete [] temp_hb_energy;
     delete [] temp_hb_type;
     for(i=0;i<=nhb;i++)
	{
	  delete [] *(temp_hbond+i);
	}
     delete [] temp_hbond;

// --- SN 02:2006
// ---  Update "nhb"; j should be equal to "khb" at this point

     nhb = khb;

}
