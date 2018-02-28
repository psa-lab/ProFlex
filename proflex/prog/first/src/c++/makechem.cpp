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
/* makechem.cpp                                                       */
/*                                                                    */
/* These routines are used to create the proflexdataset file using the  */
/* original pdb file of the protein being analyzed. The lines are read*/
/* in, sorted into a linked list of class node, and output to the     */
/* *_proflexdataset file. If multiple conformations of a given amino-   */
/* acid exist, they are dealt with in this file.                      */
/**********************************************************************/

/**********************************************************************/
/* getrecord reads the individual lines of the pdb file, and parses   */
/* the data from the ATOM and HETATM lines into a data structure of   */
/* type record. records are a subclass of the class node. The data is */
/* read into a linked list of type node.                              */
/**********************************************************************/
#include <map>
void record::getrecord(void) {

  int 
    i, l,
    is_atom = 0, is_hetatm = 0, is_ter = 0;
    
  char line[90];
  char sr_mod_num[6];
  int modelno,linelengthflag=0,shortlineflag=0;


  //if( !ifil.getline(line, 90) ) // if the next line is the EOF, return
  //  return;
  
  ifil.getline( line, 90 );
  line[89] = '\0';
  

  is_atom   = !strncmp( line, "ATOM",   4 );
  is_hetatm = !strncmp( line, "HETATM", 6 );
  is_ter    = !strncmp( line, "TER",    3 );

  /* This code seems to keep reading the lines and skip CONECT/MASTER/END
   * records until either a ATOM/HETATM/TER record appears. Meanwhile if it  
   * gets the model number, it sets the model no and returns  - Sameer Mar 4 2004
   */

  while( !is_atom && !is_hetatm && !is_ter && !ifil.eof() ){
    // print to file all lines that aren't CONECT, MASTER, and END. 
    // _______ SN 2007:01:30	--- Also skip ANISOU records ______
    if( strncmp(line, "CONECT", 6) && 
	strncmp(line, "ANISOU", 6) &&
	strncmp(line, "MASTER", 6) && 
	strncmp(line, "END",    3)){
// __trying to cope with MODEL & LINK input lines_________AJR 05.08.02
      if( !strncmp(line, "MODEL",5)) { 
	    strncpy( sr_mod_num,  line+11, 5 );
	    modelno = atoi(sr_mod_num);
	    modelflag=modelno;
	    if(modelflag >1) return;
      }
      ofil << line << "\n";
    } 
    
    // The following code checks to see if the getline operation on ifil failed.
    // This will occur if the input buffer exceeds the buffer size of the variable
    // line. BMH 5.24.02. 
    if( ifil.fail() ){
      cout << "Error reading from input file in record::getrecord.\n" << endl;
      exit(14);
    }
   

    if( ifil.getline( line, 90 ) ) {
      is_atom   = !strncmp( line, "ATOM",   4 );
      is_hetatm = !strncmp( line, "HETATM", 6 );
      is_ter    = !strncmp( line, "TER",    3 );
    }
  }
  
 
  if( is_atom || is_hetatm || is_ter ){
// AJR 07.31.02 these lines and below deals with case when input line lacks
//              the occupancy and temperature data, i.e. line is too short.    
    if( !linelengthflag) {
      if(strlen(line)<66) shortlineflag = 1;
    }
 
	// re-initialze the buffer on finding the  TER record
    if( !strncmp(line, "TER", 3) ) {
      l = strlen(line);
      for( i = 0; i < (85-l); i++ ) {
	line[l+i]=' ';
      }
      line[84]='\0';
    }
    
    strncpy( field1, line,    6 );
    strncpy( sr_no,  line+6,  5 );
    strncpy( aname,  line+11, 5 );
    strncpy( altloc, line+16, 1 );
    strncpy( rname,  line+17, 3 );
    strncpy( chain,  line+20, 2 );
    strncpy( res_no, line+22, 4 );
    strncpy( code,   line+26, 4 );
    strncpy( strx,   line+30, 8 );
    strncpy( stry,   line+38, 8 );
    strncpy( strz,   line+46, 8 );
// AJR 07.31.02 fix the case when occupancy and temp aren't provided.
    if(shortlineflag) {
      strcpy( occupancy, "  0.00");
      strcpy( temp_string, "  0.00");
    }
    else {
      strncpy( occupancy, line+54, 6 );
      strncpy( temp_string,   line+60, 6 );
    } 
 
    coord[1]    = atof(strx);
    coord[2]    = atof(stry);
    coord[3]    = atof(strz);
    occ         = atof(occupancy);
    temperature = atof(temp_string);
    ri_sr_no    = atoi(sr_no);

// AJR 07.31.02 temporary fix of nonstandard hydrogen names
    if( !strncmp(aname," H",2) && !strcmp(field1,"ATOM  ")) {
      strcpy(aname, "  H  ");
    }
// AJR 07.31.02 end temporary fix of nonstandard hydrogen names

// _____ SN 2007:02:08 	--- Change 'D'uetirium atoms to 'H'    
// _____ SN 2007:02:20	--- EXIT and output error to user asking to 
//			    replace all 'D'uetirium atoms to 'H'
    if( !strncmp(&aname[2],"D",1)) {
	strncpy(&aname[2],"H",1);
    }
    if( !strcmp(rname,"DOD")) {
        strcpy(rname,"HOH");
    }    
// ______________ End of fix ______________

    if( !strncmp(line,"ATOM",  4) || 
	!strncmp(line,"HETATM",6) ) {
      no_atoms += 1;
    }
    else {
      ri_sr_no = 0;
    }
  }
  return;
}
/**********************************************************************/
/**********************************************************************/

/**********************************************************************/
/* To read a line containing CONECT fron the .pdb file                */
/* The PDB file format allows crystallographers to identify explicitly*/
/* when specific atoms are to be connected. Search for these records  */
/* and add a connect record to the dataset file for each CONECT record*/
/**********************************************************************/
void conect::getconect(void) {

  char sr_no[6],conect1_sr_no[6],conect2_sr_no[6];
  char conect3_sr_no[6],conect4_sr_no[6];
  char line[90];
  int i,num[4], flag=0;

  ifil.getline(line,90);

  while( flag == 0 ) {
    
    if(strncmp(line,"CONECT",6)==0) {
      strncpy(sr_no,line+6,5);
      sr_no[5]='\0';
      strncpy(conect1_sr_no,line+11,5);
      conect1_sr_no[5]='\0';
      strncpy(conect2_sr_no,line+16,5);
      conect2_sr_no[5]='\0';
      strncpy(conect3_sr_no,line+21,5);
      conect3_sr_no[5]='\0';
      strncpy(conect4_sr_no,line+26,5);
      conect4_sr_no[5]='\0';
      
      ci_sr_no=atoi(sr_no);
      num[0]=atoi(conect1_sr_no);
      num[1]=atoi(conect2_sr_no);
      num[2]=atoi(conect3_sr_no);
      num[3]=atoi(conect4_sr_no);
      
      ci_sr_no=map_array[ci_sr_no];
      for(i=0;i<4;i++) {
	i_conect_sr_no[i]=map_array[num[i]];
      }
      
      /*i_conect1_sr_no=map_array[i_conect1_sr_no];
	i_conect2_sr_no=map_array[i_conect2_sr_no];
	i_conect3_sr_no=map_array[i_conect3_sr_no];
	i_conect4_sr_no=map_array[i_conect4_sr_no];*/
      
      if(ci_sr_no!=0) {
	ofil<<"CONECT"<<setw(5)<<ci_sr_no;
	for(i=0;i<4;i++) {
	  if(i_conect_sr_no[i]!=0) {
	    ofil<<setw(5)<<i_conect_sr_no[i];
	  }
	  else {
	    // ofil<<"     ";
	    break;
	  }
	}
	ofil<<"\n";
	
	/*	   ofil<<"CONECT"<<setw(5)<<i_sr_no<<setw(5)
		   <<i_conect1_sr_no<<setw(5)
		   <<i_conect2_sr_no<<setw(5)
		   <<i_conect3_sr_no
		   <<setw(5)<<i_conect4_sr_no<<"\n";*/
      }
      flag = 1;
    }
    else {
      ifil.getline(line,90);

      if( ifil.fail() ) {
#ifdef DEBUG
	cout << "end of file" << endl;
#endif
	ifil.close();
	ofil.close();
	break;
      }
    }

  }
  return;
} // end of getconnect method
/**********************************************************************/
/**********************************************************************/


/**********************************************************************/
//void node::add_rec_node(int mflag){ 
  void node::add_rec_node(void){ 

	r1.getrecord();
	//if(r1.modelflag){ cout<<"model is :"<<r1.modelflag<<endl;}
	return;
}
/**********************************************************************/


/**********************************************************************/
void conect_node::add_conect_node(void)
{
	c1.getconect();
	return;
}
/**********************************************************************/


/**********************************************************************/
/* Adds a record in the sorted linked list l1. Performs a check to    */
/* sort the incoming pdb data in order of increasing chain ID (a, b,..*/
/* residue number, and atom number.                                   */
/**********************************************************************/
void list::add_record(int usage) {

  int  flag2 = 1, file_open_flag = 0, flag1 = 0, i = 0, 
    length = 0, k=0, l=0, m=0, flag4 = 0;
  
  char buffer[3], chainletter[50], line[100];
  // choice is a class variable -- use a different variable name
  std::string my_choice;
  
  ofstream nh_file;

  ifstream fp;  

  if( file_open_flag == 0 ) {

    length = strlen(outputfile);
    
    if( length < 13 ) {
      cout << endl << "File name ERROR in makechem.cpp " << endl;
      exit(-1);
    }
    
    ifil.open(inputfile, ios::in);
	/* I don't understand why first line is being read here, as no
	 * further operation seems to happen on the buffer - Sameer, Mar 5, 04
	 * Commenting following two lines :
    	ifil.getline( line, 90 );
	    line[89] = '\0';
	 */

    if( ifil.fail() ) {
      cout << endl << "File ["<< inputfile << "] not found!" << endl;
      exit(-1);
    }

    // Read in previous conformations, if any
    fp.open(outputfile,ios::in);
    if( fp ) {
      i = 0;
      while(1) {
	fp.getline(line,100);
	if( fp == NULL) 
	  break;
	if( !strncmp(line, "REMARK:w:conformation_user", 26) ) {
	  fp.getline(line,80);
	  do {
	    k = strlen(line);
	    prev_sel[i]=line[k-1];
	    i++;
	    fp.getline(line,80);
	  } while(strncmp(line,"REMARK:W:conformation_user",26)==0);
	  prev_sel_flag=1;
	  break;
	}
      }
    }
    fp.close();
    ofil.open(outputfile,ios::out);
    file_open_flag=1;
  }
  
  map_array = new int[100000];
  std::fill(map_array, map_array + 100000, 0);

  //----Get a new line from the input file
  node *fresh, *temp, *prev;
  fresh = temp = prev = NULL;

  while( !ifil.eof() ) {

    fresh = new node;
    fresh->add_rec_node(); // read a line from the pdb file    
    if( fresh->r1.modelflag >= 2){
      model_flag = fresh->r1.modelflag;
      break;
    }

    fresh->next = fresh->prior = NULL;

    if( ifil.eof() ) // This line is necessary, as a new node is allocated before
      break;         // it is known whether the end-of-file has been reached. 

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  AJR 05.07.02 testing enum idea %%% */
    //cout<< fresh->r1.aname <<"\t"<<a1<<"\t";
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

    if( !strcmp(fresh->r1.chain,"  ") )	
      strcpy(fresh->r1.chain,"~~");
    
    /* Initialize the linked list */
    if( start == NULL )	{ 
      start = fresh;
      last  = fresh;
    }

    else {
      temp = start;
      while( flag2 != 0 ) {

	// DEBUG -ve residue indices
        //if( strncmp( temp->r1.res_no, "  -", 3 ) == 0 )
        //{
        //   printf("\n fresh: %s temp: %s strcmp: %d",fresh->r1.res_no,
	//			 temp->r1.res_no,
	//			 strcmp(fresh->r1.res_no,temp->r1.res_no));
        //}

	/************************************************************/
	/* If the last record and this record have different chain  */
	/* IDs.                                                     */
	/************************************************************/
	if( strcmp(fresh->r1.chain,temp->r1.chain) > 0 ) {
	  if( temp == start ) {
	    fresh->next=start;
	    start->prior=fresh;
	    start=fresh;
	    flag2=0;
	  }
	  else {
	    fresh->next=temp;
	    fresh->prior=temp->prior;
	    prev=temp->prior;
	    prev->next=fresh;
	    temp->prior=fresh;
	    flag2=0;
	  }
	}

	/************************************************************/
	/* If the last record and this record have THE SAME chain   */
	/* IDs.                                                     */
	/************************************************************/
	else {
	  if( strcmp(fresh->r1.chain,temp->r1.chain) == 0 ) {

	    /* SN 2008:07 	bugfix to handle negative rsd indices:
		 		Change '>' to '!=' as: " 1" > "-1" is FALSE */
	    if( strcmp(fresh->r1.res_no,temp->r1.res_no) != 0){
	      if( temp == start ) {
		fresh->next=start;
		start->prior=fresh;
		start=fresh;
		flag2=0;
	      }
	      else {
		fresh->next=temp;
		fresh->prior=temp->prior;
		prev=temp->prior;
		prev->next=fresh;
		temp->prior=fresh;
		flag2=0;
	      }
	    }
	    else {
	      if( strcmp(fresh->r1.res_no,temp->r1.res_no) == 0 ) {

		if( strcmp( fresh->r1.sr_no,temp->r1.sr_no) != 0 ) {
		  if( temp == start ) {
		    fresh->next=start;
		    start->prior=fresh;
		    start=fresh;
		    flag2=0;
		  }
		  else {
		    fresh->next=temp;
		    fresh->prior=temp->prior;
		    prev=temp->prior;
		    prev->next=fresh;
		    temp->prior=fresh;
		    flag2=0;
		  }
		}
		else {
		  if( temp == last ) {
		    last->next   = fresh;
		    fresh->prior = last;
		    last = fresh;
		    flag2 = 0;
		  }
		  else {
		    temp = temp->next;
		  }
		}
	      }
	      else {
		if( temp == last ) {
		  last->next=fresh;
		  fresh->prior=last;
		  last=fresh;
		  flag2=0;
		}
		else {
		  temp = temp->next;
		}
	      }
	    }
	  }
	  else {
	    if( temp == last ) {
	      last->next   = fresh;
	      fresh->prior = last;
	      last  = fresh;
	      flag2 = 0;
	    }
	    else {
	      temp = temp->next;
	    }
	  }
	}
      }
      flag2 = 1;
    }
  }

  /************************************************************
  the error produced by the SUN compiler (CC) orrur somewhere
    before the following error check. Seems to be including a 
    ghost line with all fields empty except ri_sr_no = 0. BMH 
  test_chem.open("test_chem", ios::out);

// A possible reason for the error is the class used to have a class variable 
// with the name temp and we have a local variable in this function with the 
// name temp.  Compiler warnings about shadowing the "temp" class variable 
// lend support to this hypothesis.

  temp = last;
  while( temp ) {
    test_chem << "check" << temp->r1.field1<<setw(5)<<temp->r1.ri_sr_no
	      << temp->r1.aname<<temp->r1.altloc
	      << temp->r1.rname<<temp->r1.chain
	      << temp->r1.res_no<<temp->r1.code
	      << temp->r1.strx<<temp->r1.stry
	      << temp->r1.strz<<temp->r1.occupancy
	      << temp->r1.temp_string<<setw(5)<<temp->r1.mod_res_no
	      << " "<<temp->r1.DAH_type
	      << endl;	
    temp = temp->prior;
  } 
  test_chem.close();
  ************************************************************/

  //--------------------------Check if any records are present
  if( !start && !last)	{
    cout<< endl << "\a\tThere were no atom records in the PDB file!"
	<< endl << endl;
    answer[0] = 'r';
    answer[1] = 'm';
    answer[2] = ' ';
    answer[3] = '\0';
    strcat(answer,outputfile);
    system( answer );
    exit(11);
  } 
  
  //--------------Marking the start of a new chain
  temp = last;
  chainptr[0] = temp;
  chainletter[0] = chainptr[0]->r1.chain[1];
  while( temp != start ) {
    temp=temp->prior;
    if(strcmp(temp->r1.chain,chainptr[chaincount]->r1.chain)!=0) {
      chainptr[chaincount+1]=temp;
      ++chaincount;
      chainletter[chaincount]=chainptr[chaincount]->r1.chain[1];
      
    }
  }

  //--------------------------------------------check for chain IDs that are U,V,..etc
  chain_change_count = 0;
  ++chaincount;
  m = 0;

  for( i = 0; i < chaincount; i++ ) {
    //cout<<endl<<chainletter[i];
    if( chainletter[i]=='U'||chainletter[i]=='V'||chainletter[i]=='W' ||
	chainletter[i]=='X'||chainletter[i]=='Y'||chainletter[i]=='Z') {
      for( k = 0; k < chaincount; k++ ) {
	flag4 = 0;
	for( l = 0; l < chaincount; l++ ) {
	  if( chainletter[l] == ('A'+m) ) {
	    flag4 = 1;
	    break;
	  }
	}
	if( flag4 != 1 ) {
	  
	  //----------------------replace chainletter[i] with 'A'+m
	  chain_change_count++;
	  chain_change[chain_change_count][0]=chainptr[i]->r1.chain[1];
	  chain_change[chain_change_count][1]=('A'+m);
	  
	  temp = chainptr[i];
	  strcpy( buffer, temp->r1.chain );
	  temp->r1.chain[1]=('A'+m);

	  /* rewrote, see below. BMH 3.16.00
	  while(1) {
	    temp = temp->prior;
	    if( !strcmp( temp->r1.chain, buffer ) ) {
	      temp->r1.chain[1]=('A'+m);
	      if( temp == start ) {
		break;
	      }
	    }
	    else {
	      break;
	    }
	  }
	  */
	  
	  temp = temp->prior;
	  while( (temp != NULL) && !strcmp( temp->r1.chain, buffer ) ) {
	    temp->r1.chain[1]=('A'+m);
	    temp = temp->prior;
	  }
	  
	  m++;
	  break;
	}
	else {
	  m++;
	}
      }
    }
  }
  
  i = 0;
  temp = chainptr[i];
  temp->r1.mod_res_no = 1;
  prev = temp;
  temp = temp->prior;

  // Modify the residue numbers	
  while( i < chaincount ) {

    if( !strcmp(temp->r1.chain,chainptr[i]->r1.chain) )	{
      if( !strcmp(temp->r1.res_no,prev->r1.res_no) ) {
	if( !strcmp(temp->r1.code,prev->r1.code) ) {
	  temp->r1.mod_res_no=prev->r1.mod_res_no;
	  if( temp == start ) {
	    break;
	  }
	  prev = temp;
	  temp = temp->prior;
	}
	else {
	  temp->r1.mod_res_no = prev->r1.mod_res_no+1;
	  if( temp == start ) {
	    break;
	  }
	  prev = temp;
	  temp=temp->prior;
	}
      }
      else {
	temp->r1.mod_res_no=prev->r1.mod_res_no+1;
	if( temp == start ) {
	  break;
	}
	prev = temp;
	temp = temp->prior;
      }
    }
    else {
      if( temp == start ) {
	break;
      }
      i++;
      temp=chainptr[i];
      temp->r1.mod_res_no=1;
      prev=temp;
      temp=temp->prior;
    }
  }
  temp = last;
  if( !strncmp(temp->r1.field1,"TER",3) ) {
    temp = temp->prior;
  }
  if( strcmp(temp->r1.altloc," ") ) {
    altlocptr[0]=temp;
    no_conf[altloccount]++;
    altloccount++;
    flag1 = 1;
  }	
  while( temp != start ) {
    prev = temp;
    temp = temp->prior;
    if( !strcmp(temp->r1.field1,"TER   ") ) {
      if( temp == start ) {
	break;
      }
      temp = temp->prior;
    }
    if( strcmp(temp->r1.altloc," ") ) {
      if( !flag1 ) {
	if( !strcmp(temp->prior->r1.aname,temp->r1.aname) ) {
	  altlocptr[0]=temp;
	  no_conf[altloccount]++;
	  altloccount++;
	  flag1=1;
	}
      }
      else {
	if( temp->r1.mod_res_no != altlocptr[altloccount-1]->r1.mod_res_no ) {
	  if( !strcmp(temp->prior->r1.aname,temp->r1.aname) ) {
	    altlocptr[altloccount]=temp;
	    no_conf[altloccount]++;
	    altloccount++;
	  }
	}
	else {
	  if( !strcmp(altlocptr[altloccount-1]->r1.aname,temp->r1.aname) ) {
	    no_conf[altloccount-1]++;
	  }
	}
      }
    }
  }
  
  if(altloccount != 0){

//---------------------------- AJR 04.29.02 adding default, "-non", option <-> usage = 1.
    cout<<endl<<"Multiple conformations were found for certain residues."<<endl;
    cout<<endl<<"Choose one of the following: "<<endl<<endl;
    cout<<"    1. Default (selects lowest mobility)"<<endl;
    cout<<"    2. New selection"<<endl;
    if( prev_sel_flag != 0 ) {
      cout<<"    3. Previous selection"<<endl;
    }
    if(usage == 0){
      cout<<endl<<"Enter Choice: ";
      cin>> my_choice;
    }
    else if(usage) {
      if( prev_sel_flag == 1) my_choice = "3";
      if( prev_sel_flag == 0) my_choice = "1";
    }

    if( prev_sel_flag == 0 ) {
      while(my_choice != "1" && my_choice != "2"){
	//	fflush(cin);
	//	cout<<endl<<endl<<"\t\t\tInvalid Choice!!"<<endl;
	cout<<endl<<"Enter choice: ";
	cin>> my_choice;
      }
      conf_chosen = new char[altloccount+1]; 
      del_conf(my_choice,usage);
    }
    else {
      while(my_choice != "1" && my_choice != "2" && my_choice != "3"){
	//	fflush(cin);
	//	cout<<endl<<endl<<"\t\t\tInvalid Choice!!"<<endl;
	cout<<endl<<"Enter choice: ";
	cin>> my_choice;
	//ccccccccccccccccccccccccccccccccccc
	}
      conf_chosen = new char[altloccount+1]; 
      del_conf(my_choice,usage);
    }
  }

/* _____ SN 2007:02 --- Check if the input file is missing even polar 'H' atoms
 			Note that all ATOM recs are stored in chainptr[0] with
			consecutive records linked by "prior" pointer and not 
			"next"!
	--- It is sufficient to check just the first residue's amide nitrogen */			
/* _____ SN 2008:02 --- Check all the chains!
			Also, make sure that only standard residue types are 
			tested for hydrogens as ligand atoms form a different 
			chain and usually are assigned non-standard rsd names!*/

  // The variable that used to be used here shadows the global variable
  // "residue *res" in ../../include/class.h

  // Cheat by using a map
  std::map<std::string, bool> std_res_names;
  std_res_names["GLY"] = true;
  std_res_names["ALA"] = true;
  std_res_names["VAL"] = true;
  std_res_names["LEU"] = true;
  std_res_names["ILE"] = true;
  std_res_names["SER"] = true;
  std_res_names["THR"] = true;
  std_res_names["CYS"] = true;
  std_res_names["MET"] = true;
  std_res_names["PRO"] = true;
  std_res_names["ASP"] = true;
  std_res_names["ASN"] = true;
  std_res_names["GLU"] = true;
  std_res_names["GLN"] = true;
  std_res_names["LYS"] = true;
  std_res_names["ARG"] = true;
  std_res_names["HIS"] = true;
  std_res_names["PHE"] = true;
  std_res_names["TYR"] = true;
  std_res_names["TRP"] = true;

  flag1 = 0;	// Count all chains that have polar Hs on their atoms
  while(flag1 < chaincount)
  {
    temp = chainptr[flag1]; 
    std::string curr_resName = temp->r1.rname;

    // Check if the residue name is a standard one
    if(std_res_names.find(curr_resName) != std_res_names.end()){
      while(temp != NULL && curr_resName == temp->r1.rname) {	
        if(temp->r1.aname[2] == 'H') {
           flag1++;	
	   break;		/* Break at the first occurrence of an 'H' */
        }
        temp = temp->prior;
      }
    }
    //Non-standard residue; goto next chain
    else flag1++;
  }

  if(flag1 != chaincount) {
     cout<<"\n################################################################";
     cout<<"\n#                                                              #";
     cout<<"\n# WARNING: Input file is missing hydrogens on residues.        #";
     cout<<"\n#          Please add and equilibrate hydrogens using software #";
     cout<<"\n#          such as AMBER, CHARMM, Reduce or WhatIf,            #";
     cout<<"\n#          delete the proflexdataset file, and re-run ProFlex. #";
     cout<<"\n#                                                              #";
     cout<<"\n#          Make sure the ligand and other co-factor atoms are  #";
     cout<<"\n#          appropriately protonated.                           #";
     cout<<"\n#                                                              #";
     cout<<"\n################################################################\n";
	//cout<<"\n\natom #: "<<temp->r1.ri_sr_no;
     exit(1);	   
   }

// _____ End of fix _____



  temp = last;
  flag1 = 0;

  while( temp != NULL ) {

// _____ SN 2007:02 --- fix to catch water molecules with no 'H's in input PDB file

     if(!strcmp(temp->r1.rname,"HOH") && temp->r1.aname[2] == 'O') {
        flag1++;
        if((temp->prior)->r1.aname[2] != 'H') {
     	  cout<<"\n################################################################";
     	  cout<<"\n#                                                              #";
     	  cout<<"\n# WARNING: Input file is missing hydrogens on water molecules  #";
     	  cout<<"\n#          Please add and equilibrate hydrogens using software #";
     	  cout<<"\n#          such as AMBER, CHARMM, Reduce or WhatIf,            #";
     	  cout<<"\n#          delete the proflexdataset file, and re-run ProFlex  #";
     	  cout<<"\n#                                                              #";
          cout<<"\n#          Make sure the ligand and other co-factor atoms are  #";
          cout<<"\n#          appropriately protonated.                           #";
          cout<<"\n#                                                              #";
     	  cout<<"\n################################################################\n";
          exit(1);
	}
     }
     else if(flag1 > 1) // Two water molecules checked and found 'H' atoms
			// Assuming the rest of the waters have 'H's we break out of check
	break;

    temp = temp->prior;
  }
// _____ End of fix _______


  int count = 1;
  temp = last;

  nh = 0;
  while( temp != NULL ) {

    if(temp->r1.aname[2]=='H') {
      nh++;
    }
    if(temp->r1.aname[2]=='S') {
      ns++;
    }
    
    if( strncmp( temp->r1.field1,"TER",3 ) ) {
      map_array[temp->r1.ri_sr_no] = count;
      temp->r1.ri_sr_no = count;
      if( temp != start )
	count++;
    }
    temp = temp->prior;
  }

  ofil.close();	
  ofil.clear();
}
/**********************************************************************/
/*               End of  void list::add_record(int )                  */
/**********************************************************************/


/**********************************************************************/
void list::del_conf(std::string choice_in, int usage) {

  int i, j, k, f = 0, occ_flag, ich;
//  int i, j, k, f = 0, flag = 0, occ_flag, temp_flag=0, ich;
  
  float temp_sum[20], max_occ, low_temp;

  char ans;
  std::string my_choice;
  
  node *temp, *prev, *t;
  
  ich = atoi(choice_in.c_str());
  switch(ich) {
  case 1:		//---------------------------------Default
    {
      // TESTING 	cout << "altloccount = " << altloccount << endl;
      for(i=0;i<altloccount;i++)
	{
	  occ_flag=0;                        //symmetry
	  temp=altlocptr[i];
	  
	  max_occ=temp->r1.occ;
	  t=temp;
				//occu[0]=max_occ;
	  for(j=1;j<no_conf[i];j++)	// Determine the conformation that has max. occupancy
	    {
	      temp=temp->prior;
	      /* TESTING
		 cout << "altloc = " << i << " conf. = " 
		 << j << " --> " << temp->r1.altloc << " on Res: " 
		 << temp->r1.rname << " Res. # = " << temp->r1.res_no
		 << " atom # = " << temp->r1.sr_no 
		 << " occ = " << temp->r1.occ << " max = " 
		 << max_occ << endl;
		 END OF TESTING */
	      //occu[i]=temp->r1.occ;
	      if(max_occ < temp->r1.occ)
		{
		  max_occ=temp->r1.occ;
		  t=temp;
		  temp=temp->prior;
		  // TESTING                      cout << "temp pointer changed!!" << endl;
		  occ_flag=1;          //asymmetry
		}
	      else
		{
		  // TESTING                      cout << "CONSTANT pointer." << endl;
		  if(max_occ!=temp->r1.occ)
		    {
		      // TESTING                      cout << "asymmetry." << endl;
		      occ_flag=1;  //asymmetry
		    }
		}
	    }
	  for(j=0;j<20;j++)
	    {
	      temp_sum[j]=0;
	    }
	  if(occ_flag==0)	// Determine the conformation that has highest mobility
	    {
	      temp=altlocptr[i];          //extra line
	      for(j=0;j<no_conf[i];j++)
		{
		  while(temp->r1.mod_res_no == altlocptr[i]->r1.mod_res_no)
		    {
		      if(strcmp(temp->r1.altloc,t->r1.altloc)==0)
			{
			  temp_sum[j]+=temp->r1.temperature;
                                /* TESTING
				   cout << "conf. = " << j << " --> " 
				   << temp->r1.altloc << " on Res: " 
				   << temp->r1.rname << " Res # = "
				   << temp->r1.res_no << " atom # = "
				   << temp->r1.sr_no << endl;
				   END OF TESTING */
			} 
		      temp=temp->prior;
		    }
		  t=t->prior;
		  temp=t;
		}
	      //aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
                                /* TESTING
				   for(j=0;j<no_conf[i];j++)
				   {
				   cout<<endl<<"Residue "<<i<<"conformation "<<j<<": sum of temp.= "<<temp_sum[j]<< endl;
				   }
				   END OF TESTING */
	      
	      low_temp=temp_sum[0];
	      k=0;
	      for(j=1;j<no_conf[i];j++)
		{
		  if(low_temp > temp_sum[j])
		    {
		      low_temp=temp_sum[j];
		      k=j;
//		      temp_flag=1;
		    }
		}
	      t=altlocptr[i];
	      for(j=0;j<k;j++)
		{
		  t=t->prior;
		}
	      temp=altlocptr[i];
	      prev=temp->prior;
	      while(temp->r1.mod_res_no == altlocptr[i]->r1.mod_res_no)
		{
		  if(strcmp(temp->r1.altloc,t->r1.altloc)!=0 && strcmp(temp->r1.altloc," ")!=0)
		    {
		      if(temp==altlocptr[i])
			{
			  altlocptr[i]=t;
			}
		      for(j=0;j<chaincount;j++)
			{
			  if(chainptr[j]==temp)
			    {
			      chainptr[j]=prev;
			      break;
			    }
			}
		      prev->next=temp->next;
		      temp->next->prior=prev;
		      delete temp;        
		    }
		  temp=prev;
		  prev=temp->prior;
		}
	    }
	  else
	    {
	      temp=altlocptr[i];
	      prev=temp->prior;
	      while(temp->r1.mod_res_no == altlocptr[i]->r1.mod_res_no)
		{
		  if(strcmp(temp->r1.altloc,t->r1.altloc)!=0 && strcmp(temp->r1.altloc," ")!=0)
		    {
		      if(temp==altlocptr[i])
			{
			  altlocptr[i]=t;
			}
		      for(j=0;j<chaincount;j++)
			{
			  if(chainptr[j]==temp)
			    {
			      chainptr[j]=prev;
			      break;
			    }
			}
		      prev->next=temp->next;
		      temp->next->prior=prev;
		      delete temp;
		    }
		  temp=prev;
		  prev=temp->prior;
		}
	    }
	  conf_chosen[i]=t->r1.altloc[0];
	  conf_option=1;
	}
      break;
    }
  
  case 3:		//-----------------------------------------------Previous selection
    {
      
      if(prev_sel_flag==0)
	{
	  //system("clear");
	  cout<<endl<<"No previous selections were made !!"<<endl;	
	  cout<<endl<<"Choose one of the following: "<<endl;
	  cout<<endl<<"    1. Default (selects lowest mobility)"<<endl;
	  cout<<"    2. New Selection"<<endl;
	  cout<<endl<<"Enter Choice: ";
	  cin>>my_choice;
	  if(usage) my_choice = "1";
	  while(my_choice != "1" && my_choice != "2")
	    {
	      //	cout<<endl<<endl<<"\t\t\tInvalid Choice!!"<<endl;
	      cout<<endl<<"Enter choice: ";
	      cin>>my_choice;
	    }
	  if(my_choice == "1" || my_choice == "2")
	    {
	      del_conf(my_choice,usage);
	      //flag = 1;
	      break;
	    }
	}
      else
	{
	  //system("clear");
	  cout<<endl<<endl<<"Your previous selections were: "<<endl;
	  cout<<endl<<"RESIDUE "<<"RESIDUE NO. "<<"CHAIN "
	      <<"CONFORMATIONS"<<" SELECTION"<<endl<<endl;
	  for(i=0;i<altloccount;i++)
	    {
	      //pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
	      /*
		cout<<"  "<<altlocptr[i]->r1.rname<<setw(10)
		<<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code<<"  "
		<<altlocptr[i]->r1.chain<<"       ";
		*/
	      
	      if(strcmp(altlocptr[i]->r1.chain,"~~")!=0)
		{
		  cout<<"  "<<altlocptr[i]->r1.rname<<setw(10)
		      <<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code<<"  "
		      <<altlocptr[i]->r1.chain<<"       ";
		}
	      else
		{
		  cout<<"  "<<altlocptr[i]->r1.rname<<setw(10)
		      <<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code<<"  "
		      <<"  "<<"       ";
		}
	      
	      temp=altlocptr[i];
	      for(j=0;j<no_conf[i];j++)
		{
		  cout<<" "<<temp->r1.altloc;
		  //   <<"          "
		  //   <<conf_ans[i];
		  temp=temp->prior;
		}
	      cout<<"          "<<prev_sel[i]<<endl;
	    }
	  cout<<endl<<endl<<"Type  n  for new selection: ";
	  cin >> answer;
	  if( strlen(answer) == 1 ) ans = answer[0];
	  else ans = 'Y';
	  if(ans!='n' && ans!='N')
	    {		
	      for(i=0;i<altloccount;i++)
		{
		  temp=altlocptr[i];
		  prev=temp->prior;
		  while(temp->r1.mod_res_no == altlocptr[i]->r1.mod_res_no)
		    {
		      if(temp->r1.altloc[0]!=prev_sel[i] && strcmp(temp->r1.altloc," ")!=0)
			{
			  if(temp==altlocptr[i])
			    {
			      altlocptr[i]=temp->prior;
			    }
			  for(j=0;j<chaincount;j++)
			    {
			      if(chainptr[j]==temp)
				{
				  chainptr[j]=prev;
				  break;
				}
			    }
			  prev->next=temp->next;
			  temp->next->prior=prev;
			  delete temp;
			}
		      temp=prev;
		      prev=temp->prior;
		    }
		  conf_chosen[i]=prev_sel[i];
		}
	      
	    }
	  else
	    {
	      //system("clear");
	      cout<<endl<<"Choose one of the following: "<<endl;
	      cout<<endl<<"    1. Default (selects lowest mobility)"<<endl;
	      cout<<"    2. New Selection"<<endl;
	      cout<<endl<<"Enter Choice: ";
	      cin>>my_choice;
	      if(usage) my_choice = "1";
	      while(my_choice != "1" && my_choice != "2")
		{
		  //	cout<<endl<<endl<<"\t\t\tInvalid Choice!!"<<endl;
		  cout<<endl<<"Enter choice: ";
		  cin>>my_choice;
		}
	      //	if(my_choice==1 || my_choice==2)
	      //	{
	      del_conf(my_choice,usage);
	      //flag=1;
	      break;
	      //	}
	    }
	}
      break;
    }
  
  case 2:		//----------------------------------------------New selection
    {
				//system("clear");
      
      cout<<endl<<"RESIDUE "<<"RESIDUE NO. "<<"CHAIN "<<"CONFORMATIONS"<<endl<<endl;
      for(i=0;i<altloccount;i++)
	{
	  /*
	    cout<<"  "<<altlocptr[i]->r1.rname<<setw(10)
	    <<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code<<"  "
	    <<altlocptr[i]->r1.chain<<"       ";
	    */
	  
	  if(strcmp(altlocptr[i]->r1.chain,"~~")==0)
	    {
	      cout<<"  "<<altlocptr[i]->r1.rname<<setw(10)
		  <<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code<<"  "
		  <<"  "<<"       ";
	    }
	  else
	    {
	      cout<<"  "<<altlocptr[i]->r1.rname<<setw(10)
		  <<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code<<"  "
		  <<altlocptr[i]->r1.chain<<"       ";
	    }
	  
	  temp=altlocptr[i];
	  for(j=0;j<no_conf[i];j++)
	    {
	      cout<<" "<<temp->r1.altloc;
	      //	strcpy(string[i],temp->r1.altloc);
	      temp=temp->prior;
	    }
	  cout<<endl;
	}
      
      for(i=0;i<altloccount;i++)
	{
	  temp=altlocptr[i];
	  f=0;
	  while(f==0)
	    {	
	      cout<<endl<<"Choose the conformation for residue "
		  <<altlocptr[i]->r1.res_no<<altlocptr[i]->r1.code<<" ("
		  <<altlocptr[i]->r1.rname<<"): ";
	      cin>>conf_chosen[i];
	      for(j=0;j<no_conf[i];j++)
		{
		  if(temp->r1.altloc[0]==conf_chosen[i])
		    {
		      f=1;
		      break;
		    }
		  temp=temp->prior;
		}
	      if(f==0)
		{
		  cout<<endl<<"\t Invalid selection!!"<<endl;
		}
	    }
	  conf_option=2;
	  f=0;
	}
      for(i=0;i<altloccount;i++)
	{
	  /*
	    temp=altlocptr[i];
	    prev=temp->prior;
	    while(temp->r1.mod_res_no == altlocptr[i]->r1.mod_res_no)
	    {
	    cout << endl << "temp = " << temp << " prior = " << temp->prior << " next = " << temp->next << endl;
	    cout << "temp:  " << "Res = " << temp->r1.rname << " Res # = " << temp->r1.res_no << " conf = " << temp->r1.altloc << " Atom = " << temp->r1.aname << " Atom # = " << temp->r1.sr_no << endl;
	    prev = temp->prior;
	    cout << "prior: " << "Res = " << prev->r1.rname << " Res # = " << prev->r1.res_no << " conf = " << prev->r1.altloc << " Atom = " << prev->r1.aname << " Atom # = " << prev->r1.sr_no << endl;
	    prev = temp->next;
	    cout << "next:  " << "Res = " << prev->r1.rname << " Res # = " << prev->r1.res_no << " conf = " << prev->r1.altloc << " Atom = " << prev->r1.aname << " Atom # = " << prev->r1.sr_no << endl;
	    prev=temp->prior;
	    temp=prev;
	    }
	    cout << endl << endl;
	    cin >> k;
	    */
	  temp=altlocptr[i];
	  prev=temp->prior;
	  while(temp->r1.mod_res_no == altlocptr[i]->r1.mod_res_no)
	    {
	      if(temp->r1.altloc[0]!=conf_chosen[i] && strcmp(temp->r1.altloc," ")!=0)
		{
		  if(temp==altlocptr[i])
		    {
		      altlocptr[i]=temp->prior;
		    }
		  for(j=0;j<chaincount;j++)
		    {
		      if(chainptr[j]==temp)
			{
			  chainptr[j]=prev;
			  break;
			}
		    }
		  prev->next=temp->next;
		  temp->next->prior=prev;
		  delete temp; 
		}
	      temp=prev;
	      prev=temp->prior;
	    }
	}
      break;
    }
  }
  /*
    cout << endl << endl;
    for(i=0;i<altloccount;i++)
    {
    cout << conf_chosen[i] << "  " << altlocptr[i]->r1.altloc << endl;
    }
    cin >> i;
    */
  return;
  
}
/**********************************************************************/
/*            End of void list::del_conf( char ch[] int usage)        */
/**********************************************************************/


/**********************************************************************/
/*  Add a CONECT record in the link list l2                           */
/**********************************************************************/
void conect_list::add_conect(void) {

  conect_node *fresh;
  
  /* BMH the code below is pointless. what the hell?
     if( file_open_flag == 0 ) {
     ifil.open(inputfile,ios::in);
     ofil.open(outputfile,ios::app);
     file_open_flag=1;
     }
  */
  
  ifil.open(inputfile,ios::in);
  ofil.open(outputfile,ios::app);

  while( ifil != NULL ) {
    
    fresh = new conect_node;
    
    fresh->add_conect_node();
    
    if( ifil == NULL )
      break;
    
    if( start == NULL ) {
      start = fresh;
      last  = fresh;
    }
    
    else {
      fresh->next=start;
      start->prior=fresh;
      start=fresh;
    }
    
  }

  ifil.close();
  ifil.clear();
  ofil.close();	
  ofil.clear();
  return;
}
/**********************************************************************/
/*          End of   void  connect_list::add_connect(void)            */
/**********************************************************************/


/**********************************************************************/
void list::make_chem_rec(void) {

  int flag3 = 0;
  
  node *tem;
	
  tem = last;
  ofstream xofil;
  
  xofil.open(outputfile,ios::app);

  if( xofil == NULL ) {
    cout<<"Error in opening outputfile"<<endl;
    exit(-1);
  }
  while( flag3 != 1 ) {
    if( tem == start ) {
      
      xofil << tem->r1.field1 << setw(5) << tem->r1.ri_sr_no
	    << tem->r1.aname  << tem->r1.altloc
	    << tem->r1.rname  << tem->r1.chain
	    << tem->r1.res_no << tem->r1.code
	    << tem->r1.strx   << tem->r1.stry
	    << tem->r1.strz   << tem->r1.occupancy
	    << tem->r1.temp_string << setw(5) << tem->r1.mod_res_no
	    << " " << tem->r1.DAH_type
	    << endl;	
      
      xofil.close();
      flag3 = 1;	
    }
    else {
      if(strncmp(tem->r1.field1,"TER",3) ){ // AJR 05.13.02 omit TER lines
      xofil << tem->r1.field1 << setw(5) << tem->r1.ri_sr_no
	    << tem->r1.aname  << tem->r1.altloc
	    << tem->r1.rname  << tem->r1.chain
	    << tem->r1.res_no << tem->r1.code
	    << tem->r1.strx   << tem->r1.stry
	    << tem->r1.strz   << tem->r1.occupancy
	    << tem->r1.temp_string << setw(5) << tem->r1.mod_res_no
	    << " " << tem->r1.DAH_type
	    << endl;	
      }
      tem = tem->prior;
    }
  }
}
/**********************************************************************/
/*              End of void list::make_chem_rec(void)                 */
/**********************************************************************/

		
/**********************************************************************/
/* Determine the largest and smallest residue number, and x, y, and z */
/* coordinate value.                                                  */
/**********************************************************************/
void list::cal_maxmin(void) {

  int flag3 = 0, max_res_count=1;

  node *tem;
  
  tem = last; // last is the first ATOM in the protein
  no_atoms = 0;

  while( flag3 != 1 ) {
    if( !strcmp(tem->r1.chain,"~~") ) {
      strcpy(tem->r1.chain,"  ");
    }
    if( tem == start ) {
      if(tem->r1.field1[0]!='T'){
	no_atoms++;
	//cout << tem->r1.sr_no << endl;
      }
      if(tem->r1.coord[1]<xmin)
	xmin=tem->r1.coord[1];
      if(tem->r1.coord[2]<ymin)
	ymin=tem->r1.coord[2];
      if(tem->r1.coord[3]<zmin)
	zmin=tem->r1.coord[3];
      
      if( largest_res < max_res_count ) {
	largest_res = max_res_count;
      }
      flag3 = 1;	
    }
    
    else {
      if( tem->r1.field1[0] != 'T' ){  // if the current list node is not a TER record in the pdb file
	no_atoms++;
	//cout << tem->r1.sr_no << endl;
      }
      if( tem->r1.coord[1] < xmin )  
	xmin = tem->r1.coord[1];
      if( tem->r1.coord[2] < ymin )
	ymin = tem->r1.coord[2];
      if( tem->r1.coord[3] < zmin )
	zmin = tem->r1.coord[3];
      
      if( tem->r1.mod_res_no == tem->prior->r1.mod_res_no )	{
	max_res_count++;
      }
      else {
	if( largest_res < max_res_count ) {
	  largest_res = max_res_count;
	}
	max_res_count = 1;
      }
      
      tem = tem->prior;
    }
  }
  //cout << "number of atoms" << no_atoms << endl;
}
/**********************************************************************/
/*             End of   void list::cal_maxmin(void)                   */
/**********************************************************************/



