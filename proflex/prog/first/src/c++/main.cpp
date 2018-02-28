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
/* main.cpp                                                           */
/*                                                                    */
/* This is the front end for the flexibility analysis algorithm FIRST */
/* The program reads in the coordinates of a protein structure, with  */
/* explicit hydrogen atoms, and computes the rigid and flexible reg-  */
/* ions in the structure. The initial data parsing is perfermed by    */
/* several C++ programs called by this routine, main.cpp. The actual  */
/* algorithm for performing the analysis is written in FORTRAN, and   */
/* may be found in ../fortran.                                        */
/*                                                                    */
/* This program has been developed by Don Jacobs, et al.              */
/* Additions and bugs fixes have been performed by                    */
/* AJ Rader, BM Hespenheide, and Sandeep Namilikonda at Michigan      */
/* State University.                                                  */
/*                    DISTRIBUTION VERSION  4.0                       */
/*                         revision: 05.13.02 by AJR                  */
/*                         revision: 02.07.03 by BMH                  */
/*                         revision: 02.13.04 by Sameer               */
/*								      */
/* 		V5.0  Last revision: 04.  .08 by Sandeep Namilikonda  */
/**********************************************************************/

using namespace std;
#include <vector>
#include <string>
#include <cstring>

#include"../../include/class.h"  // Contains the class descriptions and headers
#include"makechem.cpp"           // Make the proflexdataset file given a protein 
                                 //    structure in pdb format
#include"fixbabel.cpp"           // Creates the linked list describing the 
                                 //    connectivity of the protein atoms
#include"find_hbonds.cpp"        // Given the explicit hydrogen atoms, identify 
                                 //    all potential hydrogen bonds
#include"transitions.cpp"        // Read in the data from a proflexdataset file
#include"pick_hbonds.cpp"        // Selection of desired H-bonds
#include"hydrophobic.cpp"        // Finds and adds Hydrophobic tethers




/* ------- Usage Msg Function -----------------------------------------*/
void display_usage(char * program_name){
    //system( "clear" );
    cout<<"\n    Usage : "<< program_name<<" [-r] [-e] [-nonh|nonf] <-h|{-p|-pw}> <filename>"<<endl;
    // Adding hydrogens through FIRST is not yet supported - Sameer Arora Jan 7 '04
    //cout<<"\n    Usage is:  first [-ah][-h][-p][-pw][-res][filename]"<<endl;
    cout<<"\n    <filename> Expected to have a .pdb or .PDB file extension."<<endl;
    cout<<  "\t       With the -p or -pw option the extension may be _proflexdataset."<<endl;
    cout<<endl;
    cout<<"\n\t\aOption definitions:"<<endl;
    cout<<"\n\t-nonh  Use for non-interactive hydrogen bond dilution."<<endl;
    cout<<"\n\t-nonf  Use for non-interactive flexibility analysis."<<endl;
    cout<<"\n\t       **NOTE: Any filtering using stereochemical, energy,"; 
    cout<<"\n\t		hydrophobic, side chain, and water criteria from the"; 
    cout<<"\n\t		previous ProFlex run generating the proflex_dataset file"; 
    cout<<"\n\t		are not applied in -non mode.  Only the default stereochemical"; 
    cout<<"\n\t		criteria are applied. Be sure to only include buried"; 
    cout<<"\n\t		water molecules in the input PDB file.**"<<endl;
    //cout<<"\n\t-ah    Use to Add Hydrogen atoms explicitly to the PDB formated file."<<endl;
    //cout<<  "\t       This is the DEFAULT assumption."<<endl;
    //cout<<  "\t     > Hydrogen atoms that are currently listed are IGNORED, unless"<<endl;
    //cout<<  "\t     > the -h option is also used. (see -h option)"<<endl;
    cout<<"\n\t-h     Use with a PDB file having hydrogens."<<endl;
    cout<<  "\t       Hydrogens can be added explicitly by programs like WhatIf."<<endl;
    cout<<"\n\t-r     Distance between van der Wal's surfaces of atoms participating"<< endl;
    cout<<  "\t       in a hydrophobic interaction. 0.5 is default."<<endl;
    cout<<"\n\t-e     Max energy a hydrogen bond can have to be included in the analysis."<< endl;
    //cout<<  "\t       ( E.g. -0.1 for including only H-bonds with energy < -0.1 kcal/mol."<<endl;
    //cout<<  "\t        Default = -0.1 kcal/mol)"<<endl;
    //cout<<  "\t     > When used with the -ah option, the listed hydrogen atoms are"<<endl;
    //cout<<  "\t     > preserved with their coordinates unchanged in the optimization."<<endl;
    cout<<"\n\t-p     Use when the original PDB file has been previously processed."<<endl;
    cout<<  "\t       *proflexdataset* file from previous run must exist in the same directory."<<endl;
    cout<<"\n\t-pw    Same as [-p] except previous warning messages are re-displayed."<<endl;
    cout<< endl << endl << endl;
}
/* ------- Main Function ----------------------------------------------*/
int main(int argc, char *argv[]) {
  
  // Define function prototypes
  void getfiles(void);
  
  // Define variables	
  int 
    i, j, k, l, opt_sum, command_option[13];

  int 
    length, 
    correct_args_flag = 0,
    modusage = 0;
  
  float percent;
   
  char temp_res[4], ans;
  char line[90], *buffer;
  list l1;

  // Set the path variable to the root of the FIRST programs distribution. 
  // The user needs to setenv the variable, FIRSTPTB to point the
  // location of the /prog/ directory. 
  //path = getenv("FIRSTPTB");
  //if( !path ){
  //  cout << "The environment variable $FIRSTPTB is not set properly." << path << endl;
  //  return(1);
  //}

  // Set the path variable to the root of the ProFlex programs distribution.
  // The user needs to setenv the variable, PROFLEX_HOME to point the
  // location of the /prog/ directory.
  path = getenv("PROFLEX_HOME");
  if( !path ){
    cout<< "The environment variable $PROFLEX_HOME is not set properly."<< path << endl;
    return(1);
  }

  // define option utilities, display usage check
  if( argc == 1 ) {
	display_usage(argv[0]);
    exit( 1 );
  }
  
  // Initialize all command line settings to 0
  // Need to start from 0!
  for( i = 0; i < 13; i++) {
    command_option[i] = 0;
  }
  int usage = 0; // flag set 1 for 'best-effort' non-interactive execution

  // read and record the command line arguments 
  for( i = 1; i < argc; i++) {
    if(strcmp(argv[i],"-nonh") == 0) {
      command_option[0] =1;
      usage = 1;
      continue;
    }

    //the following option allows for automated flexibility analysis
    if(strcmp(argv[i],"-nonf") == 0) {
      command_option[0] =1;
      usage = 2;
      continue;
    }   
 
    if(strcmp(argv[i],"-h")  == 0 ) { 
      command_option[1] = 1; 
      continue; 
    }
    
    if(strcmp(argv[i],"-p")  == 0 ) { 
      command_option[2] = 1; 
      continue; 
    }

    if(strcmp(argv[i],"-pw") == 0 ) {
      command_option[2] = 1;
      command_option[3] = 1;
      continue;
    }
    if(!strncmp( argv[i], "-r", 2 ) ) {
      	// SN 2008:02 --- We never checked for input before!
	//		  Also, no validation of input is performed!
        if(*(argv[i]+2) == '\0')
	{
	  cout<<"\n\n\tERROR: Please enter a valid R_factor value!\n";
	  cout<<"\t       Note that the value should immediately follow '-r'\n\n";
	  exit(1);
        }
      R_factor = atof( (argv[i])+2 );
      continue;
    }

    if(!strncmp( argv[i], "-e", 2 ) ) {
      	// SN 2008:02 --- We never checked for input before!
	//		  Also, no validation of input is performed!
        if(*(argv[i]+2) == '\0')
	{
	  cout<<"\n\n\tERROR: Please enter a valid energy value!\n";
	  cout<<"\t       Note that the value should immediately follow '-e'\n\n";
	  exit(1);
        }
      max_energy = atof( (argv[i])+2 );
      continue;
    }

    length = strlen(argv[i]);
    
    if( length > 3 ) {
      j = length - 4;
      if( !strncmp(argv[i]+j, ".pdb", 4) || \
	  !strncmp(argv[i]+j,".PDB",4)){

          for( j = 0; j <= length; j++)
            inputfile[j] = *(argv[i]+j);
	
	/*
	 * If the input to proflex is in a different directory then
	 * the input argument has to be stripped of any path information
	 */
/*
	for( j = length-5; j >= 0; j-- )
	  if( *(argv[i]+j) == '/' )
	    break;
	
	if(j == -1)
	  for(j=0 ; j <= length; j++ )
	    inputfile[j] = *(argv[i]+j);
	else
	  for( k=j+1; k <= length; k++ )
	    inputfile[k-(j+1)] = *(argv[i]+k);
*/

	
	command_option[11] += 1; 
	correct_args_flag = 1;
	continue;
      }
    }
   
    /* 
     * 2006:02:28 SN	- Change the length to 15 = length(_proflexdataset)
     */ 
    if( length > 15) {
      j = length - 15;
      if( !strncmp(argv[i]+j, "_proflexdataset", 15) ) {

       for( j = 0; j <= length; j++)
          inputfile[j] = *(argv[i]+j);

	command_option[12] += 1;
	correct_args_flag = 1;
	continue;
      }
    }
    if( correct_args_flag == 0 )  {
    /* looks like incorrect args - Sameer, 13 Feb 2004*/
	display_usage(argv[0]);
    /* SN 2008:04 
	change the Warning message to indicate what is a correct input! */
	printf("ERROR: Invalid arguments! Please see above for a description of valid input");
		exit(-1);
    }
  }
  
  // check expected input file name
  opt_sum = command_option[11] + command_option[12];
  
  // if the used the -p flag, but supplied no file name. 
  if(opt_sum == 0) {
    if( command_option[2] == 1 ) {
      cout<<"\n\tNOTE: Use extension .pdb instead of _proflexdataset" << endl;
      cout<<"\t      even if the original .pdb file is not present.";
    }
    cout<<"\n\n\aPlease enter filename (with extension .pdb): ";
    cin >> inputfile;
    l = strlen(inputfile);
    if( strncmp( inputfile+(l-4), ".pdb", 4) ) {
      cout << "\n\n\aExpected file name must have a .pdb extension.\n\n" << endl;
      exit(2);
    }
    command_option[11] = 1;
  }
  
  else if(opt_sum > 1) {
    if( command_option[2] == 1 ) {
      cout<<"\n\tNOTE: Use extension .pdb instead of _proflexdataset" << endl;
      cout<<"\t      even if the original .pdb file is not present.";
    }
    cout<<"\n\n\aPlease enter ONLY ONE filename (with extension .pdb): ";
    cin >> inputfile;
    l = strlen(inputfile);
    if( strncmp( inputfile+(l-4), ".pdb", 4) ) {
      cout<<"\n\n\aExpected file name must have a .pdb extension.\n\n"<<endl;
      exit(2);
    }
    command_option[11] = 1;
    command_option[12] = 0;
  } 
  
  // If they supply the proflexdataset file but DID NOT use the -p option flag. 
  // Prompts for a file with extension .pdb . 
  else if( command_option[12] == 1 && command_option[2] == 0 ) {
    cout<<"\n\n\a\tThe extension _proflexdataset works only with the" << endl;
    cout<<"\t-p or -pw options after FIRST has previously been run.";
    cout<<"\n\nPlease enter filename (with extension .pdb):";
    exit(2);
  } 
  
  // process the -p and -pw options
  if(command_option[2] == 1) {
    command_option[1] = 1;
    command_option[10] = 0;
    for( i = 4; i <= 9; i++ )  
      command_option[i] = 0;
  }
  
  // process the -h and -ah options 
  if(command_option[1]  == 0 )
    command_option[10] = 1;  // the -ah is Default
  if(command_option[10] == 1 ) {
    cout<<"\n\n\aSorry; Currently this version cannot place hydrogen atoms."<<endl;
    cout<<    "       Place hydrogen atoms externally before running ProFlex."<<endl;
    cout<<    "       One such program that can be used is WhatIF."<<endl;
    cout<<    "       "<<endl;
    cout<<    "       Then run  ProFlex  again using the -h option.\n\n"<<endl;
    exit(3);
  }

  opt_sum = 0;
  percent = 0;        // Structure resolution info. query --- deprecated (v4.0)
  l1.resolution = percent;
  
  // define the outputfile. Append 'proflexdataset' to the root of the pdb file.

   /*
    * If the input PDB file is in a different directory then
    * the input argument has to be stripped of any path information.
    * Although, with '-p' option the proflexdataset has to be in
    * current working directory!
    */
   length = strlen(inputfile);
   for( j = length-1; j >= 0; j-- )
    if( inputfile[j] == '/' )			// Strip of path info.
    {
       strcpy(outputfile,(inputfile+j+1));
       break;
    }

   if(j == -1)
    strcpy(outputfile,inputfile);

  if( command_option[12] != 1 || command_option[2] != 1 ) {
    l = strlen(outputfile);
    outputfile[l-4]='\0';
    strcat(outputfile,"_proflexdataset");
  }
  
  //-----------------------------------------------------------------------------  
  //-----------------------------------------------------------------------------  
  // display the introduction screen, FIRST
  //system( "clear" );
  cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
  //system( "clear" );

  cout << "\n"
       <<  "\t               #    #   ####   #    #\n"
       <<  "\t               ##  ##  #       #    #\n"
       <<  "\t               # ## #   ####   #    #\n"
       <<  "\t               #    #       #  #    #\n"
       <<  "\t               #    #  #    #  #    #\n"
       <<  "\t               #    #   ####    ####\n"
	   <<  "\n"
       <<  "\t ######                  #######\n"
       <<  "\t #     #  #####    ####  #        #       ######  #    #\n"
       <<  "\t #     #  #    #  #    # #        #       #        #  #\n"
       <<  "\t ######   #    #  #    # #####    #       #####     ##\n"
       <<  "\t #        #####   #    # #        #       #         ##\n"
       <<  "\t #        #   #   #    # #        #       #        #  #\n"
       <<  "\t #        #    #   ####  #        ######  ######  #    #\n"
       <<  "\t----------------------------------------------------------\n"
       <<  "\t Software for Protein Flexibility Prediction and Analysis\n"
       <<  "\t                   version 5.0\n"
       <<  "\t    Copyright (C) 1997 - 2008, Michigan State University\n"
       <<  "\t                (formerly called FIRST)\n"
       <<  "\t     Protein Structural Analysis and Design Laboratory\n"
       <<  "\t                Michigan State University\n"                
       <<  "\t                   East Lansing MI, USA\n" << endl;
  
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------   
  
  
  // The -p option. Work with the previous proflexdataset
  if( command_option[2] == 1 ) {
    ifil.open( outputfile, ios::in );
    if( ifil == NULL ) {
      cout<< "\n\n\aPrevious data set ["<<outputfile<<"] not found!";
      cout<< "\n\aPlease make sure that the dataset file is in the "; 
      cout<< "current working directory.\n\n";
      ifil.close(); 
      exit(4);
    }
    else
      cout<<"\n\n\t Working with previous file: " << outputfile << endl;
    
    // The -pw option. Show previous warnings to the screen
    if( command_option[3] == 1 ) {
      i = -1;
      ifil.getline( line, 90 );
      while( strncmp(line,"END",3) != 0 ) {
	if( strncmp( line, "REMARK:w:", 9 ) == 0) {
	  cout << "\n\t\a WARNING MESSAGES PRESENT:\n" << endl;
	  i = 1;
	  goto Warning;
	}
	ifil.getline(line,90);
      }
    Warning:
      if( i < 0 ) {
	cout << "\n\n\t NO COVALENT BOND WARNING MESSAGES PRESENT.\n" << endl;
      }
      else {
	while( strncmp(line,"END",3) !=0) {
        /*
         * 2007:07:07   SN      --- Add warning msg to indicate Hbonds filtered
         *
          if(strncmp(line+6,":w:Hbnd_filtered",16)==0) {
            cout << "\t Some H-bonds were filtered out by the preferred"<< endl
 << "\t    length and angular checks. See 'filtered_Hbonds.log' for details."<<endl;
          }
	*/
	  if(strncmp(line+6,":w:check_bonding",16)==0) {
	    cout << "\t Some covalent bonds were identified by distance" << endl 
		 << "\t    criteria: See check_bonding in *proflexdataset* for details."<< endl;
	  }
	  if(strncmp(line+6,":w:isolated_atoms",17)==0) {
	    l = strlen(line);
	    for( j=l; j>25 ; j-- ) {
	      if( strncmp(line+j,"=",1) == 0 ) {
		i = j+1;
		break;
	      }
	    }
	    j = atoi(line+i);
	    cout << "\t Number of isolated atoms found = " 
		 << j << "." << endl;
	  }
	  if(strncmp(line+6,":w:isolated_water",17)==0) {
	    l = strlen(line);
	    for( j=l; j>25 ; j-- ) {
	      if( strncmp(line+j,"=",1) == 0 ) {
		i = j+1;
		break;
	      }
	    }
	    j = atoi(line+i);
	    cout << "\t Number of isolated oxygen atoms found in water = " 
		 << j <<  "." <<endl;
	  }
	  if(strncmp(line+6,":w:isolated_hydro",17)==0) {
	    l = strlen(line);
	    for( j=l; j>25 ; j-- ) {
	      if( strncmp(line+j,"=",1) == 0 ) {
		i = j+1;
		break;
	      }
	    }
	    j = atoi(line+i);
	    cout << "\t Number of isolated hydrogen atoms found = " 
		 << j << "." << endl;
	  }
	  if(strncmp(line+6,":w:disulfide_bond",17)==0) {
	    l = strlen(line);
	    for( j=l; j>25 ; j-- ) {
	      if( strncmp(line+j,"=",1) == 0 ) {
		i = j+1;
		break;
	      }
	    }
	    j = atoi(line+i);
	    cout << "\t Number of disulfide bonds found = " 
		 << j << "." << endl;
	  }
	  if(strncmp(line+6,":w:skip_nonhydrog",17)==0) {
	    cout << "\t Some heavy atoms were not connected to anything:" << endl;
	    cout << "\t    See skip_nonhydrogen in *proflexdataset* for details." 
		 << endl;
	  }
	  if(strncmp(line+6,":w:poor_bond",12)==0) {
	    cout << "\t Some [Expected] covalent bonds were too poor" << endl;
	    cout << "\t    to connect: See poor_bond in *proflexdataset* for details." 
		 << endl;
	  }
	  if(strncmp(line+6,":w:skip_hydrogen",16)==0) {
	    cout << "\t Some hydrogen atoms were not connected to anything:" << endl;
	    cout << "\t    See skip_hydrogen in *proflexdataset* for details." 
		 << endl;
	  }
	  if(strncmp(line+6,":w:missing_atoms",16)==0) {
	    cout << "\t Some [Expected] heavy atoms were missing from" << endl;
	    cout << "\t    known groups: See missing_atoms in *proflexdataset* for details." 
		 << endl;
	  }
	  if(strncmp(line+6,":w:chain_id",11)==0) {
	    cout << "\t Some chain IDs in the original PDB file were" << endl;
	    cout << "\t    re-labeled: See chain_id in *proflexdataset* for details" << endl;
	  }
	  if(strncmp(line+6,":w:conformation_auto",20)==0) {
	    cout << "\t Multiple conformations were present in the original" << endl;
	    cout << "\t    PDB file: Default option selected conformations" << endl;
	    cout << "\t    with maximum occupancy or minimum mobility." << endl;
	  }
	  if(strncmp(line+6,":w:conformation_user",20)==0) {
	    cout << "\t Multiple conformations were present in the original" << endl;
	    cout << "\t    PDB file: Each conformation was selected by user" << endl;
	    cout << "\t    See conformation_user in *proflexdataset* for details." << endl;
	  }
	  ifil.getline(line,90);
	}
      }
    }
    else {
      cout << "\n\tWARNING MESSAGES SUPPRESSED with -p option.\n" << endl;
    }

    ifil.close();
    
    l1.read_proflexdataset();
    
    l1.check();	
    
    getfiles();
    
  } // end of the code for the -p option
  

  // The -h option. Work with the pdb file in the absence of a _proflexdataset file
  else {
    // check to see if a _proflexdataset file exists
    ifil.open( outputfile, ios::in );
    if( !ifil ) {
      cout<< "\n\nProcessing the new data set ["<<inputfile<<"] --> ["
	  <<                             outputfile<<"]\n\n";
    }
    else { // if a *_proflexdataset file exists in this directory
      cout<<"\n\n\aCAUTION: Previous file [" << outputfile 
	  <<"] already exists!\n" << endl;
      cout<<"         This file will be over-written, but not its" << endl; 
      cout<<"         associated FIRST analysis output files. It is" << endl; 
      cout<<"         recommended to either rename " << inputfile << endl; 
      cout<<"         or remove the existing data files, or rename the " << endl; 
      cout<<"         existing data files. This is because the file with" << endl; 
      cout<<"         an extension  _proflexdataset serves as a record." << endl; 
      //cout<<"\n\a      ** The .ijkl files will be incremented while the" << endl; 
      //cout<<"         original input file will be over-written!\n" << endl; 
      cout<<"\n\a         Use [-p] to keep intact: [" << outputfile << "]" << endl << endl; 
      
      if( command_option[0] == 0 ) {
	cout <<"         Type  I  (default --> stop) to IGNORE this WARNING: ";
	ans = cin.get();
	if( ans == '\n' ){
	  cout << endl;
	  ifil.close(); 
	  exit(5);
	}
	else
	  cin.ignore(80, '\n');
      }	
      
    }
    
    ifil.close();
    ifil.clear();
    
    l1.add_record(usage);
    
    l1.cal_maxmin();

    ifil.close();
    ifil.clear();

    //------------------------------------------Read in the standard residue library
    residue_lib[0] = '\0';
    strcat( residue_lib,path );
    strcat( residue_lib, "/first/lib/residue.lib" );
    ifil.open( residue_lib, ios::in );

    if( ifil.fail() ) {
      cout<<endl<<"\aStandard Residue Library not found!"<< endl;
      cout<<"Expected location:  Path = " << residue_lib << endl;
      cout<<"Define new Path in class.h" << endl << endl << endl;
      ifil.close(); exit(6);
    }

    buffer = new char[201];
    while( !ifil.eof() ) {
      ifil.getline(buffer,200);
      res_count++;
    }

    delete [] buffer;

    res = new residue[res_count+1];	
    res_count=0;

    ifil.close();
    ifil.clear();
    ifil.open( residue_lib, ios::in );
    ifil.getline(line,90);
    
    ifil >> temp_res;
    while( strcmp(temp_res,"end") != 0 ) {
      res[res_count].get_residue(temp_res);
      res_count++;
      ifil >> temp_res;
    }

    ifil.close();
    ifil.clear();
    res_count++;

    // The following 50 or so lines replace the old code that tried to 
    // swap even if the locations were the same.  One issue that may have
    // caused problems is the classes do not have copy constructor or assignment
    // operators defined and copying an object to itself could be undefined.
    std::vector<std::string> residue_order;
    residue_order.push_back("GLY");
    residue_order.push_back("ALA");
    residue_order.push_back("VAL");
    residue_order.push_back("LEU");
    residue_order.push_back("ILE");
    residue_order.push_back("SER");
    residue_order.push_back("THR");
    residue_order.push_back("CYS");
    residue_order.push_back("MET");
    residue_order.push_back("PRO");
    residue_order.push_back("ASP");
    residue_order.push_back("ASN");
    residue_order.push_back("GLU");
    residue_order.push_back("GLN");
    residue_order.push_back("LYS");
    residue_order.push_back("ARG");
    residue_order.push_back("HIS");
    residue_order.push_back("PHE");
    residue_order.push_back("TYR");
    residue_order.push_back("TRP");

    k=0;
    for(size_t III = 0; III < residue_order.size(); ++III){
      for(size_t JJJ = III; JJJ < res_count; ++JJJ){
        // The jth residue in the res array does not have the same name as
        // the ith (current) residue in our prescribed order
        if(std::string(res[JJJ].rname) != residue_order[III]) continue;
   
        // The current residue (jth residue) is not in the correct slot in the
        // res array
        if(JJJ != III){
	  res[res_count] = res[III];
	  res[III] = res[JJJ]; 
	  res[JJJ] = res[res_count];
        }
        // If the current residue (jth residue) is in the correct slot in 
        // the res array, there is no need to move/copy anything
        ++k;
        break;
      }
    }

    if( k != 20 ) {
      cout << endl <<endl<< "\aERROR in reading the residue library: "<<residue_lib
	   << endl << "The standard residues were found to be incomplete!" 
	   << endl << endl;
      exit(7);
    }
    res_count--;

    /************************************************************/
    /* Read in the polar_lookup library.  Contains information  */
    /* on common HETATM groups, such as water. It lists the     */
    /* Hbonding properties of the atoms in each molecule listed.*/
    /************************************************************/
    polar_lookup[0] = '\0';
    strcat(polar_lookup,path);
    strcat(polar_lookup,"/first/lib/polar_lookup.lib");
    ifil.open(polar_lookup,ios::in);
    if( ifil.fail() ) {
      cout<<endl<<"\aPolar_lookup Library not found!"<<endl;
      cout<<"Expected location:  Path = " << path <<endl;
      cout<<"Define new Path in class.h"<<endl<<endl<<endl;
      ifil.close(); exit(6);
    }

    /************************************************************/
    /* Find out how many line are in the file polar_lookup.lib  */
    /* and create an array of that many objects of type DAHspecs*/
    buffer = new char[201];
    while( !ifil.eof() ) {
      ifil.getline( buffer, 200 );
      polar_count++;
    }
    delete [] buffer;
    polar = new DAHspecs[polar_count];	
    polar_count = 0;
    ifil.close();
    ifil.clear();
    /************************************************************/

    /************************************************************/
    /* For each line describing a molecule in the polar_lookup  */
    /* file, store the information.                             */
    ifil.open(polar_lookup,ios::in);
    ifil.getline(line,90);
    ifil >> temp_res;
    while( strcmp(temp_res,"end") != 0 ) {
      polar[polar_count].get_polar_record(temp_res);
      polar_count++;
      ifil >> temp_res;
    }
    ifil.close();
    ifil.clear();
    /************************************************************/
    
    l1.make_connectivity(usage); //--- 	SN 2006:03 --- This also finds 
				 //	and processes tethers
//cout<<"\nTEST1:"<<h_list_count;
    l1.find_Hbonds(usage);
//cout<<"\nTEST2\n";
    l1.check();	
    
    // Show warning messages
    i = -1;
    ifil.open( outputfile, ios::in );
    ifil.getline(line,90);
    while( strncmp(line,"END",3) ) {
      if(strncmp(line,"REMARK:w:",9)==0) {
	cout << "\n\t\a WARNING MESSAGES PRESENT:\n" << endl;
	i = 1;
	goto Warning2;
      }
      ifil.getline(line,90);
    }
    
  Warning2:
    if( i < 0 ) {
      cout << "\n\n\t NO WARNING MESSAGES PRESENT.\n" << endl;
    }
    else {
      while( strncmp(line,"END",3) !=0) {
	if(strncmp(line+6,":w:check_bonding",16)==0) {
	  cout << "\t Some covalent bonds were identified by distance" << endl 
	       << "\t    criteria: See check_bonding in *proflexdataset* for details." << endl;
	}
	if(strncmp(line+6,":w:isolated_atoms",17)==0) {
	  l = strlen(line);
	  for( j=l; j>25 ; j-- ) {
   	    if( strncmp(line+j,"=",1) == 0 ) {
	      i = j+1;
	      break;
	    }
	  }
	  j = atoi(line+i);
	  cout << "\t Number of isolated atoms found = " 
	       << j << "." << endl;
	}
	if(strncmp(line+6,":w:isolated_water",17)==0) {
	  l = strlen(line);
	  for( j=l; j>25 ; j-- ) {
	    if( strncmp(line+j,"=",1) == 0 ) {
	      i = j+1;
	      break;
	    }
	  }
	  j = atoi(line+i);
	  cout << "\t Number of isolated oxygen atoms found in water = " 
	       << j <<  "." <<endl;
	}
	if(strncmp(line+6,":w:isolated_hydro",17)==0) {
	  l = strlen(line);
	  for( j=l; j>25 ; j-- ) {
	    if( strncmp(line+j,"=",1) == 0 ) {
	      i = j+1;
	      break;
	    }
	  }
	  j = atoi(line+i);
	  cout << "\t Number of isolated hydrogen atoms found = " 
	       << j << "." << endl;
	}
	if(strncmp(line+6,":w:disulfide_bond",17)==0) {
	  l = strlen(line);
	  for( j=l; j>25 ; j-- ) {
	    if( strncmp(line+j,"=",1) == 0 ) {
	      i = j+1;
	      break;
	    }
	  }
	  j = atoi(line+i);
	  cout << "\t Number of disulfide bonds found = " 
	       << j << "." << endl;
	}
	if(strncmp(line+6,":w:skip_nonhydrog",17)==0) {
	   cout << "\t Some heavy atoms were not connected to anything:" << endl;
	   cout << "\t    See skip_nonhydrogen in *proflexdataset* for details." 
		<< endl;
	}
	if(strncmp(line+6,":w:poor_bond",12)==0) {
	  cout << "\t Some [Expected] covalent bonds were too poor" << endl;
	  cout << "\t    to connect: See poor_bond in *proflexdataset* for details." 
	       << endl;
	}
	if(strncmp(line+6,":w:skip_hydrogen",16)==0) {
	  cout << "\t Some hydrogen atoms were not connected to anything:" << endl;
	  cout << "\t    See skip_hydrogen in *proflexdataset* for details." 
	       << endl;
	}
	if(strncmp(line+6,":w:missing_atoms",16)==0) {
	  cout << "\t Some [Expected] heavy atoms were missing from" << endl;
	  cout << "\t    known groups: See missing_atoms in *proflexdataset* for details." 
	       << endl;
	}
	if(strncmp(line+6,":w:chain_id",11)==0) {
	  cout << "\t Some chain IDs in the original PDB file were" << endl;
	  cout << "\t    re-labeled: See chain_id in *proflexdataset* for details" << endl;
	}
	if(strncmp(line+6,":w:conformation_auto",20)==0) {
	  cout << "\t Multiple conformations were present in the original" << endl;
	  cout << "\t    PDB file: Default option selected conformations" << endl;
	  cout << "\t    with maximum occupancy or minimum mobility." << endl;
	}
	if(strncmp(line+6,":w:conformation_user",20)==0) {
	  cout << "\t Multiple conformations were present in the original" << endl;
	  cout << "\t    PDB file: Each conformation was selected by user" << endl;
	  cout << "\t    See conformation_user in *proflexdataset* for details." << endl;
	}
	ifil.getline(line,90);
      }
    }
    getfiles();
  }

  if( !usage ){
    // Mandatory pause
    cout << endl << "Type  \"s\"  to stop or any other key to continue: ";
    ans = cin.get();
    if( ans == 's' || ans == 'S') {
      cout << endl << endl;
      exit(5);
    }
    if( ans != '\n' )
      cin.ignore(80,'\n');
  }
  
  summon_fortran[0] = '\0';
  strcat( summon_fortran, path);
  strcat( summon_fortran, "/first/bin/first" );

  // 03.22.02 AJR allows correct hbond list output for -p or -h with -non option
  if( usage ) 
    modusage = usage + command_option[2];

  l1.pick_hbonds(modusage,command_option[2]); 	//--- 2006:03:15 SN 
						// preacptr_info will be read only is
						// '-p' is selected and user wants to do
						// stringent stereochemical filtering

  //------------------------------------------------------ summon the FORTRAN code 
  ofil.open( "qXyZaB.proflexdataset", ios::out );
  ofil << outputfile << endl;
  ofil << nfile << endl;
  ofil << usage <<endl;
  ofil.close();
 
  system( summon_fortran );
 // cout<<"Back to CPP";	SN --- DEBUGING
  
  qfil.open( "qXyZaB.proflexdataset", ios::in );
  qfil >> answer;
  
  if( strcmp(answer, "qXyZaB.proflexdataset") ) {
    strcpy( outfile[0], answer );
    for( i = 1; i < 7; i++ )
      qfil >> outfile[i];
    
    qfil >> nfile;
  }
  else {
    answer[0] = '\0';
    strcat(answer,"rm ");
    strcat(answer,outfile[0]);
    system( answer );
  }
  qfil.close();
}
