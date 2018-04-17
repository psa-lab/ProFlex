# ProFlex


MSU ProFlex (formerly called FIRST) is a computational tool for identifying rigid and flexible regions in protein structures and protein-ligand complexes. Analysis of a single, static three-dimensional protein structure can capture the essential conformational flexibility of the protein.

Hydrogen bonds, salt bridges and hydrophobic contacts are identified by geometric and energetic criteria once hydrogens have been placed using an external program (e.g., WHAT IF (using the HB2NET command), AMBER, YASARA, or GROMACS). Using a constraint counting algorithm, all-atom calculations on proteins of over 1000 residues can be completed within seconds.

This software is available as source code and has been applied to many proteins and their complexes (see publications below).

- K. S. Keating, S. C. Flores, M. B. Gerstein, and L. A. Kuhn (2008) StoneHinge: Hinge Prediction by Network Analysis of Individual Protein Structures ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/StoneHinge_ProteinScience_in_press.pdf)), Protein Science.

- M. I. Zavodszky, M. Lei, A. R. Day, M. F. Thorpe, and L. A. Kuhn (2004) Modeling Correlated Main-chain Motions in Proteins for Flexible Molecular Recognition ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/Zavodszky_etal_Proteins04.pdf)) Proteins: Struct. Funct. Bioinf., 57(2), 243-261.

- H. Gohlke, L. A. Kuhn, and D. A. Case (2004) Change in Protein Flexibility Upon Complex Formation: Analysis of Ras-Raf Using Molecular Dynamics and a Molecular Framework Approach ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/Gohlke_etal_Proteins04.pdf)) Proteins: Struct. Funct. Bioinf., 56, 322-337.

- Brandon M. Hespenheide, A.J. Rader, M.F. Thorpe, and Leslie A. Kuhn (2002) Identifying protein folding cores from the evolution of flexible regions during unfolding ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/hespenheidejmgm2002.pdf)) J. Molec. Graphics and Modelling, 21, 195-207.

- A. J. Rader, B. M. Hespenheide, L. A. Kuhn, and M. F. Thorpe (2002) Protein Unfolding: Rigidity Lost ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/raderpnas2002.pdf)) Proceedings of the National Academy of Sciences USA 99, 3540-3545.

- D. J. Jacobs, A. J. Rader, L. A. Kuhn, and M. F. Thorpe (2001) Protein Flexibility Predictions Using Graph Theory ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/jacobsproteins2001.pdf)) Proteins: Structure, Function, and Genetics 44, 150-165.

## Supported Platforms

ProFlex is very resource-efficient and can analyze a protein complex within approximately a second. While other operating systems might be compatible, we recommened installing and running ProFlex on macOS or CentOS Linux. 


## Quick Installation Guide

This quick installation guide walks you through an example for how to install ProFlex on macOS 10.13 (High Sierra). Similar steps should apply to various Linux distributions like RedHat or CentOS but might be less straightforward. For general installation information, please refer to the documents linked in the **User Manual** section below.

1. Download ProFlex from GitHub by clicking the "Clone or Download" and then "Download ZIP" button in the upper right corner of this repository.
2. Unzip the downloaded `ProFlex-master.zip` file, open a new Terminal window, and navigate (`cd`) into the unzipped `ProFlex-master` directory and then `cd` into the `proflex/prog` subdirectory.
3. Execute `pwd` in the Terminal window and note the exact path of the `proflex/prog` subdirectory, e.g., `/Users/sebastian/Desktop/ProFlex-Master/proflex/prog`. This path is required for setting the `PROFLEX_HOME` environment variable of your Terminal's bash shell. You can do this by executing e.g., `export PROFLEX_HOME=/Users/sebastian/Desktop/ProFlex-Master/proflex/prog`.
4. In addition to setting the `PROFLEX_HOME` environment variable, it is recommended to also add the line (`export PROFLEX_HOME=/Users/sebastian/Desktop/ProFlex-Master/proflex/prog`) to your `~/.bash_profile` file so that you don't need to execute it each time you open a new Terminal window and want to run ProFlex.
5. Now, navigate from the `ProFlex-Master/proflex/prog` subdirectory back into the `ProFlex-Master/proflex` subdirectory and execute the command `make`. You might be prompted to install Apple's Developer Tools in case you haven't compiled code before. If you are prompted to do so, please go ahead and install these developer tools.
6. If you had to install the developer tools in step 5, attempt the compilation again by running `make`. In case you see an error regarding missing F77 compiler. If you don't have a Fortran compiler like `gfortran` installed, it's now time to do so (you can check if `gfortran` installed by running `which gfortran` in the Terminal). Binary installers for different versions of macOS are available form this website: https://gcc.gnu.org/wiki/GFortranBinaries.
7. If you execute `which gfortran` and see that a valid path is returned, this indicates that you have a version of GFortran installed. What you need to do next is to execute `export F77=gfortran` in the Terminal and execute the `make` command from the `ProFlex-Master/proflex` directory command one more time.
8. After the compilation succcessfully completed, the ProFlex executable will be available as `ProFlex-Master/proflex/bin/proflex`.

**Debugging Notes**

- If you successfully installed ProFlex but later get an `Segmentation Fault: 11` error when you try to execute it in a new Terminal window, this typically means that the `PROFLEX_HOME` isn't set in the current bash session. See steps 3 & 4 above to fix this issue.

## User Manual

- [A Quick Guide to ProFlex 5.1](proflex/docs/QuickGuideToProFlex.txt)

## License Agreements

ProFlex is now available for licensing to academic and commercial researchers.  
For academic use, please refer to the GNU/GPLv2 license that is available in this repository ([LICENSE.txt](LICENSE.txt)).

To arrange a commercial license or for scientific inquiries, please contact: 

Leslie Kuhn. 
502C Biochemistry Building   
MSU, East Lansing, MI 48824  

Telephone: (517) 353-8745  
E-mail: kuhnlab@msu.edu

## Release Notes

**Version 5.1 – Latest Release : January 2009.**

- ProFlex version 5.1 includes the implementation of a new rainbow color scheme for showing flexibility index results in PyMol. Please see [release_notes/version_5.1.txt](release_notes/version_5.1.txt) for a detailed list of changes made to the code for versions 5.0 and 5.1.

**Version 5.0 – Released : June 2008.**

- ProFlex version, 5.0, included significant enhancements over the previous release of ProFlex.

**Version 4.0 – Initial Release : Mar 2004.**

- Version 4.0 was the first public source code release of ProFlex.
