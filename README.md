# ProFlex


MSU ProFlex (formerly called FIRST) is a computational tool for identifying rigid and flexible regions in protein structures and protein-ligand complexes. Analysis of a single, static three-dimensional protein structure can capture the essential conformational flexibility of the protein.

Hydrogen bonds, salt bridges and hydrophobic contacts are identified by geometric and energetic criteria once polar hydrogens have been placed using an external program (such as What If). Using a constraint counting algorithm, all-atom calculations on proteins of over 1000 residues can be completed within seconds.

This software is available as source code and has been applied to many proteins and their complexes (see publications below).

- K. S. Keating, S. C. Flores, M. B. Gerstein, and L. A. Kuhn (2008) StoneHinge: Hinge Prediction by Network Analysis of Individual Protein Structures ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/StoneHinge_ProteinScience_in_press.pdf)), Protein Science.

- M. I. Zavodszky, M. Lei, A. R. Day, M. F. Thorpe, and L. A. Kuhn (2004) Modeling Correlated Main-chain Motions in Proteins for Flexible Molecular Recognition ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/Zavodszky_etal_Proteins04.pdf)) Proteins: Struct. Funct. Bioinf., 57(2), 243-261.

- H. Gohlke, L. A. Kuhn, and D. A. Case (2004) Change in Protein Flexibility Upon Complex Formation: Analysis of Ras-Raf Using Molecular Dynamics and a Molecular Framework Approach ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/Gohlke_etal_Proteins04.pdf)) Proteins: Struct. Funct. Bioinf., 56, 322-337.

- Brandon M. Hespenheide, A.J. Rader, M.F. Thorpe, and Leslie A. Kuhn (2002) Identifying protein folding cores from the evolution of flexible regions during unfolding ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/hespenheidejmgm2002.pdf)) J. Molec. Graphics and Modelling, 21, 195-207.

- A. J. Rader, B. M. Hespenheide, L. A. Kuhn, and M. F. Thorpe (2002) Protein Unfolding: Rigidity Lost ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/raderpnas2002.pdf)) Proceedings of the National Academy of Sciences USA 99, 3540-3545.

- D. J. Jacobs, A. J. Rader, L. A. Kuhn, and M. F. Thorpe (2001) Protein Flexibility Predictions Using Graph Theory ([pdf](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/jacobsproteins2001.pdf)) Proteins: Structure, Function, and Genetics 44, 150-165.

## User Manual

An updated version (v5.0) of the manual will be available shortly. Meanwhile, version 4.0 of the manual can be downloaded by clicking the link below. Note that the background on ProFlex in this manual will be useful, but many of the aspects of running the software, including output file formats and visualization, have changed. For details of running the code, please refer to the Quick Guide (for v. 5.1) distributed with the software (and also made available for download below).

- [ProFlex / FIRST Manual 4.0](http://www.kuhnlab.bmb.msu.edu/projects/first/docs/FIRST_Manual.pdf) (pdf format)
- [A Quick Guide to ProFlex 5.1](http://www.kuhnlab.bmb.msu.edu/projects/first/docs/A_Quick_Guide_to_ProFlex.pdf) (pdf format)

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
