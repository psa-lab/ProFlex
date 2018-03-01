Before running proflex, please pre-process the protein input file to 
contain polar hydrogen atoms (see Quick Guide).  This has already been 
done for the 1bck.pdb file provided in the INPUT directory of this example.

Next, run the proflex command:

	proflex -h 1bck.pdb

In the interactive proflex run, you are presented the default 
stereochemical validation criteria with a choice to further tighten the 
constraints. In the example run (see transcript file), the user chose to 
enter more stringent constraints. Once the contraints are updated, the H-
bond list is accordingly filtered and the new count of H-bonds is 
displayed on the screen.

Next, proflex presents the NON-COVALENT BOND SCREENING OPTIONS MENU. (See 
the Quick Guide and README file in the top level of the examples 
directory for guidance here.)

In the fourth step, proflex implements an energy-based filter. The user 
is presented with the range of energies at which H-bonds or salt bridges 
were formed in the input molecule. To choose an energy cutoff, a user has 
to simply enter the magnitude of the energy, e.g., -2.3 (see the run-
transcript).

The fifth menu provides the user a choice to perform flexibility analysis
or to perform an H-bond dilution procedure.


			OUTPUT DESCRIPTION
			------------------

As described in the QuickGuide, the ProFlex output consists of various
textfiles and PyMol scripts generated based on the contents of text 
files. Here we describe the three PyMol scripts in detail.

The output displayed by PyMol upon running the three ".pml" scripts is 
stored as the following images:

1bck_FC_tubeview.png			- 1bck_FC_0001.pml
1bck_RC_ribbonview.png			- 1bck_RC_0001.pml
1bck_flexindex_pymol_cartoonview.png	- 1bck_flex_0001.pml 
												
(1) 1bck_flex_0001.pml uses the 1bck_flex.pdb as the input. For each 
atom, the input file contains a relative flexibility index on the scale 
of 0 to 99, with 0 representing maximum rigidity and 99 for maximum 
flexibility in the B-value column of the *flex.pdb file. The script 
defines a color code that maps the B-values over a range of 0 to 99 onto 
11 bins: 

0-42 	colored blue for max rigidity
42-49 	four bins with decreasing intensity of blue
49-50 	colored grey to represent isostatic regions
50-58 	four bins colored starting with yellow and progressing 	
	towards red as the flexibility increases
58-99 	colored red for max flexibility

(2) 1bck_RC_0001.pml again uses 1bck_flex.pdb as its input.  The atom 
number column contains a re-ordered atom index such that atoms that 
belong to the same flexible or rigid cluster have contiguous atom 
indices. This allows the script to specify the cluster boundaries easily 
to PyMol. In this script, each rigid cluster with 7 or more atoms is 
distinctly named as "RC<number>" although the number does not reflect the 
rank of the cluster with respect to the cluster size. Hence, when the 
script is run in PyMol, the user may notice non-contiguous RC<number> 
objects. All the rigid clusters that have fewer than 7 atoms are grouped 
into objects named "smallRCs<number>" where each "smallRCs" object may 
consist of at most 50 small, individual clusters.
(The size 50 is used as the limit to prevent inadvertent segmentation 
fault that may result due to PyMol's implicit limit on the string data 
type.) Since the emphasis here is on displaying the rigid clusters, all 
the flexible clusters are grouped into one single "FC" object and colored 
white, whereas each RC is given a distinct color until the total number 
of RCs exceeds 20, in which case, colors are reused. Hydrophobic 
interactions are represented as tiny spheres and grouped under the object 
name "HPHOB<number>".

(3) 1bck_FC_0001.pml is similar to the _RC_ script except that the
emphasis is on flexible clusters. Each flexible cluster with a size of 
at least 7 atoms is named "FC<number>" and colored with a distinct pastel 
color until there are 20 clusters and then colors are reused. The smaller
cluster groups are called "smallFCs<number>". The rigid clusters are 
grouped under one object called "RC" and colored blue. Hydrophobic 
interactions are represented as tiny spheres and grouped under the object 
name "HPHOB<number>". In addition to all the above objects, this script
creates another object called "dangle" to represent the dangling bonds 
that bear a zero cluster label in the "allbond" file. These bonds are 
flexible and not coupled to other regions in the structure, and are 
colored white. 

NOTE: 
	Information in the "allbond" file is used to generate the 
"_flex.pdb" file from the "proflexdataset" file.  The contents of the 
allbond file are explained in the Quick Guide, and may be analyzed by 
users independent of the PyMol interface presented here to display rigid 
clusters, flexible clusters, and relative flexibility indices.  Remember 
that atoms in the flex.pdb file have been renumbered, and thus the atom 
numbers in the PDB file that was input to proflex should be used for 
interpreting the allbonds file. 

