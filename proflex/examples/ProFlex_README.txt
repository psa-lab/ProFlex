We have included ProFlex runs for two proteins in this directory to guide
the user through the process of running ProFlex and familiarizing with the 
input-output files and their contents. Here is how each directory has been 
organized:

				<Protein>
		    		    |
	 	 ---------------------------------------
		|					|
	Flexibility-Analysis			Standard-HBond-Dilution
		|					|
    	    ---------				   ----------
   	   |	     |		    	      	  |	     |
	 Input	   Output			Input      Output

The two directories under each protein correspond to two different
modes of analysis possible using ProFlex: 
	(a) flexibility analysis at a particular energy cutoff followed by 
a post-processing step that generates visualization scripts for viewing the 
protein color-coded with respect to rigid cluster decomposition, flexible 
cluster decomposition, or by flexibility index; 
	(b) HBdilute flexibility analysis done through gradual dissolution 
of H-bonds from the weakest to the strongest of the H-bonds that remain after 
the initial user-specified (or default) energy cutoff is applied. This is 
done via the hbdilute option in ProFlex, providing output as a postscript 
file that contains a color coded graph depicting how the protein residues 
map onto distinct rigid clusters.

The workflow of ProFlex is as follows:

	(a) An input PDB file containing polar hydrogen atoms is provided 
to ProFlex, only including those water molecules and bound ligands or metals 
the user intends to include in the flexibility analysis.

	(b) ProFlex queries the user about the hydrogen-bond potential of 
non-protein atom (e.g., metals or ligand atoms) to assign donor, acceptor, 
donor/acceptor (e.g., hydroxyl group), or "none" as the H-bonding character.

	(c) If an atom appears to have improper covalent or non-covalent 
valence, ProFlex will query about which atoms to be considered as bonded to 
an atom. Note that a bond (including hydrogen-bond) is either present or 
absent in ProFlex, rather than having partial-bond character. Stereochemical 
and energy criteria are intended to only include H-bonds and salt bridges 
that act as constraints in the structure. 

	(d) stereochemical filtering of H-bonds to ensure that only 
good-geometry hydrogen bonds and salt bridges are included.

	(e) filtering to only include certain classes of interactions 
(hydrophobic tethers, H-bonds involving waters, side-chains, etc.), as 
desired on a case-by-case basis.

	(f) energy-based filtering of H-bonds, using the default of 
-1.0 Kcal/mol as the weakest hydrogen bond or salt bridge to include, or 
more stringent filtering as required for that particular protein (e.g., 
if the protein appears entirely rigid at the default cutoff, and relative 
flexibility information is desired).

Under each protein directory, the input PDB file to ProFlex is specified in 
the directory named "INPUT" and the corresponding output from the ProFlex run 
is found in the directory named "OUTPUT". We also provide "README.txt" files in
each directory for further guidance, and "run-transcript.txt" as a transcript 
of the actual ProFlex run used to create this data.  This may help the user 
understand what information is required at each step the run, and provide an 
example of an appropriate response.



