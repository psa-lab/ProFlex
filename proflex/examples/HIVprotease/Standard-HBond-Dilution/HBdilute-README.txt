This example shows the output of the hydrogen bond dilution (hbdilute) module 
of ProFlex. H-bond dilution module is executed only after the initial four
steps as detailed in "Flexibility-Analysis/README.txt". So, it can be run
with "-h" as well as "-p" options in which case the input would be 
<protein>.pdb and <protein>_proflexdataset, respectively. The output
comprises of a text file called "decomp_list", which encodes the rigid 
cluster decomposition of the input molecule. This text file is then processed
to generate a postscript file called the "<protein>_h-bonds.ps", which
graphically depicts the change in the flexibility profile of the protein as 
H-bonds that are stronger than the chosen energy cut-off are dissolved
one at a time, starting with the weakest one first. 

Apart from the above two files, a list of filtered H-bonds and a list of
hydrophobic tethers is output in two separate files.

When the fifth menu called the "ANALYSIS" menu in the run-transcript
is presented, the user should select option 2. The user is next
presented two modes of generating H-bond dilution plot. The recommended
option is option 1 or the standard H-bond dilution where the H-bonds are
dissolved in the increasing order of their strngth with the weakest bond
first. If the second option is chosen then the H-bonds are dissloved one
at a time in a random order based on an initial random number (seed) input
by the user. 

