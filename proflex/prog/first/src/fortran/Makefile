################################################################################
#  MSU ProFlex, formerly called FIRST, is a software developed to predict and  #
#  analyze protein flexibility.                                                #
#  This source file is a part of MSU ProFlex.                                  #
#                                                                              #
#  Copyright (C) 1997 - 2008, Michigan State University.                       #
#                                                                              #
#  This program is free software; you can redistribute to academic users only, #
#  it and/or modify it under the terms of the GNU General Public License,      #
#  version 2, as published by the Free Software Foundation.                    #
#                                                                              #
#  This program is distributed in the hope that it will be useful,             #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of              # 
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
#  GNU General Public License for more details.                                #
#                                                                              #
#  You should have received a copy of the GNU General Public License           #
#  along with this program; if not, write to the Free Software                 #
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA,  #
#  or see http://www.gnu.org/licenses/gpl.txt                                  #
################################################################################

# define your fortran compiler here.
ifndef F77
F77 = g77
endif

#path of binary directory
BINPATH = ../../bin
EXECUTABLE = first

# any flags?
FFLAGS = -g

$(BINPATH)/$(EXECUTABLE): clstrlabels.o colorclustr.o decomp_R3dbb.o decomp_H3dbb.o decomp_L3dbb.o \
	decomp_C3dbb.o first.o get_files.o read_data.o out_chime.o \
        mapmolecule.o name_change.o out_PDBfile.o out_decomp.o out_script.o \
        outputfirst.o pebgam_3dbb.o place_Hbond.o pruning.o read_hb_list.o \
        relabelatms.o molflex_3dbb.o hbdilute.o decomp_list.o hbfunction.o \
	hbpruning.o output_HBdilute.o output_bondwts.o set_parameter
	${F77} -O3 -o $(BINPATH)/$(EXECUTABLE) clstrlabels.o colorclustr.o decomp_R3dbb.o decomp_H3dbb.o decomp_L3dbb.o \
	decomp_C3dbb.o first.o get_files.o read_data.o out_chime.o \
        mapmolecule.o name_change.o out_PDBfile.o out_decomp.o out_script.o \
        outputfirst.o pebgam_3dbb.o place_Hbond.o pruning.o read_hb_list.o \
        relabelatms.o molflex_3dbb.o hbdilute.o decomp_list.o hbfunction.o \
	hbpruning.o output_HBdilute.o output_bondwts.o ${FFLAGS}

clstrlabels.o: clstrlabels.f set_parameter
	${F77} -O3 -c clstrlabels.f ${FFLAGS}

colorclustr.o: colorclustr.f set_parameter
	${F77} -O3 -c -Iset_parameter colorclustr.f ${FFLAGS}

decomp_R3dbb.o: decomp_R3dbb.f  set_parameter
	${F77} -O3 -c decomp_R3dbb.f ${FFLAGS}

decomp_H3dbb.o: decomp_H3dbb.f  set_parameter
	${F77} -O3 -c decomp_H3dbb.f ${FFLAGS}

decomp_L3dbb.o: decomp_L3dbb.f  set_parameter
	${F77} -O3 -c decomp_L3dbb.f ${FFLAGS}

decomp_C3dbb.o: decomp_C3dbb.f  set_parameter
	${F77} -O3 -c decomp_C3dbb.f ${FFLAGS}

first.o: first.f  set_parameter
	${F77} -O3 -c first.f ${FFLAGS}

get_files.o: get_files.f set_parameter
	${F77} -O3 -c get_files.f ${FFLAGS}

read_data.o: read_data.f set_parameter
	${F77} -O3 -c read_data.f ${FFLAGS}

mapmolecule.o: mapmolecule.f  set_parameter
	${F77} -O3 -c mapmolecule.f ${FFLAGS}

name_change.o: name_change.f set_parameter
	${F77} -O3 -c name_change.f ${FFLAGS}

out_PDBfile.o: out_PDBfile.f set_parameter
	${F77} -O3 -c out_PDBfile.f ${FFLAGS}

out_decomp.o: out_decomp.f set_parameter
	${F77} -O3 -c out_decomp.f ${FFLAGS}

out_chime.o: out_chime.f set_parameter
	${F77} -O3 -c out_chime.f ${FFLAGS}

out_script.o: out_script.f set_parameter
	${F77} -O3 -c out_script.f ${FFLAGS}

outputfirst.o: outputfirst.f  set_parameter
	${F77} -O3 -c outputfirst.f ${FFLAGS}

output_bondwts.o: output_bondwts.f  set_parameter
	${F77} -O3 -c output_bondwts.f ${FFLAGS} 

output_HBdilute.o: output_HBdilute.f  set_parameter
	${F77} -O3 -c output_HBdilute.f ${FFLAGS}

pebgam_3dbb.o: pebgam_3dbb.f  set_parameter
	${F77} -O3 -c pebgam_3dbb.f ${FFLAGS}

place_Hbond.o: place_Hbond.f  set_parameter
	${F77} -O3 -c place_Hbond.f ${FFLAGS}

pruning.o: pruning.f  set_parameter
	${F77} -O3 -c pruning.f ${FFLAGS}

hbpruning.o: hbpruning.f  set_parameter
	${F77} -O3 -c hbpruning.f ${FFLAGS}

read_hb_list.o: read_hb_list.f  set_parameter
	${F77} -O3 -c read_hb_list.f ${FFLAGS}

relabelatms.o: relabelatms.f  set_parameter
	${F77} -O3 -c relabelatms.f ${FFLAGS}

molflex_3dbb.o: molflex_3dbb.f set_parameter
	${F77} -O3 -c molflex_3dbb.f ${FFLAGS}

hbdilute.o: hbdilute.f set_parameter
	${F77} -O3 -c hbdilute.f ${FFLAGS}

decomp_list.o: decomp_list.f set_parameter
	${F77} -O3 -c decomp_list.f ${FFLAGS}

hbfunction.o: hbfunction.f set_parameter
	${F77} -O3 -c hbfunction.f ${FFLAGS}

clean:
	rm -f *.o $(BINPATH)/$(EXECUTABLE)
