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

#This is a top level makefile to execute the makefiles in the
# subdirectories.

all: bindir
	(cd ./prog/first/src/c++; make)
	(cd ./prog/first/src/fortran; make)
	(cd ./prog/hbdilute/src; make )
	(cd ./util/; make )
clean:
	(cd ./prog/first/src/c++; make clean)
	(cd ./prog/first/src/fortran; make clean)
	(cd ./prog/hbdilute/src; make clean)
	(cd ./util/; make clean )
	rm -rf ./bin

bindir:
	#creates a bin folder in the proflex directory that contains links to all executables
	mkdir -p bin
	ln -s ${PROFLEX_HOME}/first/bin/proflex bin/proflex
	ln -s ${PROFLEX_HOME}/hbdilute/bin/hbdilute bin/hbdilute
	ln -s ${PROFLEX_HOME}/../util/flex_index bin/flex_index

 
