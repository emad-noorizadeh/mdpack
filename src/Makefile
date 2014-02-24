#  Makefile to compile mdpack
# 
#  Copyright (C) 2012  Emad Noorizadeh.
# 
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#   


#-----------------------------------------------------------------------------
# Shell
#-----------------------------------------------------------------------------
SHELL = /bin/sh

#-----------------------------------------------------------------------------
#    Directories 
#-----------------------------------------------------------------------------
OUT = ../lib/libmdpack.a
CC = g++
ODIR = ../obj
INC =  ../include
#-----------------------------------------------------------------------------
# Compilers
#-----------------------------------------------------------------------------
CC = g++
#-----------------------------------------------------------------------------
# Flags for Compilers
#-----------------------------------------------------------------------------
DEBUG = -g
OFLAG = -O3
CFLAGS = -Wall $(OFLAG) $(DEBUG)
#-----------------------------------------------------------------------------
# Source files
#-----------------------------------------------------------------------------
SRCS = 	vector3.cpp \
	atom.cpp \
	bond.cpp \
	angle.cpp \
	dihedral.cpp \
	impDihedral.cpp \
	utilities.cpp \
	ljTable.cpp  \
	random.cpp  \
	methods.cpp \
	md.cpp \
	
#-----------------------------------------------------------------------------
# Objects to be created
#-----------------------------------------------------------------------------
OBJS = $(SRCS:.cpp=.o)



libmdpack.a: $(OBJS)
	ar cru libmdpack.a $(OBJS)


# adp: $(OBJS)
# 	$(CC) $(CFLAGS) -o  $(PROG) $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS)  -c $*.cpp -I$(INC)
# $(SDIR)/%.o:
# 	$(CC) $(CFLAGS) -I$(INC) -c $*.cpp
	


	
	


# Dependencies are supposed to be in a file ``make.depend''
# which is inclduded by make.
include make.depend
# Because of the following rule, ``make'' will attempt to 
# create ``make.depend'' if it does not exist or if one
# of the files in $(SRCS) is more recent than ``make.depend''
make.depend:	$(SRCS)
		g++ -MM -I$(INC) $^ >$@	
		
		
clean:
	rm -f $(OBJS) $(PROG) make.depend
	rm -r ../obj
	rm  ../lib/libmdpack.a
	
install: libmdpack.a
	mv libmdpack.a ../lib/
	mkdir ../obj
	mv $(OBJS) ../obj/	

# depend:
# 	makedepend $(CFLAGS) $(SRCS)      
#        


# $(ODIR)/%.o: $(SDIR)/%.cpp 
#  
# 	$(CC) $(CFLAGS) -I$(INC) $(SRCS) 
# 
# $(OUT): $(OBJS) 
# 	ar rvs $(OUT) $^
# 
# .PHONY: clean
# 
# clean:
# 	rm -f $(ODIR)/*.o $(OUT)


	