#  Makefile to compile adp.cpp
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
#  To compile type
#  $ make
#  To run type
#  $ ./adp nSteps  (nSteps is the number of step the you want to run the simulation)


#-----------------------------------------------------------------------------
# Shell
#-----------------------------------------------------------------------------
SHELL = /bin/sh

#-----------------------------------------------------------------------------
#    Directories 
#-----------------------------------------------------------------------------
CC = g++
# where you installed mdpack library libmdpack.a
LIBDIR = ../../lib       
INC =  ../../include
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
SRCS =  adp.cpp
#-----------------------------------------------------------------------------
OBJS = $(SRCS:.cpp=.o)

# name of executable
PROG = adp

$(PROG):
	$(CC) $(CFLAGS)  -o $(PROG) $(SRCS) -I$(INC) -L$(LIBDIR) -lmdpack

clean:
	rm -f $(OBJS) $(PROG)
	
	


	