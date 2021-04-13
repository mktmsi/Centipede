# makefile

#
# Definitions of main compilation tools
#

#For Mac?
CC            = gcc
CPPFLAGS      = -g -O3
#For Linux
#CXX            = g++
#CFLAGS        = -g -O3

#For Linux
#LDFLAGS		= 	-L/usr/X11R6/lib -L/usr/X11/lib -L/usr/lib/X11R6 -L/usr/lib/X11
#LDLIBS		= -lpng -lglut -lGLU -lGL -lXmu -lXi -lXext -lX11 -lm 

#For Mac
LDFLAGS     =   -I/usr/X11R6/include -L/usr/X11R6/lib -framework OpenGL -framework GLUT -framework Foundation
LDLIBS      =   -lpng -lGLU -lGL -lXmu -lXi -lXext -lX11 -lm -lstdc++ 

PROGRAM1      = main
#OBJECTS       = 

#.SUFFIXES   : .o .c .cpp .h

$(PROGRAM1) : $(PROGRAM1).o Monitor.o

clean      :

	-rm -rf obj
	-rm -f *.o core*
	-rm -f *~
	@echo "all cleared"
