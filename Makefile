#  MinGW
ifeq "$(OS)" "Windows_NT"
V_CFLG= -O3 -Wall
V_LIBS=-lglut32cu -lglu32 -lopengl32
CLEAN=del *.exe *.o *.a
else
#  OSX
ifeq ("$(shell uname)","Darwin")
V_CFLG=-g -O3 -Wall -Wno-deprecated-declarations $(shell sdl2-config --cflags)
V_LIBS=-framework OpenGL $(shell sdl2-config --libs)
S_CFLG=-g -O3 -Wall -Wno-deprecated-declarations
S_LIBS=-lm
#  Linux/Unix/Solaris
else
V_CFLG=-g -O3 -Wall $(shell sdl2-config --cflags) -DGL_GLEXT_PROTOTYPES
V_LIBS=-lGLU -lGL -lm $(shell sdl2-config --libs)
S_CFLG=-g -O3 -Wall 
S_LIBS=-lm 
endif
#  OSX/Linux/Unix/Solaris
CLEAN=rm -rf sim simg vis *.o *.a *.dSYM
endif

all:sim simg vis

#  Compile
#.c.o:
#	gcc -std=c99 -c $(V_CFLG) $<
#.cpp.o:
#	g++ -c $(V_CFLG) $<
vis.o:vis.cpp objects.h pixlight.vert pixlight.frag
	g++ -std=c++11 -c $(V_CFLG) $<

sim.o:sim.cpp
	g++ -std=c++11 -c $(S_CFLG) $< -fopenmp

simg.o:sim.cu
	nvcc -c -O3 -Xcompiler "-std=c++11 -c $(S_CFLG)" -o $@ $<

objects.o: objects.cpp objects.h
	g++ -c $(V_CFLG) $<

#  link
vis:vis.o objects.o
	g++ -g -O3 -o $@ $^ $(V_LIBS)

sim:sim.o
	g++ -g -O3 -o $@ $^ $(S_LIBS) -fopenmp

simg:simg.o
	nvcc -o $@ $^ -g -O3 $(S_LIBS)

#  Clean
clean:
	$(CLEAN)
