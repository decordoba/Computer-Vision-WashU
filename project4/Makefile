# Makefile for project 2

PROJ=svm
PROJ_OBJS=$(shell ls *.cpp | awk '{sub(".cpp", ".o"); print}') svmCMath.o
IMAGELIB=ImageLib/libImage.a

CC=g++

UNAME := $(shell uname)
# Flags for Linux
CPPFLAGS=-Wall -O3  -I /usr/include
LIB_PATH=-L/user/lib/x86_64-linux-gnu
LIBS= -lpng -ljpeg -lGL -lglut -lGLU -lX11 -lXft -lXinerama -lXext -ldl -lfontconfig -lfltk -lfltk_cairo -lfltk_forms -lfltk_gl -lfltk_images

all: $(PROJ)

$(IMAGELIB): 
	make -C ImageLib

$(PROJ): $(PROJ_OBJS) $(IMAGELIB)
	$(CC) -o $@ $(LIB_PATH) $(PROJ_OBJS) $(IMAGELIB) $(LIBS) 

clean:
	make -C ImageLib clean
	rm -f *.o *~ $(PROJ)
