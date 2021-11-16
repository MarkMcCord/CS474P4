all: image.o ReadImage.o WriteImage.o  main
CFLAGS = -g -Wno-deprecated 


ReadImage.o:	image.h ReadImage.cpp
	g++ -c $(CFLAGS) ReadImage.cpp

WriteImage.o:	image.h WriteImage.cpp
	g++ -c $(CFLAGS) WriteImage.cpp

image.o:	image.h image.cpp
	g++ -c $(CFLAGS) image.cpp

main:	image.h image.o ReadImage.o WriteImage.o \
	main.cpp
	g++ -o main $(CFLAGS)  image.o ReadImage.o \
					WriteImage.o main.cpp

