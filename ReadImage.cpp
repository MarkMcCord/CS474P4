#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "image.h"

void readImage(char fname[], ImageType& image)
{
 int i, j;
 int N, M, Q;
 unsigned char *charImage;
 char header [100], *ptr;
 ifstream ifp;

 ifp.open(fname, ios::in | ios::binary);

 if (!ifp) {
   cout << "Can't read image: " << fname << endl;
   exit(1);
 }

 // read header

 ifp.getline(header,100,'\n'); //P5 thing
 if ( (header[0]!='P') ||    /* 'P' */
      (header[1]!='5') ) {   /* '5' */
      cout << "Image " << fname << " is not PGM" << endl;
      exit(1);
 }
 //cout << header << endl;

 ifp.getline(header,100,'\n'); //comment
 //cout << header << endl;

 ifp.getline(header,100,'\n'); //size
 //cout << header << endl;
 N=strtol(header ,&ptr,0);
 M=strtol(ptr,&ptr,0);

 ifp.getline(header,100,'\n'); //range
 //cout << header << endl;
 Q=strtol(header,&ptr,0);

 //cout << M << " " << N << " " << Q << endl;
 charImage = (unsigned char *) new unsigned char [M*N];

 ifp.read( reinterpret_cast<char *>(charImage), (M*N)*sizeof(unsigned char));

 if (ifp.fail()) {
   cout << "Image " << fname << " has wrong size" << endl;
   exit(1);
 }

 ifp.close();

 //
 // Convert the unsigned characters to integers
 //

 int val;

 for(i=0; i<M; i++)
   for(j=0; j<N; j++) {
     val = (int)charImage[i*N+j];
     image.setPixelVal(i, j, val);     
   }

 delete [] charImage;

}
