#include "log.h"
#include "format.h"
#include "mio.h"
#include <stdio.h>
#include <unistd.h>



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  Stringa buffer;
  FILE *fp;

  if (argc != 3) {
    usage ("%s <samples.txt> <pathOutputDirectory>",argv[0]);
  }
  buffer = stringCreate (100);
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    stringPrintf (buffer,"mkdir %s/%s",argv[2],currMatrix->name);
    hlr_system (string (buffer),0);
    stringPrintf (buffer,"%s/%s/%s.matrix",argv[2],currMatrix->name,currMatrix->name);
    fp = fopen (string (buffer),"w");
    if (fp == NULL) {
      die ("Unable to open file: %s",string (buffer));
    }
    fprintf (fp,"%s\n",mio_writeMatrix (currMatrix,0));
    fclose (fp);
    stringPrintf (buffer,"%s/%s",argv[2],currMatrix->name);
    chdir (string (buffer));
    stringPrintf (buffer,"mio2images ../%s < %s/%s/%s.matrix",argv[1],argv[2],currMatrix->name,currMatrix->name);
    hlr_system (string (buffer),0);
    stringPrintf (buffer,"mio2bed ../%s < %s/%s/%s.matrix",argv[1],argv[2],currMatrix->name,currMatrix->name);
    hlr_system (string (buffer),0);
  }
  mio_deInit ();
  return 0;
}
