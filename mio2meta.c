#include "log.h"
#include "format.h"
#include "mio.h"



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  
  if (argc != 2) {
    usage ("%s <samples.txt>",argv[0]);
  }
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    puts (mio_writeMatrix (currMatrix,1));
  }
  mio_deInit ();
  return 0;
}
