#include "log.h"
#include "format.h"
#include "mio.h"



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  
  if (argc != 3) {
    usage ("%s <samples.txt> <name>",argv[0]);
  }
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    if (strEqual (currMatrix->name,argv[2])) {
      puts (mio_writeMatrix (currMatrix,0));
      break;
    }
  }
  mio_deInit ();
  return 0;
}
