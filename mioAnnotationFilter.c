#include "log.h"
#include "format.h"
#include "mio.h"



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  int i;
  Annotation *currAnnotation;
  int mod;

  if (argc != 4) {
    usage ("%s <samples.txt> <first|second> <nameAnnotationSet>",argv[0]);
  }
  mio_init ("-",argv[1]);
  if (strCaseEqual (argv[2],"first")) {
    mod = 1;
  }
  else if (strCaseEqual (argv[2],"second")) {
    mod = 0;
  }
  else {
    usage ("%s <samples.txt> <first|second> <nameAnnotationSet>",argv[0]);
  }
  while (currMatrix = mio_getNextMatrix ()) {
    i = 0; 
    while (i < arrayMax (currMatrix->annotations)) {
      currAnnotation = arrp (currMatrix->annotations,i,Annotation);
      if ((currAnnotation->intervalNumber % 2) == mod &&
          strCaseEqual (argv[3],currAnnotation->nameAnnotationSet)) {
        break;
      }
      i++;
    }
    if (i < arrayMax (currMatrix->annotations)) {
      puts (mio_writeMatrix (currMatrix,0));
    }
  }
  mio_deInit ();
  return 0;
}
