#include "log.h"
#include "format.h"
#include "mio.h"
#include <stdio.h>



static void writeBedFile (FILE *fp, Interval *currInterval)
{
  int i;
  SubInterval *currSubInterval;
  
  for (i = 0; i < arrayMax (currInterval->subIntervals); i++) {
    currSubInterval = arrp (currInterval->subIntervals,i,SubInterval);
    fprintf (fp,"%s\t%d\t%d\n",currInterval->chromosome,currSubInterval->start,currSubInterval->end);
  }
}



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  int i;
  Interval *intervalPtr1,*intervalPtr2;
  char *name1 = NULL;
  char *name2 = NULL;
  int pairNumber1,pairNumber2;
  int intervalNumber1,intervalNumber2;
  FILE *fp;
  Stringa buffer;

  if (argc != 2) {
    usage ("%s <samples.txt>",argv[0]);
  }
  buffer = stringCreate (100);
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    for (i = 0; i < arrayMax (currMatrix->intervals); i = i + 2) {
      intervalPtr1 = arrp (currMatrix->intervals,i,Interval);
      intervalPtr2 = arrp (currMatrix->intervals,i + 1,Interval);
      mio_getIntervalInfo (intervalPtr1->name,&name1,&pairNumber1,&intervalNumber1);
      mio_getIntervalInfo (intervalPtr2->name,&name2,&pairNumber2,&intervalNumber2);
      if (pairNumber1 != pairNumber2) {
        die ("Expected same pair!");
      }
      stringPrintf (buffer,"%s_%d_%d.bed",currMatrix->name,pairNumber1,intervalNumber1);
      fp = fopen (string (buffer),"w");
      if (fp == NULL) {
        die ("Unable to open file: %s",string (buffer));
      }
      writeBedFile (fp,intervalPtr1);
      fclose (fp);
      stringPrintf (buffer,"%s_%d_%d.bed",currMatrix->name,pairNumber2,intervalNumber2);
      fp = fopen (string (buffer),"w");
      if (fp == NULL) {
        die ("Unable to open file: %s",string (buffer));
      }
      writeBedFile (fp,intervalPtr2);
      fclose (fp);
    }
  }
  mio_deInit ();
  stringDestroy (buffer);
  return 0;
}
