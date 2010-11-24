#include "log.h"
#include "format.h"
#include "mio.h"




int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  Interval *intervalPtr1,*intervalPtr2;
  char *name1 = NULL;
  char *name2 = NULL;
  int pairNumber1,pairNumber2;
  int intervalNumber1,intervalNumber2;
  int length;
  double percentIdentity;
  int i;
  Statistic *currStatistic1,*currStatistic2;

  if (argc != 2) {
    usage ("%s <samples.txt>",argv[0]);
  }
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    if (arrayMax (currMatrix->intervals) != 2) {
      continue;
    }
    for (i = 0; i < arrayMax (currMatrix->intervals); i = i + 2) {
      intervalPtr1 = arrp (currMatrix->intervals,i,Interval);
      intervalPtr2 = arrp (currMatrix->intervals,i + 1,Interval);
      mio_getIntervalInfo (intervalPtr1->name,&name1,&pairNumber1,&intervalNumber1);
      mio_getIntervalInfo (intervalPtr2->name,&name2,&pairNumber2,&intervalNumber2);
      if (pairNumber1 != pairNumber2) {
        die ("Expected same pair!");
      }
      mio_getMetaInfo (currMatrix,pairNumber1,&length,&percentIdentity);
      currStatistic1 = mio_getStatisticPointer (currMatrix,intervalNumber1);
      currStatistic2 = mio_getStatisticPointer (currMatrix,intervalNumber2);
      printf ("%s\t%f\t%f\t%f\t%d\n",currMatrix->name,currStatistic1->overallAverage,currStatistic2->overallAverage,percentIdentity,length);
    }
  }
  mio_deInit ();
  return 0;
}

