#include "log.h"
#include "format.h"
#include "mio.h"



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  int i;
  Statistic *currStatistic;
  double minAverageIntervalExpressionLevel;
  int mod;

  if (argc != 4) {
    usage ("%s <samples.txt> <first|second> <minAverageIntervalExpressionLevel>",argv[0]);
  }
  mio_init ("-",argv[1]);
  if (strCaseEqual (argv[2],"first")) {
    mod = 1;
  }
  else if (strCaseEqual (argv[2],"second")) {
    mod = 0;
  }
  else {
    usage ("%s <samples.txt> <first|second> <minAverageIntervalExpressionLevel>",argv[0]);
  }
  minAverageIntervalExpressionLevel = atof (argv[3]);
  while (currMatrix = mio_getNextMatrix ()) {
    i = 0; 
    while (i < arrayMax (currMatrix->statistics)) {
      currStatistic = arrp (currMatrix->statistics,i,Statistic);
      if ((currStatistic->intervalNumber % 2) == mod &&
          currStatistic->overallAverage > minAverageIntervalExpressionLevel) {
        break;
      }
      i++;
    }
    if (i < arrayMax (currMatrix->statistics)) {
      puts (mio_writeMatrix (currMatrix,0));
    }
  }
  mio_deInit ();
  return 0;
}
