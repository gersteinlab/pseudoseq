#include "log.h"
#include "format.h"
#include "mio.h"



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  int i;
  Interval *intervalPtr1,*intervalPtr2;
  char *name1 = NULL;
  char *name2 = NULL;
  int pairNumber1,pairNumber2;
  int intervalNumber1,intervalNumber2;
  Score *currScore;
  Statistic *currStatistic1,*currStatistic2;
  PairCorrelation *currPairCorrelation;

  if (argc != 2) {
    usage ("%s <samples.txt>",argv[0]);
  }
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    if (arrayMax (currMatrix->statistics) == 0) {
      die ("No statistics found! Run 'mioAddStatistics'");
    }
    if (arrayMax (currMatrix->pairCorrelations) == 0) {
      die ("No correlations found! Run 'mioAddCorrelations'");
    }
    for (i = 0; i < arrayMax (currMatrix->intervals); i = i + 2) {
      intervalPtr1 = arrp (currMatrix->intervals,i,Interval);
      intervalPtr2 = arrp (currMatrix->intervals,i + 1,Interval);
      mio_getIntervalInfo (intervalPtr1->name,&name1,&pairNumber1,&intervalNumber1);
      mio_getIntervalInfo (intervalPtr2->name,&name2,&pairNumber2,&intervalNumber2);
      if (pairNumber1 != pairNumber2) {
        die ("Expected same pair!");
      }
      currStatistic1 = mio_getStatisticPointer (currMatrix,intervalNumber1);
      currStatistic2 = mio_getStatisticPointer (currMatrix,intervalNumber2);
      currPairCorrelation = mio_getPairCorrelationPointer (currMatrix,pairNumber1);
      currScore = arrayp (currMatrix->scores,arrayMax (currMatrix->scores),Score);
      currScore->pairNumber = pairNumber1;
      if (currStatistic1->overallAverage > currStatistic2->overallAverage) {
        currScore->type = SCORE_TYPE_TRANSCRIPTION;
      }
      else if (currPairCorrelation->averageCorrelation < 0) {
        currScore->type = SCORE_TYPE_DISCORDANT;
      }
      else {
        currScore->type = SCORE_TYPE_CONCORDANT;
      }   
      currScore->score = -1.0 * currStatistic1->overallAverage * currPairCorrelation->averageCorrelation * currStatistic1->overallAverage / currStatistic2->overallAverage;
    }
    puts (mio_writeMatrix (currMatrix,0));
  }
  mio_deInit ();
  return 0;
}

