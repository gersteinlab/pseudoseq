#include "log.h"
#include "format.h"
#include "numUtil.h"
#include "intervalFind.h"
#include "mio.h"
#include <gsl/gsl_statistics.h>
#include <math.h>



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  int i,j,k;
  Array rowIndicesInterval1,rowIndicesInterval2;
  Interval *intervalPtr1,*intervalPtr2;
  char *name1 = NULL;
  char *name2 = NULL;
  int pairNumber1,pairNumber2;
  int intervalNumber1,intervalNumber2;
  PairCorrelation *currPairCorrelation;
  double correlation;
  double *values1,*values2;
  int count;

  if (argc != 2) {
    usage ("%s <samples.txt>",argv[0]);
  }
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    if (arrayMax (currMatrix->pairCorrelations) != 0) {
      puts (mio_writeMatrix (currMatrix,0));
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
      rowIndicesInterval1 = mio_getRowIndicesForInterval (currMatrix,intervalNumber1);
      rowIndicesInterval2 = mio_getRowIndicesForInterval (currMatrix,intervalNumber2);
      if (arrayMax (rowIndicesInterval1) != arrayMax (rowIndicesInterval2)) {
        warn ("Unequal length: %d %d",arrayMax (rowIndicesInterval1),arrayMax (rowIndicesInterval2));
      }
      currPairCorrelation = arrayp (currMatrix->pairCorrelations,arrayMax (currMatrix->pairCorrelations),PairCorrelation);
      currPairCorrelation->pairNumber = pairNumber1;
      currPairCorrelation->averageCorrelation = 0.0;
      currPairCorrelation->correlations = arrayCreate (mio_getNumSamples (),double);
      count = 0;
      for (j = 0; j < mio_getNumSamples (); j++) {
        values1 = (double*)malloc (arrayMax (rowIndicesInterval1)*sizeof(double));
        values2 = (double*)malloc (arrayMax (rowIndicesInterval1)*sizeof(double));
        for (k = 0; k < arrayMax (rowIndicesInterval1); k++) {
          values1[k] = mio_getValue (currMatrix,arru (rowIndicesInterval1,k,int),j);
          values2[k] = mio_getValue (currMatrix,arru (rowIndicesInterval2,k,int),j);
        }
        correlation = gsl_stats_correlation (values1,1,values2,1,arrayMax (rowIndicesInterval1));
        if (!isnan (correlation)) {
          currPairCorrelation->averageCorrelation += correlation;
          count++;
        }
        array (currPairCorrelation->correlations,arrayMax (currPairCorrelation->correlations),double) = correlation;
        free (values1);
        free (values2);
      }
      currPairCorrelation->averageCorrelation = currPairCorrelation->averageCorrelation / count; 
      arrayDestroy (rowIndicesInterval1);
      arrayDestroy (rowIndicesInterval2);
    }
    puts (mio_writeMatrix (currMatrix,0));
  }
  mio_deInit ();
  return 0;
}
