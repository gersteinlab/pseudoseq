#include "log.h"
#include "format.h"
#include "mio.h"



typedef struct {
  int sampleNumber;
  int intervalNumber;
  double value;
} Entry;



static int sortEntries (Entry *a, Entry *b)
{
  int diff;

  diff = a->intervalNumber - b->intervalNumber;
  if (diff != 0) {
    return diff;
  }
  return a->sampleNumber - b->sampleNumber;
}



static void calculateStatistics (Array statistics, Array entryPointers, int intervalNumber) 
{
  int i,j;
  Entry *currEntry,*nextEntry;
  double sum,min,max;
  Statistic *currStatistic;
  int count;

  currStatistic = arrayp (statistics,arrayMax (statistics),Statistic);
  currStatistic->intervalNumber= intervalNumber;
  currStatistic->averages = arrayCreate (mio_getNumSamples (),double);
  currStatistic->mins = arrayCreate (mio_getNumSamples (),double);
  currStatistic->maxs = arrayCreate (mio_getNumSamples (),double);
  i = 0;
  while (i < arrayMax (entryPointers)) {
    currEntry = arru (entryPointers,i,Entry*);
    min = currEntry->value;
    max = currEntry->value;
    sum = currEntry->value;
    count = 1;
    j = i + 1;
    while (j < arrayMax (entryPointers)) {
      nextEntry = arru (entryPointers,j,Entry*);
      if (currEntry->sampleNumber != nextEntry->sampleNumber) {
        break;
      }
      count++;
      sum += nextEntry->value;
      if (nextEntry->value < min) {
        min = nextEntry->value;
      } 
      if (nextEntry->value > max) {
        max = nextEntry->value;
      } 
      j++;
    }
    array (currStatistic->averages,arrayMax (currStatistic->averages),double) = sum / count;
    array (currStatistic->mins,arrayMax (currStatistic->mins),double) = min;
    array (currStatistic->maxs,arrayMax (currStatistic->maxs),double) = max; 
    i = j;
  }   
  currStatistic->overallAverage = arru (currStatistic->averages,0,double);
  currStatistic->overallMax = arru (currStatistic->maxs,0,double);
  currStatistic->overallMin = arru (currStatistic->mins,0,double);
  for (i = 1; i < mio_getNumSamples (); i++) {
    currStatistic->overallAverage += arru (currStatistic->averages,i,double);
    if (arru (currStatistic->maxs,i,double) > currStatistic->overallMax) {
      currStatistic->overallMax = arru (currStatistic->maxs,i,double);
    }
    if (arru (currStatistic->mins,i,double) < currStatistic->overallMin) {
      currStatistic->overallMin = arru (currStatistic->mins,i,double);
    }
  }
  currStatistic->overallAverage = currStatistic->overallAverage / mio_getNumSamples ();
}



static void getIntervalStatistics (Matrix *currMatrix)
{
  int i,j;
  Row *currRow;
  Array entries;
  Entry *currEntry,*nextEntry;
  static Array entryPointers = NULL;

  entries = arrayCreate (1000,Entry);
  for (i = 0; i < arrayMax (currMatrix->rows); i++) {
    currRow =  arrp (currMatrix->rows,i,Row);
    for (j = 0; j < mio_getNumSamples (); j++) {
      currEntry = arrayp (entries,arrayMax (entries),Entry);
      currEntry->sampleNumber = j;
      currEntry->intervalNumber = currRow->intervalNumber;
      currEntry->value = arru (currRow->values,j,double);
    }
  }
  arraySort (entries,(ARRAYORDERF)sortEntries);

  if (entryPointers == NULL) {
    entryPointers = arrayCreate (1000,Entry*);
  }
  i = 0;
  while (i < arrayMax (entries)) { 
    arrayClear (entryPointers);
    currEntry = arrp (entries,i,Entry);
    array (entryPointers,arrayMax (entryPointers),Entry*) = currEntry;
    j = i + 1;
    while (j < arrayMax (entries)) {
      nextEntry = arrp (entries,j,Entry);
      if (currEntry->intervalNumber != nextEntry->intervalNumber) {
        break;
      }
      array (entryPointers,arrayMax (entryPointers),Entry*) = nextEntry;
      j++;
    }
    calculateStatistics (currMatrix->statistics,entryPointers,currEntry->intervalNumber);  
    i = j;
  }
  arrayDestroy (entries);
}



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  
  if (argc != 2) {
    usage ("%s <samples.txt>",argv[0]);
  }
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    if (arrayMax (currMatrix->statistics) == 0) {
      getIntervalStatistics (currMatrix);
    }
    puts (mio_writeMatrix (currMatrix,0));
  }
  mio_deInit ();
  return 0;
}

