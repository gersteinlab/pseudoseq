#include "log.h"
#include "format.h"
#include "numUtil.h"
#include "intervalFind.h"
#include "common.h"
#include "mio.h"




static int getSize (Interval *currInterval) 
{
  SubInterval *currSubInterval;
  int i;
  int sum;

  sum = 0;
  for (i = 0; i < arrayMax (currInterval->subIntervals); i++) {
    currSubInterval = arrp (currInterval->subIntervals,i,SubInterval);
    sum += (currSubInterval->end - currSubInterval->start);
  }
  return sum;
}



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  Interval *currInterval;
  Array overlapIntervals;
  Interval *thisInterval;
  int i,j,k,l;
  char *sampleName = NULL;
  int pairNumber,intervalNumber;
  int overlap;
  SubInterval *currIntervalSubInterval,*thisIntervalSubInterval;
  Annotation *currAnnotation;

  if (argc != 4) {
    usage ("%s <samples.txt> <annotation.interval> <nameAnnotationSet>",argv[0]);
  }
  mio_init ("-",argv[1]);
  intervalFind_addIntervalsToSearchSpace (argv[2],0);
  while (currMatrix = mio_getNextMatrix ()) {
    for (i = 0; i < arrayMax (currMatrix->intervals); i++) {
      currInterval = arrp (currMatrix->intervals,i,Interval);
      mio_getIntervalInfo (currInterval->name,&sampleName,&pairNumber,&intervalNumber);
      overlapIntervals = intervalFind_getOverlappingIntervals (currInterval->chromosome,currInterval->start,currInterval->end);
      for (j = 0; j < arrayMax (overlapIntervals); j++) {
        thisInterval = arru (overlapIntervals,j,Interval*);
        overlap = 0;
        for (k = 0; k < arrayMax (currInterval->subIntervals); k++) {
          currIntervalSubInterval = arrp (currInterval->subIntervals,k,SubInterval);
          for (l = 0; l < arrayMax (thisInterval->subIntervals); l++) {
            thisIntervalSubInterval = arrp (thisInterval->subIntervals,l,SubInterval);
            overlap += positiveRangeIntersection (currIntervalSubInterval->start,currIntervalSubInterval->end,
                                                  thisIntervalSubInterval->start,thisIntervalSubInterval->end);
          }
        }
        if (overlap > 0) {
          currAnnotation = arrayp (currMatrix->annotations,arrayMax (currMatrix->annotations),Annotation);
          currAnnotation->nameAnnotationSet = hlr_strdup (argv[3]);
          AllocVar (currAnnotation->annotation);
          currAnnotation->annotation->name = hlr_strdup (thisInterval->name);
          currAnnotation->annotation->chromosome = hlr_strdup (thisInterval->chromosome);
          currAnnotation->annotation->source = thisInterval->source;
          currAnnotation->annotation->strand = thisInterval->strand;
          currAnnotation->annotation->start = thisInterval->start;
          currAnnotation->annotation->end = thisInterval->end;
          currAnnotation->annotation->subIntervalCount = thisInterval->subIntervalCount;
          currAnnotation->annotation->subIntervals = arrayCopy (thisInterval->subIntervals);
          currAnnotation->intervalNumber = intervalNumber;
          currAnnotation->percentIntervalOverlap = 1.0 * overlap / getSize (currInterval) * 100;
          currAnnotation->percentAnnotationOverlap = 1.0 * overlap / getSize (thisInterval) * 100;
        }
      }
    }
    puts (mio_writeMatrix (currMatrix,0));
  }
  mio_deInit ();
  return 0;
}
