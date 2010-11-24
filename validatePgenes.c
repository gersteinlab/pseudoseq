#include "log.h"
#include "format.h"
#include "numUtil.h"
#include "intervalFind.h"



static int overlapIntervals (Interval *currPgene, Interval *currCds) 
{
  SubInterval *currPgeneSubInterval,*currCdsSubInterval;
  int k,l;
  int overlap;

  overlap = 0;
  for (k = 0; k < arrayMax (currPgene->subIntervals); k++) {
    currPgeneSubInterval = arrp (currPgene->subIntervals,k,SubInterval);
    for (l = 0; l < arrayMax (currCds->subIntervals); l++) {
      currCdsSubInterval = arrp (currCds->subIntervals,l,SubInterval);
      overlap += positiveRangeIntersection (currPgeneSubInterval->start,currPgeneSubInterval->end,currCdsSubInterval->start,currCdsSubInterval->end);  
    }
  }
  if (overlap > 0) {
    warn ("%d\t%s\t%s",overlap,currPgene->name,currCds->name);
    return 1;
  }
  return 0;
}



int main (int argc, char *argv[])
{ 
  int i,j;
  Array pgenes,cdss;
  Interval *currPgene,*currCds;
  int numOverlaps;
  
  if (argc != 3) {
    usage ("%s <cdsAnnotation.interval> <pgeneAnnotation.interval>",argv[0]);
  }
   
  intervalFind_addIntervalsToSearchSpace (argv[1],0);
  pgenes = intervalFind_parseFile (argv[2],0);
  for (i = 0; i < arrayMax (pgenes); i++) {
    currPgene = arrp (pgenes,i,Interval);
    cdss = intervalFind_getOverlappingIntervals (currPgene->chromosome,currPgene->start,currPgene->end);
    numOverlaps = 0;
    for (j = 0; j < arrayMax (cdss); j++) {
      currCds = arru (cdss,j,Interval*);
      numOverlaps += overlapIntervals (currPgene,currCds);
    }
    if (numOverlaps == 0) {
      puts (intervalFind_writeInterval (currPgene));
    }
  }
  return 0;
}

