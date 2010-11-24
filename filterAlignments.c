#include "log.h"
#include "format.h"
#include "intervalFind.h"
#include "numUtil.h"
#include "common.h"



static int hasAnnotationOverlap (Interval *currInterval)
{
  Array overlapIntervals;
  Interval *thisInterval;
  int i,j,k;
  SubInterval *currIntervalSubInterval,*thisIntervalSubInterval;

  overlapIntervals = intervalFind_getOverlappingIntervals (currInterval->chromosome,currInterval->start,currInterval->end);
  for (i = 0; i < arrayMax (overlapIntervals); i++) {
    thisInterval = arru (overlapIntervals,i,Interval*);
    for (j = 0; j < arrayMax (currInterval->subIntervals); j++) {
      currIntervalSubInterval = arrp (currInterval->subIntervals,j,SubInterval);
      for (k = 0; k < arrayMax (thisInterval->subIntervals); k++) {
        thisIntervalSubInterval = arrp (thisInterval->subIntervals,k,SubInterval);
        if (rangeIntersection (currIntervalSubInterval->start,currIntervalSubInterval->end,
                               thisIntervalSubInterval->start,thisIntervalSubInterval->end) > 0) {
          return 1;
        }
      }
    }
  }
  return 0;
}



int main (int argc, char *argv[])
{ 
  Array alignments;
  Texta alignmentInfos;
  int i;
  Interval *currAlignmentItem1,*currAlignmentItem2;
  Stringa buffer;
  FILE *fp1,*fp2;
  int index;
  int mode;

  if (argc != 4) {
    usage ("%s <prefixAlignmentFiles> <first|second> <annotation.interval>",argv[0]);
  }
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s.interval",argv[1]);
  alignments = intervalFind_parseFile (string (buffer),0);
  stringPrintf (buffer,"%s.meta",argv[1]);
  alignmentInfos = readList (string (buffer));
  if (2 * arrayMax (alignmentInfos) != arrayMax (alignments)) {
    die ("Unexpected counts");
  }
  if (strCaseEqual (argv[2],"first")) {
    mode = 1;
  }
  else if (strCaseEqual (argv[2],"second")) {
    mode = 2;
  }
  else {
    usage ("%s <prefixAlignmentFiles> <first|second> <annotation.interval>",argv[0]);
  }
  intervalFind_addIntervalsToSearchSpace (argv[3],0);
  stringPrintf (buffer,"%s.filtered.interval",argv[1]);
  fp1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s.filtered.meta",argv[1]);
  fp2 = fopen (string (buffer),"w");
  if (fp1 == NULL || fp2 == NULL) {
    die ("Unable to open output files!");
  }
  index = 0;
  for (i = 1; i < arrayMax (alignments); i = i + 2) {
    currAlignmentItem1 = arrp (alignments,i - 1,Interval);
    currAlignmentItem2 = arrp (alignments,i,Interval);
    if ((mode == 1 && hasAnnotationOverlap (currAlignmentItem1)) || (mode == 2 && hasAnnotationOverlap (currAlignmentItem2))) {
      fprintf (fp1,"%s\n",intervalFind_writeInterval (currAlignmentItem1));
      fprintf (fp1,"%s\n",intervalFind_writeInterval (currAlignmentItem2));
      fprintf (fp2,"%s\n",textItem (alignmentInfos,index));
    }
    index++;
  }
  fclose (fp1);
  fclose (fp2);
  return 0;
}
