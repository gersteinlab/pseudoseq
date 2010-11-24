#include "log.h"
#include "format.h"
#include "intervalFind.h"
#include "numUtil.h"
#include "common.h"
#include "mio.h"
  


static int printOutput (FILE *fp1, FILE *fp2, Array alignmentPtrs, Array alignmentInfoPtrs, int *intervalLineCount, int maxNumIntervalsPerEntry)
{
  int i;

  if (arrayMax (alignmentPtrs) != (2 * arrayMax (alignmentInfoPtrs))) {
    die ("Unexpected counts");
  }
  if (arrayMax (alignmentPtrs) > maxNumIntervalsPerEntry) {
    return 0;
  }
  for (i = 0; i < arrayMax (alignmentPtrs); i++) {
    fprintf (fp1,"%s\n",intervalFind_writeInterval (arru (alignmentPtrs,i,Interval*)));
  }
  *intervalLineCount += arrayMax (alignmentPtrs);
  for (i = 0; i < arrayMax (alignmentInfoPtrs); i++) {
    fprintf (fp2,"%s\n",*arru (alignmentInfoPtrs,i,char**));
  }  
  return 1;
}



int main (int argc, char *argv[])
{ 
  Array alignments;
  Texta alignmentInfos;
  int i,j;
  Interval *currAlignment,*nextAlignment;;
  Stringa buffer;
  FILE *fp1 = NULL;
  FILE *fp2 = NULL;
  char *name1 = NULL;
  char *name2 = NULL;
  int pairNumber1,pairNumber2,intervalNumber1,intervalNumber2;
  int subsetCount;
  int maxNumIntervalsPerFile;
  int intervalLineCount;
  Array alignmentPtrs;
  Array alignmentInfoPtrs;
  int maxNumIntervalsPerEntry;
  int writeSuccess;

  if (argc != 4) {
    usage ("%s <prefixAlignmentFiles> <maxNumIntervalsPerFile> <maxNumIntervalsPerEntry>",argv[0]);
  }
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s.interval",argv[1]);
  alignments = intervalFind_parseFile (string (buffer),0);
  stringPrintf (buffer,"%s.meta",argv[1]);
  alignmentInfos = readList (string (buffer));
  if (2 * arrayMax (alignmentInfos) != arrayMax (alignments)) {
    die ("Unexpected counts");
  }
  maxNumIntervalsPerFile = atoi (argv[2]);
  maxNumIntervalsPerEntry = atoi (argv[3]);
  subsetCount = 0;
  intervalLineCount = 0;
  writeSuccess = 1;
  alignmentPtrs = arrayCreate (1000,Interval*);
  alignmentInfoPtrs = arrayCreate (1000,char**);
  i = 0; 
  while (i < arrayMax (alignments)) {
    arrayClear (alignmentPtrs);
    arrayClear (alignmentInfoPtrs);
    currAlignment = arrp (alignments,i,Interval);
    mio_getIntervalInfo (currAlignment->name,&name1,&pairNumber1,&intervalNumber1);
    if (intervalLineCount == 0 && writeSuccess == 1) {
      if (fp1 != NULL) {
        fclose (fp1);
      }
      if (fp2 != NULL) {
        fclose (fp2);
      }
      subsetCount++;
      stringPrintf (buffer,"%s.subset%04d.interval",argv[1],subsetCount);
      fp1 = fopen (string (buffer),"w");
      stringPrintf (buffer,"%s.subset%04d.meta",argv[1],subsetCount);
      fp2 = fopen (string (buffer),"w");
      if (fp1 == NULL || fp2 == NULL) {
        die ("Unable to open output files");
      }
    }
    array (alignmentPtrs,arrayMax (alignmentPtrs),Interval*) = currAlignment;
    j = i + 1;
    while (j < arrayMax (alignments)) {
      nextAlignment = arrp (alignments,j,Interval);
      mio_getIntervalInfo (nextAlignment->name,&name2,&pairNumber2,&intervalNumber2);
      if (!strEqual (name1,name2)) {
        break;
      }
      if ((j % 2) == 1) {
        array (alignmentInfoPtrs,arrayMax (alignmentInfoPtrs),char**) = &textItem (alignmentInfos,j / 2);
      }
      array (alignmentPtrs,arrayMax (alignmentPtrs),Interval*) = nextAlignment;
      j++;
    }
    i = j;
    writeSuccess = printOutput (fp1,fp2,alignmentPtrs,alignmentInfoPtrs,&intervalLineCount,maxNumIntervalsPerEntry);
    if (intervalLineCount > maxNumIntervalsPerFile) {
      intervalLineCount = 0;
    }
  }
  fclose (fp1);
  fclose (fp2);
  stringDestroy (buffer);
  arrayDestroy (alignmentPtrs);
  arrayDestroy (alignmentInfoPtrs);
  return 0;
}
