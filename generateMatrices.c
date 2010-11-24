#include "log.h"
#include "format.h"
#include "common.h"
#include "linestream.h"
#include "intervalFind.h"
#include "mio.h"



typedef struct {
  char* name;
  int length;
  double percentIdentity;
} AlignmentInfo;



typedef struct {
  char *chromosome;
  int position;
  double value;
  int sampleNumber;
  int intervalNumber;
  int exonNumber;
  int intervalGroupNumber;
} Entry;



typedef struct {
  char *chromosome;
  int start;
  int end;
  double value;
} BedGraph;



typedef struct {
  char *name;
  int intervalIndex;
  int intervalNumber;
  int intervalGroupNumber;
} IntervalMap;



typedef struct {
  char *name;
  int intervalGroupNumber;
  Array intervalIndices;
} SuperIntervalMap;



static int sortSuperIntervalMapsByIntervalGroupNumber (SuperIntervalMap *a, SuperIntervalMap *b) 
{
  return a->intervalGroupNumber - b->intervalGroupNumber;
}



static int sortIntervalMapsByIntervalIndex (IntervalMap *a, IntervalMap *b) 
{
  return a->intervalIndex - b->intervalIndex;
}



static int sortIntervalMapsByIntervalGroupNumber (IntervalMap *a, IntervalMap *b) 
{
  int diff;

  diff = a->intervalGroupNumber - b->intervalGroupNumber;
  if (diff != 0) {
    return diff;
  }
  return a->intervalIndex - b->intervalIndex;
}



static void addIntervalMap (Array intervalMaps, char *name, int intervalIndex, int intervalNumber, int intervalGroupNumber)
{
  IntervalMap *currIM;

  currIM = arrayp (intervalMaps,arrayMax (intervalMaps),IntervalMap);
  currIM->name = hlr_strdup (name);
  currIM->intervalIndex = intervalIndex;
  currIM->intervalNumber = intervalNumber;
  currIM->intervalGroupNumber = intervalGroupNumber;
}



static Array generateIntervalMap (Array intervals)
{
  Array intervalMaps;
  Interval *currInterval,*nextInterval;
  int i,j;
  static char *name1 = NULL;
  static char *name2 = NULL;
  int pairNumber1,pairNumber2,intervalNumber1,intervalNumber2;
  int intervalGroupNumber;
  
  intervalMaps = arrayCreate (100,IntervalMap);
  intervalGroupNumber = 0;
  i = 0; 
  while (i < arrayMax (intervals)) {
    currInterval = arrp (intervals,i,Interval);
    mio_getIntervalInfo (currInterval->name,&name1,&pairNumber1,&intervalNumber1);
    addIntervalMap (intervalMaps,name1,i,intervalNumber1,intervalGroupNumber);
    j = i + 1;
    while (j < arrayMax (intervals)) {
      nextInterval = arrp (intervals,j,Interval);
      mio_getIntervalInfo (nextInterval->name,&name2,&pairNumber2,&intervalNumber2);
      if (!strEqual (name1,name2)) {
        intervalGroupNumber++;
        break;
      }
      addIntervalMap (intervalMaps,name2,j,intervalNumber2,intervalGroupNumber);
      j++;
    }
    i = j;
  }
  return intervalMaps;
}



static Array readBedGraphFile (char *fileName)
{
  LineStream ls;
  char *line;
  Array bedGraphs;
  BedGraph *currBedGraph;
  WordIter w;
  
  bedGraphs = arrayCreate (1000000,BedGraph);
  ls = ls_createFromFile (fileName);
  ls_nextLine (ls); 
  while (line = ls_nextLine (ls)) {
    if (strStartsWithC (line,"track")) {
      continue;
    }
    currBedGraph = arrayp (bedGraphs,arrayMax (bedGraphs),BedGraph);
    w = wordIterCreate (line," \t",1);
    currBedGraph->chromosome = hlr_strdup (wordNext (w));
    currBedGraph->start = atoi (wordNext (w));
    currBedGraph->end = atoi (wordNext (w));
    currBedGraph->value = atof (wordNext (w));
    wordIterDestroy (w);
  }
  ls_destroy (ls);
  return bedGraphs;
}
                            


static int sortBedGraphs (BedGraph *a, BedGraph *b)
{
  int diff;

  diff = strcmp (a->chromosome,b->chromosome);
  if (diff != 0) {
    return diff;
  }
  diff = a->start - b->start;
  if (diff != 0) {
    return diff;
  }
  return b->end - a->end;
}



static void freeBedGraphs (Array bedGraphs)
{
  BedGraph *currBedGraph;
  int i;
  
  for (i = 0; i < arrayMax (bedGraphs); i++) {
    currBedGraph = arrp (bedGraphs,i,BedGraph);
    hlr_free (currBedGraph->chromosome);
  }
  arrayDestroy (bedGraphs);
}



static void addEntry (Array entries, char *chromosome, int position,
                      int intervalGroupNumber, int sampleNumber, int intervalNumber, 
                      int exonNumber, double value) 
{
  Entry *currEntry;
 
  currEntry = arrayp (entries,arrayMax (entries),Entry);
  currEntry->chromosome = hlr_strdup (chromosome);
  currEntry->position = position;
  currEntry->sampleNumber = sampleNumber;
  currEntry->intervalNumber = intervalNumber;
  currEntry->exonNumber = exonNumber;
  currEntry->intervalGroupNumber = intervalGroupNumber;
  currEntry->value = value;
}



static void getValuesForRegion (Array entries, Array bedGraphs, char *chromosome, char strand, int start, int end,
                                int intervalGroupNumber, int sampleNumber, int intervalNumber, int exonNumber)
{
  BedGraph testBedGraph;
  int index;
  int i,j;
  BedGraph *currBedGraph;
  static Array bedGraphPtrs = NULL;
  int numOccurances;
  double value;

  if (bedGraphPtrs == NULL) {
    bedGraphPtrs = arrayCreate (100,BedGraph*);
  }
  else {
    arrayClear (bedGraphPtrs);
  }
  testBedGraph.chromosome = hlr_strdup (chromosome);
  testBedGraph.start = start;
  testBedGraph.end = end;
  arrayFind (bedGraphs,&testBedGraph,&index,(ARRAYORDERF)sortBedGraphs); 
  i = index;
  while (i >= 0) {
    currBedGraph = arrp (bedGraphs,i,BedGraph);
    if (!strEqual (chromosome,currBedGraph->chromosome) || currBedGraph->end < start) {
      break;
    }
    array (bedGraphPtrs,arrayMax (bedGraphPtrs),BedGraph*) = currBedGraph;
    i--;
  }
  i = index + 1;
  while (i < arrayMax (bedGraphs)) {
    currBedGraph = arrp (bedGraphs,i,BedGraph);
    if (!strEqual (chromosome,currBedGraph->chromosome) || currBedGraph->start > end) {
      break;
    }
    array (bedGraphPtrs,arrayMax (bedGraphPtrs),BedGraph*) = currBedGraph;
    i++;
  }
  for (i = start; i < end; i++) {
    numOccurances = 0;
    for (j = 0; j < arrayMax (bedGraphPtrs); j++) {
      currBedGraph = arru (bedGraphPtrs,j,BedGraph*);
      if (currBedGraph->start <= i && i < currBedGraph->end) {
        numOccurances++;
        value = currBedGraph->value;
      } 
    }
    if (numOccurances == 1) {
      addEntry (entries,chromosome,i,intervalGroupNumber,sampleNumber,intervalNumber,exonNumber,value); 
    }
    else if (numOccurances == 0) {
      addEntry (entries,chromosome,i,intervalGroupNumber,sampleNumber,intervalNumber,exonNumber,0.0); 
    }
    else {
      die ("Overlapping BEDGRAPH entries!");
    }
  }
  hlr_free (testBedGraph.chromosome);
}



static int sortEntries (Entry *a, Entry *b)
{
  int diff;

  diff = a->intervalGroupNumber - b->intervalGroupNumber;
  if (diff != 0) {
    return diff;
  }
  diff = a->intervalNumber - b->intervalNumber;
  if (diff != 0) {
    return diff;
  }
  diff = a->position - b->position;
  if (diff != 0) {
    return diff;
  }
  return a->sampleNumber - b->sampleNumber;
}



static Array readAlignmentInfoFile (char *fileName)
{
  LineStream ls;
  char *line;
  Array alignmentInfos;
  AlignmentInfo *currAM;
  WordIter w;
  
  alignmentInfos = arrayCreate (10000,AlignmentInfo);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    currAM = arrayp (alignmentInfos,arrayMax (alignmentInfos),AlignmentInfo);
    w = wordIterCreate (line,"\t",1);
    currAM->name = hlr_strdup (wordNext (w));
    currAM->length = atoi (wordNext (w));
    currAM->percentIdentity = atof (wordNext (w));
    wordIterDestroy (w);
  }
  ls_destroy (ls);
  return alignmentInfos;
}



static int sortAlignmentInfos (AlignmentInfo *a, AlignmentInfo *b)
{
  return  strcmp (a->name,b->name);
}



static Texta generatePairNames (Texta names)
{
  static Texta pairNames = NULL;
  char *pos;
  int i;

  textCreateClear (pairNames,100);
  for (i = 0; i < arrayMax (names); i++) {
    pos = strrchr (textItem (names,i),'|');
    if (pos == NULL) {
      die ("Expected '|' in name: %s",textItem (names,i));
    }
    *pos = '\0';
    textAdd (pairNames,textItem (names,i));
  }
  textUniqKeepOrder (pairNames);
  return pairNames;
}



int main (int argc, char *argv[]) 
{
  Array bedGraphs;
  int i,j,k,m;
  Texta samples;
  Stringa buffer;
  Array intervals;
  Interval *currInterval;
  SubInterval *currSubInterval;
  Array entries;
  Entry *currEntry,*nextEntry;
  char *pos;
  Array intervalMaps;
  IntervalMap *currIM,*nextIM,testIM;
  int index;
  SuperIntervalMap *currSIM,testSIM;
  Array superIntervalMaps;
  Array alignmentInfos;
  AlignmentInfo *currAM,testAM;
  Texta names,pairNames;
  
  if (argc != 5) {
    usage ("%s <samples.txt> <coordinates.interval> <file.meta> <pathToBedGraphDirectory>",argv[0]);
  }
  samples = readList (argv[1]);
  intervals = intervalFind_parseFile (argv[2],0);
  intervalMaps = generateIntervalMap (intervals);
  arraySort (intervalMaps,(ARRAYORDERF)sortIntervalMapsByIntervalIndex);
  alignmentInfos = readAlignmentInfoFile (argv[3]);
  arraySort (alignmentInfos,(ARRAYORDERF)sortAlignmentInfos);
  buffer = stringCreate (100);
  entries = arrayCreate (100000000,Entry);
  for (i = 0; i < arrayMax (samples); i++) {
    pos = strchr (textItem (samples,i),'\t');
    if (pos == NULL) {
      die ("Expected tab between prefix and short name: %s",textItem (samples,i));
    }
    *pos = '\0';
    stringPrintf (buffer,"%s/%s",argv[4],textItem (samples,i));
    bedGraphs = readBedGraphFile (string (buffer));
    arraySort (bedGraphs,(ARRAYORDERF)sortBedGraphs); 
    for (j = 0; j < arrayMax (intervals); j++) {
      currInterval = arrp (intervals,j,Interval);
      testIM.intervalIndex = j;
      if (!arrayFind (intervalMaps,&testIM,&index,(ARRAYORDERF)sortIntervalMapsByIntervalIndex)) {
        die ("Expected to find interval map for %d:",j);
      }
      currIM = arrp (intervalMaps,index,IntervalMap);
      for (k = 0; k < arrayMax (currInterval->subIntervals); k++) {
        currSubInterval = arrp (currInterval->subIntervals,k,SubInterval);
        getValuesForRegion (entries,
                            bedGraphs,
                            currInterval->chromosome,
                            currInterval->strand,
                            currSubInterval->start,
                            currSubInterval->end,
                            currIM->intervalGroupNumber,
                            i + 1,
                            currIM->intervalNumber,
                            k + 1);
      }
    }
    freeBedGraphs (bedGraphs);
    warn ("Done processing: %s",textItem (samples,i));
  } 
  arraySort (intervalMaps,(ARRAYORDERF)sortIntervalMapsByIntervalGroupNumber);
  superIntervalMaps = arrayCreate (1000,SuperIntervalMap);
  i = 0;
  while (i < arrayMax (intervalMaps)) {
    currIM = arrp (intervalMaps,i,IntervalMap);
    currSIM = arrayp (superIntervalMaps,arrayMax (superIntervalMaps),SuperIntervalMap);
    currSIM->name = hlr_strdup (currIM->name);
    currSIM->intervalGroupNumber = currIM->intervalGroupNumber;
    currSIM->intervalIndices = arrayCreate (100,int);
    array (currSIM->intervalIndices,arrayMax (currSIM->intervalIndices),int) = currIM->intervalIndex;
    j = i + 1;
    while (j < arrayMax (intervalMaps)) {
      nextIM = arrp (intervalMaps,j,IntervalMap);
      if (currIM->intervalGroupNumber != nextIM->intervalGroupNumber) {
        break;
      }
      array (currSIM->intervalIndices,arrayMax (currSIM->intervalIndices),int) = nextIM->intervalIndex;
      j++;
    }
    i = j;
  }
  arraySort (entries,(ARRAYORDERF)sortEntries);

  names = textCreate (1000);
  i = 0;
  while (i < arrayMax (entries)) {
    currEntry = arrp (entries,i,Entry);
    testSIM.intervalGroupNumber = currEntry->intervalGroupNumber;
    if (!arrayFind (superIntervalMaps,&testSIM,&index,(ARRAYORDERF)sortSuperIntervalMapsByIntervalGroupNumber)) {
      die ("Expected to find super interval map for %d:",currEntry->intervalGroupNumber);
    }
    currSIM = arrp (superIntervalMaps,index,SuperIntervalMap);
    printf ("Name:\t%s\n",currSIM->name);
    textClear (names);
    for (k = 0; k < arrayMax (currSIM->intervalIndices); k++) {
      currInterval = arrp (intervals,arru (currSIM->intervalIndices,k,int),Interval);
      textAdd (names,currInterval->name);
      printf ("AlignmentInterval:\t%s\t%s\t%c\t%d\t%d\t%d\t",
              currInterval->name,
              currInterval->chromosome,
              currInterval->strand,
              currInterval->start,
              currInterval->end,
              arrayMax (currInterval->subIntervals));
      for (m = 0; m < arrayMax (currInterval->subIntervals); m++) {
        currSubInterval = arrp (currInterval->subIntervals,m,SubInterval);
        printf ("%d%s",currSubInterval->start,m < arrayMax (currInterval->subIntervals) - 1 ? "," : "\t");
      }
      for (m = 0; m < arrayMax (currInterval->subIntervals); m++) {
        currSubInterval = arrp (currInterval->subIntervals,m,SubInterval);
        printf ("%d%s",currSubInterval->end,m < arrayMax (currInterval->subIntervals) - 1 ? "," : "\n");
      }
    }
    pairNames = generatePairNames (names);
    for (k = 0; k < arrayMax (pairNames); k++) {
      testAM.name = hlr_strdup (textItem (pairNames,k));
      if (!arrayFind (alignmentInfos,&testAM,&index,(ARRAYORDERF)sortAlignmentInfos)) {
        die ("Expected to find alignment info: %s",textItem (pairNames,k));
      }
      currAM = arrp (alignmentInfos,index,AlignmentInfo);
      pos = strchr (currAM->name,'|');
      printf ("AlignmentMeta:\t%s\t%d\t%f\n",pos + 1,currAM->length,currAM->percentIdentity);
      hlr_free (testAM.name);
    }
    printf ("%s\t%d\t%d\t%d\t%.4f\t",currEntry->chromosome,currEntry->position,currEntry->intervalNumber,currEntry->exonNumber,currEntry->value);
    j = i + 1;
    while (j < arrayMax (entries)) {
      nextEntry = arrp (entries,j,Entry);
      if (currEntry->intervalGroupNumber != nextEntry->intervalGroupNumber) {
        break;
      }
      if (nextEntry->sampleNumber == 1) {
        printf ("%s\t%d\t%d\t%d\t",nextEntry->chromosome,nextEntry->position,nextEntry->intervalNumber,nextEntry->exonNumber);
      }
      printf ("%.4f",nextEntry->value);
      if (nextEntry->sampleNumber == arrayMax (samples)) {
        printf ("\n");
      }
      else {
        printf ("\t");
      }
      j++;
    }
    i = j;
    puts ("###");
  }
  return 0;
}


