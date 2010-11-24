#include "log.h"
#include "format.h"
#include "linestream.h"
#include "intervalFind.h"
#include "common.h"
#include "stringUtil.h"
#include "mio.h"



static LineStream ls = NULL;
static Texta sampleNames = NULL;
static Texta samplePrefixes = NULL;
static int numSamples;



void mio_init (char *fileNameMatrices, char *fileNameSamples)
{
  int i;
  char *pos;
  Texta samples;

  ls = ls_createFromFile (fileNameMatrices);;
  samples = readList (fileNameSamples);
  numSamples = arrayMax (samples);
  sampleNames = textCreate (arrayMax (samples));
  samplePrefixes = textCreate (arrayMax (samples));
  for (i = 0; i < arrayMax (samples); i++) {
    pos = strchr (textItem (samples,i),'\t');
    if (pos == NULL) {
      die ("Expected tab between prefix and name: %s",textItem (samples,i));
    }
    *pos = '\0';
    textAdd (samplePrefixes,textItem (samples,i));
    textAdd (sampleNames,pos + 1);
  }
  textDestroy (samples);
}



void mio_deInit (void)
{
  ls_destroy (ls);
  textDestroy (sampleNames);
  textDestroy (samplePrefixes);
}



int mio_getNumSamples (void)
{
  return numSamples;
}



Texta mio_getSamplesNames (void)
{
  return sampleNames;
}



static void mio_freeMatrix (Matrix *currMatrix)
{
  int i;
  Row *currRow;
  Column *currColumn;
  Statistic *currStatistic;
  Interval *currInterval;
  PairCorrelation *currPairCorrelation;
  Annotation *currAnnotation;

  if (currMatrix == NULL) {
    return;
  }
  hlr_free (currMatrix->name);
  for (i = 0; i < arrayMax (currMatrix->rows); i++) {
    currRow = arrp (currMatrix->rows,i,Row);
    hlr_free (currRow->chromosome);
    arrayDestroy (currRow->values);
  }
  arrayDestroy (currMatrix->rows);
  for (i = 0; i < arrayMax (currMatrix->columns); i++) {
    currColumn = arrp (currMatrix->columns,i,Column);
    hlr_free (currColumn->samplePrefix);
    hlr_free (currColumn->sampleName);
  }
  arrayDestroy (currMatrix->columns);
  for (i = 0; i < arrayMax (currMatrix->statistics); i++) {
    currStatistic = arrp (currMatrix->statistics,i,Statistic);
    arrayDestroy (currStatistic->averages);
    arrayDestroy (currStatistic->mins);
    arrayDestroy (currStatistic->maxs);
  }
  arrayDestroy (currMatrix->statistics);
  arrayDestroy (currMatrix->scores);
  for (i = 0; i < arrayMax (currMatrix->intervals); i++) {
    currInterval = arrp (currMatrix->intervals,i,Interval);
    hlr_free (currInterval->name);
    hlr_free (currInterval->chromosome);
    arrayDestroy (currInterval->subIntervals);
  }
  arrayDestroy (currMatrix->intervals);
  arrayDestroy (currMatrix->alignmentMetas);
  for (i = 0; i < arrayMax (currMatrix->pairCorrelations); i++) {
    currPairCorrelation = arrp (currMatrix->pairCorrelations,i,PairCorrelation);
    arrayDestroy (currPairCorrelation->correlations);
  }
  arrayDestroy (currMatrix->pairCorrelations);
  for (i = 0; i < arrayMax (currMatrix->annotations); i++) {
    currAnnotation = arrp (currMatrix->annotations,i,Annotation);
    hlr_free (currAnnotation->nameAnnotationSet);
    hlr_free (currAnnotation->annotation->name);
    hlr_free (currAnnotation->annotation->chromosome);
    arrayDestroy (currAnnotation->annotation->subIntervals);
    freeMem (currAnnotation->annotation);
  }
  arrayDestroy (currMatrix->annotations);
  freeMem (currMatrix);
  currMatrix = NULL;
}



static void mio_processNextMatrix (Matrix **thisMatrix)
{
  char *line;
  Row *currRow;
  Column *currColumn;
  WordIter w,w2;
  char *token;
  int i;
  char *pos;
  Interval *currInterval;
  int proceed;
  AlignmentMeta *currAlignmentMeta;
  PairCorrelation *currPairCorrelation;
  Statistic *currStatistic;
  Annotation *currAnnotation;
  static Stringa buffer = NULL;
  Matrix *currMatrix;
  Score *currScore;

  currMatrix = *thisMatrix;
  proceed = 0;
  while (line = ls_nextLine (ls)) {
    if (strEqual (line,"###")) {
      break;
    }
    if (strStartsWithC (line,"Name:")) {
      proceed = 1;
      pos = strchr (line,'\t');
      if (pos == NULL) {
        die ("Unexpected line: %s",line);
      }
      currMatrix->name = hlr_strdup (pos + 1);
      currMatrix->rows = arrayCreate (1000,Row);
      currMatrix->columns = arrayCreate (mio_getNumSamples (),Column);
      currMatrix->intervals = arrayCreate (10,Interval);
      currMatrix->statistics = arrayCreate (10,Statistic);
      currMatrix->alignmentMetas = arrayCreate (10,AlignmentMeta);
      currMatrix->pairCorrelations = arrayCreate (10,PairCorrelation);
      currMatrix->annotations = arrayCreate (10,Annotation);
      currMatrix->scores = arrayCreate (10,Score);
      continue;
    }
    if (strStartsWithC (line,"AlignmentInterval:")) {
      pos = strchr (line,'\t');
      if (pos == NULL) {
        die ("Unexpected line: %s",line);
      }
      currInterval = arrayp (currMatrix->intervals,arrayMax (currMatrix->intervals),Interval);
      intervalFind_parseLine (currInterval,pos + 1,0);
      continue;
    }
    if (strStartsWithC (line,"AlignmentMeta:")) {
      currAlignmentMeta = arrayp (currMatrix->alignmentMetas,arrayMax (currMatrix->alignmentMetas),AlignmentMeta);
      w = wordIterCreate (line,"\t",1);
      wordNext (w); // get rid off prefix
      currAlignmentMeta->pairNumber = atoi (wordNext (w));
      currAlignmentMeta->length = atoi (wordNext (w));
      currAlignmentMeta->percentIdentity = atof (wordNext (w));
      wordIterDestroy (w);
      continue;
    }
    if (strStartsWithC (line,"PairCorrelation:")) {
      currPairCorrelation = arrayp (currMatrix->pairCorrelations,arrayMax (currMatrix->pairCorrelations),PairCorrelation);
      w = wordIterCreate (line,"\t",1);
      wordNext (w); // get rid off prefix
      currPairCorrelation->pairNumber = atoi (wordNext (w));
      currPairCorrelation->averageCorrelation = atof (wordNext (w));
      currPairCorrelation->correlations = arrayCreate (mio_getNumSamples (),double);
      w2 = wordIterCreate (wordNext (w),",",1);
      while (token = wordNext (w2)) {
        array (currPairCorrelation->correlations,arrayMax (currPairCorrelation->correlations),double) = atof (token);
      }
      wordIterDestroy (w2);
      wordIterDestroy (w);
      continue;
    }
    if (strStartsWithC (line,"Statistic:")) {
      currStatistic = arrayp (currMatrix->statistics,arrayMax (currMatrix->statistics),Statistic);
      w = wordIterCreate (line,"\t",1);
      wordNext (w); // get rid off prefix
      currStatistic->intervalNumber = atoi (wordNext (w));
      currStatistic->overallAverage = atof (wordNext (w));
      currStatistic->overallMin = atof (wordNext (w));
      currStatistic->overallMax = atof (wordNext (w));
      currStatistic->averages = arrayCreate (mio_getNumSamples (),double);
      currStatistic->mins = arrayCreate (mio_getNumSamples (),double);
      currStatistic->maxs = arrayCreate (mio_getNumSamples (),double);
      w2 = wordIterCreate (wordNext (w),",",1);
      while (token = wordNext (w2)) {
        array (currStatistic->averages,arrayMax (currStatistic->averages),double) = atof (token);
      }
      wordIterDestroy (w2);
      w2 = wordIterCreate (wordNext (w),",",1);
      while (token = wordNext (w2)) {
        array (currStatistic->mins,arrayMax (currStatistic->mins),double) = atof (token);
      }
      wordIterDestroy (w2);
      w2 = wordIterCreate (wordNext (w),",",1);
      while (token = wordNext (w2)) {
        array (currStatistic->maxs,arrayMax (currStatistic->maxs),double) = atof (token);
      }
      wordIterDestroy (w2);
      wordIterDestroy (w);
      continue;
    }
    if (strStartsWithC (line,"Score:")) {
      currScore = arrayp (currMatrix->scores,arrayMax (currMatrix->scores),Score);
      w = wordIterCreate (line,"\t",1);
      wordNext (w); // get rid off prefix
      currScore->pairNumber = atoi (wordNext (w));
      currScore->type = atoi (wordNext (w));
      currScore->score = atof (wordNext (w));
      wordIterDestroy (w);
      continue;
    }
    if (strStartsWithC (line,"Annotation:")) {
      currAnnotation = arrayp (currMatrix->annotations,arrayMax (currMatrix->annotations),Annotation);
      w = wordIterCreate (line,"\t",1);
      wordNext (w); // get rid off prefix
      currAnnotation->nameAnnotationSet = hlr_strdup (wordNext (w));
      currAnnotation->intervalNumber = atoi (wordNext (w));
      currAnnotation->percentIntervalOverlap = atof (wordNext (w));
      currAnnotation->percentAnnotationOverlap = atof (wordNext (w));
      AllocVar (currAnnotation->annotation);
      stringCreateClear (buffer,100);
      while (token = wordNext (w)) {
        if (stringLen (buffer) > 0) {
          stringCat (buffer,"\t");
        }
        stringCat (buffer,token);
      }
      intervalFind_parseLine (currAnnotation->annotation,string (buffer),0);
      wordIterDestroy (w);
      continue;
    }
    currRow = arrayp (currMatrix->rows,arrayMax (currMatrix->rows),Row);
    w = wordIterCreate (line,"\t",0);
    currRow->chromosome = hlr_strdup (wordNext (w));
    currRow->position = atoi (wordNext (w));
    currRow->intervalNumber = atoi (wordNext (w));
    currRow->exonNumber = atoi (wordNext (w));
    currRow->values = arrayCreate (100,double);
    while (token = wordNext (w)) {
      array (currRow->values,arrayMax (currRow->values),double) = atof (token);
    }
    wordIterDestroy (w);
  }
  if (proceed == 0) {
    *thisMatrix = NULL;
    return;
  }
  for (i = 0; i < numSamples; i++) {
    currColumn =  arrayp (currMatrix->columns,arrayMax (currMatrix->columns),Column);
    currColumn->sampleNumber = i + 1;
    currColumn->samplePrefix = hlr_strdup (textItem (samplePrefixes,i));
    currColumn->sampleName = hlr_strdup (textItem (sampleNames,i));
  }
}



Matrix* mio_getNextMatrix (void)
{
  static Matrix *currMatrix = NULL; 
 
  mio_freeMatrix (currMatrix);
  AllocVar (currMatrix);  
  mio_processNextMatrix (&currMatrix); 
  return currMatrix;
}



Array mio_parse (void)
{
  Array matrices;
  Matrix *currMatrix;
  int count;
  
  count = 0;
  matrices = arrayCreate (100000,Matrix);
  while (1) {
    currMatrix = arrayp (matrices,arrayMax (matrices),Matrix);
    mio_processNextMatrix (&currMatrix); 
    if (currMatrix == NULL) {
      arraySetMax (matrices,count);
      break;
    }
    count++;
  }
  return matrices;
}



double mio_getValue (Matrix *currMatrix, int rowIndex, int columnIndex)
{
  Row *currRow;
  
  currRow = arrp (currMatrix->rows,rowIndex,Row);
  return arru (currRow->values,columnIndex,double);
}



char* mio_writeMatrix (Matrix *currMatrix, int writeMetaInformationOnly)
{
  Row *currRow;
  int i,j;
  Interval *currInterval;
  AlignmentMeta *currAlignmentMeta;
  Statistic *currStatistic;
  PairCorrelation *currPairCorrelation;
  Annotation *currAnnotation;
  Score *currScore;
  static Stringa buffer = NULL;

  stringCreateClear (buffer,100000);
  stringPrintf (buffer,"Name:\t%s\n",currMatrix->name);
  for (i = 0; i < arrayMax (currMatrix->intervals); i++) {
    currInterval = arrp (currMatrix->intervals,i,Interval);
    stringAppendf (buffer,"AlignmentInterval:\t");
    stringAppendf (buffer,"%s\n",intervalFind_writeInterval (currInterval));
  }
  for (i = 0; i < arrayMax (currMatrix->alignmentMetas); i++) {
    currAlignmentMeta = arrp (currMatrix->alignmentMetas,i,AlignmentMeta);
    stringAppendf (buffer,"AlignmentMeta:\t%d\t%d\t%f\n",currAlignmentMeta->pairNumber,currAlignmentMeta->length,currAlignmentMeta->percentIdentity);
  }
  for (i = 0; i < arrayMax (currMatrix->pairCorrelations); i++) {
    currPairCorrelation = arrp (currMatrix->pairCorrelations,i,PairCorrelation);
    stringAppendf (buffer,"PairCorrelation:\t%d\t%f\t",currPairCorrelation->pairNumber,currPairCorrelation->averageCorrelation);
    for (j = 0; j < arrayMax (currPairCorrelation->correlations); j++) {
      stringAppendf (buffer,"%f%s",arru (currPairCorrelation->correlations,j,double),j < arrayMax (currPairCorrelation->correlations) - 1 ? "," : "\n");
    }
  }
  for (i = 0; i < arrayMax (currMatrix->statistics); i++) {
    currStatistic = arrp (currMatrix->statistics,i,Statistic);
    stringAppendf (buffer,"Statistic:\t%d\t%f\t%f\t%f\t",currStatistic->intervalNumber,currStatistic->overallAverage,currStatistic->overallMin,currStatistic->overallMax);
    for (j = 0; j < arrayMax (currStatistic->averages); j++) {
      stringAppendf (buffer,"%f%s",arru (currStatistic->averages,j,double),j < arrayMax (currStatistic->averages) - 1 ? "," : "\t");
    }
    for (j = 0; j < arrayMax (currStatistic->mins); j++) {
      stringAppendf (buffer,"%f%s",arru (currStatistic->mins,j,double),j < arrayMax (currStatistic->mins) - 1 ? "," : "\t");
    }
    for (j = 0; j < arrayMax (currStatistic->maxs); j++) {
      stringAppendf (buffer,"%f%s",arru (currStatistic->maxs,j,double),j < arrayMax (currStatistic->maxs) - 1 ? "," : "\n");
    }
  }
  for (i = 0; i < arrayMax (currMatrix->scores); i++) {
    currScore = arrp (currMatrix->scores,i,Score);
    stringAppendf (buffer,"Score:\t%d\t%d\t%f\n",currScore->pairNumber,currScore->type,currScore->score);
  }
  for (i = 0; i < arrayMax (currMatrix->annotations); i++) {
    currAnnotation = arrp (currMatrix->annotations,i,Annotation);
    stringAppendf (buffer,"Annotation:\t%s\t%d\t%.1f\t%.1f\t",currAnnotation->nameAnnotationSet,currAnnotation->intervalNumber,currAnnotation->percentIntervalOverlap,currAnnotation->percentAnnotationOverlap);
    stringAppendf (buffer,"%s\n",intervalFind_writeInterval (currAnnotation->annotation));
  }
  if (writeMetaInformationOnly == 0) {
    for (i = 0; i < arrayMax (currMatrix->rows); i++) {
      currRow = arrp (currMatrix->rows,i,Row);
      stringAppendf (buffer,"%s\t%d\t%d\t%d\t",currRow->chromosome,currRow->position,currRow->intervalNumber,currRow->exonNumber);
      for (j = 0; j < arrayMax (currRow->values); j++) {
        stringAppendf (buffer,"%.4f%s",arru (currRow->values,j,double),j < arrayMax (currRow->values) - 1 ? "\t" : "\n");
      }
    }
  }
  stringAppendf (buffer,"%s","###");
  return string (buffer);
}



Array mio_getRowIndicesForInterval (Matrix *currMatrix, int intervalNumber)
{
  Row *currRow;
  int i;
  Array rowIndices;
 
  rowIndices = arrayCreate (10000,int);
  for (i = 0; i < arrayMax (currMatrix->rows); i++) {
    currRow = arrp (currMatrix->rows,i,Row);
    if (currRow->intervalNumber == intervalNumber) {
      array (rowIndices,arrayMax (rowIndices),int) = i;
    }
  }
  return rowIndices;
}



void mio_getIntervalInfo (char *intervalName, char** name, int *pairNumber, int *intervalNumber)
{
  Texta tokens;

  tokens = textFieldtokP (intervalName,"|");
  if (arrayMax (tokens) != 3) {
    die ("Unexpected interval name format: %s",intervalName);
  }
  strReplace (name,textItem (tokens,0));
  *pairNumber = atoi (textItem (tokens,1));
  *intervalNumber = atoi (textItem (tokens,2));
  textDestroy (tokens);
}



void mio_getMetaInfo (Matrix *currMatrix, int pairNumber, int *length, double *percentIdentity)
{
  int i;
  AlignmentMeta *currAlignmentMeta;

  for (i = 0; i < arrayMax (currMatrix->alignmentMetas); i++) {
    currAlignmentMeta = arrp (currMatrix->alignmentMetas,i,AlignmentMeta);
    if (currAlignmentMeta->pairNumber == pairNumber) {
      *length = currAlignmentMeta->length;
      *percentIdentity = currAlignmentMeta->percentIdentity;
      return;
    }
  }
}



PairCorrelation* mio_getPairCorrelationPointer (Matrix *currMatrix, int pairNumber) 
{
  int i;
  PairCorrelation *currPairCorrelation;

  for (i = 0; i < arrayMax (currMatrix->pairCorrelations); i++) {
    currPairCorrelation = arrp (currMatrix->pairCorrelations,i,PairCorrelation);
    if (currPairCorrelation->pairNumber == pairNumber) {
      return currPairCorrelation;
    }
  }
  return NULL;
}



Statistic* mio_getStatisticPointer (Matrix *currMatrix, int intervalNumber)
{
  int i;
  Statistic *currStatistic;

  for (i = 0; i < arrayMax (currMatrix->statistics); i++) {
    currStatistic = arrp (currMatrix->statistics,i,Statistic);
    if (currStatistic->intervalNumber == intervalNumber) {
      return currStatistic;
    }
  }
  return NULL;
}



Array mio_getAnnotationPointers (Matrix *currMatrix, int intervalNumber)
{
  int i;
  Annotation *currAnnotation;
  Array annotationPointers;

  annotationPointers = arrayCreate (100,Annotation*);
  for (i = 0; i < arrayMax (currMatrix->annotations); i++) {
    currAnnotation = arrp (currMatrix->annotations,i,Annotation);
    if (currAnnotation->intervalNumber == intervalNumber) {
      array (annotationPointers,arrayMax (annotationPointers),Annotation*) = currAnnotation;
    }
  }
  return annotationPointers;
}



Score* mio_getScorePointer (Matrix *currMatrix, int pairNumber)
{
  int i;
  Score *currScore;

  for (i = 0; i < arrayMax (currMatrix->scores); i++) {
    currScore = arrp (currMatrix->scores,i,Score);
    if (currScore->pairNumber == pairNumber) {
      return currScore;
    }
  }
  return NULL;
}
