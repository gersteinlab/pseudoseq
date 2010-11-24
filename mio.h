#ifndef DEF_MIO_H
#define DEF_MIO_H



#include "intervalFind.h"



#define SCORE_TYPE_ALL 0
#define SCORE_TYPE_TRANSCRIPTION 1
#define SCORE_TYPE_DISCORDANT 2
#define SCORE_TYPE_CONCORDANT 3



typedef struct {
  char *name;
  Array intervals;
  Array rows;
  Array columns;
  Array statistics;
  Array alignmentMetas;
  Array pairCorrelations;
  Array annotations;
  Array scores;
} Matrix;



typedef struct {
  char *chromosome;
  int position;
  int intervalNumber;
  int exonNumber;
  Array values;
} Row;



typedef struct {
  int sampleNumber;
  char *samplePrefix;
  char *sampleName;
} Column;



typedef struct {
  int intervalNumber;
  double overallAverage;
  double overallMin;
  double overallMax;
  Array averages;
  Array mins;
  Array maxs;
} Statistic;



typedef struct {
  int pairNumber;
  int length;
  double percentIdentity;
} AlignmentMeta;



typedef struct {
  int pairNumber;
  double averageCorrelation;
  Array correlations;
} PairCorrelation;



typedef struct {
  char *nameAnnotationSet;
  int intervalNumber;
  Interval *annotation;
  double percentIntervalOverlap;
  double percentAnnotationOverlap;
} Annotation;



typedef struct {
  int pairNumber;
  int type;
  double score;
} Score;



extern void mio_init (char *fileNameMatrices, char *fileNameSamples);
extern void mio_deInit (void);
extern int mio_getNumSamples (void);
extern Texta mio_getSamplesNames (void);
extern Matrix* mio_getNextMatrix (void);
extern Array mio_parse (void);
extern char* mio_writeMatrix (Matrix *currMatrix, int writeMetaInformationOnly);
extern double mio_getValue (Matrix *currMatrix, int rowIndex, int columnIndex);
extern Array mio_getRowIndicesForInterval (Matrix *currMatrix, int intervalNumber);
extern void mio_getIntervalInfo (char *intervalName, char** sampleName, int *pairNumber, int *intervalNumber);
extern void mio_getMetaInfo (Matrix *currMatrix, int pairNumber, int *length, double *percentIdentity);
extern PairCorrelation* mio_getPairCorrelationPointer (Matrix *currMatrix, int pairNumber);
extern Statistic* mio_getStatisticPointer (Matrix *currMatrix, int intervalNumber);
extern Array mio_getAnnotationPointers (Matrix *currMatrix, int intervalNumber);
extern Score* mio_getScorePointer (Matrix *currMatrix, int pairNumber);



#endif
