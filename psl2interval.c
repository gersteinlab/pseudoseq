#include "log.h"
#include "format.h"
#include "numUtil.h"
#include "blatParser.h"
#include "linestream.h"



typedef struct {
  int queryOffset;
  int genomicOffset;
} Offset;



typedef struct {
  int start;
  int end;
} QueryExon;



typedef struct {
  char *id;
  char *name;
} SequenceId;



static int sortOffsets (Offset *a, Offset *b)
{
  return a->queryOffset - b->queryOffset;
}



static int queryOffset2genomicOffset (Array offsets, int queryOffset)
{
  int index;
  Offset testOffset;

  testOffset.queryOffset = queryOffset;
  if (!arrayFind (offsets,&testOffset,&index,(ARRAYORDERF)sortOffsets)) {
    die ("Expected to find queryOffset: %d",queryOffset);
  }
  return arrp (offsets,index,Offset)->genomicOffset;
}



static void extractQueryInformation (Array offsets, Array queryExons, char *sequenceName, char **name, char **chromosome, char *strand)
{
  Texta tokens;
  int i,j;
  Offset *currOffset;
  QueryExon *currQueryExon;
  int index;
  static Stringa buffer = NULL;

  stringCreateClear (buffer,100);
  tokens = textFieldtokP (sequenceName,"|");
  stringPrintf (buffer,"%s_%s",textItem (tokens,0),textItem (tokens,1));
  strReplace (name,string (buffer));
  strReplace (chromosome,textItem (tokens,4));
  *strand = textItem (tokens,5)[0];
  index = 0;
  for (i = 6; i < arrayMax (tokens); i = i + 2) {
    currQueryExon = arrayp (queryExons,arrayMax (queryExons),QueryExon);
    currQueryExon->start = atoi (textItem (tokens,i));
    currQueryExon->end = atoi (textItem (tokens,i + 1));
    for (j = currQueryExon->start; j < currQueryExon->end; j++) {
      currOffset = arrayp (offsets,arrayMax (offsets),Offset);
      currOffset->queryOffset = index;
      currOffset->genomicOffset = j;
      index++;
    }
  }
  textDestroy (tokens);
}



static void printInterval (FILE *fp, char *name, int pairNumber, int pairId, char *chromosome, char strand, Array starts, Array ends)
{
  int i;

  fprintf (fp,"%s|%d|%d\t",name,pairNumber,pairId);
  fprintf (fp,"%s\t",chromosome);
  fprintf (fp,"%c\t",strand);
  fprintf (fp,"%d\t",arru (starts,0,int));
  fprintf (fp,"%d\t",arru (ends,arrayMax (ends) - 1,int));
  fprintf (fp,"%d\t",arrayMax (starts));
  for (i = 0; i < arrayMax (starts); i++) {
    fprintf (fp,"%d%s",arru (starts,i,int),i < arrayMax (starts) - 1 ? "," : "\t");
  }
  for (i = 0; i < arrayMax (ends); i++) {
    fprintf (fp,"%d%s",arru (ends,i,int),i < arrayMax (ends) - 1 ? "," : "\n");
  }
}



static void printQueryInterval (FILE *fp, Array queryExons, char *name, int pairNumber, int pairId, char *chromosome, char strand, Array starts, Array ends)
{
  int i,j;
  Array newStarts;
  Array newEnds;
  QueryExon *currQueryExon;
  int overlap;
  
  if (arrayMax (starts) !=  arrayMax (ends)) {
    die ("Expected same number of starts and ends");
  }
  newStarts = arrayCreate (100,int);
  newEnds = arrayCreate (100,int);
  for (i = 0; i < arrayMax (starts); i++) {
    for (j = 0; j < arrayMax (queryExons); j++) {
      currQueryExon = arrp (queryExons,j,QueryExon);
      overlap = rangeIntersection (arru (starts,i,int),arru (ends,i,int),currQueryExon->start,currQueryExon->end);
      if (overlap > 0) {
        if (currQueryExon->start <= arru (starts,i,int) && arru (ends,i,int) <= currQueryExon->end) {
          array (newStarts,arrayMax (newStarts),int) = arru (starts,i,int);
          array (newEnds,arrayMax (newEnds),int) = arru (ends,i,int);
        }
        else if (arru (starts,i,int) < currQueryExon->start && currQueryExon->end < arru (ends,i,int)) {
          array (newStarts,arrayMax (newStarts),int) = currQueryExon->start;
          array (newEnds,arrayMax (newEnds),int) = currQueryExon->end;
        }
        else if (arru (ends,i,int) > currQueryExon->end) {
          array (newStarts,arrayMax (newStarts),int) = arru (starts,i,int);
          array (newEnds,arrayMax (newEnds),int) = currQueryExon->end;
        }
        else if (arru (starts,i,int) < currQueryExon->start) {
          array (newStarts,arrayMax (newStarts),int) = currQueryExon->start;
          array (newEnds,arrayMax (newEnds),int) = arru (ends,i,int);
        }
        else {
          die ("Unexpected outcome");
        }
      }
    }
  }
  if (arrayMax (newStarts) !=  arrayMax (newEnds)) {
    die ("Expected same number of starts and ends");
  }
  printInterval (fp,name,pairNumber,pairId,chromosome,strand,newStarts,newEnds);
  arrayDestroy (newStarts);
  arrayDestroy (newEnds);
}



static void printTargetInterval (FILE *fp, char *name, int pairNumber, int pairId, char *chromosome, char strand, Array starts, Array ends)
{
  if (arrayMax (starts) !=  arrayMax (ends)) {
    die ("Expected same number of starts and ends");
  }
  printInterval (fp,name,pairNumber,pairId,chromosome,strand,starts,ends);
}



static int targetSameAsQuery (Array queryExons, char* queryChromosome, PslEntry *currPE)
{
  int i;
  QueryExon *currQueryExon;

  for (i = 0; i < arrayMax (queryExons); i++) {
   currQueryExon = arrp (queryExons,i,QueryExon);
   if (rangeIntersection (currPE->tStart,currPE->tEnd,currQueryExon->start,currQueryExon->end) > 0) {
     return 1;
   }
  }
  return 0;
}



static Array readSequenceIds (char *fileName)
{
  LineStream ls;
  char *line;
  char *pos;
  SequenceId *currSI;
  Array sequenceIds;
  
  sequenceIds = arrayCreate (10000,SequenceId);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    pos = strchr (line,'\t');
    if (pos == NULL) {
      die ("Unexpected line: %s",line);
    }
    *pos = '\0';
    currSI = arrayp (sequenceIds,arrayMax (sequenceIds),SequenceId);
    currSI->id = hlr_strdup (line);
    currSI->name = hlr_strdup (pos + 1);
  }
  ls_destroy (ls);
  return sequenceIds;
}



static int sortSequenceIds (SequenceId *a, SequenceId *b)
{
  return strcmp (a->id,b->id);
}



static int sortPslEntries (PslEntry **a, PslEntry **b) 
{
  double aPercentIdentity,bPercentIdentity;

  aPercentIdentity = (double)(*a)->matches / ((*a)->matches + (*a)->misMatches);
  bPercentIdentity = (double)(*b)->matches / ((*b)->matches + (*b)->misMatches);
  if (aPercentIdentity < bPercentIdentity) {
    return 1;
  }
  if (aPercentIdentity > bPercentIdentity) {
    return -1;
  }
  return 0;
}



int main (int argc, char *argv[])
{ 
  BlatQuery *currBQ;
  PslEntry *currPE;
  int i,j;
  int minBlockSize;
  char *name = NULL;
  char *chromosome = NULL;
  char strand;
  Array offsets;
  Array queryExons;
  Array starts,ends;
  Array validHits;
  Array sequenceIds;
  SequenceId *currSI,testId;
  int index;
  FILE *fp1,*fp2;
  Stringa buffer;
  
  if (argc != 3) {
    usage ("%s <prefixAlignmentFiles> <minBlockSize>",argv[0]);
  }  
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s.psl",argv[1]);
  blatParser_initFromFile (string (buffer));
  stringPrintf (buffer,"%s.sequenceIds",argv[1]);
  sequenceIds = readSequenceIds (string (buffer));
  arraySort (sequenceIds,(ARRAYORDERF)sortSequenceIds);
  minBlockSize = atoi (argv[2]);
  offsets = arrayCreate (1000,Offset);
  queryExons = arrayCreate (100,QueryExon);
  starts = arrayCreate (100,int);
  ends = arrayCreate (100,int); 
  validHits = arrayCreate (100,PslEntry*);
  stringPrintf (buffer,"%s.minBlock%s.interval",argv[1],argv[2]);
  fp1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s.minBlock%s.meta",argv[1],argv[2]);
  fp2 = fopen (string (buffer),"w");
  if (fp1 == NULL || fp2 == NULL) {
    die ("Unable to open output files");
  }
  while (currBQ = blatParser_nextQuery ()) {
    arrayClear (offsets);
    arrayClear (queryExons);
    arrayClear (validHits);
    testId.id = hlr_strdup (currBQ->qName);
    if (!arrayFind (sequenceIds,&testId,&index,(ARRAYORDERF)sortSequenceIds)) {
      die ("Expected to find sequenceId: %s",testId.id);
    }
    hlr_free (testId.id);
    currSI = arrp (sequenceIds,index,SequenceId);
    extractQueryInformation (offsets,queryExons,currSI->name,&name,&chromosome,&strand);
    for (i = 0; i < arrayMax (currBQ->entries); i++) {
      currPE = arrp (currBQ->entries,i,PslEntry);
      if (targetSameAsQuery (queryExons,chromosome,currPE) == 1) {
        continue;
      }
      j = 0; 
      while (j < arrayMax (currPE->blockSizes)) {
        if (arru (currPE->blockSizes,j,int) > minBlockSize) {
          break;
        }
        j++;
      }
      if (j == arrayMax (currPE->blockSizes)) {
        continue;
      }
      array (validHits,arrayMax (validHits),PslEntry*) = currPE; 
    }
    arraySort (validHits,(ARRAYORDERF)sortPslEntries);
    for (i = 0; i < arrayMax (validHits); i++) {
      currPE = arru (validHits,i,PslEntry*);  
      if (currPE->strand == '+') {
        arrayClear (starts);
        arrayClear (ends);
        for (j = 0; j < arrayMax (currPE->blockSizes); j++) {
          array (starts,arrayMax (starts),int) = queryOffset2genomicOffset (offsets,arru (currPE->qStarts,j,int));
        }
        for (j = 0; j < arrayMax (currPE->blockSizes); j++) {
          array (ends,arrayMax (ends),int) = queryOffset2genomicOffset (offsets,arru (currPE->qStarts,j,int) + arru (currPE->blockSizes,j,int) - 1) + 1;
        }
        printQueryInterval (fp1,queryExons,name,i + 1,2 * i + 1,chromosome,strand,starts,ends);
        arrayClear (starts);
        arrayClear (ends);
        for (j = 0; j < arrayMax (currPE->blockSizes); j++) {
          array (starts,arrayMax (starts),int) = arru (currPE->tStarts,j,int);
        }
        for (j = 0; j < arrayMax (currPE->blockSizes); j++) {
          array (ends,arrayMax (ends),int) = arru (currPE->tStarts,j,int) + arru (currPE->blockSizes,j,int);
        }
        printTargetInterval (fp1,name,i + 1,2 * i + 2,currPE->tName,currPE->strand,starts,ends);
      }
      else if (currPE->strand == '-') {
        arrayClear (starts);
        arrayClear (ends);
        for (j = arrayMax (currPE->blockSizes) - 1; j >= 0; j--) {
          array (starts,arrayMax (starts),int) = queryOffset2genomicOffset (offsets,currPE->qSize - arru (currPE->qStarts,j,int) - arru (currPE->blockSizes,j,int));
        }
        for (j = arrayMax (currPE->blockSizes) - 1; j >= 0; j--) {
          array (ends,arrayMax (ends),int) = queryOffset2genomicOffset (offsets,currPE->qSize - arru (currPE->qStarts,j,int) - 1) + 1;
        }
        printQueryInterval (fp1,queryExons,name,i + 1,2 * i + 1,chromosome,strand,starts,ends);
        arrayClear (starts);
        arrayClear (ends);
        for (j = 0; j < arrayMax (currPE->blockSizes); j++) {
          array (starts,arrayMax (starts),int) = arru (currPE->tStarts,j,int);
        }
        for (j = 0; j < arrayMax (currPE->blockSizes); j++) {
          array (ends,arrayMax (ends),int) = arru (currPE->tStarts,j,int) + arru (currPE->blockSizes,j,int);
        }
        printTargetInterval (fp1,name,i + 1,2 * i + 2,currPE->tName,currPE->strand,starts,ends);
      }
      fprintf (fp2,"%s|%d\t%d\t%f\n",name,i + 1,currPE->matches + currPE->misMatches,(double)currPE->matches / (currPE->matches + currPE->misMatches));
    }
  }
  blatParser_deInit ();
  fclose (fp1);
  fclose (fp2);
  stringDestroy (buffer);
  arrayDestroy (offsets);
  arrayDestroy (queryExons);
  arrayDestroy (starts);
  arrayDestroy (ends);
  arrayDestroy (validHits);
  return 0;
}

 

