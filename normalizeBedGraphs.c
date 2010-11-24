#include "log.h"
#include "format.h"
#include "linestream.h"
#include <stdio.h>



int main (int argc, char *argv[])
{ 
  Stringa buffer;
  LineStream ls1,ls2;
  char *line1,*line2;
  char *pos;
  char *file = NULL;
  int numReadsPerMillion;
  WordIter w;
  char *chromosome = NULL;
  int start,end;
  double value;
  char *prefix = NULL;
  FILE *fp;
  
  if (argc != 2) {
    usage ("%s <file.txt>",argv[0]);
  }
  buffer = stringCreate (100);
  ls1 = ls_createFromFile (argv[1]);
  while (line1 = ls_nextLine (ls1)) {
    pos = strchr (line1,'\t');
    *pos = '\0';
    if (pos == NULL) {
      die ("Unexpected event",line1);
    }
    numReadsPerMillion = atoi (pos + 1);
    strReplace (&file,line1);
    pos = strstr (line1,".nonNormalized.bgr");
    if (pos == NULL) {
      die ("Expected to find '.nonNormalized.bgr' in file name: %s",line1);
    }
    *pos = '\0';
    strReplace (&prefix,line1);
    stringPrintf (buffer,"%s.bgr",prefix);
    fp = fopen (string (buffer),"w");
    if (fp == NULL) {
      die ("Unable to open file: %s",string (buffer));
    }
    ls2 = ls_createFromFile (file);
    while (line2 = ls_nextLine (ls2)) {
      if (strStartsWithC (line2,"track")) {
        continue;
      }
      w = wordIterCreate (line2,"\t",1);
      strReplace (&chromosome,wordNext (w));
      start = atoi (wordNext (w));
      end = atoi (wordNext (w));
      value = atoi (wordNext (w));
      fprintf (fp,"%s\t%d\t%d\t%.4f\n",chromosome,start,end,value / numReadsPerMillion);
      wordIterDestroy (w);
    }
    ls_destroy (ls2);
    fclose (fp);
  }
  ls_destroy (ls1);
  return 0;
}



