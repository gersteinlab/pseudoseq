#include "log.h"
#include "format.h"
#include "fasta.h"
#include "stdio.h"



int main (int argc, char *argv[])
{ 
  Seq *currSeq;
  Array sequences;
  int i;
  FILE *fp1,*fp2;
  Stringa buffer;

  if (argc != 2) {
    usage ("%s <fileNamePrefix>",argv[0]);
  }
  seq_init ();
  fasta_initFromFile ("-");
  sequences = fasta_readAllSequences (0);
  fasta_deInit ();
  buffer = stringCreate (100);
  stringPrintf (buffer,"%s.fas",argv[1]);
  fp1 = fopen (string (buffer),"w");
  stringPrintf (buffer,"%s.sequenceIds",argv[1]);
  fp2 = fopen (string (buffer),"w");
  if (fp1 == NULL || fp2 == NULL) {
    die ("Unable to open files");
  }
  for (i = 0; i < arrayMax (sequences); i++) {
    currSeq = arrp (sequences,i,Seq);
    fprintf (fp1,">Sequence_%d\n%s\n",i + 1,currSeq->sequence);
    fprintf (fp2,"Sequence_%d\t%s\n",i + 1,currSeq->name);
  }
  fclose (fp1);
  fclose (fp2);
  stringDestroy (buffer);
  return 0;
}
