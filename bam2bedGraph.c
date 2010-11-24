#include "log.h"
#include "format.h"
#include "linestream.h"



int main (int argc, char *argv[])
{ 
  Stringa buffer;
  LineStream ls;
  char *line;

  if (argc != 3) {
    usage ("%s <file.bam> <outputPrefix>",argv[0]);
  }
  buffer = stringCreate (100);
  stringPrintf (buffer,"samtools view %s | sam2mrf | cut -f1 > %s.mrf",argv[1],argv[2]);
  warn ("Running: %s",string (buffer));
  hlr_system (string (buffer),1);
  stringPrintf (buffer,"mrfSubsetByTargetName %s < %s.mrf",argv[2],argv[2]);
  warn ("Running: %s",string (buffer));
  hlr_system (string (buffer),1);
  stringPrintf (buffer,"ls -1 %s_chr*.mrf",argv[2]);
  warn ("Running: %s",string (buffer));
  ls = ls_createFromPipe (string (buffer));
  while (line = ls_nextLine (ls)) {
    stringPrintf (buffer,"mrf2bedGraph %s doNotNormalize < %s",argv[2],line);
    warn ("Running: %s",string (buffer));
    hlr_system (string (buffer),1);
    stringPrintf (buffer,"rm -rf %s",line);
    warn ("Running: %s",string (buffer));
    hlr_system (string (buffer),1);
  }
  ls_destroy (ls);
  stringPrintf (buffer,"cat %s_chr*.bgr > %s.nonNormalized.bgr",argv[2],argv[2]);
  warn ("Running: %s",string (buffer));
  hlr_system (string (buffer),1);
  stringPrintf (buffer,"rm -rf %s_chr*.bgr",argv[2]);
  warn ("Running: %s",string (buffer));
  hlr_system (string (buffer),1);
  stringDestroy (buffer);
  return 0;
}


