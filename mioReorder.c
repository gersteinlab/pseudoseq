#include "log.h"
#include "format.h"
#include "mio.h"
#include "intervalFind.h"



static void reorderInterval (Matrix *currMatrix, Array rowIndices) 
{
  Row thisRow;
  int i;

  for (i = 0; i < arrayMax (rowIndices) / 2; i++) {
    thisRow = arru (currMatrix->rows,arru (rowIndices,i,int),Row);
    arru (currMatrix->rows,arru (rowIndices,i,int),Row) = arru (currMatrix->rows,arru (rowIndices,arrayMax (rowIndices) - 1 - i,int),Row);
    arru (currMatrix->rows,arru (rowIndices,i,int),Row).intervalNumber = thisRow.intervalNumber;
    arru (currMatrix->rows,arru (rowIndices,arrayMax (rowIndices) - 1 - i,int),Row) = thisRow;
    arru (currMatrix->rows,arru (rowIndices,arrayMax (rowIndices) - 1 - i,int),Row).intervalNumber = arru (currMatrix->rows,arru (rowIndices,arrayMax (rowIndices) - 1 - i,int),Row).intervalNumber;
  }
}



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  int i;
  Interval *currInterval,*prevInterval;
  int pairNumber1,pairNumber2;
  int intervalNumber1,intervalNumber2;
  char *name1 = NULL;
  char *name2 = NULL;

  if (argc != 2) {
    usage ("%s <samples.txt>",argv[0]);
  }
  mio_init ("-",argv[1]);
  while (currMatrix = mio_getNextMatrix ()) {
    for (i = 1; i < arrayMax (currMatrix->intervals); i = i + 2) {
      prevInterval = arrp (currMatrix->intervals,i - 1,Interval);
      currInterval = arrp (currMatrix->intervals,i,Interval);
      mio_getIntervalInfo (prevInterval->name,&name1,&pairNumber1,&intervalNumber1);
      mio_getIntervalInfo (currInterval->name,&name2,&pairNumber2,&intervalNumber2);
      if (prevInterval->strand == '+' && currInterval->strand == '+') {
        ; // do nothing
      }
      else if (prevInterval->strand == '+' && currInterval->strand == '-') {
        reorderInterval (currMatrix,mio_getRowIndicesForInterval (currMatrix,intervalNumber2));
      }
      else if (prevInterval->strand == '-' && currInterval->strand == '+') {
        reorderInterval (currMatrix,mio_getRowIndicesForInterval (currMatrix,intervalNumber1));
        reorderInterval (currMatrix,mio_getRowIndicesForInterval (currMatrix,intervalNumber2));
      }
      else if (prevInterval->strand == '-' && currInterval->strand == '-') {
        reorderInterval (currMatrix,mio_getRowIndicesForInterval (currMatrix,intervalNumber1));
      }
      else {
        die ("Unexpected outcome!");
      }
    }
    puts (mio_writeMatrix (currMatrix,0));
  }
  mio_deInit ();
  return 0;
}
