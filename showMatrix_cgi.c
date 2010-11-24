#include "log.h"
#include "format.h"
#include "html.h"
#include "htmlLinker.h"
#include "common.h"
#include "mio.h"
#include <math.h>



#define WEB_URL_CGI "http://dynamic.gersteinlab.org/people/lh372"
#define WEB_DATA_DIR "/home/lh372/public_html/PseudoSeq"
#define WEB_DATA_URL "http://homes.gersteinlab.org/people/lh372/PseudoSeq"
#define WEB_DATA_SET_ENCODE "ENCODE"
#define WEB_DATA_SET_HBM "HumanBodyMap"
#define WEB_DATA_META_FILE_ENCODE "ENCODE.minBlock100.filtered.all.final.minExpression0.1.score.meta"
#define WEB_DATA_META_FILE_HBM "HumanBodyMap.minBlock100.filtered.all.final.minExpression0.1.score.meta"
#define WEB_DATA_SAMPLE_FILE_ENCODE "samples.txt"
#define WEB_DATA_SAMPLE_FILE_HBM "samples.txt"



static void printJavaScript (int numPairs)
{
  int i;

  puts ("<meta charset=\"utf-8\">");
  puts ("<link rel=\"stylesheet\" href=\"http://jqueryui.com/themes/base/jquery.ui.all.css\">");
  puts ("<Script src=\"http://code.jquery.com/jquery-1.4.3.js\"></script>");
  puts ("<script src=\"http://jqueryui.com/Ui/jquery.ui.core.js\"></script>");
  puts ("<script src=\"http://jqueryui.com/ui/jquery.ui.widget.js\"></Script>");
  puts ("<script src=\"http://jqueryui.com/ui/jquery.ui.mouse.js\"></script>");
  puts ("<script src=\"http://jqueryui.com/ui/jquery.ui.tabs.js\"></script>");
  puts ("<script>");
  puts ("");
  for (i = 0; i < numPairs; i++) {
    puts ("$(function() {");
    printf ("	$( \"#tabs-%d\" ).tabs()\n",i + 1);
    puts ("});");
  }
  puts ("</script>");
}



static void processData (char *dataSet, char *matrixName, int pairNumber)
{
  Matrix *currMatrix;
  static Stringa buffer1 = NULL;
  static Stringa buffer2 = NULL;
  char *sampleFile = NULL;
  int i,j;
  Interval *intervalPtr1,*intervalPtr2;
  char *name1 = NULL;
  char *name2 = NULL;
  int pairNumber1,pairNumber2;
  int intervalNumber1,intervalNumber2;
  int length;
  double percentIdentity;
  Array annotationPointers1,annotationPointers2;
  Statistic *currStatistic1,*currStatistic2;
  PairCorrelation *currPairCorrelation;
  Score *currScore;
  int index;
  Texta sampleNames;
  Annotation *currAnnotation;

  stringCreateClear (buffer1,100);
  stringCreateClear (buffer2,100);
  if (strEqual (dataSet,WEB_DATA_SET_ENCODE)) {
    strReplace (&sampleFile,WEB_DATA_SAMPLE_FILE_ENCODE);
  }
  else if (strEqual (dataSet,WEB_DATA_SET_HBM)) {
    strReplace (&sampleFile,WEB_DATA_SAMPLE_FILE_HBM);
  }
  else {
    die ("Unknown data set: %s",dataSet);
  }
  stringPrintf (buffer1,"%s/%s/%s/%s.matrix",WEB_DATA_DIR,dataSet,matrixName,matrixName);
  stringPrintf (buffer2,"%s/%s/%s",WEB_DATA_DIR,dataSet,sampleFile);

  mio_init (string (buffer1),string (buffer2));
  while (currMatrix = mio_getNextMatrix ()) {
    if (arrayMax (currMatrix->statistics) == 0) {
      die ("No statistics found! Run 'mioAddStatistics'");
    }
    if (arrayMax (currMatrix->pairCorrelations) == 0) {
      die ("No correlations found! Run 'mioAddCorrelations'");
    }
    if (arrayMax (currMatrix->annotations) == 0) {
      die ("No annotations found! Run 'mioAddAnnotations'");
    }
    if (arrayMax (currMatrix->scores) == 0) {
      die ("No scores found! Run 'mioAddScores'");
    }
    puts ("<html>");
    puts ("<head>");
    html_printGenericStyleSheet (12);
    puts ("<title>PseudoSeq - Details</title>\n");
    printJavaScript (arrayMax (currMatrix->intervals) / 2);
    puts ("</head>");
    puts ("<body>");
    puts ("<a name=top>&nbsp;</a>");
    printf ("<h1><center>Results for %s (%s)</center></h1><br>\n",matrixName,dataSet);
    puts ("<br><br><center><h2>Overview</h2></center><br><br>");
    puts ("<table border=1 cellpadding=2 width=80% align=center>");
    puts ("<tr align=left>");
    puts ("<th>Pair number</th>");
    puts ("<th>Alignment length</th>");
    puts ("<th>Percent identity</th>");
    puts ("<th>Average pair-wise correlation</th>");
    puts ("<th>Average expression PGENE</th>");
    puts ("<th>Average expression HIT</th>");
    puts ("<th>Score type</th>");
    puts ("<th>Score</th>");
    puts ("</tr>");
    for (i = 0; i < arrayMax (currMatrix->intervals); i = i + 2) {
      intervalPtr1 = arrp (currMatrix->intervals,i,Interval);
      intervalPtr2 = arrp (currMatrix->intervals,i + 1,Interval);
      mio_getIntervalInfo (intervalPtr1->name,&name1,&pairNumber1,&intervalNumber1);
      mio_getIntervalInfo (intervalPtr2->name,&name2,&pairNumber2,&intervalNumber2);
      if (pairNumber1 != pairNumber2) {
        die ("Expected same pair!");
      }
      mio_getMetaInfo (currMatrix,pairNumber1,&length,&percentIdentity);
      currStatistic1 = mio_getStatisticPointer (currMatrix,intervalNumber1);
      currStatistic2 = mio_getStatisticPointer (currMatrix,intervalNumber2);
      annotationPointers1 = mio_getAnnotationPointers (currMatrix,intervalNumber1);
      annotationPointers2 = mio_getAnnotationPointers (currMatrix,intervalNumber2);
      currPairCorrelation = mio_getPairCorrelationPointer (currMatrix,pairNumber1);
      currScore = mio_getScorePointer (currMatrix,pairNumber1);
      printf ("<tr bgcolor=%s>\n",pairNumber == pairNumber1 ? "E8E8E8" : "FFFFFF");
      printf ("<td><a href=#pairNumber%d>%d</a></td>\n",pairNumber1,pairNumber1);
      printf ("<td>%d</td>\n",length);
      printf ("<td>%.2f</td>\n",percentIdentity);
      printf ("<td>%.4f</td>\n",currPairCorrelation->averageCorrelation);
      printf ("<td>%.4f</td>\n",currStatistic1->overallAverage);
      printf ("<td>%.4f</td>\n",currStatistic2->overallAverage);
      printf ("<td>%d</td>\n",currScore->type);
      printf ("<td>%.2f</td>\n",currScore->score);
      puts ("</tr>");
      arrayDestroy (annotationPointers1);
      arrayDestroy (annotationPointers2);
      fflush (stdout);
    }
    puts ("</table>");
    
    index = 1;
    for (i = 0; i < arrayMax (currMatrix->intervals); i = i + 2) {
      intervalPtr1 = arrp (currMatrix->intervals,i,Interval);
      intervalPtr2 = arrp (currMatrix->intervals,i + 1,Interval);
      mio_getIntervalInfo (intervalPtr1->name,&name1,&pairNumber1,&intervalNumber1);
      mio_getIntervalInfo (intervalPtr2->name,&name2,&pairNumber2,&intervalNumber2);
      if (pairNumber1 != pairNumber2) {
        die ("Expected same pair!");
      }
      mio_getMetaInfo (currMatrix,pairNumber1,&length,&percentIdentity);
      currStatistic1 = mio_getStatisticPointer (currMatrix,intervalNumber1);
      currStatistic2 = mio_getStatisticPointer (currMatrix,intervalNumber2);
      annotationPointers1 = mio_getAnnotationPointers (currMatrix,intervalNumber1);
      annotationPointers2 = mio_getAnnotationPointers (currMatrix,intervalNumber2);
      currPairCorrelation = mio_getPairCorrelationPointer (currMatrix,pairNumber1);
      currScore = mio_getScorePointer (currMatrix,pairNumber1);
      puts ("<br><br><br><br>");
      puts ("<hr>");
      puts ("<br><br><br>");
      printf ("<a name=pairNumber%d>&nbsp;</a>\n",pairNumber1);
      puts ("<br><center>[<a href=#top>Top</a>]<center><br>");
      printf ("<center><h2>Summary for alignment pair number %d</h2></center>\n",pairNumber1);
      puts ("<br>");
      puts ("<center><h3>Alignment summary</h3></center>");      
      puts ("<br>");
      printf ("<center><b>Alignment length</b>: %d &nbsp;&nbsp;<b>Percent sequence identity</b>: %.2f<br></center>\n",length,percentIdentity * 100);
      puts ("<br>");
      puts ("<table border=1 cellpadding=2 width=80% align=center>");
      puts ("<tr align=left>");
      puts ("<td><b>Feature</b></td>");
      puts ("<td><b>PGENE</b></td>");
      puts ("<td><b>HIT</b></td>");
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td><b>Number of alignment blocks</b></td>");
      printf ("<td>%d</td>\n",intervalPtr1->subIntervalCount);
      printf ("<td>%d</td>\n",intervalPtr2->subIntervalCount);
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td><b>Genomic coordinates</b></td>");
      printf ("<td>%s:%d-%d</td>\n",intervalPtr1->chromosome,intervalPtr1->start,intervalPtr1->end);
      printf ("<td>%s:%d-%d</td>\n",intervalPtr2->chromosome,intervalPtr2->start,intervalPtr2->end);
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td><b>Strand</b></td>");
      printf ("<td>%c</td>\n",intervalPtr1->strand);
      printf ("<td>%c</td>\n",intervalPtr2->strand);
      puts ("</tr>");
      puts ("</table>");
      puts ("<br>");
      puts ("<br>");
      puts ("<br>");
      puts ("<center><h3>Annotation summary</h3></center>");
      puts ("<br>");
      puts ("<center><h4>PGENE annotations</h4></center>");
      puts ("<table border=1 cellpadding=2 width=80% align=center>");
      puts ("<tr align=left>");
      puts ("<td><b>Annotation set</b></td>");
      puts ("<td><b>Name</b></td>");
      puts ("<td><b>Genomic coordinates</b></td>");
      puts ("<td><b>Strand</b></td>");
      puts ("<td><b>Percent alignment overlap</b></td>");
      puts ("<td><b>Percent annotation overlap</b></td>");
      puts ("</tr>");
      for (j = 0; j < arrayMax (annotationPointers1); j++) {
        currAnnotation = arru (annotationPointers1,j,Annotation*);
        puts ("<tr align=left>");
        printf ("<td>%s</td>\n",currAnnotation->nameAnnotationSet);
        printf ("<td>%s</td>\n",currAnnotation->annotation->name);
        printf ("<td>%s:%d-%d</td>\n",
                currAnnotation->annotation->chromosome,
                currAnnotation->annotation->start,
                currAnnotation->annotation->end);
        printf ("<td>%c</td>\n",currAnnotation->annotation->strand);
        printf ("<td>%.2f</td>\n",currAnnotation->percentIntervalOverlap);
        printf ("<td>%.2f</td>\n",currAnnotation->percentAnnotationOverlap);
        puts ("</tr>");
      }
      puts ("</table>");
      puts ("<br>");
      puts ("<center><h4>HIT annotations</h4></center>");
      puts ("<table border=1 cellpadding=2 width=80% align=center>");
      puts ("<tr align=left>");
      puts ("<td><b>Annotation set</b></td>");
      puts ("<td><b>Name</b></td>");
      puts ("<td><b>Genomic coordinates</b></td>");
      puts ("<td><b>Strand</b></td>");
      puts ("<td><b>Percent alignment overlap</b></td>");
      puts ("<td><b>Percent annotation overlap</b></td>");
      puts ("</tr>");
      for (j = 0; j < arrayMax (annotationPointers2); j++) {
        currAnnotation = arru (annotationPointers2,j,Annotation*);
        puts ("<tr align=left>");
        printf ("<td>%s</td>\n",currAnnotation->nameAnnotationSet);
        printf ("<td>%s</td>\n",currAnnotation->annotation->name);
        printf ("<td>%s:%d-%d</td>\n",
                currAnnotation->annotation->chromosome,
                currAnnotation->annotation->start,
                currAnnotation->annotation->end);
        printf ("<td>%c</td>\n",currAnnotation->annotation->strand);
        printf ("<td>%.2f</td>\n",currAnnotation->percentIntervalOverlap);
        printf ("<td>%.2f</td>\n",currAnnotation->percentAnnotationOverlap);
        puts ("</tr>");
      }
      puts ("</table>");
      puts ("<br>");
      puts ("<br>");
      puts ("<br>");
      puts ("<center><h3>Sample statistics</h3></center>");
      puts ("<br>");
      sampleNames = mio_getSamplesNames ();
      puts ("<table border=1 cellpadding=2 width=80% align=center>");
      puts ("<tr align=left>");
      puts ("<td rowspan=2><b>Feature</b></td>");
      puts ("<td colspan=3><b>PGENE expression</b></td>");
      puts ("<td colspan=3><b>HIT expression</b></td>");
      puts ("<td rowspan=2><b>Pair-wise correlation</b></td>");
      puts ("</tr>");
      puts ("<tr align=left>");
      puts ("<td><b>Average</b></td>");
      puts ("<td><b>Minimum</b></td>");
      puts ("<td><b>Maximum</b></td>");
      puts ("<td><b>Average</b></td>");
      puts ("<td><b>Minimum</b></td>");
      puts ("<td><b>Maximum</b></td>");
      puts ("</tr>");
      for (j = 0; j < arrayMax (sampleNames); j++) {
        puts ("<tr align=left>");
        printf ("<td>%s</td>\n",textItem (sampleNames,j));
        printf ("<td>%.3f</td><td>%.3f</td><td>%.3f</td>\n",
                arru (currStatistic1->averages,j,double),
                arru (currStatistic1->mins,j,double),
                arru (currStatistic1->maxs,j,double));
        printf ("<td>%.3f</td><td>%.3f</td><td>%.3f</td>\n",
                arru (currStatistic2->averages,j,double),
                arru (currStatistic2->mins,j,double),
                arru (currStatistic2->maxs,j,double));
        printf ("<td>%.3f</td>\n",arru (currPairCorrelation->correlations,j,double));
        puts ("</tr>");
      }
      puts ("</table>");
      puts ("<br>");
      puts ("<br>");
      puts ("<br>");
      puts ("<center><h3>Graphic representation</h3></center>");
      puts ("<br>");
      printf ("<div id=tabs-%d style=\"width:90%%;\">\n",index);
      puts ("<ul>");
      printf ("<li><a href=#tabs-%d-1>Unscaled image</a></li>\n",index);
      printf ("<li><a href=#tabs-%d-2>Scaled image</a></li>\n",index);
      puts ("</ul>");
      printf ("<div id=tabs-%d-1>\n",index);
      stringPrintf (buffer1,"%s/%s/%s/%s_%d.unscaled.png",WEB_DATA_URL,dataSet,matrixName,matrixName,pairNumber1);
      printf ("<img src=%s />\n",string (buffer1));
      puts ("</div>");
      printf ("<div id=tabs-%d-2>\n",index);
      stringPrintf (buffer1,"%s/%s/%s/%s_%d.scaled.png",WEB_DATA_URL,dataSet,matrixName,matrixName,pairNumber1);
      printf ("<img src=%s />\n",string (buffer1));
      puts ("</div>");
      puts ("</div>");
      arrayDestroy (annotationPointers1);
      arrayDestroy (annotationPointers2);
      fflush (stdout);
      index++;
    }
  }
  for (i = 0; i < 100; i++) {
    puts ("<br>");
  }
  mio_deInit ();
  puts ("</body>");
  puts ("</html>");
  fflush (stdout);  

}



int main (int argc, char *argv[]) 
{
  char *queryString;
  int first;
  Stringa item = stringCreate (20);
  Stringa value = stringCreate (20);
  char *iPtr,*vPtr;
  char *matrixName = NULL;
  int pairNumber;
  char *dataSet = NULL;
  
  cgiInit();
  cgiHeader("text/html");
  queryString = cgiGet2Post();
  first = 1;
  cgiGetInit ();
  while (cgiGetNextPair (&first,item,value)) {
    iPtr = string (item);
    vPtr = string (value);
    if (strEqual (iPtr,"matrixName")) {
      strReplace (&matrixName,vPtr);
    }
    else if (strEqual (iPtr,"pairNumber")) {
      pairNumber = atoi (vPtr);
    }
    else if (strEqual (iPtr,"dataSet")) {
      strReplace (&dataSet,vPtr);
    }
    else {
      die ("Unexpected inputs: '%s' '%s'\n",iPtr,vPtr); 
    }
  }
  processData (dataSet,matrixName,pairNumber);
  fflush (stdout);
  return 0;
}
