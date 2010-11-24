#include "log.h"
#include "format.h"
#include "html.h"
#include "htmlLinker.h"
#include "common.h"
#include "mio.h"
#include <math.h>



#define WEB_URL_CGI "http://dynamic.gersteinlab.org/people/lh372"
#define WEB_DATA_DIR "/home/lh372/public_html/PseudoSeq"
#define WEB_DATA_SET_ENCODE "ENCODE"
#define WEB_DATA_SET_HBM "HumanBodyMap"
#define WEB_DATA_META_FILE_ENCODE "ENCODE.minBlock100.filtered.all.final.minExpression0.1.score.meta"
#define WEB_DATA_META_FILE_HBM "HumanBodyMap.minBlock100.filtered.all.final.minExpression0.1.score.meta"
#define WEB_DATA_SAMPLE_FILE_ENCODE "samples.txt"
#define WEB_DATA_SAMPLE_FILE_HBM "samples.txt"



typedef struct {
  char *name;
  int pairNumber;
  Interval *alignment1;
  Interval *alignment2;
  Statistic *statistic1;
  Statistic *statistic2;
  int length;
  double percentIdentity;
  PairCorrelation* pairCorrelation;
  Score *score;
  int numberOfAlignmentPairs;
} Pair;



static void printJavaScript (void)
{
  puts ("<meta charset=\"utf-8\">");
  puts ("<link rel=\"stylesheet\" href=\"http://jqueryui.com/themes/base/jquery.ui.all.css\">");
  puts ("<Script src=\"http://code.jquery.com/jquery-1.4.3.js\"></script>");
  puts ("<script src=\"http://jqueryui.com/Ui/jquery.ui.core.js\"></script>");
  puts ("<script src=\"http://jqueryui.com/ui/jquery.ui.widget.js\"></Script>");
  puts ("<script src=\"http://jqueryui.com/ui/jquery.ui.mouse.js\"></script>");
  puts ("<script src=\"http://jqueryui.com/ui/jquery.ui.slider.js\"></script>");
  puts ("<script>");

  puts ("");

  puts ("$(function() {");
  puts ("	$( \"#slider-range\" ).slider({");
  puts ("		range: true,");
  puts ("		min: -100,");
  puts ("		max: 100,");
  puts ("		values: [ -100, 0 ],");
  puts ("		slide: function( event, ui ) {");
  puts ("			$( \"#correlation\" ).val( \"[\" + ui.values[ 0 ] / 100 + \" \" + ui.values[ 1 ] / 100 + \"]\")");
  puts ("		}");
  puts ("	});");
  puts ("});");
  puts ("</script>");
  
  puts ("");

  puts ("<script type=\"text/javascript\">");
  puts ("function hideAll(div_id) {");
  puts ("    changeObjectVisibility(\"scoreField\",\"hidden\");");
  puts ("    changeObjectVisibility(\"expressionField\",\"hidden\");");
  puts ("    $( \"#slider-range\" ).slider (\"values\",0,-100);");
  puts ("    $( \"#slider-range\" ).slider (\"values\",1,0);");
  puts ("}");

  puts ("");

  puts ("function getStyleObject(objectId) {");
  puts ("    // cross-browser function to get an object's style object given its id");
  puts ("    if(document.getElementById && document.getElementById(objectId)) {");
  puts ("	// W3C DOM");
  puts ("	return document.getElementById(objectId).style;");
  puts ("    } else if (document.all && document.all(objectId)) {");
  puts ("	// MSIE 4 DOM");
  puts ("	return document.all(objectId).style;");
  puts ("    } else if (document.layers && document.layers[objectId]) {");
  puts ("	// NN 4 DOM.. note: this won't find nested layers");
  puts ("	return document.layers[objectId];");
  puts ("    } else {");
  puts ("	return false;");
  puts ("    }");
  puts ("} // getStyleObject");

  puts ("");

  puts ("function changeObjectVisibility(objectId, newVisibility) {");
  puts ("    // get a reference to the cross-browser style object and make sure the object exists");
  puts ("    var styleObject = getStyleObject(objectId);");
  puts ("    if(styleObject) {");
  puts ("	styleObject.visibility = newVisibility;");
  puts ("	return true;");
  puts ("    }"); 
  puts ("    else {");
  puts ("	// we couldn't find the object, so we can't change its visibility");
  puts ("	return false;");
  puts ("    }");
  puts ("} // changeObjectVisibility");

  puts ("");

  puts ("</script>");
}



static void printIntroForm (void) 
{
  puts ("<html>");
  puts ("<head>");
  html_printGenericStyleSheet (12);
  puts ("<title>PseudoSeq</title>\n");
  puts ("</head>");
  puts ("<body>");
  puts ("<h1>PseudoSeq: Analysis of transcribed pseudogenes using multiple RNA-Seq samples.</h1><br><br>");
  printf ("<form action=%s/pseudoSeq_cgi method=get>\n",WEB_URL_CGI);
  puts ("<b>Data set</b>:&nbsp;");
  puts ("<select name=mode>");
  puts ("<option value=dataENCODE>ENCODE");
  puts ("<option value=dataHBM selected>Human Body Map");
  puts ("</select>");
  puts ("<br><br><br>");
  puts ("<input type=submit value=Submit>");
  puts ("<input type=reset value=Reset>");
  puts ("</form>");
  puts ("</body>");
  puts ("</html>");
  fflush (stdout);  
}



static void printMainForm (char *dataSet, char *sampleFile) 
{
  int i;
  Texta samples;
  static Stringa buffer = NULL;
  char *pos;

  stringCreateClear (buffer,100);
  puts ("<html>");
  puts ("<head>");
  html_printGenericStyleSheet (12);
  printJavaScript ();
  puts ("<title>PseudoSeq</title>\n");
  puts ("</head>");
  puts ("<body>");
  puts ("<h1>PseudoSeq: Analysis of transcribed pseudogenes using multiple RNA-Seq samples.</h1><br>");
  printf ("<h2>Data set: %s</h2><br><br>\n",dataSet);
  puts ("<table border=1 cellpadding=20 width=100%>");
  puts ("<tr>");
  puts ("<td valign=top>");
  printf ("<form action=%s/pseudoSeq_cgi method=get>\n",WEB_URL_CGI);
  puts ("<b>Select type:</b>&nbsp;");
  puts ("<select name=type>");
  printf ("<option value=%d selected>Expression (PGENE) > Expression (HIT)\n",SCORE_TYPE_TRANSCRIPTION);
  printf ("<option value=%d>Discordant expression patterns\n",SCORE_TYPE_DISCORDANT);
  printf ("<option value=%d>Concordant expression patterns\n",SCORE_TYPE_CONCORDANT);
  printf ("<option value=%d>Select all\n",SCORE_TYPE_ALL);
  puts ("</select>");
  puts ("<br><br><br>");
  puts ("<b>Minimum score:</b>&nbsp;");
  puts ("<input type=radio name=scoreType value=all checked onClick=changeObjectVisibility('scoreField','hidden');>All scores&nbsp;&nbsp;&nbsp;");
  puts ("<input type=radio name=scoreType value=user onClick=changeObjectVisibility('scoreField','visible');>User-defined score");
  puts ("<div id=scoreField style=\"visibility:hidden; display:inline;\">&nbsp;&nbsp;&nbsp;<input type=text name=userScore value=1.0 size=10></div>");
  puts ("<br><br><br>");
  puts ("<b>Minimum alignment length:</b>&nbsp;");
  puts ("<select name=minAlignmentLength>");
  for (i = 100; i <= 1000; i = i + 100) {
    printf ("<option value=%d%s>%d\n",i,i == 100 ? " selected" : "",i);
  }
  puts ("</select>");
  puts ("<br><br><br>");
  puts ("<b>Minimum percent identity:</b>&nbsp;");
  puts ("<select name=minPercentIdentiy>");
  for (i = 90; i <= 100; i++) {
    printf ("<option value=%d%s>%d\n",i,i == 90 ? " selected" : "",i);
  }
  puts ("</select>");
  puts ("<br><br><br>");
  puts ("<b>Range average pair-wise Pearson correlation coefficient:</b>&nbsp;");
  puts ("<input type=text name=correlationValues id=correlation style=\"border:0; color:#f6931f; font-weight:bold;\" value=\"[-1 0]\" size=12/>&nbsp;");
  puts ("<br><br>");
  puts ("<div id=\"slider-range\" style=\"width:75%;\"></div>"); 
  puts ("</td>");
  puts ("<td valign=top>");
  puts ("<b>Maximum number of alignment pairs:</b>&nbsp;");
  puts ("<select name=maxNumAlignmentPairs>");
  for (i = 1; i <= 25; i++) {
    printf ("<option value=%d%s>%d\n",i,i == 25 ? " selected" : "",i);
  }
  puts ("</select>");
  puts ("<br><br><br>");
  puts ("<b>Minimum expression requirement:</b>&nbsp;");
  puts ("<input type=radio name=expressionRequirement value=none checked onClick=changeObjectVisibility('expressionField','hidden');>None&nbsp;&nbsp;&nbsp;");
  puts ("<input type=radio name=expressionRequirement value=user onClick=changeObjectVisibility('expressionField','visible');>User-defined<br>");
  puts ("<div id=expressionField style=visibility:hidden;>");
  puts ("<br>");
  stringPrintf (buffer,"%s/%s/%s",WEB_DATA_DIR,dataSet,sampleFile);
  samples = readList (string (buffer));
  puts ("<table border=0 cellpadding=10 cellspacing=10 width=100%>");
  puts ("<tr>");
  puts ("<td>");
  for (i = 0; i < arrayMax (samples); i++) {
    pos = strchr (textItem (samples,i),'\t');
    if (pos == NULL) {
      die ("Expected '\t' in file with samples!");
    }
    if ((i % 5) == 0 && i > 0) {
      puts ("</td><td>");
    }
    printf ("<input type=checkbox name=samples value=%d>%s<br>\n",i,pos + 1);
  }
  puts ("</td>");
  puts ("</tr>");
  puts ("</table>"); 
  puts ("<br>");
  puts ("<b>Minimum expression value for selected samples:</b>&nbsp;");
  puts ("<input type=text name=minExpressionValue value=\"0.1\">");
  puts ("</div>");
  puts ("</td>");
  puts ("</tr>");
  puts ("</table>"); 
  puts ("<br><br><br>");
  if (strEqual (dataSet,WEB_DATA_SET_ENCODE)) {
    puts ("<input type=hidden name=mode value=analysisENCODE>");
  }
  else  if (strEqual (dataSet,WEB_DATA_SET_HBM)) {
    puts ("<input type=hidden name=mode value=analysisHBM>");
  }
  else {
    die ("Unknown data set: %s",dataSet);
  }
  puts ("<input type=submit value=Submit>");
  puts ("<input type=reset value=Reset onClick=\"hideAll()\";>");
  puts ("</form>");
  puts ("</body>");
  puts ("</html>");
  fflush (stdout);  
}

 

static int isValidScore (Score *currScore, int type, char *scoreType, double userScore) 
{
  if (type == SCORE_TYPE_ALL) {
    return 1;
  }
  if (currScore->type == type && strEqual (scoreType,"all")) {
    return 1;
  }
  if (currScore->type == type && strEqual (scoreType,"user") && currScore->score > userScore) {
    return 1;
  }
  return 0; 
}



static int isValidExpression (char *expressionRequirement, double minExpressionValue, Array samplesIndices, Statistic *currStatistic)
{
  int i;

  if (strEqual (expressionRequirement,"none")) {
    return 1;
  }
  for (i = 0; i < arrayMax (samplesIndices); i++) {
    if (arru (currStatistic->averages,arru (samplesIndices,i,int),double) < minExpressionValue) {
      return 0;
    }
  }
  return 1;
}



static int sortPairs (Pair *a, Pair *b) 
{
  int diff;

  diff = a->score->type - b->score->type;
  if (diff != 0) {
    return diff;
  }
  if (isnan (a->score->score)) {
    return 1;
  }
  if (isnan (b->score->score)) {
    return -1;
  }
  if (a->score->score < b->score->score) {
    return 1;
  }
  if (a->score->score > b->score->score) {
    return -1;
  }
  if (a->statistic1->overallAverage < b->statistic1->overallAverage) {
    return 1;
  }
  if (a->statistic1->overallAverage > b->statistic1->overallAverage) {
    return -1;
  }
  return 0;
}



static void processData (char *dataSet, char *metaFile, char *sampleFile, int type, char *scoreType,
                         double userScore, int minAlignmentLength, double minPercentIdentiy, double minCorrelation,
                         double maxCorrelation, int maxNumAlignmentPairs, char *expressionRequirement, 
                         double minExpressionValue, Array samplesIndices)
{
  Matrix *currMatrix;
  static Stringa buffer1 = NULL;
  static Stringa buffer2 = NULL;
  Array matrices;
  int i,j;
  Array pairs;
  Pair *currPair;
  Interval *intervalPtr1,*intervalPtr2;
  char *name1 = NULL;
  char *name2 = NULL;
  int pairNumber1,pairNumber2;
  int intervalNumber1,intervalNumber2;
  int length;
  double percentIdentity;
  Statistic *currStatistic1;
  Statistic *currStatistic2;
  PairCorrelation *currPairCorrelation;
  Score *currScore;

  puts ("<html>");
  puts ("<head>");
  html_printGenericStyleSheet (12);
  puts ("<title>PseudoSeq</title>\n");
  puts ("</head>");
  puts ("<body>");
  stringCreateClear (buffer1,100);
  stringCreateClear (buffer2,100);
  stringPrintf (buffer1,"%s/%s/%s",WEB_DATA_DIR,dataSet,metaFile);
  stringPrintf (buffer2,"%s/%s/%s",WEB_DATA_DIR,dataSet,sampleFile);  
  pairs = arrayCreate (10000,Pair);
  mio_init (string (buffer1),string (buffer2));
  matrices = mio_parse ();
  for (i = 0; i < arrayMax (matrices); i++) {
    currMatrix = arrp (matrices,i,Matrix);
    for (j = 0; j < arrayMax (currMatrix->intervals); j = j + 2) {
      intervalPtr1 = arrp (currMatrix->intervals,j,Interval);
      intervalPtr2 = arrp (currMatrix->intervals,j + 1,Interval);
      mio_getIntervalInfo (intervalPtr1->name,&name1,&pairNumber1,&intervalNumber1);
      mio_getIntervalInfo (intervalPtr2->name,&name2,&pairNumber2,&intervalNumber2);
      if (pairNumber1 != pairNumber2) {
        die ("Expected same pair!");
      }
      mio_getMetaInfo (currMatrix,pairNumber1,&length,&percentIdentity);
      currStatistic1 = mio_getStatisticPointer (currMatrix,intervalNumber1);
      currStatistic2 = mio_getStatisticPointer (currMatrix,intervalNumber2);
      currPairCorrelation = mio_getPairCorrelationPointer (currMatrix,pairNumber1);
      currScore = mio_getScorePointer (currMatrix,pairNumber1);
      if (isValidScore (currScore,type,scoreType,userScore) == 0) {
        continue;
      }
      if (length < minAlignmentLength || percentIdentity < minPercentIdentiy ||
          currPairCorrelation->averageCorrelation < minCorrelation || currPairCorrelation->averageCorrelation > maxCorrelation) {
        continue;
      }
      if ((arrayMax (currMatrix->intervals) / 2) > maxNumAlignmentPairs) {
        continue;
      }
      if (isValidExpression (expressionRequirement,minExpressionValue,samplesIndices,currStatistic1) == 0) {
        continue;
      }
      currPair = arrayp (pairs,arrayMax (pairs),Pair);
      currPair->name = currMatrix->name;
      currPair->pairNumber = pairNumber1;
      currPair->alignment1 = intervalPtr1;
      currPair->alignment2 = intervalPtr2;   
      currPair->length = length;
      currPair->percentIdentity = percentIdentity;
      currPair->statistic1 = currStatistic1;
      currPair->statistic2 = currStatistic2;
      currPair->pairCorrelation = currPairCorrelation;
      currPair->score = currScore;
      currPair->numberOfAlignmentPairs = arrayMax (currMatrix->intervals) / 2;
    }
  }
  arraySort (pairs,(ARRAYORDERF)sortPairs);
  printf ("<h1><center>Results for data set: %s</center></h1><br>\n",dataSet);
  puts ("<br>");
  puts ("<table border=1 cellpadding=2 width=90% align=center>");
  puts ("<tr align=left>");
  puts ("<th>Name</th>");
  puts ("<th>Pair number</th>");
  puts ("<th>Alignment length</th>");
  puts ("<th>Percent identity</th>");
  puts ("<th>Average pair-wise correlation</th>");
  puts ("<th>Average expression PGENE</th>");
  puts ("<th>Average expression HIT</th>");
  puts ("<th>Number of alignments</th>");
  puts ("<th>Score type</th>");
  puts ("<th>Score</th>");
  puts ("<th>Details</th>");
  puts ("</tr>");
  for (i = 0; i < arrayMax (pairs); i++) {
    currPair = arrp (pairs,i,Pair);
    puts ("<tr align=left>");
    printf ("<td>%s</td>\n",currPair->name);
    printf ("<td>%d</td>\n",currPair->pairNumber);
    printf ("<td>%d</td>\n",currPair->length);
    printf ("<td>%.2f</td>\n",currPair->percentIdentity);
    printf ("<td>%.4f</td>\n",currPair->pairCorrelation->averageCorrelation);
    printf ("<td>%.4f</td>\n",currPair->statistic1->overallAverage);
    printf ("<td>%.4f</td>\n",currPair->statistic2->overallAverage);
    printf ("<td>%d</td>\n",currPair->numberOfAlignmentPairs);
    printf ("<td>%d</td>\n",currPair->score->type);
    printf ("<td>%.2f</td>\n",currPair->score->score);
    printf ("<td><a href=%s/showMatrix_cgi?matrixName=%s&dataSet=%s&pairNumber=%d>Link</a></td>\n",
            WEB_URL_CGI,currPair->name,dataSet,currPair->pairNumber);
    puts ("</tr>");
  }
  puts ("</table>");
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
  char *mode = NULL;
  int type;
  char *scoreType = NULL;
  double userScore;
  int minAlignmentLength;
  double minPercentIdentiy;
  double minCorrelation;
  double maxCorrelation;
  int maxNumAlignmentPairs;
  char *expressionRequirement = NULL;
  double minExpressionValue;
  char *copy = NULL;
  char *pos1,*pos2,*pos3;
  Array samplesIndices;

  cgiInit();
  cgiHeader("text/html");
  queryString = cgiGet2Post();
  if (queryString[0] == '\0') {
    warn ("Wrong URL, use the following <a href=%s/pseudoSeq_cgi?mode=intro>link</a>!",WEB_URL_CGI);
    fflush (stdout);
    return 0;
  }
  samplesIndices = arrayCreate (100,int);
  first = 1;
  cgiGetInit ();
  while (cgiGetNextPair (&first,item,value)) {
    iPtr = string (item);
    vPtr = string (value);
    if (strEqual (iPtr,"mode")) {
      strReplace (&mode,vPtr);
    }
    else if (strEqual (iPtr,"type")) {
      type = atoi (vPtr);
    }
    else if (strEqual (iPtr,"scoreType")) {
      strReplace (&scoreType,vPtr);
    } 
    else if (strEqual (iPtr,"userScore")) {
      userScore = atof (vPtr);
    }
    else if (strEqual (iPtr,"minAlignmentLength")) {
      minAlignmentLength = atoi (vPtr);
    }
    else if (strEqual (iPtr,"minPercentIdentiy")) {
      minPercentIdentiy = atof (vPtr) / 100;
    }
    else if (strEqual (iPtr,"correlationValues")) {
      strReplace (&copy,vPtr);
      pos1 = strchr (copy,'[');
      pos2 = strchr (copy,' ');
      pos3 = strchr (copy,']');
      if (pos1 == NULL || pos2 == NULL || pos3 == NULL) {
        die ("Unexpected correlation values: %p %p %p",pos1,pos2,pos3);
      }
      *pos2 = '\0';
      *pos3 = '\0';
      minCorrelation = atof (pos1 + 1);
      maxCorrelation = atof (pos2 + 1);
    }
    else if (strEqual (iPtr,"maxNumAlignmentPairs")) {
      maxNumAlignmentPairs = atoi (vPtr);
    }
    else if (strEqual (iPtr,"expressionRequirement")) {
      strReplace (&expressionRequirement,vPtr);
    }
    else if (strEqual (iPtr,"minExpressionValue")) {
      minExpressionValue = atof (vPtr);
    }
    else if (strEqual (iPtr,"samples")) {
      array (samplesIndices,arrayMax (samplesIndices),int) = atoi (vPtr);
    }
    else {
      die ("Unexpected inputs: '%s' '%s'\n",iPtr,vPtr); 
    }
  }
  if (strEqual (mode,"intro")) {
    printIntroForm ();
  }
  else if (strEqual (mode,"dataENCODE")) {
    printMainForm (WEB_DATA_SET_ENCODE,WEB_DATA_SAMPLE_FILE_ENCODE);
  }
  else if (strEqual (mode,"dataHBM")) {
    printMainForm (WEB_DATA_SET_HBM,WEB_DATA_SAMPLE_FILE_HBM);
  }
  else if (strEqual (mode,"analysisENCODE")) {
  
  }
  else if (strEqual (mode,"analysisHBM")) {
    processData (WEB_DATA_SET_HBM,WEB_DATA_META_FILE_HBM,WEB_DATA_SAMPLE_FILE_HBM,
                 type,scoreType,userScore,minAlignmentLength,minPercentIdentiy,minCorrelation,maxCorrelation,
                 maxNumAlignmentPairs,expressionRequirement,minExpressionValue,samplesIndices);
  }
  else {
    die ("Unknown Mode: '%s'",mode);
  }
  fflush (stdout);
  arrayDestroy (samplesIndices);
  return 0;
}
