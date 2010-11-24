#include "log.h"
#include "format.h"
#include "numUtil.h"
#include "intervalFind.h"
#include "bits.h"
#include "mio.h"
#include "gd.h"
#include "gdfontmb.h"
#include "gdfonts.h"
#include "gdfontt.h"
#include <gsl/gsl_statistics.h>



#define IMAGE_HEIGHT 1200
#define IMAGE_WIDTH 1000
#define MARGIN_TOP 20
#define MARGIN_BOTTOM 20
#define MARGIN_LEFT 20
#define MARGIN_RIGHT 20
#define SIZE_LEFT_LABEL 110
#define SIZE_RIGHT_LABEL 85
#define SIZE_TOP_AREA 50
#define VERTICAL_TOP_LABEL_OFFSET 15
#define LEFT_LABEL_OFFSET 50
#define RIGHT_LABEL_OFFSET 10
#define SIZE_ANNOTATION 6
#define MAX_NUM_ANNOTATIONS_PER_INTERVAL 10
#define FRACTION_INTRONIC 0.25



typedef struct {
  int genomic;
  int pixel;
} Coordinate; 



typedef struct {
  int start;
  int end;
  int size;
  int color;
} Block;



static double min3 (float a, float b, float c) 
{
  float min;

  min = a;
  if (b < min) {
    min = b;
  }
  if (c < min) {
    min = c;
  }
  return min;
}



static double max3 (float a, float b, float c) 
{
  float max;

  max = a;
  if (b > max) {
    max = b;
  }
  if (c > max) {
    max = c;
  }
  return max;
}



// r,g,b values are from 0 to 255
// h = [0,360], s = [0,1], v = [0,1]
//		if s == 0, then h = -1 (undefined)
void rgb2hsv (float r, float g, float b, float *h, float *s, float *v )
{
  float min, max, delta;

  r = r / 255;
  g = g / 255;
  b = b / 255;
  min = min3 ( r, g, b );
  max = max3 ( r, g, b );
  *v = max;			// v
  delta = max - min;
  if( max != 0 )
    *s = delta / max;		// s
  else {
    // r = g = b = 0		// s = 0, v is undefined
    *s = 0;
    *h = -1;
    return;
  }
  if( r == max )
    *h = ( g - b ) / delta;	// between yellow & magenta
  else if( g == max )
    *h = 2 + ( b - r ) / delta;	// between cyan & yellow
  else
    *h = 4 + ( r - g ) / delta;	// between magenta & cyan
  *h *= 60;			// degrees
  if( *h < 0 )
    *h += 360;
}



void hsv2rgb (float *r, float *g, float *b, float h, float s, float v)
{
  int i;
  float f, p, q, t;
	
  if( s == 0 ) {
    // achromatic (grey)
    *r = *g = *b = v;
    return;
  }
  h /= 60;			// sector 0 to 5
  i = floor( h );
  f = h - i;			// factorial part of h
  p = v * ( 1 - s );
  q = v * ( 1 - s * f );
  t = v * ( 1 - s * ( 1 - f ) );
  switch( i ) {
  case 0:
    *r = v;
    *g = t;
    *b = p;
    break;
  case 1:
    *r = q;
    *g = v;
    *b = p;
    break;
  case 2:
    *r = p;
    *g = v;
    *b = t;
    break;
  case 3:
    *r = p;
    *g = q;
    *b = v;
    break;
  case 4:
    *r = t;
    *g = p;
    *b = v;
    break;
  default:		// case 5:
    *r = v;
    *g = p;
    *b = q;
    break;
  }
  *r = *r * 255;
  *g = *g * 255;
  *b = *b * 255;
}



static void drawTrack (gdImagePtr im, Matrix *currMatrix, Array rowIndices, int columnIndex, 
                       GraphCoordTrans gct, double min, double max,
                       int trackNumber, int trackHeight, int heightAnnotations, int color)
{
  int i,j;
  int count;
  int peakHeight;
  double sumValues;
  double value;

  i = 0; 
  while (i < arrayMax (rowIndices)) {
    value = mio_getValue (currMatrix,arru (rowIndices,i,int),columnIndex);
    if (min != max) {
      value = (value - min) / (max - min);
    }
    sumValues = value;
    count = 1;
    j = i + 1;
    while (j < arrayMax (rowIndices)) {
      if (gr_ct_toPix (gct,i) == gr_ct_toPix (gct,j)) {
        value = mio_getValue (currMatrix,arru (rowIndices,j,int),columnIndex);
        if (min != max) {
          value = (value - min) / (max - min);
        }
        sumValues += value;
        count++;
      }
      else {
        break;
      }
      j++;
    }
    peakHeight = trackHeight * (sumValues / count);
    gdImageLine (im,
                 MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (trackNumber * trackHeight),
                 MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (trackNumber * trackHeight) - peakHeight,
                 color);
    i = j;
  }
}



static double calculateAverage (Array values)
{
  int i;
  double sum;

  sum = 0.0;
  for (i = 0; i < arrayMax (values); i++) {
    sum = sum + arru (values,i,double);
  }
  return sum / arrayMax (values);
}



static int sortCoordinates (Coordinate *a, Coordinate *b)
{
  return a->genomic - b->genomic;
}



static Array generateCoordinateSystem (Array annotationPointers, int numAnnotations, int width)
{
  Annotation *currAnnotation;
  Interval *currTranscript;
  int min,max;
  int k,l,m;
  Bits *bits;
  int bitCount;
  SubInterval *currExon;
  Array coordinates;
  Coordinate *currCoordinate;
  Array exonStarts,exonEnds;
  int compositeTranscriptLength;
  int intronSize;
  int totalExonSize;
  int exonSize;
  GraphCoordTrans gct;
  int start;

  exonStarts = arrayCreate (100,int);
  exonEnds = arrayCreate (100,int);
  currAnnotation = arru (annotationPointers,0,Annotation*);
  currTranscript = currAnnotation->annotation;
  min = currTranscript->start;
  max = currTranscript->end;
  for (k = 1; k < numAnnotations; k++) {
    currAnnotation = arru (annotationPointers,k,Annotation*);
    currTranscript = currAnnotation->annotation;
    if (currTranscript->start < min) {
      min = currTranscript->start;
    }
    if (currTranscript->end > max) {
      max = currTranscript->end;
    }
  }
  bitCount = max - min + 1;
  bits = bitAlloc (bitCount);
  for (k = 0; k < numAnnotations; k++) {
    currAnnotation = arru (annotationPointers,k,Annotation*);
    currTranscript = currAnnotation->annotation;
    for (l = 0; l < arrayMax (currTranscript->subIntervals); l++) {
      currExon = arrp (currTranscript->subIntervals,l,SubInterval);
      for (m = currExon->start; m <= currExon->end; m++) {
        bitSetOne (bits,m - min);
      }
    }
  }
  if (bitReadOne (bits,0) == 1) {
    array (exonStarts,arrayMax (exonStarts),int) = min;
  }
  for (k = min; k < max; k++) {
    if (bitReadOne (bits,k - min) == 1 &&
        bitReadOne (bits,k - min + 1) == 0) {
      array (exonEnds,arrayMax (exonEnds),int) = k;
    }
    if (bitReadOne (bits,k - min) == 0 &&
        bitReadOne (bits,k - min + 1) == 1) {
      array (exonStarts,arrayMax (exonStarts),int) = k + 1;
    }
  }
  if (bitReadOne (bits,bitCount - 1) == 1) {
    array (exonEnds,arrayMax (exonEnds),int) = max;
  }
  bitFree (&bits);
  if (arrayMax (exonStarts) != arrayMax (exonEnds)) {
    die ("exonStarts and exonEnds must have the same number of elements");
  }
  
  compositeTranscriptLength = 0;
  for (k = 0; k < arrayMax (exonStarts); k++) {
    compositeTranscriptLength += (arru (exonEnds,k,int) - arru (exonStarts,k,int));
  }
  
  intronSize = (double)width * FRACTION_INTRONIC / (arrayMax (exonStarts) - 1);
  totalExonSize = width - (intronSize * (arrayMax (exonStarts) - 1));
  coordinates = arrayCreate (10000,Coordinate);
  start = 0;
  for (k = 0; k < arrayMax (exonStarts); k++) {
    exonSize = (double)(arru (exonEnds,k,int) - arru (exonStarts,k,int)) / compositeTranscriptLength * totalExonSize;
    gct = gr_ct_create ((double)arru (exonStarts,k,int),(double)arru (exonEnds,k,int),start,start + exonSize);
    for (l = arru (exonStarts,k,int); l < arru (exonEnds,k,int); l++) {     
      currCoordinate = arrayp (coordinates,arrayMax (coordinates),Coordinate);
      currCoordinate->genomic = l;
      currCoordinate->pixel = gr_ct_toPix (gct,l);
    }
    gr_ct_destroy (gct);
    start += (exonSize + intronSize);
  }
  arraySort (coordinates,(ARRAYORDERF)sortCoordinates);
  arrayDestroy (exonStarts);
  arrayDestroy (exonEnds);
  return coordinates;
}



static int genomicPosition2pixel (Array coordinates, int genomicCoordinate) 
{
  Coordinate testCoordinate;
  int index;
 
  testCoordinate.genomic = genomicCoordinate;
  if (!arrayFind (coordinates,&testCoordinate,&index,(ARRAYORDERF)sortCoordinates)) {
    die ("Expected to find coordinate: %d",genomicCoordinate); 
  }
  return arrp (coordinates,index,Coordinate)->pixel;
}



static void drawAnnotations (gdImagePtr im, Matrix *currMatrix, Array blocks, Array annotationPointers, 
                             int numAnnotations, Array coordinates, int width, int baseColor, int *iter)
{
  int i,j,k;
  Annotation *currAnnotation;
  Interval *currInterval;
  SubInterval *prevSubInterval,*currSubInterval;
  int yCoordinate;
  int thisIter;
  int overlap;
  static Stringa buffer = NULL;
  char *pos;
  Block *currBlock;

  stringCreateClear (buffer,100);
  thisIter = *iter;
  for (i = 0; i < numAnnotations; i++) {
    currAnnotation = arru (annotationPointers,i,Annotation*);
    currInterval = currAnnotation->annotation;
    yCoordinate = thisIter * 2 * SIZE_ANNOTATION;
    for (j = 0; j < arrayMax (currInterval->subIntervals); j++) {
      currSubInterval = arrp (currInterval->subIntervals,j,SubInterval);
      gdImageFilledRectangle (im,
                              MARGIN_LEFT + SIZE_LEFT_LABEL + genomicPosition2pixel (coordinates,currSubInterval->start),
                              MARGIN_TOP + SIZE_TOP_AREA + yCoordinate - SIZE_ANNOTATION / 2,
                              MARGIN_LEFT + SIZE_LEFT_LABEL + genomicPosition2pixel (coordinates,currSubInterval->end - 1),
                              MARGIN_TOP + SIZE_TOP_AREA + yCoordinate + SIZE_ANNOTATION / 2,
                              baseColor);
      for (k = 0; k < arrayMax (blocks); k++) {
        currBlock =  arrp (blocks,k,Block);
        overlap = rangeIntersection (currSubInterval->start,currSubInterval->end - 1,currBlock->start,currBlock->end - 1);
        if (overlap > 0) {
          gdImageFilledRectangle (im,
                                  MARGIN_LEFT + SIZE_LEFT_LABEL + genomicPosition2pixel (coordinates,MAX (currSubInterval->start,currBlock->start)),
                                  MARGIN_TOP + SIZE_TOP_AREA + yCoordinate - SIZE_ANNOTATION / 2,
                                  MARGIN_LEFT + SIZE_LEFT_LABEL + genomicPosition2pixel (coordinates,MIN (currSubInterval->end - 1,currBlock->end - 1)),
                                  MARGIN_TOP + SIZE_TOP_AREA + yCoordinate + SIZE_ANNOTATION / 2,
                                  currBlock->color);
        }
      }
    }
    pos = strchr (currMatrix->name,'_');
    if (pos == NULL) {
      die ("Expected to find '_' in matrix name!");
    }
    if (strstr (currInterval->name,pos + 1)) {
      stringPrintf (buffer,"<%s>",currAnnotation->nameAnnotationSet);
    }
    else {
      stringPrintf (buffer,"%s",currAnnotation->nameAnnotationSet);
    }
    gdImageString (im,
                   gdFontGetMediumBold (),
                   MARGIN_LEFT,
                   MARGIN_TOP + SIZE_TOP_AREA + yCoordinate - gdFontGetMediumBold ()->h / 2,
                   (unsigned char*)string (buffer),
                   baseColor);
    stringPrintf (buffer,"(%.0f,%.0f)",currAnnotation->percentIntervalOverlap,currAnnotation->percentAnnotationOverlap);
    gdImageString (im,
                   gdFontGetMediumBold (),
                   MARGIN_LEFT + SIZE_LEFT_LABEL + width + RIGHT_LABEL_OFFSET,
                   MARGIN_TOP + SIZE_TOP_AREA + yCoordinate - gdFontGetMediumBold ()->h / 2,
                   (unsigned char*)string (buffer),
                   baseColor);
    thisIter++;
  }
  thisIter = *iter;
  for (i = 0; i < numAnnotations; i++) {
    currAnnotation = arru (annotationPointers,i,Annotation*);
    currInterval = currAnnotation->annotation;
    yCoordinate = thisIter * 2 * SIZE_ANNOTATION;
    for (j = 1; j < arrayMax (currInterval->subIntervals); j++) {
      prevSubInterval = arrp (currInterval->subIntervals,j - 1,SubInterval);
      currSubInterval = arrp (currInterval->subIntervals,j,SubInterval);
      gdImageLine (im,
                   MARGIN_LEFT + SIZE_LEFT_LABEL + genomicPosition2pixel (coordinates,prevSubInterval->end - 1),
                   MARGIN_TOP + SIZE_TOP_AREA + yCoordinate,
                   MARGIN_LEFT + SIZE_LEFT_LABEL + genomicPosition2pixel (coordinates,currSubInterval->start),
                   MARGIN_TOP + SIZE_TOP_AREA + yCoordinate,
                   baseColor);
    }
    thisIter++;
  }
  *iter = thisIter + 1;
}



static Array getAlignmentBlocks (Interval *alignmentInterval)
{
  Array blocks;
  int i;
  Block *currBlock;
  SubInterval *alignmentSubInterval;

  blocks = arrayCreate (arrayMax (alignmentInterval->subIntervals),Block);
  for (i = 0; i < arrayMax (alignmentInterval->subIntervals); i++) {
    alignmentSubInterval =  arrp (alignmentInterval->subIntervals,i,SubInterval);
    currBlock = arrayp (blocks,arrayMax (blocks),Block);
    currBlock->start = alignmentSubInterval->start;
    currBlock->end = alignmentSubInterval->end;
    currBlock->size = (alignmentSubInterval->end - alignmentSubInterval->start);
  }
  return blocks;
}



static void reverseAlignmentBlocks (Array blocks)
{
  Block temp;
  int i;

  for (i = 0; i < arrayMax (blocks) / 2; i++) {
    temp = arru (blocks,i,Block);
    arru (blocks,i,Block) =  arru (blocks,arrayMax (blocks) - 1 - i,Block);
    arru (blocks,arrayMax (blocks) - 1 - i,Block) = temp;
  }
}



static void assignColorsToBlocks (Array blocks, Array colors)
{
  int i;
  Block *currBlock;
  
  if (arrayMax (blocks) != arrayMax (colors)) {
    die ("Expected same number of blocks and colors");
  }
  for (i = 0; i < arrayMax (blocks); i++) {
    currBlock = arrp (blocks,i,Block);
    currBlock->color = arru (colors,i,int);
  }
}



static void drawOneAlignmentBar (gdImagePtr im, Array blocks, int yCoordinate, int length, 
                                 int trackWidth, int baseColor)
{
  int prevPosition;
  int size;
  int i;
  Block *currBlock;

  prevPosition = 0;
  for (i = 0; i < arrayMax (blocks); i++) {
    currBlock = arrp (blocks,i,Block);
    size = round ((double)currBlock->size / length * trackWidth);
    gdImageFilledRectangle (im,
                            MARGIN_LEFT + SIZE_LEFT_LABEL + prevPosition,
                            MARGIN_TOP + SIZE_TOP_AREA + yCoordinate - SIZE_ANNOTATION / 2,
                            MARGIN_LEFT + SIZE_LEFT_LABEL + prevPosition + size,
                            MARGIN_TOP + SIZE_TOP_AREA + yCoordinate + SIZE_ANNOTATION / 2,
                            currBlock->color);
    gdImageRectangle (im,
                      MARGIN_LEFT + SIZE_LEFT_LABEL + prevPosition,
                      MARGIN_TOP + SIZE_TOP_AREA + yCoordinate - SIZE_ANNOTATION / 2,
                      MARGIN_LEFT + SIZE_LEFT_LABEL + prevPosition + size,
                      MARGIN_TOP + SIZE_TOP_AREA + yCoordinate + SIZE_ANNOTATION / 2,
                      baseColor);
    prevPosition += size;
  }
}



static void drawAlignmentBars (gdImagePtr im, Array blocks1, Array blocks2, int length, 
                               int trackWidth, int baseColor, int *iter)
{
  int yCoordinate;
  static Stringa buffer = NULL;

  stringCreateClear (buffer,100);
  yCoordinate = *iter * 2 * SIZE_ANNOTATION;
  drawOneAlignmentBar (im,blocks1,yCoordinate,length,trackWidth,baseColor);
  yCoordinate = ((*iter + 1) * 2 * SIZE_ANNOTATION) - SIZE_ANNOTATION;
  drawOneAlignmentBar (im,blocks2,yCoordinate,length,trackWidth,baseColor);
  yCoordinate = ((*iter + 1) * 2 * SIZE_ANNOTATION) - 2 * SIZE_ANNOTATION;
  stringPrintf (buffer,"%6s","5'");
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT + SIZE_LEFT_LABEL - LEFT_LABEL_OFFSET,
                 MARGIN_TOP + SIZE_TOP_AREA + yCoordinate - (gdFontGetMediumBold ()->h / 2),
                 (unsigned char*)string(buffer),
                 baseColor);
  stringPrintf (buffer,"%s","3'");
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth + RIGHT_LABEL_OFFSET,
                 MARGIN_TOP + SIZE_TOP_AREA + yCoordinate - (gdFontGetMediumBold ()->h / 2),
                 (unsigned char*)string(buffer),
                 baseColor);
}



static Array generateColors (gdImagePtr im, int numColors, float lowH, float highH, float s, float v)
{
  Array colors;
  float step;
  int i;
  float r,g,b;
  float currentH;

  colors = arrayCreate (numColors,int);
  step = (highH - lowH) / numColors;
  currentH = lowH;
  for (i = 0; i < numColors; i++) {
    hsv2rgb (&r,&g,&b,currentH,s,v);
    array (colors,arrayMax (colors),int) = gdImageColorAllocate (im,r,g,b);
    currentH += step;
  }
  return colors;
}



static void drawCorrelationPlot (gdImagePtr im, Matrix *currMatrix, Array rowIndicesInterval1, Array rowIndicesInterval2, int heightAnnotations,
                                 GraphCoordTrans gct, int trackNumber, int trackHeight, int trackWidth, int baseColor)
{
  int baseLine;
  int i,j;
  double *values1,*values2;
  double correlation;
  static Stringa buffer = NULL;

  stringCreateClear (buffer,100);
  baseLine = trackNumber + 1;
  for (i = 0; i < arrayMax (rowIndicesInterval1); i++) { 
    values1 = (double*)malloc (mio_getNumSamples ()*sizeof(double));
    values2 = (double*)malloc (mio_getNumSamples ()*sizeof(double));
    for (j = 0; j < mio_getNumSamples (); j++) {
      values1[j] = mio_getValue (currMatrix,arru (rowIndicesInterval1,i,int),j);
      values2[j] = mio_getValue (currMatrix,arru (rowIndicesInterval2,i,int),j);
    }
    correlation = gsl_stats_correlation (values1,1,values2,1,mio_getNumSamples ());
    free (values1);
    free (values2);
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine * trackHeight) - (trackHeight * correlation),
                     baseColor);
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine * trackHeight) - (trackHeight * correlation) + 1,
                     baseColor);
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine * trackHeight) - (trackHeight * correlation) - 1,
                     baseColor);
  } 
  gdImageLine (im,
                   MARGIN_LEFT + SIZE_LEFT_LABEL,
               MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine - 1) * trackHeight,
               MARGIN_LEFT + SIZE_LEFT_LABEL,
               MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine + 1) * trackHeight,
               baseColor);
  gdImageLine (im,
               MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth,
               MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine - 1) * trackHeight,
               MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth,
               MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine + 1) * trackHeight,
               baseColor);
  gdImageDashedLine (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL,
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + baseLine * trackHeight,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth,
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + baseLine * trackHeight,
                     baseColor);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT,
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + baseLine * trackHeight - gdFontGetMediumBold ()->h / 2,
                 (unsigned char*)"Correlation",
                 baseColor);
  stringPrintf (buffer,"%6d",1);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT + SIZE_LEFT_LABEL - LEFT_LABEL_OFFSET,
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine - 1) * trackHeight - gdFontGetMediumBold ()->h / 2,
                 (unsigned char*)string(buffer),
                 baseColor);
  stringPrintf (buffer,"%6d",-1);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT + SIZE_LEFT_LABEL - LEFT_LABEL_OFFSET,
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine + 1) * trackHeight - gdFontGetMediumBold ()->h / 2,
                 (unsigned char*)string (buffer),
                 baseColor);
}



static void drawExpressionPlot (gdImagePtr im, Matrix *currMatrix, Array rowIndicesInterval1, Array rowIndicesInterval2, int heightAnnotations,
                                GraphCoordTrans gct, int trackNumber, int trackHeight, int trackWidth, int baseColor, int color1, int color2, int isScaled)
{
  int i;
  Row *r1,*r2;
  int baseLine;
  double v1,v2;
  double min1,max1,min2,max2;
  static Stringa buffer = NULL;

  stringCreateClear (buffer,100);
  baseLine = trackNumber + 4;
  for (i = 0; i < arrayMax (rowIndicesInterval1); i++) { 
    r1 = arrp (currMatrix->rows,arru (rowIndicesInterval1,i,int),Row);
    r2 = arrp (currMatrix->rows,arru (rowIndicesInterval2,i,int),Row);
    v1 = calculateAverage (r1->values);
    v2 = calculateAverage (r2->values);
    if (i == 0) {
      min1 = v1;
      max1 = v1;  
      min2 = v2;
      max2 = v2;
      continue;
    }
    else {
      if (v1 < min1) {
        min1 = v1;
      }
      if (v1 > max1) {
        max1 = v1;
      }
      if (v2 < min2) {
        min2 = v2;
      }
      if (v2 > max2) {
        max2 = v2;
      }
    }
  }
  if (isScaled == 0) {
    min1 = min2 = MIN (min1,min2);
    max1 = max2 = MAX (max1,max2);
  }
  for (i = 0; i < arrayMax (rowIndicesInterval1); i++) { 
    r1 = arrp (currMatrix->rows,arru (rowIndicesInterval1,i,int),Row);
    v1 = calculateAverage (r1->values);
    if (min1 != max1) {
      v1 = (v1 - min1) / (max1 - min1);
    }
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + ((baseLine + 1) * trackHeight) - (2 * trackHeight * v1),
                     color1);
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + ((baseLine + 1) * trackHeight) - (2 * trackHeight * v1) + 1,
                     color1);
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + ((baseLine + 1) * trackHeight) - (2 * trackHeight * v1) - 1,
                     color1);
    r2 = arrp (currMatrix->rows,arru (rowIndicesInterval2,i,int),Row);
    v2 = calculateAverage (r2->values);
    if (min2 != max2) {
      v2 = (v2 - min2) / (max2 - min2);
    }
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + ((baseLine + 1) * trackHeight) - (2 * trackHeight * v2),
                     color2);
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + ((baseLine + 1) * trackHeight) - (2 * trackHeight * v2) + 1,
                     color2);
    gdImageSetPixel (im,
                     MARGIN_LEFT + SIZE_LEFT_LABEL + gr_ct_toPix (gct,i),
                     MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + ((baseLine + 1) * trackHeight) - (2 * trackHeight * v2) - 1,
                     color2);
  } 
  gdImageLine (im,
               MARGIN_LEFT + SIZE_LEFT_LABEL,
               MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine - 1) * trackHeight,
               MARGIN_LEFT + SIZE_LEFT_LABEL,
               MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine + 1) * trackHeight,
               baseColor);
  gdImageLine (im,
               MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth,
               MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine - 1) * trackHeight,
               MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth,
               MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine + 1) * trackHeight,
               baseColor);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT,
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + baseLine * trackHeight - gdFontGetMediumBold ()->h / 2,
                 (unsigned char*)"Expression",
                 baseColor);
  stringPrintf (buffer,"%6.1f",max1);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT + SIZE_LEFT_LABEL - LEFT_LABEL_OFFSET,
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine - 1) * trackHeight - gdFontGetMediumBold ()->h / 2,
                 (unsigned char*)string (buffer),
                 color1);
  stringPrintf (buffer,"%6.1f",min1);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT + SIZE_LEFT_LABEL - LEFT_LABEL_OFFSET,
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine + 1) * trackHeight - gdFontGetMediumBold ()->h / 2,
                 (unsigned char*)string (buffer),
                 color1);
  stringPrintf (buffer,"%.1f",max2);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth + RIGHT_LABEL_OFFSET,
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine - 1) * trackHeight - gdFontGetMediumBold ()->h / 2,
                 (unsigned char*)string (buffer),
                 color2);
  stringPrintf (buffer,"%.1f",min2);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth + RIGHT_LABEL_OFFSET,
                 MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (baseLine + 1) * trackHeight - gdFontGetMediumBold ()->h / 2,
                 (unsigned char*)string (buffer),
                 color2);
}



static void  drawBoxesAroundSampleTracks (gdImagePtr im, int heightAnnotations, int trackHeight, int trackWidth, int baseColor)
{
  int i;

  for (i = 0; i < mio_getNumSamples () * 2; i++) {
    gdImageRectangle (im,
                      MARGIN_LEFT + SIZE_LEFT_LABEL,
                      MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + i * trackHeight,
                      MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth,
                      MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (i + 1) * trackHeight,
                      baseColor);
  }
}



static void generateSampleLabels (gdImagePtr im, Matrix *currMatrix, int heightAnnotations, int trackHeight, int trackWidth, int baseColor)
{
  int i;
  Column *currColumn;

  for (i = 0; i < mio_getNumSamples (); i++) {
    if (i < mio_getNumSamples () - 1) {
      gdImageDashedLine (im,
                         MARGIN_LEFT,
                         MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (2 * i + 1) * trackHeight + trackHeight,
                         MARGIN_LEFT + SIZE_LEFT_LABEL,
                         MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (2 * i + 1) * trackHeight + trackHeight,
                         baseColor);
    }
  }
  for (i = 0; i < mio_getNumSamples (); i++) {
    currColumn = arrp (currMatrix->columns,i,Column);
    gdImageString (im,
                   gdFontGetMediumBold (),
                   MARGIN_LEFT,
                   MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (2 * i + 1) * trackHeight - gdFontGetMediumBold ()->h / 2,
                   (unsigned char*)currColumn->sampleName,
                   baseColor);
  }
}



static void  generatePairwiseCorrelationsLabels (gdImagePtr im, Matrix *currMatrix, int pairNumber, int heightAnnotations, 
                                                 int trackHeight, int trackWidth, int baseColor)
{
  int i;
  PairCorrelation *currPairCorrelation;
  static Stringa buffer = NULL;

  stringCreateClear (buffer,100);
  for (i = 0; i < mio_getNumSamples (); i++) {
    if (i < mio_getNumSamples () - 1) {
      gdImageDashedLine (im,
                         MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth,
                         MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (2 * i + 1) * trackHeight + trackHeight,
                         MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth + SIZE_RIGHT_LABEL,
                         MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (2 * i + 1) * trackHeight + trackHeight,
                         baseColor);
    }
  }
  currPairCorrelation = mio_getPairCorrelationPointer (currMatrix,pairNumber);
  for (i = 0; i < arrayMax (currPairCorrelation->correlations); i++) {
    stringPrintf (buffer,"%.3f",arru(currPairCorrelation->correlations,i,double));
    gdImageString (im,
                   gdFontGetMediumBold (),
                   MARGIN_LEFT + SIZE_LEFT_LABEL + trackWidth + RIGHT_LABEL_OFFSET,
                   MARGIN_TOP + SIZE_TOP_AREA + heightAnnotations + (2 * i + 1) * trackHeight - gdFontGetMediumBold ()->h / 2,
                   (unsigned char*)string (buffer),
                   baseColor);
  }
}



static void generateTopLabels (gdImagePtr im, Matrix *currMatrix, int pairNumber, Statistic *currStatistic1, Statistic *currStatistic2,
                               int length, double percentIdentity, int baseColor, int color1, int color2, int isScaled)
{
  PairCorrelation *currPairCorrelation;
  static Stringa buffer = NULL;

  stringCreateClear (buffer,100);
  currPairCorrelation = mio_getPairCorrelationPointer (currMatrix,pairNumber);
  stringPrintf (buffer,
                "Name: %s_%d, Length: %d, %%ID: %.1f, APC: %.3f",
                currMatrix->name,pairNumber,length,percentIdentity * 100,currPairCorrelation->averageCorrelation);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT,
                 MARGIN_TOP,
                 (unsigned char*)string (buffer),
                 baseColor);
 
  stringPrintf (buffer,"Scale [%.2f %.2f], Average: %.2f",
                isScaled == 1 ? currStatistic1->overallMin : MIN (currStatistic1->overallMin,currStatistic2->overallMin),
                isScaled == 1 ? currStatistic1->overallMax : MAX (currStatistic1->overallMax,currStatistic2->overallMax),
                currStatistic1->overallAverage);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT,
                 MARGIN_TOP + VERTICAL_TOP_LABEL_OFFSET,
                 (unsigned char*)string (buffer),
                 color1);
  stringPrintf (buffer,"Scale [%.2f %.2f], Average: %.2f",
                isScaled == 1 ? currStatistic2->overallMin : MIN (currStatistic1->overallMin,currStatistic2->overallMin),
                isScaled == 1 ? currStatistic2->overallMax : MAX (currStatistic1->overallMax,currStatistic2->overallMax),
                currStatistic2->overallAverage);
  gdImageString (im,
                 gdFontGetMediumBold (),
                 MARGIN_LEFT,
                 MARGIN_TOP + 2 * VERTICAL_TOP_LABEL_OFFSET,
                 (unsigned char*)string (buffer),
                 color2);
}



int main (int argc, char *argv[])
{ 
  Matrix *currMatrix;
  gdImagePtr im1,im2;
  FILE *pngout;
  int black1,white1,red1,blue1;
  int black2,white2,red2,blue2;
  int trackHeight,trackWidth;
  int i,j;
  GraphCoordTrans gct;
  Stringa buffer;
  int trackNumber;
  Array rowIndicesInterval1,rowIndicesInterval2;
  Interval *intervalPtr1,*intervalPtr2;
  char *name1 = NULL;
  char *name2 = NULL;
  int pairNumber1,pairNumber2;
  int intervalNumber1,intervalNumber2;
  int length;
  double percentIdentity;
  Array annotationPointers1,annotationPointers2;
  int heightAnnotations;
  int numAnnotations1,numAnnotations2;
  int iter;
  Array colors11,colors12,colors21,colors22;
  Array coordinates1,coordinates2;
  Statistic *currStatistic1,*currStatistic2;
  Array blocks1,blocks2;

  if (argc != 2) {
    usage ("%s <samples.txt>",argv[0]);
  }

  buffer = stringCreate (100);
  mio_init ("-",argv[1]);
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
    for (i = 0; i < arrayMax (currMatrix->intervals); i = i + 2) {
      
      im1 = gdImageCreate (IMAGE_WIDTH,IMAGE_HEIGHT);
      white1 = gdImageColorAllocate (im1,255,255,255);
      black1 = gdImageColorAllocate (im1,0,0,0);
      red1 = gdImageColorAllocate (im1,255,0,0);
      blue1 = gdImageColorAllocate (im1,0,0,255);
      
      im2 = gdImageCreate (IMAGE_WIDTH,IMAGE_HEIGHT);
      white2 = gdImageColorAllocate (im2,255,255,255);
      black2 = gdImageColorAllocate (im2,0,0,0);
      red2 = gdImageColorAllocate (im2,255,0,0);
      blue2 = gdImageColorAllocate (im2,0,0,255);

      intervalPtr1 = arrp (currMatrix->intervals,i,Interval);
      intervalPtr2 = arrp (currMatrix->intervals,i + 1,Interval);
      mio_getIntervalInfo (intervalPtr1->name,&name1,&pairNumber1,&intervalNumber1);
      mio_getIntervalInfo (intervalPtr2->name,&name2,&pairNumber2,&intervalNumber2);
      if (pairNumber1 != pairNumber2) {
        die ("Expected same pair!");
      }
      mio_getMetaInfo (currMatrix,pairNumber1,&length,&percentIdentity);
      rowIndicesInterval1 = mio_getRowIndicesForInterval (currMatrix,intervalNumber1);
      rowIndicesInterval2 = mio_getRowIndicesForInterval (currMatrix,intervalNumber2);
      if (arrayMax (rowIndicesInterval1) != arrayMax (rowIndicesInterval2)) {
        die ("Unequal length: %d %d",arrayMax (rowIndicesInterval1),arrayMax (rowIndicesInterval2));
      }   

      colors11 = generateColors (im1,intervalPtr1->subIntervalCount,0,60,1,1);
      colors12 = generateColors (im1,intervalPtr2->subIntervalCount,240,300,1,1);
      colors21 = generateColors (im2,intervalPtr1->subIntervalCount,0,60,1,1);
      colors22 = generateColors (im2,intervalPtr2->subIntervalCount,240,300,1,1);
     
      annotationPointers1 = mio_getAnnotationPointers (currMatrix,intervalNumber1);
      annotationPointers2 = mio_getAnnotationPointers (currMatrix,intervalNumber2);
      numAnnotations1 = MIN (arrayMax (annotationPointers1),MAX_NUM_ANNOTATIONS_PER_INTERVAL);
      numAnnotations2 = MIN (arrayMax (annotationPointers2),MAX_NUM_ANNOTATIONS_PER_INTERVAL);
      heightAnnotations = (2 * (numAnnotations1 + numAnnotations2) + 9) * SIZE_ANNOTATION; 
      trackWidth = IMAGE_WIDTH - MARGIN_LEFT - MARGIN_RIGHT - SIZE_LEFT_LABEL - SIZE_RIGHT_LABEL;
      // 6 tracks are reserved for correlation and expression plots
      trackHeight = (IMAGE_HEIGHT - MARGIN_TOP - MARGIN_BOTTOM - SIZE_TOP_AREA - heightAnnotations) / ((mio_getNumSamples () * 2) + 6); 
    
      coordinates1 = generateCoordinateSystem (annotationPointers1,numAnnotations1,trackWidth); 
      coordinates2 = generateCoordinateSystem (annotationPointers2,numAnnotations2,trackWidth); 

      blocks1 = getAlignmentBlocks (intervalPtr1);
      blocks2 = getAlignmentBlocks (intervalPtr2);
      if (intervalPtr1->strand == '+' && intervalPtr2->strand == '+') {
        ; // do nothing
      }
      else if (intervalPtr1->strand == '+' && intervalPtr2->strand == '-') {
        reverseAlignmentBlocks (blocks2);
      }
      else if (intervalPtr1->strand == '-' && intervalPtr2->strand == '+') {
        reverseAlignmentBlocks (blocks1);
        reverseAlignmentBlocks (blocks2);
      }
      else if (intervalPtr1->strand == '-' && intervalPtr2->strand == '-') {
        reverseAlignmentBlocks (blocks1);
      }
      else {
        die ("Unexpected outcome!");
      }
      assignColorsToBlocks (blocks1,colors11);
      assignColorsToBlocks (blocks2,colors12);

      iter = 1;
      drawAnnotations (im1,currMatrix,blocks1,annotationPointers1,numAnnotations1,coordinates1,trackWidth,black1,&iter);
      drawAnnotations (im1,currMatrix,blocks2,annotationPointers2,numAnnotations2,coordinates2,trackWidth,black1,&iter);
      drawAlignmentBars (im1,blocks1,blocks2,length,trackWidth,black1,&iter);
     
      iter = 1;  
      drawAnnotations (im2,currMatrix,blocks1,annotationPointers1,numAnnotations1,coordinates1,trackWidth,black2,&iter);
      drawAnnotations (im2,currMatrix,blocks2,annotationPointers2,numAnnotations2,coordinates2,trackWidth,black2,&iter);
      drawAlignmentBars (im2,blocks1,blocks2,length,trackWidth,black2,&iter);
      
      gct = gr_ct_create (0.0,(double)length,0,trackWidth);
      currStatistic1 = mio_getStatisticPointer (currMatrix,intervalNumber1);
      currStatistic2 = mio_getStatisticPointer (currMatrix,intervalNumber2);
      trackNumber = 1;
      for (j = 0; j < mio_getNumSamples (); j++) {   
        drawTrack (im1,currMatrix,rowIndicesInterval1,j,gct,currStatistic1->overallMin,currStatistic1->overallMax,trackNumber,trackHeight,heightAnnotations,red1);
        drawTrack (im2,currMatrix,rowIndicesInterval1,j,gct,MIN (currStatistic1->overallMin,currStatistic2->overallMin),MAX (currStatistic1->overallMax,currStatistic2->overallMax),trackNumber,trackHeight,heightAnnotations,red2);
        trackNumber++;
        drawTrack (im1,currMatrix,rowIndicesInterval2,j,gct,currStatistic2->overallMin,currStatistic2->overallMax,trackNumber,trackHeight,heightAnnotations,blue1);
        drawTrack (im2,currMatrix,rowIndicesInterval2,j,gct,MIN (currStatistic1->overallMin,currStatistic2->overallMin),MAX (currStatistic1->overallMax,currStatistic2->overallMax),trackNumber,trackHeight,heightAnnotations,blue2);
        trackNumber++;
      }
      
      drawCorrelationPlot (im1,currMatrix,rowIndicesInterval1,rowIndicesInterval2,heightAnnotations,gct,trackNumber,trackHeight,trackWidth,black1);
      drawExpressionPlot (im1,currMatrix,rowIndicesInterval1,rowIndicesInterval2,heightAnnotations,gct,trackNumber,trackHeight,trackWidth,black1,red1,blue1,1);
      drawBoxesAroundSampleTracks (im1,heightAnnotations,trackHeight,trackWidth,black1);
      generateSampleLabels (im1,currMatrix,heightAnnotations,trackHeight,trackWidth,black1);
      generatePairwiseCorrelationsLabels (im1,currMatrix,pairNumber1,heightAnnotations,trackHeight,trackWidth,black1);
      generateTopLabels (im1,currMatrix,pairNumber1,currStatistic1,currStatistic2,length,percentIdentity,black1,red1,blue1,1);

      drawCorrelationPlot (im2,currMatrix,rowIndicesInterval1,rowIndicesInterval2,heightAnnotations,gct,trackNumber,trackHeight,trackWidth,black2);
      drawExpressionPlot (im2,currMatrix,rowIndicesInterval1,rowIndicesInterval2,heightAnnotations,gct,trackNumber,trackHeight,trackWidth,black2,red2,blue2,0);
      drawBoxesAroundSampleTracks (im2,heightAnnotations,trackHeight,trackWidth,black2);
      generateSampleLabels (im2,currMatrix,heightAnnotations,trackHeight,trackWidth,black2);
      generatePairwiseCorrelationsLabels (im2,currMatrix,pairNumber1,heightAnnotations,trackHeight,trackWidth,black2);
      generateTopLabels (im2,currMatrix,pairNumber1,currStatistic1,currStatistic2,length,percentIdentity,black2,red2,blue2,0);
      
      stringPrintf (buffer,"%s_%d.scaled.png",currMatrix->name,pairNumber1);
      pngout = fopen (string (buffer),"wb");
      gdImagePng (im1,pngout);
      fclose (pngout);
      warn ("Done... %s",string (buffer));
      
      stringPrintf (buffer,"%s_%d.unscaled.png",currMatrix->name,pairNumber1);
      pngout = fopen (string (buffer),"wb");
      gdImagePng (im2,pngout);
      fclose (pngout);
      warn ("Done... %s",string (buffer));

      arrayDestroy (rowIndicesInterval1);
      arrayDestroy (rowIndicesInterval2);
      arrayDestroy (annotationPointers1);
      arrayDestroy (annotationPointers2);
      arrayDestroy (colors11);
      arrayDestroy (colors12);
      arrayDestroy (colors21);
      arrayDestroy (colors22);
      arrayDestroy (coordinates1);
      arrayDestroy (coordinates2);
      arrayDestroy (blocks1);
      arrayDestroy (blocks2);
      gr_ct_destroy (gct);
      gdImageDestroy (im1); 
      gdImageDestroy (im2); 
    }
  }
  mio_deInit ();
  stringDestroy (buffer);
  return 0;
}




