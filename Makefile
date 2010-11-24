# required input from environment:
# $(BIOINFOCONFDIR) and $(BICOSN)

# ---------------- setup symbols needed in the rest of the file --------

DBMS = oracle
include $(BIOINFOCONFDIR)/biosdefs.make
.SUFFIXES:
SHELL = /bin/sh



ROOTFLAGS = -D_REENTRANT -pthread -m64
ROOTLIBS  = -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf \
            -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix \
            -lPhysics -lMathCore -lThread -lfreetype -pthread -lm -ldl
ROOTINC = -I$(ROOTSYS)/include


GDDIR       = /home1/lh372/SOFT/gd-2.0.35
GDINC       = -I$(GDDIR)/include
GDLIB       = -L$(GDDIR)/lib



# ----------------------- entry points --------------

PROGRAMS = substituteSequenceNames psl2interval generateMatrices mio2images mioExpressionFilter mioAnnotationFilter mioReorder mioSelectMatrix mioAddStatistics mioAddAnnotation bam2bedGraph normalizeBedGraphs filterAlignments subsetAlignments mioAddCorrelations validatePgenes mioAddScores mio2meta mioSubsetMatrices mio2bed testMio

MODULES = mio.o

CGIS= pseudoSeq_cgi showMatrix_cgi

all: allprogs


allprogs: $(MODULES) $(PROGRAMS) 


cgi: $(CGIS)
	@echo "use 'make deploy' to install the CGIs"

clean: 
	/bin/rm -f $(PROGRAMS) $(MODULES) $(CGIS)

deploy:
	rsync -a --delete $(CGIS) lh372@gw.gersteinlab.org:/nfs/web/lh372/cgi-bin




# ----------------------------------------------------

substituteSequenceNames: substituteSequenceNames.c $(BIOSLIB)
	-@/bin/rm -f substituteSequenceNames
	$(CC) $(CFLAGSO) $(BIOSINC) substituteSequenceNames.c -o substituteSequenceNames $(BIOSLNK) -lm

psl2interval: psl2interval.c $(BIOSLIB)
	-@/bin/rm -f psl2interval
	$(CC) $(CFLAGSO) $(BIOSINC) psl2interval.c -o psl2interval $(BIOSLNK) -lm

generateMatrices: generateMatrices.c mio.o $(BIOSLIB)
	-@/bin/rm -f generateMatrices
	$(CC) $(CFLAGSO) $(BIOSINC) generateMatrices.c mio.o -o generateMatrices $(BIOSLNK) -lm

mio2images: mio2images.c mio.o $(GDDIR)/include/gd.h $(GDDIR)/include/gdfontmb.h $(BIOSLIB)
	-@/bin/rm -f $O/mio2images
	$(CC) $(CFLAGSO) $(BIOSINC) $(GDINC) -I$(BIOINFOGSLDIR)/include mio2images.c mio.o -o mio2images $(BIOSLNK) $(BIOINFOGSLDIR)/lib/libgsl.a $(BIOINFOGSLDIR)/lib/libgslcblas.a $(GDLIB) -lgd -lpng -lm

mioExpressionFilter: mioExpressionFilter.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioExpressionFilter
	$(CC) $(CFLAGSO) $(BIOSINC) mioExpressionFilter.c mio.o -o mioExpressionFilter $(BIOSLNK) -lm

mioAnnotationFilter: mioAnnotationFilter.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioAnnotationFilter
	$(CC) $(CFLAGSO) $(BIOSINC) mioAnnotationFilter.c mio.o -o mioAnnotationFilter $(BIOSLNK) -lm

mioReorder: mioReorder.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioReorder 
	$(CC) $(CFLAGSO) $(BIOSINC) mioReorder.c mio.o -o mioReorder $(BIOSLNK) -lm

mioSelectMatrix: mioSelectMatrix.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioSelectMatrix 
	$(CC) $(CFLAGSO) $(BIOSINC) mioSelectMatrix.c mio.o -o mioSelectMatrix $(BIOSLNK) -lm

mioAddCorrelations: mioAddCorrelations.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioAddCorrelations
	$(CC) $(CFLAGSO) $(BIOSINC) -I$(BIOINFOGSLDIR)/include mioAddCorrelations.c mio.o -o mioAddCorrelations $(BIOSLNK) $(BIOINFOGSLDIR)/lib/libgsl.a $(BIOINFOGSLDIR)/lib/libgslcblas.a -lm

mioAddStatistics: mioAddStatistics.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioAddStatistics 
	$(CC) $(CFLAGSO) $(BIOSINC) mioAddStatistics.c mio.o -o mioAddStatistics $(BIOSLNK) -lm

mioAddAnnotation: mioAddAnnotation.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioAddAnnotation 
	$(CC) $(CFLAGSO) $(BIOSINC) mioAddAnnotation.c mio.o -o mioAddAnnotation $(BIOSLNK) -lm

bam2bedGraph: bam2bedGraph.c $(BIOSLIB)
	-@/bin/rm -f bam2bedGraph
	$(CC) $(CFLAGSO) $(BIOSINC) bam2bedGraph.c -o bam2bedGraph $(BIOSLNK) -lm

normalizeBedGraphs: normalizeBedGraphs.c $(BIOSLIB)
	-@/bin/rm -f normalizeBedGraphs
	$(CC) $(CFLAGSO) $(BIOSINC) normalizeBedGraphs.c -o normalizeBedGraphs $(BIOSLNK) -lm

filterAlignments: filterAlignments.c $(BIOSLIB)
	-@/bin/rm -f filterAlignments
	$(CC) $(CFLAGSO) $(BIOSINC) filterAlignments.c -o filterAlignments $(BIOSLNK) -lm

subsetAlignments: subsetAlignments.c mio.o $(BIOSLIB)
	-@/bin/rm -f subsetAlignments 
	$(CC) $(CFLAGSO) $(BIOSINC) subsetAlignments.c mio.o -o subsetAlignments $(BIOSLNK) -lm

validatePgenes: validatePgenes.c $(BIOSLIB)
	-@/bin/rm -f validatePgenes
	$(CC) $(CFLAGSO) $(BIOSINC) validatePgenes.c -o validatePgenes $(BIOSLNK) -lm

mioAddScores: mioAddScores.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioAddScores 
	$(CC) $(CFLAGSO) $(BIOSINC) mioAddScores.c mio.o -o mioAddScores $(BIOSLNK) -lm

mio2meta: mio2meta.c mio.o $(BIOSLIB)
	-@/bin/rm -f mio2meta 
	$(CC) $(CFLAGSO) $(BIOSINC) mio2meta.c mio.o -o mio2meta $(BIOSLNK) -lm

mioSubsetMatrices: mioSubsetMatrices.c mio.o $(BIOSLIB)
	-@/bin/rm -f mioSubsetMatrices 
	$(CC) $(CFLAGSO) $(BIOSINC) mioSubsetMatrices.c mio.o -o mioSubsetMatrices $(BIOSLNK) -lm

mio2bed: mio2bed.c mio.o $(BIOSLIB)
	-@/bin/rm -f mio2bed 
	$(CC) $(CFLAGSO) $(BIOSINC) mio2bed.c mio.o -o mio2bed $(BIOSLNK) -lm

testMio: testMio.c mio.o $(BIOSLIB)
	-@/bin/rm -f testMio
	$(CC) $(CFLAGSO) $(BIOSINC) testMio.c mio.o -o testMio $(BIOSLNK) -lm



pseudoSeq_cgi: pseudoSeq_cgi.c mio.o $(BIOSLIB)
	-@/bin/rm -f pseudoSeq_cgi 
	$(CC) $(CFLAGSO) $(BIOSINC) pseudoSeq_cgi.c mio.o -o pseudoSeq_cgi $(BIOSLNK) -lm

showMatrix_cgi: showMatrix_cgi.c mio.o $(BIOSLIB)
	-@/bin/rm -f showMatrix_cgi 
	$(CC) $(CFLAGSO) $(BIOSINC) showMatrix_cgi.c mio.o -o showMatrix_cgi $(BIOSLNK) -lm


mio.o: mio.c mio.h $(BIOSLIB)  
	-@/bin/rm -f $O/mio.o
	$(CC) $(CFLAGSO) $(BIOSINC) mio.c -c -o mio.o





