/*++


Module Name:

    AlignerContext.cpp

Abstract:

    Common parameters for running single & paired alignment.

Authors:

    Ravi Pandya, May, 2012

Environment:
`
    User mode service.

Revision History:

    Integrated from SingleAligner.cpp & PairedAligner.cpp

--*/

#include "stdafx.h"
#include "Compat.h"
#include "options.h"
#include "AlignerOptions.h"
#include "AlignerContext.h"
#include "AlignerStats.h"
#include "BaseAligner.h"
#include "FileFormat.h"
#include "exit.h"
#include "PairedAligner.h"
#include "Error.h"
#include "Util.h"
#include "CommandProcessor.h"

using std::max;
using std::min;

//
// Save the index & index directory globally so that we don't need to reload them on multiple runs.
//
GenomeIndex *g_index = NULL;
char *g_indexDirectory = NULL;

AlignerContext::AlignerContext(int i_argc, const char **i_argv, const char *i_version, AlignerExtension* i_extension)
    :
    index(NULL),
    writerSupplier(NULL),
    options(NULL),
    stats(NULL),
    extension(i_extension != NULL ? i_extension : new AlignerExtension()),
    readWriter(NULL),
    argc(i_argc),
    argv(i_argv),
    version(i_version),
    perfFile(NULL)
{
}

AlignerContext::~AlignerContext()
{
    delete extension;
    if (NULL != perfFile) {
        fclose(perfFile);
    }
}

void AlignerContext::runAlignment(int argc, const char **argv, const char *version, unsigned *argsConsumed)
{
    options = parseOptions(argc, argv, version, argsConsumed, isPaired());

	if (NULL == options) {	// Didn't parse correctly
		*argsConsumed = argc;
		return;
	}

#ifdef _MSC_VER
	useTimingBarrier = options->useTimingBarrier;
#endif
	
	if (!initialize()) {
		return;
	}
	_int64 alignIterStart = timeInMillis();
	WriteStatusMessage("Beginning aligner timer.\n");
    extension->initialize();
    
    if (! extension->skipAlignment()) {
        WriteStatusMessage("Aligning.\n");

        beginIteration();

        runTask();
            
        finishIteration();

        printStatsHeader();

        printStats();

        nextIteration();    // This probably should get rolled into something else; it's really cleanup code, not "next iteration"
    }

    extension->finishAlignment();
    PrintBigAllocProfile();
    PrintWaitProfile();
    _int64 alignIterTime = timeInMillis() - alignIterStart;
    WriteStatusMessage("Aligning took: %d seconds %.2f minutes", alignIterTime / 1000, alignIterTime / (1000 * 60));
}

    void
AlignerContext::initializeThread()
{
    stats = newStats(); // separate copy per thread
    stats->extra = extension->extraStats();
    readWriter = writerSupplier != NULL ? writerSupplier->getWriter() : NULL;
    extension = extension->copy();
}

    void
AlignerContext::runThread()
{
    extension->beginThread();
    runIterationThread();
    if (readWriter != NULL) {
        readWriter->close();
        delete readWriter;
    }
    extension->finishThread();
}
    
    void
AlignerContext::finishThread(AlignerContext* common)
{
    common->stats->add(stats);
    delete stats;
    stats = NULL;
    delete extension;
    extension = NULL;
}

    bool
AlignerContext::initialize()
{
    if (g_indexDirectory == NULL || strcmp(g_indexDirectory, options->indexDir) != 0) {
        delete g_index;
        g_index = NULL;
        delete g_indexDirectory;
        g_indexDirectory = new char [strlen(options->indexDir) + 1];
        strcpy(g_indexDirectory, options->indexDir);

        if (strcmp(options->indexDir, "-") != 0) {
            WriteStatusMessage("Loading index from directory... ");
 
            fflush(stdout);
            _int64 loadStart = timeInMillis();
            index = GenomeIndex::loadFromDirectory((char*) options->indexDir, options->mapIndex, options->prefetchIndex);
            if (index == NULL) {
                WriteErrorMessage("Index load failed, aborting.\n");
				return false;
            }
            g_index = index;

            _int64 loadTime = timeInMillis() - loadStart;
             WriteStatusMessage("%llds.  %u bases, seed size %d\n",
                    loadTime / 1000, index->getGenome()->getCountOfBases(), index->getSeedLength());
         } else {
            WriteStatusMessage("no alignment, input/output only\n");
        }
    } else {
        index = g_index;
    }

    maxHits_ = options->maxHits;
    maxDist_ = options->maxDist;
    extraSearchDepth = options->extraSearchDepth;
    noUkkonen = options->noUkkonen;
    noOrderedEvaluation = options->noOrderedEvaluation;
	noTruncation = options->noTruncation;
    maxSecondaryAlignmentAdditionalEditDistance = options->maxSecondaryAlignmentAdditionalEditDistance;
	maxSecondaryAlignments = options->maxSecondaryAlignments;
	minReadLength = options->minReadLength;

	if (index != NULL && (int)minReadLength < index->getSeedLength()) {
		WriteErrorMessage("The min read length (%d) must be at least the seed length (%d), or there's no hope of aligning reads that short.\n", minReadLength, index->getSeedLength());
		return false;
	}

    if (options->perfFileName != NULL) {
        perfFile = fopen(options->perfFileName,"a");
        if (NULL == perfFile) {
            WriteErrorMessage("Unable to open perf file '%s'\n", options->perfFileName);
			return false;
        }
    }

    DataSupplier::ThreadCount = options->numThreads;

    return true;
}

    void
AlignerContext::printStatsHeader()
{
	WriteStatusMessage("Total Reads    Aligned, MAPQ >= %2d    Aligned, MAPQ < %2d     Unaligned              Too Short/Too Many Ns  %%Pairs\tReads/s   Time in Aligner (s)\n", MAPQ_LIMIT_FOR_SINGLE_HIT, MAPQ_LIMIT_FOR_SINGLE_HIT);
}

    void
AlignerContext::beginIteration()
{
    writerSupplier = NULL;
    alignStart = timeInMillis();
    clipping = options->clipping;
    totalThreads = options->numThreads;
    bindToProcessors = options->bindToProcessors;
    maxDist = maxDist_;
    maxHits = maxHits_;
    numSeedsFromCommandLine = options->numSeedsFromCommandLine;
    seedCoverage = options->seedCoverage;
    minWeightToCheck = options->minWeightToCheck;
    if (stats != NULL) {
        delete stats;
    }
    stats = newStats();
    stats->extra = extension->extraStats();
    extension->beginIteration();
    
    memset(&readerContext, 0, sizeof(readerContext));
    readerContext.clipping = options->clipping;
    readerContext.defaultReadGroup = options->defaultReadGroup;
    readerContext.genome = index != NULL ? index->getGenome() : NULL;
    readerContext.ignoreSecondaryAlignments = options->ignoreSecondaryAlignments;
    readerContext.ignoreSupplementaryAlignments = options->ignoreSecondaryAlignments;   // Maybe we should split them out
    strncpy(readerContext.junctionSeq, options->junctionSeq, MAX_JUNCTION_TRIM);
    DataSupplier::ExpansionFactor = options->expansionFactor;

    typeSpecificBeginIteration();

    if (UnknownFileType != options->outputFile.fileType) {
        const FileFormat* format;
        if (SAMFile == options->outputFile.fileType) {
            format = FileFormat::SAM[options->useM];
        } else if (BAMFile == options->outputFile.fileType) {
            format = FileFormat::BAM[options->useM];
        } else {
            //
            // This shouldn't happen, because the command line parser should catch it.  Perhaps you've added a new output file format and just
            // forgoten to add it here.
            //
            WriteErrorMessage("AlignerContext::beginIteration(): unknown file type %d for '%s'\n", options->outputFile.fileType, options->outputFile.fileName);
            soft_exit(1);
        }
        format->setupReaderContext(options, &readerContext);

        writerSupplier = format->getWriterSupplier(options, readerContext.genome);
        ReadWriter* headerWriter = writerSupplier->getWriter();
        headerWriter->writeHeader(readerContext, options->sortOutput, argc, argv, version, options->rgLineContents, options->outputFile.omitSQLines);
        headerWriter->close();
        delete headerWriter;
    }
}

    void
AlignerContext::finishIteration()
{
    extension->finishIteration();

    if (NULL != writerSupplier) {
        writerSupplier->close();
        delete writerSupplier;
        writerSupplier = NULL;
    }

    alignTime = /*timeInMillis() - alignStart -- use the time from ParallelTask.h, that may exclude memory allocation time*/ time;
}

    bool
AlignerContext::nextIteration()
{
    //
    // This thing is a vestage of when we used to allow parameter ranges.
    //
    typeSpecificNextIteration();
    return false;
}

extern char *FormatUIntWithCommas(_uint64 val, char *outputBuffer, size_t outputBufferSize);	// Relying on the one in Util.h results in an "internal compiler error" for Visual Studio.

//
// Take an integer and a percentage, and turn it into a string of the form "number (percentage%)<padding>" where
// number has commas and the whole thing is padded out with spaces to a specific length.
//
char *numPctAndPad(char *buffer, _uint64 num, double pct, size_t desiredWidth, size_t bufferLen)
{
	_ASSERT(desiredWidth < bufferLen);	// < to leave room for trailing null.

	FormatUIntWithCommas(num, buffer, bufferLen);
	const size_t percentageBufferSize = 100;	// Plenty big enough for any value
	char percentageBuffer[percentageBufferSize];

	sprintf(percentageBuffer, " (%.02f%%)", pct);
	if (strlen(percentageBuffer) + strlen(buffer) >= bufferLen) {
		WriteErrorMessage("numPctAndPad: overflowed output buffer\n");
		soft_exit(1);
	}

	strcat(buffer, percentageBuffer);
	for (size_t x = strlen(buffer); x < desiredWidth; x++) {
		strcat(buffer, " ");
	}

	return buffer;
}

    void
AlignerContext::printStats()
{
    double usefulReads = max((double) stats->usefulReads, 1.0);

	const size_t strBufLen = 50;	// Way more than enough for 64 bit numbers with commas
	char tooShort[strBufLen];
	char single[strBufLen];
	char multi[strBufLen];
	char unaligned[strBufLen];
	char numReads[strBufLen];
	char readsPerSecond[strBufLen];
	char alignTimeString[strBufLen];

	WriteStatusMessage("%-14s %s %s %s %s %.02f%%%\t%-9s %s\n",
		FormatUIntWithCommas(stats->totalReads, numReads, strBufLen),
		numPctAndPad(single, stats->singleHits, 100.0 * stats->singleHits / stats->totalReads, 22, strBufLen),
		numPctAndPad(multi, stats->multiHits, 100.0 * stats->multiHits / stats->totalReads, 22, strBufLen),
		numPctAndPad(unaligned, stats->notFound, 100.0 * stats->notFound / stats->totalReads, 22, strBufLen),
		numPctAndPad(tooShort, stats->totalReads - stats->usefulReads, 100.0 * (stats->totalReads - stats->usefulReads) / max(stats->totalReads, (_int64)1), 22, strBufLen),
		100.0 * stats->alignedAsPairs / stats->totalReads,
		FormatUIntWithCommas((unsigned _int64)(1000 * stats->totalReads / max(alignTime, (_int64)1)), readsPerSecond, strBufLen),	// Aligntime is in ms
		FormatUIntWithCommas((alignTime + 500) / 1000, alignTimeString, strBufLen)
		);

    if (NULL != perfFile) {
        fprintf(perfFile, "%d\t%d\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%lld\t%lld\tt%.0f\n",
                maxHits_, maxDist_, 
                100.0 * usefulReads / max(stats->totalReads, (_int64) 1),
				100.0 * stats->singleHits / stats->totalReads,
				100.0 * stats->multiHits / stats->totalReads,
				100.0 * stats->notFound / stats->totalReads,
                stats->lvCalls,
				100.0 * stats->alignedAsPairs / stats->totalReads,
                stats->totalReads,
                (1000.0 * usefulReads) / max(alignTime, (_int64) 1));

        fprintf(perfFile,"\n");
    }


#if TIME_HISTOGRAM
    WriteStatusMessage("Per-read alignment time histogram:\nlog2(ns)\tcount\ttotal time (ns)\n");
    for (int i = 0; i < 31; i++) {
        WriteStatusMessage("%d\t%lld\t%lld\n", i, stats->countByTimeBucket[i], stats->nanosByTimeBucket[i]);
    }
#endif // TIME_HISTOGRAM


    stats->printHistograms(stdout);

#ifdef  TIME_STRING_DISTANCE
    WriteStatusMessage("%llds, %lld calls in BSD noneClose, not -1\n",  stats->nanosTimeInBSD[0][1]/1000000000, stats->BSDCounts[0][1]);
    WriteStatusMessage("%llds, %lld calls in BSD noneClose, -1\n",      stats->nanosTimeInBSD[0][0]/1000000000, stats->BSDCounts[0][0]);
    WriteStatusMessage("%llds, %lld calls in BSD close, not -1\n",      stats->nanosTimeInBSD[1][1]/1000000000, stats->BSDCounts[1][1]);
    WriteStatusMessage("%llds, %lld calls in BSD close, -1\n",          stats->nanosTimeInBSD[1][0]/1000000000, stats->BSDCounts[1][0]);
    WriteStatusMessage("%llds, %lld calls in Hamming\n",                stats->hammingNanos/1000000000,         stats->hammingCount);
#endif  // TIME_STRING_DISTANCE

    extension->printStats();
}



        AlignerOptions*
AlignerContext::parseOptions(
    int i_argc,
    const char **i_argv,
    const char *i_version,
    unsigned *argsConsumed,
    bool      paired)
{
    argc = i_argc;
    argv = i_argv;
    version = i_version;

    AlignerOptions *options;

    if (paired) {
        options = new PairedAlignerOptions("snap paired <index-dir> <inputFile(s)> [<options>] where <input file(s)> is a list of files to process.\n");
    } else {
        options = new AlignerOptions("snap single <index-dir> <inputFile(s)> [<options>] where <input file(s)> is a list of files to process.\n");
    }

    options->extra = extension->extraOptions();
    if (argc < 3) {
        WriteErrorMessage("Too few parameters\n");
        options->usage();
		delete options;
		return NULL;
    }

    options->indexDir = argv[1];
    struct InputList {
        SNAPFile    input;
        InputList*  next;
    } *inputList = NULL;

    //
    // Now build the input array and parse options.
    //

    bool inputFromStdio = false;

    int i;
    int nInputs = 0;
    for (i = 2; i < argc; i++) {	// Starting at 2 skips single/paired and the index

        if (',' == argv[i][0]  && '\0' == argv[i][1]) {
            i++;    // Consume the comma
            break;
        }

        int argsConsumed;
        SNAPFile input;
        if (SNAPFile::generateFromCommandLine(argv+i, argc-i, &argsConsumed, &input, paired, true)) {
            if (input.isStdio) {
				if (CommandPipe != NULL) {
					WriteErrorMessage("You may not use stdin/stdout in daemon mode\n");
					delete options;
					return NULL;
				}

                if (inputFromStdio) {
                    WriteErrorMessage("You specified stdin ('-') specified for more than one input, which isn't permitted.\n");
					delete options;
					return NULL;
                } else {
                    inputFromStdio = true;
                }
            }

            InputList *listEntry = new InputList;
            listEntry->input = input;
            listEntry->next = inputList;
            inputList = listEntry;      // Yes, this puts them in backwards.  a) We reverse them at the end and b) it doesn't matter anyway

            nInputs++;
            i += argsConsumed - 1;
            continue;
        }

        bool done;
        int oldI = i;
        if (!options->parse(argv, argc, i, &done)) {
            WriteErrorMessage("Didn't understand options starting at %s\n", argv[oldI]);
            options->usage();
			delete options;
			return NULL;
        }

        if (done) {
            i++;    // For the ',' arg
            break;
        }
    }

    if (0 == nInputs) {
        WriteErrorMessage("No input files specified.\n");
		delete options;
		return NULL;
    }

    if (options->maxDist + options->extraSearchDepth >= MAX_K) {
        WriteErrorMessage("You specified too large of a maximum edit distance combined with extra search depth.  The must add up to less than %d.\n", MAX_K);
        WriteErrorMessage("Either reduce their sum, or change MAX_K in LandauVishkin.h and recompile.\n");
		delete options;
		return NULL;
    }

    if (options->maxSecondaryAlignmentAdditionalEditDistance > (int)options->extraSearchDepth) {
        WriteErrorMessage("You can't have the max edit distance for secondary alignments (-om) be bigger than the max search depth (-D)\n");
		delete options;
		return NULL;
    }

    options->nInputs = nInputs;
    options->inputs = new SNAPFile[nInputs];
    for (int j = nInputs - 1; j >= 0; j --) {
        // The loop runs backwards so that we reverse the reversing that we did when we built it.  Not that it matters anyway.
        _ASSERT(NULL != inputList);
        options->inputs[j] = inputList->input;
        InputList *dying = inputList;
        inputList = inputList->next;
        delete dying;
    }
    _ASSERT(NULL == inputList);

    *argsConsumed = i;
    return options;
}
