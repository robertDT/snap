/*++

Module Name:

    SeparatePairedEndAligner.h

Abstract:

    A paired-end aligner that always aligns each read signly and doesn't penalize chimeric
    reads mapq score. Useful for mate pair libraries that are known to have large and variable
    insert sizes aligning to denovo assemblies.

Authors:

    Chuck Sugnet, December 2014

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "PairedEndAligner.h"
#include "ChimericPairedEndAligner.h"
#include "BaseAligner.h"
#include "BigAlloc.h"

class SeparatePairedEndAligner : public ChimericPairedEndAligner {
public:
    SeparatePairedEndAligner(
        GenomeIndex         *index_,
        unsigned            maxReadSize,
        unsigned            maxHits,
        unsigned            maxK,
        unsigned            maxSeedsFromCommandLine,
        double              seedCoverage,
	unsigned            minWeightToCheck,
        bool                forceSpacing_,
        unsigned            extraSearchDepth,
        bool                noUkkonen,
        bool                noOrderedEvaluation,
		bool				noTruncation,
		unsigned			minReadLength_,
        BigAllocator        *allocator);
    
    virtual ~SeparatePairedEndAligner();
    static size_t getBigAllocatorReservation(GenomeIndex * index, unsigned maxReadSize, unsigned maxHits, unsigned seedLen, unsigned maxSeedsFromCommandLine, 
                                             double seedCoverage, unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize, int             maxSecondaryAlignmentsPerContig);
    
    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(SeparatePairedEndAligner)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/* do nothing.  Memory gets cleaned up when the allocator is deleted.*/}

    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        int                    secondaryResultBufferSize,
        int                   *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
        int                    singleSecondaryBufferSize,
        int                   *nSingleEndSecondaryResultsForFirstRead,
        int                   *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        );

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    virtual _int64 getLocationsScored() const {
        return singleAligner->getLocationsScored();
    }
};
