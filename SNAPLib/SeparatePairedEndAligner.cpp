/*++

  Module Name:

  SeparatePairedEndAligner.cpp

  Abstract:

  A paired-end aligner calls into a different paired-end aligner, and if
  it fails to find an alignment, aligns each of the reads singly.  This handles
  chimeric reads that would otherwise be unalignable.

  Authors:

  Bill Bolosky, June, 2013

  Environment:

  User mode service.

  Revision History:

  --*/


#include "stdafx.h"
#include "SeparatePairedEndAligner.h"
#include "mapq.h"
#include "directions.h"
#include "BigAlloc.h"
#include "Util.h"

using namespace std;

#ifdef TRACE_PAIRED_ALIGNER
#define TRACE printf
#else
#define TRACE(...) {}
#endif

SeparatePairedEndAligner::SeparatePairedEndAligner(
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
   bool noTruncation,
   unsigned	minReadLength_,
   BigAllocator        *allocator)
: ChimericPairedEndAligner(index_, maxReadSize, maxHits, maxK, maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, forceSpacing_, extraSearchDepth, noUkkonen, noOrderedEvaluation, noTruncation, NULL, minReadLength_, allocator)
{
  //  std::cout << "SeparatePairedEndAligner();" << std::endl;
}

size_t 
SeparatePairedEndAligner::getBigAllocatorReservation(
					     GenomeIndex *   index, 
						     unsigned        maxReadSize, 
						     unsigned        maxHits, 
						     unsigned        seedLen, 
						     unsigned        maxSeedsFromCommandLine, 
						     double          seedCoverage,
						     unsigned        maxEditDistanceToConsider, 
						     unsigned        maxExtraSearchDepth, 
					     unsigned        maxCandidatePoolSize,
					     int             maxSecondaryAlignmentsPerContig)
{
  return BaseAligner::getBigAllocatorReservation(false, maxHits, maxReadSize, seedLen, maxSeedsFromCommandLine, seedCoverage) + sizeof(SeparatePairedEndAligner) + sizeof(_uint64);
}


SeparatePairedEndAligner::~SeparatePairedEndAligner()
{
  // All handled by parent
}

#ifdef _DEBUG
extern bool _DumpAlignments;
#endif // _DEBUG


void SeparatePairedEndAligner::align(
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
				     )
{
  result->status[0] = result->status[1] = NotFound;
  *nSecondaryResults = 0;
  *nSingleEndSecondaryResultsForFirstRead = 0;
  *nSingleEndSecondaryResultsForSecondRead = 0;

  if (read0->getDataLength() < minReadLength && read1->getDataLength() < minReadLength) {
    TRACE("Reads are both too short -- returning");
    for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
      result->location[whichRead] = 0;
      result->mapq[whichRead] = 0;
      result->score[whichRead] = 0;
      result->status[whichRead] = NotFound;
    }
    result->alignedAsPair = false;
    result->fromAlignTogether = false;
    result->nanosInAlignTogether = 0;
    result->nLVCalls = 0;
    result->nSmallHits = 0;
    return;
  }

  _int64 start = timeInNanos();
  Read *read[NUM_READS_PER_PAIR] = {read0, read1};
  int *resultCount[2] = {nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead};
  bool allOk = true;
  SingleAlignmentResult singleResult[2];
  int tried_aligning[2] = {0,0};
  for (int r = 0; r < NUM_READS_PER_PAIR; r++) {

    int singleEndSecondaryResultsThisTime = 0;
    if (read[r]->getDataLength() <= minReadLength) {
      result->status[r] = NotFound;
      result->mapq[r] = 0;
      result->direction[r] = FORWARD;
      result->location[r] = 0;
      result->score[r] = 0;
    } else {
      tried_aligning[r] = 1;
    // We're using *nSingleEndSecondaryResultsForFirstRead because it's either 0 or what all we've seen (i.e., we know NUM_READS_PER_PAIR is 2)
      singleAligner->AlignRead(read[r], &singleResult[r], 
			       maxEditDistanceForSecondaryResults,
			       singleSecondaryBufferSize - *nSingleEndSecondaryResultsForFirstRead, 
			       &singleEndSecondaryResultsThisTime,
			       singleEndSecondaryResults + *nSingleEndSecondaryResultsForFirstRead);
      *(resultCount[r]) = singleEndSecondaryResultsThisTime;
      if (singleResult[r].location != InvalidGenomeLocation) {
	result->status[r] = singleResult[r].status;
	result->mapq[r] = singleResult[r].mapq; // no penalty for chimeric reads with this aligner
	result->direction[r] = singleResult[r].direction;
	result->location[r] = singleResult[r].location;
	result->score[r] = singleResult[r].score;
      }
      else {
      	result->status[r] = NotFound;
      	result->mapq[r] = 0;
      	result->direction[r] = FORWARD;
      	result->location[r] = 0;
      	result->score[r] = 0;
      }
    }
  }
  if (result->status[0] == NotFound || result->status[1] == NotFound) {
    result->fromAlignTogether = false;
    result->alignedAsPair = true;
  }
  else {
    result->fromAlignTogether = true;
    result->alignedAsPair = true;
  }
  for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
    if (result->mapq[r] > 1000) {
      //std::cout << "Bad mapq?" << std::endl;
      printf("bad mapq - SeparatePairedEndAligner: (%u, %u) score (%d, %d), MAPQ (%d, %d) datalength (%d, %d) tried aligning (%d, %d)\n",result->location[0], result->location[1],
	     result->score[0], result->score[1], result->mapq[0], result->mapq[1], read0->getDataLength(), read1->getDataLength(), tried_aligning[0], tried_aligning[1]);
      //      assert(0);
    }
  }


#ifdef _DEBUG
  if (_DumpAlignments) {
    printf("SeparatePairedEndAligner: (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n",result->location[0], result->location[1],
	   result->score[0], result->score[1], result->mapq[0], result->mapq[1]);
  }
#endif // _DEBUG
                    
}
