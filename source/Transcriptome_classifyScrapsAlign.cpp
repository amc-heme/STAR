#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"
#include "AlignVsTranscript.h"
#include "ReadAnnotations.h"
#include <bitset>

void Transcriptome::classifyScrapsAlign (Transcript **alignG, uint64 nAlignG, ReadAnnotations &readAnnot) 
{
    // Like Gene feature, but only requires the 5' end of the alignment to overlap an exon
    // Positive strand (Str==0): 5' end is the leftmost position (start of first aligned block)
    // Negative strand (Str==1): 5' end is the rightmost position (end of last aligned block)

    ReadAnnotFeature &annFeat = readAnnot.annotFeatures[SoloFeatureTypes::Scraps];

    annFeat.fAlign.resize(nAlignG);    
       
    for (uint iag=0; iag<nAlignG; iag++) {
        
        Transcript &aG=*alignG[iag];

        // Determine the 5' end position based on strand
        uint64 pos5p;
        if (aG.Str == 0) {
            // Positive strand: 5' end is the leftmost position (start of first aligned block)
            pos5p = aG.exons[0][EX_G];
        } else {
            // Negative strand: 5' end is the rightmost position (end of last aligned block)
            pos5p = aG.exons[aG.nExons-1][EX_G] + aG.exons[aG.nExons-1][EX_L] - 1;
        }

        // Binary search through transcript starts using the 5' position
        uint32 tr1=binarySearch1a<uint>(pos5p, trS, nTr);
        if (tr1==(uint32) -1) 
            continue; // The 5' position is before all transcripts

        ++tr1;
        do {
            --tr1;
            
            // Check strand compatibility (if strand-specific mode is enabled)
            if (P.pSolo.strand >= 0 && (trStr[tr1]==1 ? aG.Str : 1-aG.Str) != (uint32)P.pSolo.strand)
                continue;
            
            // Check if the 5' position is within this transcript's range
            if (pos5p < trS[tr1] || pos5p > trE[tr1])
                continue;
            
            // Check if the 5' position overlaps an exon in this transcript
            bool overlapsExon = false;
            uint32 exI = trExI[tr1]; // Starting index for this transcript's exons in the exSE array
            for (uint32 iex = 0; iex < trExN[tr1]; iex++) {
                uint32 exStart = exSE[2*(exI + iex)];     // Exon start position
                uint32 exEnd = exSE[2*(exI + iex) + 1];   // Exon end position
                
                if (pos5p >= exStart && pos5p <= exEnd) {
                    // The 5' end of the read overlaps this exon
                    overlapsExon = true;
                    break;
                }
            }
            
            if (overlapsExon) {
                annFeat.fSet.insert(trGene[tr1]);
                annFeat.fAlign[iag].insert(trGene[tr1]);
            }
            
        } while (tr1 > 0 && trS[tr1] <= pos5p);  // Continue while transcript starts are before or at the 5' position
    }
    
    if (annFeat.fSet.size() > 0)
        annFeat.ovType = ReadAnnotFeature::overlapTypes::exonic;
}
