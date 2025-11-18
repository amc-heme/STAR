#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"
#include "AlignVsTranscript.h"
#include "ReadAnnotations.h"
#include <bitset>

void Transcriptome::classifyScrapsAlign (Transcript **alignG, uint64 nAlignG, ReadAnnotations &readAnnot) 
{
    // Only checks if the 3' end of the alignment overlaps an exon (strand-specific)
    // Positive strand (Str==0): 3' end is the rightmost position (end of last exon)
    // Negative strand (Str==1): 3' end is the leftmost position (start of first exon)
    // For paired-end reads, only considers mate 1 (iFrag==0)

    ReadAnnotFeature &annFeat = readAnnot.annotFeatures[SoloFeatureTypes::Scraps];

    annFeat.fAlign.resize(nAlignG);    
       
    for (uint iag=0; iag<nAlignG; iag++) {
        
        Transcript &aG=*alignG[iag];
        
        // Skip mate 2 in paired-end reads (only process mate 1, iFrag==0)
        if (aG.iFrag > 0)
            continue;

        // Determine the 3' end position based on strand
        uint64 pos3p;
        if (aG.Str == 0) {
            // Positive strand: 3' end is the rightmost position (end of last exon)
            pos3p = aG.exons[aG.nExons-1][EX_G] + aG.exons[aG.nExons-1][EX_L] - 1;
        } else {
            // Negative strand: 3' end is the leftmost position (start of first exon)
            pos3p = aG.exons[0][EX_G];
        }

        // Binary search through transcript starts to find transcripts that might contain this position
        uint32 tr1=binarySearch1a<uint>(pos3p, trS, nTr);
        if (tr1==(uint32) -1) 
            continue; // This position is before all transcripts

        ++tr1;
        do {
            --tr1;
            
            // Check if this transcript contains the 3' position
            if (pos3p < trS[tr1] || pos3p > trE[tr1])
                continue;
            
            // Check strand compatibility
            if (P.pSolo.strand >= 0 && (trStr[tr1]==1 ? aG.Str : 1-aG.Str) != (uint32)P.pSolo.strand)
                continue;
            
            // Now check if the 3' position overlaps an exon in this transcript
            bool overlapsExon = false;
            uint32 exI = trExI[tr1]; // Starting index for this transcript's exons
            for (uint32 iex = 0; iex < trExN[tr1]; iex++) {
                uint32 exStart = exSE[2*(exI + iex)];
                uint32 exEnd = exSE[2*(exI + iex) + 1];
                
                if (pos3p >= exStart && pos3p <= exEnd) {
                    overlapsExon = true;
                    break;
                }
            }
            
            if (overlapsExon) {
                // The 3' end overlaps an exon in this transcript
                annFeat.fSet.insert(trGene[tr1]);
                annFeat.fAlign[iag].insert(trGene[tr1]);
            }
            
        } while (trEmax[tr1] >= pos3p && tr1 > 0);
    }
    
    if (annFeat.fSet.size() > 0)
        annFeat.ovType = ReadAnnotFeature::overlapTypes::exonic;
}
