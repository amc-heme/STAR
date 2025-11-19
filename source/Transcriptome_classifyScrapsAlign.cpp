#include "Transcriptome.h"
#include "ReadAlign.h"
#include "serviceFuns.cpp"
#include "AlignVsTranscript.h"
#include "ReadAnnotations.h"
#include <bitset>
#include <iostream>

void Transcriptome::classifyScrapsAlign (Transcript **alignG, uint64 nAlignG, ReadAnnotations &readAnnot) 
{
    // Like Gene feature, but only requires the 5' end of the alignment to overlap an exon
    // Positive strand (Str==0): 5' end is the leftmost position (start of first aligned block)
    // Negative strand (Str==1): 5' end is the rightmost position (end of last aligned block)

    ReadAnnotFeature &annFeat = readAnnot.annotFeatures[SoloFeatureTypes::Scraps];

    annFeat.fAlign.resize(nAlignG);    
    
    static int debugCount = 0;
    bool doDebug = (debugCount < 10);
    if (doDebug) debugCount++;
       
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

        if (doDebug) {
            std::cerr << "DEBUG Scraps: align " << iag << " nExons=" << aG.nExons 
                      << " Str=" << aG.Str << " pos5p=" << pos5p 
                      << " readStart=" << aG.exons[0][EX_G] 
                      << " readEnd=" << (aG.exons[aG.nExons-1][EX_G]+aG.exons[aG.nExons-1][EX_L]-1) << std::endl;
        }

        // Binary search through transcript starts using the 5' position
        uint32 tr1=binarySearch1a<uint>(pos5p, trS, nTr);
        if (tr1==(uint32) -1) {
            if (doDebug) std::cerr << "  Binary search returned -1 (before all transcripts)" << std::endl;
            continue; // The 5' position is before all transcripts
        }

        if (doDebug) std::cerr << "  Binary search returned tr1=" << tr1 << " nTr=" << nTr << std::endl;

        ++tr1;
        int transcriptsChecked = 0;
        do {
            --tr1;
            transcriptsChecked++;
            
            if (doDebug && transcriptsChecked <= 5) {
                std::cerr << "    Checking transcript " << tr1 << " trS=" << trS[tr1] 
                          << " trE=" << trE[tr1] << " trStr=" << (int)trStr[tr1] 
                          << " trGene=" << trGene[tr1] << " trExN=" << trExN[tr1] << std::endl;
            }
            
            // Check strand compatibility (if strand-specific mode is enabled)
            if (P.pSolo.strand >= 0 && (trStr[tr1]==1 ? aG.Str : 1-aG.Str) != (uint32)P.pSolo.strand) {
                if (doDebug && transcriptsChecked <= 5) std::cerr << "      Skipped: strand mismatch" << std::endl;
                continue;
            }
            
            // Check if the 5' position is within this transcript's range
            if (pos5p < trS[tr1] || pos5p > trE[tr1]) {
                if (doDebug && transcriptsChecked <= 5) std::cerr << "      Skipped: pos5p outside transcript range" << std::endl;
                continue;
            }
            
            // Check if the 5' position overlaps an exon in this transcript
            bool overlapsExon = false;
            uint32 exI = trExI[tr1]; // Starting index for this transcript's exons in the exSE array
            for (uint32 iex = 0; iex < trExN[tr1]; iex++) {
                uint32 exStart = exSE[2*(exI + iex)];     // Exon start position
                uint32 exEnd = exSE[2*(exI + iex) + 1];   // Exon end position
                
                if (doDebug && transcriptsChecked <= 5 && iex < 3) {
                    std::cerr << "        Exon " << iex << ": " << exStart << "-" << exEnd;
                }
                
                if (pos5p >= exStart && pos5p <= exEnd) {
                    // The 5' end of the read overlaps this exon
                    overlapsExon = true;
                    if (doDebug && transcriptsChecked <= 5) std::cerr << " MATCH!" << std::endl;
                    break;
                } else if (doDebug && transcriptsChecked <= 5 && iex < 3) {
                    std::cerr << " no" << std::endl;
                }
            }
            
            if (overlapsExon) {
                if (doDebug) std::cerr << "      ASSIGNED to gene " << trGene[tr1] << std::endl;
                annFeat.fSet.insert(trGene[tr1]);
                annFeat.fAlign[iag].insert(trGene[tr1]);
            }
            
        } while (tr1 > 0 && trS[tr1] <= pos5p);  // Continue while transcript starts are before or at the 5' position
        
        if (doDebug) std::cerr << "  Checked " << transcriptsChecked << " transcripts total" << std::endl;
    }
    
    if (annFeat.fSet.size() > 0)
        annFeat.ovType = ReadAnnotFeature::overlapTypes::exonic;
}
