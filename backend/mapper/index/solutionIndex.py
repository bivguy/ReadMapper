from typing import IO, List

from dataclasses import dataclass, field
from extend.extender import Alignment


@dataclass
class SolutionIndex:
    start:      int 
    end:        int

@dataclass 
class Metrics:
    FP:         int 
    TP:         int 
    FN:         int 
    TN:         int 

    Precision:  float
    Recall:     float
    Accuracy:   float


class MetricAccumulator:
    def __init__(self, solution_map: dict):
        self.solution_map = solution_map
        self.tp = 0
        self.fp = 0
        self.total_mapped_reads = 0
        self.seen_reads = set() # Track IDs to calculate FN/TN later
        self.correct_reads = set() # Track IDs that were TPs

    def update(self, batch_alignments: List[Alignment]):
        for a in batch_alignments:
            self.seen_reads.add(a.readId)
            
            if not a.mapped:
                continue
                
            self.total_mapped_reads += 1
            true_mapping : SolutionIndex = self.solution_map.get(a.readId)

            if true_mapping:
                # Read is in ground truth
                start_diff = abs(a.ref_start - true_mapping.start)
                end_diff = abs(a.ref_end - true_mapping.end)
                
                # Check threshold (using 5 as per your previous code)
                if start_diff > 5 or end_diff > 5:
                    self.fp += 1 # Mapped, but to wrong place
                else:
                    self.tp += 1 # Mapped to correct place
                    self.correct_reads.add(a.readId)
            else:
                # Read is NOT in ground truth, but we mapped it -> False Positive
                self.fp += 1

    def compute_final_metrics(self, total_reads_processed: int) -> Metrics:
        # FN: Reads in truth that were NOT True Positives
        # (This covers unmapped truth reads AND truth reads mapped to wrong place)
        fn = 0
        for read_id in self.solution_map.keys():
            if read_id not in self.correct_reads:
                fn += 1
        
        # TN: Reads NOT in truth that were NOT mapped
        # Total Negatives in Input = (Total Reads) - (Reads in Truth)
        # TN = (Total Negatives in Input) - (False Positives on Negatives)
        # However, easier logic: 
        # TN = (Total Reads Processed) - (TP + FP + FN) 
        # (Note: This assumes every read is exclusively TP, FP, FN, or TN)
        
        # True Total Negative Reads (Reads that shouldn't map)
        true_total_negatives = total_reads_processed - len(self.solution_map)
        
        # Reads that shouldn't map, but we mapped them anywa
        
        # Simplest valid TN calc:
        tn = true_total_negatives - (self.fp - (len(self.seen_reads) - len(self.correct_reads) - fn))
        
        tn = total_reads_processed - self.tp - self.fp - fn
        if tn < 0: tn = 0

        precision = self.tp / (self.tp + self.fp) if (self.tp + self.fp) > 0 else 0.0
        recall = self.tp / (self.tp + fn) if (self.tp + fn) > 0 else 0.0
        accuracy = (self.tp + tn) / (self.tp + tn + self.fp + fn) if (self.tp + tn + self.fp + fn) > 0 else 0.0

        return Metrics(FP=self.fp, TP=self.tp, FN=fn, TN=tn, 
                      Precision=precision, Recall=recall, Accuracy=accuracy)



class SolutionIndexBuilder:
    def __init__(self):
        # create a dictionary mapping the ground truth readID to its start and end position
        pass
    
    def getSolutionMap(self, solutionFile: IO) -> dict:
        solutionMap = {}
        while True:
            item  = self.getItem(solutionFile)
            if not item:
                break
            solutionMap[item[0]] = SolutionIndex(start = int(item[1]), end= int(item[2]))
        return solutionMap
    
    def getItem(sefl, solutionFile: IO)-> List:
        item = []
        line = solutionFile.readline().rstrip('\n')
        if not line:
            return []
        item = line.split("\t")
        return item

    def computeMetrics(self, alignments: List[Alignment], solutionMap: dict) -> Metrics:
        fp = tp = fn = tn = correct = 0
        trueTotalNegativeReads = len(alignments) - len(solutionMap)
        
        seenReads = set()          # reads for which we produced ANY mapping
        correctForRead = set()     # reads we counted as TP (avoid double TPs if duplicates)
        negMappedReads = set()     # reads that are negatives in truth but we mapped anyway
        
        for a in alignments:
            if not a.mapped:
                continue
            seenReads.add(a.readId)
            trueMapping: SolutionIndex = solutionMap.get(a.readId)
            # read was mapped
            if trueMapping:
                startDifference = abs(a.ref_start - trueMapping.start)
                endDifference = abs(a.ref_end - trueMapping.end)

                # false positive
                if startDifference > 5 or endDifference > 5:
                    fp += 1
                # true positive
                else:
                    correctForRead.add(a.readId)
                    tp += 1

            # we thought this was a valid read mapping when ti was not
            else:
                 # truth says this read is negative, but we mapped it -> FP on negatives
                fp += 1
                negMappedReads.add(a.readId)

        
        # FN: truth-positive reads we never mapped at all (a was None or absent)
        for readId in solutionMap.keys():
            if readId not in seenReads and readId not in correctForRead:
                fn += 1

        # TN: negatives in truth that we also left unmapped
        tn = trueTotalNegativeReads - len(negMappedReads)
        if tn < 0:
            tn = 0  # safety guard if inputs ever violate assumptions

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        accuracy = (tp + tn) / (tp + tn + fp + fn) 

        return Metrics(
            FP=fp,
            TP=tp,
            FN=fn,
            TN=tn,
            Precision=precision,
            Recall=recall,
            Accuracy=accuracy
        )
        


                
            

        
