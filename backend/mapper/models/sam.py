from typing import IO
from dataclasses import dataclass

@dataclass
class SAMInput:
    QNAME:  str     = ""            # Query template NAME
    FLAG:   int     = -1            # Bitwise FLAG
    RNAME:  str     = ""            # Reference sequence NAME
    POS:    int     = -1            # 1-based leftmost mapping position
    MAPQ:   int     = -1            # Mapping quality
    CIGAR:  str     = ""            # CIGAR string
    RNEXT:  str     = "*"           # Reference name of the mate/next read
    PNEXT:  int     = -1            # Position of the mate/next read
    TLEN:   int     = -1            # Observed template length
    SEQ:    str     = "*"           # Segment sequence
    QUAL:   str     = "*"           # ASCII of Phred-scaled base quality +33

class SAM:
    """
    A simple writer for SAM (Sequence Alignment/Map) format files.

    This class handles writing SAM alignment records to an output file.
    It enforces the correct field order for the 11 required fields
    (QNAME through QUAL) and fills in missing fields with '*'.
    """
    def __init__(self, referenceName: str, referenceSize: int, outputFile: IO):
        """
        Initialize a SAM writer.

        Parameters
        ----------
        outputFile : IO
            An open writable file handle (e.g., from `open("out.sam", "w")`)
            where SAM records will be written.
        """     
        self.outputFile = outputFile
        # generate the SAM header

        # TODO: figure out what we need to put in the header
        self.outputFile.write("@HD VN:1.7 SO:unsorted\n")
        self.outputFile.write(f"@SQ\tSN:{referenceName}\tLN:{referenceSize}\n")

    def WriteReadToSam(self, input: SAMInput):
        """
        Write a single alignment record to the SAM file.

        Parameters
        ----------
        samMappings : dict
            A dictionary mapping SAM field names to values.
            Expected keys are the 11 required SAM fields:
            QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT,
            PNEXT, TLEN, SEQ, QUAL.
            
            Missing fields will be written as '*'.
            Extra fields (optional tags like "NM:i:1") are ignored
            in this version.
        """
        output = []

        # generate a list that represents each field in the correct order
        for f in input.__dataclass_fields__:
                    output.append(str(getattr(input, f)))
        # create the line we will write to the file
        outputLine = "\t".join(output) + "\n"
        # write the line to the output file
        self.outputFile.write(outputLine)