from dataclasses import dataclass
from typing import IO

@dataclass 
class ReadMapperInput:
    # required fields
    readsOne: IO 
    readsTwo: IO 
    referenceGenome: IO 

    outputLocation: str

    # optional fields 
    groundTruth: IO = None
    kmerSize: int = 15 
    windowSize: int = 30

@dataclass 
class ReadMapperOutput:
    samOutput: IO
    numberOfMappedReads: int = -1


