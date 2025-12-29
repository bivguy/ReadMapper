from ..constants.constants import KMERSIZE
from typing import IO
from ..models.read import Read
# from ..models.read import Read


class Parser:
    def __init__(self, readFile: IO, referenceFile: IO = None): 
        self.readFile = readFile
        '''self.referenceFile = referenceFile
        self.referenceDescription = self.referenceFile.readline().strip('\n')'''
        self.kmerSize = KMERSIZE

        return

    def getNextReferenceKmer(self) -> str:
        res = []

        i = 0
        while i < self.kmerSize:
            char = self.referenceFile.read(1)
            if not char:
                break 
            if char == '\n' or char == '\r':
                continue 
            res.append(char)
            i += 1

        return ''.join(res)

    # def getReferenceString(self) -> str:
    #     return "".join(line.strip() for line in self.referenceFile)

    # TODO: handle both cases of the reads
    def parseNextRead(self) -> Read:
        identifier = ""
        sequence = ""
        qualityScore = ""
        isFront = True
        for i in range(4):
            line = self.readFile.readline()
            if not line:
                return None

            # removing the newline characters if applicable
            line = line.strip('\n')
            line = line.strip('\r\n')

            if i == 0:
                identifier = line[1:]
                # this is to check if the read is from the front or the back
                isFront = False if line[-1] == '2' else True
                # reversing it to left to right order if it's from the
            if i == 1:
                # if not isFront:
                #     sequence = line[::-1]
                # else:
                sequence = line
            if i == 2:
                continue
            if i == 3:
                qualityScore = line
        
        return Read(identifier=identifier, sequence=sequence, qualityScore=qualityScore, isFront=isFront)