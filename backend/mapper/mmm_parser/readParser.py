from ..mmm_parser.parser import Parser
from ..models.read import FullRead
from ..models.read import Read
from typing import List


class ReadParser:
    def __init__(self, parserFront: Parser, parserBack: Parser):
        self.parserFront = parserFront
        self.parserBack = parserBack 

    def parseFullRead(self) -> FullRead:
        frontRead : Read = self.parserFront.parseNextRead()
        backRead : Read = self.parserBack.parseNextRead()

        # if either parser is done, stop
        if frontRead is None or backRead is None:
            return None

        return FullRead(frontRead, backRead)

    def parseReadPair(self) -> List[Read]:
        frontRead : Read = self.parserFront.parseNextRead()
        backRead : Read = self.parserBack.parseNextRead()

        if frontRead is None or backRead is None:
            return []
        return [frontRead, backRead]

    def parseAllReads(self) -> List[FullRead]:
        reads: List[FullRead] = []

        while True:
            fullRead: FullRead = self.parseFullRead()
            if not fullRead:
                break
            reads.append(fullRead)

        return reads

    def parseAllReadPairs(self) -> List[List[FullRead]]:
        readPairs: List[List[FullRead]] = []

        while True:
            readPair: List[FullRead] = self.parseReadPair()
            if not readPair:
                break
            readPairs.append(readPair)

        return readPairs
