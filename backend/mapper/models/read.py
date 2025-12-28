class Read:
    def __init__(self, identifier: str, sequence: str, qualityScore: str, isFront: bool):
        self.identifier = identifier
        self.sequence = sequence
        self.qualityScore = qualityScore
        self.isFront = isFront
    
    def getIdentifier(self) -> str:
        return self.identifier

    def getSequence(self) -> str:
        return self.sequence

    def getQualityScore(self) -> str:
        return self.qualityScore    

class FullRead:
    def __init__(self, frontRead: Read, backRead: Read):
        self.frontRead = frontRead
        self.backRead = backRead
    
    def getFullRead(self) -> str:
        return self.getFrontReadString() + self.getBackReadString()

    def getReadID(self) -> str:
        return self.frontRead.getIdentifier()

    def getQualityScore(self) -> str:
        # TODO: back read quality score might be the backwards order
        return self.frontRead.getQualityScore() + self.backRead.getQualityScore()

    def getFrontReadString(self) -> str:
        return self.frontRead.getSequence()
    
    def getBackReadString(self) -> str:
        return self.backRead.getSequence()