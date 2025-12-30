from fastapi import FastAPI, UploadFile, Form
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Annotated
from mapper.readMapper.readMapper import ReadMapper, AReadMapper
from mapper.models.readMapper import ReadMapperInput, ReadMapperOutput

class ReadMapperApiInput(BaseModel):
    # files to uplod
    fastOne: UploadFile
    fastTwo: UploadFile
    referenceGenome: UploadFile

app = FastAPI()
readMapper : AReadMapper = ReadMapper()


@app.get("/")
async def root():
    return {"message": "Test"}

@app.post("/mapper/mapReads")
async def mapReads(
    fastOne: UploadFile,
    fastTwo: UploadFile,
    referenceGenome: UploadFile
    ):
    outputFileName : str = "SAMOutputFile.SAM"
    outputFolderLocation: str = "mapper/io/outputs/"
    outputFileLocation : str = outputFolderLocation + outputFileName
    readMapperInput = ReadMapperInput(
        readsOne        = fastOne.file,
        readsTwo        = fastTwo.file,
        referenceGenome = referenceGenome.file,
        outputLocation  = outputFileLocation,
        kmerSize=15,
        windowSize=30,
    ) 

    readMapperOutput : ReadMapperOutput = readMapper.mapReads(readMapperInput)

    headers = {
        "Mapped-Reads-Count": str(readMapperOutput.numberOfMappedReads),
        "Access-Control-Expose-Headers": "Mapped-Reads-Count" 
    }

    return FileResponse(
        path=outputFileLocation,
        filename=outputFileName,
        headers = headers
    )