from fastapi import FastAPI, UploadFile, Form
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
    # apiInput: Annotated[ReadMapperApiInput, Form()]
    ):
    readMapperInput = ReadMapperInput(
        readsOne        = fastOne.file,
        readsTwo        = fastTwo.file,
        referenceGenome = referenceGenome.file
    )

    ReadMapperOutput = readMapper.mapReads(readMapperInput)
    return ReadMapperOutput