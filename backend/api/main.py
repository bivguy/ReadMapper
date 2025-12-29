from fastapi import FastAPI
from pydantic import BaseModel

app = FastAPI()



@app.get("/")
async def root():
    return {"message": "Test"}

@app.post("/mapper/mapReads")
async def mapReads():
    return {"status": "Temp"}