"""Explicit design-preview endpoint, independent of discovery activity ranking."""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, ConfigDict, Field
from starlette.concurrency import run_in_threadpool

router = APIRouter(prefix='/api/motor', tags=['ALK5 design'])


class AnalogueDesignRequest(BaseModel):
    model_config = ConfigDict(extra='forbid')
    count: int = Field(default=30, ge=1, le=100, strict=True)


@router.post('/alk5/analogues')
async def alk5_analogues(body: AnalogueDesignRequest):
    from app.services.alk5_analogue_design import design_analogues
    try:
        return await run_in_threadpool(design_analogues, body.count)
    except (ImportError, FileNotFoundError) as exc:
        raise HTTPException(status_code=503, detail='ALK5 design requires RDKit and the curated reference asset') from exc


async def predict_profiles(smiles):
    from app.services.admet_client import predict_admet
    return await predict_admet(smiles)


@router.post('/alk5/analogues/admet')
async def alk5_analogue_admet(body: AnalogueDesignRequest):
    from app.services.alk5_design_admet import compare_profiles, profile_inputs
    designs = await alk5_analogues(body)
    predictions = await predict_profiles(profile_inputs(designs))
    try:
        return compare_profiles(designs, predictions)
    except (ValueError, KeyError, TypeError) as exc:
        raise HTTPException(status_code=502, detail='Ambiguous or invalid ADMET service response') from exc
