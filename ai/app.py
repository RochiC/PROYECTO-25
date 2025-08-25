# app.py  — Servicio FastAPI de IA que genera SMILES nuevos con tu modelo
from fastapi import FastAPI
from pydantic import BaseModel
from transformers import AutoTokenizer, AutoModelForCausalLM
import torch

# RDKit es opcional: si no está instalado, seguimos sin validación/propiedades
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski
    RDKit_OK = True
except Exception:
    RDKit_OK = False

app = FastAPI()

# === Carga tu modelo entrenado ===
# Asegurate que la carpeta ./CHEMBert_modelo exista al lado de este archivo.
tokenizer = AutoTokenizer.from_pretrained("./CHEMBert_modelo")
model = AutoModelForCausalLM.from_pretrained("./CHEMBert_modelo")
model.eval()

class InferInput(BaseModel):
    smiles: str | None = None   # SMILES base para condicionar la generación
    prompt: str | None = None   # opcional

def es_valido(smiles: str) -> bool:
    if not RDKit_OK:
        # Sin RDKit no podemos validar de verdad: damos por válido
        return bool(smiles) and len(smiles) > 0
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def generar_molecula(smiles_base: str = "C", max_len: int = 40) -> str:
    """
    Usa tu modelo causal para generar un SMILES nuevo
    condicionado por el SMILES base. Reintenta algunas veces
    hasta que sea distinto/valido si RDKit está disponible.
    """
    input_ids = tokenizer.encode(smiles_base, return_tensors="pt")

    # Probamos algunos intentos por si sale inválido o igual al base
    for _ in range(8):
        outputs = model.generate(
            input_ids,
            max_length=max_len,
            do_sample=True,
            top_k=50,
            top_p=0.95,
            temperature=1.0,
            pad_token_id=tokenizer.eos_token_id
            if tokenizer.eos_token_id is not None else None,
        )
        new_smiles = tokenizer.decode(outputs[0], skip_special_tokens=True).strip()
        if new_smiles and new_smiles != smiles_base and es_valido(new_smiles):
            return new_smiles

    # Si no hubo suerte, devolvemos aunque sea el último (o el base como fallback)
    return new_smiles if new_smiles else smiles_base

def props_desde_smiles(smiles: str) -> dict:
    """
    Calcula propiedades básicas (si RDKit está). Si no, devuelve None.
    """
    if not RDKit_OK:
        return {
            "peso_molecular": None,
            "prediccion_bioactiva": None,
            "lipinski_ok": None,
            "toxicidad_potencial": None,
        }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "peso_molecular": None,
            "prediccion_bioactiva": None,
            "lipinski_ok": None,
            "toxicidad_potencial": None,
        }

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hba = Lipinski.NumHAcceptors(mol)
    hbd = Lipinski.NumHDonors(mol)

    # Regla de Lipinski (sencilla)
    lip_ok = (mw < 500) and (logp <= 5) and (hba <= 10) and (hbd <= 5)

    # Una “toxicidad” muy básica a modo demo
    tox = "baja" if (logp <= 3.5 and mw < 600) else "media"

    # “predicción bioactiva” placeholder 0..1 (demo basada en logP)
    pred_bio = max(0.0, min(1.0, 1.0 - (abs(logp) / 6.0)))

    return {
        "peso_molecular": round(float(mw), 2),
        "prediccion_bioactiva": round(float(pred_bio), 3),
        "lipinski_ok": bool(lip_ok),
        "toxicidad_potencial": tox,
    }

@app.post("/infer")
def infer(body: InferInput):
    base = (body.smiles or body.prompt or "C").strip()
    new_smiles = generar_molecula(base)

    # Propiedades (si RDKit está, se calculan de verdad; si no, vienen como None)
    props = props_desde_smiles(new_smiles)

    return {
        "smiles": new_smiles,
        **props
    }
