# app.py — Servicio FastAPI para generar SMILES nuevos
from fastapi import FastAPI
from pydantic import BaseModel
import os, random

# Intentar usar modelo de Hugging Face o local si existe
USE_MODEL = False
tokenizer = None
model = None

HUGGINGFACE_ID = os.getenv("HUGGINGFACE_ID", "").strip()
LOCAL_DIR = os.getenv("LOCAL_MODEL_DIR", "").strip()  # opcional: ruta absoluta o relativa

def _try_load_model():
    global USE_MODEL, tokenizer, model
    try:
        if HUGGINGFACE_ID:
            from transformers import AutoTokenizer, AutoModelForCausalLM
            tokenizer = AutoTokenizer.from_pretrained(HUGGINGFACE_ID)
            model = AutoModelForCausalLM.from_pretrained(HUGGINGFACE_ID)
            model.eval()
            USE_MODEL = True
            print(f"[IA] Cargado desde HuggingFace: {HUGGINGFACE_ID}")
            return
    except Exception as e:
        print(f"[IA] No pude cargar {HUGGINGFACE_ID}: {e}")

    try:
        if LOCAL_DIR and os.path.isdir(LOCAL_DIR):
            from transformers import AutoTokenizer, AutoModelForCausalLM
            tokenizer = AutoTokenizer.from_pretrained(LOCAL_DIR)
            model = AutoModelForCausalLM.from_pretrained(LOCAL_DIR)
            model.eval()
            USE_MODEL = True
            print(f"[IA] Cargado modelo local: {LOCAL_DIR}")
            return
    except Exception as e:
        print(f"[IA] No pude cargar modelo local {LOCAL_DIR}: {e}")

    # sin modelo → fallback
    USE_MODEL = False
    print("[IA] Sin modelo. Uso fallback de mutación de SMILES.")

_try_load_model()

# RDKit opcional para validar/calcular props (si no está instalado, seguimos)
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski
    RDKit_OK = True
except Exception:
    RDKit_OK = False

app = FastAPI()

class InferInput(BaseModel):
    smiles: str | None = None
    prompt: str | None = None

def mutate_smiles(smiles: str) -> str:
    """
    Fallback: cambia el SMILES de forma no trivial.
    - inserta, reemplaza o borra un símbolo común.
    """
    if not smiles:
        smiles = "C"
    tokens = list(smiles)
    ops = ["insert", "replace", "delete"]
    op = random.choice(ops)

    # Algunos símbolos “orgánicos” frecuentes en SMILES:
    pool = ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "B", "=", "#", "(", ")", "[O-]", "[NH+]", "c", "n", "o"]

    if op == "insert":
        t = random.choice(pool)
        pos = random.randrange(0, len(tokens)+1)
        if len(t) == 2:  # ej Cl, Br
            tokens[pos:pos] = list(t)
        else:
            tokens.insert(pos, t)
    elif op == "replace" and len(tokens) > 0:
        t = random.choice(pool)
        pos = random.randrange(0, len(tokens))
        if len(t) == 1:
            tokens[pos] = t
        else:
            # reemplazo por varios chars
            tokens[pos:pos+1] = list(t)
    elif op == "delete" and len(tokens) > 1:
        pos = random.randrange(0, len(tokens))
        del tokens[pos]

    candidate = "".join(tokens)
    candidate = candidate.replace("Cl", "Cl").replace("Br", "Br")  # no-op para mantener dígrafos
    if candidate == smiles:
        return smiles + "O"  # asegurar cambio visible
    return candidate

def generar_con_modelo(base: str, max_len: int = 60) -> str:
    """
    Si hay modelo causal cargado, genera desde el SMILES base.
    """
    from transformers import StoppingCriteria, StoppingCriteriaList
    input_ids = tokenizer.encode(base, return_tensors="pt")
    out_ids = model.generate(
        input_ids,
        max_length=max_len,
        do_sample=True,
        top_k=50,
        top_p=0.95,
        temperature=1.0,
        pad_token_id=tokenizer.eos_token_id if tokenizer.eos_token_id else None,
    )
    return tokenizer.decode(out_ids[0], skip_special_tokens=True).strip()

def es_valido(smiles: str) -> bool:
    if not RDKit_OK:
        return bool(smiles)
    return Chem.MolFromSmiles(smiles) is not None

def props_desde_smiles(smiles: str) -> dict:
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
    lip_ok = (mw < 500) and (logp <= 5) and (hba <= 10) and (hbd <= 5)
    tox = "baja" if (logp <= 3.5 and mw < 600) else "media"
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

    if USE_MODEL:
        # intentar algunas veces que sea distinto y válido
        new_smiles = None
        for _ in range(6):
            cand = generar_con_modelo(base)
            if cand and cand != base and es_valido(cand):
                new_smiles = cand
                break
        if not new_smiles:
            # si el modelo no dio algo útil, uso fallback
            new_smiles = mutate_smiles(base)
    else:
        new_smiles = mutate_smiles(base)

    props = props_desde_smiles(new_smiles)
    return {"smiles": new_smiles, **props}
