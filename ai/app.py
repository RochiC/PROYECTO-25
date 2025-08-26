# app.py — Servicio FastAPI de IA (mock generativo "decente")
from fastapi import FastAPI
from pydantic import BaseModel
import random, re, hashlib

app = FastAPI()

class InferInput(BaseModel):
    smiles: str | None = None
    prompt: str | None = None

# Fragmentos "químicos" básicos para armar mutaciones simpáticas
FRAGS = [
    "C", "N", "O", "S", "F", "Cl", "Br",
    "CC", "CN", "OC", "C=O", "C#N", "C=C",
    "c1ccccc1",       # fenilo
    "C(C)C",          # isopropilo
    "CCO", "NC=O"
]

TWO_LETTER = {"Cl","Br"}  # para no romper halógenos al cortar

def tokenize_simple(s: str) -> list[str]:
    """Tokenizador muy simple que respeta Cl/Br como tokens de 2 letras."""
    out, i = [], 0
    while i < len(s):
        if i+1 < len(s) and s[i:i+2] in TWO_LETTER:
            out.append(s[i:i+2]); i += 2
        else:
            out.append(s[i]); i += 1
    return out

def det_rand(smiles: str) -> random.Random:
    """RNG determinístico por SMILES (para propiedades estables)."""
    h = hashlib.sha1(smiles.encode("utf-8")).hexdigest()
    seed = int(h[:16], 16)
    return random.Random(seed)

def generar_mutacion(base: str) -> str:
    """Aplica 2–3 mutaciones al SMILES base sin devolver el mismo texto."""
    if not base or len(base) < 1:
        base = "C"

    toks = tokenize_simple(base)
    rng = random.Random()  # mutación sí debe ser aleatoria en cada request
    n_ops = rng.randint(2, 3)

    for _ in range(n_ops):
        if not toks:
            toks = tokenize_simple("C")

        op = rng.choice(["insert", "replace", "delete"])

        if op == "insert":
            frag = rng.choice(FRAGS)
            pos = rng.randrange(0, len(toks)+1)
            toks[pos:pos] = tokenize_simple(frag)

        elif op == "replace" and len(toks) > 0:
            frag = rng.choice(FRAGS)
            pos = rng.randrange(0, len(toks))
            toks[pos:pos+1] = tokenize_simple(frag)

        elif op == "delete" and len(toks) > 2:
            i = rng.randrange(0, len(toks))
            j = min(len(toks), i + rng.randint(1, 2))
            del toks[i:j]

        # acotar tamaño
        if len(toks) > 60:
            toks = toks[:60]

    out = "".join(toks).strip()
    # evitar devolver igual
    if out == base or len(out) < 2:
        out = base + rng.choice(["C","N","O","Cl","Br"])
    return out

def props_mock(smiles: str) -> dict:
    """Propiedades plausibles (mock) derivadas de un RNG determinístico por SMILES."""
    rng = det_rand(smiles)

    # Peso entre 120 y 520 aprox.
    peso = round(120 + rng.random() * 400, 2)

    # Predicción bioactiva entre 0..1
    p_bio = round(rng.random(), 3)

    # Conteos simples para 'regla' de Lipinski aproximada
    n_N = len(re.findall(r"N", smiles))
    n_O = len(re.findall(r"O", smiles))
    n_hal = len(re.findall(r"Cl|Br|F", smiles))

    # Heurística tosca de aceptores/donores
    hbd = min(5, n_N)               # donores ~ N
    hba = min(10, n_N + n_O)        # aceptores ~ N+O

    # Lipinski (mock): peso < 500 y límites HBD/HBA
    lip_ok = (peso < 500) and (hbd <= 5) and (hba <= 10)

    # Toxicidad a ojímetro (más halógenos/peso -> sube)
    tox_score = (n_hal * 0.7) + (0 if peso < 300 else 0.8) + (0.5 if not lip_ok else 0)
    if tox_score < 0.8:
        tox = "baja"
    elif tox_score < 1.6:
        tox = "media"
    else:
        tox = "alta"

    return {
        "peso_molecular": peso,
        "prediccion_bioactiva": p_bio,
        "lipinski_ok": lip_ok,
        "toxicidad_potencial": tox,
    }

@app.post("/infer")
def infer(body: InferInput):
    base = (body.smiles or body.prompt or "C").strip()
    new_smiles = generar_mutacion(base)
    props = props_mock(new_smiles)
    return {"smiles": new_smiles, **props}
@app.get("/")
def root():
    return {"mensaje": "IA viva"}
