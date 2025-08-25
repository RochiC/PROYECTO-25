// ai/httpClient.js
// Con Node 22 podés usar fetch global. Si preferís node-fetch, descomentá la 1ª línea:
// import fetch from "node-fetch";

export class IAHttpClient {
  constructor({ url, timeoutMs = 15000 }) {
    if (!url) throw new Error("Falta IA_URL");
    this.url = url;
    this.timeoutMs = timeoutMs;
  }

  async inferir({ prompt = null, smiles = null }) {
    // Timeout con AbortController
    const ctrl = new AbortController();
    const timer = setTimeout(() => ctrl.abort(), this.timeoutMs);

    try {
      const res = await fetch(this.url, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ prompt, smiles }),
        signal: ctrl.signal,
      });

      if (!res.ok) {
        const txt = await res.text().catch(() => "");
        throw new Error(`IA ${res.status}: ${txt}`);
      }

      const j = await res.json();

      // Normalización estricta al contrato que usa tu back
      return {
        smiles: String(j.smiles ?? "").trim(),
        peso_molecular: Number(j.peso_molecular),
        prediccion_bioactiva: Number(j.prediccion_bioactiva),
        lipinski_ok: Boolean(j.lipinski_ok),
        toxicidad_potencial: String(j.toxicidad_potencial ?? "").trim(),
      };
    } finally {
      clearTimeout(timer);
    }
  }
}

