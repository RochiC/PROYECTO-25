import { PrismaClient } from "./generated/prisma/index.js";

const prisma = new PrismaClient();

export const allMolechules = async (req, res) => {
  try {
    const data = await prisma.molecula.findMany({});
    return res.json({ ok: true, data, mensaje: "ok" });
  } catch (error) {
    console.error(error);
    return res.status(500).json({ ok: false, mensaje: "error interno" });
  }
};

export const createMolechules = async (req, res) => {
  try {
    let {
      lipinski_ok,
      peso_molecular,
      prediccion_bioactiva,
      smiles,
      toxicidad_potencial,
      prompt,
    } = req.body ?? {};

    const faltan =
      lipinski_ok === undefined ||
      peso_molecular === undefined ||
      prediccion_bioactiva === undefined ||
      !toxicidad_potencial ||
      !smiles;

    if (faltan) {
      const ia = req.app.locals.ia;
      const out = await ia.inferir({ prompt, smiles });
      smiles = out.smiles || smiles;
      lipinski_ok = out.lipinski_ok;
      peso_molecular = out.peso_molecular;
      prediccion_bioactiva = out.prediccion_bioactiva;
      toxicidad_potencial = out.toxicidad_potencial;
    }

    if (!smiles) {
      return res.status(422).json({ ok: false, mensaje: "SMILES faltante" });
    }

    const existe = await prisma.molecula.findFirst({ where: { smiles } });
    if (existe) {
      return res.status(400).json({ ok: false, mensaje: "La molécula ya existe" });
    }

    const created = await prisma.molecula.create({
      data: {
        smiles,
        peso_molecular: parseFloat(peso_molecular),
        prediccion_bioactiva: parseFloat(prediccion_bioactiva),
        lipinski_ok: Boolean(lipinski_ok),
        toxicidad_potencial: String(toxicidad_potencial),
      },
    });

    return res.status(201).json({
      ok: true,
      data: created,
      mensaje: "Molécula creada correctamente",
    });
  } catch (error) {
    console.error(error);
    return res.status(500).json({ ok: false, mensaje: error.message ?? "error interno" });
  }
};

export const suggestNewMolecule = async (req, res) => {
  try {
    const { base_smiles, save = false, max_attempts = 5 } = req.body ?? {};
    if (!base_smiles || typeof base_smiles !== "string") {
      return res.status(400).json({ ok: false, mensaje: "Falta base_smiles" });
    }

    const ia = req.app.locals.ia;
    let generated = null;
    let attempts = 0;

    while (attempts < Number(max_attempts)) {
      attempts += 1;
      const out = await ia.inferir({ smiles: base_smiles });
      const newSmiles = String(out.smiles || "").trim();
      if (!newSmiles || newSmiles === base_smiles) continue;

      const dup = await prisma.molecula.findFirst({ where: { smiles: newSmiles } });
      if (dup) continue;

      generated = {
        smiles: newSmiles,
        peso_molecular: Number(out.peso_molecular),
        prediccion_bioactiva: Number(out.prediccion_bioactiva),
        lipinski_ok: Boolean(out.lipinski_ok),
        toxicidad_potencial: String(out.toxicidad_potencial),
      };
      break;
    }

    if (!generated) {
      return res.status(409).json({
        ok: false,
        mensaje: `No se pudo generar una molécula nueva y única en ${max_attempts} intentos`,
      });
    }

    if (!save) {
      return res.json({
        ok: true,
        data: generated,
        mensaje: "Molécula generada (no guardada)",
      });
    }

    const created = await prisma.molecula.create({ data: generated });
    return res.status(201).json({
      ok: true,
      data: created,
      mensaje: "Molécula generada y guardada",
    });
  } catch (error) {
    console.error(error);
    return res.status(500).json({ ok: false, mensaje: error.message ?? "error interno" });
  }
};

// --- NUEVO: endpoint simple que conecta directo con la IA (NO usa DB) ---
export const suggestRaw = async (req, res) => {
  try {
    const { smiles = null, prompt = null } = req.body ?? {};
    if (!smiles && !prompt) {
      return res.status(400).json({ ok: false, mensaje: "Mandá smiles o prompt" });
    }
    const ia = req.app.locals.ia;
    const out = await ia.inferir({ smiles, prompt });
    return res.json({ ok: true, data: out, mensaje: "OK (sin DB)" });
  } catch (e) {
    console.error("suggestRaw error:", e);
    return res.status(500).json({ ok: false, mensaje: e.message || "error llamando a IA" });
  }
};
