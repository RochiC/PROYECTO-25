import 'dotenv/config';
import express from "express";
import cors from "cors";

const app = express();
const PORT = process.env.PORT || 3000;

app.use(express.json());
app.use(cors());

// ---- función para mutar un SMILES de forma más variada ----
function mutarSmiles(s) {
  if (!s || typeof s !== "string") return s;

  const frags = ["N", "O", "Cl", "Br", "C=O", "C", "S"];
  const chars = s.split("");
  const op = Math.floor(Math.random() * 3); // 0=insertar, 1=reemplazar, 2=borrar

  if (op === 0) {
    // insertar fragmento
    const frag = frags[Math.floor(Math.random() * frags.length)];
    const pos = Math.floor(Math.random() * (chars.length + 1));
    chars.splice(pos, 0, frag);
  } else if (op === 1 && chars.length > 0) {
    // reemplazar
    const frag = frags[Math.floor(Math.random() * frags.length)];
    const pos = Math.floor(Math.random() * chars.length);
    chars[pos] = frag;
  } else if (op === 2 && chars.length > 1) {
    // borrar
    const pos = Math.floor(Math.random() * chars.length);
    chars.splice(pos, 1);
  }

  return chars.join("");
}

// ---- endpoint que devuelve mutación + propiedades falsas ----
app.post("/mutar", (req, res) => {
  const { smiles } = req.body ?? {};
  if (!smiles) {
    return res.status(400).json({ ok: false, mensaje: "Falta 'smiles' en el body" });
  }

  const nuevo = mutarSmiles(smiles);

  // Valores mock (falsos pero variados)
  const peso = (100 + Math.random() * 400).toFixed(2);
  const prediccion = Math.random().toFixed(3);
  const lipinski = Math.random() > 0.3; // true o false al azar
  const tox = ["baja", "media", "alta"][Math.floor(Math.random() * 3)];

  return res.json({
    ok: true,
    data: {
      smiles: nuevo,
      peso_molecular: Number(peso),
      prediccion_bioactiva: Number(prediccion),
      lipinski_ok: lipinski,
      toxicidad_potencial: tox
    },
    mensaje: "mutación mock aplicada"
  });
});

app.get("/api", (_req, res) => res.json({ mensaje: "¡Hola desde el back-end!" }));

app.listen(PORT, () => {
  console.log(`Servidor escuchando en http://localhost:${PORT}`);
});
