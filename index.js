import 'dotenv/config';
import express from "express";
import cors from "cors";

const app = express();
const PORT = process.env.PORT || 3000;

app.use(express.json());
app.use(cors());

// ---- función ultra simple para "cambiar" el SMILES ----
function mutarSmiles(s) {
  if (!s || typeof s !== "string") return s;

  // Si ya termina en 'N', agrego 'O'; si no, agrego 'N'.
  // Así garantizamos que sea distinto.
  return s.endsWith("N") ? s + "O" : s + "N";
}

// Endpoint mínimo: recibe {smiles} y devuelve {smiles:nuevo}
app.post("/mutar", (req, res) => {
  const { smiles } = req.body ?? {};
  if (!smiles) {
    return res.status(400).json({ ok: false, mensaje: "Falta 'smiles' en el body" });
  }
  const nuevo = mutarSmiles(smiles);
  return res.json({
    ok: true,
    data: { smiles: nuevo },
    mensaje: "mutación simple aplicada",
  });
});

app.get("/api", (_req, res) => res.json({ mensaje: "¡Hola desde el back-end!" }));

app.listen(PORT, () => {
  console.log(`Servidor escuchando en http://localhost:${PORT}`);
});
