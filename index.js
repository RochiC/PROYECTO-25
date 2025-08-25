import 'dotenv/config';
import express from "express";
import cors from "cors";
import {
  allMolechules,
  createMolechules,
  suggestNewMolecule,
  suggestRaw,          // 👈 import del nuevo endpoint
} from "./molecula_controller.js";
import { IAHttpClient } from "./ai/httpClient.js";

const app = express();
const PORT = process.env.PORT || 3000;

app.use(express.json());
app.use(cors());

// ===== Inyección del cliente de IA =====
const IA_URL = process.env.IA_URL || "http://localhost:8000/infer"; 
const ia = new IAHttpClient({ url: IA_URL, timeoutMs: 15000 });
app.locals.ia = ia;
// =======================================

app.get("/api", (_req, res) => {
  res.json({ mensaje: "¡Hola desde el back-end!" });
});

// Rutas existentes
app.get("/all", async (req, res) => { await allMolechules(req, res); });
app.post("/create", async (req, res) => { await createMolechules(req, res); });

// Ruta que genera y guarda (si querés usar DB)
app.post("/suggest", async (req, res) => { await suggestNewMolecule(req, res); });

// NUEVA ruta: llamada directa a la IA (sin DB)
app.post("/suggest-raw", async (req, res) => { await suggestRaw(req, res); });

app.listen(PORT, () => {
  console.log(`Servidor escuchando en http://localhost:${PORT}`);
});
