import express from "express"
import cors from "cors"
import { allMolechules, createMolechules } from "./molecula_controller.js";
const app = express();
const PORT = 3000;
app.use(express.json());
app.use(cors())

app.get('/api', (req, res) => {
    res.json({ mensaje: '¡Hola desde el back-end!' });
    });

app.listen(PORT, () => {
    console.log(`Servidor escuchando en http://localhost:${PORT}`);
      });

app.get("/all",async(req,res)=>{await allMolechules(req,res)})
app.post("/create",async(req,res)=>{await createMolechules(req,res)})

