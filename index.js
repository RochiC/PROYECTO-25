import express from "express"
import cors from "cors"
import { allMolechules } from "./molecula_controller.js";
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

app.post("/todo",(req,res)=>{allMolechules(req,res)})

