import express from "express"
const app = express();
const PORT = 3000;
app.use(express.json());

app.get('/api', (req, res) => {
    res.json({ mensaje: '¡Hola desde el back-end!' });
    });

app.listen(PORT, () => {
    console.log(`Servidor escuchando en http://localhost:${PORT}`);
      });


