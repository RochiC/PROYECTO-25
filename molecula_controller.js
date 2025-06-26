import { PrismaClient } from './generated/prisma/index.js'

const prisma = new PrismaClient()

export const allMolechules = async(req,res)=>{
    try {
        let molechules = await prisma.molecula.findMany({})
        res.send(molechules)

    } catch (error) {
        res.status(500).send(error)
    }
}

export const createMolechules = async(req,res)=>{
    const {lipinski_ok,peso_molecular,prediccion_bioactiva,smiles,toxicidad_potencial} = req.body
   try {
    let molechules = await prisma.molecula.create({
        data:{
            lipinski_ok,
            peso_molecular,
            prediccion_bioactiva,
            smiles,
            toxicidad_potencial,
        }
    })
    res.status(201).send("Creado correctamente")
   } catch (error) {
    res.status(500).send(error)
   }
}