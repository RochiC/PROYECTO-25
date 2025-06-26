import { PrismaClient } from './generated/prisma/index.js'

const prisma = new PrismaClient()

export const allMolechules = async(req,res)=>{
    let molechules = await prisma.molecula.findMany({})
    res.send(molechules)
}