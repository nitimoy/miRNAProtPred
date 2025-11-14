import { Pool } from 'pg'

let pool: Pool | null = null

export function getPool() {
  if (!pool) {
    pool = new Pool({
      user: process.env.DB_USERNAME,
      password: process.env.DB_PASSWORD,
      host: process.env.DB_HOST,
      port: parseInt(process.env.DB_PORT || '5432'),
      database: process.env.DB_NAME,
    })
  }
  return pool
}

export async function query(text: string, params?: any[]) {
  const pool = getPool()
  const result = await pool.query(text, params)
  return result
}

export interface Visit {
  id: number
  ip_address: string
  city: string
  country: string
  visit_count: number
  visit_time: Date
}

export interface MiRNARow {
  r_id: number
  Description: string
  HumanmiRNAID: string
  Accession: string
  Sequence: string
  seed_1: string
  seed_2: string
  seed_3: string
}
