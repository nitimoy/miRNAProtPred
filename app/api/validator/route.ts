import { NextRequest, NextResponse } from 'next/server'
import { query, MiRNARow } from '@/lib/db'
import { transcriptionTranslation } from '@/lib/ncbi'

export async function POST(request: NextRequest) {
  try {
    const { sequence, miRNA_ids } = await request.json()

    if (!sequence || !miRNA_ids) {
      return NextResponse.json(
        { error: 'Sequence and miRNA IDs are required' },
        { status: 400 }
      )
    }

    // Clean and process sequence
    const cleanedSequence = sequence.replace(/\s+/g, '').toUpperCase()

    // Convert to DNA sequence
    let dnaSequence: string
    try {
      dnaSequence = await transcriptionTranslation(cleanedSequence)
    } catch (error: any) {
      return NextResponse.json({ error: error.message }, { status: 400 })
    }

    // Parse miRNA IDs
    const ids = miRNA_ids
      .split(',')
      .map((id: string) => id.trim())
      .map((id: string) => (id.startsWith('hsa-') ? id : `hsa-${id}`))

    // Validate each miRNA ID
    const results = []

    for (const miRNAId of ids) {
      const result = await query('SELECT * FROM mirdb WHERE mir_id = $1', [miRNAId])

      if (result.rows.length === 0) {
        results.push({
          id: miRNAId,
          found: false,
        })
      } else {
        const miRNA: MiRNARow = result.rows[0]
        const seeds = [miRNA.seed_1, miRNA.seed_2, miRNA.seed_3].filter(Boolean)

        // Check if any seed matches in the sequence
        const seedMatch = seeds.some((seed) => dnaSequence.includes(seed))

        results.push({
          id: miRNAId,
          found: true,
          seed_match: seedMatch,
        })
      }
    }

    return NextResponse.json({ results })
  } catch (error) {
    console.error('Error in validator:', error)
    return NextResponse.json(
      { error: 'An error occurred while processing your request' },
      { status: 500 }
    )
  }
}
