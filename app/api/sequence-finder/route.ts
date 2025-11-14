import { NextRequest, NextResponse } from 'next/server'
import { query, MiRNARow } from '@/lib/db'
import { boyerMooreSearch } from '@/lib/algorithms'
import { transcriptionTranslation } from '@/lib/ncbi'
import ExcelJS from 'exceljs'
import { writeFile, mkdir } from 'fs/promises'
import path from 'path'

export async function POST(request: NextRequest) {
  try {
    const { sequence, job_id } = await request.json()

    if (!sequence) {
      return NextResponse.json({ error: 'Sequence is required' }, { status: 400 })
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

    // Fetch miRNA database
    const result = await query('SELECT * FROM mirdb')
    const rows: MiRNARow[] = result.rows

    // Extract target sequences from seed columns
    const targetSequences: string[] = []
    rows.forEach((row) => {
      if (row.seed_1) targetSequences.push(row.seed_1)
      if (row.seed_2) targetSequences.push(row.seed_2)
      if (row.seed_3) targetSequences.push(row.seed_3)
    })

    // Find matches
    interface Match {
      sequence: string
      occurrences: number[]
    }

    const matches: Match[] = []
    targetSequences.forEach((targetSeq) => {
      const occurrences = boyerMooreSearch(dnaSequence, targetSeq)
      if (occurrences.length > 0) {
        matches.push({ sequence: targetSeq, occurrences })
      }
    })

    if (matches.length === 0) {
      return NextResponse.json({
        no_results: true,
        message: 'No matching rows found based on your input.',
      })
    }

    // Build results
    interface ResultRow extends MiRNARow {
      Seed: string
      Position: number
    }

    const matchingRows: ResultRow[] = []
    const addedRowIds = new Set<number>()

    matches.forEach(({ sequence: targetSeq, occurrences }) => {
      const matchingMiRNAs = rows.filter(
        (row) =>
          (row.seed_1 === targetSeq ||
            row.seed_2 === targetSeq ||
            row.seed_3 === targetSeq) &&
          !addedRowIds.has(row.r_id)
      )

      matchingMiRNAs.forEach((miRNA) => {
        occurrences.forEach((position) => {
          matchingRows.push({
            ...miRNA,
            Seed: targetSeq,
            Position: position + 1,
          })
        })
        addedRowIds.add(miRNA.r_id)
      })
    })

    // Sort by position
    matchingRows.sort((a, b) => a.Position - b.Position)

    // Get top 20% for display
    const displayRows =
      matchingRows.length > 10
        ? matchingRows.slice(0, Math.floor(matchingRows.length * 0.2))
        : matchingRows

    // Generate Excel file
    const workbook = new ExcelJS.Workbook()
    const worksheet = workbook.addWorksheet('Results')

    // Add headers
    worksheet.columns = [
      { header: 'Description', key: 'Description', width: 30 },
      { header: 'Human miRNA ID', key: 'HumanmiRNAID', width: 20 },
      { header: 'Accession', key: 'Accession', width: 15 },
      { header: 'Sequence', key: 'Sequence', width: 30 },
      { header: 'Seed', key: 'Seed', width: 15 },
      { header: 'Position', key: 'Position', width: 10 },
    ]

    // Add data
    matchingRows.forEach((row) => {
      worksheet.addRow({
        Description: row.Description,
        HumanmiRNAID: row.HumanmiRNAID,
        Accession: row.Accession,
        Sequence: row.Sequence,
        Seed: row.Seed,
        Position: row.Position,
      })
    })

    // Save Excel file
    const timestamp = Date.now()
    const randomStr = Math.random().toString(36).substring(7)
    const filename = `${job_id || `matching_rows_${timestamp}_${randomStr}`}.xlsx`
    const filepath = path.join(process.cwd(), 'public', 'downloads', filename)

    // Ensure downloads directory exists
    await mkdir(path.join(process.cwd(), 'public', 'downloads'), { recursive: true })

    await workbook.xlsx.writeFile(filepath)

    return NextResponse.json({
      results: displayRows.map((row) => ({
        Description: row.Description,
        HumanmiRNAID: row.HumanmiRNAID,
        Accession: row.Accession,
        Sequence: row.Sequence,
        Seed: row.Seed,
        Position: row.Position,
      })),
      excel_path: `downloads/${filename}`,
    })
  } catch (error) {
    console.error('Error in sequence-finder:', error)
    return NextResponse.json(
      { error: 'An error occurred while processing your request' },
      { status: 500 }
    )
  }
}
