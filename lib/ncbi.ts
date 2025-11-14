import axios from 'axios'

const ENTREZ_BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
const BLAST_BASE_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'

export interface BlastResult {
  ncbi_gi: string
  ncbi_id: string
  start_pos: number
  end_pos: number
}

/**
 * Perform BLAST search using NCBI's BLAST API
 * Note: This is a simplified version. For production, consider using a more robust library
 * or implementing the full BLAST workflow with polling
 */
export async function performBlast(sequence: string): Promise<BlastResult> {
  try {
    // Submit BLAST request
    const submitResponse = await axios.get(BLAST_BASE_URL, {
      params: {
        CMD: 'Put',
        PROGRAM: 'tblastn',
        DATABASE: 'nt',
        QUERY: sequence,
        FORMAT_TYPE: 'XML',
      },
    })

    // Extract RID (Request ID)
    const ridMatch = submitResponse.data.match(/RID = (\w+)/)
    if (!ridMatch) {
      throw new Error('Failed to submit BLAST request')
    }
    const rid = ridMatch[1]

    // Poll for results (simplified - in production, add proper polling with delays)
    await new Promise((resolve) => setTimeout(resolve, 30000)) // Wait 30 seconds

    const resultsResponse = await axios.get(BLAST_BASE_URL, {
      params: {
        CMD: 'Get',
        RID: rid,
        FORMAT_TYPE: 'XML',
      },
    })

    // Parse XML results (simplified - in production, use a proper XML parser)
    // This is a placeholder - actual parsing would be more complex
    const xmlData = resultsResponse.data

    // For now, return a placeholder structure
    // In a real implementation, you'd parse the XML properly
    return {
      ncbi_gi: 'placeholder_gi',
      ncbi_id: 'placeholder_id',
      start_pos: 1,
      end_pos: 100,
    }
  } catch (error) {
    console.error('BLAST error:', error)
    throw new Error('BLAST search failed')
  }
}

/**
 * Retrieve sequence from NCBI using Entrez
 */
export async function retrieveSequence(
  ncbiId: string,
  startPos: number,
  endPos: number
): Promise<string> {
  try {
    const response = await axios.get(`${ENTREZ_BASE_URL}/efetch.fcgi`, {
      params: {
        db: 'nucleotide',
        id: ncbiId,
        rettype: 'fasta',
        retmode: 'text',
        seq_start: startPos,
        seq_stop: endPos,
      },
    })

    // Extract sequence from FASTA format
    const lines = response.data.split('\n')
    const sequence = lines.slice(1).join('').replace(/\s/g, '')
    return sequence
  } catch (error) {
    console.error('Entrez fetch error:', error)
    throw new Error('Failed to retrieve sequence from NCBI')
  }
}

/**
 * Detect input format and convert to DNA sequence
 */
export async function transcriptionTranslation(sequence: string): Promise<string> {
  // Detect format
  const isDNA = /^[ATGC]+$/i.test(sequence)
  const isRNA = /^[AUGC]+$/i.test(sequence)
  const isProtein = /^[ARDNCEQGHILKMFPSTWYV-]+$/i.test(sequence)

  if (isDNA) {
    return sequence.toLowerCase()
  } else if (isRNA) {
    return sequence.replace(/U/gi, 'T').toLowerCase()
  } else if (isProtein) {
    // For protein sequences, we would need to perform BLAST
    // This is complex and time-consuming, so for now we'll throw an error
    // In production, you'd implement the full BLAST workflow
    throw new Error(
      'Protein sequence input requires BLAST search, which is not fully implemented in this version. Please provide DNA or RNA sequence.'
    )
  } else {
    throw new Error('Invalid input sequence format')
  }
}
