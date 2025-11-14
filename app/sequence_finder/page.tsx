'use client'

import { useState } from 'react'
import HeroSection from '../components/HeroSection'

interface Result {
  Description: string
  HumanmiRNAID: string
  Accession: string
  Sequence: string
  Seed: string
  Position: number
}

export default function SequenceFinder() {
  const [jobId, setJobId] = useState('')
  const [sequence, setSequence] = useState('')
  const [results, setResults] = useState<Result[] | null>(null)
  const [excelPath, setExcelPath] = useState('')
  const [error, setError] = useState('')
  const [noResults, setNoResults] = useState('')
  const [isLoading, setIsLoading] = useState(false)

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    setIsLoading(true)
    setError('')
    setNoResults('')
    setResults(null)
    setExcelPath('')

    window.dispatchEvent(new Event('loadingStart'))

    try {
      const response = await fetch('/api/sequence-finder', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequence, job_id: jobId }),
      })

      const data = await response.json()

      if (!response.ok) {
        setError(data.error || 'An error occurred')
      } else if (data.no_results) {
        setNoResults(data.message)
      } else {
        setResults(data.results)
        setExcelPath(data.excel_path)
      }
    } catch (err) {
      setError('Failed to process request. Please try again.')
    } finally {
      setIsLoading(false)
      window.dispatchEvent(new Event('loadingEnd'))
    }
  }

  return (
    <>
      <HeroSection title="miRNA-SeqFinder" />

      <div className="container mt-5">
        <h6 className="mt-4">
          • This tool helps find potential human miRNA binding sites in any DNA, RNA, or protein sequence
          from pathogens.
        </h6>
        <h6 className="mt-4">
          • Users can input any sequence type, and the tool will identify matching miRNA seed sequences.
        </h6>

        <form onSubmit={handleSubmit} className="mt-4">
          <div className="form-group">
            <label htmlFor="job_id">Job Id (Optional):</label>
            <input
              type="text"
              className="form-control"
              id="job_id"
              value={jobId}
              onChange={(e) => setJobId(e.target.value)}
            />
          </div>
          <div className="form-group">
            <label htmlFor="sequence">Enter Sequence (DNA/RNA/Protein):</label>
            <textarea
              className="form-control"
              id="sequence"
              rows={4}
              required
              value={sequence}
              onChange={(e) => setSequence(e.target.value)}
            />
          </div>
          <br />
          <button type="submit" className="btn btn-primary btn-lg" disabled={isLoading} style={{ width: '100%' }}>
            {isLoading ? 'Processing...' : 'Search'}
          </button>
        </form>

        {error && (
          <div className="alert alert-danger mt-3" role="alert">
            {error}
          </div>
        )}

        {noResults && (
          <div className="alert alert-info mt-3" role="alert">
            {noResults}
          </div>
        )}

        {results && results.length > 0 && (
          <>
            <h2 className="mt-4">Matching Rows: {jobId || ''}</h2>
            <div style={{ overflowX: 'auto' }}>
              <table className="table table-striped">
                <thead>
                  <tr>
                    <th>Description</th>
                    <th>Human miRNA ID</th>
                    <th>Accession Number</th>
                    <th>Sequence 5&apos; to 3&apos;</th>
                    <th>Seed</th>
                    <th>Position</th>
                  </tr>
                </thead>
                <tbody>
                  {results.map((row, index) => (
                    <tr key={index}>
                      <td>{row.Description}</td>
                      <td>{row.HumanmiRNAID}</td>
                      <td>{row.Accession}</td>
                      <td>{row.Sequence}</td>
                      <td>{row.Seed}</td>
                      <td>{row.Position}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            {excelPath && (
              <a href={`/${excelPath}`} download className="btn btn-secondary">
                Download Complete Data
              </a>
            )}
          </>
        )}
      </div>
      <br /><br />
    </>
  )
}
