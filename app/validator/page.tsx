'use client'

import { useState } from 'react'
import HeroSection from '../components/HeroSection'

interface ValidationResult {
  id: string
  found: boolean
  seed_match?: boolean
}

export default function Validator() {
  const [sequence, setSequence] = useState('')
  const [miRNAIds, setMiRNAIds] = useState('')
  const [results, setResults] = useState<ValidationResult[] | null>(null)
  const [error, setError] = useState('')
  const [isLoading, setIsLoading] = useState(false)

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    setIsLoading(true)
    setError('')
    setResults(null)

    window.dispatchEvent(new Event('loadingStart'))

    try {
      const response = await fetch('/api/validator', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequence, miRNA_ids: miRNAIds }),
      })

      const data = await response.json()

      if (!response.ok) {
        setError(data.error || 'An error occurred')
      } else {
        setResults(data.results)
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
      <HeroSection title="miRNA-Validator" />

      <div className="container mt-5">
        <h6 className="mt-4">
          • This tool validates specific miRNA IDs against a given sequence to check if they have
          potential binding sites.
        </h6>
        <h6 className="mt-4">
          • Enter a comma-separated list of miRNA IDs (e.g., miR-21, miR-155, miR-122) and your
          target sequence.
        </h6>

        <form onSubmit={handleSubmit} className="mt-4">
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
          <div className="form-group">
            <label htmlFor="miRNA_ids">Enter miRNA IDs (comma-separated):</label>
            <input
              type="text"
              className="form-control"
              id="miRNA_ids"
              required
              placeholder="e.g., miR-21, miR-155, miR-122"
              value={miRNAIds}
              onChange={(e) => setMiRNAIds(e.target.value)}
            />
          </div>
          <br />
          <button type="submit" className="btn btn-primary btn-lg" disabled={isLoading} style={{ width: '100%' }}>
            {isLoading ? 'Validating...' : 'Validate'}
          </button>
        </form>

        {error && (
          <div className="alert alert-danger mt-3" role="alert">
            {error}
          </div>
        )}

        {results && results.length > 0 && (
          <>
            <h2 className="mt-4">Validation Results</h2>
            <div style={{ overflowX: 'auto' }}>
              <table className="table table-striped">
                <thead>
                  <tr>
                    <th>miRNA ID</th>
                    <th>Found in Database</th>
                    <th>Seed Match in Sequence</th>
                  </tr>
                </thead>
                <tbody>
                  {results.map((result, index) => (
                    <tr key={index}>
                      <td>{result.id}</td>
                      <td>{result.found ? '✓ Yes' : '✗ No'}</td>
                      <td>
                        {result.found
                          ? result.seed_match
                            ? '✓ Yes'
                            : '✗ No'
                          : 'N/A'}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </>
        )}
      </div>
      <br /><br />
    </>
  )
}
