'use client'

import { useEffect, useState } from 'react'

export default function LoadingModal() {
  const [isLoading, setIsLoading] = useState(false)

  useEffect(() => {
    // Listen for custom loading events
    const handleLoadingStart = () => setIsLoading(true)
    const handleLoadingEnd = () => setIsLoading(false)

    window.addEventListener('loadingStart', handleLoadingStart)
    window.addEventListener('loadingEnd', handleLoadingEnd)

    return () => {
      window.removeEventListener('loadingStart', handleLoadingStart)
      window.removeEventListener('loadingEnd', handleLoadingEnd)
    }
  }, [])

  if (!isLoading) return null

  return (
    <div className="modal show" style={{ display: 'flex' }}>
      <div className="modal-content">
        <h5>Processing Request</h5>
        <div className="text-center" style={{ marginTop: '20px' }}>
          <p>Please wait while we process your request. This may take a few moments.</p>
          <div className="spinner-border" role="status">
            <span style={{ position: 'absolute', width: '1px', height: '1px', overflow: 'hidden' }}>
              Loading...
            </span>
          </div>
        </div>
      </div>
    </div>
  )
}
