'use client'

import { useEffect, useState } from 'react'
import Link from 'next/link'

export default function Footer() {
  const [visitStats, setVisitStats] = useState({ visitCount: 0, uniqueVisitors: 0 })

  useEffect(() => {
    // Fetch visit statistics
    fetch('/api/visit-stats')
      .then(res => res.json())
      .then(data => setVisitStats(data))
      .catch(err => console.error('Error fetching visit stats:', err))
  }, [])

  const navItems = [
    { href: '/', label: 'Home' },
    { href: '/utr_prime', label: "3'UTRmiRNA-Pred" },
    { href: '/sequence_finder', label: 'miRNA-SeqFinder' },
    { href: '/validator', label: 'miRNA-Validator' },
    { href: 'https://mirnaprotpred.com/downloads/miRNAProtPred_Help.pdf', label: 'Help', external: true },
    { href: '/about_us', label: 'Team' },
    { href: '/contact_us', label: 'Contact Us' },
  ]

  return (
    <footer className="bg-cyan text-white text-center py-3">
      <div className="container">
        <ul className="footer-nav">
          {navItems.map((item) => (
            <li key={item.href} style={{ margin: '0 10px' }}>
              {item.external ? (
                <a href={item.href} target="_blank" rel="noopener noreferrer">
                  {item.label}
                </a>
              ) : (
                <Link href={item.href}>{item.label}</Link>
              )}
            </li>
          ))}
        </ul>
        <p className="text-center">© 2023 miRNA-ProtPred. All rights reserved.</p>
        <p>Unique Visitors: {visitStats.uniqueVisitors}</p>
        <p>Total Visits: {visitStats.visitCount}</p>
      </div>
    </footer>
  )
}
