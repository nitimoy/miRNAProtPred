'use client'

import Link from 'next/link'
import { usePathname } from 'next/navigation'
import Image from 'next/image'
import { useState } from 'react'

export default function Navbar() {
  const pathname = usePathname()
  const [isOpen, setIsOpen] = useState(false)

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
    <nav className="navbar">
      <div className="container" style={{ marginBottom: 0 }}>
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', width: '100%' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '1rem' }}>
            <Link href="/">
              <Image src="/logo.png" alt="miRNA-ProtPred Logo" width={90} height={90} />
            </Link>
            <Link href="/" className="navbar-brand">
              miRNAProtPred
            </Link>
          </div>

          <button
            className="navbar-toggler"
            onClick={() => setIsOpen(!isOpen)}
            style={{
              display: 'none',
              background: 'none',
              border: 'none',
              color: 'white',
              fontSize: '1.5rem',
              cursor: 'pointer'
            }}
          >
            ☰
          </button>

          <ul className="navbar-nav" style={{ display: isOpen ? 'flex' : 'flex' }}>
            {navItems.map((item) => (
              <li key={item.href} className="nav-item">
                {item.external ? (
                  <a
                    href={item.href}
                    className="nav-link"
                    target="_blank"
                    rel="noopener noreferrer"
                  >
                    {item.label}
                  </a>
                ) : (
                  <Link
                    href={item.href}
                    className={`nav-link ${pathname === item.href ? 'active' : ''}`}
                  >
                    {item.label}
                  </Link>
                )}
              </li>
            ))}
          </ul>
        </div>
      </div>

      <style jsx>{`
        .navbar-toggler {
          display: none;
        }

        @media (max-width: 768px) {
          .navbar-toggler {
            display: block !important;
          }

          .navbar-nav {
            display: ${isOpen ? 'flex' : 'none'} !important;
            flex-direction: column;
            width: 100%;
            margin-top: 1rem;
          }
        }
      `}</style>
    </nav>
  )
}
