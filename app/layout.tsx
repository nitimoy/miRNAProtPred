import type { Metadata } from 'next'
import './globals.css'
import Navbar from './components/Navbar'
import Footer from './components/Footer'
import LoadingModal from './components/LoadingModal'

export const metadata: Metadata = {
  title: 'miRNAProtPred',
  description: 'A Web Tool for Prediction of Human miRNA Associated with Pathogens',
  icons: {
    icon: '/favicon.ico',
  },
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang="en">
      <head>
        <link
          rel="stylesheet"
          href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"
        />
      </head>
      <body>
        <Navbar />
        {children}
        <Footer />
        <LoadingModal />
      </body>
    </html>
  )
}
