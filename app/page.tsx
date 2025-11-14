import HeroSection from './components/HeroSection'
import Image from 'next/image'

export default function Home() {
  return (
    <>
      <HeroSection
        title="miRNAProtPred"
        subtitle="A Web Tool for Prediction of Human miRNA Associated with Pathogens."
      />

      <div className="container mt-5">
        <h2 className="text-center mb-4">Why miRNA-ProtPred?</h2>
        <p>
          miRNA-ProtPred predicts the human miRNAs that have potential to inhibit the expression of
          gene/protein of different pathogens.
        </p>
        <p>
          Most techniques predict miRNA against human proteins; however, miRNA-ProtPred is the first
          online server to predict human miRNA against 3&apos;UTR, protein, RNA, or DNA of human pathogens.
        </p>
        <div className="text-center my-4">
          <Image
            src="https://s11.gifyu.com/images/S43TY.gif"
            alt="miRNA-ProtPred Tool Diagram"
            width={800}
            height={400}
            className="img-fluid"
            unoptimized
          />
        </div>
        <h3 className="mt-4"><u>Why miRNA?</u></h3>
        <ul className="list-unstyled">
          <li className="mb-2">
            <strong>&bull;</strong> It plays an important role in gene regulation and protein expression.
          </li>
          <li className="mb-2">
            <strong>&bull;</strong> During pathogenic infection, the miRNA of the human host can bind with
            3&apos; UTR of mRNA of pathogen and can inhibit the gene/protein expression and ultimately helps
            in disease control.
          </li>
          <li className="mb-2">
            <strong>&bull;</strong> The FDA has not yet approved any miRNA-based therapeutics (as of 2023).
            However, some miRNA mimics are under different phases of clinical trials.
          </li>
        </ul>
      </div>
    </>
  )
}
