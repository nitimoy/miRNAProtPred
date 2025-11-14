interface HeroSectionProps {
  title: string
  subtitle?: string
}

export default function HeroSection({ title, subtitle }: HeroSectionProps) {
  return (
    <section className="hero-section text-center text-white bg-primary py-5">
      <div className="container">
        <h1 className="text-white font-weight-bold">{title}</h1>
        {subtitle && <h3 className="text-white font-weight-bold">{subtitle}</h3>}
      </div>
    </section>
  )
}
