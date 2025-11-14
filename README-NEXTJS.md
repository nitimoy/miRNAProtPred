# miRNAProtPred - Next.js Version

A Web Tool for Prediction of Human miRNA Associated with Pathogens.

This is the Next.js conversion of the original Flask application, providing a modern, performant web interface for miRNA prediction and analysis.

## Features

- **3'UTR miRNA-Pred**: Predict human miRNAs that can inhibit pathogen 3'UTR sequences
- **miRNA-SeqFinder**: Find potential human miRNA binding sites in any DNA, RNA, or protein sequence
- **miRNA-Validator**: Validate specific miRNA IDs against target sequences
- **Team Page**: Meet the researchers behind miRNAProtPred
- **Contact Form**: Get in touch with the team

## Technology Stack

- **Frontend**: Next.js 14 (App Router), React 18, TypeScript
- **Backend**: Next.js API Routes
- **Database**: PostgreSQL
- **Email**: Nodemailer (Gmail SMTP)
- **Bioinformatics**: Custom Boyer-Moore algorithm, NCBI BLAST/Entrez integration
- **Excel Export**: ExcelJS

## Prerequisites

- Node.js 18+ and npm
- PostgreSQL database
- Gmail account for SMTP (or other email service)
- NCBI Entrez email registration
- Geolocation API token (ipinfo.io)

## Installation

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd miRNAProtPred
   ```

2. **Install dependencies**
   ```bash
   npm install
   ```

3. **Set up environment variables**

   Copy `.env.example` to `.env.local`:
   ```bash
   cp .env.example .env.local
   ```

   Edit `.env.local` and fill in your configuration:
   ```env
   # Database Configuration
   DB_USERNAME=your_db_username
   DB_PASSWORD=your_db_password
   DB_HOST=localhost
   DB_PORT=5432
   DB_NAME=mirnaprotpred

   # Email Configuration
   EMAIL_USER=your_email@gmail.com
   EMAIL_PASS=your_app_password

   # Geolocation API
   GEOLOCATION_API_TOKEN=your_ipinfo_token

   # NCBI Entrez API
   NCBI_EMAIL=your_email@example.com
   ```

4. **Set up the database**

   Create the PostgreSQL database and tables:
   ```sql
   CREATE DATABASE mirnaprotpred;

   \c mirnaprotpred

   -- Visits table
   CREATE TABLE visits (
       id SERIAL PRIMARY KEY,
       ip_address VARCHAR(50) UNIQUE,
       city VARCHAR(100),
       country VARCHAR(100),
       visit_count INTEGER DEFAULT 1,
       visit_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
   );

   -- miRNA database table (import your existing mirdb data)
   CREATE TABLE mirdb (
       r_id SERIAL PRIMARY KEY,
       Description TEXT,
       HumanmiRNAID VARCHAR(50),
       mir_id VARCHAR(50),
       Accession VARCHAR(50),
       Sequence TEXT,
       seed_1 VARCHAR(50),
       seed_2 VARCHAR(50),
       seed_3 VARCHAR(50)
   );
   ```

   Import your miRNA database data into the `mirdb` table.

5. **Run the development server**
   ```bash
   npm run dev
   ```

   Open [http://localhost:3000](http://localhost:3000) in your browser.

## Building for Production

1. **Build the application**
   ```bash
   npm run build
   ```

2. **Start the production server**
   ```bash
   npm start
   ```

## Project Structure

```
miRNAProtPred/
├── app/                        # Next.js app directory
│   ├── api/                    # API routes
│   │   ├── contact/           # Contact form endpoint
│   │   ├── sequence-finder/   # Sequence finder endpoint
│   │   ├── utr-prime/         # UTR prime endpoint
│   │   ├── validator/         # Validator endpoint
│   │   └── visit-stats/       # Visit statistics endpoint
│   ├── about_us/              # About us page
│   ├── components/            # Shared React components
│   │   ├── Footer.tsx
│   │   ├── HeroSection.tsx
│   │   ├── LoadingModal.tsx
│   │   └── Navbar.tsx
│   ├── contact_us/            # Contact us page
│   ├── sequence_finder/       # Sequence finder page
│   ├── utr_prime/             # UTR prime page
│   ├── validator/             # Validator page
│   ├── globals.css            # Global styles
│   ├── layout.tsx             # Root layout
│   └── page.tsx               # Home page
├── lib/                       # Utility libraries
│   ├── algorithms.ts          # Boyer-Moore algorithm
│   ├── db.ts                  # Database connection
│   └── ncbi.ts                # NCBI API integration
├── public/                    # Static assets
│   ├── downloads/             # Generated Excel files
│   ├── img/                   # Images
│   ├── favicon.ico
│   └── logo.png
├── .env.example               # Environment variables template
├── .env.local                 # Your local environment variables (not in git)
├── next.config.js             # Next.js configuration
├── package.json               # Dependencies
├── tsconfig.json              # TypeScript configuration
└── README-NEXTJS.md           # This file
```

## Key Differences from Flask Version

1. **Server-Side Rendering**: Next.js provides automatic server-side rendering and static site generation for better performance

2. **API Routes**: All backend logic is now in Next.js API routes instead of Flask routes

3. **TypeScript**: The entire application is now type-safe with TypeScript

4. **React Components**: The UI is built with reusable React components instead of Jinja2 templates

5. **Modern Build System**: Uses Next.js's optimized build system with automatic code splitting

6. **Image Optimization**: Next.js automatically optimizes images for better performance

## Notes on BLAST Integration

The NCBI BLAST integration is simplified in this version. For protein sequences that require BLAST search:

- The current implementation provides a basic BLAST workflow
- For production use, consider implementing proper polling with delays
- Alternatively, you may want to run BLAST searches in a background job queue
- DNA and RNA sequences work without BLAST and are fully functional

## Email Configuration

For Gmail SMTP:
1. Enable 2-factor authentication on your Google account
2. Generate an App Password: https://myaccount.google.com/apppasswords
3. Use the App Password in your `.env.local` file

## Deployment

### Vercel (Recommended)

1. Push your code to GitHub
2. Import your repository in Vercel
3. Add environment variables in Vercel dashboard
4. Deploy

### Other Platforms

The application can be deployed to any platform that supports Node.js:
- AWS (EC2, Elastic Beanstalk, or Lambda with API Gateway)
- Google Cloud Platform
- Azure
- DigitalOcean
- Heroku

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

© 2023 miRNA-ProtPred. All rights reserved.

## Team

- **Dr. Manisha Pritam** - Principal Investigator
- **Somenath Dutta** - Lead Investigator
- **Nitimoy Mondal** - Research Software Developer

## Citation

If you use miRNAProtPred in your research, please cite:

[Add citation information here]

## Support

For questions or issues, please contact us through the Contact Us page or open an issue on GitHub.
