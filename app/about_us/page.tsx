import HeroSection from '../components/HeroSection'
import Image from 'next/image'

export default function AboutUs() {
  return (
    <>
      <HeroSection
        title="Meet Our Team"
        subtitle="Get to know the people behind miRNAProtPred"
      />

      <div className="container mt-5">
        <table className="table">
          <tbody>
            <tr>
              <td>
                <p>
                  <strong>Dr. Manisha Pritam (Principal Investigator)</strong>
                  &nbsp;
                  <a
                    href="https://www.linkedin.com/in/manisha-pritam-041a64148/"
                    target="_blank"
                    rel="noopener noreferrer"
                  >
                    <i className="fa fa-linkedin-square" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                  &nbsp;
                  <a href="https://www.manishapritam.me/" target="_blank" rel="noopener noreferrer">
                    <i className="fa fa-globe" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                  &nbsp;
                  <a href="mailto:manishapritam12@gmail.com">
                    <i className="fa fa-envelope" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                </p>
                <p>Post Doctoral Fellow</p>
                <p>National Institute of Allergy and Infectious Diseases</p>
                <p>National Institutes of Health, United States.</p>
              </td>
              <td className="text-center">
                <Image src="/img/manisha.jpg" width={180} height={200} alt="Dr. Manisha Pritam" />
              </td>
            </tr>
            <tr>
              <td>
                <p>
                  <strong>Somenath Dutta (Lead Investigator)</strong>
                  &nbsp;
                  <a
                    href="https://www.linkedin.com/in/somenath-dutta-a59181199/"
                    target="_blank"
                    rel="noopener noreferrer"
                  >
                    <i className="fa fa-linkedin-square" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                  &nbsp;
                  <a href="https://www.sduttamedicine.works/" target="_blank" rel="noopener noreferrer">
                    <i className="fa fa-globe" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                  &nbsp;
                  <a href="mailto:sduttabiotech@gmail.com">
                    <i className="fa fa-envelope" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                </p>
                <p>Doctoral Fellow</p>
                <p>Department of Biomolecular Engineering,</p>
                <p>Pusan National University, Busan, South Korea.</p>
              </td>
              <td className="text-center">
                <Image src="/img/somenath.jpg" width={180} height={200} alt="Somenath Dutta" />
              </td>
            </tr>
            <tr>
              <td>
                <p>
                  <strong>Nitimoy Mondal (Research Software Developer)</strong>
                  &nbsp;
                  <a href="https://www.linkedin.com/in/nitimoy/" target="_blank" rel="noopener noreferrer">
                    <i className="fa fa-linkedin-square" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                  &nbsp;
                  <a href="https://nitimoy.netlify.app/" target="_blank" rel="noopener noreferrer">
                    <i className="fa fa-globe" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                  &nbsp;
                  <a href="mailto:nitimoy@gmail.com">
                    <i className="fa fa-envelope" style={{ fontSize: '24px', color: 'blue' }}></i>
                  </a>
                </p>
                <p>Masters of Technology (MTech)</p>
                <p>Department of Computer Science,</p>
                <p>Rajasthan Technical University, Rajasthan, India</p>
              </td>
              <td className="text-center">
                <Image src="/img/nitimoy.jpg" width={180} height={200} alt="Nitimoy Mondal" />
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </>
  )
}
