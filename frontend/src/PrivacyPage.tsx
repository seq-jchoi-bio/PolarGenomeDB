import { SectionIndex } from './SectionIndex';
import './privacy.css';

const sections = [
  ['overview', 'Overview'], ['records', 'Operational records'],
  ['analysis', 'Analysis data'], ['browser', 'Browser storage'],
  ['retention', 'Retention & deletion'], ['providers', 'Service providers'],
  ['contact', 'Your rights & contact'],
] as const;

export function PrivacyPage() {
  return <section className="page-section about-page privacy-page" aria-label="Privacy Policy">
    <div className="about-layout">
      <SectionIndex sections={sections} prefix="privacy-" className="about-index" label="Privacy policy sections" heading="On this page" />
      <div className="about-document">
        <section className="about-item" id="privacy-overview">
          <div className="privacy-version">Version 1.0 <span>Updated 5 October 2026</span></div>
          <h1>Privacy Policy</h1>
          <p>Polar Genome DB is a research platform for gene and genome analysis. We process analysis inputs and limited technical information to provide the service, resolve errors, and maintain security.</p>
        </section>
        <section className="about-item" id="privacy-records">
          <h2>Operational records</h2>
          <p>We record IP addresses, access times, service routes, request methods, response status, response size and processing time. Monitoring also uses browser and operating-system information, referring website domains, and approximate IP-based location information.</p>
          <p>These records help us investigate errors, protect the service, and understand usage. Our web access logs exclude search query strings, submitted sequences, cookies, and authentication credentials. Routes are recorded without individual gene or job identifiers.</p>
        </section>
        <section className="about-item" id="privacy-analysis">
          <h2>Analysis data</h2>
          <p>Gene identifiers, sequences, uploaded files, and analysis settings are processed to carry out your requests. Inputs, temporary files, job status, and results may be stored to run analyses and provide results. These are separate from web access logs.</p>
          <p>The tools are intended for genomic research data, not personal or patient information. Please avoid including personal information in submitted sequences, file names, or sequence headers.</p>
        </section>
        <section className="about-item" id="privacy-browser">
          <h2>Browser storage & authentication</h2>
          <p>The gene basket saves selected species and gene identifiers in your browser. You can remove entries from the basket or clear them through your browser's site-data settings. Entries are submitted to the server when you use them in an analysis.</p>
          <p>Security and administrative sign-in services may use authentication cookies. You can manage cookies in your browser settings; blocking them may affect protected features.</p>
        </section>
        <section className="about-item" id="privacy-retention">
          <h2>Retention & deletion</h2>
          {/* Publication requires verification of source-log and monitoring-database deletion,
              provider disclosures, and completion of the controller/contact details. */}
          <p>Our retention policy limits detailed visitor access records to a maximum of 90 days, after which they are deleted. Records may be removed sooner when no longer needed. These detailed records are excluded from routine backups.</p>
          <p>Aggregate usage statistics without IP addresses, submitted inputs, or individual identifiers are retained for long-term service reporting. Analysis files and results are kept separately for job processing and result retrieval, and are removed through routine cleanup rather than included in long-term backups.</p>
          <p>Expired records are removed from the relevant log files and monitoring databases. Where retention is required by law, only the necessary records are kept separately for the required period.</p>
        </section>
        <section className="about-item" id="privacy-providers">
          <h2>Service providers & safeguards</h2>
          <p>We use Oracle Cloud Infrastructure for hosting and Cloudflare for web delivery, security, and administrative authentication. Hosting is in Japan, and provider services may process information in other countries. Their processing and retention practices are separate from our own access-log policy.</p>
          <p>Connections use HTTPS, and access to operational records is restricted to authorized administrators. For information about provider practices, see the <a href="https://www.oracle.com/legal/privacy/" target="_blank" rel="noopener noreferrer">Oracle privacy information</a> and <a href="https://www.cloudflare.com/privacypolicy/" target="_blank" rel="noopener noreferrer">Cloudflare Privacy Policy</a>.</p>
        </section>
        <section className="about-item" id="privacy-contact">
          <h2>Your rights & contact</h2>
          <p>You may request access to, correction or deletion of, or restrictions on the processing of your personal information, subject to applicable law. The responsible institution and privacy contact will be listed below.</p>
          <dl className="privacy-contact-fields">
            <div><dt>Responsible institution</dt><dd>To be confirmed</dd></div>
            <div><dt>Privacy contact</dt><dd>Department and contact details to be confirmed</dd></div>
          </dl>
          <p>Updates to this policy will be posted on this page with the revision date.</p>
        </section>
      </div>
    </div>
  </section>;
}
