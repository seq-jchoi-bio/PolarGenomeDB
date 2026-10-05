import { ArrowUpRight, Check, LockKeyhole, Files, Globe2, Dna, Network, Database, PanelsTopLeft } from 'lucide-react';
import './resources.css';
import { SectionIndex } from './SectionIndex';

const resourceSections = [
  ['data-availability', 'Data Availability'],
  ['downloadable-data', 'Downloadable Data'],
  ['platform-tools', 'Platform Tools'],
] as const;

const privateDatasets = [
  { name: 'Genome & protein sequences', detail: 'In-house genome assemblies and annotated protein sequences for Cq and Da.', format: 'FASTA' },
  { name: 'Gene annotations', detail: 'Gene, transcript, exon and CDS coordinates for the in-house assemblies.', format: 'GTF / GFF3' },
  { name: 'RNA-seq alignments & coverage', detail: 'Aligned reads, alignment indexes and derived coverage tracks.', format: 'BAM / BAI / CSI / bigWig' },
  { name: 'Expression datasets', detail: 'Raw count tables and sample metadata describing the RNA-seq experiments.', format: 'Count tables / Sample metadata' },
];

const referenceDatasets = [
  { name: 'Arabidopsis thaliana · TAIR10', detail: 'Genome FASTA, protein FASTA and gene annotation. Local annotation source: Ensembl Plants release 62.', format: 'FASTA / GTF', href: 'https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index', source: 'Ensembl Plants' },
  { name: 'Oryza sativa · IRGSP-1.0', detail: 'Genome, representative transcripts, proteins and functional annotation. Local annotation source: RAP-DB, 2026-02-05.', format: 'FASTA / GFF / Annotation tables', href: 'https://rapdb.dna.affrc.go.jp/download/irgsp1.html', source: 'RAP-DB' },
  { name: 'Pfam', detail: 'Protein family profiles used for domain annotation and comparison.', format: 'HMM profiles', href: 'https://www.ebi.ac.uk/interpro/download/', source: 'InterPro / Pfam' },
  { name: 'Gene Ontology', detail: 'Ontology terms and relationships used for GO enrichment and semantic analysis.', format: 'OBO', href: 'https://geneontology.org/docs/download-ontology/', source: 'Gene Ontology' },
];

const software = [
  { group: 'Sequence analysis', icon: Dna, tools: [
    ['BLAST+', 'Nucleotide and protein sequence searches.', 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'],
    ['DIAMOND', 'Fast protein searches and reference homolog discovery.', 'https://github.com/bbuchfink/diamond'],
    ['HMMER', 'Protein domain detection against Pfam profiles.', 'https://github.com/EddyRivasLab/hmmer'],
    ['MAFFT', 'Multiple protein sequence alignment.', 'https://mafft.cbrc.jp/alignment/software/'],
    ['IQ-TREE', 'Maximum-likelihood protein phylogenies.', 'https://iqtree.github.io/'],
    ['Primer3', 'PCR primer candidate design.', 'https://primer3.org/'],
    ['primer3-py', 'Primer thermodynamics and secondary-structure screening.', 'https://github.com/libnano/primer3-py'],
    ['Bowtie2', 'Genome-based exact-match primer specificity screening.', 'https://bowtie-bio.sourceforge.net/bowtie2/index.shtml'],
  ] },
  { group: 'GO analysis & visualization', icon: Network, tools: [
    ['GoFigure!', 'Semantic GO term grouping and summary visualization.', 'https://gitlab.com/evogenlab/GO-Figure'],
    ['scikit-learn', 'Multidimensional scaling within GoFigure!.', 'https://scikit-learn.org/'],
    ['Matplotlib', 'Semantic GO figure rendering through GoFigure!.', 'https://matplotlib.org/'],
    ['seaborn', 'Plot styling and palettes within GoFigure!.', 'https://seaborn.pydata.org/'],
  ] },
  { group: 'Data preparation', icon: Database, tools: [
    ['SAMtools', 'Genome FASTA and alignment file indexing.', 'https://www.htslib.org/'],
    ['HTSlib', 'bgzip compression and tabix annotation indexing.', 'https://www.htslib.org/'],
    ['deepTools', 'CPM-normalized RNA-seq bigWig generation.', 'https://github.com/deeptools/deepTools'],
    ['Biopython', 'Sequence records and annotated GenBank exports.', 'https://biopython.org/'],
    ['SQLite', 'Species metadata, annotations and analysis records.', 'https://www.sqlite.org/'],
    ['PyYAML', 'Manifest-driven dataset configuration.', 'https://pyyaml.org/'],
  ] },
  { group: 'Web platform', icon: PanelsTopLeft, tools: [
    ['JBrowse 2', 'Interactive genome and RNA-seq track visualization.', 'https://jbrowse.org/jb2/'],
    ['FastAPI', 'Backend APIs and analysis request handling.', 'https://fastapi.tiangolo.com/'],
    ['Uvicorn', 'ASGI application server for the backend.', 'https://github.com/Kludex/uvicorn'],
    ['React', 'Interactive multi-species analysis interface.', 'https://react.dev/'],
    ['TypeScript', 'Typed frontend application code.', 'https://www.typescriptlang.org/'],
    ['Vite', 'Frontend development and production builds.', 'https://vite.dev/'],
    ['Lucide', 'Interface icons.', 'https://lucide.dev/'],
  ] },
];

export function ResourcesPage() {
  return <section className="page-section resources-page" aria-label="Resources">
    <div className="resources-inner">
      <header className="resources-heading"><h1>Resources</h1><p>Genome datasets, reference resources and the tools behind Polar Genome DB.</p></header>
      <div className="resources-layout">
        <SectionIndex sections={resourceSections} className="resources-index" label="Resources sections" />
        <div className="resources-document">
          <section id="data-availability" className="resources-section">
            <h2>Data Availability</h2>
            <p>The in-house genome and RNA-seq datasets for <em>Colobanthus quitensis</em> and <em>Deschampsia antarctica</em> are not yet publicly released. They support browsing and analysis in this platform, but the underlying datasets are not offered for download here.</p>
            <dl className="resources-access-list">
              <div><dt>In-house genome datasets</dt><dd><LockKeyhole size={16} />Bulk download not available</dd></div>
              <div><dt>RNA-seq & expression datasets</dt><dd><LockKeyhole size={16} />Not publicly released</dd></div>
              <div><dt>Public reference datasets</dt><dd className="resources-available"><Check size={16} />Available from original providers</dd></div>
            </dl>
          </section>
          <section id="downloadable-data" className="resources-section">
            <h2>Downloadable Data</h2>
            <p>Source datasets used to build the platform, grouped by release status. Locked entries are not available for download.</p>
            <h3 className="resources-dataset-heading"><Files size={20} aria-hidden="true" />In-house datasets <span>Cq / Da</span></h3>
            <div className="resources-catalog-labels" aria-hidden="true"><span>Dataset</span><span>File format</span><span>Access</span></div>
            <div className="resources-download-list resources-private-datasets">
              {privateDatasets.map(item => <article className="resources-download-row" key={item.name}>
                <div><h4>{item.name}</h4><p>{item.detail}</p></div>
                <span className="resources-formats">{item.format}</span>
                <span className="resources-locked"><LockKeyhole size={17} aria-hidden="true" />Not yet released</span>
              </article>)}
            </div>
            <h3 className="resources-dataset-heading"><Globe2 size={20} aria-hidden="true" />Public reference datasets <span>External sources</span></h3>
            <div className="resources-catalog-labels" aria-hidden="true"><span>Dataset</span><span>File format</span><span>Provider</span></div>
            <div className="resources-download-list resources-public-datasets">
              {referenceDatasets.map(item => <article className="resources-download-row" key={item.name}>
                <div><h4>{item.name}</h4><p>{item.detail}</p></div>
                <span className="resources-formats">{item.format}</span>
                <a className="resources-provider-link" href={item.href} target="_blank" rel="noreferrer">{item.source}<ArrowUpRight size={17} aria-hidden="true" /></a>
              </article>)}
            </div>
          </section>
          <section id="platform-tools" className="resources-section">
            <h2>Platform Tools</h2>
            <p>Established bioinformatics software supports the analysis pipelines, with a shared API and multi-species web interface.</p>
            {software.map(group => <div className="resources-software-group" key={group.group}>
              <h3><group.icon size={23} strokeWidth={1.6} aria-hidden="true" />{group.group}</h3>
              <dl>{group.tools.map(([name, role, href]) => <div key={name}>
                <dt><a href={href} target="_blank" rel="noreferrer">{name}<ArrowUpRight size={15} aria-hidden="true" /></a></dt>
                <dd>{role}</dd>
              </div>)}</dl>
            </div>)}
            <div className="resources-reference-data"><h3>Reference data</h3><p>Reference-based analyses use Arabidopsis TAIR10 for dicots and rice RAP-DB IRGSP-1.0 for monocots, together with Pfam protein families and the Gene Ontology. These are reference datasets, distinct from the in-house polar genome assemblies.</p></div>
          </section>
        </div>
      </div>
    </div>
  </section>;
}
