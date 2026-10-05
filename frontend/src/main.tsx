import React, { useEffect, useMemo, useRef, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { createViewState, JBrowseLinearGenomeView } from '@jbrowse/react-linear-genome-view2';
import '@fontsource/geist/400.css';
import '@fontsource/geist/500.css';
import '@fontsource/geist/600.css';
import '@fontsource/geist/700.css';
import '@fontsource/geist/latin-400-italic.css';
import '@fontsource/ibm-plex-mono/400.css';
import '@fontsource/ibm-plex-mono/500.css';
import '@fontsource/ibm-plex-mono/700.css';
import { BarChart3, BookOpen, Check, CheckCheck, ChevronDown, ChevronLeft, ChevronRight, CircleHelp, ClipboardPaste, Copy, Database, Dna, Download, FileDown, FileText, FlaskConical, GitCompareArrows, Github, GitBranch, Layers3, ListX, LoaderCircle, Minus, Network, Play, Plus, RotateCcw, Search, ShoppingBasket, Upload } from 'lucide-react';
import { API_BASE, apiGet, apiPost, SpeciesSummary } from './api/client';
import './styles.css';
import './fonts.css';
import { exportTree, exportSvgPlot } from './exportTree';
import { goTermsTsv, wrapGoLabel } from './goExport';
import { spliceExons, projectSplicedFeature } from './transcriptGeometry';
import { classifyExpressionPattern, EXPRESSION_PATTERN_POLICY } from './expressionPattern';
import { normalizePrimerRequest, annotatedTranscriptLength, primerHitLabel, primerExport, primerEmptyReason } from './primerPresentation';
import { SpeciesPage } from './SpeciesPage';
import { ResourcesPage } from './ResourcesPage';
import { PrivacyPage } from './PrivacyPage';
import { SectionIndex } from './SectionIndex';
import { BasketProvider, BasketButton, BasketMenu } from './Basket';

type SpeciesIndex = {
  species: SpeciesSummary[];
  species_count: number;
};

type Module = {
  id: string;
  title: string;
  scope: string;
  status: string;
  icon: React.ReactNode;
};

type Page = 'home' | 'about' | 'species' | 'tools' | 'resources' | 'privacy';
type NavMenu = 'species' | 'tools' | null;
type ToolId = 'genome' | 'gene' | 'blast' | 'primer' | 'similar' | 'go';
type SearchProgram = 'blastn' | 'blastp' | 'diamond';

const maxSequenceQueriesForProgram = (program: SearchProgram) => (program === 'diamond' ? 100 : 10);
const countSequenceQueries = (query: string) => {
  const normalized = query.trim();
  if (!normalized) return 0;
  if (!normalized.startsWith('>')) return 1;
  return normalized.split(/\r?\n/).filter((line) => line.startsWith('>')).length;
};
const sequenceQueryIds = (query: string) => {
  const normalized = query.trim();
  if (!normalized) return [];
  if (!normalized.startsWith('>')) return ['query'];
  return normalized.split(/\r?\n/)
    .filter((line) => line.startsWith('>'))
    .map((line) => line.slice(1).trim().split(/\s+/, 1)[0])
    .filter(Boolean);
};
const sequenceQueryRecords = (query: string) => {
  const normalized = query.trim();
  if (!normalized) return [] as Array<{ id: string; sequence: string }>;
  if (!normalized.startsWith('>')) return [{ id: 'query', sequence: normalized.replace(/\s+/g, '') }];
  const records: Array<{ id: string; sequence: string }> = [];
  let id = '';
  let sequence = '';
  for (const line of normalized.split(/\r?\n/)) {
    if (line.startsWith('>')) {
      if (id) records.push({ id, sequence: sequence.replace(/\s+/g, '') });
      id = line.slice(1).trim().split(/\s+/, 1)[0] ?? '';
      sequence = '';
    } else {
      sequence += line.trim();
    }
  }
  if (id) records.push({ id, sequence: sequence.replace(/\s+/g, '') });
  return records.filter((record) => record.id);
};
const alphabeticalLabel = (index: number) => {
  let value = index + 1;
  let label = '';
  while (value > 0) {
    value -= 1;
    label = String.fromCharCode(65 + (value % 26)) + label;
    value = Math.floor(value / 26);
  }
  return label;
};
const sequenceExampleCache = new Map<string, ApiRecord[]>();
const sequenceExampleRequests = new Map<string, Promise<ApiRecord[]>>();
const loadSequenceSearchExamples = (speciesId: string, program: SearchProgram) => {
  const key = `${speciesId}:${program}`;
  const cached = sequenceExampleCache.get(key);
  if (cached) return Promise.resolve(cached);
  const pending = sequenceExampleRequests.get(key);
  if (pending) return pending;
  const request = apiGet<ApiRecord>(`/api/species/${speciesId}/search-examples?program=${program}`)
    .then((response) => {
      const records = Array.isArray(response.records) ? response.records as ApiRecord[] : [];
      if (!records.length) throw new Error('No example sequences are available for this species.');
      sequenceExampleCache.set(key, records);
      return records;
    })
    .finally(() => sequenceExampleRequests.delete(key));
  sequenceExampleRequests.set(key, request);
  return request;
};

const modules: Module[] = [
  { id: 'genome', title: 'Genome Browser', scope: 'FASTA / GFF3 / bigWig', status: 'Runtime ready', icon: <Layers3 size={34} /> },
  { id: 'gene', title: 'Gene Search', scope: 'integrated gene report', status: 'Integrated report', icon: <Search size={34} /> },
  { id: 'expression', title: 'Expression Matrix', scope: 'raw count / CPM', status: 'Optional data', icon: <BarChart3 size={34} /> },
  { id: 'similar', title: 'Similar Gene', scope: 'dicot / monocot reference', status: 'Precomputed', icon: <GitBranch size={34} /> },
  { id: 'go', title: 'GO Annotation', scope: 'GO enrichment through reference homologs', status: 'GO inference', icon: <Network size={34} /> },
  { id: 'primer', title: 'PCR Primer', scope: 'primer3 / BLAST specificity', status: 'Job API', icon: <FlaskConical size={34} /> },
  { id: 'blast', title: 'Sequence Search', scope: 'BLAST / DIAMOND', status: 'Job API', icon: <Dna size={34} /> },
  { id: 'assets', title: 'Build Assets', scope: 'SQLite / indexes / JBrowse', status: 'Manifest driven', icon: <Database size={34} /> },
];

const entryModules = modules.filter((module) => ['gene', 'genome', 'blast', 'primer', 'similar', 'go'].includes(module.id));
const moduleDescriptions: Record<string, string> = {
  genome: 'Explore genome assemblies with gene annotation and RNA-seq coverage tracks through an embedded JBrowse viewer.',
  gene: 'Search one gene for an integrated report, or paste multiple gene IDs to compare annotation, expression, and reference evidence.',
  similar: 'Resolve reference protein isoforms and identify homologous genes in the selected polar genome.',
  go: 'Identify enriched biological processes, molecular functions, and cellular components from reference-inferred GO annotations.',
  primer: 'Design genomic PCR primer pairs for a selected gene, with configurable product size and genome specificity screening.',
  blast: 'Search nucleotide or protein sequences against genome and protein resources prepared for the selected species.',
};
const speciesImages: Record<string, string> = {
  Cq: '/images/withheld.svg',
  Da: '/images/withheld.svg',
};

const formatSpeciesShortName = (scientificName: string) => {
  const [genus, speciesName] = scientificName.split(' ');
  if (!genus || !speciesName) return scientificName;
  return `${genus.charAt(0)}. ${speciesName}`;
};
const homeSlides = [
  '/images/withheld.svg',
  '/images/withheld.svg',
  '/images/withheld.svg',
  '/images/withheld.svg',
];

const pathToPage = (pathname: string): Page => {
  const first = pathname.split('/').filter(Boolean)[0];
  if (first === 'about' || first === 'species' || first === 'tools' || first === 'resources' || first === 'privacy') return first;
  if (first === 'modules' || first === 'workspace') return 'tools';
  if (first === 'system') return 'resources';
  return 'home';
};

const pageToPath = (page: Page) => (page === 'home' ? '/' : `/${page}`);

const exampleGenesForSpecies = (speciesId: string) => (
  speciesId === 'Da'
    ? ['DaV2_01g000010', 'DaV2_01g000020', 'DaV2_01g000030']
    : ['Cq01g000010', 'Cq01g000020', 'Cq01g000030']
);
const primerExampleGeneForSpecies = (speciesId: string) => (
  speciesId === 'Da' ? 'DaV2_01g032610' : 'Cq06g010610'
);
const geneSearchExampleGenesForSpecies = (speciesId: string) => (
  speciesId === 'Da'
    ? ['DaV2_01g032610', 'DaV2_12g031450', 'DaV2_02g043450']
    : ['Cq01g010020', 'Cq06g010680', 'Cq06g010610']
);
const goEnrichmentExampleGenesForSpecies = (speciesId: string) => (
  speciesId === 'Da'
    ? ['DaV2_01g006350', 'DaV2_01g008980', 'DaV2_01g018370', 'DaV2_01g019210', 'DaV2_01g021630', 'DaV2_01g022880', 'DaV2_01g022890', 'DaV2_01g022910', 'DaV2_01g022920', 'DaV2_01g022940', 'DaV2_01g022950', 'DaV2_01g022970', 'DaV2_01g022990', 'DaV2_01g024940', 'DaV2_01g024980', 'DaV2_01g024990', 'DaV2_01g025000', 'DaV2_01g025010', 'DaV2_01g025030', 'DaV2_01g025040']
    : ['Cq10g005380', 'Cq19g002780', 'Cq21g005260', 'Cq17g006080', 'Cq26g007460', 'Cq07g004470', 'Cq12g002450', 'Cq12g000260', 'Cq34g007590', 'Cq32g004110', 'Cq12g000270', 'Cq23g007510', 'Cq29g006350', 'Cq07g003910', 'Cq21g003120', 'Cq26g006710', 'Cq30g001470', 'Cq20g007400', 'Cq21g004750', 'Cq02g010370']
);
const defaultGeneForSpecies = (speciesId: string) => (speciesId === 'Da' ? 'DaV2_01g000020' : 'Cq01g000010');
const defaultGeneForTool = (toolId: ToolId, speciesId: string) => ((toolId === 'gene' || toolId === 'similar' || toolId === 'primer') ? '' : defaultGeneForSpecies(speciesId));
const referenceContextForTaxonGroup = (taxonGroup: string) => (
  taxonGroup === 'monocot'
    ? { group: 'Monocot', scientificName: 'Oryza sativa Japonica Group', database: 'RAP-DB IRGSP-1.0', examples: ['Os04g0628100', 'RPS5', 'Os09t0251500-01'] }
    : { group: 'Dicot', scientificName: 'Arabidopsis thaliana', database: 'TAIR10', examples: ['AT5G66030', 'AT3G51390', 'RIN2'] }
);
const isToolId = (value: string): value is ToolId => ['genome', 'gene', 'blast', 'primer', 'similar', 'go'].includes(value);
const toolFromSearch = (search: string): ToolId => {
  const tool = new URLSearchParams(search).get('tool');
  return tool && isToolId(tool) ? tool : 'genome';
};

const demoSpecies: SpeciesSummary[] = [
  {
    species_id: 'Cq',
    scientific_name: 'Colobanthus quitensis',
    display_name: 'Cq',
    taxon_group: 'dicot',
    assembly_version: 'pilot assembly',
    counts: { genes: 46492, proteins: 46492, bigwig_tracks: 12 },
    expression_value_types: [],
  },
  {
    species_id: 'Da',
    scientific_name: 'Deschampsia antarctica',
    display_name: 'Da',
    taxon_group: 'monocot',
    assembly_version: 'large-genome test assembly',
    counts: { genes: 50848, proteins: 50848, bigwig_tracks: 8 },
    expression_value_types: [],
  },
];

function App() {
  const [species, setSpecies] = useState<SpeciesSummary[]>([]);
  const [selectedSpeciesId, setSelectedSpeciesId] = useState(() => new URLSearchParams(window.location.search).get('species') || 'Cq');
  const [error, setError] = useState('');
  const [page, setPage] = useState<Page>(() => pathToPage(window.location.pathname));
  const [showReleaseNotes, setShowReleaseNotes] = useState(false);
  const [activeSlide, setActiveSlide] = useState(0);
  const [openNavMenu, setOpenNavMenu] = useState<NavMenu>(null);
  const [activeTool, setActiveTool] = useState<ToolId>(() => toolFromSearch(window.location.search));
  const [geneTarget, setGeneTarget] = useState(() => new URLSearchParams(window.location.search).get('gene') ?? '');
  const [retainedTools, setRetainedTools] = useState<ToolId[]>([]);
  const [basketInput, setBasketInput] = useState<{ speciesId: string; geneIds: string[] } | null>(null);
  const resultScroll = useRef<Partial<Record<ToolId, number>>>({});
  const releaseNotesRef = useRef<HTMLDivElement | null>(null);
  const headerNavRef = useRef<HTMLElement | null>(null);

  useEffect(() => {
    const id = window.location.hash.slice(1);
    if (!id) return;
    const frame = requestAnimationFrame(() => document.getElementById(id)?.scrollIntoView());
    return () => cancelAnimationFrame(frame);
  }, [page]);

  useEffect(() => {
    const onPopState = () => {
      const nextPage = pathToPage(window.location.pathname);
      setPage(nextPage);
      if (nextPage === 'tools') {
        const tool = toolFromSearch(window.location.search);
        setActiveTool(tool);
        const params = new URLSearchParams(window.location.search);
        setGeneTarget(params.get('gene') ?? '');
        if (params.get('species')) setSelectedSpeciesId(params.get('species')!);
        if (resultScroll.current[tool] != null) requestAnimationFrame(() => requestAnimationFrame(() => window.scrollTo({ top: resultScroll.current[tool] ?? 0, behavior: 'instant' })));
      }
    };
    window.addEventListener('popstate', onPopState);
    return () => window.removeEventListener('popstate', onPopState);
  }, []);

  useEffect(() => {
    if (!showReleaseNotes) return undefined;

    const onPointerDown = (event: PointerEvent) => {
      if (!releaseNotesRef.current?.contains(event.target as Node)) {
        setShowReleaseNotes(false);
      }
    };

    document.addEventListener('pointerdown', onPointerDown);
    return () => document.removeEventListener('pointerdown', onPointerDown);
  }, [showReleaseNotes]);

  useEffect(() => {
    if (!openNavMenu) return undefined;

    const onPointerDown = (event: PointerEvent) => {
      if (!headerNavRef.current?.contains(event.target as Node)) {
        setOpenNavMenu(null);
      }
    };

    document.addEventListener('pointerdown', onPointerDown);
    return () => document.removeEventListener('pointerdown', onPointerDown);
  }, [openNavMenu]);

  useEffect(() => {
    const timer = window.setInterval(() => {
      setActiveSlide((current) => (current + 1) % homeSlides.length);
    }, 6000);
    return () => window.clearInterval(timer);
  }, []);

  function navigate(nextPage: Page, nextTool?: ToolId) {
    setPage(nextPage);
    setOpenNavMenu(null);
    const nextPath = nextPage === 'tools' && nextTool ? `/tools?tool=${nextTool}` : pageToPath(nextPage);
    window.history.pushState(null, '', nextPath);
    window.scrollTo({ top: 0, behavior: 'smooth' });
  }

  function selectSpeciesFromNav(speciesId: string) {
    setSelectedSpeciesId(speciesId);
    navigate('species');
  }

  function openTool(toolId: string) {
    if (isToolId(toolId)) {
      setGeneTarget('');
      setActiveTool(toolId);
      navigate('tools', toolId);
      return;
    }
  }

  function openCandidateGene(speciesId: string, geneId: string) {
    resultScroll.current[activeTool] = window.scrollY;
    setRetainedTools((current) => [...new Set([...current, activeTool])]);
    window.history.replaceState(window.history.state, '', `/tools?tool=${activeTool}&species=${encodeURIComponent(speciesId)}`);
    window.history.pushState(null, '', `/tools?tool=gene&species=${encodeURIComponent(speciesId)}&gene=${encodeURIComponent(geneId)}`);
    setSelectedSpeciesId(speciesId);
    setGeneTarget(geneId);
    setActiveTool('gene');
    window.scrollTo({ top: 0, behavior: 'instant' });
  }

  useEffect(() => {
    apiGet<SpeciesIndex>('/api/species')
      .then((speciesIndex) => {
        setSpecies(speciesIndex.species);
        if (speciesIndex.species.length && !speciesIndex.species.some((item) => item.species_id === selectedSpeciesId)) {
          setSelectedSpeciesId(speciesIndex.species[0].species_id);
        }
      })
      .catch(() => {
        setSpecies(demoSpecies);
        setError('Preview mode: backend API is not connected');
      });
  }, []);

  const selectedSpecies = useMemo(
    () => species.find((item) => item.species_id === selectedSpeciesId) ?? species[0],
    [species, selectedSpeciesId],
  );

  return (
    <main className="app-shell">
      <header className="global-header">
        <button type="button" className="brand-block" onClick={() => navigate('home')} aria-label="Go to home">
          <span className="brand-mark" aria-hidden="true">
            <img src="/images/withheld.svg" alt="" />
          </span>
          <div>
            <h1>Polar Genome DB</h1>
            <p>Genomics platform for polar ecosystems</p>
          </div>
        </button>
        <nav className="header-nav" aria-label="Primary navigation" ref={headerNavRef}>
          <button type="button" className={page === 'about' ? 'active' : ''} onClick={() => navigate('about')}>About</button>
          <div className="nav-menu">
            <button
              type="button"
              className={page === 'species' ? 'active' : ''}
              onClick={() => setOpenNavMenu((current) => (current === 'species' ? null : 'species'))}
              aria-expanded={openNavMenu === 'species'}
            >
              Species
              <ChevronDown size={17} />
            </button>
            {openNavMenu === 'species' && (
              <div className="nav-popover species-popover">
                {species.map((item) => (
                  <button type="button" key={item.species_id} onClick={() => selectSpeciesFromNav(item.species_id)}>
                    <img src={speciesImages[item.species_id] ?? '/images/withheld.svg'} alt="" />
                    <span>
                      <strong>{item.scientific_name}</strong>
                      <small><em>{formatSpeciesShortName(item.scientific_name)}</em> · {item.taxon_group}</small>
                    </span>
                  </button>
                ))}
              </div>
            )}
          </div>
          <div className="nav-menu">
            <button
              type="button"
              className={page === 'tools' ? 'active' : ''}
              onClick={() => setOpenNavMenu((current) => (current === 'tools' ? null : 'tools'))}
              aria-expanded={openNavMenu === 'tools'}
            >
              Tools
              <ChevronDown size={17} />
            </button>
            {openNavMenu === 'tools' && (
              <div className="nav-popover tools-popover">
                {entryModules.map((module) => (
                  <button type="button" key={module.id} onClick={() => openTool(module.id)}>
                    {module.icon}
                    <span>
                      <strong>{module.title}</strong>
                      <small>{module.scope}</small>
                    </span>
                  </button>
                ))}
              </div>
            )}
          </div>
          <button type="button" className={page === 'resources' ? 'active' : ''} onClick={() => navigate('resources')}>Resources</button>
          <button type="button" className="icon-nav-button" aria-label="Tutorial" data-tooltip="Tutorial">
            <BookOpen size={21} />
          </button>
          <BasketMenu species={species} onSearch={(speciesId, geneIds) => {
            setSelectedSpeciesId(speciesId);
            setBasketInput({ speciesId, geneIds });
            openTool('gene');
          }} />
          <a className="icon-nav-button" href="https://github.com/" target="_blank" rel="noreferrer" aria-label="Creator GitHub" data-tooltip="Creator GitHub">
            <Github size={21} />
          </a>
        </nav>
      </header>

      {page === 'home' && (
        <section className="landing-gate">
          <div className="image-slider" aria-label="Antarctic plant field image slideshow">
            {homeSlides.map((src, index) => (
              <img
                key={src}
                src={src}
                alt={`Antarctic plant landscape slide ${index + 1}`}
                className={index === activeSlide ? 'active' : ''}
              />
            ))}
            <button type="button" className="slide-control previous" onClick={() => setActiveSlide((current) => (current + homeSlides.length - 1) % homeSlides.length)} aria-label="Previous image">
              <ChevronLeft size={24} />
            </button>
            <button type="button" className="slide-control next" onClick={() => setActiveSlide((current) => (current + 1) % homeSlides.length)} aria-label="Next image">
              <ChevronRight size={24} />
            </button>
            <div className="slide-dots" aria-label="Select image">
              {homeSlides.map((src, index) => (
                <button
                  type="button"
                  key={src}
                  className={index === activeSlide ? 'active' : ''}
                  onClick={() => setActiveSlide(index)}
                  aria-label={`Show image ${index + 1}`}
                />
              ))}
            </div>
            <span className="image-credit">© Polar Genome DB</span>
          </div>

          <div className="landing-copy">
            <div className="version-badge-wrap" ref={releaseNotesRef}>
              <button type="button" className="version-badge" aria-expanded={showReleaseNotes} aria-controls="home-release-notes" onClick={() => setShowReleaseNotes((value) => !value)}>
                <span className="version-badge-label">Beta</span>
                <span className="version-badge-number">1.0.1</span>
                <ChevronDown size={13} strokeWidth={1.8} aria-hidden="true" />
              </button>
              {showReleaseNotes && (
                <div className="release-notes" id="home-release-notes">
                  <strong>Beta 1.0.1 · October 5, 2026</strong>
                  <span>Updated privacy policy.</span>
                  <span>Prepared source code for public release.</span>
                  <span>Added system architecture documentation.</span>
                </div>
              )}
            </div>
            <h2>Polar Genome DB</h2>
            <p>
              Polar Genome DB is an integrated genomics platform for organisms from polar and cold-environment ecosystems. It provides genome browsing, gene and sequence search, functional annotation, comparative genomics, expression visualization, and primer design tools to support the study of genetic diversity, adaptation, and evolution in cold-adapted species.
            </p>

            <div className="beta-disclaimer">
              <strong>Beta service notice</strong>
              <ul>
                <li>Polar Genome DB is currently provided as a beta research database prototype.</li>
                <li>Some functions may be limited, revised, or temporarily unavailable during active development.</li>
                <li>The current service covers two species only: Cq and Da.</li>
                <li>Data, annotations, inferred functions, and primer predictions should be independently verified before publication or experimental use.</li>
              </ul>
            </div>
          </div>

          <div className="entry-grid" aria-label="Main entry points">
            {entryModules.map((module) => (
              <button type="button" key={module.id} className="entry-link" onClick={() => openTool(module.id)}>
                <span>{module.icon}</span>
                <strong>{module.title}</strong>
                <small>{module.scope}</small>
                <ChevronRight size={22} />
              </button>
            ))}
          </div>

        </section>
      )}

      {page !== 'home' && page !== 'privacy' && error && <div className="status-line muted-line">{error}</div>}

      {page === 'about' && (
        <section className="page-section about-page" aria-label="About Polar Genome DB">
          <div className="about-layout">
            <AboutIndex />
            <div className="about-document">
              <AboutItem id="about-purpose" title="Purpose">
                <p>
                  Polar Genome DB is an integrated genomics platform for organisms from polar and cold-environment ecosystems. It provides genome
                  browsing, gene and sequence search, functional annotation, comparative genomics, expression data visualization, and primer design
                  tools. The database supports the exploration of genetic diversity, adaptation mechanisms, and evolutionary relationships in Arctic,
                  Antarctic, and other cold-adapted species. Reference-based annotations and comparative analyses are provided to facilitate functional
                  interpretation across species.
                </p>
              </AboutItem>
              <AboutItem id="about-reference" title="Reference Genome and Annotation">
                <p>
                  The current species resources are based on in-house genome assemblies and annotations prepared for this project. These reference
                  genomes are not currently public NCBI reference genomes, and sequence identifiers, gene models, and annotation versions should be cited
                  or interpreted according to the project-provided release information.
                </p>
              </AboutItem>
              <AboutItem id="about-scope" title="Current Scope">
                <p>
                  The beta release currently supports two species, <em>Colobanthus quitensis</em> and <em>Deschampsia antarctica</em>. Available modules
                  include genome browsing, gene search, expression track visualization, BLAST/DIAMOND sequence search, similar gene lookup, GO annotation,
                  semantic GO grouping, and PCR primer design.
                </p>
                <p>
                  During the beta period, batch sizes and concurrent analysis jobs are limited, and available functions may change as the service is
                  refined.
                </p>
              </AboutItem>
              <AboutItem id="about-disclaimer" title="Disclaimer">
                <p>
                  Polar Genome DB is developed to support transparent and reproducible research on polar and cold-adapted organisms.
                  We are committed to improving data quality, analysis reliability, and documentation as the platform evolves.
                </p>
                <p>
                  This platform is provided for research use. Reference-based annotations, inferred functions, and computational predictions
                  represent supporting evidence rather than experimental validation. Users should independently assess results before drawing
                  biological conclusions, designing experiments, or publishing findings.
                </p>
                <p>
                  As a beta service, datasets, analysis methods, and results may change between releases. Please record the platform version
                  and relevant analysis settings when reporting results, and contact the team to report suspected errors.
                </p>
              </AboutItem>
              <AboutItem id="about-privacy" title="Privacy">
                <p>Operational access records support service reliability and security. Analysis inputs and results are processed separately to provide the requested tools.</p>
                <p><a href="/privacy">Read the Privacy Policy</a> for information about access records, analysis data, browser storage, and retention.</p>
              </AboutItem>
              <AboutItem id="about-team" title="Team & Contributors">
                <p className="about-team-intro">Bringing together platform development, polar genome resources, and scientific guidance.</p>
                <dl className="about-team-list">
                  {[
                    ['Janghyun Choi', 'Inha University', 'Platform design & development', 'Platform architecture, database and analysis pipelines, and frontend and backend development.'],
                    ['Seung Chul Shin', 'KOPRI', 'Genome resources & annotation', 'Genome resource preparation and annotation.'],
                    ['Jungeun Lee', 'KOPRI', 'Transcriptomics & interpretation', 'RNA-seq datasets and biological interpretation.'],
                    ['Jinmi Yoon', 'Inha University', 'Academic supervision', 'Research guidance and scientific advice.'],
                    ['Hyoungseok Lee', 'KOPRI', 'Project supervision', 'Overall project oversight and scientific direction.'],
                  ].map(([name, affiliation, role, contribution]) => (
                    <div className="about-team-member" key={name}>
                      <dt><span className="about-team-name">{name}<span className="about-team-degree">, Ph.D.</span></span><span className="about-team-affiliation">{affiliation}</span></dt>
                      <dd><span className="about-team-role">{role}</span><span className="about-team-contribution">{contribution}</span></dd>
                    </div>
                  ))}
                </dl>
              </AboutItem>
              <AboutItem id="about-source" title="Source Code">
                <p>
                  The Polar Genome DB source code is available through the project GitHub repository, bringing together the web application,
                  backend analysis services, and manifest-driven database build scripts. These components document how genome resources are
                  organized, processed, and presented within the platform.
                </p>
                <p>
                  The code supports inspection of the analysis workflows and the development of compatible multi-species resources.
                  Source code availability is separate from access to the in-house genome and RNA-seq datasets; dataset availability
                  and reference sources are described on the Resources page.
                </p>
                <a className="about-source-link" role="link" aria-disabled="true" tabIndex={0} title="GitHub repository">
                  <Github size={18} aria-hidden="true" />GitHub repository
                </a>
              </AboutItem>
            </div>
          </div>
        </section>
      )}

      {page === 'species' && (
        <SpeciesPage species={species} selected={selectedSpecies} onSelect={setSelectedSpeciesId}
          modules={modules} onOpenTool={(id) => openTool(id as ToolId)} speciesImages={speciesImages} />
      )}

      {page === 'tools' && (
        <section className="page-section tools-page-section">
          {/* Retain navigation state within a species, but isolate late responses after a species change. */}
          {Array.from(new Set<ToolId>([...retainedTools, activeTool])).map((workspaceTool) => <div key={`${workspaceTool}:${selectedSpeciesId}`} style={{ display: workspaceTool === activeTool ? undefined : 'none' }}>
          <ToolWorkspace
            activeTool={workspaceTool}
            modules={entryModules}
            onSelectTool={openTool}
            geneTarget={workspaceTool === 'gene' ? geneTarget : ''}
            basketInput={basketInput}
            onBasketInputApplied={() => setBasketInput(null)}
            onOpenCandidateGene={openCandidateGene}
            selectedSpecies={selectedSpecies}
            selectedSpeciesId={selectedSpeciesId}
            species={species}
            onSelectSpecies={(nextSpeciesId) => { setGeneTarget(''); setSelectedSpeciesId(nextSpeciesId); }}
          />
          </div>)}
        </section>
      )}

      {page === 'resources' && (
        <ResourcesPage />
      )}
      {page === 'privacy' && <PrivacyPage />}

      <footer className="site-footer">
        <span>© 2026 Polar Genome DB. All rights reserved.</span>
        <nav className="footer-links" aria-label="Footer navigation">
          <a href="/about">About</a>
          <span aria-hidden="true">·</span>
          <a href="/resources#data-availability">Data Availability</a>
          <span aria-hidden="true">·</span>
          <a href="/privacy">Privacy Policy</a>
          <span aria-hidden="true">·</span>
          <a href="/about#about-team">Contact</a>
        </nav>
      </footer>
    </main>
  );
}

const aboutSections = [
  ['purpose', 'Purpose'], ['reference', 'Reference data'],
  ['scope', 'Current scope'], ['disclaimer', 'Disclaimer'],
  ['privacy', 'Privacy'],
  ['team', 'Team & contributors'], ['source', 'Source code'],
] as const;

function AboutIndex() {
  return <SectionIndex sections={aboutSections} prefix="about-" className="about-index" label="About page sections" heading="On this page" />;
}

function AboutItem({ id, title, children }: { id: string; title: string; children: React.ReactNode }) {
  return (
    <section className="about-item" id={id}>
      <h2>{title}</h2>
      {children}
    </section>
  );
}

type ApiRecord = Record<string, any>;
type GeneSearchResult = {
  query: string;
  found: boolean;
  gene?: ApiRecord;
  transcripts?: ApiRecord[];
  proteins?: ApiRecord[];
  expression?: ApiRecord | null;
  similar?: ApiRecord | null;
  domains?: ApiRecord | null;
  error?: string;
};

function parseGeneIds(value: string) {
  return [...new Set(value.split(/[\s,;]+/).map((item) => item.trim()).filter(Boolean))];
}

async function optionalApiGet(path: string) {
  try {
    return await apiGet<ApiRecord>(path);
  } catch {
    return null;
  }
}

function ToolWorkspace({
  activeTool,
  geneTarget,
  basketInput,
  onBasketInputApplied,
  onOpenCandidateGene,
  modules,
  onSelectTool,
  selectedSpecies,
  selectedSpeciesId,
  species,
  onSelectSpecies,
}: {
  activeTool: ToolId;
  geneTarget: string;
  basketInput: { speciesId: string; geneIds: string[] } | null;
  onBasketInputApplied: () => void;
  onOpenCandidateGene: (speciesId: string, geneId: string) => void;
  modules: Module[];
  onSelectTool: (toolId: ToolId) => void;
  selectedSpecies: SpeciesSummary | undefined;
  selectedSpeciesId: string;
  species: SpeciesSummary[];
  onSelectSpecies: (speciesId: string) => void;
}) {
  const sequenceFileInputRef = useRef<HTMLInputElement | null>(null);
  const speciesOptions = species.length ? species : demoSpecies;
  const activeSpecies = selectedSpecies ?? speciesOptions.find((item) => item.species_id === selectedSpeciesId) ?? speciesOptions[0];
  const speciesId = activeSpecies?.species_id ?? selectedSpeciesId;
  const [geneInput, setGeneInput] = useState(defaultGeneForTool(activeTool, speciesId));
  const [geneSearchResults, setGeneSearchResults] = useState<GeneSearchResult[] | null>(null);
  const [geneData, setGeneData] = useState<ApiRecord | null>(null);
  const [geneExpressionData, setGeneExpressionData] = useState<ApiRecord | null>(null);
  const [geneSequenceData, setGeneSequenceData] = useState<ApiRecord | null>(null);
  const [geneDomainData, setGeneDomainData] = useState<ApiRecord | null>(null);
  const [geneMatrixData, setGeneMatrixData] = useState<ApiRecord | null>(null);
  const [geneSimilarById, setGeneSimilarById] = useState<Record<string, ApiRecord>>({});
  const [geneGoById, setGeneGoById] = useState<Record<string, ApiRecord>>({});
  const [viewerLocation, setViewerLocation] = useState<ApiRecord | null>(null);
  const [jbrowseData, setJbrowseData] = useState<ApiRecord | null>(null);
  const [similarData, setSimilarData] = useState<ApiRecord | null>(null);
  const [similarLookupData, setSimilarLookupData] = useState<ApiRecord | null>(null);
  const [similarResolutionData, setSimilarResolutionData] = useState<ApiRecord | null>(null);
  const [selectedSimilarProteinIds, setSelectedSimilarProteinIds] = useState<string[]>([]);
  const [similarJob, setSimilarJob] = useState<ApiRecord | null>(null);
  const [pairwiseJobs, setPairwiseJobs] = useState<Record<string, ApiRecord>>({});
  const [activePairwiseKey, setActivePairwiseKey] = useState('');
  const [similarProteinId, setSimilarProteinId] = useState('');
  const [goData, setGoData] = useState<ApiRecord | null>(null);
  const [semanticData, setSemanticData] = useState<ApiRecord | null>(null);
  const [goEnrichmentData, setGoEnrichmentData] = useState<ApiRecord | null>(null);
  const [goGeneInput, setGoGeneInput] = useState('');
  const [goBackgroundInput, setGoBackgroundInput] = useState('');
  const [goMinCount, setGoMinCount] = useState(3);
  const [goCorrection, setGoCorrection] = useState('benjamini-hochberg');
  const [goConfidence, setGoConfidence] = useState<Set<string>>(new Set(['high', 'medium']));
  const [searchProgram, setSearchProgram] = useState<SearchProgram>('blastn');
  const [sequenceQuery, setSequenceQuery] = useState('');
  const [searchJob, setSearchJob] = useState<ApiRecord | null>(null);
  const [primerFlank, setPrimerFlank] = useState(0);
  const [primerGenomicRegion, setPrimerGenomicRegion] = useState<'gene_body' | 'tss_window'>('gene_body');
  const [primerTssUpstream, setPrimerTssUpstream] = useState(1000);
  const [primerTssDownstream, setPrimerTssDownstream] = useState(300);
  const [primerProductRange, setPrimerProductRange] = useState('100-300');
  const [standardPrimerSettings, setStandardPrimerSettings] = useState({ lengthMin: 18, lengthMax: 24, tmMin: 57, tmMax: 63, gcMin: 55, gcMax: 65 });
  const [cloningPrimerSettings, setCloningPrimerSettings] = useState({ lengthMin: 20, lengthMax: 40, tmMin: 57, tmMax: 65, gcMin: 30, gcMax: 70 });
  const [primerSource, setPrimerSource] = useState<'gene' | 'sequence'>('gene');
  const [primerTemplateType, setPrimerTemplateType] = useState<'gdna' | 'cdna' | 'cdna_cloning'>('gdna');
  const isCloningInput = primerSource === 'gene' && primerTemplateType === 'cdna_cloning';
  const { lengthMin: primerLengthMin, lengthMax: primerLengthMax, tmMin: primerTmMin, tmMax: primerTmMax, gcMin: primerGcMin, gcMax: primerGcMax } = isCloningInput ? cloningPrimerSettings : standardPrimerSettings;
  const setPrimerSettings = isCloningInput ? setCloningPrimerSettings : setStandardPrimerSettings;
  const setPrimerLengthMin = (value: number) => setPrimerSettings((current) => ({ ...current, lengthMin: value }));
  const setPrimerLengthMax = (value: number) => setPrimerSettings((current) => ({ ...current, lengthMax: value }));
  const setPrimerTmMin = (value: number) => setPrimerSettings((current) => ({ ...current, tmMin: value }));
  const setPrimerTmMax = (value: number) => setPrimerSettings((current) => ({ ...current, tmMax: value }));
  const setPrimerGcMin = (value: number) => setPrimerSettings((current) => ({ ...current, gcMin: value }));
  const setPrimerGcMax = (value: number) => setPrimerSettings((current) => ({ ...current, gcMax: value }));
  const [primerSequenceInput, setPrimerSequenceInput] = useState('');
  const [primerSequenceError, setPrimerSequenceError] = useState('');
  const [primerJunctionMode, setPrimerJunctionMode] = useState<'no_junction' | 'allow_junction'>('no_junction');
  const [primerStopCodonMode, setPrimerStopCodonMode] = useState<'keep' | 'remove'>('keep');
  const [primerForwardTail, setPrimerForwardTail] = useState('');
  const [primerReverseTail, setPrimerReverseTail] = useState('');
  const [checkSpecificity, setCheckSpecificity] = useState(true);
  const [primerJob, setPrimerJob] = useState<ApiRecord | null>(null);
  const [isToolNavExpanded, setIsToolNavExpanded] = useState(false);
  const [busy, setBusy] = useState<string | null>(null);
  const [toolError, setToolError] = useState('');
  const geneSearchGeneration = useRef(0);

  useEffect(() => {
    geneSearchGeneration.current += 1;
    setBusy((current) => current === 'gene' ? null : current);
    return () => { geneSearchGeneration.current += 1; };
  }, [geneInput]);

  useEffect(() => {
    // Preserve an incoming basket selection when StrictMode replays mount effects.
    setGeneInput(activeTool === 'gene' && basketInput?.speciesId === speciesId
      ? basketInput.geneIds.join('\n')
      : defaultGeneForTool(activeTool, speciesId));
    setGeneSearchResults(null);
    setGeneData(null);
    setGeneExpressionData(null);
    setGeneSequenceData(null);
    setGeneDomainData(null);
    setGeneMatrixData(null);
    setGeneSimilarById({});
    setGeneGoById({});
    setViewerLocation(null);
    setJbrowseData(null);
    setSimilarData(null);
    setSimilarLookupData(null);
    setSimilarResolutionData(null);
    setSelectedSimilarProteinIds([]);
    setSimilarJob(null);
    setPairwiseJobs({});
    setSimilarProteinId('');
    setGoData(null);
    setSemanticData(null);
    setGoEnrichmentData(null);
    setGoGeneInput('');
    setGoBackgroundInput('');
    setGoMinCount(3);
    setGoCorrection('benjamini-hochberg');
    setGoConfidence(new Set(['high', 'medium']));
    setSearchJob(null);
    setPrimerJob(null);
    setPrimerFlank(0);
    setPrimerGenomicRegion('gene_body');
    setPrimerTssUpstream(1000);
    setPrimerTssDownstream(300);
    setPrimerProductRange('100-300');
    setPrimerLengthMin(18);
    setPrimerLengthMax(24);
    setPrimerTmMin(57);
    setPrimerTmMax(63);
    setPrimerGcMin(55);
    setPrimerGcMax(65);
    setPrimerSource('gene');
    setPrimerTemplateType('gdna');
    setPrimerSequenceInput('');
    setPrimerSequenceError('');
    setPrimerJunctionMode('no_junction');
    setPrimerStopCodonMode('keep');
    setPrimerForwardTail('');
    setPrimerReverseTail('');
    setCheckSpecificity(true);
    setToolError('');
  }, [activeTool, speciesId]);

  useEffect(() => {
    if (activeTool !== 'gene' || !geneTarget) return;
    setGeneInput(geneTarget);
    openGeneReport(geneTarget, true);
  }, [activeTool, speciesId, geneTarget]);

  const appliedBasketInput = useRef<typeof basketInput>(null);
  useEffect(() => {
    if (activeTool !== 'gene' || !basketInput || basketInput.speciesId !== speciesId || appliedBasketInput.current === basketInput) return;
    appliedBasketInput.current = basketInput;
    setGeneInput(basketInput.geneIds.join('\n'));
    setGeneSearchResults(null);
    setGeneData(null);
    onBasketInputApplied();
  }, [activeTool, speciesId, basketInput, onBasketInputApplied]);

  useEffect(() => {
    if (activeTool !== 'genome') return;
    const params = new URLSearchParams(window.location.search);
    const defaultGene = params.get('species') === speciesId && params.get('gene') ? params.get('gene')! : defaultGeneForTool('genome', speciesId);
    setGeneInput(defaultGene);
    let cancelled = false;
    setBusy('genome');
    setToolError('');
    Promise.all([
      apiGet<ApiRecord>(`/api/species/${speciesId}/viewer-location/${encodeURIComponent(defaultGene)}`),
      apiGet<ApiRecord>(`/api/species/${speciesId}/jbrowse-config`).catch(() => null),
    ])
      .then(([location, config]) => {
        if (cancelled) return;
        setViewerLocation(location);
        setJbrowseData(config);
      })
      .catch((error) => {
        if (cancelled) return;
        setToolError(error instanceof Error ? error.message : String(error));
      })
      .finally(() => {
        if (!cancelled) setBusy(null);
      });
    return () => {
      cancelled = true;
    };
  }, [activeTool, speciesId]);

  const selectedModule = modules.find((module) => module.id === activeTool) ?? modules[0];
  const geneId = geneInput.trim();

  async function withBusy(label: string, task: () => Promise<void>, isCurrent = () => true) {
    setBusy(label);
    setToolError('');
    try {
      await task();
    } catch (error) {
      if (isCurrent()) setToolError(error instanceof Error ? error.message : String(error));
    } finally {
      if (isCurrent()) setBusy(null);
    }
  }

  async function pollJob(path: string, setter: (payload: ApiRecord) => void) {
    for (let attempt = 0; attempt < 35; attempt += 1) {
      await new Promise((resolve) => window.setTimeout(resolve, attempt === 0 ? 250 : 900));
      const current = await apiGet<ApiRecord>(path);
      setter(current);
      if (current.status === 'completed' || current.status === 'failed') return current;
    }
    return null;
  }

  function requireGene({ allowMultiple = false }: { allowMultiple?: boolean } = {}) {
    const geneIds = parseGeneIds(geneInput);
    if (!geneIds.length) {
      throw new Error('Gene ID is required.');
    }
    if (!allowMultiple && geneIds.length > 1) {
      throw new Error('Use one Gene ID for this tool. Gene Search supports multi-gene reports.');
    }
    return encodeURIComponent(geneIds[0]);
  }

  function runGeneSearch() {
    const generation = ++geneSearchGeneration.current;
    const isCurrent = () => generation === geneSearchGeneration.current;
    void withBusy('gene', async () => {
      const geneIds = parseGeneIds(geneInput);
      if (!geneIds.length) {
        throw new Error('Gene ID is required.');
      }

      setGeneData(null);
      setGeneExpressionData(null);
      setGeneSequenceData(null);
      setGeneDomainData(null);
      setGeneMatrixData(null);
      setGeneSimilarById({});
      setGeneGoById({});
      setSimilarData(null);
      setSimilarLookupData(null);
      setGoData(null);
      setSemanticData(null);

      const entries = await Promise.all(geneIds.map(async (query): Promise<GeneSearchResult> => {
        const encodedGene = encodeURIComponent(query);
        try {
          const detail = await apiGet<ApiRecord>(`/api/species/${speciesId}/genes/${encodedGene}`);
          const [expression, similar, domains] = await Promise.all([
            optionalApiGet(`/api/species/${speciesId}/expression/${encodedGene}`),
            optionalApiGet(`/api/species/${speciesId}/similar/${encodedGene}?limit=1`),
            optionalApiGet(`/api/species/${speciesId}/domains/${encodedGene}?cached_only=true`),
          ]);
          return {
            query,
            found: true,
            gene: detail.gene,
            transcripts: Array.isArray(detail.transcripts) ? detail.transcripts : [],
            proteins: Array.isArray(detail.proteins) ? detail.proteins : [],
            expression,
            similar,
            domains,
          };
        } catch (error) {
          return {
            query,
            found: false,
            error: error instanceof Error ? error.message : String(error),
          };
        }
      }));

      if (!isCurrent()) return;
      setGeneSearchResults(entries);
      setGeneSimilarById(Object.fromEntries(entries.filter((entry) => entry.found && entry.gene && entry.similar).map((entry) => [entry.gene!.gene_id, entry.similar as ApiRecord])));
    }, isCurrent);
  }

  function fillGeneSearchExample() {
    setGeneInput(geneSearchExampleGenesForSpecies(speciesId).join('\n'));
    setGeneSearchResults(null);
    setGeneData(null);
    setGeneExpressionData(null);
    setGeneSequenceData(null);
    setGeneDomainData(null);
    setGeneMatrixData(null);
    setSimilarData(null);
    setGoData(null);
    setSemanticData(null);
    setToolError('');
  }

  function fillSimilarGeneExample() {
    setGeneInput(referenceContextForTaxonGroup(activeSpecies?.taxon_group ?? 'dicot').examples.join('\n'));
    setGeneData(null);
    setSimilarData(null);
    setSimilarLookupData(null);
    setSimilarResolutionData(null);
    setSelectedSimilarProteinIds([]);
    setSimilarJob(null);
    setPairwiseJobs({});
    setSimilarProteinId('');
    setToolError('');
  }

  function fillPrimerExample() {
    setGeneInput(primerExampleGeneForSpecies(speciesId));
    setPrimerSource('gene');
    setPrimerTemplateType('gdna');
    setPrimerSequenceInput('');
    setPrimerJunctionMode('no_junction');
    setPrimerStopCodonMode('keep');
    setPrimerForwardTail('');
    setPrimerReverseTail('');
    setPrimerFlank(0);
    setPrimerGenomicRegion('gene_body');
    setPrimerTssUpstream(1000);
    setPrimerTssDownstream(300);
    setPrimerProductRange('100-300');
    setPrimerLengthMin(18);
    setPrimerLengthMax(24);
    setPrimerTmMin(57);
    setPrimerTmMax(63);
    setPrimerGcMin(55);
    setPrimerGcMax(65);
    setCheckSpecificity(true);
    setPrimerJob(null);
    setToolError('');
  }

  function openGeneReport(geneIdToOpen: string, direct = false) {
    void withBusy('gene report', async () => {
      const encodedGene = encodeURIComponent(geneIdToOpen);
      const [detail, expression, sequence, domains, similar] = await Promise.all([
        apiGet<ApiRecord>(`/api/species/${speciesId}/genes/${encodedGene}`),
        optionalApiGet(`/api/species/${speciesId}/expression/${encodedGene}`),
        optionalApiGet(`/api/species/${speciesId}/genes/${encodedGene}/sequence`),
        optionalApiGet(`/api/species/${speciesId}/domains/${encodedGene}?cached_only=true`),
        optionalApiGet(`/api/species/${speciesId}/similar/${encodedGene}?limit=10`),
      ]);
      setGeneData(detail);
      setGeneExpressionData(expression);
      setGeneSequenceData(sequence);
      setGeneDomainData(domains);
      setSimilarData(similar);
      if (direct) setGeneSearchResults([{
        query: geneIdToOpen, found: true, gene: detail.gene,
        transcripts: detail.transcripts ?? [], proteins: detail.proteins ?? [],
        expression, similar, domains,
      }]);
      window.setTimeout(() => {
        document.querySelector('.selected-report')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }, 30);
    });
  }

  function runGenomeLookup() {
    void withBusy('genome', async () => {
      const encodedGene = requireGene({ allowMultiple: true });
      const [location, config] = await Promise.all([
        apiGet<ApiRecord>(`/api/species/${speciesId}/viewer-location/${encodedGene}`),
        apiGet<ApiRecord>(`/api/species/${speciesId}/jbrowse-config`).catch(() => null),
      ]);
      setViewerLocation(location);
      setJbrowseData(config);
    });
  }

  function runSimilarLookup() {
    void withBusy('similar', async () => {
      const queries = parseGeneIds(geneInput);
      if (!queries.length) throw new Error('Enter at least one reference gene, transcript ID, or symbol.');
      const result = await apiPost<ApiRecord>(`/api/species/${speciesId}/similar/resolve`, { queries });
      const groups = Array.isArray(result.queries) ? result.queries as ApiRecord[] : [];
      const representativeIsoformPattern = activeSpecies?.taxon_group === 'monocot' ? /-01$/ : /\.1$/;
      const representativeProteinIds = [...new Set(groups.flatMap((group) => {
        const matches = Array.isArray(group.matches) ? group.matches as ApiRecord[] : [];
        const representative = matches.find((item) => representativeIsoformPattern.test(String(item.reference_transcript_id ?? item.reference_protein_id ?? ''))) ?? matches[0];
        const proteinId = String(representative?.reference_protein_id ?? '');
        return proteinId ? [proteinId] : [];
      }))];
      setSimilarResolutionData(result);
      setSelectedSimilarProteinIds(representativeProteinIds);
      setSimilarLookupData(null);
      setSimilarJob(null);
      setPairwiseJobs({});
      setActivePairwiseKey('');
    });
  }

  function toggleSimilarProteinSelection(proteinId: string) {
    setSelectedSimilarProteinIds((current) => (
      current.includes(proteinId) ? current.filter((item) => item !== proteinId) : [...current, proteinId]
    ));
  }

  function selectAllSimilarProteins(proteinIds: string[]) {
    setSelectedSimilarProteinIds([...new Set(proteinIds)]);
  }

  function clearSimilarProteinSelection() {
    setSelectedSimilarProteinIds([]);
  }

  function runSimilarCandidateSearch() {
    void withBusy('similar candidates', async () => {
      if (!selectedSimilarProteinIds.length) {
        throw new Error('Select at least one reference protein isoform.');
      }
      setSimilarJob(null);
      setActivePairwiseKey('');
      const job = await apiPost<ApiRecord>(`/api/species/${speciesId}/similar/jobs`, {
        species_id: speciesId,
        reference_protein_ids: selectedSimilarProteinIds,
      });
      setSimilarJob(job);
      const completed = await pollJob(`/api/species/${speciesId}/similar/jobs/${job.job_id}`, setSimilarJob);
      if (completed?.status === 'failed') {
        throw new Error(String(completed.error ?? 'Similar gene search failed.'));
      }
    });
  }

  function runPairwiseDomainComparison(referenceProteinId: string, targetProteinId: string) {
    const key = `${referenceProteinId}:${targetProteinId}`;
    setActivePairwiseKey(key);
    void (async () => {
      try {
        setToolError('');
        const job = await apiPost<ApiRecord>(`/api/species/${speciesId}/similar/pairwise`, {
          species_id: speciesId,
          reference_protein_id: referenceProteinId,
          target_protein_id: targetProteinId,
        });
        setPairwiseJobs((current) => ({ ...current, [key]: job }));
        await pollJob(`/api/species/${speciesId}/similar/pairwise/${job.job_id}`, (current) => {
          setPairwiseJobs((jobs) => ({ ...jobs, [key]: current }));
        });
      } catch (error) {
        setToolError(error instanceof Error ? error.message : String(error));
      }
    })();
  }

  function selectSimilarityProtein(proteinId: string) {
    setSimilarProteinId(proteinId);
    void withBusy('similar', async () => {
      const encodedGene = requireGene();
      const result = await apiGet<ApiRecord>(`/api/species/${speciesId}/similar/${encodedGene}?protein_id=${encodeURIComponent(proteinId)}&limit=5`);
      setSimilarData(result);
    });
  }

  function runGoLookup() {
    void withBusy('go', async () => {
      const encodedGene = requireGene();
      const [goResult, semanticResult] = await Promise.all([
        apiGet<ApiRecord>(`/api/species/${speciesId}/go/${encodedGene}`),
        apiGet<ApiRecord>(`/api/species/${speciesId}/go/${encodedGene}/semantic`),
      ]);
      setGoData(goResult);
      setSemanticData(semanticResult);
    });
  }

  function runGoEnrichment() {
    void withBusy('GO enrichment', async () => {
      const geneIds = parseGeneIds(goGeneInput);
      const backgroundIds = parseGeneIds(goBackgroundInput);
      if (!geneIds.length) {
        throw new Error('Enter at least one gene ID.');
      }
      if (geneIds.length > 1000) {
        throw new Error('Gene list accepts up to 1,000 unique IDs.');
      }
      if (backgroundIds.length > 100000) {
        throw new Error('Background list accepts up to 100,000 unique IDs.');
      }
      if (!goConfidence.size) {
        throw new Error('Select at least one homolog confidence level.');
      }
      setGoEnrichmentData(null);
      const result = await apiPost<ApiRecord>(`/api/species/${speciesId}/go/enrichment`, {
        genes: goGeneInput,
        background: goBackgroundInput,
        min_count: goMinCount,
        correction: goCorrection,
        confidence: Array.from(goConfidence),
        max_similarity_rank: 5,
      });
      setGoEnrichmentData(result);
    });
  }

  function runSequenceSearch() {
    void withBusy('blast', async () => {
      if (!sequenceQuery.trim()) {
        throw new Error('Query sequence is required.');
      }
      const queryCount = countSequenceQueries(sequenceQuery);
      const queryLimit = maxSequenceQueriesForProgram(searchProgram);
      if (queryCount > queryLimit) {
        throw new Error(`${searchProgram.toUpperCase()} accepts up to ${queryLimit} FASTA sequences; ${queryCount} were provided.`);
      }
      setSearchJob(null);
      const endpoint = searchProgram === 'diamond' ? '/api/search/diamond' : '/api/search/blast';
      const queued = await apiPost<ApiRecord>(endpoint, {
        species_id: speciesId,
        program: searchProgram,
        query: sequenceQuery,
        max_hits: searchProgram === 'diamond' ? 1 : 10,
      });
      setSearchJob(queued);
      await pollJob(`/api/search/jobs/${queued.job_id}`, setSearchJob);
    });
  }

  function fillSequenceSearchExample() {
    void (async () => {
      try {
        const records = await loadSequenceSearchExamples(speciesId, searchProgram);
        const labelPrefix = searchProgram === 'blastn' ? 'Gene' : 'Protein';
        setSequenceQuery(records.map((record, index) => (
          `>${labelPrefix}_${alphabeticalLabel(index)}\n${String(record.sequence ?? '').replace(/\s+/g, '').slice(0, searchProgram === 'blastn' ? 400 : 100)}`
        )).join('\n'));
        setSearchJob(null);
        setToolError('');
      } catch (error) {
        setToolError(error instanceof Error ? error.message : 'Example sequences could not be loaded.');
      }
    })();
  }

  useEffect(() => {
    if (activeTool !== 'blast') return;
    void loadSequenceSearchExamples(speciesId, searchProgram).catch(() => undefined);
  }, [activeTool, searchProgram, speciesId]);

  function uploadSequenceFasta(event: React.ChangeEvent<HTMLInputElement>) {
    const file = event.currentTarget.files?.[0];
    event.currentTarget.value = '';
    if (!file) return;
    if (!/\.(fa|fasta|fna|faa)$/i.test(file.name)) {
      setToolError('Upload a FASTA file with a .fa, .fasta, .fna, or .faa extension.');
      return;
    }
    void file.text().then((contents) => {
      if (!contents.trim()) {
        setToolError('The uploaded FASTA file is empty.');
        return;
      }
      setSequenceQuery(contents.trim());
      setSearchJob(null);
      setToolError('');
    }).catch(() => setToolError('The FASTA file could not be read.'));
  }

  function changeSequenceProgram(program: SearchProgram) {
    const isProteinProgram = (value: SearchProgram) => value !== 'blastn';
    if (isProteinProgram(program) !== isProteinProgram(searchProgram)) {
      setSequenceQuery('');
    }
    setSearchProgram(program);
    setSearchJob(null);
    setToolError('');
  }

  function runPrimerDesign(transcriptId?: string) {
    if (primerSource === 'sequence') {
      const sequenceBody = primerSequenceInput
        .split(/\r?\n/)
        .filter((line) => !line.trim().startsWith('>'))
        .join('')
        .replace(/\s+/g, '');
      if (sequenceBody && /[^ACGTacgt]/.test(sequenceBody)) {
        setPrimerSequenceError('Sequence must contain only A, C, G, and T.');
        setToolError('');
        return;
      }
    }
    setPrimerSequenceError('');
    void withBusy('primer', async () => {
      if (primerSource === 'gene' && !geneInput.trim()) {
        throw new Error('Enter a gene ID and select Resolve before designing primers.');
      }
      if (primerSource === 'sequence' && !primerSequenceInput.trim()) {
        throw new Error('Paste a DNA sequence before designing primers.');
      }
      const queued = await apiPost<ApiRecord>('/api/primer/design', normalizePrimerRequest({
        species_id: speciesId,
        source: primerSource,
        template_type: primerTemplateType,
        gene_id: primerSource === 'gene' ? geneInput.trim() : undefined,
        transcript_id: primerSource === 'gene' && primerTemplateType !== 'gdna' ? transcriptId : undefined,
        sequence: primerSource === 'sequence' ? primerSequenceInput : undefined,
        flank: primerFlank,
        genomic_region: primerGenomicRegion,
        tss_upstream: primerTssUpstream,
        tss_downstream: primerTssDownstream,
        junction_mode: primerJunctionMode,
        stop_codon_mode: primerStopCodonMode,
        forward_tail: primerForwardTail,
        reverse_tail: primerReverseTail,
        product_size_range: primerProductRange,
        primer_min_size: primerLengthMin,
        primer_max_size: primerLengthMax,
        primer_min_tm: primerTmMin,
        primer_max_tm: primerTmMax,
        primer_min_gc: primerGcMin,
        primer_max_gc: primerGcMax,
        check_specificity: checkSpecificity,
      }));
      setPrimerJob(queued);
      await pollJob(`/api/primer/jobs/${queued.job_id}`, setPrimerJob);
    });
  }

  const gene = geneData?.gene;
  const transcripts = Array.isArray(geneData?.transcripts) ? geneData.transcripts : [];
  const transcriptModels = Array.isArray(geneData?.transcript_models) ? geneData.transcript_models : [];
  const proteins = Array.isArray(geneData?.proteins) ? geneData.proteins : [];
  const features = Array.isArray(geneData?.features) ? geneData.features : [];
  const searchHits = Array.isArray(searchJob?.results) ? searchJob.results : [];
  const sequenceQueryLimit = maxSequenceQueriesForProgram(searchProgram);
  const sequenceQueryCount = countSequenceQueries(sequenceQuery);
  const sequenceQueryLimitExceeded = sequenceQueryCount > sequenceQueryLimit;
  const similarHits = Array.isArray(similarData?.hits) ? similarData.hits : [];
  const similarResolutionGroups = Array.isArray(similarResolutionData?.queries) ? similarResolutionData.queries as ApiRecord[] : [];
  const similarProteinOptions = [...proteins].sort((a, b) => String(a.transcript_id ?? '').localeCompare(String(b.transcript_id ?? ''), undefined, { numeric: true, sensitivity: 'base' }));
  const activeSimilarProtein = similarProteinOptions.find((protein) => protein.protein_id === similarProteinId) ?? similarProteinOptions[0];
  const goTerms = Array.isArray(goData?.go_terms) ? goData.go_terms : [];
  const semanticClusters = Array.isArray(semanticData?.clusters) ? semanticData.clusters : [];
  const primers = Array.isArray(primerJob?.results?.primers) ? primerJob.results.primers : [];

  return (
    <div className="tools-workspace">
      <div className={isToolNavExpanded ? 'tool-shell nav-expanded' : 'tool-shell nav-collapsed'}>
        <aside
          className={isToolNavExpanded ? 'tool-sidebar expanded' : 'tool-sidebar collapsed'}
          aria-label="Tool modules"
          onMouseEnter={() => setIsToolNavExpanded(true)}
          onMouseLeave={() => setIsToolNavExpanded(false)}
          onFocus={() => setIsToolNavExpanded(true)}
          onBlur={(event) => {
            if (!event.currentTarget.contains(event.relatedTarget)) {
              setIsToolNavExpanded(false);
            }
          }}
        >
          {modules.map((module) => (
            <button
              type="button"
              key={module.id}
              className={module.id === activeTool ? 'active' : ''}
              onClick={() => {
                onSelectTool(module.id as ToolId);
              }}
              title={module.title}
            >
              {module.icon}
              <span>
                <strong>{module.title}</strong>
                <small>{module.scope}</small>
              </span>
            </button>
          ))}
        </aside>

        <section className={activeTool === 'gene' ? 'tool-panel gene-tool-panel' : activeTool === 'similar' ? 'tool-panel similar-tool-panel' : activeTool === 'go' ? 'tool-panel go-tool-panel' : activeTool === 'primer' ? 'tool-panel primer-tool-panel' : activeTool === 'blast' ? 'tool-panel sequence-tool-panel' : 'tool-panel'}>
          <div className="tool-panel-head">
            <div>
              <h3>{selectedModule.title}</h3>
              <p>{moduleDescriptions[selectedModule.id] ?? selectedModule.scope}</p>
            </div>
            <strong>{speciesId}</strong>
          </div>

          {toolError && activeTool !== 'go' && activeTool !== 'blast' && <div className="tool-message error">{toolError}</div>}
          {busy && activeTool !== 'gene' && activeTool !== 'similar' && activeTool !== 'go' && activeTool !== 'primer' && activeTool !== 'blast' && <div className="tool-message">Running {busy} request...</div>}

          {activeTool === 'genome' && (
            <div className="tool-body wide">
              <GenomeControlPanel
                geneInput={geneInput}
                onGeneInput={setGeneInput}
                onRun={runGenomeLookup}
                onSelectSpecies={onSelectSpecies}
                speciesId={speciesId}
                speciesOptions={speciesOptions}
              />
              {viewerLocation ? (
                <EmbeddedGenomeBrowser config={jbrowseData?.config} location={viewerLocation.location} />
              ) : !toolError ? (
                <GenomeLoading />
              ) : (
                <div className="browser-empty">
                  <GenomePreview />
                  <EmptyResult text="Select a species and resolve a gene to load JBrowse with annotation and bigWig tracks." />
                </div>
              )}
            </div>
          )}

          {activeTool === 'gene' && (
            <div className="tool-body constrained gene-tool-body">
              <GeneSearchPanel geneInput={geneInput} onGeneInput={setGeneInput} speciesId={speciesId} speciesOptions={speciesOptions} onSelectSpecies={onSelectSpecies} />
              <div className="tool-actions">
                <button type="button" className="gene-search-submit" onClick={runGeneSearch}>Search</button>
                <button type="button" className="secondary gene-example-button" onClick={fillGeneSearchExample}>Example</button>
              </div>
              {busy === 'gene' ? (
                <div className="gene-request-status" role="status" aria-live="polite">
                  <LoaderCircle size={17} aria-hidden="true" />
                  <span>Searching gene records and integrated evidence...</span>
                </div>
              ) : geneSearchResults?.some((item) => !item.found) ? (
                <div className="gene-request-status gene-request-notice" role="status" aria-live="polite">
                  {geneSearchResults.filter((item) => !item.found).map((item) => (
                    <span key={item.query}>{/^404\b/.test(item.error ?? '') ? `Gene ${item.query} could not be found.` : `Could not retrieve gene ${item.query}. Please try again.`}</span>
                  ))}
                </div>
              ) : null}
              {geneSearchResults ? (
                <GeneSearchResultsPanel
                  speciesId={speciesId}
                  results={geneSearchResults}
                  selectedGeneId={gene?.gene_id}
                  onOpenReport={openGeneReport}
                />
              ) : null}
              {gene ? (
                <section className="selected-report">
                  <GeneReport
                    domainData={geneDomainData}
                    expressionData={geneExpressionData}
                    features={features}
                    gene={gene}
                    proteins={proteins}
                    sequenceData={geneSequenceData}
                    similarData={similarData}
                    speciesId={speciesId}
                    transcriptModels={transcriptModels}
                    transcripts={transcripts}
                  />
                </section>
              ) : null}
            </div>
          )}

          {activeTool === 'blast' && (
            <div className="tool-body constrained sequence-tool-body">
              <SequenceSearchPanel
                speciesId={speciesId}
                speciesOptions={speciesOptions}
                onSelectSpecies={onSelectSpecies}
                program={searchProgram}
                query={sequenceQuery}
                inputError={toolError}
                onProgramChange={changeSequenceProgram}
                onQueryChange={(query) => { setSequenceQuery(query); setSearchJob(null); setToolError(''); }}
              />
              <div className="tool-actions">
                <button type="button" className="gene-search-submit" onClick={runSequenceSearch} disabled={busy === 'blast' || sequenceQueryLimitExceeded}>{busy === 'blast' ? <><LoaderCircle size={18} aria-hidden="true" /> Searching</> : 'Search'}</button>
                <input ref={sequenceFileInputRef} className="sequence-file-input" type="file" accept=".fa,.fasta,.fna,.faa,text/x-fasta" onChange={uploadSequenceFasta} />
                <button type="button" className="secondary sequence-upload-button" onClick={() => sequenceFileInputRef.current?.click()} disabled={busy === 'blast'}><Upload size={16} aria-hidden="true" /> Upload FASTA</button>
                <button type="button" className="secondary gene-example-button" onClick={fillSequenceSearchExample} disabled={busy === 'blast'}>Example</button>
              </div>
              {busy === 'blast' ? <div className="gene-request-status" role="status" aria-live="polite"><LoaderCircle size={17} aria-hidden="true" /><span>Searching the selected species runtime index...</span></div> : null}
              {searchJob ? (
                <SequenceSearchResults speciesId={String(searchJob.species_id ?? speciesId)} onOpenGene={(geneId) => onOpenCandidateGene(String(searchJob.species_id ?? speciesId), geneId)} jobId={String(searchJob.job_id ?? '')} program={searchProgram} status={String(searchJob.status ?? '')} userError={String(searchJob.user_error ?? '')} submittedQueries={sequenceQueryCount} queryIds={sequenceQueryIds(sequenceQuery)} queryRecords={sequenceQueryRecords(sequenceQuery)} hits={searchHits} />
              ) : null}
            </div>
          )}

          {activeTool === 'primer' && (
            <div className="tool-body constrained primer-tool-body">
              <PrimerDesignPanel
                geneInput={geneInput}
                speciesId={speciesId}
                speciesOptions={speciesOptions}
                flank={primerFlank}
                genomicRegion={primerGenomicRegion}
                tssUpstream={primerTssUpstream}
                tssDownstream={primerTssDownstream}
                productRange={primerProductRange}
                primerLengthMin={primerLengthMin}
                primerLengthMax={primerLengthMax}
                primerTmMin={primerTmMin}
                primerTmMax={primerTmMax}
                primerGcMin={primerGcMin}
                primerGcMax={primerGcMax}
                source={primerSource}
                templateType={primerTemplateType}
                sequenceInput={primerSequenceInput}
                sequenceError={primerSequenceError}
                junctionMode={primerJunctionMode}
                stopCodonMode={primerStopCodonMode}
                forwardTail={primerForwardTail}
                reverseTail={primerReverseTail}
                checkSpecificity={checkSpecificity}
                isRunning={busy === 'primer'}
                onGeneInput={setGeneInput}
                onSelectSpecies={onSelectSpecies}
                onFlank={setPrimerFlank}
                onGenomicRegion={setPrimerGenomicRegion}
                onTssUpstream={setPrimerTssUpstream}
                onTssDownstream={setPrimerTssDownstream}
                onProductRange={setPrimerProductRange}
                onPrimerLengthMin={setPrimerLengthMin}
                onPrimerLengthMax={setPrimerLengthMax}
                onPrimerTmMin={setPrimerTmMin}
                onPrimerTmMax={setPrimerTmMax}
                onPrimerGcMin={setPrimerGcMin}
                onPrimerGcMax={setPrimerGcMax}
                onSource={setPrimerSource}
                onTemplateType={setPrimerTemplateType}
                onSequenceInput={(value) => {
                  setPrimerSequenceInput(value);
                  setPrimerSequenceError('');
                }}
                onJunctionMode={setPrimerJunctionMode}
                onStopCodonMode={setPrimerStopCodonMode}
                onForwardTail={setPrimerForwardTail}
                onReverseTail={setPrimerReverseTail}
                onCheckSpecificity={setCheckSpecificity}
                onExample={fillPrimerExample}
                onRun={runPrimerDesign}
              />
              {primerJob ? (
                <PrimerTable job={primerJob} primers={primers} />
              ) : null}
            </div>
          )}

          {activeTool === 'similar' && (
            <div className="tool-body constrained similar-tool-body">
              <SimilarGeneSearchPanel
                geneInput={geneInput}
                onGeneInput={setGeneInput}
                speciesId={speciesId}
                taxonGroup={activeSpecies?.taxon_group ?? 'dicot'}
                speciesOptions={speciesOptions}
                onSelectSpecies={onSelectSpecies}
              />
              <div className="tool-actions">
                <button type="button" className="gene-search-submit" onClick={runSimilarLookup}>Find similar genes</button>
                <button type="button" className="secondary gene-example-button" onClick={fillSimilarGeneExample}>Example</button>
              </div>
              {similarResolutionData ? (
                <ReferenceIsoformSelectionPanel
                  groups={similarResolutionGroups}
                  selectedProteinIds={selectedSimilarProteinIds}
                  onToggleProtein={toggleSimilarProteinSelection}
                  onSelectAll={selectAllSimilarProteins}
                  onClearSelection={clearSimilarProteinSelection}
                  onRun={runSimilarCandidateSearch}
                  isRunning={busy === 'similar candidates'}
                />
              ) : null}
              {similarJob && (similarJob.status === 'queued' || similarJob.status === 'running') ? (
                <div className="similar-run-status" role="status" aria-live="polite">
                  <LoaderCircle size={18} aria-hidden="true" />
                  <span>Running initial DIAMOND rank search...</span>
                </div>
              ) : null}
              {similarJob?.status === 'completed' && similarJob.results ? <SimilarCandidateResultsPanel speciesId={speciesId} similarJobId={String(similarJob.job_id)} results={similarJob.results} pairwiseJobs={pairwiseJobs} activePairwiseKey={activePairwiseKey} onRunPairwise={runPairwiseDomainComparison} onOpenGene={(geneId) => onOpenCandidateGene(speciesId, geneId)} /> : null}
            </div>
          )}

          {activeTool === 'go' && (
            <div className="tool-body constrained go-tool-body">
              <GoAnalysisInputPanel
                geneInput={goGeneInput}
                backgroundInput={goBackgroundInput}
                minCount={goMinCount}
                correction={goCorrection}
                confidence={goConfidence}
                speciesId={speciesId}
                taxonGroup={activeSpecies?.taxon_group ?? 'dicot'}
                speciesOptions={speciesOptions}
                onGeneInput={setGoGeneInput}
                onBackgroundInput={setGoBackgroundInput}
                onMinCount={setGoMinCount}
                onCorrection={setGoCorrection}
                onConfidence={setGoConfidence}
                onSelectSpecies={onSelectSpecies}
                onRun={runGoEnrichment}
                isRunning={busy === 'GO enrichment'}
              />
              {toolError ? <p className="go-input-error" role="alert">{toolError}</p> : null}
              {busy === 'GO enrichment' ? (
                <div className="gene-request-status go-request-status" role="status" aria-live="polite">
                  <LoaderCircle size={17} aria-hidden="true" />
                  <span>Resolving reference-derived GO annotations and testing enriched terms...</span>
                </div>
              ) : null}
              {goEnrichmentData ? <GoEnrichmentResults result={goEnrichmentData} /> : null}
            </div>
          )}
        </section>
      </div>
    </div>
  );
}

function SequenceSearchPanel({
  speciesId,
  speciesOptions,
  onSelectSpecies,
  program,
  query,
  inputError,
  onProgramChange,
  onQueryChange,
}: {
  speciesId: string;
  speciesOptions: SpeciesSummary[];
  onSelectSpecies: (speciesId: string) => void;
  program: SearchProgram;
  query: string;
  inputError: string;
  onProgramChange: (program: SearchProgram) => void;
  onQueryChange: (query: string) => void;
}) {
  const programDetails = {
    blastn: {
      title: 'BLASTN',
      subtitle: 'Nucleotide sequence',
      placeholder: 'Paste a nucleotide sequence (A, C, G, T, or IUPAC ambiguity codes).',
      notes: ['Search the selected genome assembly with the current BLASTN default parameters.'],
    },
    blastp: {
      title: 'BLASTP',
      subtitle: 'Protein sequence',
      placeholder: 'Paste an amino-acid sequence in one-letter code.',
      notes: ['Search the selected protein set with the current BLASTP default parameters.'],
    },
    diamond: {
      title: 'DIAMOND',
      subtitle: 'Batch protein search',
      placeholder: 'Paste an amino-acid sequence in one-letter code.',
      notes: ['High-throughput protein sequence analysis using DIAMOND default parameters.'],
    },
  } as const;
  const current = programDetails[program];
  const queryLimit = maxSequenceQueriesForProgram(program);
  const queryCount = countSequenceQueries(query);
  const isLimitExceeded = queryCount > queryLimit;

  return (
    <div className="sequence-search-panel">
      <CustomSpeciesSelect speciesId={speciesId} speciesOptions={speciesOptions} onSelectSpecies={onSelectSpecies} />
      <div className="sequence-program-field">
        <span>Search type</span>
        <div className="sequence-program-tabs" role="group" aria-label="Sequence search type">
          {(Object.keys(programDetails) as Array<keyof typeof programDetails>).map((item) => (
            <button key={item} type="button" className={program === item ? 'active' : ''} onClick={() => onProgramChange(item)}>
              <strong>{programDetails[item].title}</strong>
              <small>{programDetails[item].subtitle}</small>
            </button>
          ))}
        </div>
      </div>
      <ul className="sequence-search-notes">
        {current.notes.map((note) => <li key={note}>{note}</li>)}
        <li><strong>Up to {queryLimit} FASTA queries.</strong></li>
        <li>Use one <code>&gt;</code> header per sequence for a multi-query FASTA input.</li>
      </ul>
      <label className="sequence-query-field">
        <span>Query sequence</span>
        <textarea
          value={query}
          onChange={(event) => onQueryChange(event.target.value)}
          placeholder={current.placeholder}
          rows={10}
          spellCheck={false}
          aria-invalid={Boolean(inputError || isLimitExceeded)}
          aria-describedby={inputError ? 'sequence-input-error' : undefined}
        />
        {inputError && <p id="sequence-input-error" className="sequence-input-error" role="alert">{inputError}</p>}
        {isLimitExceeded && (
          <p className="sequence-query-limit-error" role="alert">
            {queryCount} FASTA headers detected. {current.title} accepts at most {queryLimit} sequences.
          </p>
        )}
      </label>
    </div>
  );
}

function SequenceGeneLink({ speciesId, geneId, onOpenGene }: { speciesId: string; geneId: string; onOpenGene: (geneId: string) => void }) {
  return <span className="sequence-gene-actions">
    <a href={`/tools?tool=gene&species=${encodeURIComponent(speciesId)}&gene=${encodeURIComponent(geneId)}`} onClick={(event) => {
      if (event.button !== 0 || event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
      event.preventDefault();
      onOpenGene(geneId);
    }}>{geneId}</a>
    <BasketButton speciesId={speciesId} geneId={geneId} className="sequence-gene-basket" />
  </span>;
}

function SequenceSearchResults({ speciesId, onOpenGene, jobId, program, status, userError, submittedQueries, queryIds, queryRecords, hits }: { speciesId: string; onOpenGene: (geneId: string) => void; jobId: string; program: SearchProgram; status: string; userError: string; submittedQueries: number; queryIds: string[]; queryRecords: Array<{ id: string; sequence: string }>; hits: ApiRecord[] }) {
  const groupedHits = new Map<string, ApiRecord[]>();
  const perQueryLimit = program === 'diamond' ? 1 : 10;
  for (const hit of hits) {
    const queryId = String(hit.qseqid ?? 'query');
    const current = groupedHits.get(queryId) ?? [];
    if (current.length < perQueryLimit) current.push(hit);
    groupedHits.set(queryId, current);
  }
  const groups = Array.from(groupedHits.entries());
  const queryLengths = new Map(queryRecords.map((record) => [record.id, record.sequence.length]));

  if (status === 'failed') {
    return <section className="sequence-results" aria-live="polite">
      <div className="sequence-results-heading"><h3>{program === 'diamond' ? 'DIAMOND Batch Results' : 'Sequence Search Results'}</h3></div>
      <p className="sequence-input-error" role="alert">{userError || 'Search could not be completed. Please retry; if it persists, contact the platform administrator.'}</p>
    </section>;
  }

  if (program === 'diamond') {
    return <DiamondBatchResults speciesId={speciesId} onOpenGene={onOpenGene} jobId={jobId} status={status} submittedQueries={submittedQueries} queryIds={queryIds} rawHits={hits} hits={groups.map(([queryId, queryHits]) => ({ queryId, hit: queryHits[0] }))} />;
  }

  return (
    <section className="sequence-results" aria-live="polite">
      <div className="sequence-results-heading">
        <div><h3>Sequence Search Results</h3><p>BLAST displays up to 10 hits per query at E-value &le; 1e-5.</p></div>
        <SequenceResultsExport jobId={jobId} program={program} queryIds={queryIds} hits={hits} />
      </div>
      {status !== 'completed' ? <EmptyResult text={`Search job status: ${status || 'queued'}.`} /> : null}
      {status === 'completed' && !groups.length ? <EmptyResult text="No matching genes were returned for this query set." /> : null}
      {status === 'completed' && groups.map(([queryId, queryHits], index) => (
        <SequenceQueryResults key={queryId} speciesId={speciesId} onOpenGene={onOpenGene} jobId={jobId} program={program} queryId={queryId} queryLength={queryLengths.get(queryId) ?? 0} queryNumber={index + 1} hits={queryHits} />
      ))}
    </section>
  );
}

function DiamondBatchResults({ speciesId, onOpenGene, jobId, status, submittedQueries, queryIds, rawHits, hits }: { speciesId: string; onOpenGene: (geneId: string) => void; jobId: string; status: string; submittedQueries: number; queryIds: string[]; rawHits: ApiRecord[]; hits: Array<{ queryId: string; hit: ApiRecord }> }) {
  const [page, setPage] = useState(0);
  const pageSize = 25;
  const pageCount = Math.max(1, Math.ceil(hits.length / pageSize));
  const activePage = Math.min(page, pageCount - 1);
  const visibleHits = hits.slice(activePage * pageSize, (activePage + 1) * pageSize);
  const noHitCount = Math.max(0, submittedQueries - hits.length);

  useEffect(() => setPage(0), [hits.length, submittedQueries]);

  return (
    <section className="sequence-results diamond-batch-results" aria-live="polite">
      <div className="sequence-results-heading">
        <div><h3>DIAMOND Batch Results</h3><p>Top matching gene per protein query using DIAMOND default parameters.</p></div>
        <SequenceResultsExport jobId={jobId} program="diamond" queryIds={queryIds} hits={rawHits} />
      </div>
      {status !== 'completed' ? <EmptyResult text={`Search job status: ${status || 'queued'}.`} /> : null}
      {status === 'completed' ? (
        <>
          <div className="diamond-batch-summary" aria-label="DIAMOND batch summary">
            <span><strong>{submittedQueries}</strong> submitted</span>
            <span><strong>{hits.length}</strong> matched</span>
            <span><strong>{noHitCount}</strong> no hit</span>
          </div>
          {hits.length ? (
            <>
              <div className="table-wrap sequence-result-table-wrap">
                <table className="multi-table sequence-result-table diamond-result-table" data-species={speciesId}>
                  <colgroup>
                    <col className="diamond-col-query" />
                    <col className="diamond-col-id" />
                    <col className="diamond-col-id" />
                    <col className="diamond-col-id" />
                    <col className="diamond-col-identity" />
                    <col className="diamond-col-alignment" />
                    <col className="diamond-col-evalue" />
                    <col className="diamond-col-score" />
                    <col className="diamond-col-raw" />
                  </colgroup>
                  <thead>
                    <tr>
                      <th>Query</th>
                      <th>Gene</th>
                      <th>Transcript</th>
                      <th>Protein</th>
                      <th>Identity (%)</th>
                      <th>Alignment</th>
                      <th>E-value</th>
                      <th>Bit score</th>
                      <th className="sequence-raw-header">Raw</th>
                    </tr>
                  </thead>
                  <tbody>
                    {visibleHits.map(({ queryId, hit }) => {
                      const annotation = (hit.annotation ?? {}) as ApiRecord;
                      const geneId = annotation.gene_id ? String(annotation.gene_id) : '';
                      return (
                        <tr key={queryId}>
                          <td className="diamond-query-id">{queryId}</td>
                          <td className={geneId ? 'sequence-result-gene' : ''}>{geneId ? <SequenceGeneLink speciesId={speciesId} geneId={geneId} onOpenGene={onOpenGene} /> : '-'}</td>
                          <td>{annotation.transcript_id ?? '-'}</td>
                          <td>{annotation.protein_id ?? hit.sseqid ?? '-'}</td>
                          <td>{Number(hit.pident ?? 0).toFixed(1)}</td>
                          <td>{Number(hit.length ?? 0).toLocaleString()} aa</td>
                          <td>{formatCell(hit.evalue, 'evalue')}</td>
                          <td>{hit.bitscore != null && Number.isFinite(Number(hit.bitscore)) ? Math.round(Number(hit.bitscore)).toLocaleString() : '-'}</td>
                          <td className="sequence-raw-cell"><SequenceQueryResultsExport jobId={jobId} program="diamond" queryId={queryId} /></td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>
              {pageCount > 1 ? (
                <div className="compact-table-pagination">
                  <span>{activePage * pageSize + 1}-{Math.min((activePage + 1) * pageSize, hits.length)} of {hits.length}</span>
                  <div>
                    <button type="button" aria-label="Previous DIAMOND result page" disabled={activePage === 0} onClick={() => setPage((value) => Math.max(0, value - 1))}><ChevronLeft size={16} /></button>
                    <span>{activePage + 1} / {pageCount}</span>
                    <button type="button" aria-label="Next DIAMOND result page" disabled={activePage >= pageCount - 1} onClick={() => setPage((value) => Math.min(pageCount - 1, value + 1))}><ChevronRight size={16} /></button>
                  </div>
                </div>
              ) : null}
            </>
          ) : <EmptyResult text="No matching genes were returned for this protein query set." />}
        </>
      ) : null}
    </section>
  );
}

function SequenceResultsExport({ jobId, program, queryIds, hits }: { jobId: string; program: SearchProgram; queryIds: string[]; hits: ApiRecord[] }) {
  const downloadAllResults = () => {
    if (jobId) {
      window.location.assign(`${API_BASE}/api/search/jobs/${encodeURIComponent(jobId)}/raw-results`);
      return;
    }
    const fields = ['query_id', 'rank', 'subject_id', 'gene_id', 'transcript_id', 'protein_id', 'genome_coordinate', 'identity_percent', 'alignment_length', 'alignment_unit', 'evalue', 'bit_score'];
    const hitsByQuery = new Map<string, ApiRecord[]>();
    for (const hit of hits) {
      const queryId = String(hit.qseqid ?? 'query');
      hitsByQuery.set(queryId, [...(hitsByQuery.get(queryId) ?? []), hit]);
    }
    const orderedQueryIds = queryIds.length ? queryIds : Array.from(hitsByQuery.keys());
    const sanitize = (value: unknown) => String(value ?? '').replace(/[\t\r\n]+/g, ' ').trim();
    const lines = orderedQueryIds.flatMap((queryId) => {
      const queryHits = hitsByQuery.get(queryId) ?? [];
      if (!queryHits.length) return [[queryId, '', '', '', '', '', '', '', '', program === 'blastn' ? 'bp' : 'aa', '', '']];
      return queryHits.map((hit, index) => {
        const annotation = (hit.annotation ?? {}) as ApiRecord;
        const coordinate = annotation.seqid && annotation.start && annotation.end
          ? `${annotation.seqid}:${annotation.start}-${annotation.end}${annotation.strand ? ` (${annotation.strand})` : ''}`
          : '';
        return [
          queryId,
          index + 1,
          hit.sseqid,
          annotation.gene_id,
          annotation.transcript_id,
          annotation.protein_id,
          coordinate,
          hit.pident,
          hit.length,
          program === 'blastn' ? 'bp' : 'aa',
          hit.evalue,
          hit.bitscore,
        ];
      });
    });
    const tsv = [fields, ...lines].map((row) => row.map(sanitize).join('\t')).join('\n');
    const blob = new Blob([`${tsv}\n`], { type: 'text/tab-separated-values;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = url;
    anchor.download = `polar-genome-db-${program}-all-queries.tsv`;
    anchor.click();
    URL.revokeObjectURL(url);
  };

  return (
    <span className="sequence-export-tooltip" data-tooltip={`Download raw ${program === 'diamond' ? 'DIAMOND' : 'BLAST'} hits for all queries as TSV`}>
      <button type="button" className="sequence-result-export" aria-label={`Download raw ${program === 'diamond' ? 'DIAMOND' : 'BLAST'} hits for all queries as TSV`} onClick={downloadAllResults}>
        <FileDown size={18} aria-hidden="true" />
        <span>TSV</span>
      </button>
    </span>
  );
}

function SequenceQueryResultsExport({ jobId, program, queryId }: { jobId: string; program: SearchProgram; queryId: string }) {
  const downloadQueryResults = (event: React.MouseEvent<HTMLButtonElement>) => {
    event.preventDefault();
    event.stopPropagation();
    if (!jobId) return;
    window.location.assign(`${API_BASE}/api/search/jobs/${encodeURIComponent(jobId)}/raw-results?query_id=${encodeURIComponent(queryId)}`);
  };

  return (
    <span className="sequence-export-tooltip sequence-query-export-tooltip" data-tooltip={`Download raw ${program === 'diamond' ? 'DIAMOND' : 'BLAST'} hits for this query as TSV`}>
      <button type="button" className={`sequence-query-export${program === 'diamond' ? ' sequence-table-raw-export' : ''}`} aria-label={`Download raw hits for ${queryId} as TSV`} onClick={downloadQueryResults}>
        <FileDown size={17} aria-hidden="true" />
      </button>
    </span>
  );
}

function SequenceQueryResults({ speciesId, onOpenGene, jobId, program, queryId, queryLength, queryNumber, hits }: { speciesId: string; onOpenGene: (geneId: string) => void; jobId: string; program: SearchProgram; queryId: string; queryLength: number; queryNumber: number; hits: ApiRecord[] }) {
  const [page, setPage] = useState(0);
  const pageSize = 5;
  const pageCount = Math.max(1, Math.ceil(hits.length / pageSize));
  const activePage = Math.min(page, pageCount - 1);
  const visibleHits = hits.slice(activePage * pageSize, (activePage + 1) * pageSize);

  useEffect(() => setPage(0), [queryId, hits.length]);

  return (
    <details className="candidate-group sequence-query-results" open={queryNumber === 1}>
      <summary className="query-summary">
        <ChevronRight className="query-toggle" size={20} aria-hidden="true" />
        <span className="query-index">Query {queryNumber}</span>
        <span className="query-title">{queryId}</span>
        <div className="query-metadata">
          <span className="query-pill">{hits.length} {hits.length === 1 ? 'match' : 'matches'}</span>
          {queryLength ? <span className="query-pill">{queryLength.toLocaleString()} {program === 'blastn' ? 'bp' : 'aa'}</span> : null}
          <SequenceQueryResultsExport jobId={jobId} program={program} queryId={queryId} />
        </div>
      </summary>
      <section className="sequence-query-content">
        <div className="table-wrap sequence-result-table-wrap">
        <table className="multi-table sequence-result-table">
          <thead>
            <tr>
              <th>Rank</th>
              <th>Gene</th>
              {program === 'blastn' ? <th>Genome coordinate</th> : <><th>Transcript</th><th>Protein</th></>}
              <th>Identity (%)</th>
              <th>Alignment</th>
              <th>E-value</th>
              <th>Bit score</th>
            </tr>
          </thead>
          <tbody>
            {visibleHits.map((hit, index) => {
              const annotation = (hit.annotation ?? {}) as ApiRecord;
              const geneId = annotation.gene_id ? String(annotation.gene_id) : '';
              const coordinate = annotation.seqid && annotation.start && annotation.end
                ? `${annotation.seqid}:${Number(annotation.start).toLocaleString()}-${Number(annotation.end).toLocaleString()}${annotation.strand ? ` (${annotation.strand})` : ''}`
                : '-';
              return (
                <tr key={`${queryId}-${activePage}-${index}`}>
                  <td>{activePage * pageSize + index + 1}</td>
                  <td className={geneId ? 'sequence-result-gene' : ''}>{geneId ? <SequenceGeneLink speciesId={speciesId} geneId={geneId} onOpenGene={onOpenGene} /> : (annotation.target_type === 'genomic_region' ? 'Intergenic region' : '-')}</td>
                  {program === 'blastn' ? <td>{coordinate}</td> : <><td>{annotation.transcript_id ?? '-'}</td><td>{annotation.protein_id ?? (annotation.target_type === 'protein' ? hit.sseqid : '-')}</td></>}
                  <td>{Number(hit.pident ?? 0).toFixed(1)}</td>
                  <td>{Number(hit.length ?? 0).toLocaleString()} {program === 'blastn' ? 'bp' : 'aa'}</td>
                  <td>{formatCell(hit.evalue, 'evalue')}</td>
                  <td>{formatCell(hit.bitscore, 'bitscore')}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
        </div>
        {pageCount > 1 ? (
        <div className="compact-table-pagination">
          <span>{activePage * pageSize + 1}-{Math.min((activePage + 1) * pageSize, hits.length)} of {hits.length}</span>
          <div>
            <button type="button" aria-label={`Previous results for ${queryId}`} disabled={activePage === 0} onClick={() => setPage((value) => Math.max(0, value - 1))}><ChevronLeft size={16} /></button>
            <span>{activePage + 1} / {pageCount}</span>
            <button type="button" aria-label={`Next results for ${queryId}`} disabled={activePage >= pageCount - 1} onClick={() => setPage((value) => Math.min(pageCount - 1, value + 1))}><ChevronRight size={16} /></button>
          </div>
        </div>
        ) : null}
      </section>
    </details>
  );
}

function GeneSearchPanel({
  geneInput,
  onGeneInput,
  speciesId,
  speciesOptions,
  onSelectSpecies,
}: {
  geneInput: string;
  onGeneInput: (value: string) => void;
  speciesId: string;
  speciesOptions: SpeciesSummary[];
  onSelectSpecies: (speciesId: string) => void;
}) {
  const exampleGenes = geneSearchExampleGenesForSpecies(speciesId);

  return (
    <div className="gene-search-panel">
      <CustomSpeciesSelect speciesId={speciesId} speciesOptions={speciesOptions} onSelectSpecies={onSelectSpecies} />
      <label className="gene-query-field">
        <span>Gene IDs</span>
        <textarea
          className="gene-query-input"
          value={geneInput}
          onChange={(event) => onGeneInput(event.target.value)}
          placeholder={`Enter gene IDs for ${speciesOptions.find((item) => item.species_id === speciesId)?.scientific_name ?? speciesId}\n\n${exampleGenes.join('\n')}`}
          rows={3}
        />
        <p className="gene-query-help">Search one or more genes to inspect integrated annotation, expression, and reference evidence.</p>
      </label>
    </div>
  );
}

function SimilarGeneSearchPanel({
  geneInput,
  onGeneInput,
  speciesId,
  taxonGroup,
  speciesOptions,
  onSelectSpecies,
}: {
  geneInput: string;
  onGeneInput: (value: string) => void;
  speciesId: string;
  taxonGroup: string;
  speciesOptions: SpeciesSummary[];
  onSelectSpecies: (speciesId: string) => void;
}) {
  const reference = referenceContextForTaxonGroup(taxonGroup);

  return (
    <div className="gene-search-panel similar-gene-search-panel">
      <div className="similar-species-context" aria-label="Species and reference configuration">
        <CustomSpeciesSelect speciesId={speciesId} speciesOptions={speciesOptions} onSelectSpecies={onSelectSpecies} />
        <p className="similar-reference-inline">
          <span className="similar-taxon-group">{reference.group}</span>
          <span>Functional reference:</span>
          <strong><em>{reference.scientificName}</em> {reference.database}</strong>
        </p>
      </div>
      <label className="gene-query-field">
        <span>Reference gene IDs, transcript IDs, or symbols</span>
        <textarea
          className="gene-query-input"
          value={geneInput}
          onChange={(event) => onGeneInput(event.target.value)}
          placeholder={`Enter one or more ${reference.scientificName} gene IDs, transcript IDs, or symbols\n\n${reference.examples.join('\n')}`}
          rows={6}
        />
        <p className="gene-query-help">Each query is resolved against the selected reference genome before you choose the protein isoforms to compare.</p>
      </label>
    </div>
  );
}

function GoAnalysisInputPanel({
  geneInput,
  backgroundInput,
  minCount,
  correction,
  confidence,
  speciesId,
  taxonGroup,
  speciesOptions,
  onGeneInput,
  onBackgroundInput,
  onMinCount,
  onCorrection,
  onConfidence,
  onSelectSpecies,
  onRun,
  isRunning,
}: {
  geneInput: string;
  backgroundInput: string;
  minCount: number;
  correction: string;
  confidence: Set<string>;
  speciesId: string;
  taxonGroup: string;
  speciesOptions: SpeciesSummary[];
  onGeneInput: (value: string) => void;
  onBackgroundInput: (value: string) => void;
  onMinCount: (value: number) => void;
  onCorrection: (value: string) => void;
  onConfidence: (value: Set<string>) => void;
  onSelectSpecies: (speciesId: string) => void;
  onRun: () => void;
  isRunning: boolean;
}) {
  const reference = referenceContextForTaxonGroup(taxonGroup);
  const exampleGenes = goEnrichmentExampleGenesForSpecies(speciesId).join('\n');
  const [openHint, setOpenHint] = useState<string | null>(null);
  useEffect(() => {
    if (!openHint) return;
    const closeHint = (event: MouseEvent) => {
      const target = event.target;
      if (target instanceof Element && target.closest('.go-option-help')) return;
      setOpenHint(null);
    };
    document.addEventListener('click', closeHint);
    return () => document.removeEventListener('click', closeHint);
  }, [openHint]);
  const readIds = async (event: React.ChangeEvent<HTMLInputElement>, onRead: (value: string) => void) => {
    const file = event.target.files?.[0];
    if (!file) return;
    onRead(await file.text());
    event.target.value = '';
  };
  const toggleConfidence = (value: string) => {
    const next = new Set(confidence);
    if (next.has(value)) next.delete(value);
    else next.add(value);
    onConfidence(next);
  };

  return (
    <div className="go-analysis-workspace">
      <div className="go-species-context" aria-label="Species and reference configuration">
        <CustomSpeciesSelect speciesId={speciesId} speciesOptions={speciesOptions} onSelectSpecies={onSelectSpecies} />
        <p className="similar-reference-inline">
          <span className="similar-taxon-group">{reference.group}</span>
          <span>Functional reference:</span>
          <strong><em>{reference.scientificName}</em> {reference.database}</strong>
        </p>
      </div>

      <section className="go-analysis-fields" aria-label="GO enrichment input">
        <div className="go-input-field">
          <label className="go-input-label">
            <span>Gene list</span>
            <textarea
              value={geneInput}
              spellCheck={false}
              placeholder={exampleGenes}
              onChange={(event) => onGeneInput(event.target.value)}
            />
          </label>
          <div className="go-input-foot">
            <small>One gene ID per line, or IDs separated by commas, tabs, or spaces.</small>
            <div className="go-upload-meta">
              <label className="go-upload-button">
                <Upload size={15} aria-hidden="true" />
                <span>Upload IDs</span>
                <input type="file" accept=".txt,.csv,.tsv,text/plain,text/csv" onChange={(event) => void readIds(event, onGeneInput)} />
              </label>
              <small className="go-input-limit">Up to 1,000 IDs</small>
            </div>
          </div>
        </div>
        <div className="go-input-field">
          <label className="go-input-label">
            <span>Background list <em>optional</em></span>
            <textarea
              value={backgroundInput}
              spellCheck={false}
              placeholder="Leave blank to use all genes with inferred GO annotation."
              onChange={(event) => onBackgroundInput(event.target.value)}
            />
          </label>
          <div className="go-input-foot">
            <small>Use a custom background when the experiment tested a defined subset of genes.</small>
            <div className="go-upload-meta">
              <label className="go-upload-button">
                <Upload size={15} aria-hidden="true" />
                <span>Upload background</span>
                <input type="file" accept=".txt,.csv,.tsv,text/plain,text/csv" onChange={(event) => void readIds(event, onBackgroundInput)} />
              </label>
              <small className="go-input-limit">Up to 100,000 IDs</small>
            </div>
          </div>
        </div>
      </section>

      <section className="go-analysis-options" aria-label="GO analysis options">
        <label>
          <span className="go-option-heading">Minimum term count <GoOptionHelp id="minimum-count" text="Only GO terms supported by at least this many submitted genes are tested." openHint={openHint} onToggle={setOpenHint} /></span>
          <input type="number" min={1} max={1000} value={minCount} onChange={(event) => onMinCount(Number(event.target.value) || 1)} />
        </label>
        <label>
          <span className="go-option-heading">Multiple testing correction <GoOptionHelp id="multiple-testing" text="Adjusts enrichment p-values for the many GO terms tested in one analysis." openHint={openHint} onToggle={setOpenHint} /></span>
          <CompactSelect value={correction} onChange={onCorrection} ariaLabel="Multiple testing correction" options={[
            { value: 'benjamini-hochberg', label: 'Benjamini-Hochberg' },
            { value: 'bonferroni', label: 'Bonferroni' },
            { value: 'holm', label: 'Holm' },
          ]} />
        </label>
        <fieldset>
          <legend className="go-option-heading">Homolog confidence <GoOptionHelp id="homolog-confidence" text="Select which similarity confidence tiers may contribute reference-derived GO annotations." openHint={openHint} onToggle={setOpenHint} /></legend>
          <div className="go-confidence-options">
            {['high', 'medium', 'low'].map((value) => (
              <label key={value}>
                <input type="checkbox" checked={confidence.has(value)} onChange={() => toggleConfidence(value)} />
                <span>{value.charAt(0).toUpperCase() + value.slice(1)}</span>
              </label>
            ))}
          </div>
        </fieldset>
      </section>
      <div className="tool-actions go-input-actions">
        <button type="button" className="gene-search-submit" onClick={onRun} disabled={isRunning}>{isRunning ? 'Running...' : 'Run GO analysis'}</button>
        <button type="button" className="secondary gene-example-button" onClick={() => onGeneInput(exampleGenes)}>Example</button>
      </div>
    </div>
  );
}

function GoOptionHelp({
  id,
  text,
  openHint,
  onToggle,
}: {
  id: string;
  text: string;
  openHint: string | null;
  onToggle: (id: string | null) => void;
}) {
  const isOpen = openHint === id;
  return (
    <span className="go-option-help">
      <button type="button" aria-label={`About ${id.replace('-', ' ')}`} aria-expanded={isOpen} onClick={() => onToggle(isOpen ? null : id)}>
        <CircleHelp size={15} aria-hidden="true" />
      </button>
      {isOpen ? <span className="go-option-tooltip" role="tooltip">{text}</span> : null}
    </span>
  );
}

type CompactSelectOption = { value: string; label: string };

function CompactSelect({
  value,
  options,
  onChange,
  disabled = false,
  ariaLabel,
  variant = 'compact',
}: {
  value: string;
  options: CompactSelectOption[];
  onChange: (value: string) => void;
  disabled?: boolean;
  ariaLabel: string;
  variant?: 'compact' | 'primer';
}) {
  const [isOpen, setIsOpen] = useState(false);
  const selectRef = useRef<HTMLDivElement>(null);
  const selected = options.find((option) => option.value === value) ?? options[0];

  useEffect(() => {
    const closeMenu = (event: MouseEvent) => {
      if (!selectRef.current?.contains(event.target as Node)) setIsOpen(false);
    };
    document.addEventListener('mousedown', closeMenu);
    return () => document.removeEventListener('mousedown', closeMenu);
  }, []);

  useEffect(() => {
    if (disabled) setIsOpen(false);
  }, [disabled]);

  return (
    <div ref={selectRef} className={`${variant === 'primer' ? 'primer-transcript-field report-select' : 'compact-select'}${isOpen ? ' open' : ''}${disabled ? ' disabled' : ''}`} onKeyDown={(event) => { if (event.key === 'Escape') { setIsOpen(false); selectRef.current?.querySelector('button')?.focus(); } }}>
      <button type="button" className={variant === 'primer' ? 'primer-transcript-trigger' : 'compact-select-trigger'} aria-haspopup="listbox" aria-expanded={isOpen} aria-label={ariaLabel} disabled={disabled} onClick={() => setIsOpen((open) => !open)}>
        <span>{selected?.label ?? ''}</span>
        <ChevronDown size={16} aria-hidden="true" />
      </button>
      {isOpen ? (
        <div className={variant === 'primer' ? 'primer-transcript-menu' : 'compact-select-menu'} role="listbox" aria-label={ariaLabel}>
          {options.map((option) => (
            <button type="button" key={option.value} className={option.value === value ? 'active' : ''} role="option" aria-selected={option.value === value} onClick={() => { onChange(option.value); setIsOpen(false); }}>
              {option.label}
            </button>
          ))}
        </div>
      ) : null}
    </div>
  );
}

function GoEnrichmentResults({ result }: { result: ApiRecord }) {
  const categories = Array.isArray(result.categories) ? result.categories as ApiRecord[] : [];
  const invalid = Array.isArray(result.invalid) ? result.invalid : [];
  const [activeTab, setActiveTab] = useState<'standard' | 'semantic'>('standard');
  const [enrichmentMode, setEnrichmentMode] = useState<'all' | 'direct'>('all');
  const resultSets = result.result_sets as ApiRecord | undefined;
  const activeCategories = Array.isArray(resultSets?.[enrichmentMode]) ? resultSets?.[enrichmentMode] as ApiRecord[] : categories;
  return (
    <section className="go-enrichment-results multi-panel" aria-label="GO enrichment results">
      <header className="multi-head go-enrichment-head">
        <div>
          <h2>GO Functional Annotation Chart</h2>
          <p>Reference homolog-inferred annotations · one-sided hypergeometric test</p>
        </div>
      </header>
      {invalid.length ? <p className="go-result-note">{invalid.length.toLocaleString()} submitted ID{invalid.length === 1 ? '' : 's'} could not be resolved and were excluded.</p> : null}
      {Number(result.query_outside_background_count ?? 0) ? <p className="go-result-note">{Number(result.query_outside_background_count).toLocaleString()} query gene{Number(result.query_outside_background_count) === 1 ? '' : 's'} were outside the custom background and were excluded.</p> : null}
      <nav className="go-result-tabs" aria-label="GO result views">
        <button type="button" aria-pressed={activeTab === 'standard'} className={activeTab === 'standard' ? 'active' : ''} onClick={() => setActiveTab('standard')}><BarChart3 size={20} aria-hidden="true" /><span>Standard Enrichment</span></button>
        <button type="button" aria-pressed={activeTab === 'semantic'} className={activeTab === 'semantic' ? 'active' : ''} onClick={() => setActiveTab('semantic')}><Network size={20} aria-hidden="true" /><span>Semantic GO</span></button>
      </nav>
      {activeTab === 'standard' ? (
        <div className="go-enrichment-categories">
          {activeCategories.map((category) => <GoEnrichmentCategory key={`${enrichmentMode}-${String(category.namespace)}`} category={category} minCountFloor={Number(result.min_count ?? 1)} enrichmentMode={enrichmentMode} onEnrichmentMode={setEnrichmentMode} />)}
        </div>
      ) : (
        <GoSemanticPanel categories={Array.isArray(resultSets?.all) ? resultSets?.all as ApiRecord[] : categories} minCountFloor={Number(result.min_count ?? 1)} />
      )}
    </section>
  );
}

function GoEnrichmentCategory({
  category,
  minCountFloor,
  enrichmentMode,
  onEnrichmentMode,
}: {
  category: ApiRecord;
  minCountFloor: number;
  enrichmentMode: 'direct' | 'all';
  onEnrichmentMode: (mode: 'direct' | 'all') => void;
}) {
  const terms = Array.isArray(category.terms) ? category.terms as ApiRecord[] : [];
  const [filter, setFilter] = useState({ level: 'all', minCount: minCountFloor, maxPvalue: 0.05, maxFdr: 1 });
  const [page, setPage] = useState(0);
  const defaultSelected = useMemo(() => new Set<string>(terms
    .filter((term) => Number(term.count ?? 0) >= minCountFloor
      && Number(term.pvalue ?? 1) <= 0.05
      && Number(term.adjusted_pvalue ?? 1) <= 1)
    .slice(0, 10)
    .map((term) => String(term.go_id))), [category.terms, minCountFloor]);
  const [selected, setSelected] = useState<Set<string>>(() => new Set(defaultSelected));
  const [bubbleOpen, setBubbleOpen] = useState(false);
  const [bubbleOptions, setBubbleOptions] = useState({ layout: 'auto', xMetric: 'gene_ratio', colorMetric: 'adjusted_pvalue', palette: 'viridis' });
  const pageSize = 10;
  const filteredTerms = terms.filter((term) => (
    goLevelAllowed(term.level, filter.level)
    &&
    Number(term.count ?? 0) >= filter.minCount
    && Number(term.pvalue ?? 1) <= filter.maxPvalue
    && Number(term.adjusted_pvalue ?? 1) <= filter.maxFdr
  ));
  const totalPages = Math.max(1, Math.ceil(filteredTerms.length / pageSize));
  const currentPage = Math.min(page, totalPages - 1);
  const visibleTerms = filteredTerms.slice(currentPage * pageSize, (currentPage + 1) * pageSize);

  useEffect(() => {
    setFilter({ level: 'all', minCount: minCountFloor, maxPvalue: 0.05, maxFdr: 1 });
    setPage(0);
    setSelected(new Set(defaultSelected));
    setBubbleOpen(false);
    setBubbleOptions({ layout: 'auto', xMetric: 'gene_ratio', colorMetric: 'adjusted_pvalue', palette: 'viridis' });
  }, [defaultSelected, minCountFloor]);

  const updateFilter = (patch: Partial<typeof filter>) => {
    setFilter((current) => ({
      ...current,
      ...patch,
      minCount: Math.max(minCountFloor, Number(patch.minCount ?? current.minCount) || minCountFloor),
    }));
    setPage(0);
  };
  const resetFilter = () => {
    setFilter({ level: 'all', minCount: minCountFloor, maxPvalue: 0.05, maxFdr: 1 });
    setPage(0);
  };
  const toggleSelected = (goId: string) => {
    setSelected((current) => {
      const next = new Set(current);
      if (next.has(goId)) next.delete(goId);
      else next.add(goId);
      return next;
    });
  };
  const selectedTerms = filteredTerms.filter((term) => selected.has(String(term.go_id))).slice(0, 15);
  const exportSelectedTerms = terms.filter((term) => selected.has(String(term.go_id)));
  const exportTerms = (scope: 'all' | 'selected') => downloadTsv(
    goTermsTsv(scope === 'all' ? terms : exportSelectedTerms, String(category.namespace ?? category.short_label), enrichmentMode),
    `${String(category.short_label).toLowerCase()}_${enrichmentMode}_enrichment_${scope}.tsv`,
  );

  return (
    <details className="candidate-group go-enrichment-category" open={category.short_label === 'BP'}>
      <summary className="query-summary">
        <ChevronRight size={19} className="query-toggle" aria-hidden="true" />
        <span className="query-index">{category.short_label}</span>
        <span className="query-title">{category.label}</span>
        <span className="query-metadata">
          <span className="query-pill">{filteredTerms.length.toLocaleString()} visible terms</span>
          <span className="query-pill">{Number(category.query_annotated ?? 0).toLocaleString()} annotated query genes</span>
        </span>
      </summary>
      <div className="go-enrichment-table-section">
        <p className="go-category-note">{Number(category.term_count ?? 0).toLocaleString()} tested terms · minimum count cannot be lower than the submitted analysis count ({minCountFloor.toLocaleString()}).</p>
        <div className="go-filter-panel" aria-label={`${category.label} display filters`}>
          <label>Annotation mode
            <CompactSelect value={enrichmentMode} onChange={(value) => onEnrichmentMode(value as 'direct' | 'all')} ariaLabel="Annotation mode" options={[
              { value: 'all', label: 'All ancestors' },
              { value: 'direct', label: 'Direct' },
            ]} />
          </label>
          <label>Level
            <CompactSelect value={filter.level} disabled={enrichmentMode === 'direct'} onChange={(value) => updateFilter({ level: value })} ariaLabel="GO hierarchy level" options={[
              { value: 'all', label: 'All levels' },
              { value: '2', label: 'Level 2' },
              { value: '3', label: 'Level 3' },
              { value: '4', label: 'Level 4' },
              { value: '5', label: 'Level 5' },
              { value: '6', label: 'Level 6' },
              { value: '7plus', label: 'Level 7+' },
            ]} />
          </label>
          <label>Count &gt;=
            <input type="number" min={minCountFloor} value={filter.minCount} onChange={(event) => updateFilter({ minCount: Number(event.target.value) || minCountFloor })} />
          </label>
          <label>p-value &lt;=
            <input type="number" min={0} max={1} step={0.001} value={filter.maxPvalue} onChange={(event) => updateFilter({ maxPvalue: Number(event.target.value) || 0 })} />
          </label>
          <label>Adjusted p-value &lt;=
            <input type="number" min={0} max={1} step={0.001} value={filter.maxFdr} onChange={(event) => updateFilter({ maxFdr: Number(event.target.value) || 0 })} />
          </label>
          <div className="go-filter-actions">
            <span className="go-filter-tooltip" data-tooltip="Reset filters">
              <button type="button" className="go-filter-icon-button" onClick={resetFilter} aria-label="Reset filters">
                <RotateCcw size={17} aria-hidden="true" />
              </button>
            </span>
            <span className="go-filter-tooltip" data-tooltip="Clear selected terms">
              <button type="button" className="go-filter-icon-button" disabled={!selected.size} onClick={() => setSelected(new Set())} aria-label="Clear selected terms">
                <ListX size={18} aria-hidden="true" />
              </button>
            </span>
            <span className="go-filter-tooltip" data-tooltip={'All terms TSV\nCurrent annotation mode'}>
              <button type="button" className="go-filter-icon-button" disabled={!terms.length} onClick={() => exportTerms('all')} aria-label="Export all terms">
                <Download size={17} aria-hidden="true" />
              </button>
            </span>
            <span className="go-filter-tooltip" data-tooltip="Selected terms TSV">
              <button type="button" className="go-filter-icon-button" disabled={!exportSelectedTerms.length} onClick={() => exportTerms('selected')} aria-label="Export selected terms">
                <FileDown size={17} aria-hidden="true" />
              </button>
            </span>
            <span className="go-filter-tooltip" data-tooltip={bubbleOpen ? 'Hide bubble plot' : 'Show bubble plot'}>
              <button type="button" className="go-filter-icon-button go-bubble-toggle" disabled={!selected.size} onClick={() => setBubbleOpen((value) => !value)} aria-label={bubbleOpen ? 'Hide bubble plot' : 'Show bubble plot'}>
                <BarChart3 size={18} aria-hidden="true" />
              </button>
            </span>
          </div>
        </div>
        {bubbleOpen ? <GoBubblePlot section={String(category.short_label)} rows={selectedTerms} selectedCount={selected.size} options={bubbleOptions} onOptions={setBubbleOptions} /> : null}
        {filteredTerms.length ? (
          <>
            <div className="table-wrap">
              <table className="compact-table go-enrichment-table">
                <thead>
                  <tr>
                    <th>Select</th>
                    <th>GO term</th>
                    <th>Term</th>
                    <th>Count</th>
                    <th>Percent</th>
                    <th>p-value</th>
                    <th>Adjusted<br />p-value</th>
                    <th>Fold<br />enrichment</th>
                    <th>Population</th>
                  </tr>
                </thead>
                <tbody>
                  {visibleTerms.map((term) => (
                    <tr key={String(term.go_id)}>
                      <td><input className="go-term-select" type="checkbox" checked={selected.has(String(term.go_id))} onChange={() => toggleSelected(String(term.go_id))} /></td>
                      <td className="go-term-id">{term.go_id}</td>
                      <td className="go-term-name">{term.term}</td>
                      <td>{Number(term.count ?? 0).toLocaleString()}</td>
                      <td>{formatGoPercent(Number(term.percent ?? 0))}</td>
                      <td>{formatGoPValue(Number(term.pvalue ?? 1))}</td>
                      <td>{formatGoPValue(Number(term.adjusted_pvalue ?? 1))}</td>
                      <td>{Number(term.fold_enrichment ?? 0).toFixed(2)}</td>
                      <td>{Number(term.pop_hits ?? 0).toLocaleString()} / {Number(term.pop_total ?? 0).toLocaleString()}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            {totalPages > 1 ? (
              <div className="compact-table-pagination">
                <span>{currentPage * pageSize + 1}-{Math.min((currentPage + 1) * pageSize, filteredTerms.length)} of {filteredTerms.length}</span>
                <div>
                  <button type="button" aria-label={`Previous ${category.label} terms`} disabled={currentPage === 0} onClick={() => setPage((value) => Math.max(0, value - 1))}><ChevronLeft size={16} /></button>
                  <span>{currentPage + 1} / {totalPages}</span>
                  <button type="button" aria-label={`Next ${category.label} terms`} disabled={currentPage >= totalPages - 1} onClick={() => setPage((value) => Math.min(totalPages - 1, value + 1))}><ChevronRight size={16} /></button>
                </div>
              </div>
            ) : null}
          </>
        ) : <p className="go-empty-category">No {category.label} terms matched the current display filters.</p>}
      </div>
    </details>
  );
}

type SemanticOptions = {
  palette: 'plasma' | 'viridis' | 'cividis' | 'magma' | 'Set2';
  si: number;
  maxTerms: number;
  labelMode: 'numbered' | 'description' | 'description-numbered';
  bubbleSize: 'members' | 'frequency';
};

const defaultSemanticOptions = (): SemanticOptions => ({
  palette: 'plasma',
  si: 0.5,
  maxTerms: 15,
  labelMode: 'numbered',
  bubbleSize: 'members',
});

function GoSemanticPanel({ categories, minCountFloor }: { categories: ApiRecord[]; minCountFloor: number }) {
  return (
    <div className="go-semantic-categories">
      <p className="go-semantic-intro">Semantic GO uses All ancestors terms only. Filtered enrichment terms are reduced with GO-Figure! using semantic similarity.</p>
      {categories.map((category) => <GoSemanticCategory key={String(category.namespace)} category={category} minCountFloor={minCountFloor} />)}
    </div>
  );
}

function GoSemanticCategory({ category, minCountFloor }: { category: ApiRecord; minCountFloor: number }) {
  const plotRef = useRef<HTMLDivElement>(null);
  const terms = Array.isArray(category.terms) ? category.terms as ApiRecord[] : [];
  const filter = { minCount: minCountFloor, maxPvalue: 0.05, maxFdr: 1 };
  const [options, setOptions] = useState<SemanticOptions>(defaultSemanticOptions);
  const [semanticHelp, setSemanticHelp] = useState<string | null>(null);
  const [semanticResult, setSemanticResult] = useState<ApiRecord | null>(null);
  const [semanticError, setSemanticError] = useState('');
  const [isRunning, setIsRunning] = useState(false);
  const semanticOptionsRef = useRef<HTMLDivElement>(null);
  const filteredTerms = terms
    .filter((term) => (
      Number(term.count ?? 0) >= filter.minCount
      && Number(term.pvalue ?? 1) <= filter.maxPvalue
      && Number(term.adjusted_pvalue ?? 1) <= filter.maxFdr
    ))
    .sort((left, right) => Number(left.pvalue ?? 1) - Number(right.pvalue ?? 1));
  const selectedTerms = filteredTerms.slice(0, options.maxTerms);
  useEffect(() => {
    if (!semanticHelp) return undefined;
    const closeHelp = (event: MouseEvent) => {
      if (!semanticOptionsRef.current?.contains(event.target as Node)) setSemanticHelp(null);
    };
    document.addEventListener('mousedown', closeHelp);
    return () => document.removeEventListener('mousedown', closeHelp);
  }, [semanticHelp]);
  const runSemantic = async () => {
    if (selectedTerms.length < 2) {
      setSemanticError('At least two filtered All ancestors GO terms are required for semantic plotting.');
      return;
    }
    setIsRunning(true);
    setSemanticError('');
    try {
      const result = await apiPost<ApiRecord>('/api/go/semantic', {
        section: category.short_label,
        palette: options.palette,
        si: options.si,
        max_terms: options.maxTerms,
        label_mode: options.labelMode,
        bubble_size: options.bubbleSize,
        font_size: 'small',
        legend_size: 'small',
        terms: selectedTerms.map((term) => ({ go_id: term.go_id, pvalue: term.pvalue })),
      });
      setSemanticResult(result);
    } catch (error) {
      setSemanticError(error instanceof Error ? error.message : 'Semantic GO failed.');
    } finally {
      setIsRunning(false);
    }
  };
  const tsvRows = parseSemanticTsv(String(semanticResult?.semantic_tsv ?? ''));

  return (
    <details className="candidate-group go-enrichment-category go-semantic-category" open={category.short_label === 'BP'}>
      <summary className="query-summary">
        <ChevronRight size={19} className="query-toggle" aria-hidden="true" />
        <span className="query-index">{category.short_label}</span>
        <span className="query-title">{category.label}</span>
        <span className="query-metadata">
          <span className="query-pill">{filteredTerms.length.toLocaleString()} eligible terms</span>
          <span className="query-pill">{selectedTerms.length.toLocaleString()} GO-Figure terms</span>
        </span>
      </summary>
      <div className="go-enrichment-table-section">
        <p className="go-category-note">All ancestors only · GO-Figure receives the top {options.maxTerms.toLocaleString()} eligible terms ranked by raw p-value.</p>
        <div ref={semanticOptionsRef} className="semantic-options-panel" aria-label={`${category.label} semantic GO options`}>
          <label><span className="semantic-option-heading">Palette <GoOptionHelp id="palette" text="Color palette used by GO-Figure for semantic clusters." openHint={semanticHelp} onToggle={setSemanticHelp} /></span>
            <CompactSelect value={options.palette} onChange={(value) => setOptions((current) => ({ ...current, palette: value as SemanticOptions['palette'] }))} ariaLabel="Semantic GO palette" options={[
              { value: 'plasma', label: 'Plasma' }, { value: 'viridis', label: 'Viridis' }, { value: 'cividis', label: 'Cividis' }, { value: 'magma', label: 'Magma' }, { value: 'Set2', label: 'Set2' },
            ]} />
          </label>
          <label><span className="semantic-option-heading">SI <GoOptionHelp id="si" text="Similarity cutoff. Lower values merge related GO terms; higher values keep groups more separate." openHint={semanticHelp} onToggle={setSemanticHelp} /></span>
            <input type="number" min={0.1} max={0.9} step={0.05} value={options.si} onChange={(event) => setOptions((current) => ({ ...current, si: Math.max(0.1, Math.min(0.9, Number(event.target.value) || 0.5)) }))} />
          </label>
          <label><span className="semantic-option-heading">Terms <GoOptionHelp id="terms" text="Maximum enriched GO terms sent to GO-Figure after ranking by p-value." openHint={semanticHelp} onToggle={setSemanticHelp} /></span>
            <input type="number" min={5} max={30} step={1} value={options.maxTerms} onChange={(event) => setOptions((current) => ({ ...current, maxTerms: Math.max(5, Math.min(30, Math.floor(Number(event.target.value) || 15))) }))} />
          </label>
          <label><span className="semantic-option-heading">Plot label <GoOptionHelp id="plot-label" text="Label style shown inside the semantic GO plot. Numbered labels are mapped in the summary table." openHint={semanticHelp} onToggle={setSemanticHelp} /></span>
            <CompactSelect value={options.labelMode} onChange={(value) => setOptions((current) => ({ ...current, labelMode: value as SemanticOptions['labelMode'] }))} ariaLabel="Semantic GO plot label" options={[
              { value: 'numbered', label: 'Numbered' }, { value: 'description', label: 'Description' }, { value: 'description-numbered', label: 'Description + number' },
            ]} />
          </label>
          <label><span className="semantic-option-heading">Size <GoOptionHelp id="size" text="Selects whether bubble size represents cluster members or GO term frequency." openHint={semanticHelp} onToggle={setSemanticHelp} /></span>
            <CompactSelect value={options.bubbleSize} onChange={(value) => setOptions((current) => ({ ...current, bubbleSize: value as SemanticOptions['bubbleSize'] }))} ariaLabel="Semantic GO bubble size" options={[
              { value: 'members', label: 'Members' }, { value: 'frequency', label: 'Frequency' },
            ]} />
          </label>
          <button type="button" className="semantic-run-button" disabled={selectedTerms.length < 2 || isRunning} onClick={() => void runSemantic()}>
            {isRunning ? <><LoaderCircle size={17} className="semantic-spinner" aria-hidden="true" />Running</> : 'Run semantic GO'}
          </button>
        </div>
        {semanticError ? <p className="go-semantic-error">{semanticError}</p> : null}
        {semanticResult ? (
          <div className="gofigure-panel">
            <div className="gofigure-head">
              <span>{Number(semanticResult.terms_used ?? selectedTerms.length).toLocaleString()} GO-Figure terms{semanticResult.truncated ? ' · top terms used' : ''}</span>
              <div className="go-bubble-export-actions" aria-label="Semantic plot export controls">
                <span className="go-filter-tooltip" data-tooltip="Download summary TSV"><button type="button" className="go-filter-icon-button" aria-label="Download semantic summary TSV" onClick={() => downloadTsv(String(semanticResult.semantic_tsv ?? ''), `${category.short_label.toLowerCase()}_semantic_go_summary.tsv`)}><FileText size={17} aria-hidden="true" /></button></span>
                <span className="go-filter-tooltip" data-tooltip="Download full TSV"><button type="button" className="go-filter-icon-button" aria-label="Download semantic full TSV" onClick={() => downloadTsv(String(semanticResult.semantic_full_tsv ?? ''), `${category.short_label.toLowerCase()}_semantic_go_full.tsv`)}><FileDown size={17} aria-hidden="true" /></button></span>
                <GoPlotExportButtons plotRef={plotRef} filename={`${String(category.short_label).toLowerCase()}_semantic_go`} label="semantic plot" />
              </div>
            </div>
            {tsvRows.length ? <SemanticSummaryTable rows={tsvRows} /> : null}
            <div ref={plotRef} className="gofigure-svg-wrap" dangerouslySetInnerHTML={{ __html: String(semanticResult.svg ?? '') }} />
          </div>
        ) : null}
      </div>
    </details>
  );
}

function parseSemanticTsv(tsv: string) {
  const lines = tsv.trim().split(/\r?\n/).filter(Boolean);
  if (lines.length < 2) return [] as Array<Record<string, string>>;
  const headers = lines[0].split('\t');
  return lines.slice(1).map((line) => headers.reduce((record, header, index) => {
    record[header] = line.split('\t')[index] ?? '';
    return record;
  }, {} as Record<string, string>));
}

function downloadTsv(contents: string, filename: string) {
  if (!contents) return;
  const href = URL.createObjectURL(new Blob([contents], { type: 'text/tab-separated-values;charset=utf-8' }));
  const anchor = document.createElement('a');
  anchor.href = href;
  anchor.download = filename;
  anchor.click();
  window.setTimeout(() => URL.revokeObjectURL(href), 1000);
}

function GoPlotExportButtons({ plotRef, filename, label }: { plotRef: React.RefObject<HTMLDivElement | null>; filename: string; label: string }) {
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState('');
  async function download(format: 'png' | 'svg') {
    setBusy(true);
    setError('');
    try {
      const svg = plotRef.current?.querySelector('svg');
      if (!svg) throw new Error('No plot is available to export.');
      await exportSvgPlot(svg, format, filename, 'Arial, Helvetica, sans-serif');
    } catch (error) {
      setError(error instanceof Error ? error.message : 'Plot export failed.');
    } finally {
      setBusy(false);
    }
  }
  return <>
    {(['png', 'svg'] as const).map((format) => <span key={format} className="go-filter-tooltip" data-tooltip={`${format.toUpperCase()} export`}><button type="button" className="go-bubble-export-button" disabled={busy} aria-label={`Export ${label} as ${format.toUpperCase()}`} onClick={() => void download(format)}>{format.toUpperCase()}</button></span>)}
    {error ? <span className="go-semantic-error" role="alert">{error}</span> : null}
  </>;
}

function SemanticSummaryTable({ rows }: { rows: Array<Record<string, string>> }) {
  return (
    <div className="table-wrap go-semantic-summary-wrap">
      <table className="compact-table go-semantic-summary-table">
        <thead><tr><th>Plot No.</th><th>Representative</th><th>Description</th><th>Members</th><th>p-value</th><th>Frequency</th></tr></thead>
        <tbody>{rows.map((row, index) => <tr key={`${row.representative ?? row.GO ?? index}`}><td>{row.dotCount ?? String(index + 1)}</td><td>{row.representative ?? row.GO ?? ''}</td><td>{row.description ?? row.legend ?? ''}</td><td>{row.members ?? ''}</td><td>{row.pval ?? row.PValue ?? ''}</td><td>{row.frequency ?? ''}</td></tr>)}</tbody>
      </table>
    </div>
  );
}

function goLevelAllowed(level: unknown, selector: string) {
  if (selector === 'all') return true;
  const value = Number(level);
  if (!Number.isFinite(value)) return false;
  if (selector === '7plus') return value >= 7;
  return value === Number(selector);
}

function GoBubblePlot({
  section,
  rows,
  selectedCount,
  options,
  onOptions,
}: {
  section: string;
  rows: ApiRecord[];
  selectedCount: number;
  options: { layout: string; xMetric: string; colorMetric: string; palette: string };
  onOptions: (value: { layout: string; xMetric: string; colorMetric: string; palette: string }) => void;
}) {
  const plotRef = useRef<HTMLDivElement>(null);
  if (!rows.length) {
    return <p className="go-empty-category">Select one or more visible GO terms to draw the bubble plot.</p>;
  }
  const metric = (term: ApiRecord) => options.xMetric === 'fold_enrichment'
    ? Number(term.fold_enrichment ?? 0)
    : Number(term.percent ?? 0) / 100;
  const colorMetric = (term: ApiRecord) => Number(options.colorMetric === 'pvalue' ? term.pvalue : term.adjusted_pvalue) || 1;
  const ordered = [...rows].sort((a, b) => metric(b) - metric(a));
  const layout = options.layout === 'auto' ? (ordered.length > 10 ? 'square' : 'rectangle') : options.layout;
  const left = 318;
  const right = 162;
  const top = 24;
  const labelLines = ordered.map((term) => wrapGoLabel(String(term.term ?? term.go_id), left - 28));
  const lineHeight = 19;
  let contentHeight = 0;
  const rowCenters = labelLines.map((lines) => {
    const rowHeight = Math.max(40, lines.length * lineHeight + 16);
    const center = contentHeight + rowHeight / 2;
    contentHeight += rowHeight;
    return top + center;
  });
  const plotWidth = layout === 'square' ? 350 : 330;
  const plotHeight = Math.max(layout === 'square' ? 280 : 250, contentHeight);
  const width = left + plotWidth + right;
  const height = top + plotHeight + 58;
  const values = ordered.map(metric);
  const maxValue = Math.max(...values, options.xMetric === 'fold_enrichment' ? 1 : 0.05) * 1.12;
  const counts = ordered.map((term) => Number(term.count ?? 0));
  const minCount = Math.min(...counts);
  const maxCount = Math.max(...counts);
  const countLegend = bubbleCountLegend(maxCount);
  const scores = ordered.map((term) => -Math.log10(Math.max(colorMetric(term), 1e-300)));
  const minScore = Math.min(...scores);
  const maxScore = Math.max(...scores);
  const xTicks = nicePlotTicks(maxValue);
  const xAxisMax = xTicks[xTicks.length - 1] ?? maxValue;
  const radius = (count: number) => {
    if (maxCount === minCount) return 8;
    const normalized = Math.max(0, Math.min(1, (count - minCount) / (maxCount - minCount)));
    return 6 + Math.sqrt(normalized) * 7;
  };
  const fill = (score: number) => {
    const t = maxScore === minScore ? 1 : (score - minScore) / (maxScore - minScore);
    return goBubblePalette(options.palette, t);
  };
  const xLabel = options.xMetric === 'fold_enrichment' ? 'Fold enrichment' : 'Gene ratio (count / query)';
  const colorLabel = options.colorMetric === 'pvalue' ? 'p-value' : 'Adjusted p-value';
  const gradientId = `go-bubble-gradient-${section}`;
  const countLegendStartY = top + 42;
  const colorLegendLabelY = countLegendStartY + countLegend.length * 31 + 22;
  const colorLegendBarY = colorLegendLabelY + 12;
  return (
    <section className="go-bubble-panel" aria-label={`${section} GO bubble plot`}>
      <div className="go-bubble-heading">
        <div><span>Selected Bubble Plot</span><small>{Math.min(selectedCount, 15)} of {selectedCount} selected terms</small></div>
      </div>
      <div className="go-bubble-controls">
        <label>Layout
          <CompactSelect value={options.layout} onChange={(value) => onOptions({ ...options, layout: value })} ariaLabel="Bubble plot layout" options={[
            { value: 'auto', label: 'Auto' }, { value: 'rectangle', label: 'Rectangle' }, { value: 'square', label: 'Square' },
          ]} />
        </label>
        <label>X-axis
          <CompactSelect value={options.xMetric} onChange={(value) => onOptions({ ...options, xMetric: value })} ariaLabel="Bubble plot x axis" options={[
            { value: 'gene_ratio', label: 'Gene ratio' }, { value: 'fold_enrichment', label: 'Fold enrichment' },
          ]} />
        </label>
        <label>Color
          <CompactSelect value={options.colorMetric} onChange={(value) => onOptions({ ...options, colorMetric: value })} ariaLabel="Bubble plot color" options={[
            { value: 'adjusted_pvalue', label: 'Adjusted p-value' }, { value: 'pvalue', label: 'p-value' },
          ]} />
        </label>
        <label>Palette
          <CompactSelect value={options.palette} onChange={(value) => onOptions({ ...options, palette: value })} ariaLabel="Bubble plot palette" options={[
            { value: 'viridis', label: 'Viridis' }, { value: 'cividis', label: 'Cividis' }, { value: 'plasma', label: 'Plasma' }, { value: 'blue_red', label: 'Blue-red' },
          ]} />
        </label>
        <div className="go-bubble-export-actions" aria-label="Bubble plot export controls">
          <GoPlotExportButtons plotRef={plotRef} filename={`${section.toLowerCase()}_standard_go_bubble`} label="bubble plot" />
        </div>
      </div>
      <div ref={plotRef} className="go-bubble-svg-wrap">
        <svg className="go-bubble-svg" viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${section} GO bubble plot, ${xLabel} by ${colorLabel}`}>
          <defs><linearGradient id={gradientId} x1="0%" x2="0%" y1="100%" y2="0%"><stop offset="0%" stopColor={goBubblePalette(options.palette, 0)} /><stop offset="100%" stopColor={goBubblePalette(options.palette, 1)} /></linearGradient></defs>
          <rect x={left} y={top} width={plotWidth} height={plotHeight} className="go-chart-frame" />
          {xTicks.map((tick, index) => {
            const x = left + (tick / xAxisMax) * plotWidth;
            return <g key={index}><line x1={x} x2={x} y1={top} y2={top + plotHeight} className="go-chart-grid" /><text x={x} y={top + plotHeight + 22} textAnchor="middle" className="go-chart-axis">{formatGoPlotTick(tick)}</text></g>;
          })}
          {ordered.map((term, index) => {
            const value = values[index];
            const x = left + (value / xAxisMax) * plotWidth;
            const y = rowCenters[index];
            const label = String(term.term ?? term.go_id);
            return <g key={String(term.go_id)}><line x1={left} x2={left + plotWidth} y1={y} y2={y} className="go-chart-row-grid" /><text x={left - 10} textAnchor="end" className="go-chart-term">{labelLines[index].map((line, lineIndex) => <tspan key={lineIndex} x={left - 10} y={y + 4 + (lineIndex - (labelLines[index].length - 1) / 2) * lineHeight}>{line}</tspan>)}</text><circle cx={x} cy={y} r={radius(Number(term.count ?? 0))} fill={fill(scores[index])} className="go-chart-bubble"><title>{`${term.go_id} · ${label}\n${xLabel}: ${value.toFixed(options.xMetric === 'fold_enrichment' ? 2 : 3)}\nCount: ${term.count}\n${colorLabel}: ${formatGoPValue(colorMetric(term))}`}</title></circle></g>;
          })}
          <text x={left + plotWidth / 2} y={height - 8} textAnchor="middle" className="go-chart-label">{xLabel}</text>
          <text x={left + plotWidth + 28} y={top + 18} className="go-chart-count-heading">Count</text>
          {countLegend.map((count, index) => <g key={count}><circle cx={left + plotWidth + 44} cy={countLegendStartY + index * 31} r={radius(count)} className="go-chart-size" /><text x={left + plotWidth + 66} y={countLegendStartY + 5 + index * 31} className="go-chart-count-label">{count}</text></g>)}
          <text x={left + plotWidth + 28} y={colorLegendLabelY} className="go-chart-label">{colorLabel}</text>
          <rect x={left + plotWidth + 40} y={colorLegendBarY} width="15" height="86" fill={`url(#${gradientId})`} className="go-chart-gradient" />
          <text x={left + plotWidth + 64} y={colorLegendBarY + 10} className="go-chart-axis">{formatGoPValue(10 ** -maxScore)}</text>
          <text x={left + plotWidth + 64} y={colorLegendBarY + 86} className="go-chart-axis">{formatGoPValue(10 ** -minScore)}</text>
        </svg>
      </div>
    </section>
  );
}

function nicePlotTicks(maxValue: number, desiredTicks = 4) {
  if (!Number.isFinite(maxValue) || maxValue <= 0) return [0, 1];
  const step = niceStep(maxValue / desiredTicks);
  const upper = step * Math.ceil(maxValue / step);
  const ticks: number[] = [];
  for (let tick = 0; tick <= upper + step * 0.001; tick += step) ticks.push(Math.round(tick * 1e10) / 1e10);
  return ticks;
}

function bubbleCountLegend(maximum: number) {
  if (maximum <= 1) return [Math.max(0, maximum)];
  const step = Math.max(1, niceStep(maximum / 4));
  const values: number[] = [];
  for (let value = step; value < maximum && values.length < 3; value += step) values.push(value);
  return [...new Set([...values, maximum])].slice(-4);
}

function formatGoPlotTick(value: number) {
  if (Number.isInteger(value)) return value.toLocaleString();
  if (Math.abs(value) < 0.1) return value.toFixed(2);
  return value.toFixed(1).replace(/\.0$/, '');
}

function goBubblePalette(name: string, value: number) {
  const palettes: Record<string, [number, number, number][]> = {
    viridis: [[68, 1, 84], [59, 82, 139], [33, 145, 140], [253, 231, 37]],
    cividis: [[0, 34, 78], [80, 103, 111], [169, 151, 92], [253, 234, 69]],
    plasma: [[13, 8, 135], [126, 3, 168], [203, 71, 119], [240, 249, 33]],
    blue_red: [[37, 99, 235], [220, 38, 38]],
  };
  const stops = palettes[name] ?? palettes.viridis;
  const scaled = Math.max(0, Math.min(1, value)) * (stops.length - 1);
  const index = Math.min(stops.length - 2, Math.floor(scaled));
  const fraction = scaled - index;
  const color = stops[index].map((channel, channelIndex) => Math.round(channel + (stops[index + 1][channelIndex] - channel) * fraction));
  return `rgb(${color.join(', ')})`;
}

function formatGoPValue(value: number) {
  if (!Number.isFinite(value)) return '-';
  if (value === 0) return '< 1e-300';
  if (value < 0.001) return value.toExponential(2);
  return value.toFixed(3);
}

function formatGoPercent(value: number) {
  if (!Number.isFinite(value)) return '-';
  return `${value.toFixed(value < 10 ? 1 : 0)}%`;
}

function ReferenceIsoformSelectionPanel({
  groups,
  selectedProteinIds,
  onToggleProtein,
  onSelectAll,
  onClearSelection,
  onRun,
  isRunning,
}: {
  groups: ApiRecord[];
  selectedProteinIds: string[];
  onToggleProtein: (proteinId: string) => void;
  onSelectAll: (proteinIds: string[]) => void;
  onClearSelection: () => void;
  onRun: () => void;
  isRunning: boolean;
}) {
  const matchedGroups = groups.filter((group) => Array.isArray(group.matches) && group.matches.length);
  const unresolvedGroups = groups.filter((group) => !Array.isArray(group.matches) || !group.matches.length);
  const proteinIds = [...new Set(matchedGroups.flatMap((group) => (
    (group.matches as ApiRecord[]).map((match) => String(match.reference_protein_id ?? '')).filter(Boolean)
  )))];
  return (
    <section className="similar-lookup-results">
      <article className="multi-panel similar-isoform-panel">
        <header className="multi-head">
          <div>
            <h2>Reference Protein Isoforms</h2>
            <p>{selectedProteinIds.length.toLocaleString()} selected</p>
          </div>
          <div className="isoform-selection-actions" aria-label="Isoform selection controls">
            <button type="button" className="table-icon-button" onClick={() => onSelectAll(proteinIds)} aria-label="Select all resolved isoforms" title="Select all isoforms">
              <CheckCheck size={19} aria-hidden="true" />
            </button>
            <button type="button" className="table-icon-button" onClick={onClearSelection} aria-label="Clear isoform selection" title="Clear selection">
              <ListX size={19} aria-hidden="true" />
            </button>
          </div>
        </header>
        {matchedGroups.length ? (
          <div className="table-wrap">
            <table className="multi-table similar-isoform-table">
              <thead>
                <tr>
                  <th>Select</th>
                  <th>Input query</th>
                  <th>Reference transcript</th>
                  <th>Gene</th>
                  <th>Gene symbol</th>
                  <th>Length</th>
                </tr>
              </thead>
              <tbody>
                {matchedGroups.flatMap((group) => (group.matches as ApiRecord[]).map((match) => {
                  const proteinId = String(match.reference_protein_id ?? '');
                  return (
                    <tr key={`${group.query}-${proteinId}`} className={selectedProteinIds.includes(proteinId) ? 'active' : ''}>
                      <td><input type="checkbox" checked={selectedProteinIds.includes(proteinId)} onChange={() => onToggleProtein(proteinId)} aria-label={`Select ${match.reference_transcript_id ?? proteinId}`} /></td>
                      <td>{String(group.query)}</td>
                      <td>{match.reference_transcript_id ?? proteinId}</td>
                      <td>{match.reference_gene_id ?? '-'}</td>
                      <td>{match.reference_gene_symbol ?? '-'}</td>
                      <td>{Number(match.length ?? 0).toLocaleString()} aa</td>
                    </tr>
                  );
                }))}
              </tbody>
            </table>
          </div>
        ) : <EmptyResult text="No reference protein isoforms were resolved." />}
        {matchedGroups.length ? (
          <div className="similar-run-row">
            <div>
              <button type="button" className="similar-run-button" onClick={onRun} disabled={!selectedProteinIds.length || isRunning}>
                {isRunning ? <LoaderCircle size={16} className="semantic-spinner" aria-hidden="true" /> : <Play size={16} fill="currentColor" aria-hidden="true" />}
                <span>{isRunning ? 'Running...' : 'Run selected isoforms'}</span>
              </button>
            </div>
          </div>
        ) : null}
        {unresolvedGroups.length ? <div className="multi-missing">{unresolvedGroups.map((group) => <span className="error-inline" key={String(group.query)}>{group.query}: unresolved</span>)}</div> : null}
      </article>
    </section>
  );
}

function SimilarQueryExport({ speciesId, sourceJobId, proteinId }: { speciesId: string; sourceJobId: string; proteinId: string }) {
  const [jobId, setJobId] = useState('');
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState('');
  useEffect(() => {
    if (!jobId) return;
    let cancelled = false;
    let timer: number | undefined;
    const poll = async () => {
      try {
        const result = await apiGet<ApiRecord>(`/api/species/${speciesId}/similar/exports/${jobId}`);
        if (cancelled) return;
        if (result.status === 'failed') throw new Error(String(result.error ?? 'Export failed.'));
        if (result.status === 'completed') {
          window.location.assign(`${API_BASE}/api/species/${speciesId}/similar/exports/${jobId}/download`);
          setBusy(false);
          setJobId('');
          return;
        }
        timer = window.setTimeout(() => void poll(), 1000);
      } catch (error) {
        if (!cancelled) {
          setError(error instanceof Error ? error.message : 'Export failed.');
          setBusy(false);
          setJobId('');
        }
      }
    };
    void poll();
    return () => { cancelled = true; window.clearTimeout(timer); };
  }, [jobId, speciesId]);
  async function startExport() {
    setBusy(true);
    setError('');
    try {
      const result = await apiPost<ApiRecord>(`/api/species/${speciesId}/similar/jobs/${sourceJobId}/export`, { reference_protein_id: proteinId });
      setJobId(String(result.job_id));
    } catch (error) {
      setError(error instanceof Error ? error.message : 'Export failed.');
      setBusy(false);
    }
  }
  return <span className="query-export-control" onClick={(event) => { event.stopPropagation(); }}>
    <span className="go-filter-tooltip" data-tooltip={busy ? 'Preparing TSV...' : 'TSV export\nCandidates + pairwise'}>
      <button type="button" className="query-export-button" aria-label={`Export ${proteinId} candidate and pairwise TSV`} disabled={busy} onClick={(event) => { event.preventDefault(); void startExport(); }}>
        {busy ? <LoaderCircle size={19} className="semantic-spinner" aria-hidden="true" /> : <FileDown size={19} aria-hidden="true" />}
      </button>
    </span>
    {error ? <span role="alert" className="query-export-error">{error}</span> : null}
  </span>;
}

function SimilarCandidateResultsPanel({
  speciesId,
  similarJobId,
  results,
  pairwiseJobs,
  activePairwiseKey,
  onRunPairwise,
  onOpenGene,
}: {
  speciesId: string;
  similarJobId: string;
  results: ApiRecord;
  pairwiseJobs: Record<string, ApiRecord>;
  activePairwiseKey: string;
  onRunPairwise: (referenceProteinId: string, targetProteinId: string) => void;
  onOpenGene: (geneId: string) => void;
}) {
  const groups = Array.isArray(results.groups) ? results.groups as ApiRecord[] : [];
  const candidateCount = Number(results.candidate_count ?? 0);
  return (
    <section className="similar-lookup-results">
      <article className="multi-panel similar-candidate-panel">
        <header className="multi-head">
          <div>
            <h2>Candidate Rank Results</h2>
            <p>{candidateCount.toLocaleString()} target gene candidates</p>
          </div>
        </header>
        {groups.length ? groups.map((group, groupIndex) => {
          const query = group.query_protein as ApiRecord | undefined;
          const candidates = (Array.isArray(group.candidates) ? group.candidates as ApiRecord[] : []).slice(0, 5);
          const queryName = String(query?.gene_symbol ?? query?.gene_id ?? query?.protein_id ?? '-');
          const referenceProteinId = String(query?.protein_id ?? '');
          return (
            <details className="candidate-group" key={String(query?.protein_id ?? groupIndex)} open={groupIndex === 0}>
              <summary className="query-summary">
                <ChevronRight className="query-toggle" size={20} aria-hidden="true" />
                <span className="query-index">Query {groupIndex + 1}</span>
                <span className="query-title">{queryName}</span>
                <div className="query-metadata" onClick={(event) => event.stopPropagation()}>
                  <span className="query-pill">Gene {query?.gene_id ?? '-'}</span>
                  <span className="query-pill">Isoform {query?.transcript_id ?? query?.protein_id ?? '-'}</span>
                  <span className="query-pill">{Number(query?.length ?? 0).toLocaleString()} aa</span>
                  <span className="query-pill">{candidates.length} candidates</span>
                </div>
                <SimilarQueryExport speciesId={speciesId} sourceJobId={similarJobId} proteinId={referenceProteinId} />
              </summary>
              <section className="candidate-rank-section">
                <h4>Candidate Rank Table</h4>
                {candidates.length ? (
                  <div className="table-wrap">
                    <table className="multi-table similar-candidate-table">
                      <thead>
                        <tr><th>Rank</th><th aria-label="Basket"></th><th>Target gene</th><th>Transcript</th><th>Identity</th><th>Coverage</th><th>E-value</th><th>Bitscore</th><th>Pairwise</th></tr>
                      </thead>
                      <tbody>
                        {candidates.map((candidate) => (
                          <tr key={`${query?.protein_id}-${candidate.gene_id}-${candidate.protein_id}`}>
                            <td>{candidate.rank}</td>
                            <td>
                              {candidate.gene_id && <BasketButton speciesId={speciesId} geneId={String(candidate.gene_id)} className="secondary table-icon-button basket-gene-button" />}
                            </td>
                            <td className="similar-candidate-gene"><a className="candidate-gene-link" href={`/tools?tool=gene&species=${encodeURIComponent(speciesId)}&gene=${encodeURIComponent(String(candidate.gene_id))}`} onClick={(event) => { if (event.button !== 0 || event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return; event.preventDefault(); onOpenGene(String(candidate.gene_id)); }}>{candidate.gene_id ?? '-'}</a></td>
                            <td>{candidate.transcript_id ?? '-'}</td>
                            <td>{Number(candidate.pident ?? 0).toFixed(1)}%</td>
                            <td>{candidate.alignment_length ?? 0}/{candidate.query_length ?? 0} aa</td>
                            <td>{Number(candidate.evalue ?? 0).toExponential(1)}</td>
                            <td>{Number(candidate.bitscore ?? 0).toFixed(1)}</td>
                            <td>
                              <PairwiseDomainButton
                                job={pairwiseJobs[`${referenceProteinId}:${candidate.protein_id}`]}
                                onRun={() => onRunPairwise(referenceProteinId, String(candidate.protein_id ?? ''))}
                              />
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                ) : <EmptyResult text="No target gene candidates passed the DIAMOND search threshold." />}
              </section>
              {candidates.map((candidate) => {
                const pairwiseKey = `${referenceProteinId}:${candidate.protein_id}`;
                const job = pairwiseJobs[pairwiseKey];
                return activePairwiseKey === pairwiseKey && job?.status === 'completed' && job.results ? <PairwiseDomainPanel key={pairwiseKey} result={job.results} candidate={candidate} /> : null;
              })}
              <CandidateExpressionHeatmap candidates={candidates} expression={results.expression as ApiRecord | undefined} />
            </details>
          );
        }) : <EmptyResult text="No candidate result groups were returned." />}
      </article>
      {groups.length ? (
        <article className="multi-panel similar-tree-panel">
          <SimilarInputPhylogenyPanel speciesId={speciesId} similarJobId={similarJobId} groups={groups} />
        </article>
      ) : null}
    </section>
  );
}

function SimilarInputPhylogenyPanel({ speciesId, similarJobId, groups }: { speciesId: string; similarJobId: string; groups: ApiRecord[] }) {
  const treeRef = useRef<HTMLDivElement>(null);
  const [exporting, setExporting] = useState(false);
  const [exportError, setExportError] = useState('');
  async function downloadTree(format: 'png' | 'svg') {
    const svg = treeRef.current?.querySelector<SVGSVGElement>('svg.reference-tree');
    if (!svg) return;
    setExporting(true);
    setExportError('');
    try { await exportTree(svg, format, `${speciesId}_similar_${similarJobId}_phylogeny`); }
    catch (error) { setExportError(error instanceof Error ? error.message : 'Tree export failed.'); }
    finally { setExporting(false); }
  }
  const inputs = useMemo(() => groups.map((group) => group.query_protein as ApiRecord | undefined).filter((item): item is ApiRecord => Boolean(item?.protein_id)), [groups]);
  const inputIds = useMemo(() => inputs.map((item) => String(item.protein_id)), [inputs]);
  const [selectedIds, setSelectedIds] = useState<string[]>(() => inputIds.slice(0, 1));
  const [model, setModel] = useState('LG+G');
  const [job, setJob] = useState<ApiRecord | null>(null);
  const [error, setError] = useState('');

  useEffect(() => {
    setSelectedIds(inputIds.slice(0, 1));
    setJob(null);
    setError('');
  }, [similarJobId]);

  const running = job?.status === 'queued' || job?.status === 'running';
  const result = job?.status === 'completed' && job.results ? job.results as ApiRecord : null;
  const toggle = (proteinId: string) => {
    setSelectedIds((current) => current.includes(proteinId) ? current.filter((item) => item !== proteinId) : [...current, proteinId]);
  };
  const buildTree = async () => {
    if (!selectedIds.length) {
      setError('Select at least one input isoform.');
      return;
    }
    if (selectedIds.length > 8) {
      setError('A maximum of 8 input isoforms can be included in one tree.');
      return;
    }
    setError('');
    setJob({ status: 'queued' });
    try {
      const queued = await apiPost<ApiRecord>(`/api/species/${speciesId}/similar/tree`, {
        similar_job_id: similarJobId,
        reference_protein_ids: selectedIds,
        model,
      });
      setJob(queued);
      for (let attempt = 0; attempt < 120; attempt += 1) {
        await new Promise((resolve) => window.setTimeout(resolve, attempt === 0 ? 250 : 700));
        const current = await apiGet<ApiRecord>(`/api/phylogeny/jobs/${queued.job_id}`);
        setJob(current);
        if (current.status === 'completed' || current.status === 'failed') {
          if (current.status === 'failed') setError(String(current.error ?? 'Phylogenetic tree calculation failed.'));
          return;
        }
      }
      setError('Phylogenetic tree calculation timed out.');
    } catch (requestError) {
      setError(requestError instanceof Error ? requestError.message : String(requestError));
    }
  };

  return (
    <section className="similar-tree-section">
      <header className="multi-head">
        <div>
          <h2>Phylogenetic Tree</h2>
          <p>Selected reference isoforms and their top target candidates</p>
        </div>
      </header>
      <div className="similar-tree-controls">
        <label className="phylogeny-model-select similar-tree-model-select">
          <span>Protein model</span>
          <CompactSelect value={model} disabled={running} onChange={setModel} ariaLabel="Protein substitution model" options={[
            { value: 'LG+G', label: 'LG+G' },
            { value: 'WAG+G', label: 'WAG+G' },
            { value: 'JTT+G', label: 'JTT+G' },
            { value: 'Auto', label: 'Auto (LG / WAG / JTT)' },
          ]} />
        </label>
        <div className="similar-tree-inputs" role="group" aria-label="Input isoforms included in the tree">
          {inputs.map((input) => {
            const proteinId = String(input.protein_id);
            const label = String(input.gene_symbol ?? input.gene_id ?? proteinId);
            const transcript = String(input.transcript_id ?? proteinId);
            return <label className="similar-tree-isoform" key={proteinId}><input type="checkbox" checked={selectedIds.includes(proteinId)} onChange={() => toggle(proteinId)} /><span>{label} · {transcript}</span></label>;
          })}
        </div>
        <button type="button" className="similar-tree-run" disabled={running || !selectedIds.length} onClick={() => void buildTree()}>
          {running ? <LoaderCircle size={16} className="semantic-spinner" aria-hidden="true" /> : <GitBranch size={16} aria-hidden="true" />}
          <span>{running ? 'Building tree...' : 'Build phylogenetic tree'}</span>
        </button>
        <div className="tree-export-actions">
          {(['png', 'svg'] as const).map(format => <span key={format} className="go-filter-tooltip" data-tooltip={`${format.toUpperCase()} export`}>
            <button type="button" className="go-bubble-export-button" aria-label={`Export similar tree as ${format.toUpperCase()}`} disabled={!result || exporting} onClick={() => void downloadTree(format)}>{format.toUpperCase()}</button>
          </span>)}
        </div>
      </div>
      {error ? <div className="tool-message error">{error}</div> : null}
      {!error && running ? <p className="phylogeny-status" role="status"><LoaderCircle size={17} className="semantic-spinner" aria-hidden="true" />Aligning proteins and inferring tree...</p> : null}
      {exportError && <p role="alert" className="tool-message error">{exportError}</p>}
      <div ref={treeRef}>{result ? <ProteinPhylogenyTree gene={{ gene_id: 'Selected similar-gene inputs' }} result={result} /> : null}</div>
    </section>
  );
}

function PairwiseDomainButton({ job, onRun }: { job?: ApiRecord; onRun: () => void }) {
  const running = job?.status === 'queued' || job?.status === 'running';
  return (
    <button type="button" className="secondary pairwise-button" onClick={onRun} disabled={running} title="Compare shared Pfam domains">
      <GitCompareArrows size={16} aria-hidden="true" />
      <span>{running ? 'Checking...' : 'Pairwise'}</span>
    </button>
  );
}

function CandidateExpressionHeatmap({ candidates, expression }: { candidates: ApiRecord[]; expression?: ApiRecord }) {
  const samples = Array.isArray(expression?.samples) ? expression!.samples as ApiRecord[] : [];
  const values = expression?.values as Record<string, Record<string, number>> | undefined;
  if (!expression?.available || !samples.length || !values || !candidates.length) return null;
  const relevantCandidates = candidates.slice(0, 5);
  const contextForSample = (sample: ApiRecord) => {
    const metadata = (sample.metadata ?? {}) as ApiRecord;
    const temperature = String(metadata.temperature ?? '');
    const infection = String(metadata.infection ?? '');
    const tissue = String(metadata.tissue ?? '');
    const condition = String(sample.condition ?? metadata.condition ?? '');
    const replicate = String(sample.replicate ?? String(sample.sample_id ?? '').match(/_(\d+)$/)?.[1] ?? '');
    return {
      primary: temperature ? `${temperature} C` : tissue || condition || String(sample.display_name ?? sample.sample_id),
      secondary: infection ? (infection === '+' ? 'Infected' : 'Control') : temperature ? condition : tissue ? condition : '',
      replicate,
      temperature: temperature ? Number(temperature) : Number.NaN,
      infectionOrder: infection === '-' ? 0 : infection === '+' ? 1 : 2,
    };
  };
  const orderedSamples = [...samples].sort((left, right) => {
    const a = contextForSample(left);
    const b = contextForSample(right);
    return (Number.isFinite(a.temperature) ? a.temperature : Number.MAX_SAFE_INTEGER) - (Number.isFinite(b.temperature) ? b.temperature : Number.MAX_SAFE_INTEGER)
      || a.infectionOrder - b.infectionOrder
      || a.primary.localeCompare(b.primary)
      || Number(a.replicate || 0) - Number(b.replicate || 0);
  });
  const transformed = relevantCandidates.flatMap((candidate) => samples.map((sample) => Math.log2(Number(values[String(candidate.gene_id)]?.[String(sample.sample_id)] ?? 0) + 1)));
  const maximum = Math.max(1, ...transformed);
  return (
    <section className="candidate-expression-section">
      <header>
        <div>
          <h4>Candidate Expression Heatmap</h4>
          <p>Top target candidates · log2(CPM + 1)</p>
        </div>
      </header>
      <div className="candidate-expression-scroll">
        <div className="candidate-expression-grid" style={{ gridTemplateColumns: `minmax(178px, 1.3fr) repeat(${orderedSamples.length}, minmax(72px, 1fr))` }}>
          <span className="candidate-expression-corner">Target gene</span>
          {orderedSamples.map((sample) => {
            const context = contextForSample(sample);
            return <span className="candidate-expression-sample" key={String(sample.sample_id)}><b>{context.primary}</b>{context.secondary ? <b className="sample-context">{context.secondary}</b> : null}{context.replicate ? <small>Rep {context.replicate}</small> : null}</span>;
          })}
          {relevantCandidates.map((candidate) => (
            <React.Fragment key={String(candidate.gene_id)}>
              <span className="candidate-expression-gene">{String(candidate.gene_id)}</span>
              {orderedSamples.map((sample) => {
                const rawValue = Number(values[String(candidate.gene_id)]?.[String(sample.sample_id)] ?? 0);
                const logValue = Math.log2(rawValue + 1);
                const intensity = Math.max(0, Math.min(1, logValue / maximum));
                const background = `rgb(${Math.round(246 - 209 * intensity)}, ${Math.round(248 - 135 * intensity)}, ${Math.round(247 - 145 * intensity)})`;
                return <span className="candidate-expression-cell" key={String(sample.sample_id)} style={{ backgroundColor: background, color: intensity > 0.57 ? '#ffffff' : '#24302d' }} title={`${candidate.gene_id} · ${sample.display_name ?? sample.sample_id}: ${rawValue.toFixed(2)} CPM`}>{logValue.toFixed(1)}</span>;
              })}
            </React.Fragment>
          ))}
        </div>
      </div>
      <div className="candidate-expression-legend"><span>Low</span><i /><span>High</span></div>
    </section>
  );
}

function PairwiseDomainPanel({ result, candidate }: { result: ApiRecord; candidate: ApiRecord }) {
  const domains = Array.isArray(result.shared_domains) ? result.shared_domains as ApiRecord[] : [];
  return (
    <section className="pair-domain-section">
      <h4>Pairwise Domain Comparison: {candidate.gene_id}</h4>
      {!result.available ? <EmptyResult text={String(result.reason ?? 'Pairwise Pfam comparison is unavailable.')} /> : domains.length ? (
        <>
          <PairwiseDomainSchematic result={result} domains={domains} />
          <div className="table-wrap">
            <table className="multi-table pairwise-domain-table">
              <thead><tr><th>Pfam domain</th><th>Accession</th><th>Query range</th><th>Target range</th><th>Identity</th><th>Aligned positions</th></tr></thead>
              <tbody>{domains.map((domain) => <tr key={`${domain.pfam_id}-${domain.query_range}-${domain.target_range}`}><td>{domain.domain_name ?? domain.pfam_id}</td><td>{domain.pfam_id ?? '-'}</td><td>{domain.query_range}</td><td>{domain.target_range}</td><td>{domain.identity == null ? '-' : `${Number(domain.identity).toFixed(1)}%`}</td><td>{domain.aligned_positions ?? '-'}</td></tr>)}</tbody>
            </table>
          </div>
        </>
      ) : <EmptyResult text="No shared Pfam domains were found for this candidate." />}
    </section>
  );
}

function PairwiseDomainSchematic({ result, domains }: { result: ApiRecord; domains: ApiRecord[] }) {
  return (
    <div className="pairwise-domain-architectures">
      {domains.map((domain, index) => <PairwiseDomainArchitecture key={`${domain.pfam_id}-${domain.query_range}-${domain.target_range}`} result={result} domain={domain} index={index} />)}
    </div>
  );
}

function PairwiseDomainArchitecture({ result, domain, index }: { result: ApiRecord; domain: ApiRecord; index: number }) {
  const width = 920;
  const height = 86;
  const trackX = 245;
  const trackW = 590;
  const queryLength = Math.max(1, Number(result.query_length ?? 1));
  const targetLength = Math.max(1, Number(result.target_length ?? 1));
  const colors = ['#126b5b', '#b14b36', '#5b6f95', '#8a6d1d', '#6d5a7b', '#44736a', '#9a5a78', '#3f7d9b'];
  const color = colors[index % colors.length];
  const parseRange = (value: unknown) => String(value ?? '').split('-').map((item) => Number(item));
  const position = (aa: number, length: number) => trackX + ((Math.max(1, aa) - 1) / length) * trackW;
  const [queryStart, queryEnd] = parseRange(domain.query_range);
  const [targetStart, targetEnd] = parseRange(domain.target_range);
  const queryX = position(queryStart, queryLength);
  const targetX = position(targetStart, targetLength);
  const queryWidth = Math.max(4, position(queryEnd, queryLength) - queryX);
  const targetWidth = Math.max(4, position(targetEnd, targetLength) - targetX);
  const domainName = String(domain.domain_name ?? domain.pfam_id);
  return (
    <svg className="pairwise-domain-map" viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`Shared Pfam domain ${index + 1}: ${domainName}`}>
      <title>{`${domain.pfam_id} · ${domainName}`}</title>
      <text x={trackX} y="16" className="pairwise-domain-map-title">{index + 1}. {domain.pfam_id} · {domainName}</text>
      <text x="18" y="47" className="pairwise-domain-protein-label">Reference | {result.query_transcript_id ?? result.query_protein_id}</text>
      <text x="18" y="73" className="pairwise-domain-protein-label">Target | {result.target_transcript_id ?? result.target_protein_id}</text>
      <line x1={trackX} x2={trackX + trackW} y1="46" y2="46" className="domain-backbone" />
      <line x1={trackX} x2={trackX + trackW} y1="72" y2="72" className="domain-backbone" />
      <text x={trackX} y="34" className="domain-map-axis">1</text>
      <text x={trackX + trackW} y="34" className="domain-map-axis" textAnchor="end">{queryLength.toLocaleString()} aa</text>
      <text x={trackX} y="59" className="domain-map-axis">1</text>
      <text x={trackX + trackW} y="59" className="domain-map-axis" textAnchor="end">{targetLength.toLocaleString()} aa</text>
      <rect x={queryX} y="42" width={queryWidth} height="8" rx="2" fill={color} />
      <rect x={targetX} y="68" width={targetWidth} height="8" rx="2" fill={color} />
    </svg>
  );
}

function GeneSearchResultsPanel({
  speciesId,
  results,
  selectedGeneId,
  onOpenReport,
}: {
  speciesId: string;
  results: GeneSearchResult[];
  selectedGeneId?: string;
  onOpenReport: (geneId: string) => void;
}) {
  const foundResults = results.filter((item) => item.found && item.gene);
  const missingResults = results.filter((item) => !item.found || !item.gene);

  return (
    <article className="multi-panel search-results-panel">
      <header className="multi-head">
        <div>
          <h2>Search Results</h2>
          <p>{foundResults.length.toLocaleString()} found{missingResults.length ? `, ${missingResults.length.toLocaleString()} unresolved` : ''}</p>
        </div>
      </header>
      {foundResults.length ? (
        <div className="table-wrap">
          <table className="multi-table search-gene-table">
            <thead>
              <tr>
                <th></th>
                <th></th>
                <th>Gene</th>
                <th>Location</th>
                <th>Strand</th>
                <th>Length</th>
                <th>Transcript</th>
                <th>Expression</th>
                <th>Top Pfam</th>
                <th>Top reference</th>
              </tr>
            </thead>
            <tbody>
              {foundResults.map((item) => {
                const gene = item.gene!;
                const topHit = Array.isArray(item.similar?.hits) ? item.similar.hits[0] : null;
                return (
                  <tr
                    key={`${item.query}-${gene.gene_id}`}
                    className={selectedGeneId === gene.gene_id ? 'active' : ''}
                  >
                    <td>
                      <BasketButton speciesId={speciesId} geneId={gene.gene_id} className="secondary table-icon-button basket-gene-button" />
                    </td>
                    <td>
                      <button
                        type="button"
                        className="secondary table-icon-button view-report"
                        aria-label={`Open report for ${gene.gene_id}`}
                        title="Report"
                        onClick={() => onOpenReport(gene.gene_id)}
                      >
                        <FileText size={18} strokeWidth={2.1} />
                      </button>
                    </td>
                    <td>
                      <button type="button" className="gene-result-link" onClick={() => onOpenReport(gene.gene_id)}>
                        {gene.gene_id}
                      </button>
                    </td>
                    <td>{formatSearchLocation(gene)}</td>
                    <td>{gene.strand ?? '-'}</td>
                    <td>{(Number(gene.end ?? 0) - Number(gene.start ?? 0) + 1).toLocaleString()} bp</td>
                    <td>{Number(item.transcripts?.length ?? 0).toLocaleString()}</td>
                    <td>{expressionPeakForSearch(item.expression)}</td>
                    <td className="topDomain">{topPfamForSearch(item)}</td>
                    <td className="topReference">{topHit?.reference_gene_id ?? '-'}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      ) : (
        <EmptyResult text="No matched genes were found." />
      )}
    </article>
  );
}

function GeneResultTable({
  expressionData,
  gene,
  proteins,
  similarData,
  transcripts,
}: {
  expressionData: ApiRecord | null;
  gene: ApiRecord;
  proteins: ApiRecord[];
  similarData: ApiRecord | null;
  transcripts: ApiRecord[];
}) {
  const similarHits: ApiRecord[] = Array.isArray(similarData?.hits) ? similarData.hits : [];
  const cpmRows: ApiRecord[] = Array.isArray(expressionData?.expression?.cpm) ? expressionData.expression.cpm : [];
  const topHit = similarHits[0];
  const meanCpm = cpmRows.length ? cpmRows.reduce((sum, item) => sum + Number(item.value ?? 0), 0) / cpmRows.length : null;
  return (
    <CompactTable
      title="Search results"
      columns={['gene_id', 'location', 'transcripts', 'protein_id', 'mean_cpm', 'top_reference_gene', 'identity']}
      rows={[{
        gene_id: gene.gene_id,
        location: formatLocation(gene),
        transcripts: transcripts.length,
        protein_id: proteins[0]?.protein_id,
        mean_cpm: meanCpm,
        top_reference_gene: topHit?.reference_gene_id,
        identity: topHit?.pident ? `${formatCell(topHit.pident)}%` : null,
      }]}
    />
  );
}

function GeneReport({
  domainData,
  expressionData,
  features,
  gene,
  proteins,
  sequenceData,
  similarData,
  speciesId,
  transcriptModels,
  transcripts,
}: {
  domainData: ApiRecord | null;
  expressionData: ApiRecord | null;
  features: ApiRecord[];
  gene: ApiRecord;
  proteins: ApiRecord[];
  sequenceData: ApiRecord | null;
  similarData: ApiRecord | null;
  speciesId: string;
  transcriptModels: ApiRecord[];
  transcripts: ApiRecord[];
}) {
  const similarHits = Array.isArray(similarData?.hits) ? similarData.hits : [];
  const maxFeatureCount = (featureType: string) => Math.max(
    0,
    ...transcripts.map((transcript) => features.filter(
      (feature) => feature.transcript_id === transcript.transcript_id && feature.feature_type === featureType,
    ).length),
  );
  const cdsCount = maxFeatureCount('CDS');
  const exonCount = maxFeatureCount('exon');
  const topHit = similarHits[0];

  return (
    <article className="gene-card smarlens-report" data-gene={gene.gene_id}>
      <header className="gene-head">
        <div className="gene-title">
          <h2>{gene.gene_id}</h2>
        </div>
      </header>
      <div className="sections">
      <section className="block">
        <h3>Summary</h3>
        <div className="kv">
          <span>Gene</span><span className="summary-gene-actions expression-mode-group">
            <strong>{gene.gene_id}</strong>
            <BasketButton speciesId={speciesId} geneId={gene.gene_id} label />
            <button type="button" onClick={() => window.open(`/tools?tool=genome&species=${encodeURIComponent(speciesId)}&gene=${encodeURIComponent(gene.gene_id)}`, '_blank', 'noopener,noreferrer')}><Layers3 size={15} aria-hidden="true" />Genome Browser</button>
            <button type="button" title="Download gene-oriented genomic sequence" onClick={() => window.location.assign(`${API_BASE}/api/species/${encodeURIComponent(speciesId)}/genes/${encodeURIComponent(gene.gene_id)}/export?format=fasta`)}><Dna size={15} aria-hidden="true" />Sequence FASTA</button>
            <button type="button" title="Download all transcript annotations for this gene" onClick={() => window.location.assign(`${API_BASE}/api/species/${encodeURIComponent(speciesId)}/genes/${encodeURIComponent(gene.gene_id)}/export?format=gtf`)}><FileDown size={15} aria-hidden="true" />Transcript GTF</button>
          </span>
          <span>Locus</span><span>{formatLocation(gene)}</span>
          <span>Gene length</span><span>{(Number(gene.end ?? 0) - Number(gene.start ?? 0) + 1).toLocaleString()} bp</span>
          <span>Transcript model</span><span>{transcripts.length.toLocaleString()} transcripts · {exonCount.toLocaleString()} exons · {cdsCount.toLocaleString()} CDS</span>
          <span>Expression</span><span>{expressionSummaryForReport(expressionData)}</span>
          <span>Functional hints</span><span>{functionalHintsForReport(topHit, domainData)}</span>
        </div>
      </section>

      <section className="block">
        <h3>Expression</h3>
        <ExpressionProfile expressionData={expressionData} speciesId={speciesId} />
      </section>

      <section className="block">
        <h3>Transcript Models</h3>
        <TranscriptModelFigure gene={gene} transcriptModels={transcriptModels} transcripts={transcripts} features={features} />
        <div className="legend">
          <span><i className="swatch cds-swatch" />CDS</span>
          <span><i className="swatch exon-swatch" />Exon</span>
          <span><i className="swatch intron-swatch" />Intron</span>
          <span className="legend-note">UTR features are not available in the current annotation.</span>
        </div>
      </section>

      <section className="block">
        <h3>Sequence</h3>
        <SequencePanel sequenceData={sequenceData} gene={gene} speciesId={speciesId} transcripts={transcripts} />
      </section>

      <section className="block domain-block">
        <div className="block-head-row">
          <h3>Protein Domains</h3>
        </div>
        <DomainEvidencePanel domainData={domainData} geneId={gene.gene_id} proteins={proteins} speciesId={speciesId} />
      </section>

      <section className="block similarity-block">
        <div className="block-head-row">
          <h3>Reference Similarity</h3>
        </div>
        <ProteinSimilarityPanel geneId={gene.gene_id} proteins={proteins} speciesId={speciesId} />
      </section>

      <section className="block phylogeny-block">
        <div className="block-head-row">
          <h3>Phylogenetic Tree</h3>
        </div>
        <ReferencePhylogenyPanel gene={gene} proteins={proteins} speciesId={speciesId} />
      </section>
      </div>
    </article>
  );
}

function TranscriptModelFigure({
  features,
  gene,
  transcriptModels,
  transcripts,
}: {
  features: ApiRecord[];
  gene: ApiRecord;
  transcriptModels: ApiRecord[];
  transcripts: ApiRecord[];
}) {
  const [showAllModels, setShowAllModels] = useState(false);
  const [modelMode, setModelMode] = useState<'genomic' | 'spliced'>('genomic');
  useEffect(() => {
    setShowAllModels(false);
    setModelMode('genomic');
  }, [gene.gene_id]);
  const models: ApiRecord[] = transcriptModels.length
    ? transcriptModels
    : transcripts.map((transcript) => ({
        ...transcript,
        features: features.filter((feature) => feature.transcript_id === transcript.transcript_id),
      }));
  const sortedModels = [...models].sort((a, b) => {
    const left = String(a.transcript_id ?? a.id ?? '');
    const right = String(b.transcript_id ?? b.id ?? '');
    const byId = left.localeCompare(right, undefined, { numeric: true, sensitivity: 'base' });
    if (byId !== 0) return byId;
    return Number(a.start ?? 0) - Number(b.start ?? 0);
  });
  if (!sortedModels.length) {
    return <EmptyResult text="No transcript model is available for this gene." />;
  }

  const groupedModels = sortedModels.reduce<Map<string, ApiRecord[]>>((groups, model) => {
    const signature = (Array.isArray(model.features) ? model.features : [])
      .filter((feature: ApiRecord) => feature.feature_type === 'exon' || feature.feature_type === 'CDS')
      .map((feature: ApiRecord) => `${feature.feature_type}:${feature.start}-${feature.end}`)
      .sort()
      .join('|');
    const group = groups.get(signature) ?? [];
    group.push(model);
    groups.set(signature, group);
    return groups;
  }, new Map<string, ApiRecord[]>());
  const modelGroups: ApiRecord[][] = Array.from(groupedModels.values());
  const initialModelCount = 6;
  const visibleGroups = showAllModels ? modelGroups : modelGroups.slice(0, initialModelCount);

  const orderedExons = (model: ApiRecord) => (Array.isArray(model.features) ? model.features : [])
    .filter((feature: ApiRecord) => feature.feature_type === 'exon')
    .sort((a: ApiRecord, b: ApiRecord) => {
      const left = Number(a.start ?? 0);
      const right = Number(b.start ?? 0);
      return gene.strand === '-' ? right - left : left - right;
    });

  const min = Number(gene.start ?? 1);
  const max = Number(gene.end ?? min);
  const isSpliced = modelMode === 'spliced';
  const interval = (feature: ApiRecord) => ({ start: Number(feature.start), end: Number(feature.end) });
  const maxSplicedLength = Math.max(1, ...sortedModels.map((model) =>
    spliceExons(orderedExons(model).map(interval), String(gene.strand)).length));
  const width = 1060;
  const labelWidth = 220;
  const rightPad = 42;
  const rowHeight = 84;
  const top = 64;
  const height = top + visibleGroups.length * rowHeight + 18;
  const scale = (value: number) => labelWidth + ((value - min) / Math.max(1, max - min + 1)) * (width - labelWidth - rightPad);
  const splicedScale = (offset: number) => labelWidth + offset / maxSplicedLength * (width - labelWidth - rightPad);
  const trackCenter = labelWidth + (width - labelWidth - rightPad) / 2;

  return (
    <>
    <div className="expression-mode-group" role="group" aria-label="Transcript model view">
      <button type="button" className={!isSpliced ? 'active' : ''} aria-pressed={!isSpliced} onClick={() => setModelMode('genomic')}>Genomic</button>
      <button type="button" className={isSpliced ? 'active' : ''} aria-pressed={isSpliced} onClick={() => setModelMode('spliced')}>Spliced</button>
    </div>
    <div className="model-scroll">
      <svg className="transcript-model" data-mode={modelMode} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${gene.gene_id} ${modelMode} transcript model`}>
        <text x={trackCenter} y="20" className="model-axis-title" textAnchor="middle">{isSpliced ? 'Spliced transcript length (nt), 5′ → 3′' : `${gene.seqid}:${Number(gene.start).toLocaleString()}-${Number(gene.end).toLocaleString()} (${gene.strand})`}</text>
        <line x1={labelWidth} x2={width - rightPad} y1="40" y2="40" className="model-axis" />
        <text x={labelWidth} y="54" className="model-axis-label">{isSpliced ? '0' : Number(gene.start).toLocaleString()}</text>
        <text x={width - rightPad} y="54" className="model-axis-label" textAnchor="end">{isSpliced ? maxSplicedLength.toLocaleString() : Number(gene.end).toLocaleString()}</text>
        {visibleGroups.map((group, index) => {
          const model = group[0];
          const y = top + index * rowHeight;
          const yCenter = y + 46;
          const txStart = Number(model.start ?? gene.start);
          const txEnd = Number(model.end ?? gene.end);
          const modelFeatures = Array.isArray(model.features) ? model.features : [];
          const exons = orderedExons(model);
          const spliced = spliceExons(exons.map(interval), String(gene.strand));
          const cdsFeatures = modelFeatures
            .filter((feature: ApiRecord) => feature.feature_type === 'CDS')
            .sort((a: ApiRecord, b: ApiRecord) => Number(a.start ?? 0) - Number(b.start ?? 0));
          const exonBlocks = exons.map((feature: ApiRecord, featureIndex: number) => {
            const x = isSpliced ? splicedScale(spliced.blocks[featureIndex].offsetStart) : scale(Number(feature.start));
            const blockWidth = isSpliced
              ? splicedScale(spliced.blocks[featureIndex].offsetEnd) - x
              : Math.max(2, scale(Number(feature.end) + 1) - x);
            const exonNumber = featureIndex + 1;
            const exonDetail = `${model.transcript_id}: exon ${exonNumber} of ${exons.length}. ${gene.seqid}:${feature.start}-${feature.end} (${Number(feature.end) - Number(feature.start) + 1} bp), strand ${gene.strand}.`;
            return { feature, featureIndex, x, blockWidth, exonNumber, exonDetail, labelX: x + blockWidth / 2 };
          });
          const labelLanes = new Map<number, number>();
          const laneEnds: number[] = [];
          for (const block of [...exonBlocks].sort((a, b) => a.labelX - b.labelX)) {
            const halfWidth = String(block.exonNumber).length * 4 + 3;
            let lane = laneEnds.findIndex((end) => end < block.labelX - halfWidth);
            if (lane < 0) lane = laneEnds.length;
            laneEnds[lane] = block.labelX + halfWidth;
            labelLanes.set(block.featureIndex, lane);
          }
          const transcriptIds = group.map((item) => String(item.transcript_id));
          const label = group.length === 1 ? transcriptIds[0] : `${transcriptIds[0]} +${group.length - 1}`;
          const details = group.length === 1
            ? `(${exons.length} exons, ${cdsFeatures.length} CDS)`
            : `(${group.length} isoforms · ${exons.length} exons, ${cdsFeatures.length} CDS)`;
          return (
            <g key={transcriptIds.join('|')} data-transcript={model.transcript_id} data-spliced-length={spliced.length}>
              <title>{transcriptIds.join(', ')}</title>
              <rect x={labelWidth} y={y + 3} width={width - labelWidth - rightPad} height={rowHeight - 8} className="model-row-bg" />
              <text x="12" y={yCenter - 3} className="model-transcript-label">{label}</text>
              <text x="12" y={yCenter + 14} className="model-transcript-meta">{details}</text>
              <text x="12" y={yCenter + 31} className="model-transcript-meta">{spliced.length.toLocaleString()} nt (spliced)</text>
              {!isSpliced && <line x1={scale(txStart)} x2={scale(txEnd)} y1={yCenter} y2={yCenter} className="model-intron" />}
              {exonBlocks.map(({ feature, featureIndex, x, blockWidth }) => {
                return (
                  <g key={`exon-${feature.start}-${feature.end}-${featureIndex}`}>
                    <rect
                      x={x}
                      y={yCenter - 5}
                      width={blockWidth}
                      height="10"
                      rx="1.5"
                      className="model-exon"
                    />
                  </g>
                );
              })}
              {exonBlocks.map(({ feature, featureIndex, labelX, exonNumber, exonDetail }) => (
                <text
                  key={`exon-label-${feature.start}-${feature.end}-${featureIndex}`}
                  x={labelX}
                  y={yCenter - 12 - (labelLanes.get(featureIndex) ?? 0) * 14}
                  textAnchor="middle"
                  className="model-exon-index"
                  aria-label={exonDetail}
                >
                  <title>{exonDetail}</title>
                  {exonNumber}
                </text>
              ))}
              {cdsFeatures.flatMap((feature: ApiRecord, featureIndex: number) => {
                const segments = isSpliced
                  ? projectSplicedFeature(interval(feature), spliced.blocks, String(gene.strand)).map((segment) => ({ x: splicedScale(segment.start), blockWidth: splicedScale(segment.end) - splicedScale(segment.start) }))
                  : [{ x: scale(Number(feature.start)), blockWidth: Math.max(2, scale(Number(feature.end) + 1) - scale(Number(feature.start))) }];
                return segments.map(({ x, blockWidth }, segmentIndex) => (
                  <rect
                    key={`cds-${feature.start}-${feature.end}-${featureIndex}-${segmentIndex}`}
                    x={x}
                    y={yCenter - 3.5}
                    width={blockWidth}
                    height="7"
                    rx="1"
                    className="model-cds"
                  />
                ));
              })}
              {isSpliced && exonBlocks.slice(1).map(({ x, featureIndex }) => (
                <line key={`junction-${featureIndex}`} x1={x} x2={x} y1={yCenter - 5} y2={yCenter + 5} className="model-exon-junction" />
              ))}
            </g>
          );
        })}
      </svg>
      {modelGroups.length > initialModelCount ? (
        <button type="button" className="model-expand" onClick={() => setShowAllModels((value) => !value)}>
          <ChevronDown size={16} aria-hidden="true" style={{ transform: showAllModels ? 'rotate(180deg)' : undefined }} />
          {showAllModels ? 'Show fewer transcript models' : `Show all ${modelGroups.length} transcript models`}
        </button>
      ) : null}
    </div>
    </>
  );
}

function SequencePanel({
  gene,
  sequenceData,
  speciesId,
  transcripts,
}: {
  gene: ApiRecord;
  sequenceData: ApiRecord | null;
  speciesId: string;
  transcripts: ApiRecord[];
}) {
  const [mode, setMode] = useState('gene');
  const transcriptOptions = transcripts.filter((item) => item.transcript_id);
  const initialTranscriptId = sequenceData?.views?.transcript?.transcript_id ?? transcriptOptions[0]?.transcript_id ?? '';
  const [selectedTranscriptId, setSelectedTranscriptId] = useState(initialTranscriptId);
  const [localSequenceData, setLocalSequenceData] = useState<ApiRecord | null>(sequenceData);
  const [sequenceError, setSequenceError] = useState('');
  const initialLeftFlank = Number(sequenceData?.upstream ?? 0);
  const initialRightFlank = Number(sequenceData?.downstream ?? 0);
  const [leftFlank, setLeftFlank] = useState(initialLeftFlank);
  const [rightFlank, setRightFlank] = useState(initialRightFlank);
  const [appliedLeftFlank, setAppliedLeftFlank] = useState(initialLeftFlank);
  const [appliedRightFlank, setAppliedRightFlank] = useState(initialRightFlank);
  const [copied, setCopied] = useState(false);

  useEffect(() => {
    const nextLeftFlank = Number(sequenceData?.upstream ?? 0);
    const nextRightFlank = Number(sequenceData?.downstream ?? 0);
    setLocalSequenceData(sequenceData);
    setSelectedTranscriptId(sequenceData?.views?.transcript?.transcript_id ?? transcriptOptions[0]?.transcript_id ?? '');
    setLeftFlank(nextLeftFlank);
    setRightFlank(nextRightFlank);
    setAppliedLeftFlank(nextLeftFlank);
    setAppliedRightFlank(nextRightFlank);
    setSequenceError('');
  }, [sequenceData, gene.gene_id]);

  useEffect(() => {
    const currentTranscriptId = localSequenceData?.views?.transcript?.transcript_id ?? '';
    const currentLeftFlank = Number(localSequenceData?.upstream ?? 0);
    const currentRightFlank = Number(localSequenceData?.downstream ?? 0);
    const needsFetch = selectedTranscriptId !== currentTranscriptId
      || appliedLeftFlank !== currentLeftFlank
      || appliedRightFlank !== currentRightFlank;
    if (!needsFetch) return;
    let cancelled = false;
    const params = new URLSearchParams();
    if (selectedTranscriptId) params.set('transcript_id', selectedTranscriptId);
    if (appliedLeftFlank) params.set('upstream', String(appliedLeftFlank));
    if (appliedRightFlank) params.set('downstream', String(appliedRightFlank));
    setSequenceError('');
    apiGet<ApiRecord>(`/api/species/${speciesId}/genes/${encodeURIComponent(gene.gene_id)}/sequence?${params.toString()}`)
      .then((payload) => {
        if (!cancelled) setLocalSequenceData(payload);
      })
      .catch((error) => {
        if (!cancelled) setSequenceError(error instanceof Error ? error.message : String(error));
      });
    return () => {
      cancelled = true;
    };
  }, [
    selectedTranscriptId,
    speciesId,
    gene.gene_id,
    appliedLeftFlank,
    appliedRightFlank,
    localSequenceData?.upstream,
    localSequenceData?.downstream,
    localSequenceData?.views?.transcript?.transcript_id,
  ]);

  const views = localSequenceData?.views ?? {};
  const availableModes = ['gene', 'reference', 'transcript', 'protein'].filter((item) => views[item]);
  const isReferenceRedundant = String(gene.strand ?? '') === '+' && Boolean(views.gene) && Boolean(views.reference);
  const selectableModes = availableModes.filter((item) => !(item === 'reference' && isReferenceRedundant));
  const activeMode = selectableModes.includes(mode) ? mode : selectableModes[0] ?? availableModes[0] ?? 'gene';
  const sequence = views[activeMode] ?? null;
  const text = String(sequence?.sequence ?? '');
  const canUseFlank = activeMode === 'gene' || activeMode === 'reference';
  const canSelectTranscript = transcriptOptions.length > 1 && (activeMode === 'transcript' || activeMode === 'protein');
  const normalizedLeftFlank = clampSequenceFlank(leftFlank);
  const normalizedRightFlank = clampSequenceFlank(rightFlank);
  const hasPendingFlank = normalizedLeftFlank !== appliedLeftFlank || normalizedRightFlank !== appliedRightFlank;

  async function copySequence() {
    if (!text) return;
    try {
      await navigator.clipboard.writeText(text);
      setCopied(true);
      window.setTimeout(() => setCopied(false), 1600);
    } catch {
      setSequenceError('Sequence could not be copied. Please select the sequence and copy it manually.');
    }
  }

  if (!localSequenceData || !sequence) {
    return <EmptyResult text="Sequence extraction is not available. Check genome FASTA and FAI runtime assets." />;
  }

  return (
    <div className="sequence-panel">
      {sequenceError && <div className="tool-message error">{sequenceError}</div>}
      <div className="seq-mode-group">
        {availableModes.map((item) => {
          const isDisabled = item === 'reference' && isReferenceRedundant;
          return (
            <button
              type="button"
              key={item}
              className={activeMode === item ? 'active' : ''}
              disabled={isDisabled}
              title={isDisabled ? 'Same as Gene sequence for + strand genes.' : undefined}
              onClick={() => {
                if (!isDisabled) setMode(item);
              }}
            >
              {sequenceModeLabel(item)}
            </button>
          );
        })}
      </div>
      <ul className="sequence-term-help">
        <li>Gene sequence follows the gene direction for easier reading.</li>
        <li>Reference strand follows the original genome assembly direction.</li>
      </ul>
      <div className="sequence-meta">
        <dl className="sequence-facts">
          <div>
            <dt>Position</dt>
            <dd className="sequence-position">{sequenceRangeLabel(sequence, gene)}</dd>
          </div>
          <div>
            <dt>Length</dt>
            <dd>{text.length.toLocaleString()} {activeMode === 'protein' ? 'aa' : 'bp'}</dd>
          </div>
        </dl>
        <div className="sequence-controls-export">
        {canUseFlank && (
          <div className="sequence-flank-control" aria-label="Flanking sequence controls">
            <div className="sequence-flank-heading">
              <span>Flanking sequence</span>
            </div>
            <div className="sequence-flank-field">
              <span>Start</span>
              <PrimerStepperControl label="Extend before genomic start in base pairs" value={leftFlank} min={0} max={10000} step={100} unit="bp" onChange={setLeftFlank} />
            </div>
            <div className="sequence-flank-field">
              <span>End</span>
              <PrimerStepperControl label="Extend after genomic end in base pairs" value={rightFlank} min={0} max={10000} step={100} unit="bp" onChange={setRightFlank} />
            </div>
            <button
              type="button"
              className="sequence-flank-apply"
              disabled={!hasPendingFlank}
              onClick={() => {
                setAppliedLeftFlank(normalizedLeftFlank);
                setAppliedRightFlank(normalizedRightFlank);
                setLeftFlank(normalizedLeftFlank);
                setRightFlank(normalizedRightFlank);
              }}
            >
              Apply
            </button>
            <button
              type="button"
              className="sequence-flank-reset"
              disabled={appliedLeftFlank === 0 && appliedRightFlank === 0 && leftFlank === 0 && rightFlank === 0}
              onClick={() => {
                setLeftFlank(0);
                setRightFlank(0);
                setAppliedLeftFlank(0);
                setAppliedRightFlank(0);
              }}
            >
              Reset
            </button>
          </div>
        )}
        {canSelectTranscript && (
          <div className="sequence-transcript-select">
            <span>Transcript</span>
            <CompactSelect variant="primer" ariaLabel="Sequence transcript" value={selectedTranscriptId} onChange={setSelectedTranscriptId} options={transcriptOptions.map((item) => ({ value: item.transcript_id, label: item.transcript_id }))} />
          </div>
        )}
        <span className="sequence-export-tooltip" data-tooltip="GB export">
          <button type="button" className="sequence-gb-export" aria-label="GB export" disabled={!text} onClick={() => {
            const params = new URLSearchParams({ mode: activeMode, upstream: String(localSequenceData.upstream ?? 0), downstream: String(localSequenceData.downstream ?? 0) });
            const displayedTranscript = sequence.transcript_id ?? localSequenceData.views?.transcript?.transcript_id;
            if (displayedTranscript) params.set('transcript_id', String(displayedTranscript));
            window.location.assign(`${API_BASE}/api/species/${encodeURIComponent(speciesId)}/genes/${encodeURIComponent(gene.gene_id)}/sequence.gb?${params}`);
          }}><FileDown size={18} aria-hidden="true" /></button>
        </span>
        </div>
      </div>
      <div className="seq-wrap">
        <button
          type="button"
          className={copied ? 'seq-copy copied' : 'seq-copy'}
          aria-label={copied ? 'Sequence copied' : 'Copy sequence'}
          title={copied ? 'Copied' : 'Copy sequence'}
          onClick={() => void copySequence()}
        >
          {copied ? <Check size={16} strokeWidth={2} /> : <Copy size={16} strokeWidth={1.9} />}
        </button>
        <pre className={`seq ${activeMode === 'protein' ? 'protein-seq' : activeMode === 'transcript' ? 'transcript-seq' : ''}`.trim()}>
          <SequenceText sequence={sequence} />
        </pre>
      </div>
      <ul className="seq-note"><li>{String(sequence.note ?? sequenceModeNote(activeMode))}</li></ul>
    </div>
  );
}

function SequenceText({ sequence }: { sequence: ApiRecord | null }) {
  const text = String(sequence?.sequence ?? '');
  if (!sequence || !text) return <>Sequence not available.</>;
  const marks = new Array(text.length).fill('');
  const features = Array.isArray(sequence.features) ? sequence.features : [];
  features
    .filter((feature: ApiRecord) => feature.feature_type === 'exon' || feature.feature_type === 'CDS')
    .forEach((feature: ApiRecord) => {
      const className = feature.feature_type === 'CDS' ? 'cds' : 'exon';
      const start = Math.max(0, Number(feature.display_start ?? feature.start ?? 1) - 1);
      const end = Math.min(text.length - 1, Number(feature.display_end ?? feature.end ?? 0) - 1);
      for (let index = start; index <= end; index += 1) marks[index] = className;
    });

  const chunks: React.ReactNode[] = [];
  let buffer = '';
  let current = marks[0] || '';
  const flush = (className: string, key: number) => {
    if (!buffer) return;
    chunks.push(className ? <span className={className} key={key}>{buffer}</span> : <span key={key}>{buffer}</span>);
    buffer = '';
  };
  for (let index = 0; index < text.length; index += 1) {
    const mark = marks[index] || '';
    if (mark !== current) {
      flush(current, index);
      current = mark;
    }
    buffer += text[index];
  }
  flush(current, text.length);
  return <>{chunks}</>;
}

function sequenceModeLabel(mode: string) {
  if (mode === 'gene') return 'Gene sequence';
  if (mode === 'reference') return 'Reference strand';
  if (mode === 'transcript') return 'Transcript';
  if (mode === 'protein') return 'Protein';
  return mode;
}

function clampSequenceFlank(value: number) {
  if (!Number.isFinite(value)) return 0;
  return Math.max(0, Math.min(10000, Math.round(value)));
}

function sequenceRangeLabel(sequence: ApiRecord, gene: ApiRecord): React.ReactNode {
  if (sequence.mode === 'protein') return sequence.protein_id ? `${sequence.protein_id} protein sequence` : 'Protein sequence';
  if (sequence.mode === 'transcript') return sequence.transcript_id ? `${sequence.transcript_id} spliced transcript` : 'Spliced transcript';
  const seqid = String(sequence.seqid ?? gene.seqid ?? '');
  const start = Number(sequence.start ?? gene.start).toLocaleString();
  const end = Number(sequence.end ?? gene.end).toLocaleString();
  return (
    <>
      <strong>{seqid}</strong>
      <span className="sequence-range">{start} - {end}</span>
    </>
  );
}

function sequenceModeNote(mode: string) {
  if (mode === 'gene') return 'Gene sequence is shown in gene orientation, with exon and CDS bases highlighted.';
  if (mode === 'reference') return 'Reference strand is shown as stored in the genome assembly.';
  if (mode === 'transcript') return 'Transcript sequence is exon-spliced and displayed in transcript direction.';
  if (mode === 'protein') return 'Protein sequence is loaded from the species protein FASTA build.';
  return '';
}

function DomainEvidencePanel({
  domainData,
  geneId,
  proteins,
  speciesId,
}: {
  domainData: ApiRecord | null;
  geneId: string;
  proteins: ApiRecord[];
  speciesId: string;
}) {
  const proteinOptions = [...proteins].sort((a, b) => String(a.transcript_id ?? '').localeCompare(
    String(b.transcript_id ?? ''), undefined, { numeric: true, sensitivity: 'base' },
  ));
  const defaultProteinId = String(proteinOptions[0]?.protein_id ?? '');
  const [selectedProteinId, setSelectedProteinId] = useState(defaultProteinId);
  const [selectedDomainData, setSelectedDomainData] = useState<ApiRecord | null>(domainData);
  const [loadingProteinId, setLoadingProteinId] = useState('');
  const [domainError, setDomainError] = useState('');
  const activeProteinId = proteinOptions.some((protein) => protein.protein_id === selectedProteinId)
    ? selectedProteinId
    : defaultProteinId;

  useEffect(() => {
    setSelectedProteinId(defaultProteinId);
    setSelectedDomainData(domainData);
    setLoadingProteinId('');
    setDomainError('');
  }, [geneId, speciesId]);

  useEffect(() => {
    if (!activeProteinId) return;
    if (domainData?.protein_id === activeProteinId) {
      setSelectedDomainData(domainData);
      setDomainError('');
      return;
    }
    let cancelled = false;
    setLoadingProteinId(activeProteinId);
    setDomainError('');
    void apiGet<ApiRecord>(`/api/species/${speciesId}/domains/${encodeURIComponent(geneId)}?protein_id=${encodeURIComponent(activeProteinId)}&cached_only=true`)
      .then((result) => {
        if (!cancelled) setSelectedDomainData(result);
      })
      .catch((error) => {
        if (!cancelled) setDomainError(error instanceof Error ? error.message : String(error));
      })
      .finally(() => {
        if (!cancelled) setLoadingProteinId('');
      });
    return () => {
      cancelled = true;
    };
  }, [activeProteinId, domainData, geneId, speciesId]);

  const activeDomainData = selectedDomainData?.protein_id === activeProteinId ? selectedDomainData : null;
  const domains = Array.isArray(activeDomainData?.domains) ? activeDomainData.domains : [];
  const proteinLength = Math.max(1, Number(activeDomainData?.protein_length ?? proteinOptions.find((protein) => protein.protein_id === activeProteinId)?.length ?? 1));
  if (!proteinOptions.length) {
    return <EmptyResult text="No translated protein sequence is available for this gene." />;
  }
  return (
    <>
      {proteinOptions.length > 1 ? (
        <div className="domain-protein-select">
          <span>Protein isoform</span>
          <CompactSelect variant="primer" ariaLabel="Domain protein isoform" value={activeProteinId} onChange={setSelectedProteinId} options={proteinOptions.map((protein) => ({ value: protein.protein_id, label: `${protein.transcript_id} · ${Number(protein.length ?? 0).toLocaleString()} aa` }))} />
        </div>
      ) : null}
      {loadingProteinId ? <p className="domain-loading">Loading precomputed Pfam domains...</p> : null}
      {domainError ? <EmptyResult text={domainError} /> : null}
      {!domainError && !activeDomainData && !loadingProteinId ? <EmptyResult text="Waiting for cached Pfam domains." /> : null}
      {!domainError && activeDomainData && !domains.length ? <EmptyResult text="No Pfam domains passed the current i-Evalue cutoff." /> : null}
      {!domainError && activeDomainData && domains.length ? <>
        <DomainMap domainData={activeDomainData} proteinLength={proteinLength} />
        <CompactTable
          title="Pfam domains"
          titleActions={<div className="summary-gene-actions expression-mode-group">
            <button type="button" title="Download all stored Pfam results for this protein" onClick={() => window.location.assign(`${API_BASE}/api/species/${encodeURIComponent(speciesId)}/domains/${encodeURIComponent(geneId)}/export?protein_id=${encodeURIComponent(activeProteinId)}`)}><FileDown size={15} aria-hidden="true" />TSV export</button>
          </div>}
          columns={['domain', 'accession', 'protein_range', 'hmm_range', 'i_evalue', 'domain_score', 'description']}
          pageSize={5}
          tableClassName="domain-table"
          rows={domains.map((domain: ApiRecord) => ({
            ...domain,
            protein_range: `${domain.ali_from ?? domain.start ?? '-'}-${domain.ali_to ?? domain.end ?? '-'}`,
            hmm_range: `${domain.hmm_from ?? '-'}-${domain.hmm_to ?? '-'}`,
          }))}
        />
      </> : null}
    </>
  );
}

function DomainMap({ domainData, proteinLength }: { domainData: ApiRecord; proteinLength: number }) {
  const domains = Array.isArray(domainData.domains) ? domainData.domains : [];
  const width = 920;
  const height = Math.max(112, 78 + domains.length * 26);
  const trackW = 730;
  const trackX = (width - trackW) / 2;
  const colors = ['#126b5b', '#b14b36', '#5b6f95', '#8a6d1d', '#6d5a7b', '#44736a'];
  const scale = (aa?: number) => trackX + ((Number(aa || 1) - 1) / proteinLength) * trackW;
  return (
    <svg className="domain-map" viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${domainData.query_transcript || domainData.gene_id} Pfam domain map`}>
      <text x={trackX + trackW / 2} y="17" className="domain-map-title" textAnchor="middle">{domainData.query_transcript ?? domainData.protein_id} ({proteinLength.toLocaleString()} aa)</text>
      <line x1={trackX} x2={trackX + trackW} y1="38" y2="38" className="domain-backbone" />
      <text x={trackX} y="30" className="domain-map-axis">1</text>
      <text x={trackX + trackW} y="30" className="domain-map-axis" textAnchor="end">{proteinLength.toLocaleString()}</text>
      {domains.map((domain: ApiRecord, index: number) => {
        const start = Number(domain.ali_from ?? domain.start ?? 1);
        const end = Number(domain.ali_to ?? domain.end ?? start);
        const x = scale(start);
        const blockWidth = Math.max(4, scale(end) - x);
        const y = 60 + index * 26;
        return (
          <g key={`${domain.domain ?? domain.domain_name}-${index}`}>
            <rect x={x} y={y} width={blockWidth} height="9" rx="2" fill={colors[index % colors.length]} />
            <text x={Math.min(x + blockWidth + 5, width - 12 - String(domain.domain ?? domain.domain_name).length * 6)} y={y + 8} className="domain-label">{domain.domain ?? domain.domain_name}</text>
          </g>
        );
      })}
    </svg>
  );
}

function ReferenceSimilarityPanel({ similarHits, topHit, exportUrl }: { similarHits: ApiRecord[]; topHit: ApiRecord | undefined; exportUrl: string }) {
  return (
    <>
      {similarHits.length ? (
        <CompactTable
          title={topHit?.reference_display_name ?? 'Top reference hits'}
          titleActions={<div className="summary-gene-actions expression-mode-group">
            <button type="button" title="Download all stored reference similarity results for this protein" onClick={() => window.location.assign(exportUrl)}><FileDown size={15} aria-hidden="true" />TSV export</button>
          </div>}
          columns={['rank', 'reference_gene_id', 'reference_transcript_id', 'reference_gene_symbol', 'pident', 'alignment_length', 'evalue', 'bitscore', 'reference_protein_description']}
          pageSize={5}
          tableClassName="similarity-table"
          columnLabels={{ reference_gene_id: <>Reference<br />Gene ID</>, bitscore: <>Bit<br />Score</> }}
          rows={similarHits.slice(0, 5).map((hit) => ({
            ...hit,
            reference_protein_description: compactReferenceDescription(hit.reference_protein_description),
          }))}
        />
      ) : (
        <EmptyResult text="No precomputed reference similarity hit is available for this gene." />
      )}
    </>
  );
}

function ProteinSimilarityPanel({ geneId, proteins, speciesId }: { geneId: string; proteins: ApiRecord[]; speciesId: string }) {
  const options = [...proteins].sort((a, b) => String(a.transcript_id ?? '').localeCompare(String(b.transcript_id ?? ''), undefined, { numeric: true, sensitivity: 'base' }));
  const defaultProteinId = String(options[0]?.protein_id ?? '');
  const [selectedProteinId, setSelectedProteinId] = useState(defaultProteinId);
  const [data, setData] = useState<ApiRecord | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const activeProteinId = options.some((item) => item.protein_id === selectedProteinId) ? selectedProteinId : defaultProteinId;

  useEffect(() => {
    setSelectedProteinId(defaultProteinId);
    setData(null);
    setError('');
  }, [geneId, speciesId]);

  useEffect(() => {
    if (!activeProteinId) return;
    let cancelled = false;
    setLoading(true); setError('');
    void apiGet<ApiRecord>(`/api/species/${speciesId}/similar/${encodeURIComponent(geneId)}?protein_id=${encodeURIComponent(activeProteinId)}&limit=5`)
      .then((result) => { if (!cancelled) setData(result); })
      .catch((requestError) => { if (!cancelled) setError(requestError instanceof Error ? requestError.message : String(requestError)); })
      .finally(() => { if (!cancelled) setLoading(false); });
    return () => { cancelled = true; };
  }, [activeProteinId, geneId, speciesId]);

  if (!options.length) return <EmptyResult text="No translated protein sequence is available for similarity search." />;
  const hits = Array.isArray(data?.hits) ? data.hits : [];
  return <>
    {options.length > 1 ? <div className="domain-protein-select">
      <span>Protein isoform</span>
      <CompactSelect variant="primer" ariaLabel="Similarity protein isoform" value={activeProteinId} onChange={setSelectedProteinId} options={options.map((protein) => ({ value: protein.protein_id, label: `${protein.transcript_id} · ${Number(protein.length ?? 0).toLocaleString()} aa` }))} />
    </div> : null}
    {loading ? <p className="domain-loading">Searching the reference proteome...</p> : null}
    {error ? <EmptyResult text={error} /> : null}
    {!loading && !error && data?.protein_id === activeProteinId ? <ReferenceSimilarityPanel similarHits={hits} topHit={hits[0]} exportUrl={`${API_BASE}/api/species/${encodeURIComponent(speciesId)}/similar/${encodeURIComponent(geneId)}/export?protein_id=${encodeURIComponent(activeProteinId)}`} /> : null}
  </>;
}

function ReferencePhylogenyPanel({ gene, proteins, speciesId }: { gene: ApiRecord; proteins: ApiRecord[]; speciesId: string }) {
  const panelRef = useRef<HTMLDivElement>(null);
  const [exportError, setExportError] = useState('');
  const [exporting, setExporting] = useState(false);
  async function downloadTree(format: 'png' | 'svg') {
    const svg = panelRef.current?.querySelector<SVGSVGElement>('svg.reference-tree');
    if (!svg) return;
    setExporting(true);
    setExportError('');
    try {
      await exportTree(svg, format, `${speciesId}_${gene.gene_id}_phylogeny`);
    } catch (error) {
      setExportError(error instanceof Error ? error.message : 'Tree export failed.');
    } finally {
      setExporting(false);
    }
  }
  const [model, setModel] = useState('LG+G');
  const [analysis, setAnalysis] = useState<'reference' | 'isoform'>('reference');
  const sortedProteins = [...proteins].sort((a, b) => Number(b.length ?? 0) - Number(a.length ?? 0));
  const defaultProteinId = String(sortedProteins[0]?.protein_id ?? '');
  const [proteinId, setProteinId] = useState(defaultProteinId);
  const [job, setJob] = useState<ApiRecord | null>(null);
  const [error, setError] = useState('');
  const requestToken = useRef(0);

  async function requestTree(nextModel: string, nextAnalysis = analysis, nextProteinId = proteinId) {
    const token = ++requestToken.current;
    setError('');
    setJob({ status: 'queued', model: nextModel });
    try {
      const queued = await apiPost<ApiRecord>('/api/phylogeny/tree', {
        species_id: speciesId,
        gene_id: gene.gene_id,
        model: nextModel,
        analysis: nextAnalysis,
        protein_id: nextAnalysis === 'reference' ? nextProteinId : null,
      });
      if (token !== requestToken.current) return;
      setJob(queued);
      for (let attempt = 0; attempt < 80; attempt += 1) {
        await new Promise((resolve) => window.setTimeout(resolve, attempt === 0 ? 250 : 700));
        const current = await apiGet<ApiRecord>(`/api/phylogeny/jobs/${queued.job_id}`);
        if (token !== requestToken.current) return;
        setJob(current);
        if (current.status === 'completed' || current.status === 'failed') {
          if (current.status === 'failed') setError(String(current.error ?? 'Phylogenetic tree calculation failed.'));
          return;
        }
      }
      if (token === requestToken.current) setError('Phylogenetic tree calculation timed out.');
    } catch (requestError) {
      if (token === requestToken.current) setError(requestError instanceof Error ? requestError.message : String(requestError));
    }
  }

  useEffect(() => {
    if (!proteins.length) return;
    setModel('LG+G');
    setAnalysis('reference');
    setProteinId(defaultProteinId);
    void requestTree('LG+G', 'reference', defaultProteinId);
  }, [gene.gene_id]);

  if (!proteins.length) return <EmptyResult text="No translated protein sequence is available for tree visualization." />;
  const result = job?.status === 'completed' && job.results ? job.results as ApiRecord : null;
  return (
    <div className="phylogeny-panel" ref={panelRef}>
      <div className="gene-tree-toolbar">
      <div className="phylogeny-model-select tree-analysis-select">
        <span>Analysis</span>
        <CompactSelect variant="primer" ariaLabel="Tree analysis" value={analysis} options={[{ value: 'reference', label: 'Reference similarity' }, { value: 'isoform', label: 'Isoform comparison' }]} onChange={(value) => { const next = value as 'reference' | 'isoform'; setAnalysis(next); void requestTree(model, next, proteinId); }} />
      </div>
      {analysis === 'reference' ? <div className="phylogeny-model-select">
        <span>Query transcript</span>
        <CompactSelect variant="primer" ariaLabel="Tree query transcript" value={proteinId} options={sortedProteins.map((protein) => ({ value: protein.protein_id, label: `${protein.transcript_id} · ${Number(protein.length ?? 0).toLocaleString()} aa` }))} onChange={(next) => { setProteinId(next); void requestTree(model, analysis, next); }} />
      </div> : null}
      <div className="phylogeny-model-select tree-model-select">
        <span>Model</span>
        <CompactSelect variant="primer" ariaLabel="Tree protein model"
          value={model}
          disabled={job?.status === 'queued' || job?.status === 'running'}
          onChange={(nextModel) => {
            setModel(nextModel);
            void requestTree(nextModel);
          }}
          options={['LG+G', 'WAG+G', 'JTT+G', 'Auto'].map((value) => ({ value, label: value }))}
        />
      </div>
      <div className="tree-export-actions">
        {(['png', 'svg'] as const).map((format) => <span key={format} className="go-filter-tooltip" data-tooltip={`${format.toUpperCase()} export`}><button type="button" className="go-bubble-export-button" aria-label={`Export tree as ${format.toUpperCase()}`} disabled={!result || exporting} onClick={() => void downloadTree(format)}>{format.toUpperCase()}</button></span>)}
      </div>
      </div>
      {exportError ? <p role="alert" className="tool-message error">{exportError}</p> : null}
      {error ? <div className="tool-message error">{error}</div> : null}
      {!error && !result ? <p className="phylogeny-status">{job?.status === 'running' ? 'Aligning proteins and inferring tree...' : 'Preparing protein tree...'}</p> : null}
      {result ? <ProteinPhylogenyTree gene={gene} result={result} /> : null}
    </div>
  );
}

type NewickNode = {
  children: NewickNode[];
  label: string;
  length: number;
};

type PositionedNewickNode = {
  depth: number;
  y: number;
  label: string;
  length: number;
  children: PositionedNewickNode[];
};

function parseNewick(input: string): NewickNode {
  let cursor = 0;
  const skipWhitespace = () => {
    while (/\s/.test(input[cursor] ?? '')) cursor += 1;
  };
  const readUntilDelimiter = () => {
    const start = cursor;
    while (cursor < input.length && !'():,;'.includes(input[cursor])) cursor += 1;
    return input.slice(start, cursor).trim();
  };
  const readBranch = (): NewickNode => {
    skipWhitespace();
    const children: NewickNode[] = [];
    if (input[cursor] === '(') {
      cursor += 1;
      do {
        children.push(readBranch());
        skipWhitespace();
        if (input[cursor] !== ',') break;
        cursor += 1;
      } while (cursor < input.length);
      if (input[cursor] === ')') cursor += 1;
    }
    skipWhitespace();
    const label = readUntilDelimiter();
    let length = 0;
    if (input[cursor] === ':') {
      cursor += 1;
      const parsed = Number.parseFloat(readUntilDelimiter());
      length = Number.isFinite(parsed) ? parsed : 0;
    }
    return { children, label, length };
  };
  return readBranch();
}

function ProteinPhylogenyTree({ gene, result }: { gene: ApiRecord; result: ApiRecord }) {
  let root: NewickNode;
  try {
    root = parseNewick(String(result.newick ?? ''));
  } catch {
    return <EmptyResult text="The generated Newick tree could not be rendered." />;
  }
  const leaves = Array.isArray(result.leaves) ? result.leaves as ApiRecord[] : [];
  const displayLabels = new Map(leaves.map((item) => [String(item.id), String(item.label)]));
  const leafKinds = new Map(leaves.map((item) => [String(item.id), String(item.kind)]));
  let leafIndex = 0;
  let maxDepth = 0;

  const placeTree = (node: NewickNode, parentDepth: number): PositionedNewickNode => {
    const depth = parentDepth + Math.max(0, node.length);
    maxDepth = Math.max(maxDepth, depth);
    const children = node.children.map((child) => placeTree(child, depth));
    const y = children.length
      ? children.reduce((sum, child) => sum + child.y, 0) / children.length
      : leafIndex++;
    return { ...node, depth, y, children };
  };
  const tree = placeTree(root, 0);
  const leafCount = Math.max(leafIndex, 1);
  const width = 1040;
  const treeX = 56;
  const treeW = 590;
  const labelX = 678;
  const topY = 40;
  const rowH = 44;
  const height = Math.max(180, topY + leafCount * rowH + 58);
  const normalizedDepth = Math.max(maxDepth, 0.01);
  const xForDepth = (depth: number) => treeX + (depth / normalizedDepth) * treeW;
  const yForIndex = (index: number) => topY + index * rowH;
  const scaleDistance = normalizedDepth >= 0.5 ? 0.1 : normalizedDepth >= 0.2 ? 0.05 : normalizedDepth >= 0.08 ? 0.02 : 0.01;
  const scaleWidth = Math.max(28, Math.min(treeW, (scaleDistance / normalizedDepth) * treeW));

  const renderNode = (node: PositionedNewickNode): React.ReactNode => {
    const x = xForDepth(node.depth);
    const y = yForIndex(node.y);
    if (!node.children.length) {
      const kind = leafKinds.get(node.label);
      return (
        <g key={`${node.label}-${node.depth}-${node.y}`}>
          <circle cx={x} cy={y} r="4.5" className={kind === 'query' ? 'reference-tree-query' : 'reference-tree-hit'} />
          <text x={labelX} y={y + 4} className={kind === 'query' ? 'reference-tree-query-label' : 'reference-tree-hit-label'}>
            {displayLabels.get(node.label) ?? node.label}
          </text>
        </g>
      );
    }
    const childYs = node.children.map((child) => yForIndex(child.y));
    return (
      <g key={`${node.label || 'internal'}-${node.depth}-${node.y}`}>
        <line x1={x} x2={x} y1={Math.min(...childYs)} y2={Math.max(...childYs)} className="reference-tree-spine" />
        {node.children.map((child) => {
          const childX = xForDepth(child.depth);
          const childY = yForIndex(child.y);
          return <line key={`${child.label}-${child.depth}-${child.y}-branch`} x1={x} x2={childX} y1={childY} y2={childY} className="reference-tree-branch" />;
        })}
        {node.children.map(renderNode)}
      </g>
    );
  };

  return (
    <svg className="reference-tree" viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${gene.gene_id} reference similarity tree`}>
      <text x={treeX} y="20" className="reference-tree-title">{result.selected_model ?? result.requested_model} · {Number(result.alignment_length ?? 0).toLocaleString()} aligned aa</text>
      {renderNode(tree)}
      <line x1={treeX} x2={treeX + scaleWidth} y1={height - 28} y2={height - 28} className="reference-tree-scale" />
      <line x1={treeX} x2={treeX} y1={height - 33} y2={height - 23} className="reference-tree-scale" />
      <line x1={treeX + scaleWidth} x2={treeX + scaleWidth} y1={height - 33} y2={height - 23} className="reference-tree-scale" />
      <text x={treeX} y={height - 10} className="reference-tree-scale-label">{scaleDistance} substitutions/site</text>
    </svg>
  );
}

function MultiGeneReport({
  matrixData,
  goById,
  similarById,
}: {
  matrixData: ApiRecord;
  goById: Record<string, ApiRecord>;
  similarById: Record<string, ApiRecord>;
}) {
  const genes = Array.isArray(matrixData.genes) ? matrixData.genes : [];
  const samples = Array.isArray(matrixData.samples) ? matrixData.samples : [];
  const values = matrixData.values ?? {};
  const summaryRows = genes.map((gene: ApiRecord) => {
    const geneId = gene.gene_id;
    const sampleValues = Object.values(values[geneId] ?? {}).filter((value) => typeof value === 'number') as number[];
    const meanCpm = sampleValues.length ? sampleValues.reduce((sum, value) => sum + value, 0) / sampleValues.length : null;
    const topHit = Array.isArray(similarById[geneId]?.hits) ? similarById[geneId].hits[0] : null;
    const goTerms = Array.isArray(goById[geneId]?.go_terms) ? goById[geneId].go_terms : [];
    return {
      gene_id: geneId,
      location: formatLocation(gene),
      mean_cpm: meanCpm,
      top_reference_gene: topHit?.reference_gene_id,
      pident: topHit?.pident,
      go_terms: goTerms.length,
    };
  });
  const matrixRows = genes.map((gene: ApiRecord) => {
    const row: ApiRecord = { gene_id: gene.gene_id };
    samples.forEach((sample: ApiRecord) => {
      row[sample.sample_id] = values[gene.gene_id]?.[sample.sample_id] ?? null;
    });
    return row;
  });

  return (
    <div className="gene-report">
      <section className="gene-report-section">
        <div className="gene-section-head">
          <span>Multi-gene result</span>
          <h4>{genes.length} genes · {samples.length} samples</h4>
        </div>
        {Array.isArray(matrixData.missing_genes) && matrixData.missing_genes.length > 0 && (
          <div className="tool-message error">Missing genes: {matrixData.missing_genes.join(', ')}</div>
        )}
        <CompactTable
          title="Gene summary"
          columns={['gene_id', 'location', 'mean_cpm', 'top_reference_gene', 'pident', 'go_terms']}
          rows={summaryRows}
        />
      </section>
      <section className="gene-report-section">
        <div className="gene-section-head">
          <span>Expression</span>
          <h4>CPM matrix</h4>
        </div>
        <CompactTable
          title="Expression matrix"
          columns={['gene_id', ...samples.map((sample: ApiRecord) => sample.sample_id)]}
          rows={matrixRows}
        />
      </section>
      <section className="gene-report-section">
        <div className="gene-section-head">
          <span>Domain evidence</span>
          <h4>Not yet available</h4>
        </div>
        <EmptyResult text="Protein domain comparison requires query-side domain annotation. This section is reserved for the Pfam/InterPro build step." />
      </section>
    </div>
  );
}

type ExpressionGroup = {
  key: string;
  label: string;
  mean: number;
  values: number[];
  rows: ApiRecord[];
};

function expressionSummaryForReport(expressionData: ApiRecord | null) {
  const cpmRows: ApiRecord[] = Array.isArray(expressionData?.expression?.cpm) ? expressionData.expression.cpm : [];
  if (!cpmRows.length) return 'No expression data.';
  const pattern = classifyExpressionPattern(groupExpressionRows(cpmRows));
  const cutoff = EXPRESSION_PATTERN_POLICY.minimumGroupCpm;
  const detected = cpmRows.filter((row) => Number(row.value ?? 0) >= cutoff).length;
  return `${pattern.label} (CPM group means); detected in ${detected}/${cpmRows.length} samples at CPM >= ${cutoff}.`;
}

function expressionPeakForSearch(expressionData: ApiRecord | null | undefined) {
  const cpmRows: ApiRecord[] = Array.isArray(expressionData?.expression?.cpm) ? expressionData.expression.cpm : [];
  if (!cpmRows.length) return '-';
  return classifyExpressionPattern(groupExpressionRows(cpmRows)).label;
}

function topPfamForSearch(item: GeneSearchResult) {
  const domains: ApiRecord[] = Array.isArray(item.domains?.domains) ? item.domains.domains : [];
  if (domains.length) {
    const topDomain = domains[0];
    return topDomain.domain_name ?? topDomain.domain ?? topDomain.pfam_id ?? '-';
  }
  const summary = item.gene?.cache_summary ?? item.gene?.summary ?? {};
  return summary.top_domain ?? summary.top_pfam ?? summary.pfam ?? '-';
}

function compactReferenceDescription(value: unknown) {
  const text = String(value ?? '').trim();
  if (!text) return '-';

  // Ensembl/TAIR FASTA headers carry structured metadata before the functional description.
  const descriptionMatch = text.match(/(?:^|\s)description:([\s\S]*?)(?:\s+\[Source:|$)/i);
  const displayText = (descriptionMatch?.[1] ?? text).trim();
  return displayText.length > 180 ? `${displayText.slice(0, 177).trimEnd()}...` : displayText;
}

function functionalHintsForReport(topHit: ApiRecord | undefined, domainData: ApiRecord | null) {
  const domains: ApiRecord[] = Array.isArray(domainData?.domains) ? domainData.domains : [];
  const hints = [
    domains[0]?.domain ? `Pfam: ${domains[0].domain}` : '',
    topHit?.reference_gene_id ? `Reference: ${topHit.reference_gene_id}` : '',
    topHit?.pident ? `Identity: ${formatCell(topHit.pident)}%` : '',
  ].filter(Boolean);
  return hints.length ? hints.join(' · ') : 'No domain or reference hint is available yet.';
}

function median(values: number[]) {
  if (!values.length) return 0;
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
}

function inferExpressionGroup(row: ApiRecord) {
  const sampleId = String(row.sample_id ?? row.display_name ?? '');
  const cqMatch = sampleId.match(/^(Cq|Inf)_(4|16)_(\d+)$/i);
  if (cqMatch) {
    const treatmentLabel = cqMatch[1].toLowerCase() === 'inf' ? 'Infection' : 'Control';
    const temperatureLabel = `${cqMatch[2]}°C`;
    return {
      key: `${cqMatch[1]}_${cqMatch[2]}`,
      label: `${treatmentLabel} ${temperatureLabel}`,
      order: `${cqMatch[1] === 'Cq' ? 0 : 1}-${cqMatch[2] === '4' ? 0 : 1}`,
    };
  }
  const fields = [row.tissue, row.condition, row.treatment, row.timepoint].filter(Boolean);
  const label = fields.length ? fields.join(' ') : String(row.display_name ?? row.sample_id ?? 'Sample');
  return { key: label, label, order: label };
}

function groupExpressionRows(rows: ApiRecord[]): ExpressionGroup[] {
  const grouped = new Map<string, ExpressionGroup & { order: string }>();
  rows.forEach((row) => {
    const group = inferExpressionGroup(row);
    const value = Number(row.value ?? 0);
    if (!grouped.has(group.key)) {
      grouped.set(group.key, { key: group.key, label: group.label, mean: 0, values: [], rows: [], order: group.order });
    }
    const current = grouped.get(group.key);
    if (!current) return;
    current.values.push(value);
    current.rows.push(row);
  });
  return [...grouped.values()]
    .map((group) => ({
      ...group,
      mean: group.values.length ? group.values.reduce((sum, value) => sum + value, 0) / group.values.length : 0,
    }))
    .sort((a, b) => a.order.localeCompare(b.order, undefined, { numeric: true }));
}

function expressionInterpretation(groups: ExpressionGroup[], activeMode: 'cpm' | 'median' | 'raw_count') {
  const means = groups.map((group) => Number(group.mean || 0)).filter((value) => Number.isFinite(value) && value >= 0);
  const maxMean = Math.max(...means, 0);
  const dominant = groups.reduce<ExpressionGroup | null>((best, group) => (
    !best || group.mean > best.mean ? group : best
  ), null);
  const specificity = means.length > 1 && maxMean > 0
    ? means.reduce((sum, value) => sum + (1 - value / maxMean), 0) / (means.length - 1)
    : null;
  const breadthCutoff = EXPRESSION_PATTERN_POLICY.minimumGroupCpm;
  const breadth = groups.filter((group) => activeMode === 'raw_count' ? group.mean > 0 : group.mean >= breadthCutoff).length;
  return {
    dominantLabel: dominant && dominant.mean > 0 ? dominant.label : '-',
    specificity,
    breadth,
    breadthTotal: groups.length,
    breadthCutoffLabel: activeMode === 'raw_count' ? '> 0 raw count' : `>= ${breadthCutoff} CPM`,
  };
}

function ExpressionProfile({ expressionData, speciesId }: { expressionData: ApiRecord | null; speciesId: string }) {
  const expression = expressionData?.expression ?? {};
  const hasCpm = Array.isArray(expression.cpm) && expression.cpm.length > 0;
  const hasRaw = Array.isArray(expression.raw_count) && expression.raw_count.length > 0;
  const [mode, setMode] = useState<'cpm' | 'median' | 'raw_count'>('cpm');
  const availableModes: Array<'cpm' | 'median' | 'raw_count'> = [
    ...(hasCpm ? ['cpm' as const] : []),
    ...(hasCpm || hasRaw ? ['median' as const] : []),
    ...(hasRaw ? ['raw_count' as const] : []),
  ];
  const activeMode = availableModes.includes(mode) ? mode : availableModes[0] ?? 'cpm';
  const baseRows = activeMode === 'raw_count'
    ? expression.raw_count
    : hasCpm
      ? expression.cpm
      : expression.raw_count;
  const rows = Array.isArray(baseRows) ? baseRows : [];

  if (!rows.length) {
    return <EmptyResult text="No expression count record is available for this gene." />;
  }

  const numericRows: ApiRecord[] = rows.map((item: ApiRecord) => ({ ...item, value: Number(item.value ?? 0) }));
  const medianBase = median(numericRows.map((item) => item.value).filter((value) => Number.isFinite(value) && value > 0));
  const plottedRows = numericRows.map((item) => ({
    ...item,
    value: activeMode === 'median' && medianBase > 0 ? item.value / medianBase : item.value,
  }));
  const groups = groupExpressionRows(plottedRows);
  const replicateGroups = groups.filter((group) => group.values.length > 1).length;
  const yLabel = activeMode === 'raw_count' ? 'Raw count' : activeMode === 'median' ? 'Median-scaled value' : 'CPM';
  const contextGroups = groupExpressionRows(hasCpm ? expression.cpm : expression.raw_count);
  const pattern = classifyExpressionPattern(hasCpm ? contextGroups : []);
  const interpretation = expressionInterpretation(contextGroups, hasCpm ? 'cpm' : 'raw_count');
  const specificityLabel = interpretation.specificity === null ? '-' : interpretation.specificity.toFixed(2);

  return (
    <div className="expression-profile">
      <div className="expression-main">
        <div className="expression-mode-group" role="group" aria-label="Expression value mode">
          {availableModes.map((item) => (
            <button type="button" key={item} className={activeMode === item ? 'active' : ''} onClick={() => setMode(item)}>
              {item === 'raw_count' ? 'Raw' : item === 'median' ? 'Median' : 'CPM'}
            </button>
          ))}
        </div>
        <ExpressionPlot groups={groups} yLabel={yLabel} multilineLabels={speciesId === 'Da'} />
      </div>
      <aside className="expression-context">
        <div className="expression-context-head">
          <h5>Expression context</h5>
          <span>{hasCpm ? 'CPM' : 'Raw count'}</span>
        </div>
        <div className="expression-context-primary">
          <span>Expression pattern</span>
          <strong title={pattern.description}>{pattern.label}</strong>
        </div>
        <div className="expression-context-stats">
          <div>
            <span>Specificity</span>
            <strong>{specificityLabel}</strong>
          </div>
          <div>
            <span>Breadth</span>
            <strong>{interpretation.breadth} / {interpretation.breadthTotal}</strong>
          </div>
          <div>
            <span>Groups</span>
            <strong>{groups.length}</strong>
          </div>
          <div>
            <span>Samples</span>
            <strong>{numericRows.length}</strong>
          </div>
          <div>
            <span>Replicate groups</span>
            <strong>{replicateGroups || 'None'}</strong>
          </div>
        </div>
        <div className="expression-context-note">
          <span>{pattern.description}</span>
          <span>Tau-like specificity from group means.</span>
          <span>Breadth cutoff: {interpretation.breadthCutoffLabel}.</span>
        </div>
      </aside>
    </div>
  );
}

function ExpressionPlot({ groups, yLabel, multilineLabels = false }: { groups: ExpressionGroup[]; yLabel: string; multilineLabels?: boolean }) {
  const width = 900;
  const labelLines = groups.map((group) => multilineLabels ? group.label.trim().split(/\s+/) : [group.label]);
  const extraLabelHeight = Math.max(0, ...labelLines.map((lines) => lines.length - 1)) * 22;
  const height = 430 + extraLabelHeight;
  const margin = { top: 24, right: 28, bottom: 64 + extraLabelHeight, left: 86 };
  const chartWidth = width - margin.left - margin.right;
  const chartHeight = height - margin.top - margin.bottom;
  const maxValue = Math.max(...groups.flatMap((group) => [...group.values, group.mean]), 0);
  const { ticks, yMax } = expressionAxisTicks(maxValue);
  const xStep = chartWidth / Math.max(1, groups.length);
  const barWidth = Math.min(48, Math.max(18, xStep * 0.38));
  const y = (value: number) => margin.top + chartHeight - (value / yMax) * chartHeight;

  return (
    <div className="expression-plot-wrap">
      <svg className="expression-plot" style={{ aspectRatio: `${width} / ${height}` }} viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${yLabel} expression plot`}>
        <text
          x={-(margin.top + chartHeight / 2)}
          y="22"
          className="plot-y-title"
          textAnchor="middle"
          transform="rotate(-90)"
        >
          {yLabel}
        </text>
        <line x1={margin.left} x2={margin.left} y1={margin.top} y2={margin.top + chartHeight} className="plot-axis" />
        <line x1={margin.left} x2={margin.left + chartWidth} y1={margin.top + chartHeight} y2={margin.top + chartHeight} className="plot-axis" />
        {ticks.map((tick) => (
          <g key={tick}>
            <line x1={margin.left - 6} x2={margin.left} y1={y(tick)} y2={y(tick)} className="plot-axis-tick" />
            <line x1={margin.left} x2={margin.left + chartWidth} y1={y(tick)} y2={y(tick)} className={tick === 0 ? 'plot-tick strong' : 'plot-tick'} />
            <text x={margin.left - 12} y={y(tick) + 4} textAnchor="end" className="plot-tick-label">{formatAxisTick(tick)}</text>
          </g>
        ))}
        {groups.map((group, index) => {
          const center = margin.left + xStep * index + xStep / 2;
          const barY = y(group.mean);
          const barHeight = margin.top + chartHeight - barY;
          return (
            <g key={group.key}>
              <rect x={center - barWidth / 2} y={barY} width={barWidth} height={barHeight} rx="2" className="plot-bar" />
              {group.values.length > 1 && group.values.map((value, pointIndex) => {
                const offset = ((pointIndex - (group.values.length - 1) / 2) / Math.max(1, group.values.length - 1)) * Math.min(barWidth * 0.8, 42);
                return <circle key={`${group.key}-${pointIndex}`} cx={center + offset} cy={y(value)} r="4.4" className="plot-point" />;
              })}
              <text x={center} y={margin.top + chartHeight + 24} textAnchor="middle" className="plot-x-label">
                {labelLines[index].map((line, lineIndex) => <tspan key={lineIndex} x={center} dy={lineIndex ? 22 : 0}>{line}</tspan>)}
              </text>
            </g>
          );
        })}
      </svg>
    </div>
  );
}

function niceStep(value: number) {
  if (!Number.isFinite(value) || value <= 0) return 1;
  const exponent = Math.floor(Math.log10(value));
  const base = 10 ** exponent;
  const fraction = value / base;
  const niceFraction = fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10;
  return niceFraction * base;
}

function expressionAxisTicks(maxValue: number) {
  if (!Number.isFinite(maxValue) || maxValue <= 0) {
    return { yMax: 1, ticks: [0, 0.25, 0.5, 0.75, 1] };
  }
  const step = niceStep((maxValue * 1.08) / 4);
  const yMax = step * Math.ceil((maxValue * 1.08) / step);
  return { yMax, ticks: [0, step, step * 2, step * 3, step * 4].filter((tick) => tick <= yMax + step * 0.001) };
}

function formatAxisTick(value: number) {
  if (!Number.isFinite(value)) return '-';
  if (Number.isInteger(value)) return value.toLocaleString();
  const rounded = Math.round(value * 10) / 10;
  if (Number.isInteger(rounded)) return rounded.toLocaleString();
  return rounded.toLocaleString(undefined, { minimumFractionDigits: 1, maximumFractionDigits: 1 });
}

function SpeciesInlineControl({
  speciesId,
  speciesOptions,
  onSelectSpecies,
}: {
  speciesId: string;
  speciesOptions: SpeciesSummary[];
  onSelectSpecies: (speciesId: string) => void;
}) {
  return (
    <div className="inline-control-row single">
      <label>
        Species
        <select value={speciesId} onChange={(event) => onSelectSpecies(event.target.value)}>
          {speciesOptions.map((item) => (
            <option key={item.species_id} value={item.species_id}>
              {item.display_name} · {item.scientific_name}
            </option>
          ))}
        </select>
      </label>
    </div>
  );
}

function GeneControlPanel({
  geneInput,
  onGeneInput,
  speciesId,
  speciesOptions,
  onSelectSpecies,
}: {
  geneInput: string;
  onGeneInput: (value: string) => void;
  speciesId: string;
  speciesOptions: SpeciesSummary[];
  onSelectSpecies: (speciesId: string) => void;
}) {
  const [exampleGene] = exampleGenesForSpecies(speciesId);

  return (
    <div className="inline-control-row">
      <label>
        Species
        <select value={speciesId} onChange={(event) => onSelectSpecies(event.target.value)}>
          {speciesOptions.map((item) => (
            <option key={item.species_id} value={item.species_id}>
              {item.display_name} · {item.scientific_name}
            </option>
          ))}
        </select>
      </label>
      <label>
        Gene ID
        <input value={geneInput} onChange={(event) => onGeneInput(event.target.value)} placeholder={exampleGene} />
      </label>
    </div>
  );
}

function PrimerDesignPanel({
  geneInput,
  speciesId,
  speciesOptions,
  flank,
  genomicRegion,
  tssUpstream,
  tssDownstream,
  productRange,
  primerLengthMin,
  primerLengthMax,
  primerTmMin,
  primerTmMax,
  primerGcMin,
  primerGcMax,
  source,
  templateType,
  sequenceInput,
  sequenceError,
  junctionMode,
  stopCodonMode,
  forwardTail,
  reverseTail,
  checkSpecificity,
  isRunning,
  onGeneInput,
  onSelectSpecies,
  onFlank,
  onGenomicRegion,
  onTssUpstream,
  onTssDownstream,
  onProductRange,
  onPrimerLengthMin,
  onPrimerLengthMax,
  onPrimerTmMin,
  onPrimerTmMax,
  onPrimerGcMin,
  onPrimerGcMax,
  onSource,
  onTemplateType,
  onSequenceInput,
  onJunctionMode,
  onStopCodonMode,
  onForwardTail,
  onReverseTail,
  onCheckSpecificity,
  onExample,
  onRun,
}: {
  geneInput: string;
  speciesId: string;
  speciesOptions: SpeciesSummary[];
  flank: number;
  genomicRegion: 'gene_body' | 'tss_window';
  tssUpstream: number;
  tssDownstream: number;
  productRange: string;
  primerLengthMin: number;
  primerLengthMax: number;
  primerTmMin: number;
  primerTmMax: number;
  primerGcMin: number;
  primerGcMax: number;
  source: 'gene' | 'sequence';
  templateType: 'gdna' | 'cdna' | 'cdna_cloning';
  sequenceInput: string;
  sequenceError: string;
  junctionMode: 'no_junction' | 'allow_junction';
  stopCodonMode: 'keep' | 'remove';
  forwardTail: string;
  reverseTail: string;
  checkSpecificity: boolean;
  isRunning: boolean;
  onGeneInput: (value: string) => void;
  onSelectSpecies: (speciesId: string) => void;
  onFlank: (value: number) => void;
  onGenomicRegion: (value: 'gene_body' | 'tss_window') => void;
  onTssUpstream: (value: number) => void;
  onTssDownstream: (value: number) => void;
  onProductRange: (value: string) => void;
  onPrimerLengthMin: (value: number) => void;
  onPrimerLengthMax: (value: number) => void;
  onPrimerTmMin: (value: number) => void;
  onPrimerTmMax: (value: number) => void;
  onPrimerGcMin: (value: number) => void;
  onPrimerGcMax: (value: number) => void;
  onSource: (value: 'gene' | 'sequence') => void;
  onTemplateType: (value: 'gdna' | 'cdna' | 'cdna_cloning') => void;
  onSequenceInput: (value: string) => void;
  onJunctionMode: (value: 'no_junction' | 'allow_junction') => void;
  onStopCodonMode: (value: 'keep' | 'remove') => void;
  onForwardTail: (value: string) => void;
  onReverseTail: (value: string) => void;
  onCheckSpecificity: (value: boolean) => void;
  onExample: () => void;
  onRun: (transcriptId?: string) => void;
}) {
  const selectedSpecies = speciesOptions.find((item) => item.species_id === speciesId);
  const genePlaceholder = selectedSpecies
    ? `Enter a ${selectedSpecies.display_name} gene ID, then select Resolve`
    : 'Enter a gene ID, then select Resolve';
  const [resolvedGene, setResolvedGene] = useState<ApiRecord | null>(null);
  const [resolvedGeneId, setResolvedGeneId] = useState('');
  const [selectedTranscriptId, setSelectedTranscriptId] = useState('');
  const [isResolving, setIsResolving] = useState(false);
  const [resolveError, setResolveError] = useState('');
  const [isTranscriptMenuOpen, setIsTranscriptMenuOpen] = useState(false);
  const [openParameterHelp, setOpenParameterHelp] = useState<string | null>(null);
  const transcriptMenuRef = useRef<HTMLDivElement | null>(null);
  const parameterOptionsRef = useRef<HTMLElement | null>(null);
  const resolveGeneration = useRef(0);
  const transcriptOptions = Array.isArray(resolvedGene?.transcripts) ? [...resolvedGene.transcripts].sort((a, b) => String(a.transcript_id).localeCompare(String(b.transcript_id), undefined, { numeric: true })) : [];
  const selectedTranscript = transcriptOptions.find((item) => item.transcript_id === selectedTranscriptId) ?? transcriptOptions[0];
  const [productMinValue, productMaxValue] = productRange.split('-').map((value) => Number(value));
  const productMin = Number.isFinite(productMinValue) ? productMinValue : 100;
  const productMax = Number.isFinite(productMaxValue) ? productMaxValue : 300;
  const setProductBounds = (min: number, max: number) => {
    const nextMin = Math.max(50, Math.round(min) || 50);
    const nextMax = Math.max(nextMin + 1, Math.round(max) || nextMin + 1);
    onProductRange(`${nextMin}-${nextMax}`);
  };
  const designAvailable = source === 'sequence'
    ? Boolean(sequenceInput.trim())
    : Boolean(!isResolving && resolvedGene && resolvedGeneId === geneInput.trim() && (templateType === 'gdna' || selectedTranscriptId));
  const isCloningTemplate = source === 'gene' && templateType === 'cdna_cloning';

  useEffect(() => {
    const closeTranscriptMenu = (event: MouseEvent) => {
      if (transcriptMenuRef.current && !transcriptMenuRef.current.contains(event.target as Node)) {
        setIsTranscriptMenuOpen(false);
      }
    };
    document.addEventListener('mousedown', closeTranscriptMenu);
    return () => document.removeEventListener('mousedown', closeTranscriptMenu);
  }, []);

  useEffect(() => {
    if (!openParameterHelp) return;
    const closeParameterHelp = (event: MouseEvent) => {
      if (!parameterOptionsRef.current?.contains(event.target as Node)) setOpenParameterHelp(null);
    };
    document.addEventListener('mousedown', closeParameterHelp);
    return () => document.removeEventListener('mousedown', closeParameterHelp);
  }, [openParameterHelp]);

  useEffect(() => {
    resolveGeneration.current += 1;
    setIsResolving(false);
    setResolvedGene(null);
    setResolvedGeneId('');
    setSelectedTranscriptId('');
    setResolveError('');
    setIsTranscriptMenuOpen(false);
  }, [speciesId, geneInput]);

  const resolveGene = async () => {
    const value = geneInput.trim();
    if (!value) {
      setResolveError('Enter a gene ID first.');
      return;
    }
    setIsResolving(true);
    setResolveError('');
    const generation = ++resolveGeneration.current;
    try {
      const payload = await apiGet<ApiRecord>(`/api/species/${speciesId}/genes/${encodeURIComponent(value)}`);
      if (generation !== resolveGeneration.current) return;
      const transcripts = Array.isArray(payload.transcripts) ? [...payload.transcripts].sort((a, b) => String(a.transcript_id).localeCompare(String(b.transcript_id), undefined, { numeric: true })) : [];
      setResolvedGene(payload);
      setResolvedGeneId(value);
      setSelectedTranscriptId(String(transcripts[0]?.transcript_id ?? ''));
      setIsTranscriptMenuOpen(false);
    } catch (error) {
      if (generation !== resolveGeneration.current) return;
      setResolvedGene(null);
      setResolvedGeneId('');
      setSelectedTranscriptId('');
      const message = error instanceof Error ? error.message : '';
      setResolveError(message.includes('404') ? 'No matching gene ID was found for the selected species.' : (message || 'Gene could not be resolved.'));
    } finally {
      if (generation === resolveGeneration.current) setIsResolving(false);
    }
  };

  const transcriptSummary = (transcript: ApiRecord) => {
    const featureType = templateType === 'cdna_cloning' ? 'CDS' : 'exon';
    const length = annotatedTranscriptLength(resolvedGene?.features ?? [], String(transcript.transcript_id), featureType);
    return `${transcript.transcript_id} · ${length === null ? 'Length unavailable' : `${length.toLocaleString()} nt ${featureType === 'CDS' ? 'CDS' : 'cDNA'}`}`;
  };

  const resolvedGeneSummary = resolvedGene && resolvedGeneId === geneInput.trim() ? (
    <div className="primer-resolved-gene"><span>Resolved gene</span><strong>{resolvedGene.gene?.gene_id ?? resolvedGeneId}</strong></div>
  ) : null;

  const transcriptSelector = resolvedGene && resolvedGeneId === geneInput.trim() ? (
    <section className="primer-transcript-resolution" aria-label="Resolved transcript">
      <div className="primer-transcript-field" ref={transcriptMenuRef}>
        <span>Transcript</span>
        <button type="button" className="primer-transcript-trigger" aria-haspopup="listbox" aria-expanded={isTranscriptMenuOpen} onClick={() => setIsTranscriptMenuOpen((value) => !value)}>
          <strong>{selectedTranscript ? transcriptSummary(selectedTranscript) : 'No transcript model available'}</strong>
          <ChevronDown size={17} aria-hidden="true" />
        </button>
        {isTranscriptMenuOpen && transcriptOptions.length > 0 ? (
          <div className="primer-transcript-menu" role="listbox" aria-label="Transcript">
            {transcriptOptions.map((transcript) => (
              <button type="button" key={String(transcript.transcript_id)} role="option" aria-selected={transcript.transcript_id === selectedTranscript?.transcript_id} className={transcript.transcript_id === selectedTranscript?.transcript_id ? 'active' : ''} onClick={() => {
                setSelectedTranscriptId(String(transcript.transcript_id));
                setIsTranscriptMenuOpen(false);
              }}>
                {transcriptSummary(transcript)}
              </button>
            ))}
          </div>
        ) : null}
      </div>
    </section>
  ) : (
    <p className="primer-transcript-pending">Resolve the gene ID to select an annotated transcript.</p>
  );

  return (
    <section className="primer-design-workspace" aria-label="PCR primer design input">
      <div className="primer-grid">
        <section className="primer-input-panel">
          <nav className="input-mode-tabs" aria-label="PCR primer input mode">
            <button type="button" className={source === 'gene' ? 'active' : ''} onClick={() => onSource('gene')}><Search size={16} aria-hidden="true" />Target gene</button>
            <button type="button" className={source === 'sequence' ? 'active' : ''} onClick={() => onSource('sequence')}><ClipboardPaste size={16} aria-hidden="true" />Paste sequence</button>
          </nav>
          {source === 'gene' ? (
            <>
              <CustomSpeciesSelect speciesId={speciesId} speciesOptions={speciesOptions} onSelectSpecies={onSelectSpecies} />
              <label className="field-label">
                <span className="primer-gene-label-title">Gene ID</span>
                <span className="primer-gene-resolve">
                  <input value={geneInput} onChange={(event) => {
                    onGeneInput(event.target.value);
                    if (event.target.value.trim() !== resolvedGeneId) {
                      setResolvedGene(null);
                      setSelectedTranscriptId('');
                      setResolveError('');
                    }
                  }} placeholder={genePlaceholder} spellCheck={false} onKeyDown={(event) => {
                    if (event.key === 'Enter') {
                      event.preventDefault();
                      void resolveGene();
                    }
                  }} />
                  <button type="button" onClick={() => void resolveGene()} disabled={isResolving} aria-label="Resolve gene" title="Resolve gene">
                    {isResolving ? <LoaderCircle size={16} className="semantic-spinner" aria-hidden="true" /> : <Search size={16} aria-hidden="true" />}
                  </button>
                </span>
              </label>
              {resolveError ? <p className="primer-resolve-error">{resolveError}</p> : null}
              <section className="primer-template-panel" aria-label="Gene-based primer template">
                <div className="segmented compact primer-template-tabs" role="group" aria-label="Primer template type">
                  <label className={templateType === 'gdna' ? 'active' : ''}><input type="radio" name="polarPrimerTemplate" checked={templateType === 'gdna'} onChange={() => onTemplateType('gdna')} /> Genomic DNA</label>
                  <label className={templateType === 'cdna' ? 'active' : ''}><input type="radio" name="polarPrimerTemplate" checked={templateType === 'cdna'} onChange={() => onTemplateType('cdna')} /> cDNA</label>
                  <label className={templateType === 'cdna_cloning' ? 'active' : ''}><input type="radio" name="polarPrimerTemplate" checked={templateType === 'cdna_cloning'} onChange={() => onTemplateType('cdna_cloning')} /> cDNA cloning</label>
                </div>
                <div className="primer-template-heading">
                  <strong>Template for this gene</strong>
                  <span>Choose the sequence context for primer design.</span>
                </div>
                {templateType === 'gdna' ? resolvedGeneSummary : null}
                {templateType === 'gdna' ? (
                  <div className="primer-context-options">
                    <div className="mini-section-title">Genomic region</div>
                    <div className="primer-region-options" role="group" aria-label="Genomic primer region">
                      <label className={genomicRegion === 'gene_body' ? 'active' : ''}>
                        <input type="radio" name="polarPrimerRegion" checked={genomicRegion === 'gene_body'} onChange={() => onGenomicRegion('gene_body')} />
                        <span className="primer-checkmark" aria-hidden="true" />
                        <span><strong>Gene body</strong><small>Design within the annotated gene interval.</small></span>
                      </label>
                      <label className={genomicRegion === 'tss_window' ? 'active' : ''}>
                        <input type="radio" name="polarPrimerRegion" checked={genomicRegion === 'tss_window'} onChange={() => onGenomicRegion('tss_window')} />
                        <span className="primer-checkmark" aria-hidden="true" />
                        <span><strong>TSS flanking window</strong><small>Use a strand-aware window around the transcription start site.</small></span>
                      </label>
                    </div>
                    {genomicRegion === 'tss_window' ? (
                      <div className="primer-flank-row" aria-label="TSS flanking window">
                        <label className="primer-flank-control"><span>Upstream</span>
                          <PrimerStepperControl label="Upstream base pairs" value={tssUpstream} min={0} max={5000} step={50} unit="bp" onChange={onTssUpstream} />
                        </label>
                        <label className="primer-flank-control"><span>Downstream</span>
                          <PrimerStepperControl label="Downstream base pairs" value={tssDownstream} min={0} max={5000} step={50} unit="bp" onChange={onTssDownstream} />
                        </label>
                      </div>
                    ) : null}
                  </div>
                ) : null}
                {templateType === 'cdna' ? (
                  <div className="primer-context-options">
                    {resolvedGeneSummary}
                    {transcriptSelector}
                    <div className="mini-section-title">Junction handling</div>
                    <p className="form-hint">Select the preferred primer-placement rule for this transcript.</p>
                    <div className="primer-junction-options" role="group" aria-label="cDNA junction handling">
                      <label className={junctionMode === 'no_junction' ? 'active' : ''}>
                        <input type="radio" name="polarPrimerJunction" checked={junctionMode === 'no_junction'} onChange={() => onJunctionMode('no_junction')} />
                        <span className="primer-checkmark" aria-hidden="true" />
                        <span className="primer-junction-copy"><strong>Exon only</strong><small>The entire amplicon stays within a single exon.</small></span>
                        <span className="cdna-inline-schematic" aria-hidden="true"><span className="cdna-exon" /><span className="cdna-gap" /><span className="cdna-exon" /><span className="cdna-primer exon-primer" /></span>
                      </label>
                      <label className={junctionMode === 'allow_junction' ? 'active' : ''}>
                        <input type="radio" name="polarPrimerJunction" checked={junctionMode === 'allow_junction'} onChange={() => onJunctionMode('allow_junction')} />
                        <span className="primer-checkmark" aria-hidden="true" />
                        <span className="primer-junction-copy"><strong>Multiple exons only</strong><small>The amplicon must span at least two exons. Primers need not cross a junction.</small></span>
                        <span className="cdna-inline-schematic" aria-hidden="true"><span className="cdna-exon" /><span className="cdna-gap" /><span className="cdna-exon" /><span className="cdna-primer multi-exon-product" /></span>
                      </label>
                    </div>
                  </div>
                ) : null}
                {templateType === 'cdna_cloning' ? (
                  <div className="primer-context-options">
                    {resolvedGeneSummary}
                    {transcriptSelector}
                    <div className="mini-section-title">Cloning options</div>
                    <div className="segmented compact primer-cloning-tabs" role="group" aria-label="Stop codon handling">
                      <label className={stopCodonMode === 'keep' ? 'active' : ''}><input type="radio" name="polarPrimerStop" checked={stopCodonMode === 'keep'} onChange={() => onStopCodonMode('keep')} /> Keep stop codon</label>
                      <label className={stopCodonMode === 'remove' ? 'active' : ''}><input type="radio" name="polarPrimerStop" checked={stopCodonMode === 'remove'} onChange={() => onStopCodonMode('remove')} /> Remove stop codon</label>
                    </div>
                    <div className="parameter-row primer-cloning-tails">
                      <label className="field-label"><span className="primer-cloning-label">Forward 5' addition</span>
                        <input value={forwardTail} onChange={(event) => onForwardTail(event.target.value)} placeholder="Optional tag or site" spellCheck={false} />
                      </label>
                      <label className="field-label"><span className="primer-cloning-label">Reverse 3' addition (<strong className="primer-sense-strand">sense strand</strong>)</span>
                        <input value={reverseTail} onChange={(event) => onReverseTail(event.target.value)} placeholder="Optional tag or site" spellCheck={false} />
                      </label>
                    </div>
                    <p className="form-hint">Full-length cloning primers are anchored to the selected CDS boundaries. Optional 5' additions are excluded from annealing calculations.</p>
                  </div>
                ) : null}
              </section>
            </>
          ) : (
            <>
              <CustomSpeciesSelect speciesId={speciesId} speciesOptions={speciesOptions} onSelectSpecies={onSelectSpecies} />
              <label className="field-label primer-sequence-field">DNA sequence
                <textarea value={sequenceInput} rows={10} spellCheck={false} placeholder="Paste FASTA or plain nucleotide sequence" onChange={(event) => onSequenceInput(event.target.value)} />
                {sequenceError ? <p className="primer-sequence-error">{sequenceError}</p> : null}
              </label>
            </>
          )}
        </section>

        <section ref={parameterOptionsRef} className="primer-options-panel" aria-label="Primer3 design options">
          <div className="parameter-row" aria-label={isCloningTemplate ? 'Amplicon size is fixed to the full-length CDS' : undefined}>
            <PrimerNumberStepper label="Amplicon min" value={productMin} min={50} max={5000} step={10} unit="bp" onChange={(value) => setProductBounds(value, productMax)} disabled={isCloningTemplate} openHelp={openParameterHelp} onToggleHelp={setOpenParameterHelp} />
            <PrimerNumberStepper label="Amplicon max" value={productMax} min={50} max={5000} step={10} unit="bp" onChange={(value) => setProductBounds(productMin, value)} disabled={isCloningTemplate} openHelp={openParameterHelp} onToggleHelp={setOpenParameterHelp} />
          </div>
          <div className="parameter-row">
            <PrimerNumberStepper label="Primer min" value={primerLengthMin} min={15} max={source === 'gene' && templateType === 'cdna_cloning' ? 40 : 35} step={1} unit="nt" onChange={onPrimerLengthMin} openHelp={openParameterHelp} onToggleHelp={setOpenParameterHelp} />
            <PrimerNumberStepper label="Primer max" value={primerLengthMax} min={15} max={source === 'gene' && templateType === 'cdna_cloning' ? 40 : 35} step={1} unit="nt" onChange={onPrimerLengthMax} openHelp={openParameterHelp} onToggleHelp={setOpenParameterHelp} />
          </div>
          <div className="parameter-row">
            <PrimerNumberStepper label="Tm min" value={primerTmMin} min={40} max={90} step={0.5} unit="°C" onChange={onPrimerTmMin} openHelp={openParameterHelp} onToggleHelp={setOpenParameterHelp} />
            <PrimerNumberStepper label="Tm max" value={primerTmMax} min={40} max={90} step={0.5} unit="°C" onChange={onPrimerTmMax} openHelp={openParameterHelp} onToggleHelp={setOpenParameterHelp} />
          </div>
          <div className="parameter-row">
            <PrimerNumberStepper label="GC min" value={primerGcMin} min={0} max={100} step={1} unit="%" onChange={onPrimerGcMin} openHelp={openParameterHelp} onToggleHelp={setOpenParameterHelp} />
            <PrimerNumberStepper label="GC max" value={primerGcMax} min={0} max={100} step={1} unit="%" onChange={onPrimerGcMax} openHelp={openParameterHelp} onToggleHelp={setOpenParameterHelp} />
          </div>
          {isCloningTemplate ? (
            <div className="primer-specificity-unavailable" aria-label="Genome specificity is not used for full-length cloning">
              <span>Genome specificity</span>
              <strong>Not used for full-length cloning</strong>
            </div>
          ) : (
            <label className="primer-specificity-toggle">
              <input type="checkbox" checked={checkSpecificity} onChange={(event) => onCheckSpecificity(event.target.checked)} />
              <span>
                <strong>Genome specificity</strong>
                <small>Screen exact primer matches and, where applicable, candidate pair products across the selected genome.</small>
              </span>
            </label>
          )}
          <div className="primer-design-actions">
            <button type="button" className="gene-search-submit" onClick={() => onRun(selectedTranscriptId || undefined)} disabled={isRunning || !designAvailable} title={designAvailable ? 'Design PCR primers' : source === 'gene' ? 'Resolve the selected gene before designing primers' : 'Paste a DNA sequence before designing primers'}>
              {isRunning ? <LoaderCircle size={17} className="semantic-spinner" aria-hidden="true" /> : null}
              {isRunning ? 'Designing' : 'Design primers'}
            </button>
            <button type="button" className="secondary gene-example-button" onClick={onExample} disabled={isRunning || source === 'sequence'} title={source === 'sequence' ? 'Paste a DNA sequence to design primers' : 'Load a primer-design example'}>Example</button>
          </div>
        </section>
      </div>
    </section>
  );
}

function PrimerStepperControl({ label, value, min, max, step, unit, onChange, disabled = false }: {
  label: string; value: number; min: number; max: number; step: number;
  unit: string; onChange: (value: number) => void; disabled?: boolean;
}) {
  const [draft, setDraft] = useState<string | null>(null);
  useEffect(() => { setDraft(null); }, [value, min, max, disabled]);
  const precision = step < 1 ? 1 : 0;
  const draftValue = () => {
    const text = (draft ?? String(value)).trim();
    const parsed = text === '' ? NaN : Number(text);
    return Number.isFinite(parsed) ? parsed : value;
  };
  const commit = (candidate: number) => {
    const next = Number(Math.min(max, Math.max(min, candidate)).toFixed(precision));
    setDraft(null);
    onChange(next);
  };
  return (
    <span className="primer-stepper">
      <button type="button" disabled={disabled} onMouseDown={(event) => event.preventDefault()} onClick={() => commit(draftValue() - step)} aria-label={`Decrease ${label}`}><Minus size={14} aria-hidden="true" /></button>
      <input type="text" role="spinbutton" inputMode={step < 1 ? 'decimal' : 'numeric'} disabled={disabled}
        aria-valuemin={min} aria-valuemax={max} aria-valuenow={value} aria-label={label}
        value={draft ?? String(value)} onChange={(event) => setDraft(event.target.value)}
        onBlur={() => commit(draftValue())} onKeyDown={(event) => {
          if (event.key === 'Enter') {
            event.preventDefault();
            event.currentTarget.blur();
          } else if (event.key === 'Escape') {
            event.preventDefault();
            setDraft(null);
          } else if (event.key === 'ArrowUp' || event.key === 'ArrowDown') {
            event.preventDefault();
            commit(draftValue() + (event.key === 'ArrowUp' ? step : -step));
          }
        }} />
      <button type="button" disabled={disabled} onMouseDown={(event) => event.preventDefault()} onClick={() => commit(draftValue() + step)} aria-label={`Increase ${label}`}><Plus size={14} aria-hidden="true" /></button>
      <em>{unit}</em>
    </span>
  );
}

function PrimerNumberStepper({
  label,
  value,
  min,
  max,
  step,
  unit,
  onChange,
  disabled = false,
  openHelp,
  onToggleHelp,
}: {
  label: string;
  value: number;
  min: number;
  max: number;
  step: number;
  unit: string;
  onChange: (value: number) => void;
  disabled?: boolean;
  openHelp: string | null;
  onToggleHelp: (id: string | null) => void;
}) {
  return (
    <div className={`field-label primer-parameter-control ${disabled ? 'is-disabled' : ''}`} title={disabled ? 'Fixed to full-length CDS for cloning' : undefined}>
      <PrimerParameterLabel label={label} openHelp={openHelp} onToggleHelp={onToggleHelp} />
      <PrimerStepperControl label={label} value={value} min={min} max={max} step={step} unit={unit} onChange={onChange} disabled={disabled} />
    </div>
  );
}

const primerParameterHelp: Record<string, string> = {
  'Amplicon min': 'Shortest acceptable product. Range: 50-5,000 bp.',
  'Amplicon max': 'Longest acceptable product. Range: 50-5,000 bp.',
  'Primer min': 'Shortest annealing region, excluding additions. Range: 15-35 nt; cloning up to 40 nt.',
  'Primer max': 'Longest annealing region, excluding additions. Range: 15-35 nt; cloning up to 40 nt.',
  'Tm min': 'Lowest primer melting temperature. Range: 40-90 C.',
  'Tm max': 'Highest primer melting temperature. Range: 40-90 C.',
  'GC min': 'Lowest primer GC content. Range: 0-100%.',
  'GC max': 'Highest primer GC content. Range: 0-100%.',
};

function PrimerParameterLabel({ label, openHelp, onToggleHelp }: { label: string; openHelp: string | null; onToggleHelp: (id: string | null) => void }) {
  const id = `primer-${label.toLowerCase().replace(/\s+/g, '-')}`;
  const isOpen = openHelp === id;
  return (
    <span className="primer-parameter-label">
      <span>{label}</span>
      <button type="button" aria-label={`About ${label}`} aria-expanded={isOpen} onClick={() => onToggleHelp(isOpen ? null : id)}><CircleHelp size={14} aria-hidden="true" /></button>
      {isOpen ? <span className="primer-parameter-tooltip" role="tooltip">{primerParameterHelp[label]}</span> : null}
    </span>
  );
}

function CustomSpeciesSelect({
  speciesId,
  speciesOptions,
  onSelectSpecies,
}: {
  speciesId: string;
  speciesOptions: SpeciesSummary[];
  onSelectSpecies: (speciesId: string) => void;
}) {
  const [isSpeciesMenuOpen, setIsSpeciesMenuOpen] = useState(false);
  const speciesMenuRef = useRef<HTMLDivElement | null>(null);
  const selectedSpecies = speciesOptions.find((item) => item.species_id === speciesId);

  useEffect(() => {
    const closeSpeciesMenu = (event: MouseEvent) => {
      if (speciesMenuRef.current && !speciesMenuRef.current.contains(event.target as Node)) {
        setIsSpeciesMenuOpen(false);
      }
    };
    document.addEventListener('mousedown', closeSpeciesMenu);
    return () => document.removeEventListener('mousedown', closeSpeciesMenu);
  }, []);

  return (
    <div className="genome-species-field">
      <span className="species-field-label">Species</span>
      <div className="genome-select-shell" ref={speciesMenuRef}>
        <button
          type="button"
          className="genome-select-trigger"
          aria-haspopup="listbox"
          aria-label="Select species"
          aria-expanded={isSpeciesMenuOpen}
          onClick={() => setIsSpeciesMenuOpen((value) => !value)}
        >
          <span>
            <strong>{selectedSpecies?.display_name ?? speciesId}</strong>
            <small>{selectedSpecies?.scientific_name ?? 'Select species'}</small>
          </span>
          <ChevronDown size={18} aria-hidden="true" />
        </button>
        {isSpeciesMenuOpen && (
          <div className="genome-select-menu" role="listbox" aria-label="Species">
            {speciesOptions.map((item) => (
              <button
                type="button"
                key={item.species_id}
                className={item.species_id === speciesId ? 'active' : ''}
                role="option"
                aria-selected={item.species_id === speciesId}
                onClick={() => {
                  onSelectSpecies(item.species_id);
                  setIsSpeciesMenuOpen(false);
                }}
              >
                <strong>{item.display_name}</strong>
                <small>{item.scientific_name}</small>
              </button>
            ))}
          </div>
        )}
        <select value={speciesId} onChange={(event) => onSelectSpecies(event.target.value)} aria-hidden="true" tabIndex={-1}>
          {speciesOptions.map((item) => (
            <option key={item.species_id} value={item.species_id}>
              {item.display_name} · {item.scientific_name}
            </option>
          ))}
        </select>
      </div>
    </div>
  );
}

function GenomeControlPanel({
  geneInput,
  onGeneInput,
  onRun,
  speciesId,
  speciesOptions,
  onSelectSpecies,
}: {
  geneInput: string;
  onGeneInput: (value: string) => void;
  onRun: () => void;
  speciesId: string;
  speciesOptions: SpeciesSummary[];
  onSelectSpecies: (speciesId: string) => void;
}) {
  const [exampleGene] = exampleGenesForSpecies(speciesId);

  return (
    <div className="genome-control-panel">
      <CustomSpeciesSelect speciesId={speciesId} speciesOptions={speciesOptions} onSelectSpecies={onSelectSpecies} />
      <div className="genome-locus-field">
        <label htmlFor="genome-locus-input">Gene ID or locus</label>
        <div className="genome-locus-control">
          <input id="genome-locus-input" value={geneInput} onChange={(event) => onGeneInput(event.target.value)} placeholder={exampleGene} onKeyDown={(event) => { if (event.key === 'Enter' && !event.nativeEvent.isComposing) onRun(); }} />
          <button type="button" className="genome-locus-load" onClick={onRun}>Load</button>
        </div>
      </div>
    </div>
  );
}

function EmbeddedGenomeBrowser({ config, location }: { config: ApiRecord | undefined; location: string | undefined }) {
  const assembly = config?.assemblies?.[0];
  const rawTracks = Array.isArray(config?.tracks) ? config.tracks : [];
  const tracks = useMemo(() => prepareJBrowseTracks(rawTracks), [rawTracks]);
  const jbrowseLocation = location ? normalizeJBrowseLocation(expandJBrowseLocation(location, 2)) : undefined;
  const defaultSession = useMemo(() => {
    if (!assembly || !tracks.length) return undefined;
    return buildJBrowseSession(assembly, tracks);
  }, [assembly, tracks]);
  const viewState = useMemo(() => {
    if (!assembly || !tracks.length || !defaultSession) return null;
    return createViewState({
      assembly,
      tracks,
      location: jbrowseLocation,
      defaultSession,
      disableAddTracks: false,
    });
  }, [assembly, tracks, jbrowseLocation, defaultSession]);

  if (!viewState) {
    return <EmptyResult text="JBrowse configuration is not available for this species." />;
  }

  return (
    <div className="embedded-browser">
      <JBrowseLinearGenomeView viewState={viewState} />
    </div>
  );
}

function normalizeJBrowseLocation(location: string) {
  const match = location.match(/^([^:]+):([\d,]+)-([\d,]+)$/);
  if (!match) return location;
  return `${match[1]}:${match[2]}..${match[3]}`;
}

function expandJBrowseLocation(location: string, fold: number) {
  const match = location.match(/^([^:]+):([\d,]+)-([\d,]+)$/);
  if (!match) return location;
  const start = Number(match[2].replace(/,/g, ''));
  const end = Number(match[3].replace(/,/g, ''));
  if (!Number.isFinite(start) || !Number.isFinite(end) || end <= start) return location;

  const width = end - start + 1;
  const targetWidth = Math.round(width * fold);
  const center = Math.round((start + end) / 2);
  const expandedStart = Math.max(1, center - Math.floor(targetWidth / 2));
  const expandedEnd = center + Math.ceil(targetWidth / 2);
  return `${match[1]}:${expandedStart}-${expandedEnd}`;
}

const expressionTrackOrder: Record<string, number> = {
  Cq_4: 0,
  Inf_4: 1,
  Cq_16: 2,
  Inf_16: 3,
};

const expressionTrackColors: Record<string, string> = {
  Cq_4: '#2a6f6c',
  Inf_4: '#b6782d',
  Cq_16: '#2f65a7',
  Inf_16: '#8a5aa6',
};

const sampleTrackPalette = [
  '#2a6f6c',
  '#b6782d',
  '#2f65a7',
  '#8a5aa6',
  '#687a2f',
  '#b34f64',
  '#397990',
  '#7c6232',
  '#4f6cb2',
  '#9a5a3f',
];

function getExpressionTrackGroup(track: ApiRecord) {
  const name = String(track.name ?? track.trackId ?? '');
  const match = name.match(/^(Cq|Inf)_(4|16)_/i);
  if (!match) return null;
  const prefix = match[1].toLowerCase() === 'inf' ? 'Inf' : 'Cq';
  return `${prefix}_${match[2]}`;
}

function getExpressionTrackReplicate(track: ApiRecord) {
  const name = String(track.name ?? track.trackId ?? '');
  const match = name.match(/^(?:Cq|Inf)_(?:4|16)_(\d+)/i);
  return match ? Number(match[1]) : 0;
}

function formatBigWigSourceName(track: ApiRecord) {
  return String(track.name ?? track.trackId ?? '')
    .replace(/\s*CPM bin25$/i, '')
    .replace(/^Cq_/, '')
    .replace(/^Da_/, '');
}

function getSpeciesIdFromTrack(track: ApiRecord) {
  return String(track.trackId ?? '').split('_')[0] || 'species';
}

function makeMultiQuantitativeTrack(trackId: string, name: string, sourceTracks: ApiRecord[], colors: string[]) {
  const rowHeight = 48;
  return {
    type: 'MultiQuantitativeTrack',
    trackId,
    name,
    assemblyNames: sourceTracks[0]?.assemblyNames ?? [],
    adapter: {
      type: 'MultiWiggleAdapter',
      subadapters: sourceTracks.map((track, index) => ({
        ...track.adapter,
        source: formatBigWigSourceName(track),
        color: colors[index % colors.length],
      })),
    },
    displays: [
      {
        type: 'MultiLinearWiggleDisplay',
        displayId: `${trackId}-MultiLinearWiggleDisplay`,
        defaultRendering: 'multirowxy',
        autoscale: 'local',
        minimalTicks: true,
        height: Math.max(132, sourceTracks.length * rowHeight),
      },
    ],
  };
}

function prepareJBrowseTracks(tracks: ApiRecord[]) {
  const nonQuantitativeTracks = tracks.filter((track) => track.type !== 'QuantitativeTrack');
  const quantitativeTracks = tracks.filter((track) => track.type === 'QuantitativeTrack');
  const cqGroups = quantitativeTracks.reduce<Record<string, ApiRecord[]>>((acc, track) => {
    const group = getExpressionTrackGroup(track);
    if (!group) return acc;
    acc[group] = [...(acc[group] ?? []), track];
    return acc;
  }, {});

  const groupedCqTracks = Object.entries(cqGroups)
    .sort(([groupA], [groupB]) => (expressionTrackOrder[groupA] ?? 99) - (expressionTrackOrder[groupB] ?? 99))
    .map(([group, groupTracks]) => {
      const sortedTracks = groupTracks.sort((a, b) => getExpressionTrackReplicate(a) - getExpressionTrackReplicate(b));
      const speciesId = getSpeciesIdFromTrack(sortedTracks[0]);
      return makeMultiQuantitativeTrack(
        `${speciesId}_${group}_multiwiggle`,
        `${group} CPM bin25`,
        sortedTracks,
        [expressionTrackColors[group] ?? '#2f65a7'],
      );
    });

  const ungroupedTracks = quantitativeTracks.filter((track) => !getExpressionTrackGroup(track));
  const groupedBySpecies = ungroupedTracks.reduce<Record<string, ApiRecord[]>>((acc, track) => {
    const speciesId = getSpeciesIdFromTrack(track);
    acc[speciesId] = [...(acc[speciesId] ?? []), track];
    return acc;
  }, {});
  const sampleGroupedTracks = Object.entries(groupedBySpecies).map(([speciesId, speciesTracks]) =>
    makeMultiQuantitativeTrack(
      `${speciesId}_expression_multiwiggle`,
      `${speciesId} RNA-seq CPM bin25`,
      speciesTracks,
      sampleTrackPalette,
    ),
  );

  return [...nonQuantitativeTracks, ...groupedCqTracks, ...sampleGroupedTracks];
}

function buildJBrowseSession(assembly: ApiRecord, tracks: ApiRecord[]) {
  const referenceTrackId = assembly.sequence?.trackId;
  const visibleTracks = tracks.filter((track) => ['FeatureTrack', 'QuantitativeTrack', 'MultiQuantitativeTrack'].includes(String(track.type)));
  const prioritizedTracks = [
    ...visibleTracks.filter((track) => track.type === 'FeatureTrack').slice(0, 1),
    ...visibleTracks.filter((track) => track.type === 'QuantitativeTrack' || track.type === 'MultiQuantitativeTrack'),
  ];
  return {
    name: 'Polar Genome DB browser',
    view: {
      id: 'pgdb-linear-genome-view',
      type: 'LinearGenomeView',
      hideHeader: false,
      hideNoTracksActive: true,
      trackSelectorType: 'hierarchical',
      tracks: [
        referenceTrackId
          ? {
              type: 'ReferenceSequenceTrack',
              configuration: referenceTrackId,
              displays: [
                {
                  type: 'LinearReferenceSequenceDisplay',
                  configuration: `${referenceTrackId}-LinearReferenceSequenceDisplay`,
                },
              ],
            }
          : null,
        ...prioritizedTracks.map((track) => ({
          type: track.type,
          configuration: track.trackId,
          displays: [
            {
              type: track.type === 'MultiQuantitativeTrack' ? 'MultiLinearWiggleDisplay' : track.type === 'QuantitativeTrack' ? 'LinearWiggleDisplay' : 'LinearBasicDisplay',
              configuration: `${track.trackId}-${track.type === 'MultiQuantitativeTrack' ? 'MultiLinearWiggleDisplay' : track.type === 'QuantitativeTrack' ? 'LinearWiggleDisplay' : 'LinearBasicDisplay'}`,
              ...(track.type === 'QuantitativeTrack' ? { heightPreConfig: 52 } : {}),
              ...(track.type === 'MultiQuantitativeTrack' ? { heightPreConfig: Math.max(132, ((track.adapter?.subadapters as ApiRecord[] | undefined)?.length ?? 1) * 48) } : {}),
            },
          ],
        })),
      ].filter(Boolean),
    },
  };
}

function EmptyResult({ text, multiline = false }: { text: string; multiline?: boolean }) {
  return <div className={multiline ? 'empty-result primer-empty-result' : 'empty-result'}>{text}</div>;
}

function formatLocation(item: ApiRecord | undefined) {
  if (!item) return '-';
  return `${item.seqid ?? '-'}:${Number(item.start ?? 0).toLocaleString()}-${Number(item.end ?? 0).toLocaleString()}${item.strand ? ` (${item.strand})` : ''}`;
}

function formatSearchLocation(item: ApiRecord) {
  const seqid = String(item.seqid ?? '-').replace(/^.*?((?:chr|chromosome)[_.-]?[a-z0-9]+)$/i, '$1');
  return `${seqid}:${Number(item.start ?? 0).toLocaleString()}-${Number(item.end ?? 0).toLocaleString()}`;
}

function tableColumnLabel(column: string) {
  const labels: Record<string, string> = {
    hmm_range: 'hmm range',
    i_evalue: 'i-value',
    pident: 'Identity (%)',
    evalue: 'E-value',
    reference_transcript_id: 'Reference transcript',
    reference_gene_symbol: 'Gene symbol',
  };
  const text = column.replace(/_/g, ' ');
  return labels[column] ?? text.charAt(0).toUpperCase() + text.slice(1);
}

function CompactTable({
  title,
  titleActions,
  columnLabels,
  columns,
  rows,
  pageSize,
  tableClassName = '',
}: {
  title: string;
  titleActions?: React.ReactNode;
  columnLabels?: Record<string, React.ReactNode>;
  columns: string[];
  rows: ApiRecord[];
  pageSize?: number;
  tableClassName?: string;
}) {
  const [page, setPage] = useState(0);
  const totalPages = pageSize ? Math.ceil(rows.length / pageSize) : 1;
  const currentPage = Math.min(page, Math.max(0, totalPages - 1));
  const visibleRows = pageSize ? rows.slice(currentPage * pageSize, (currentPage + 1) * pageSize) : rows;

  useEffect(() => {
    setPage(0);
  }, [title, rows.length]);

  if (!rows.length) {
    return <EmptyResult text={`${title}: no records returned.`} />;
  }
  return (
    <div className="compact-table-wrap">
      {titleActions ? <div className="compact-table-title-actions"><h4>{title}</h4>{titleActions}</div> : <h4>{title}</h4>}
      <table className={`compact-table ${tableClassName}`.trim()}>
        <thead>
          <tr>
            {columns.map((column) => <th key={column} data-column={column}>{columnLabels?.[column] ?? tableColumnLabel(column)}</th>)}
          </tr>
        </thead>
        <tbody>
          {visibleRows.map((row, index) => (
            <tr key={`${title}-${currentPage * (pageSize ?? rows.length) + index}`}>
              {columns.map((column) => <td key={column} data-column={column}>{formatCell(row[column], column)}</td>)}
            </tr>
          ))}
        </tbody>
      </table>
      {pageSize && totalPages > 1 && (
        <div className="compact-table-pagination">
          <span>{currentPage * pageSize + 1}-{Math.min((currentPage + 1) * pageSize, rows.length)} of {rows.length}</span>
          <div>
            <button
              type="button"
              aria-label="Previous page"
              title="Previous page"
              disabled={currentPage === 0}
              onClick={() => setPage((value) => Math.max(0, value - 1))}
            >
              <ChevronLeft size={16} />
            </button>
            <span>{currentPage + 1} / {totalPages}</span>
            <button
              type="button"
              aria-label="Next page"
              title="Next page"
              disabled={currentPage >= totalPages - 1}
              onClick={() => setPage((value) => Math.min(totalPages - 1, value + 1))}
            >
              <ChevronRight size={16} />
            </button>
          </div>
        </div>
      )}
    </div>
  );
}

function primerSpecificityLabel(value: unknown, specificity?: ApiRecord, exportText = false) {
  const capped = Number(value) >= 2 && (specificity?.hit_count_cap === 2 || String(specificity?.engine ?? '').includes('k2'));
  return primerHitLabel(value, capped, exportText);
}

function primerSpecificityClass(value: unknown, status?: string) {
  if (status === 'not_assessed' || value === null || value === undefined) return 'not-assessed';
  return Number(value) === 1 ? 'unique' : Number(value) === 0 ? 'not-assessed' : 'multiple';
}

function primerPairScreenNote(specificity: ApiRecord) {
  const limit = Number(specificity.screening_max_product_size ?? 1000).toLocaleString();
  const notes: Record<string, string> = {
    target_confirmed: 'Expected genomic binding positions confirmed.',
    target_not_detected: 'Expected genomic product was not detected; specificity is not confirmed.',
    target_out_of_range: `Expected genomic product exceeds the ${limit} bp screening limit.`,
    off_target_detected: 'Products outside the expected binding positions were detected.',
    genome_matches_only: 'Genome products only; the submitted sequence has no assigned target locus.',
    junction_genome_screen: 'Junction-spanning primer: genome contamination screen only, not transcriptome specificity.',
  };
  return `${notes[String(specificity.status)] ?? 'Legacy genome screen; target coordinates were not verified.'} Exact matches only, products up to ${limit} bp. Counts of 2 are lower bounds.`;
}

function primerJunctionLabel(primer: ApiRecord) {
  const junction = primer.junction;
  if (!junction) return '-';
  const left = String(junction.left_context ?? '');
  const right = String(junction.right_context ?? '');
  if (!left && !right) return '-';
  if (left === right) return left;
  return `${left} / ${right.replace(/^Exon\s+/, '')}`;
}

function PrimerAdditionsCell({ primer }: { primer: ApiRecord }) {
  const forward = String(primer.left_tail_sequence ?? '');
  const reverse = String(primer.right_sense_tail_sequence ?? '');
  return (
    <div className="primer-additions-cell">
      <div><span>Forward 5′</span><code>{forward || 'None'}</code></div>
      <div><span>Reverse 3′</span><code>{reverse || 'None'}</code></div>
    </div>
  );
}

function primerResultsTsv(job: ApiRecord, primers: ApiRecord[]) {
  const isCloning = job.results?.target?.genomic_region === 'cdna_cloning' || job.target?.genomic_region === 'cdna_cloning';
  if (isCloning) {
    const headers = ['Rank', 'Forward primer', 'Reverse primer', 'Forward annealing core', 'Reverse annealing core', "Forward 5' addition", "Reverse 3' addition (sense strand)", 'Full-length CDS product (bp)', 'CDS coverage (%)', 'Forward Tm (C)', 'Reverse Tm (C)', 'Forward GC (%)', 'Reverse GC (%)', 'Penalty', 'Final product incl additions (bp)', 'Final oligo structure Tm (JSON)', 'Tm difference (C)', 'Warnings'];
    const rows = primers.map((primer) => [
      primer.rank,
      primer.left_sequence,
      primer.right_sequence,
      primer.left_core_sequence,
      primer.right_core_sequence,
      primer.left_tail_sequence,
      primer.right_sense_tail_sequence,
      primer.product_size,
      Number(primer.full_length_coverage ?? primer.coverage ?? 100).toFixed(1),
      Number(primer.left_tm ?? 0).toFixed(2),
      Number(primer.right_tm ?? 0).toFixed(2),
      Number(primer.left_gc_percent ?? 0).toFixed(1),
      Number(primer.right_gc_percent ?? 0).toFixed(1),
      Number(primer.penalty ?? 0).toFixed(3),
      primer.final_product_size ?? '',
      primer.structure ? JSON.stringify(primer.structure) : '',
      primer.tm_difference ?? '',
      (primer.warnings ?? []).join('; '),
    ]);
    return primerExport(job, headers, rows);
  }
  const headers = ['Rank', 'Forward primer', 'Reverse primer', 'Product size (bp)', 'Forward Tm (C)', 'Reverse Tm (C)', 'Forward GC (%)', 'Reverse GC (%)', 'Forward exact hits', 'Reverse exact hits', 'Pair products', 'Penalty', 'Specificity engine', 'Screening max product (bp)', 'Expected genomic product (bp)', 'Screen status', 'Observed on-target products', 'Observed off-target products', 'Off-target count is a lower bound'];
  const rows = primers.map((primer) => [
    primer.rank,
    primer.left_sequence,
    primer.right_sequence,
    primer.product_size,
    Number(primer.left_tm ?? 0).toFixed(2),
    Number(primer.right_tm ?? 0).toFixed(2),
    Number(primer.left_gc_percent ?? 0).toFixed(1),
    Number(primer.right_gc_percent ?? 0).toFixed(1),
    primerSpecificityLabel(primer.specificity?.left_hit_count, primer.specificity, true),
    primerSpecificityLabel(primer.specificity?.right_hit_count, primer.specificity, true),
    primerSpecificityLabel(primer.specificity?.pair_amplicon_count, primer.specificity, true),
    Number(primer.penalty ?? 0).toFixed(3),
    primer.specificity?.engine ?? '',
    primer.specificity?.screening_max_product_size ?? '',
    primer.expected_genomic_product_size ?? primer.specificity?.expected_genomic_product_size ?? '',
    primer.specificity?.status ?? 'not_assessed',
    primer.specificity?.on_target_product_count ?? 'Not assessed',
    primer.specificity?.off_target_product_count ?? 'Not assessed',
    primer.specificity?.off_target_count_is_lower_bound ?? '',
  ]);
  return primerExport(job, headers, rows);
}

function PrimerProductTrack({ primer, targetLength }: { primer: ApiRecord; targetLength: number }) {
  const leftOffset = Number(primer.left_position?.offset ?? 0);
  const leftLength = Number(primer.left_position?.length ?? 0);
  const rightOffset = Number(primer.right_position?.offset ?? 0);
  const rightLength = Number(primer.right_position?.length ?? 0);
  const rightStart = Math.max(0, rightOffset - rightLength + 1);
  const denominator = Math.max(targetLength, rightOffset + 1, leftOffset + leftLength, 1);
  const productLeft = Math.max(0, Math.min(96, (leftOffset / denominator) * 100));
  const productRight = Math.max(productLeft + 2, Math.min(100, ((rightOffset + 1) / denominator) * 100));
  const leftWidth = Math.max(1.3, (leftLength / denominator) * 100);
  const rightWidth = Math.max(1.3, (rightLength / denominator) * 100);
  return (
    <div className="primer-product-track" aria-label={`${primer.product_size} base pair product`}>
      <span className="primer-track-base" />
      <span className="primer-track-product" style={{ left: `${productLeft}%`, width: `${productRight - productLeft}%` }} />
      <span className="primer-track-forward" style={{ left: `${productLeft}%`, width: `${leftWidth}%` }} />
      <span className="primer-track-reverse" style={{ left: `${Math.max(productLeft, (rightStart / denominator) * 100)}%`, width: `${rightWidth}%` }} />
    </div>
  );
}

function PrimerTable({ job, primers }: { job: ApiRecord; primers: ApiRecord[] }) {
  const [copied, setCopied] = useState('');
  const [page, setPage] = useState(0);
  const [openHeaderHelp, setOpenHeaderHelp] = useState<string | null>(null);
  const resultTableHeadRef = useRef<HTMLTableSectionElement | null>(null);
  const target = job.results?.target ?? job.target ?? {};
  const targetId = String(target.target_id ?? 'Selected target');
  const targetLength = Number(target.length ?? 1);
  const isCdna = target.genomic_region === 'cdna';
  const isCloning = target.genomic_region === 'cdna_cloning';
  const design = job.design ?? {};
  const designContext = isCloning
    ? 'Full-length CDS cloning'
    : isCdna
      ? `cDNA · ${job.results?.logic_version >= 5 ? (design.junction_mode === 'allow_junction' ? 'Multiple exons only' : 'Exon only') : (job.results?.logic_version >= 3 ? (design.junction_mode === 'allow_junction' ? 'Multiple exons preferred' : 'Exon only') : (design.junction_mode === 'allow_junction' ? 'Junction-spanning primer' : 'Exon-contained primers')) + ' (legacy constraints)'}`
      : target.genomic_region === 'submitted_sequence'
        ? 'Submitted DNA sequence'
        : design.genomic_region === 'tss_window' || target.genomic_region === 'tss_window'
          ? `Genomic DNA · TSS flanking (-${Number(design.tss_upstream ?? 1000).toLocaleString()} / +${Number(design.tss_downstream ?? 300).toLocaleString()} bp)`
          : 'Genomic DNA · Gene body';
  const hasSpecificity = primers.some((primer) => primer.specificity);
  const pageSize = isCloning ? 3 : 10;
  const totalPages = Math.max(1, Math.ceil(primers.length / pageSize));
  const currentPage = Math.min(page, totalPages - 1);
  const visiblePrimers = primers.slice(currentPage * pageSize, (currentPage + 1) * pageSize);

  useEffect(() => setPage(0), [job.job_id]);

  useEffect(() => {
    if (!openHeaderHelp) return;
    const closeHeaderHelp = (event: MouseEvent) => {
      if (!resultTableHeadRef.current?.contains(event.target as Node)) setOpenHeaderHelp(null);
    };
    document.addEventListener('mousedown', closeHeaderHelp);
    return () => document.removeEventListener('mousedown', closeHeaderHelp);
  }, [openHeaderHelp]);

  async function copyPrimer(key: string, sequence: string) {
    try {
      await navigator.clipboard.writeText(sequence);
    } catch {
      const area = document.createElement('textarea');
      area.value = sequence;
      area.style.position = 'fixed';
      area.style.left = '-9999px';
      document.body.appendChild(area);
      area.select();
      document.execCommand('copy');
      area.remove();
    }
    setCopied(key);
    window.setTimeout(() => setCopied(''), 900);
  }

  if (!primers.length) {
    return <EmptyResult multiline text={job.status === 'failed' ? (job.error ?? 'Primer design failed.') : job.status === 'completed' ? primerEmptyReason(job.results ?? {}) : 'Primer design is in progress.'} />;
  }
  const exportData = primerResultsTsv(job, primers);
  return (
    <section className="primer-results-panel">
      <header className="primer-results-head">
        <div>
          <h2>PCR Primer Design Results</h2>
          <p><strong>{targetId}</strong>{target.transcript_id ? <span>{target.transcript_id}</span> : null}<span>{isCloning ? `${primers.length} full-length pairs` : `${primers.length} candidates`}</span><span>{designContext}</span><span>{hasSpecificity ? `Exact genome screening · ≤ ${Number(primers[0]?.specificity?.screening_max_product_size ?? 1000).toLocaleString()} bp` : 'Genome specificity not assessed'}</span></p>
        </div>
        <button
          type="button"
          className="primer-export-button"
          onClick={() => downloadTsv(exportData.contents, `${exportData.target}_primer_pairs.tsv`)}
        >
          <Download size={15} /> Export TSV
        </button>
      </header>
      <div className="compact-table-wrap primer-table-wrap">
      <table className={`compact-table primer-table primer-result-table ${isCdna ? 'cdna-primer-result-table' : ''} ${isCloning ? 'cloning-primer-result-table' : ''}`}>
        <thead ref={resultTableHeadRef}>
          <tr>
            <th>No.</th>
            <th>Primer pair</th>
            <th>{isCloning ? 'Full-length CDS' : isCdna ? 'Expected product' : 'Amplicon'}</th>
            <th>Tm / GC</th>
            {isCloning ? <th>Coverage</th> : null}
            {isCloning ? <th>Additions</th> : null}
            {hasSpecificity ? <th><span className="primer-table-heading">Forward<br />hits <PrimerResultHelp id="forward-hits" text="Detected exact genome matches. ≥2 means two or more; this is not a transcriptome screen." openHelp={openHeaderHelp} onToggle={setOpenHeaderHelp} /></span></th> : null}
            {hasSpecificity ? <th><span className="primer-table-heading">Reverse<br />hits <PrimerResultHelp id="reverse-hits" text="Detected exact genome matches. ≥2 means two or more; this is not a transcriptome screen." openHelp={openHeaderHelp} onToggle={setOpenHeaderHelp} /></span></th> : null}
            {hasSpecificity ? <th><span className="primer-table-heading">Pair<br />products <PrimerResultHelp id="pair-products" text="Exact genome products within the screening limit. ≥2 means two or more. Limited means the expected genomic product is outside that limit." openHelp={openHeaderHelp} onToggle={setOpenHeaderHelp} /></span></th> : null}
            {isCdna ? <th className="primer-exon-heading" title="Annotated exon context for the forward and reverse primers on the selected spliced transcript.">Exon<br />context</th> : null}
          </tr>
        </thead>
        <tbody>
          {visiblePrimers.map((primer) => (
            <tr key={primer.rank}>
              <td>{primer.rank}</td>
              <td className="primer-sequences-cell">
                <div><span>Forward</span><code>{primer.left_sequence}</code><button type="button" title="Copy forward primer" onClick={() => void copyPrimer(`F${primer.rank}`, String(primer.left_sequence ?? ''))}>{copied === `F${primer.rank}` ? <Check size={14} /> : <Copy size={14} />}</button></div>
                <div><span>Reverse</span><code>{primer.right_sequence}</code><button type="button" title="Copy reverse primer" onClick={() => void copyPrimer(`R${primer.rank}`, String(primer.right_sequence ?? ''))}>{copied === `R${primer.rank}` ? <Check size={14} /> : <Copy size={14} />}</button></div>
              </td>
              <td className="primer-amplicon-cell"><PrimerProductTrack primer={primer} targetLength={targetLength} /><span>{isCdna ? 'cDNA ' : ''}{Number(primer.product_size ?? 0).toLocaleString()} bp</span>{isCdna ? <span title="Expected genomic product at the annotated locus, including introns. A dash means no contiguous pair can be assigned (for example, a junction-spanning primer), or no prediction is available. This is separate from genome-wide specificity screening.">gDNA {(primer.expected_genomic_product_size ?? primer.specificity?.expected_genomic_product_size) != null ? `${Number(primer.expected_genomic_product_size ?? primer.specificity?.expected_genomic_product_size).toLocaleString()} bp` : '—'}</span> : null}</td>
              <td className="primer-metrics-cell"><span>F&nbsp; {Number(primer.left_tm ?? 0).toFixed(2)} C &middot; {Number(primer.left_gc_percent ?? 0).toFixed(1)}%</span><span>R&nbsp; {Number(primer.right_tm ?? 0).toFixed(2)} C &middot; {Number(primer.right_gc_percent ?? 0).toFixed(1)}%</span>{isCloning && primer.warnings?.length ? <span className="primer-tm-warning" title={primer.warnings.join(' ')}>Tm difference {Number(primer.tm_difference).toFixed(1)} C &gt; 5 C</span> : null}</td>
              {isCloning ? <td className="primer-coverage-cell">{Number(primer.full_length_coverage ?? primer.coverage ?? 100).toFixed(0)}% CDS</td> : null}
              {isCloning ? <td><PrimerAdditionsCell primer={primer} /></td> : null}
              {hasSpecificity ? <td><span className={`primer-specificity ${primerSpecificityClass(primer.specificity?.left_hit_count)}`}>{primerSpecificityLabel(primer.specificity?.left_hit_count, primer.specificity)}</span></td> : null}
              {hasSpecificity ? <td><span className={`primer-specificity ${primerSpecificityClass(primer.specificity?.right_hit_count)}`}>{primerSpecificityLabel(primer.specificity?.right_hit_count, primer.specificity)}</span></td> : null}
              {hasSpecificity ? <td><span title={primerPairScreenNote(primer.specificity ?? {})} className={`primer-specificity ${primer.specificity?.status === 'target_confirmed' ? 'unique' : 'not-assessed'}`}>{primer.specificity?.status === 'target_out_of_range' ? 'Limited' : primerSpecificityLabel(primer.specificity?.pair_amplicon_count, primer.specificity)}</span></td> : null}
              {isCdna ? <td className="primer-exon-context">{primerJunctionLabel(primer)}</td> : null}
            </tr>
          ))}
        </tbody>
      </table>
      </div>
      {!isCloning && totalPages > 1 ? (
        <div className="compact-table-pagination primer-result-pagination">
          <span>{currentPage * pageSize + 1}-{Math.min((currentPage + 1) * pageSize, primers.length)} of {primers.length}</span>
          <div>
            <button type="button" aria-label="Previous primer page" title="Previous page" disabled={currentPage === 0} onClick={() => setPage((value) => Math.max(0, value - 1))}><ChevronLeft size={16} /></button>
            <span>{currentPage + 1} / {totalPages}</span>
            <button type="button" aria-label="Next primer page" title="Next page" disabled={currentPage >= totalPages - 1} onClick={() => setPage((value) => Math.min(totalPages - 1, value + 1))}><ChevronRight size={16} /></button>
          </div>
        </div>
      ) : null}
    </section>
  );
}

function PrimerResultHelp({ id, text, openHelp, onToggle }: { id: string; text: string; openHelp: string | null; onToggle: (id: string | null) => void }) {
  const isOpen = openHelp === id;
  return (
    <span className={`primer-result-help ${id}`}>
      <button type="button" aria-label={`About ${id.replace('-', ' ')}`} aria-expanded={isOpen} onClick={() => onToggle(isOpen ? null : id)}><CircleHelp size={14} aria-hidden="true" /></button>
      {isOpen ? <span className="primer-result-tooltip" role="tooltip">{text}</span> : null}
    </span>
  );
}

function formatCell(value: unknown, column = '') {
  if (value === null || value === undefined || value === '') return '-';
  if (typeof value === 'number') {
    const integerColumns = new Set(['rank', 'start', 'end', 'length', 'alignment_length', 'product_size', 'support_count', 'go_terms']);
    if (Number.isInteger(value) && integerColumns.has(column)) return value.toLocaleString();
    if (Math.abs(value) > 9999 || (Math.abs(value) > 0 && Math.abs(value) < 0.001)) return value.toExponential(2);
    return Number.isInteger(value) ? value.toLocaleString() : value.toFixed(3);
  }
  if (typeof value === 'object') return JSON.stringify(value);
  return String(value);
}

function GenomePreview() {
  return (
    <div className="viewer-frame">
      <div className="chromosome-ruler">
        <span />
        <span />
        <span />
        <span />
        <span />
      </div>
      <div className="track annotation-track">
        <span style={{ width: '18%', left: '9%' }} />
        <span style={{ width: '26%', left: '37%' }} />
        <span style={{ width: '12%', left: '74%' }} />
      </div>
      <div className="track signal-track">
        <i style={{ height: '28%', left: '7%' }} />
        <i style={{ height: '48%', left: '13%' }} />
        <i style={{ height: '38%', left: '21%' }} />
        <i style={{ height: '70%', left: '42%' }} />
        <i style={{ height: '44%', left: '55%' }} />
        <i style={{ height: '62%', left: '68%' }} />
        <i style={{ height: '31%', left: '81%' }} />
      </div>
      <div className="viewer-labels">
        <span>Annotation</span>
        <span>RNA-seq signal</span>
      </div>
    </div>
  );
}

function GenomeLoading() {
  return (
    <div className="browser-loading" aria-live="polite">
      <div>
        <span />
        <span />
        <span />
      </div>
      <strong>Loading genome browser</strong>
      <p>Preparing JBrowse runtime configuration and the default gene location.</p>
    </div>
  );
}

createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <BasketProvider><App /></BasketProvider>
  </React.StrictMode>,
);
