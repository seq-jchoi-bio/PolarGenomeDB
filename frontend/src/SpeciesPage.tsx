import { useEffect, useState, type ReactNode } from 'react';
import { ArrowUpRight, Check, Layers3, CircleDashed, Unplug } from 'lucide-react';
import { API_BASE, apiGet, type SpeciesSummary } from './api/client';
import './species.css';

const profiles: Record<string, { common: string; description: string; images: string[] }> = {
  Cq: {
    common: 'Antarctic pearlwort',
    description: 'Colobanthus quitensis is a cushion-forming plant and one of the two native flowering plants of Antarctica. It occurs in the maritime Antarctic, including the South Shetland Islands and the western Antarctic Peninsula.',
    images: ['/images/withheld.svg', '/images/withheld.svg'],
  },
  Da: {
    common: 'Antarctic hair grass',
    description: 'Deschampsia antarctica is one of the two native flowering plants of Antarctica. This grass occurs in the maritime Antarctic, including the South Orkney Islands, the South Shetland Islands and the western Antarctic Peninsula.',
    images: ['/images/withheld.svg', '/images/withheld.svg'],
  },
};
const number = (value: number | undefined) => value == null ? 'Not available' : value.toLocaleString('en-US');
type AssetStatus = {
  ready: boolean;
  assets: { expression: Array<{ value_type: string; genes: number; samples: number; value_count: number }> };
};

export function SpeciesPage({ species, selected, onSelect, modules, onOpenTool, speciesImages }: {
  species: SpeciesSummary[];
  selected: SpeciesSummary | undefined;
  onSelect: (id: string) => void;
  modules: Array<{ id: string; title: string; scope: string; icon: ReactNode }>;
  onOpenTool: (id: string) => void;
  speciesImages: Record<string, string>;
}) {
  const [assembly, setAssembly] = useState<{ speciesId: string; bases: number; sequences: number } | null>(null);
  const [assemblyError, setAssemblyError] = useState(false);
  const [assetCheck, setAssetCheck] = useState<{ id: string; data?: AssetStatus; failed?: boolean } | null>(null);
  const id = selected?.species_id;
  const [connection, setConnection] = useState<{ id: string; state: 'connected' | 'unavailable'; species?: SpeciesSummary } | null>(null);
  useEffect(() => {
    if (!id) return;
    let disposed = false;
    let active: AbortController | null = null;
    const check = async () => {
      if (active) return;
      const controller = new AbortController();
      active = controller;
      const timeout = window.setTimeout(() => controller.abort(), 8000);
      try {
        const response = await fetch(`${API_BASE}/api/species`, { signal: controller.signal, cache: 'no-store' });
        if (!response.ok) throw new Error('Species API unavailable');
        const data = await response.json() as { species: SpeciesSummary[] };
        const current = data.species.find(item => item.species_id === id);
        if (!current) throw new Error('Species unavailable');
        if (!disposed) setConnection({ id, state: 'connected', species: current });
        try {
          const assets = await fetch(`${API_BASE}/api/species/${encodeURIComponent(id)}/assets/status`, { signal: controller.signal, cache: 'no-store' });
          if (!assets.ok) throw new Error('Asset check unavailable');
          const data = await assets.json() as AssetStatus;
          if (!disposed) setAssetCheck({ id, data });
        } catch {
          if (!disposed) setAssetCheck({ id, failed: true });
        }
      } catch {
        if (!disposed) setConnection({ id, state: 'unavailable' });
      } finally {
        window.clearTimeout(timeout);
        active = null;
      }
    };
    setConnection(null);
    setAssetCheck(null);
    void check();
    const interval = window.setInterval(check, 30000);
    window.addEventListener('focus', check);
    return () => {
      disposed = true;
      active?.abort();
      window.clearInterval(interval);
      window.removeEventListener('focus', check);
    };
  }, [id]);
  const state = connection?.id === id ? connection?.state : undefined;
  const connected = state === 'connected';
  const liveSpecies = connected ? connection?.species : undefined;
  const statusText = connected ? 'API connected' : state === 'unavailable' ? 'Not connected' : 'Checking...';
  const StatusIcon = connected ? Check : state === 'unavailable' ? Unplug : CircleDashed;
  const resourceStatus = (moduleId: string) => {
    const check = assetCheck?.id === id ? assetCheck : null;
    if (!connected) return { text: statusText, ok: false, error: state === 'unavailable' };
    if (!check) return { text: 'Checking...', ok: false, error: false };
    if (check.failed || !check.data) return { text: 'Check unavailable', ok: false, error: true };
    if (moduleId === 'assets') return { text: check.data.ready ? 'Files ready' : 'Needs attention', ok: check.data.ready, error: !check.data.ready };
    const values = check.data.assets.expression;
    if (!values.length) return { text: 'Not provided', ok: false, error: false };
    const complete = ['raw_count', 'cpm'].every(type => values.some(value => value.value_type === type
      && value.genes > 0 && value.samples > 0 && value.value_count === value.genes * value.samples));
    return { text: complete ? 'Available' : 'Incomplete', ok: complete, error: !complete };
  };
  useEffect(() => {
    let cancelled = false;
    setAssembly(null);
    setAssemblyError(false);
    if (id) apiGet<{ chromosomes: Array<{ length: number }> }>(`/api/species/${encodeURIComponent(id)}/chromosomes`)
      .then(({ chromosomes }) => {
        if (!cancelled) setAssembly({ speciesId: id, bases: chromosomes.reduce((sum, row) => sum + row.length, 0), sequences: chromosomes.length });
      }).catch(() => { if (!cancelled) setAssemblyError(true); });
    return () => { cancelled = true; };
  }, [id, connected]);

  const profile = id ? profiles[id] : undefined;
  const currentAssembly = connected && assembly?.speciesId === id ? assembly : null;
  const facts: Array<[string, string]> = selected ? [
    ['Assembly version', selected.assembly_version && selected.assembly_version !== 'unknown' ? selected.assembly_version : 'Not specified'],
    ['Assembly span', currentAssembly ? `${(currentAssembly.bases / 1e6).toLocaleString('en-US', { maximumFractionDigits: 2 })} Mb` : assemblyError || state === 'unavailable' ? 'Unavailable' : 'Loading...'],
    ['Assembly sequences', currentAssembly ? number(currentAssembly.sequences) : assemblyError || state === 'unavailable' ? 'Unavailable' : 'Loading...'],
    ['Annotated genes', number(liveSpecies?.counts.genes)],
    ['Transcripts', number(liveSpecies?.counts.transcripts)],
    ['Protein records', number(liveSpecies?.counts.proteins)],
    ['RNA-seq samples', number(liveSpecies?.counts.samples)],
    ['Coverage tracks', number(liveSpecies?.counts.bigwig_tracks)],
  ] : [];

  return (
    <section className="species-page page-section" aria-label="Species resources">
      <div className="species-page-inner">
        <header className="species-page-heading">
          <h1>Species</h1>
        </header>
        <div className="species-page-layout">
          <nav className="species-switcher" aria-label="Select species">
            {species.map((item) => {
              const [genus, ...rest] = item.scientific_name.split(' ');
              const shortName = rest.length ? `${genus[0]}. ${rest.join(' ')}` : genus;
              return (
              <button key={item.species_id} type="button" aria-pressed={item.species_id === id} onClick={() => onSelect(item.species_id)}>
                <img src={speciesImages[item.species_id] ?? '/images/withheld.svg'} alt="" />
                <span><em>{item.scientific_name}</em><small><em>{shortName}</em> · {item.taxon_group}</small></span>
              </button>
            ); })}
          </nav>
        <div className="species-page-content">
        {!selected ? <p role="status">Species information is not available.</p> : <>
          <header className="species-introduction">
            <div className="species-naming"><h2><em>{selected.scientific_name}</em></h2>
              <span>{selected.taxon_group === 'monocot' ? 'Monocot' : selected.taxon_group === 'dicot' ? 'Dicot' : selected.taxon_group}</span>
            </div>
            {profile && <><p className="species-common-name">{profile.common}</p><p className="species-description">{profile.description}</p></>}
          </header>
          {profile && <div className="species-photographs">
            <div className="species-photo-pair">
              {profile.images.map((src, index) => <figure key={src}>
                <img src={src} alt={`${selected.scientific_name}, field photograph ${index + 1}`} width="1600" height="1067" />
                <figcaption className="image-credit">© Polar Genome DB</figcaption>
              </figure>)}
            </div>
          </div>}

          <section className="species-resource-band" aria-labelledby="species-genome-title">
            <div className="species-band-intro">
              <h2 id="species-genome-title">Genome resource</h2>
              <p>In-house genome assembly and annotation, with linked protein sequences and RNA-seq data.</p>
              <button type="button" className="species-browser-link" onClick={() => onOpenTool('genome')}><Layers3 size={18} />Open genome browser<ArrowUpRight size={17} /></button>
            </div>
            <div>
              <dl className="species-genome-facts">{facts.map(([label, value]) => <div key={label}><dt>{label}</dt><dd>{value}</dd></div>)}</dl>
              <p className="species-genome-note">{connected ? 'Counts reflect the current database release. Assembly span is the total indexed sequence length, including any unplaced sequences.' : 'Database statistics cannot be verified until the backend connection is available.'}</p>
              <p className="species-reference-note">Functional reference <span>{selected.taxon_group === 'monocot' ? <><em>Oryza sativa</em> · RAP-DB IRGSP-1.0</> : selected.taxon_group === 'dicot' ? <><em>Arabidopsis thaliana</em> · TAIR10</> : 'Not assigned'}</span></p>
            </div>
          </section>

          <section className="species-resource-band species-modules-band" aria-labelledby="species-modules-title">
            <div className="species-band-intro"><h2 id="species-modules-title">Analysis modules</h2><p>Tools and data resources for <em>{selected.scientific_name}</em>.</p><p role="status">{connected ? 'Resource status reflects matrix records and required file availability, not full integrity or analysis validation.' : state === 'unavailable' ? 'Backend connection unavailable. Analysis and data readiness cannot be verified.' : 'Checking backend connection...'}</p></div>
            <div className="species-module-list">
              {modules.map((module) => {
                const isTool = !['expression', 'assets'].includes(module.id);
                const resource = resourceStatus(module.id);
                const available = isTool ? connected : resource.ok;
                const disconnected = isTool ? state === 'unavailable' : resource.error;
                const Icon = isTool ? StatusIcon : available ? Check : disconnected ? Unplug : CircleDashed;
                return <div className="species-module-row" key={module.id}>
                  <span className="species-module-icon" aria-hidden="true">{module.icon}</span>
                  <div>{isTool ? <button type="button" onClick={() => onOpenTool(module.id)}>{module.title}<ArrowUpRight size={16} /></button> : <span className="species-module-title">{module.title}</span>}
                    <p>{module.id === 'primer' ? 'Primer3 design / Bowtie2 specificity' : module.scope}</p></div>
                  <span className={`species-module-status${available ? ' is-available' : disconnected ? ' is-disconnected' : ''}`}><Icon size={17} />{isTool ? statusText : resource.text}</span>
                </div>;
              })}
            </div>
          </section>
        </>}
        </div>
        </div>
      </div>
    </section>
  );
}
