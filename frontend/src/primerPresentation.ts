type RecordValue = Record<string, any>;

export function primerEmptyReason(results: RecordValue): string {
  const reason = String(results.empty_reason ?? 'No primer pairs met the selected constraints. Adjust the parameters and run again.');
  const diagnostics = results.diagnostics;
  const lines: string[] = [];
  const labels: Record<string, string> = {
    'GC below minimum': 'below minimum GC',
    'GC above maximum': 'above maximum GC',
    'Tm below minimum': 'below minimum Tm',
    'Tm above maximum': 'above maximum Tm',
    'homopolymer longer than 5 nt': 'with a homopolymer longer than 5 nt',
    'hairpin Tm above 47 C': 'with hairpin Tm above 47 C',
    'homodimer Tm above 47 C': 'with homodimer Tm above 47 C',
  };
  for (const [key, title] of [['forward', 'Forward'], ['reverse', 'Reverse']]) {
    const end = diagnostics?.[key];
    if (!end || end.accepted !== 0) continue;
    const failures = Object.entries(end.rejected ?? {})
      .filter(([, count]) => Number(count) > 0)
      .map(([label, count]) => `${count} candidate${Number(count) === 1 ? '' : 's'} ${labels[label] ?? label}`);
    if (failures.length) lines.push(`${title}: ${failures.join('; ')}.`);
  }
  if (!lines.length) {
    for (const [label, count] of Object.entries(diagnostics?.pair_rejections ?? {})) {
      if (Number(count) > 0) lines.push(`Pair screening: ${count} pair${Number(count) === 1 ? '' : 's'} rejected (${label}).`);
    }
  }
  return lines.length
    ? ['No full-length CDS primer pairs.', ...lines].join('\n')
    : reason.replace('No full-length CDS primer pairs. ', 'No full-length CDS primer pairs.\n');
}

export function normalizePrimerRequest(input: RecordValue): RecordValue {
  const request = { ...input };
  if (request.source === 'sequence') {
    request.template_type = 'gdna';
    request.genomic_region = 'gene_body';
    delete request.gene_id;
    delete request.transcript_id;
    request.sequence = String(request.sequence ?? '').split(/\r?\n/)
      .filter((line) => !line.trim().startsWith('>')).join('').replace(/\s+/g, '').toUpperCase();
  } else {
    delete request.sequence;
    if (request.template_type !== 'gdna') request.genomic_region = 'gene_body';
  }
  if (request.template_type === 'gdna') delete request.transcript_id;
  if (request.source !== 'gene' || request.template_type !== 'gdna') {
    delete request.flank;
    delete request.tss_upstream;
    delete request.tss_downstream;
  }
  if (request.template_type !== 'cdna') delete request.junction_mode;
  if (request.template_type === 'cdna_cloning') {
    request.check_specificity = false;
    delete request.product_size_range;
  } else {
    delete request.stop_codon_mode;
    delete request.forward_tail;
    delete request.reverse_tail;
  }
  return request;
}

export function annotatedTranscriptLength(features: RecordValue[], transcriptId: string, featureType: string): number | null {
  const intervals = new Map<string, number>();
  for (const feature of features) {
    if (feature.transcript_id !== transcriptId || feature.feature_type !== featureType) continue;
    const start = Number(feature.start), end = Number(feature.end);
    if (Number.isFinite(start) && Number.isFinite(end) && end >= start) intervals.set(`${start}:${end}`, end - start + 1);
  }
  return intervals.size ? [...intervals.values()].reduce((sum, length) => sum + length, 0) : null;
}

export function primerHitLabel(value: unknown, capped = false, exportText = false): string {
  if (value === null || value === undefined || !Number.isFinite(Number(value))) return 'Not assessed';
  return capped ? `${exportText ? '>=' : '\u2265'}${Number(value)}` : Number(value).toLocaleString('en-US');
}

export function primerExport(job: RecordValue, headers: string[], rows: unknown[][]) {
  const target = job.results?.target ?? job.target ?? {};
  const design = job.design ?? {};
  const context = {
    'Job ID': job.job_id ?? '', Species: job.species_id ?? design.species_id ?? '',
    Gene: target.target_id ?? '', Transcript: target.transcript_id ?? '',
    'Template mode': target.genomic_region ?? design.template_type ?? '',
    'Logic version': job.results?.logic_version ?? 1,
    'Design engine': job.results?.design_engine ?? 'Legacy',
    'Design parameters (JSON)': JSON.stringify(design),
  };
  const clean = (value: unknown) => String(value ?? '').replace(/[\t\r\n]/g, ' ');
  return {
    target: [context.Species, context.Gene, context.Transcript, context['Template mode']].filter(Boolean).join('_'),
    contents: [[...Object.keys(context), ...headers], ...rows.map((row) => [...Object.values(context), ...row])]
      .map((row) => row.map(clean).join('\t')).join('\n'),
  };
}
