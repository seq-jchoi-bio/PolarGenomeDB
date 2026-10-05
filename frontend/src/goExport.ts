type GoRecord = Record<string, unknown>;

export function wrapGoLabel(label: string, maxWidth: number): string[] {
  const context = document.createElement('canvas').getContext('2d');
  const fits = (text: string) => {
    if (!context) return text.length <= 30;
    return ['Geist', 'Arial'].every((font) => {
      context.font = `15px ${font}`;
      return context.measureText(text).width <= maxWidth;
    });
  };
  const lines: string[] = [];
  let line = '';
  for (const word of label.trim().split(/\s+/)) {
    const next = line ? `${line} ${word}` : word;
    if (fits(next)) { line = next; continue; }
    if (line) lines.push(line);
    line = '';
    for (const character of word) {
      if (line && !fits(line + character)) { lines.push(line); line = ''; }
      line += character;
    }
  }
  if (line) lines.push(line);
  return lines.length ? lines : [''];
}

export function goTermsTsv(terms: GoRecord[], namespace: string, annotationMode: string): string {
  const base = ['namespace', 'annotation_mode', 'go_id', 'term', 'level', 'count', 'percent', 'pvalue', 'adjusted_pvalue', 'fold_enrichment', 'pop_hits', 'pop_total'];
  const columns = [...new Set([...base, ...terms.flatMap((term) => Object.keys(term))])];
  const field = (value: unknown) => {
    const text = value == null ? '' : typeof value === 'object' ? JSON.stringify(value) : String(value);
    return /[\t\r\n"]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
  };
  const lines = terms.map((term) => {
    const record = { ...term, namespace, annotation_mode: annotationMode } as GoRecord;
    return columns.map((column) => field(record[column])).join('\t');
  });
  return [columns.join('\t'), ...lines].join('\n') + '\n';
}
