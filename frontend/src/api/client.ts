export const API_BASE = import.meta.env.VITE_API_BASE ?? '';

export type SpeciesSummary = {
  species_id: string;
  scientific_name: string;
  display_name: string;
  taxon_group: string;
  assembly_version: string | null;
  counts: Record<string, number>;
  expression_value_types: Array<Record<string, number | string>>;
};

export type Gene = {
  gene_id: string;
  seqid: string;
  start: number;
  end: number;
  strand: string | null;
  source?: string | null;
  attributes?: Record<string, unknown>;
};

export async function apiGet<T>(path: string): Promise<T> {
  const response = await fetch(`${API_BASE}${path}`);
  if (!response.ok) {
    throw new Error(`${response.status} ${response.statusText}`);
  }
  return response.json() as Promise<T>;
}

export async function apiPost<T>(path: string, body: unknown): Promise<T> {
  const response = await fetch(`${API_BASE}${path}`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body),
  });
  if (!response.ok) {
    const detail = await response.json().catch(() => ({}));
    throw new Error(detail.detail ?? `${response.status} ${response.statusText}`);
  }
  return response.json() as Promise<T>;
}
