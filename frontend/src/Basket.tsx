import { createContext, useContext, useEffect, useRef, useState, type ReactNode } from 'react';
import { createPortal } from 'react-dom';
import { Check, Search, ShoppingBasket, Trash2, X } from 'lucide-react';
import type { SpeciesSummary } from './api/client';
import './basket.css';

type Entry = { speciesId: string; geneId: string };
const storageKey = 'polargdb.gene-basket.v1';
function readBasket(): Entry[] {
  try {
    const value: unknown = JSON.parse(localStorage.getItem(storageKey) ?? '[]');
    if (!Array.isArray(value)) return [];
    const entries = value.filter((item): item is Entry => typeof item?.speciesId === 'string' && !!item.speciesId.trim() && typeof item?.geneId === 'string' && !!item.geneId.trim());
    return entries.filter((item, i) => entries.findIndex(other => other.speciesId === item.speciesId && other.geneId === item.geneId) === i);
  } catch { return []; }
}
const BasketContext = createContext<{
  entries: Entry[]; toggle: (speciesId: string, geneId: string) => void; clear: () => void; storageError: boolean;
} | null>(null);
export function BasketProvider({ children }: { children: ReactNode }) {
  const [entries, setEntries] = useState(readBasket);
  const [storageError, setStorageError] = useState(false);
  useEffect(() => {
    try { localStorage.setItem(storageKey, JSON.stringify(entries)); setStorageError(false); }
    catch { setStorageError(true); }
  }, [entries]);
  const toggle = (speciesId: string, geneId: string) => {
    if (!speciesId.trim() || !geneId.trim()) return;
    setEntries(current => current.some(item => item.speciesId === speciesId && item.geneId === geneId)
      ? current.filter(item => item.speciesId !== speciesId || item.geneId !== geneId)
      : [...current, { speciesId, geneId }]);
  };
  return <BasketContext.Provider value={{ entries, toggle, clear: () => setEntries([]), storageError }}>{children}</BasketContext.Provider>;
}
function useBasket() {
  const basket = useContext(BasketContext);
  if (!basket) throw new Error('BasketProvider is required');
  return basket;
}
export function BasketButton({ speciesId, geneId, className, label = false }: { speciesId: string; geneId: string; className?: string; label?: boolean }) {
  const { entries, toggle } = useBasket();
  const added = entries.some(item => item.speciesId === speciesId && item.geneId === geneId);
  const title = added ? 'Remove from basket' : 'Add to basket';
  return <button type="button" className={`${className ?? ''} gene-basket-action${added ? ' is-added' : ''}`} aria-pressed={added}
    aria-label={`${added ? 'Remove' : 'Add'} ${geneId} ${added ? 'from' : 'to'} basket`} title={title}
    onClick={() => toggle(speciesId, geneId)}>
    {added ? <Check size={label ? 15 : 18} aria-hidden="true" /> : <ShoppingBasket size={label ? 15 : 18} aria-hidden="true" />}
    {label && (added ? 'In basket' : 'Add basket')}
  </button>;
}

export function BasketMenu({ species, onSearch }: { species: SpeciesSummary[]; onSearch: (speciesId: string, geneIds: string[]) => void }) {
  const { entries, toggle, clear, storageError } = useBasket();
  const dialog = useRef<HTMLDialogElement>(null);
  const [confirmClear, setConfirmClear] = useState(false);
  const groups = [...new Set(entries.map(item => item.speciesId))];
  const close = () => { dialog.current?.close(); setConfirmClear(false); };
  return <>
    <button type="button" className="icon-nav-button basket-nav-button" aria-label={`Basket (${entries.length} genes)`} data-tooltip="Basket" onClick={() => dialog.current?.showModal()}>
      <ShoppingBasket size={21} />{entries.length > 0 && <span className="basket-count" aria-hidden="true">{entries.length}</span>}
    </button>
    {createPortal(<dialog ref={dialog} className="basket-dialog" aria-labelledby="basket-title" onClose={() => setConfirmClear(false)} onClick={event => { if (event.target === event.currentTarget) { const rect = event.currentTarget.getBoundingClientRect(); if (event.clientX < rect.left || event.clientX > rect.right || event.clientY < rect.top || event.clientY > rect.bottom) close(); } }}>
      <header className="basket-heading"><div><ShoppingBasket size={23} /><h2 id="basket-title">Gene basket</h2><span aria-live="polite">{entries.length} genes</span></div><button type="button" onClick={close} aria-label="Close basket" title="Close"><X size={22} /></button></header>
      {storageError && <p className="basket-warning" role="alert">Browser storage is unavailable. This basket will not be saved after reloading.</p>}
      <div className="basket-content">
        {!entries.length && <p className="basket-empty">Your basket is empty.</p>}
        {groups.map(speciesId => <section className="basket-group" key={speciesId}>
          <header><h3><em>{species.find(item => item.species_id === speciesId)?.scientific_name ?? speciesId}</em></h3>
            <button type="button" className="basket-search" onClick={() => { onSearch(speciesId, entries.filter(item => item.speciesId === speciesId).map(item => item.geneId)); close(); }}><Search size={16} />Use in Gene Search</button></header>
          <ul>{entries.filter(item => item.speciesId === speciesId).map(item => <li key={item.geneId}><span>{item.geneId}</span><button type="button" title="Remove from basket" aria-label={`Remove ${item.geneId} from basket`} onClick={() => toggle(item.speciesId, item.geneId)}><X size={17} /></button></li>)}</ul>
        </section>)}
      </div>
      {!!entries.length && <footer className="basket-footer">{confirmClear ? <><span>Remove all genes?</span><button type="button" onClick={() => { clear(); setConfirmClear(false); }}>Clear basket</button><button type="button" onClick={() => setConfirmClear(false)}>Cancel</button></> : <button type="button" onClick={() => setConfirmClear(true)}><Trash2 size={16} />Clear basket</button>}</footer>}
    </dialog>, document.body)}
  </>;
}
