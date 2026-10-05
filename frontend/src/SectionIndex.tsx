import { useEffect, useRef, useState } from 'react';

type Props = {
  sections: readonly (readonly [string, string])[];
  prefix?: string;
  className: string;
  label: string;
  heading?: string;
};

export function SectionIndex({ sections, prefix = '', className, label, heading }: Props) {
  const [active, setActive] = useState<string>(sections[0]?.[0] ?? '');
  const selectedAnchor = useRef<string | null>(null);
  useEffect(() => {
    let frame = 0;
    const update = () => {
      frame = 0;
      const scrollTop = window.scrollY;
      const maxScroll = Math.max(0, document.documentElement.scrollHeight - window.innerHeight);
      const anchors = sections.flatMap(([id]) => {
        const element = document.getElementById(`${prefix}${id}`);
        if (!element) return [];
        const offset = Number.parseFloat(getComputedStyle(element).scrollMarginTop) || 0;
        return [{ id, top: Math.max(0, Math.min(maxScroll, element.getBoundingClientRect().top + scrollTop - offset)) }];
      });
      // Short sections can share a bottom-clamped destination: respect clicks.
      const selected = anchors.find(({ id }) => id === selectedAnchor.current);
      if (selected && Math.abs(selected.top - scrollTop) <= 3) {
        setActive(selected.id);
        return;
      }
      selectedAnchor.current = null;
      if (scrollTop <= 1 || maxScroll <= 1) {
        setActive(sections[0]?.[0] ?? '');
        return;
      }
      let nearest = anchors[0];
      for (const anchor of anchors) {
        const distance = Math.abs(anchor.top - scrollTop);
        const previousDistance = nearest ? Math.abs(nearest.top - scrollTop) : Infinity;
        if (distance < previousDistance || (distance === previousDistance && scrollTop >= maxScroll - 2)) nearest = anchor;
      }
      if (nearest) setActive(nearest.id);
    };
    const schedule = () => { if (!frame) frame = requestAnimationFrame(update); };
    const syncHash = () => {
      const section = sections.find(([id]) => window.location.hash === `#${prefix}${id}`);
      if (section) {
        selectedAnchor.current = section[0];
        setActive(section[0]);
      }
      schedule();
    };
    syncHash();
    window.addEventListener('scroll', schedule, { passive: true });
    window.addEventListener('resize', schedule);
    window.addEventListener('hashchange', syncHash);
    const observer = new ResizeObserver(schedule);
    observer.observe(document.documentElement);
    return () => {
      cancelAnimationFrame(frame);
      window.removeEventListener('scroll', schedule);
      window.removeEventListener('resize', schedule);
      window.removeEventListener('hashchange', syncHash);
      observer.disconnect();
    };
  }, [sections, prefix]);
  return <nav className={className} aria-label={label}>
    {heading && <span>{heading}</span>}
    {sections.map(([id, text], index) => <a key={id} href={`#${prefix}${id}`} onClick={(event) => {
      if (event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
      selectedAnchor.current = id;
      setActive(id);
    }} aria-current={active === id ? 'location' : undefined}>
      <small aria-hidden="true">{String(index + 1).padStart(2, '0')}</small>
      {text}
      <i aria-hidden="true" />
    </a>)}
  </nav>;
}
