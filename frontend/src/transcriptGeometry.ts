type Interval = { start: number; end: number };

// Annotation coordinates are inclusive; display offsets are half-open.
export function spliceExons(exons: Interval[], strand: string) {
  let offset = 0;
  const blocks = [...exons]
    .sort((a, b) => strand === '-' ? b.start - a.start : a.start - b.start)
    .map((exon) => {
      const start = offset;
      offset += exon.end - exon.start + 1;
      return { ...exon, offsetStart: start, offsetEnd: offset };
    });
  return { blocks, length: offset };
}

export function projectSplicedFeature(
  feature: Interval,
  blocks: ReturnType<typeof spliceExons>['blocks'],
  strand: string,
) {
  return blocks.flatMap((block) => {
    const start = Math.max(feature.start, block.start);
    const end = Math.min(feature.end, block.end);
    if (start > end) return [];
    const offsetStart = block.offsetStart + (strand === '-' ? block.end - end : start - block.start);
    return [{ start: offsetStart, end: offsetStart + end - start + 1 }];
  });
}
