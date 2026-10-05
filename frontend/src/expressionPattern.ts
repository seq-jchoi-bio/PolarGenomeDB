// Display heuristics for the sampled conditions, not a tissue-specificity test.
export const EXPRESSION_PATTERN_POLICY = Object.freeze({
  minimumGroupCpm: 1,
  maximumBroadRatio: 1.5,
});

export function classifyExpressionPattern(groups: { label: string; mean: number }[]) {
  if (!groups.length || groups.some((group) => !Number.isFinite(group.mean) || group.mean < 0)) {
    return { label: '-', description: 'CPM group means are unavailable.' };
  }
  if (groups.length === 1) {
    return { label: groups[0].label, description: 'Only one condition is available; breadth cannot be assessed.' };
  }
  const maximum = Math.max(...groups.map((group) => group.mean));
  const minimum = Math.min(...groups.map((group) => group.mean));
  const { minimumGroupCpm, maximumBroadRatio } = EXPRESSION_PATTERN_POLICY;
  if (maximum < minimumGroupCpm) {
    return { label: 'Low expression', description: `All observed condition means are below ${minimumGroupCpm} CPM.` };
  }
  if (minimum >= minimumGroupCpm && maximum / minimum <= maximumBroadRatio) {
    return { label: 'Broad', description: `Across observed conditions: all means >= ${minimumGroupCpm} CPM; highest/lowest <= ${maximumBroadRatio}.` };
  }
  const peak = groups.find((group) => group.mean === maximum)!;
  return { label: peak.label, description: 'Highest mean CPM among observed conditions; not a claim of condition-specific expression.' };
}
