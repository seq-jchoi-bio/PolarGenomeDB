import geist400 from '@fontsource/geist/files/geist-latin-400-normal.woff2?url';
import geist500 from '@fontsource/geist/files/geist-latin-500-normal.woff2?url';
import geist600 from '@fontsource/geist/files/geist-latin-600-normal.woff2?url';
import geist700 from '@fontsource/geist/files/geist-latin-700-normal.woff2?url';

const fontUrls: Record<string, string> = { '400': geist400, '500': geist500, '600': geist600, '700': geist700 };
const fontData = new Map<string, Promise<string>>();
function embeddedFont(weight: string): Promise<string> {
  if (!fontData.has(weight)) {
    const pending = fetch(fontUrls[weight]).then(async (response) => {
      if (!response.ok) throw new Error('Could not embed the plot font. Please retry.');
      const blob = await response.blob();
      return new Promise<string>((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = () => resolve(String(reader.result));
        reader.onerror = () => reject(new Error('Could not read the plot font.'));
        reader.readAsDataURL(blob);
      });
    }).catch((error) => { fontData.delete(weight); throw error; });
    fontData.set(weight, pending);
  }
  return fontData.get(weight)!;
}

export async function exportSvgPlot(svg: SVGSVGElement, format: 'png' | 'svg', filename: string, fontFamily?: string) {
  await document.fonts.ready;
  const clone = svg.cloneNode(true) as SVGSVGElement;
  const originals = [svg, ...svg.querySelectorAll('*')];
  const copies = [clone, ...clone.querySelectorAll('*')];
  const weights = new Set<string>();
  // Inline computed styles so downloaded plots do not depend on application CSS.
  originals.forEach((element, index) => {
    const style = getComputedStyle(element);
    const copy = copies[index] as SVGElement;
    if (!fontFamily && element.tagName.toLowerCase() === 'text' && style.fontFamily.includes('Geist') && fontUrls[style.fontWeight]) weights.add(style.fontWeight);
    for (const property of ['fill', 'fill-opacity', 'fill-rule', 'stroke', 'stroke-width', 'stroke-opacity', 'stroke-linecap', 'stroke-linejoin', 'stroke-dasharray', 'stroke-dashoffset', 'font-family', 'font-size', 'font-weight', 'font-style', 'text-anchor', 'dominant-baseline', 'opacity', 'stop-color', 'stop-opacity', 'clip-path']) {
      // Computed paint references can become absolute document URLs; keep SVG-local IDs.
      const value = style.getPropertyValue(property).replace(/url\(["']?[^)]*#([^"')]+)["']?\)/g, 'url(#$1)');
      copy.style.setProperty(property, value);
    }
    if (fontFamily) copy.style.setProperty('font-family', fontFamily);
  });
  const view = svg.viewBox.baseVal;
  if (!(view.width > 0 && view.height > 0)) throw new Error('The plot has no valid export dimensions.');
  const bounds = svg.getBBox();
  const x = Math.min(view.x, bounds.x < view.x ? bounds.x - 12 : view.x);
  const y = Math.min(view.y, bounds.y < view.y ? bounds.y - 12 : view.y);
  const width = Math.max(view.x + view.width, bounds.x + bounds.width > view.x + view.width ? bounds.x + bounds.width + 12 : view.x + view.width) - x;
  const height = Math.max(view.y + view.height, bounds.y + bounds.height > view.y + view.height ? bounds.y + bounds.height + 12 : view.y + view.height) - y;
  clone.setAttribute('viewBox', `${x} ${y} ${width} ${height}`);
  if (weights.size) {
    const style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
    style.textContent = (await Promise.all([...weights].map(async (weight) => `@font-face{font-family:Geist;font-weight:${weight};font-style:normal;src:url("${await embeddedFont(weight)}") format("woff2");}`))).join('\n');
    clone.prepend(style);
  }
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('width', String(width));
  clone.setAttribute('height', String(height));
  const background = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  background.setAttribute('x', String(x));
  background.setAttribute('y', String(y));
  background.setAttribute('width', String(width));
  background.setAttribute('height', String(height));
  background.setAttribute('fill', '#ffffff');
  clone.prepend(background);
  let blob = new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml;charset=utf-8' });
  if (format === 'png') {
    const source = URL.createObjectURL(blob);
    try {
      const image = new Image();
      image.src = source;
      await image.decode();
      const canvas = document.createElement('canvas');
      canvas.width = Math.ceil(width * 2);
      canvas.height = Math.ceil(height * 2);
      const context = canvas.getContext('2d');
      if (!context) throw new Error('PNG export is unavailable.');
      context.drawImage(image, 0, 0, canvas.width, canvas.height);
      blob = await new Promise<Blob>((resolve, reject) => canvas.toBlob((value) => value ? resolve(value) : reject(new Error('PNG export failed.')), 'image/png'));
    } finally {
      URL.revokeObjectURL(source);
    }
  }
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = `${filename}.${format}`;
  link.click();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export const exportTree = exportSvgPlot;
