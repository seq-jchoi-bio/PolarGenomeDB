# TeX Gyre Heros

- Version: 2.004 (2009-10-30).
- Source: https://mirrors.ctan.org/fonts/tex-gyre.zip
- Project: https://www.gust.org.pl/projects/e-foundry/tex-gyre/heros
- Retrieved: 2026-09-21.
- License: GUST Font License, referring to LPPL 1.3c or later.
- Authors and copyright: see the accompanying original manifest and README.

The four non-condensed OpenType files are distributed unmodified. The
original license, manifest and README accompany them. The complete upstream
distribution, including the other files named in the manifest, is available
from the source URL above.

The application serves these files locally through CSS @font-face; no local()
lookup or external font service is used. Regular and italic have weight 400;
bold and bold italic have weight 700. Intermediate CSS weights resolve to
these real faces. Font synthesis is disabled in the application.

Geist is the active UI font. Heros is retained as a design alternative,
selected via --font-ui-alternative. Only sequence displays, sequence input,
and primer sequence columns use IBM Plex Mono. All other UI text, including
identifiers, coordinates and figure labels, uses Geist.
