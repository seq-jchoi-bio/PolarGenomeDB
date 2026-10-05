# Getting Started

Source snapshot of the multi-species genomics portal. Analysis implementations,
API routes, schema creation, build scripts, UI and synthetic unit tests are included.
This is not a copy of the production datasets or an operational server backup.

## Scope

- No genome, protein, annotation, RNA-seq, SQLite, runtime cache or index files.
- No production manifests, credentials, server addresses, private documents,
  administrator portal, monitoring, backup infrastructure or SSH/Slurm scripts.
- Photographs and branding images are replaced by an explicit placeholder.
- Example gene identifiers and public scientific descriptions remain in the UI;
  examples only work with matching datasets supplied independently.
- Original source code and technical documentation use the MIT License. See
  `LICENSE` and `NOTICE.md` for scope, image exclusions and third-party terms.

See `SOURCE_CHANGES.md` and `SOURCE_INVENTORY.json` for the release boundary.

## Local setup

Use Python 3.12 or newer and a Node.js version supported by Vite 7.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
cd frontend
npm ci
npm run build
cd ..
```

The original `environment.yml` is included, but is a minimal environment, not a
complete dependency lock. Analyses/builds additionally require the applicable
tools: BLAST+, DIAMOND, HMMER, SAMtools/HTSlib, deepTools, Primer3, Bowtie2,
MAFFT and IQ-TREE 3. Semantic GO requires GoFigure! and its Python dependencies
in a separate environment. External binaries and reference databases are not bundled.

Configure `ANTARCTIC_ORTHOLOG_BIN` to the directory containing analysis binaries.
The portable default is the active Python executable's directory. Set
`ANTARCTIC_PFAM_DB`, `ANTARCTIC_GOFIGURE_PYTHON` and `ANTARCTIC_GOFIGURE_SCRIPT`
explicitly for those modules. `.env.example` lists placeholder paths only;
export variables in your shell or supply them through your process manager.
The application does not automatically read `.env` files.

```sh
python -m uvicorn backend.main:app --host 127.0.0.1 --port 8770
```

The built UI is served by the backend. Alternatively run `npm run dev` inside
`frontend/` for the development UI at port 5174 (API proxy: port 8770).
Without a prepared database, the health endpoint and static UI are available,
but data queries and analyses are not. No live production API is used as a fallback.

## Dataset preparation

Manifest paths are relative to the manifest file unless absolute. Copy the
example manifests and supply your own authorized inputs. Runtime database
and asset paths are local configuration, not distributed production paths.

```sh
python scripts/build_species.py --manifest manifests/species.example.yaml
python scripts/build_reference.py --help
python scripts/build_go_ontology.py --help
python scripts/build_all.py --help
```

Build core species metadata first. Prepare the appropriate Arabidopsis/rice
reference, GO ontology and Pfam profiles before the full pipeline. The retained
reference IDs in the examples match the current monocot/dicot selection logic.
Inspect script `--help` for reference annotations and version-specific formats.

`build_all.py` requires an existing source DB with reference/GO resources and
real manifest inputs, even in plan-only mode. It builds into a separate staging
directory: species metadata, BLAST/DIAMOND indexes, Bowtie2 index registration,
Pfam, similarity, optional expression, browser assets and runtime indexes.
It does not download reference resources, build Bowtie2 indexes, or deploy to a server.
Bowtie2 indexes must be prepared separately. Expression is optional; raw count
inputs are imported, not recalculated from BAMs. Coverage generation needs BAMs
and deepTools. Inspect and validate the staging output before operational use.

## Tests and reproducibility

```sh
python -m pip install httpx
python -m unittest discover -s backend/tests
```

Tests use temporary or mocked records; external tool availability can affect
coverage. Exact reproduction of biological results requires the original data,
reference versions, tool versions and settings. Those data are not provided here.
This snapshot preserves scientific code, not a complete reproducibility archive.

The production-specific privacy/contact statements, institutional attributions,
species descriptions and example IDs remain for source transparency. Adapt them
before operating an independent site. This snapshot is not production deployment
guidance: configure authentication, HTTPS, rate limits and private asset access
for your own environment before exposing a service.
