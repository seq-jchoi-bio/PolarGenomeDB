# Public source boundary

This release is derived from the working source without editing that source.

## Retained

Analysis and ranking algorithms, parameter defaults, database schema creation,
API and job handling, manifest-driven build logic, public UI, input validation,
file-serving restrictions and request limits are preserved. Synthetic unit tests
are retained. Project-relative runtime path conventions remain because they are
part of the software interface; actual dataset locations/manifests are excluded.

## Minimal changes

- Private absolute tool-installation paths become environment-variable settings
  with active-Python or project-relative defaults.
- Actual manifests are replaced with separately authored examples; these examples
  are not anonymized copies of production metadata.
- Runtime image references point to a neutral placeholder. Original photographs
  and branding image files are not distributed as application assets.
- Documentation includes an approved Home screenshot showing the interface,
  photographs and logo, plus C4 Level 1 and Level 2 PNG diagrams. The README
  displays Level 1 only. Internal architecture diagrams and maintenance notes
  are excluded. Font assets retain upstream notices.
- Public-only setup documentation and a defensive `.gitignore` are added.

## Excluded

Production datasets, databases, caches, jobs, indexes, raw files, logs, reports,
personal/meeting documents, original runtime manifests, administrator and security
infrastructure, credential handling, monitoring/backup deployment, private cluster
scripts, and local UI/benchmark runners are not distributed. Two dataset-specific
smoke/benchmark runners contain actual protein sequences and are excluded entirely.
Generic API safety checks remain; removal of those checks would change behavior.

## Verification

`SOURCE_INVENTORY.json` lists copied file hashes and identifies unchanged,
sanitized and public-only files. No secret values or original absolute locations
are recorded there. Source hashes allow comparison by the maintainer without
publishing excluded source contents. An automated allowlist and scan reduce leak
risk but do not prove that every possible sensitive value has been detected.

This directory is a source snapshot, not a second development branch. Future
exports should be made into a new directory and compared before publication.
Do not recursively copy the entire working project into this directory.
