from pathlib import Path

PARTS = ('1', '2', '3', '4', 'rev.1', 'rev.2')


def index_format(prefix: Path) -> str:
    """Require one nonempty six-file Bowtie2 index, with no mixed formats."""
    formats = []
    for ext in ('bt2', 'bt2l'):
        files = [Path(f'{prefix}.{part}.{ext}') for part in PARTS]
        if any(p.exists() for p in files):
            if not all(p.is_file() and p.stat().st_size > 0 for p in files):
                raise ValueError(f'Incomplete or empty Bowtie2 {ext} index: {prefix}')
            formats.append(ext)
    if len(formats) != 1:
        raise ValueError(f'Expected exactly one complete .bt2 or .bt2l index: {prefix}')
    return formats[0]
