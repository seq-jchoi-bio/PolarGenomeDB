"""Publish complete job records so polling never reads a partial JSON file."""
import json
from pathlib import Path
import tempfile


def write_json(path: Path, payload: dict) -> None:
    content = json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + '\n'
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                         prefix=f'.{path.name}-', suffix='.tmp', delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(content)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
