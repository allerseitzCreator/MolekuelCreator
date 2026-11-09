#!/usr/bin/env python3
"""Baut die einzelne Datei index.html aus src/molekuelzeichner.html und src/chem.js.

Aufruf:  python3 build.py
"""
from pathlib import Path

root = Path(__file__).resolve().parent
src = (root / 'src' / 'molekuelzeichner.html').read_text(encoding='utf-8')
chem = (root / 'src' / 'chem.js').read_text(encoding='utf-8').replace('</script', '<\\/script')
marker = '<script>/*CHEM*/</script>'
assert src.count(marker) == 1, 'Platzhalter für chem.js nicht gefunden'
(root / 'index.html').write_text(src.replace(marker, '<script>\n' + chem + '\n</script>'), encoding='utf-8')
print('index.html erstellt')
