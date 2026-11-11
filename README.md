# MolekülCreator

Ein Molekül-Zeichenprogramm für den Chemieunterricht, das komplett im Browser läuft. Es ist eine einzige HTML-Datei ohne Installation, ohne Server und ohne externe Bibliotheken.

## Starten

- `index.html` herunterladen und im Browser öffnen, **oder**
- GitHub Pages für dieses Repository aktivieren (*Settings → Pages → Branch `main`, Ordner `/ (root)`*). Die App ist dann unter `https://<benutzer>.github.io/MolekuelCreator/` erreichbar.

## Funktionen

**Zeichnen**
- Atome, Einfach-, Doppel- und Dreifachbindungen, Ladungen, Radikale, freie Elektronenpaare, Farben
- Bindungslängen mit Arretierung: sehr kurz, kurz, mittel, lang
- CH₂-Einheit in eine markierte Bindung einschieben
- Moleküle duplizieren (Strg+D)
- Gruppieren sowie Rückgängig und Wiederholen

**Darstellungen**
- Skelettformel, Halbstrukturformel und Valenzstrichformel mit allen H-Atomen
- Fischer- und Haworth-Projektion für Zucker
- Keil-Strich-Schreibweise in der Skelettformel, in der Valenzstrichformel von Hand setzbar
- R/S-Kennzeichnung, Hybridisierung und Aromatizität

**Namen und Formeln**
- Struktur aus dem Namen (IUPAC- und Trivialnamen) oder aus SMILES erzeugen
- Summenformel und IUPAC-Name des gezeichneten Moleküls
- SMILES mit Stereoinformation kopieren

**Reaktionen beim Annähern zweier Moleküle**
- Glykosidische Bindung (Monosaccharide)
- Peptidbindung, mit dem N-Terminus links und dem C-Terminus rechts
- Polykondensation: Ester und Amid (–H₂O), Säurechlorid bzw. Phosgen (–HCl)
- Polyaddition: Urethan und Harnstoff
- Polymerisation: C=C an C=C bzw. an ein Radikal-Kettenende (Kopf-Schwanz), außerdem Kettenabbruch

**3D-Modell**
- Kugel-Stab-, Kalotten- und Stabmodell, drehbar
- Elektronendichteverteilung (elektrostatisches Potential) und Bindungspolarität mit δ+/δ− sowie Dipolmoment (Schätzung)

**Beispiele**
- Alkane, Alkohole, Carbonsäuren, Aromaten, Zucker, Aminosäuren, Peptide, Fette, Steroide
- Kunststoffe wie PE, PVC, PET, Nylon, Polycarbonat und Polyurethan

**Export**
- PNG und SVG mit frei wählbarem Dateinamen

## Aufbau des Quellcodes

| Datei | Inhalt |
|---|---|
| `src/molekuelzeichner.html` | Oberfläche, Zeichnen, Darstellungen, Reaktionen, 3D-Ansicht |
| `src/chem.js` | Chemie-Kern: SMILES-Parser, Name→Struktur, 2D-Layout, IUPAC-Namensgebung, Stereochemie |
| `build.py` | fügt beide Teile zur einzelnen Datei `index.html` zusammen |
| `index.html` | fertige, eigenständige App |

Nach Änderungen in `src/` die App neu erzeugen:

```bash
python3 build.py
```

## Hinweise

- Die Schaltfläche „IUPAC-Namen von Claude ermitteln“ steht nur zur Verfügung, wenn die App als Claude-Artefakt läuft. In der eigenständigen Version ist sie ausgeblendet. Alles andere funktioniert offline.
- Für Peptide, viele Heterocyclen und große verknüpfte Polymere berechnet die App keinen systematischen Namen.
- Das Dipolmoment und die Elektronendichte im 3D-Modell sind Näherungen (Gasteiger-Ladungen) und dienen der Veranschaulichung.
