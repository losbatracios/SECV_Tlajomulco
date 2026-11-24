# SECV_Tlajomulco

[![DOI](https://zenodo.org/badge/1457111180.svg)](https://doi.org/10.5281/zenodo.24208442)

Code and derived data for the book chapter **"Los límites de la mirada. Archivo fotogramétrico y transformación socioambiental"** (Ana Paula Montes Ruiz and Emilio José Bermeo Donoso), forthcoming in *El punto ciego*, Seminario Permanente de Estudios de Cultura Visual, CEIICH-UNAM.

The chapter reads INEGI's aerial photography archive of Tlajomulco de Zúñiga (Jalisco, Mexico) through two operations:

1. **Archive coverage.** Counts, densities and flight-line geometry from the flight index of INEGI's Sistema de Información de Índices de Vuelos (SIIV), queried on 7 July 2026 (2,074 records, six central municipalities of the Guadalajara metropolitan area).
2. **Cartographic regression.** Georeferencing of three aerial photographs (February 1971 and February 1993 over the airport corridor; February 1971 over the Presa La Playa basin), with prediction errors measured against the reference base (Esri World Imagery Wayback, version of 5 August 2026, release 26334), plus two satellite series over the basin (JRC Global Surface Water, Landsat MSS).

Every figure the chapter cites is listed, with its exact name, sample size, producing script and caveat, in [`04_Análisis/cifras_validadas_v2.csv`](04_Análisis/cifras_validadas_v2.csv). A methodological note in Spanish with the detail that did not fit in the chapter is in [`docs/nota_metodo.md`](docs/nota_metodo.md).

## Reproducing the figures

```bash
pip install -r requirements.txt
cd "04_Análisis"
python cifras_validadas_v2b.py
```

This regenerates `cifras_validadas_v2.csv` from the files in this repository. Tested on 20 November 2026 with Python 3.13 and the versions in `requirements.txt`: the regenerated file is byte-identical to the one committed here. The scripts it calls read their inputs with paths relative to the repository root, so the folder structure must be kept as is.

Each analysis script can also be run on its own from its folder (for example `python linea_base_v8.py` inside `04_Análisis/Scripts_Python/06_linea_base/`). Run from this repository, all the scripts listed under "Reproducible chain" below rewrite their tables (CSV, GeoPackage) identically. Two kinds of output can differ slightly: the transformation matrices in `_validacion_v3/matrices_pixel_a_UTM32613.json`, whose coefficients change by less than 10⁻⁵ (sub-millimetre on the ground) because of floating-point differences between library versions, and the rendered figures, at pixel level. No cited figure changes.

## Repository structure

The structure mirrors the working folder of the project, so the scripts run without edits.

```
02_Archivo_INEGI/
  Metadatos_SINFA/            SIIV flight index, query of 7 July 2026 (CSV, latin-1) and the query page (PDF)
  Imágenes/_pares/            metadata only: master table of the 98 requested frames and the copy manifest
04_Análisis/
  cifras_validadas_v2.csv     the validated figures (one row per figure, with exact name, n, script, caveat, status)
  cifras_validadas_v1.py, _v2.py, _v2b.py   generators of that table
  Scripts_Python/
    02b_…                     corpus against the official municipal area (AGEM 14 097)
    05_regresion_par_estrella/   airport-corridor pair: models, prediction error, scale; GCP and checkpoint files (_gcp_qgis/)
    06_linea_base/            municipal polygons, coverage baseline, flight lines, nominal footprints, DENUE plant
    07_regresion_humedal/     Presa La Playa photograph: fit, contour, archive coverage of the basin, JRC; GCP files
    08_mss/                   Landsat MSS contrast over the basin
    08_figuras_capitulo/      figure scripts (figures 1–3 and 3b run from this repository)
  Visualizaciones/            outputs of those scripts (CSV, JSON, GeoPackage) and figures 2, 3 and 3b
05_Datos_abiertos/
  Capas_cartograficas/        INEGI Marco Geoestadístico 2025, municipal layer for Jalisco (14mun) and derived layers
  JRC_GSW/                    monthly JRC Global Surface Water series over the basin box
  Landsat_MSS/                MSS inventory, annual chips (1972–1985) and the classification template
docs/nota_metodo.md           methodological note (Spanish)
```

### Reproducible chain

`02b_validacion_corpus_AGEM_oficial.py`, `06_linea_base/linea_base_v8.py`, `huella_v1.py`, `cifras_municipales_v12.py`, `pisa_vaso_v7.py`, `franja_norte_v12b.py`, `05_regresion_par_estrella/validacion_par_estrella_v3.py`, `escala_1971_v5d.py`, `07_regresion_humedal/cobertura_vaso_v1.py`, `08_mss/mss_contraste_v1.py`, `08_figuras_capitulo/figura_3b_huellas_v1.py`, and the three `cifras_validadas` generators.

### Steps that need the INEGI photographs (not included)

`07_regresion_humedal/validacion_humedal_v1.py` reads the georeferenced 1971 photograph of the basin; `05_regresion_par_estrella/escala_v5b_borde_fotograma.py` and `fiduciales_v5c.py` read the scanned frames; `03_inspeccion_calidad_corpus.py` reads all 98 frames; figures 4 and 5 of the chapter use the photographs and Esri imagery. Their outputs used by the chapter are included (`Visualizaciones/05_humedal_georef/`, `01_inspeccion_calidad/calidad_corpus.csv`), and the GCP and checkpoint files allow the fits to be recomputed by anyone who obtains the same frames from INEGI (frames 0237BZ L56 F41, 0713VE L1 F29 and 0249BZ L37 F5; see the master table).

### Superseded scripts kept for the record

- `07_regresion_humedal/jrc_serie_v1.py`: first reading of the JRC series. It does not run with current pandas and its root path assumed another folder layout. The citable series is the one recomputed by `cobertura_vaso_v1.py` (`Visualizaciones/05_humedal_georef/_validacion/jrc_replica_mensual.csv`).
- `07_regresion_humedal/mss_contacto_v1.py`: contact sheets for a visual classification of 132 MSS scenes. The material turned out to be 12 annual mosaics, so the 132 rows were 12 images seen several times; the citable reading is `08_mss/mss_contraste_v1.py`. The note explains this.
- `jrc_serie_v1.js` and `mss_vaso_v1.js` are the Google Earth Engine scripts that produced the JRC series and the MSS chips.

## Notes on the data

- **Requested sets.** Set A has 12 frames from 1993; set B has the 19 pre-reform frames that were requested from INEGI. The request selected B with a spherical haversine distance under 2 km; measured on the WGS84 ellipsoid, the same criterion yields 20 (frame 104AZ/4/44/1 lies at 1,998 m). The chapter reports the 19 that were requested and analysed. Set C has 67 frames of 1971.
- **Centre versus footprint.** A frame is counted in a municipality when its declared centre falls inside the municipal area (AGEM 14 097 of the Marco Geoestadístico 2025, which is a statistical area, not an official political boundary). Counting by nominal footprint instead changes some counts (see `huella_v1.py` and the note) but not the conclusions.
- **Error names.** Georeferencing errors are always reported with their exact name: RMS of residuals, corrected standard error, leave-one-out prediction error, or prediction error at independent checkpoints. All are relative to the reference base, not absolute accuracy.

## Not included

- The INEGI aerial photographs, which INEGI distributes under its own terms.
- Esri imagery, including the Wayback crops used for the chapter figures.
- Any material from the Colectiva Agua y Vida.
- The text of the chapter, until it is published.

## Licenses

- Code: MIT (see `LICENSE`).
- Data produced by the authors: CC BY 4.0 (see `LICENSE-DATA`).
- Third-party data remain under their original terms: INEGI (SIIV index, Marco Geoestadístico, DENUE; Términos de libre uso de la información del INEGI), JRC Global Surface Water (European Commission), USGS Landsat.

## Funding

Emilio José Bermeo Donoso carried out this work as a doctoral fellow of the Secretaría de Ciencia, Humanidades, Tecnología e Innovación (SECIHTI), Mexico.

## Authors

- Emilio José Bermeo Donoso, Universidad Iberoamericana. ORCID [0000-0002-1296-6713](https://orcid.org/0000-0002-1296-6713)
- Ana Paula Montes Ruiz, Instituto de Geografía, UNAM. ORCID [0000-0001-6557-7765](https://orcid.org/0000-0001-6557-7765)

## How to cite

Please cite the Zenodo record, DOI [10.5281/zenodo.24208442](https://doi.org/10.5281/zenodo.24208442) (this DOI always resolves to the latest version), and the chapter.
