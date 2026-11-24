# JRC Global Surface Water — Vaso de La Playa, Tlajomulco de Zúñiga

## Producto

**JRC Global Surface Water, Monthly Water History**  
Versión: **v1.4** (`JRC/GSW1_4/MonthlyHistory`)  
Periodo cubierto: **1984-03 a 2021-12** (38 años, hasta 456 registros mensuales posibles)  
Resolución: **30 m/px** (Landsat 5/7/8)  
Referencia: Pekel, J.-F., Cottam, A., Gorelick, N. & Belward, A. S. (2016). High-resolution mapping of global surface water and its long-term changes. *Nature*, 540, 418–422. https://doi.org/10.1038/nature20584  
Acceso vía Google Earth Engine (GEE): https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_MonthlyHistory  
Fecha de preparación del script: **2026-10-26**

## Área de análisis

**Caja del vaso de 1971** (WGS84):  
Latitud: 20.4647 – 20.4796 N  
Longitud: −103.5442 – −103.5100 O  
Área de la caja (ref. 1971): 6,538.2 píxeles de 30 m = **588 ha** (corregido por [4] el 31-10; antes decía ~568 ha). El conteo sale del propio CSV crudo: px_water + px_no_water + px_no_obs suma 6,538.2 en cada mes.

El polígono del vaso (04_Análisis/Visualizaciones/05_humedal_georef/vaso_1971.gpkg)
existe pero **no está validado** a la fecha de este script. Se usa la caja.
Cuando [4. TJ Python] entregue el polígono validado (ronda 6), actualizar el script
GEE reemplazando la variable `bbox` por ese polígono.

## Banda 'water' — valores

| Valor | Significado | Interpretación |
|-------|-------------|----------------|
| 0 | Sin observación | No hay imagen Landsat válida ese mes — **declarar como sin obs, no como sin agua** |
| 1 | Sin agua | El pixel está seco o cubierto de tierra ese mes |
| 2 | Con agua | El pixel está inundado ese mes |

## Archivos

```
05_Datos_abiertos/JRC_GSW/
├── README_JRC_GSW.md                  ← este archivo
├── jrc_gsw_monthly_bbox_v1.csv        ← CSV bruto exportado por GEE (PENDIENTE: descargar de Drive)
├── jrc_monthly_clean.csv              ← generado por jrc_serie_v1.py
├── jrc_annual_summary.csv             ← generado por jrc_serie_v1.py
├── jrc_february_series.csv            ← generado por jrc_serie_v1.py
└── jrc_agua_por_anio.png              ← gráfica, generada por jrc_serie_v1.py

04_Análisis/Scripts_Python/07_regresion_humedal/
├── jrc_serie_v1.js                    ← script GEE (pegar en Code Editor)
└── jrc_serie_v1.py                    ← procesamiento local del CSV
```

## Flujo de trabajo

1. **Ejecutar jrc_serie_v1.js en GEE** (https://code.earthengine.google.com/)  
   Exporta `jrc_gsw_monthly_bbox_v1.csv` a Google Drive, carpeta `JRC_GSW_Tlajomulco`.  
   Tiempo estimado: 2–5 min en la cola de GEE.

2. **Descargar el CSV** desde Google Drive y colocarlo en `05_Datos_abiertos/JRC_GSW/`.

3. **Ejecutar jrc_serie_v1.py** desde la raíz del proyecto:  
   ```
   python 04_Análisis/Scripts_Python/07_regresion_humedal/jrc_serie_v1.py
   ```  
   Genera los tres CSV de salida y la gráfica.

4. **Enviar a [4. TJ Python]** para replicación (cifras no citables hasta entonces).

## Notas metodológicas

- **Cronología, no causalidad.** Este análisis describe *en qué meses y años se registra
  agua* dentro de la ventana temporal del satélite (1984–2021). No establece causas.
  Verbos admitidos: *registra, no registra, coincide con, precede a, es posterior a*.
  Vetados: *secó, desecó, habilitó, acompañó, facilita*.

- **Meses sin observación.** Los registros con `px_no_obs` alto (valor 0 en la banda
  'water') son meses sin cobertura Landsat válida. Se declaran como *sin observación*,
  nunca como *sin agua*.

- **Serie de febrero.** Febrero es el mes de la fotografía aérea de 1971
  (249BZ_37_5_1971). La comparación directa 1971–satélite se limita a ese mes.

- **Área de referencia.** El denominador de `pct_bbox_area` es el área total de la caja
  de 1971 calculada en GEE. Cuando el polígono validado esté disponible, actualizar.

- **Cifras exploratorias.** Todos los valores de este conjunto son exploratorios hasta
  que [4. TJ Python] los replique desde los archivos fuente.

- **Versión v1.4 vs v1.5.** Si GEE ofrece JRC/GSW1_5/MonthlyHistory al ejecutar el
  script, verificar el ID antes de correr. La v1.4 cubre hasta 2021-12.
