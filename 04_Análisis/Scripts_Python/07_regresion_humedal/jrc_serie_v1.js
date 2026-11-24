/**
 * jrc_serie_v1.js
 * ---------------
 * Serie temporal JRC Global Surface Water — vaso de La Playa, Tlajomulco de Zúñiga
 *
 * Producto  : JRC Global Surface Water, Monthly Water History v1.4
 * Fuente    : Pekel et al. (2016), Nature 540, 418–422. doi:10.1038/nature20584
 * Colección : JRC/GSW1_4/MonthlyHistory
 * Periodo   : 1984-03 a 2021-12 (v1.4; actualizar ID si se usa v1.5+)
 * Resolución: 30 m/px (Landsat)
 * Banda     : 'water'  → 0 = sin observación, 1 = sin agua, 2 = con agua
 *
 * Área      : caja del vaso de 1971
 *   Latitud : 20.4647 – 20.4796 N
 *   Longitud: −103.5442 – −103.5100 O  (WGS84)
 *
 * Salida    : CSV exportado a Google Drive, carpeta 'JRC_GSW_Tlajomulco'
 *   Columnas: year, month, px_water, px_no_water, px_no_obs, area_ha_water,
 *             pct_bbox_area (% del área de la caja de 1971 con agua)
 *
 * Uso       : 1. Abrir https://code.earthengine.google.com/
 *             2. Crear nuevo script, pegar este código completo.
 *             3. Clic en "Run". La exportación aparece en la pestaña "Tasks".
 *             4. Clic en "Run" junto a la tarea para confirmar la exportación.
 *             5. El CSV llega a Google Drive > JRC_GSW_Tlajomulco/
 *             6. Mover a 05_Datos_abiertos/JRC_GSW/ y ejecutar jrc_serie_v1.py
 *                para generar estadísticas, resumen y gráfica.
 *
 * Tiempo estimado de ejecución: 2–5 min (depende de la cola de GEE).
 * Cifras son exploratorias hasta que [4. TJ Python] las replique.
 *
 * Fecha de preparación del script: 2026-10-26
 */

// ─── 1. Área de análisis ──────────────────────────────────────────────────
var bbox = ee.Geometry.Rectangle(
  [-103.5442, 20.4647, -103.5100, 20.4796],  // [xMin, yMin, xMax, yMax]
  'EPSG:4326',
  false  // no geodésica; caja rectangular en WGS84
);

// Área de referencia de 1971 (la caja completa, en ha)
var bboxArea_ha = bbox.area(1, 'EPSG:32613').divide(10000);

// ─── 2. Colección JRC Monthly Water History ──────────────────────────────
var gsw = ee.ImageCollection('JRC/GSW1_4/MonthlyHistory')
  .filterBounds(bbox);

// ─── 3. Función: estadísticas por imagen (un mes) ────────────────────────
function monthStats(img) {
  var date  = ee.Date(img.get('system:time_start'));
  var year  = date.get('year');
  var month = date.get('month');

  var clipped = img.clip(bbox);

  var water    = clipped.eq(2);
  var no_water = clipped.eq(1);
  var no_obs   = clipped.eq(0);

  var scale = 30;
  var proj  = 'EPSG:32613';

  var countWater = water.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: bbox, scale: scale,
    crs: proj, maxPixels: 1e9
  }).get('water');

  var countNoWater = no_water.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: bbox, scale: scale,
    crs: proj, maxPixels: 1e9
  }).get('water');

  var countNoObs = no_obs.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: bbox, scale: scale,
    crs: proj, maxPixels: 1e9
  }).get('water');

  var areaWater_ha = ee.Number(countWater).multiply(scale * scale).divide(10000);
  var pct = areaWater_ha.divide(bboxArea_ha).multiply(100);

  return ee.Feature(null, {
    'year'         : year,
    'month'        : month,
    'px_water'     : countWater,
    'px_no_water'  : countNoWater,
    'px_no_obs'    : countNoObs,
    'area_ha_water': areaWater_ha,
    'pct_bbox_area': pct
  });
}

// ─── 4. Aplicar a toda la colección ──────────────────────────────────────
var stats = gsw.map(monthStats);
var fc    = ee.FeatureCollection(stats);

// ─── 5. Imprimir primeras filas ───────────────────────────────────────────
print('Primeras 5 filas (verificación):', fc.limit(5));
print('Total de registros:', fc.size());
print('Área de la caja (ha, ref. 1971):', bboxArea_ha);

// ─── 6. Exportar a Drive ──────────────────────────────────────────────────
Export.table.toDrive({
  collection  : fc,
  description : 'jrc_gsw_vaso_la_playa_1984_2021',
  folder      : 'JRC_GSW_Tlajomulco',
  fileNamePrefix: 'jrc_gsw_monthly_bbox_v1',
  fileFormat  : 'CSV'
});

// ─── Notas metodológicas ──────────────────────────────────────────────────
// • Los meses con px_no_obs alto no se interpretan como 'sin agua'.
// • El área de referencia es la CAJA de 1971, no el polígono del vaso.
// • Verbos admitidos: registra, no registra, coincide con, precede a,
//   es posterior a. NO se infiere causa.
// • Cifras exploratorias hasta que [4. TJ Python] las replique.
// • Versión: JRC/GSW1_4/MonthlyHistory (v1.4, Pekel et al. 2016).
