/**
 * mss_vaso_v1.js — Landsat MSS sobre el vaso de La Playa
 * Proyecto: Tlajomulco · [3. TJ Datos] · despacho b1c.7
 * Fecha de creación: 2026-10-31
 *
 * INSTRUCCIONES PARA EMILIO — leer antes de Run:
 * ─────────────────────────────────────────────
 * PASO 1. Pegar TODO este contenido en Google Earth Engine
 *         (code.earthengine.google.com → botón "+" → New Script)
 * PASO 2. Hacer clic en "Run"
 * PASO 3. Panel Console → verificar el conteo de escenas encontradas
 * PASO 4. Panel Tasks → hacer "Run" en:
 *         a) "mss_inventory_vaso_la_playa"  → CSV de inventario
 *         b) Una tarea "mss_chip_YYYY" por cada año  → chips NIR anuales (mosaico)
 * PASO 5. Cuando termine, copiar desde Google Drive (carpeta Landsat_MSS_Tlajomulco):
 *         - mss_inventory.csv   →  05_Datos_abiertos/Landsat_MSS/mss_inventory.csv
 *         - mss_chip_YYYY.tif   →  05_Datos_abiertos/Landsat_MSS/chips/
 * PASO 6. Ejecutar mss_contacto_v1.py localmente para armar las hojas de contacto
 *
 * BANDAS NIR por satélite:
 *   LM01, LM02, LM03  →  B7  (NIR2, 0.80–1.10 µm)
 *   LM04, LM05         →  B4  (NIR2, 0.80–1.10 µm)
 *
 * CRITERIO DE USABILIDAD (fijado en README antes de ver imágenes):
 *   "caja sin nubes ni bandeado visible en ≥80 % del área"
 *   Se aplica sobre los chips al clasificar, no en este script.
 */

// =============================================================================
// 1. ÁREA DE INTERÉS
// =============================================================================

var bbox = ee.Geometry.Rectangle(
  [-103.5442, 20.4647, -103.5100, 20.4796],
  'EPSG:4326',
  false
);

// Buffer de 3 km para los chips de visualización
var bboxBuffer = bbox.buffer(3000);

Map.centerObject(bbox, 13);
Map.addLayer(bbox,       {color: 'red'},  'Bbox vaso La Playa');
Map.addLayer(bboxBuffer, {color: 'blue'}, 'Buffer 3 km');

// =============================================================================
// 2. CONSTRUIR COLECCIÓN UNIFICADA LM01–LM05 / T1+T2
// =============================================================================

var START_DATE = '1972-07-01';
var END_DATE   = '1985-12-31';

/**
 * Para LM01–03: la banda NIR2 se llama B7
 * Para LM04–05: la banda NIR2 se llama B4
 * Renombramos todo a 'NIR' para poder unir las colecciones.
 */
function loadAndNormalize(collId, nirBand, sat, tier) {
  return ee.ImageCollection(collId)
    .filterDate(START_DATE, END_DATE)
    .filterBounds(bbox)
    .map(function(img) {
      return img
        .select([nirBand], ['NIR'])
        .set({
          'satellite': sat,
          'tier':      tier,
          'date_str':  img.date().format('YYYY-MM-dd'),
          'year':      img.date().get('year'),
          'month':     img.date().get('month'),
          'cloud_cover_pct': ee.Number(img.get('CLOUD_COVER')).round(),
          'scene_id':  img.id(),
          'nir_band':  nirBand,
          'collection': collId
        });
    });
}

var allMSS = ee.ImageCollection([])
  .merge(loadAndNormalize('LANDSAT/LM01/C02/T1', 'B7', 'LM01', 'T1'))
  .merge(loadAndNormalize('LANDSAT/LM01/C02/T2', 'B7', 'LM01', 'T2'))
  .merge(loadAndNormalize('LANDSAT/LM02/C02/T1', 'B7', 'LM02', 'T1'))
  .merge(loadAndNormalize('LANDSAT/LM02/C02/T2', 'B7', 'LM02', 'T2'))
  .merge(loadAndNormalize('LANDSAT/LM03/C02/T1', 'B7', 'LM03', 'T1'))
  .merge(loadAndNormalize('LANDSAT/LM03/C02/T2', 'B7', 'LM03', 'T2'))
  .merge(loadAndNormalize('LANDSAT/LM04/C02/T1', 'B4', 'LM04', 'T1'))
  .merge(loadAndNormalize('LANDSAT/LM04/C02/T2', 'B4', 'LM04', 'T2'))
  .merge(loadAndNormalize('LANDSAT/LM05/C02/T1', 'B4', 'LM05', 'T1'))
  .merge(loadAndNormalize('LANDSAT/LM05/C02/T2', 'B4', 'LM05', 'T2'));

print('Total de escenas encontradas en bbox 1972–1985:', allMSS.size());
print('Primeras 20 escenas (ordenadas por fecha):', allMSS.sort('date_str').limit(20));

// =============================================================================
// 3. EXPORTAR INVENTARIO CSV
// =============================================================================

var inventoryFC = allMSS.map(function(img) {
  return ee.Feature(null, {
    'scene_id':        img.get('scene_id'),
    'satellite':       img.get('satellite'),
    'tier':            img.get('tier'),
    'date':            img.get('date_str'),
    'year':            img.get('year'),
    'month':           img.get('month'),
    'cloud_cover_pct': img.get('cloud_cover_pct'),
    'nir_band':        img.get('nir_band'),
    'collection':      img.get('collection')
  });
});

Export.table.toDrive({
  collection:     inventoryFC,
  description:    'mss_inventory_vaso_la_playa',
  folder:         'Landsat_MSS_Tlajomulco',
  fileNamePrefix: 'mss_inventory',
  fileFormat:     'CSV',
  selectors:      ['scene_id', 'satellite', 'tier', 'date', 'year', 'month',
                   'cloud_cover_pct', 'nir_band', 'collection']
});

print('✓ Tarea de inventario creada → Tasks → "mss_inventory_vaso_la_playa" → Run');

// =============================================================================
// 4. EXPORTAR CHIPS NIR — UN MOSAICO POR AÑO (mediana)
// =============================================================================
// Un mosaico anual es manejable (13 años = 13 tareas).
// Si el inventario muestra escenas de interés en meses específicos,
// se puede refinar con mss_vaso_v2.js (escenas individuales).

// Calcular estiramiento global para escala uniforme entre chips
// Usamos percentil 2–98 sobre toda la colección recortada al bbox
var nirStack = allMSS.select('NIR').median().clip(bboxBuffer);
var globalStats = nirStack.reduceRegion({
  reducer:  ee.Reducer.percentile([2, 98]),
  geometry: bboxBuffer,
  scale:    60,
  maxPixels: 1e8
});
var vizMin = ee.Number(globalStats.get('NIR_p2'));
var vizMax = ee.Number(globalStats.get('NIR_p98'));

print('Estiramiento global NIR (p2–p98):', vizMin, '–', vizMax);

// Mostrar mosaico en el mapa
Map.addLayer(
  nirStack,
  {bands: ['NIR'], min: vizMin, max: vizMax, palette: ['000000', 'ffffff']},
  'MSS NIR mediana 1972–1985'
);

// Exportar un chip por año
var years = ee.List.sequence(1972, 1985);

years.getInfo().forEach(function(yr) {
  var yearStr = yr.toString();
  var annual = allMSS.filter(ee.Filter.eq('year', yr));
  var count = annual.size().getInfo();

  if (count === 0) {
    print('  Sin escenas en ' + yearStr + ' — no se crea tarea de chip.');
    return;
  }

  print('  ' + yearStr + ': ' + count + ' escenas → creando tarea de chip...');

  var mosaic = annual.select('NIR').median().clip(bboxBuffer);

  Export.image.toDrive({
    image:          mosaic.visualize({
                      min: vizMin.getInfo(),
                      max: vizMax.getInfo(),
                      palette: ['000000', 'ffffff']
                    }),
    description:    'mss_chip_' + yearStr,
    folder:         'Landsat_MSS_Tlajomulco',
    fileNamePrefix: 'mss_chip_' + yearStr,
    region:         bboxBuffer,
    scale:          60,
    crs:            'EPSG:4326',
    maxPixels:      1e8
  });
});

print('✓ Tareas de chip anuales creadas en Tasks → hacer Run en cada una.');
print('  Después de descargar, ejecutar mss_contacto_v1.py para armar hojas de contacto.');

// =============================================================================
// 5. NOTA METODOLÓGICA (en consola)
// =============================================================================

print('');
print('═══════════════════════════════════════════════════════');
print('NOTA: el chip anual es mediana de todas las escenas del año.');
print('Si en un año hay escenas con y sin nubes, la mediana puede');
print('enmascarar parcialmente las nubes, pero no eliminarlas.');
print('La clasificación ciega se aplica igualmente: el chip es');
print('"no utilizable" si hay cobertura de nubes visible ≥20 %.');
print('La fecha real de cada escena queda en mss_inventory.csv.');
print('═══════════════════════════════════════════════════════');
