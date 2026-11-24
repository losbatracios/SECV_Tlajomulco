# Landsat MSS — Vaso de La Playa (Tlajomulco de Zúñiga)

**Producto:** Landsat Collection 2, Multispectral Scanner (MSS)  
**Fuente:** USGS / Google Earth Engine  
**Fecha de acceso al script:** 2026-10-31  
**Resolución nominal:** 60 m por píxel (resampled en Collection 2)  
**Colecciones consultadas:** LM01/C02/T1, LM01/C02/T2, LM02/C02/T1, LM02/C02/T2,  
  LM03/C02/T1, LM03/C02/T2, LM04/C02/T1, LM04/C02/T2, LM05/C02/T1, LM05/C02/T2  
**Periodo de búsqueda:** 1972-07-01 a 1985-12-31  
**Área de interés (bbox):** 20.4647–20.4796 N; −103.5442 a −103.5100 O (WGS84)  
**Polígono de referencia:** `vaso_1971.gpkg` (~416 ha; zona de tono oscuro con márgenes
  tónicos, dibujada sobre fotografía aérea CETENAL de febrero de 1971, no equivale
  a lámina de agua abierta)

---

## Criterio de escena "utilizable" — fijado antes de ver cualquier imagen

**Una escena se clasifica como "utilizable" si y solo si:**

> La banda NIR cubre la caja del vaso sin nubes ni bandeado visible
> en al menos el 80 % del área de la caja.

El juicio visual se aplica al chip NIR (banda B7 para LM01–03; B4 para LM04–05)
exportado por el script GEE, recortado sobre el bbox más un buffer de 3 km,
con escala de grises uniforme (percentil 2–98 de toda la serie).

Este criterio se aplica **antes** de ver el contenido del vaso. La clasificación
de "agua visible" se realiza sobre las escenas "utilizables" en orden aleatorio
con semilla fija, con la fecha oculta (protocolo ciego).

---

## Criterio de "agua visible"

Para la clasificación ciega, una escena se registra como **"agua visible"** si se
observa un parche oscuro en NIR de al menos ~1 km² dentro del bbox del vaso.

- "agua no visible" significa únicamente que la resolución de 60 m no registra
  ese rasgo en esa fecha. **No equivale a "no había humedal".**
- Un vaso de baja profundidad, cubierto de vegetación acuática o con suelo
  saturado puede ser invisible en NIR a 60 m incluso cuando existe.

---

## Umbral de viabilidad del análisis (b1c.7)

Si el inventario registra **menos de 3 años distintos con al menos una escena
utilizable** en el periodo 1972–1984, el análisis MSS se declara **fuera de
alcance** y se reporta como límite metodológico en el capítulo.

---

## Archivos en este directorio

| Archivo | Descripción |
|---|---|
| `README_Landsat_MSS.md` | Este documento |
| `mss_inventory.csv` | Inventario de escenas (generado desde GEE) |
| `mss_clasificacion_plantilla.csv` | Plantilla de clasificación ciega |
| `chips/` | Chips NIR por escena (exportados desde GEE a Drive) |

## Script de análisis

| Archivo | Ruta |
|---|---|
| `mss_vaso_v1.js` | `04_Análisis/Scripts_Python/07_regresion_humedal/mss_vaso_v1.js` |
| `mss_contacto_v1.py` | `04_Análisis/Scripts_Python/07_regresion_humedal/mss_contacto_v1.py` |

---

## Nota de cadena de validación

Todas las cifras de este análisis son **exploratorias** hasta que [4. TJ Python]
las tabule en ronda 6. No citar ninguna cifra de este directorio sin esa validación.

---

*Generado: 2026-10-31 · [3. TJ Datos] · despacho b1c.7*
