#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Distancia entre la planta registrada en el DENUE y la lámina de agua de 1971
=============================================================================

PROPÓSITO
---------
Medir, de forma reproducible, la distancia entre el establecimiento que el
DENUE registra en el municipio y la lámina de agua que la fotografía aérea de
1971 registra en el vaso de La Playa.

El resultado es una distancia, nada más. No se afirma ninguna relación causal
entre la planta y la desaparición de la lámina de agua: la fotografía de 1971
registra una lámina de agua y las imágenes actuales no la registran, y eso es
lo único que sostienen estos datos. Además, la lámina de agua visible en una
fotografía no es lo mismo que el humedal entendido como cuenca, que puede ser
mayor; de ahí que la distancia se dé al borde de la lámina registrada y se
diga así.

GEOMETRÍA DE REFERENCIA
-----------------------
  planta      coordenadas del DENUE, tomadas de la ficha del establecimiento.
  centro      centro de la lámina de agua de 1971 estimado en la sesión del
              25-08 sobre las fotografías del conjunto C.
  caja        envolvente de la lámina de agua de 1971, en coordenadas
              geográficas, con el error declarado en esa misma sesión.

La caja está dada por paralelos y meridianos, así que sus lados no son rectos
en una proyección métrica. Antes de proyectar se densifican los cuatro lados,
para que la distancia se mida contra el borde real de la caja y no contra la
cuerda entre sus esquinas.

Distancias en EPSG:32613 (WGS 84 / UTM 13N).

Requiere: numpy, pyproj, shapely, pyshp
"""

import sys
from pathlib import Path

import numpy as np
from shapely.geometry import Point, Polygon

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI))
from poligonos_v12 import cargar_oficial, A_UTM

# Ficha del DENUE: LABORATORIOS PISA, CLEE 14097325412000077000005803S3,
# ID 6336993, clase 325412 (fabricación de preparaciones farmacéuticas),
# Carretera libre a San Isidro Mazatepec 7000, Tlajomulco de Zúñiga;
# localidad geoestadística "El Agostadero [PISA Farmacéutica]".
# Ficha en 05_Datos_abiertos/DENUE/PISA_DENUE_6336993.pdf
PLANTA = (20.48777778, -103.55833333)
PLANTA_ETIQUETA = "establecimiento del DENUE (ID 6336993)"

# Lámina de agua de 1971 en el vaso de La Playa, estimada el 25-08.
CENTRO_VASO = (20.4723, -103.5283)
CAJA_VASO = dict(sur=20.4647, norte=20.4796, oeste=-103.5442, este=-103.5100)
ERROR_CAJA = dict(latitud=330.0, longitud=115.0)     # metros, declarado el 25-08

# Cifras exploratorias calculadas en [1.2] el 18-10, para contrastar.
EXPLORATORIO = dict(centro=3571.0, caja=1714.0)

DENSIFICAR = 200          # puntos por lado de la caja


def caja_proyectada():
    c = CAJA_VASO
    lon = np.linspace(c["oeste"], c["este"], DENSIFICAR)
    lat = np.linspace(c["sur"], c["norte"], DENSIFICAR)
    borde = (
        [(x, c["sur"]) for x in lon]
        + [(c["este"], y) for y in lat]
        + [(x, c["norte"]) for x in lon[::-1]]
        + [(c["oeste"], y) for y in lat[::-1]]
    )
    e, n = A_UTM.transform([p[0] for p in borde], [p[1] for p in borde])
    return Polygon(zip(e, n))


def punto(lat, lon):
    return Point(*A_UTM.transform(lon, lat))


if __name__ == "__main__":
    print("=" * 88)
    print("DISTANCIA ENTRE LA PLANTA DEL DENUE Y LA LÁMINA DE AGUA DE 1971")
    print("=" * 88)

    p = punto(*PLANTA)
    centro = punto(*CENTRO_VASO)
    caja = caja_proyectada()
    agem = cargar_oficial(["14097"])["14097"]

    print(f"\n  {PLANTA_ETIQUETA}")
    print(f"    coordenadas          {PLANTA[0]:.8f}, {PLANTA[1]:.8f}")
    print(f"    en UTM 13N           E {p.x:,.0f}  N {p.y:,.0f}")
    dentro = agem.contains(p)
    print(f"    ¿dentro del AGEM 14 097?   "
          f"{'sí' if dentro else 'no'}"
          f"   (distancia al borde {agem.exterior.distance(p):,.0f} m"
          f"{', por dentro' if dentro else ''})")

    d_centro = p.distance(centro)
    d_caja = p.distance(caja)
    print(f"\n  {'medida':52}{'metros':>10}{'exploratorio [1.2]':>22}")
    print("  " + "─" * 84)
    print(f"  {'al centro estimado de la lámina de agua de 1971':52}"
          f"{d_centro:>10,.0f}{EXPLORATORIO['centro']:>18,.0f} m"
          f"{'  ✓' if abs(d_centro - EXPLORATORIO['centro']) < 5 else '  ←'}")
    print(f"  {'al borde de la caja de la lámina de agua de 1971':52}"
          f"{d_caja:>10,.0f}{EXPLORATORIO['caja']:>18,.0f} m"
          f"{'  ✓' if abs(d_caja - EXPLORATORIO['caja']) < 5 else '  ←'}")

    # Descomposición: la planta queda al oeste y al norte de la esquina
    # noroeste de la caja, así que el punto más próximo de la caja es esa
    # esquina y la distancia es la diagonal entre las dos componentes.
    noroeste = punto(CAJA_VASO["norte"], CAJA_VASO["oeste"])
    print(f"\n  respecto de la esquina noroeste de la caja")
    print(f"    componente este–oeste   {abs(p.x - noroeste.x):>8,.0f} m"
          f"   ({'al oeste' if p.x < noroeste.x else 'al este'})")
    print(f"    componente norte–sur    {abs(p.y - noroeste.y):>8,.0f} m"
          f"   ({'al norte' if p.y > noroeste.y else 'al sur'})")
    print(f"    distancia a la esquina  {p.distance(noroeste):>8,.0f} m")

    print(f"\n  error declarado de la caja (25-08): {ERROR_CAJA['latitud']:,.0f} m en"
          f" latitud · {ERROR_CAJA['longitud']:,.0f} m en longitud.")
    holgura = d_caja / max(ERROR_CAJA.values())
    print(f"  La distancia al borde es {holgura:.1f} veces el mayor de esos errores.")

    print(f"\n  Formulación admitida para el texto:")
    print(f"    el establecimiento que el DENUE registra se ubica a unos"
          f" {d_caja/1000:.1f} km al oeste-noroeste")
    print(f"    del borde de la lámina de agua que registra la fotografía de 1971.")
    print(f"    No se afirma ninguna relación causal, y la lámina de agua no es")
    print(f"    el humedal entendido como cuenca.")

    print("\n" + "=" * 88)
    print("FIN DEL SCRIPT")
    print("=" * 88)
