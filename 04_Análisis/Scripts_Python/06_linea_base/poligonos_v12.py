#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Áreas geoestadísticas municipales: carga y comparación de las tres fuentes
==========================================================================

PROPÓSITO
---------
Cargar las áreas geoestadísticas municipales (AGEM) de los seis municipios de
la Zona Metropolitana de Guadalajara que intervienen en el análisis, desde la
fuente oficial, y compararlas con las dos capas que se venían usando.

Nota de precisión, que vale para todo el proyecto: el Marco Geoestadístico es
un instrumento de referencia para la captación y presentación de información
estadística; el propio INEGI advierte que su cartografía no representa el
establecimiento ni la definición oficial de límites político-administrativos.
Por eso aquí se habla de "área geoestadística municipal (AGEM)" y nunca de
"límite municipal oficial".

LAS TRES FUENTES
----------------
  oficial       Marco Geoestadístico 2025 del INEGI, Jalisco (UPC
                794551163061), capa municipal `14mun.shp`. CRS EPSG:6372
                (Mexico ITRF2008 / Lambert Conformal Conic), codificación
                ISO-8859-1 según el .cpg. Es la fuente canónica.
  derivado      `AGEM_14_ZMG6_MacWilliXD_20261013.geojson`: derivado de
                terceros del mismo Marco Geoestadístico 2025, fijado a un
                commit, con los seis municipios.
  provisional   `poligono_tlajomulco_14097.geojson`: la capa con la que se
                calcularon las cifras municipales hasta el 14-10. Solo
                Tlajomulco.

QUÉ SE COMPARA
--------------
Por municipio, entre la fuente oficial y cada una de las otras dos:

  área            de cada capa y su diferencia relativa;
  diferencia      área de la diferencia simétrica, es decir, la superficie
  simétrica       que una capa incluye y la otra no, en ambos sentidos, y su
                  peso sobre el área oficial;
  Hausdorff       la mayor distancia desde un punto de un contorno al punto
                  más próximo del otro. Es la cota del desplazamiento máximo
                  entre los dos trazos, y responde a la pregunta que importa
                  aquí: cuánto puede moverse un punto al cambiar de capa.

La distancia de Hausdorff que calcula la biblioteca es la versión discreta,
sobre los vértices de cada contorno, de modo que es una cota inferior de la
continua. Con contornos de miles de vértices la diferencia es despreciable.

Todas las medidas se hacen en EPSG:32613 (WGS 84 / UTM 13N), que es el CRS
métrico del proyecto.

USO COMO MÓDULO
---------------
    from poligonos_v12 import cargar_oficial, SEIS, NOMBRES
    agem = cargar_oficial()          # {cvegeo: polígono en EPSG:32613}
    tlaj = agem["14097"]

Requiere: numpy, pyproj, shapely, pyshp
"""

import json
from pathlib import Path

import shapefile                       # pyshp
from pyproj import Transformer
from shapely.geometry import shape, Polygon, MultiPolygon
from shapely.ops import transform as shp_transform, unary_union

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[2]                 # raíz del repositorio del proyecto

SHP_OFICIAL = (RAIZ / "05_Datos_abiertos/Capas_cartograficas/MG_INEGI_2025_Jalisco"
                    / "794551163061_s/14_jalisco/conjunto_de_datos/14mun")
GEOJSON_DERIVADO = (RAIZ / "05_Datos_abiertos/Capas_cartograficas"
                         / "AGEM_14_ZMG6_MacWilliXD_20261013.geojson")
GEOJSON_PROVISIONAL = (RAIZ / "05_Datos_abiertos/Capas_cartograficas"
                            / "poligono_tlajomulco_14097.geojson")

# El .cpg del Marco Geoestadístico dice "ISO 88591", que es ISO-8859-1 con un
# espacio de más; pyshp no reconoce esa grafía, así que se pasa explícita.
CODIFICACION_DBF = "latin-1"

CRS_OFICIAL = "EPSG:6372"              # Mexico ITRF2008 / LCC
CRS_GEOJSON = "EPSG:4326"              # coordenadas geográficas
CRS_TRABAJO = "EPSG:32613"             # WGS 84 / UTM 13N

SEIS = ["14039", "14070", "14097", "14098", "14101", "14120"]
NOMBRES = {
    "14039": "Guadalajara",
    "14070": "El Salto",
    "14097": "Tlajomulco de Zúñiga",
    "14098": "San Pedro Tlaquepaque",
    "14101": "Tonalá",
    "14120": "Zapopan",
}

_T_OFICIAL = Transformer.from_crs(CRS_OFICIAL, CRS_TRABAJO, always_xy=True)
_T_GEO = Transformer.from_crs(CRS_GEOJSON, CRS_TRABAJO, always_xy=True)
A_UTM = Transformer.from_crs(CRS_GEOJSON, CRS_TRABAJO, always_xy=True)
A_GEO = Transformer.from_crs(CRS_TRABAJO, CRS_GEOJSON, always_xy=True)


def _a_utm(geom, transformador):
    return shp_transform(lambda x, y, z=None: transformador.transform(x, y), geom)


def _poligono_de_shape(sh):
    """Convierte una geometría de pyshp en Polygon o MultiPolygon.

    Un registro del shapefile puede traer varios anillos: el primero de cada
    grupo es el contorno exterior y los que giran en sentido contrario son
    huecos. Se resuelve orientándose por el signo del área con la fórmula del
    polígono, que es el criterio del propio formato.
    """
    partes = list(sh.parts) + [len(sh.points)]
    anillos = [sh.points[partes[i]:partes[i + 1]] for i in range(len(sh.parts))]

    def area_firmada(p):
        s = 0.0
        for i in range(len(p) - 1):
            s += p[i][0] * p[i + 1][1] - p[i + 1][0] * p[i][1]
        return s / 2.0

    exteriores, huecos = [], []
    for a in anillos:
        (exteriores if area_firmada(a) < 0 else huecos).append(a)
    if not exteriores:                 # todos en el mismo sentido
        exteriores, huecos = anillos, []

    piezas = []
    for ext in exteriores:
        caja = Polygon(ext)
        propios = [h for h in huecos if caja.contains(Polygon(h).representative_point())]
        piezas.append(Polygon(ext, propios))
    geom = piezas[0] if len(piezas) == 1 else MultiPolygon(piezas)
    return geom if geom.is_valid else geom.buffer(0)


def cargar_oficial(claves=None):
    """AGEM del Marco Geoestadístico 2025, en EPSG:32613."""
    claves = set(claves or SEIS)
    r = shapefile.Reader(str(SHP_OFICIAL), encoding=CODIFICACION_DBF)
    salida = {}
    for sr in r.iterShapeRecords():
        d = sr.record.as_dict()
        cve = d.get("CVEGEO") or (d.get("CVE_ENT", "") + d.get("CVE_MUN", ""))
        if cve in claves:
            salida[cve] = _a_utm(_poligono_de_shape(sr.shape), _T_OFICIAL)
    r.close()
    faltan = claves - set(salida)
    if faltan:
        raise ValueError(f"no están en el shapefile: {sorted(faltan)}")
    return salida


def cargar_derivado(claves=None):
    """AGEM del derivado de terceros, en EPSG:32613."""
    claves = set(claves or SEIS)
    gj = json.load(open(GEOJSON_DERIVADO, encoding="utf-8"))
    salida = {}
    for f in gj["features"]:
        p = f["properties"]
        cve = str(p.get("cvegeo") or p.get("CVEGEO") or "")
        if cve in claves:
            g = shape(f["geometry"])
            salida[cve] = _a_utm(g if g.is_valid else g.buffer(0), _T_GEO)
    return salida


def cargar_provisional():
    """La capa provisional de Tlajomulco usada hasta el 14-10, en EPSG:32613."""
    gj = json.load(open(GEOJSON_PROVISIONAL, encoding="utf-8"))
    g = unary_union([shape(f["geometry"]) for f in gj["features"]])
    return _a_utm(g if g.is_valid else g.buffer(0), _T_GEO)


# ---------------------------------------------------------------------------
# Comparación
# ---------------------------------------------------------------------------

def comparar(a, b):
    """Área de cada capa, diferencia simétrica y distancia de Hausdorff."""
    dif = a.symmetric_difference(b)
    return dict(
        area_a=a.area / 1e6,
        area_b=b.area / 1e6,
        solo_a=a.difference(b).area / 1e6,
        solo_b=b.difference(a).area / 1e6,
        sim=dif.area / 1e6,
        sim_pct=100 * dif.area / a.area,
        hausdorff=a.boundary.hausdorff_distance(b.boundary),
    )


if __name__ == "__main__":
    print("=" * 92)
    print("ÁREAS GEOESTADÍSTICAS MUNICIPALES — COMPARACIÓN DE LAS TRES FUENTES")
    print("=" * 92)
    print(f"\n  oficial      {SHP_OFICIAL.name}.shp  ({CRS_OFICIAL} → {CRS_TRABAJO},"
          f" dbf en {CODIFICACION_DBF})")
    print(f"  derivado     {GEOJSON_DERIVADO.name}")
    print(f"  provisional  {GEOJSON_PROVISIONAL.name}")

    ofi = cargar_oficial()
    der = cargar_derivado()

    print("\n" + "─" * 92)
    print("OFICIAL FRENTE AL DERIVADO DE TERCEROS")
    print("─" * 92)
    print(f"\n  {'clave':>6} {'municipio':24}{'área ofi.':>12}{'área der.':>12}"
          f"{'dif.':>9}{'dif. simétrica':>17}{'Hausdorff':>12}")
    print("  " + "─" * 90)
    for c in SEIS:
        r = comparar(ofi[c], der[c])
        print(f"  {c:>6} {NOMBRES[c]:24}{r['area_a']:>10,.2f} km²"
              f"{r['area_b']:>10,.2f} km²"
              f"{100*(r['area_b']/r['area_a']-1):>+8.2f}%"
              f"{r['sim']:>11,.3f} km²{r['sim_pct']:>5.2f}%"
              f"{r['hausdorff']:>10,.1f} m")

    print("\n" + "─" * 92)
    print("OFICIAL FRENTE A LA CAPA PROVISIONAL (solo Tlajomulco)")
    print("─" * 92)
    prov = cargar_provisional()
    r = comparar(ofi["14097"], prov)
    print(f"\n  área oficial                    {r['area_a']:>12,.2f} km²")
    print(f"  área provisional                {r['area_b']:>12,.2f} km²"
          f"   ({100*(r['area_b']/r['area_a']-1):+.2f} %)")
    print(f"  solo en la oficial              {r['solo_a']:>12,.3f} km²")
    print(f"  solo en la provisional          {r['solo_b']:>12,.3f} km²")
    print(f"  diferencia simétrica            {r['sim']:>12,.3f} km²"
          f"   ({r['sim_pct']:.2f} % del área oficial)")
    print(f"  distancia de Hausdorff          {r['hausdorff']:>12,.1f} m")

    print("\n" + "=" * 92)
    print("FIN — desde aquí, el polígono canónico es el oficial")
    print("=" * 92)
