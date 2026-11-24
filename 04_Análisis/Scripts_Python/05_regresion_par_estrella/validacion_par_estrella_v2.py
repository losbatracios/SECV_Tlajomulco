"""
Validación del registro cartográfico del par comparativo 1971 / 1993.
Versión 2.

Qué añade respecto de la versión 1
----------------------------------
1. Error de predicción leave-one-out real: para cada punto de control se
   ajusta el modelo con los n−1 restantes y se predice el punto excluido.
   En la versión 1, lo que se llamaba "dejar uno fuera" era el error de
   ajuste del reajuste sin ese punto, que es otra cosa: mide cuánto baja
   el error del conjunto reducido, no cuánto se equivoca el modelo al
   predecir un punto que no vio. El nombre queda corregido.
2. Comparación de tres modelos de transformación: similitud (Helmert),
   afín y proyectiva, cada uno con sus grados de libertad declarados.
3. Geometría de los puntos de control: posición de cada uno respecto del
   centro de toma declarado en el índice de vuelos, y extensión de la
   envolvente de puntos dentro del cuadro.

Nomenclatura de los estimadores (se usa siempre el nombre completo)
-------------------------------------------------------------------
r      residuo bidimensional de un punto; S = Σ r²; n = número de puntos;
       p = número de parámetros del modelo.
- RMS de residuos            = √(S/n)
- Error estándar corregido   = √(S/(n − p/2))
- Error de predicción LOO    = RMS de los residuos de predicción, cada uno
                               calculado con el modelo ajustado sin su punto
Los tres son error relativo a la base de referencia, no exactitud.

Redundancia
-----------
Cada punto aporta 2 ecuaciones, así que el ajuste tiene 2n − p ecuaciones
redundantes. Con 2n = p el ajuste pasa por los puntos y el error de ajuste
es cero por construcción: no mide nada. El leave-one-out de un modelo
necesita que n−1 puntos dejen redundancia; si no, se declara y no se
reporta cifra.

Entrada : .points exportados por el georreferenciador de QGIS
CRS     : EPSG:32613 (WGS 84 / UTM 13N)
Requiere: numpy; pyproj es opcional (solo para imprimir latitud/longitud)
"""

import csv
import math
from pathlib import Path

import numpy as np

try:
    from pyproj import Transformer
    _A_GEO = Transformer.from_crs("EPSG:32613", "EPSG:4326", always_xy=True)
except ImportError:
    _A_GEO = None


# ---------------------------------------------------------------------------
# Rutas
# ---------------------------------------------------------------------------
# Todas las rutas se resuelven desde la ubicación de este archivo, de modo que
# el script funciona desde cualquier directorio de trabajo y lee los archivos
# canónicos del repositorio, sin copias locales.
AQUI = Path(__file__).resolve().parent          # .../04_Análisis/Scripts_Python/05_regresion_par_estrella
RAIZ = AQUI.parents[2]                          # raíz del repositorio del proyecto

RUTA_1971 = AQUI / "_gcp_qgis/PE1_0237BZ_0056_041_0200_F1312-Z15B_1971_02_0030K_270.points"
RUTA_1993 = AQUI / "_gcp_qgis/PE2_0713VE_0001_029_0116_GUADAL-JAL_1993_02_0010K_123.points"

# Inventario de calidad del corpus: una fila por fotograma, con las dimensiones
# del ráster entregado por el archivo. Es la fuente de `ancho_px` y `alto_px`.
CSV_CALIDAD = RAIZ / "04_Análisis/Visualizaciones/01_inspeccion_calidad/calidad_corpus.csv"

# Identificador de carpeta de cada toma del par, tal como aparece en el
# inventario de calidad.
CARPETA = {
    "1971": "0237BZ_0056_041_0200_F1312-Z15B_1971_02_0030K_270",
    "1993": "0713VE_0001_029_0116_GUADAL-JAL_1993_02_0010K_123",
}


def dimensiones_raster(carpeta, ruta_csv=CSV_CALIDAD):
    """Ancho y alto en píxeles de un fotograma, leídos del inventario de calidad.

    Se lee del archivo en lugar de escribirse a mano para que la cifra tenga una
    sola procedencia y cualquier recálculo del inventario se propague al ajuste.
    """
    with open(ruta_csv, newline="", encoding="utf-8") as fh:
        for fila in csv.DictReader(fh):
            if fila["carpeta"] == carpeta:
                return int(fila["ancho_px"]), int(fila["alto_px"])
    raise KeyError(f"{carpeta} no aparece en {ruta_csv}")


# Centro de toma declarado en el índice de vuelos (consulta SIIV 07-07).
FICHA = {
    "1971": dict(bote="237BZ", linea=56, foto=41, escala_nominal=30000,
                 centro_lat=20.5339, centro_lon=-103.2926),
    "1993": dict(bote="713VE", linea=1, foto=29, escala_nominal=10000,
                 centro_lat=20.5333, centro_lon=-103.2920),
}

for _fecha, _ficha in FICHA.items():
    _ficha["carpeta"] = CARPETA[_fecha]
    _ficha["ancho_px"], _ficha["alto_px"] = dimensiones_raster(CARPETA[_fecha])


# ---------------------------------------------------------------------------
# Lectura
# ---------------------------------------------------------------------------

def leer_points(ruta):
    """Devuelve (coordenadas de terreno, coordenadas de píxel, residuo QGIS).

    La primera línea define el CRS y empieza con '#'. Solo se conservan los
    puntos activos (enable = 1).
    """
    with open(ruta, encoding="utf-8") as f:
        lineas = [l for l in f if not l.startswith("#")]
    filas = [r for r in csv.DictReader(lineas) if int(float(r["enable"])) == 1]
    mapa = np.array([[float(r["mapX"]), float(r["mapY"])] for r in filas])
    fuente = np.array([[float(r["sourceX"]), float(r["sourceY"])] for r in filas])
    residuo_qgis = np.array([float(r["residual"]) for r in filas])
    return mapa, fuente, residuo_qgis


# ---------------------------------------------------------------------------
# Modelos de transformación
#
# Los tres se resuelven por mínimos cuadrados sobre coordenadas normalizadas
# (centradas en su media y escaladas a distancia media unitaria) y después se
# devuelven a las unidades originales. La normalización no cambia el modelo:
# evita que la diferencia de magnitud entre píxeles (~1e4) y coordenadas UTM
# (~1e6) degrade la solución, que es crítico en la proyectiva.
# ---------------------------------------------------------------------------

def _normalizador(pts):
    """Matriz 3x3 que centra en la media y escala a distancia media √2."""
    c = pts.mean(axis=0)
    d = np.hypot(*(pts - c).T).mean()
    s = math.sqrt(2) / d if d > 0 else 1.0
    T = np.array([[s, 0, -s * c[0]],
                  [0, s, -s * c[1]],
                  [0, 0, 1.0]])
    return T


def _homogeneas(pts):
    return np.column_stack([pts, np.ones(len(pts))])


def _aplicar_matriz(H, pts):
    q = _homogeneas(pts) @ H.T
    return q[:, :2] / q[:, 2:3]


def ajustar_similitud(fuente, mapa):
    """Helmert: rotación + escala uniforme + traslación. 4 parámetros."""
    n = len(fuente)
    A = np.zeros((2 * n, 4))
    y = np.zeros(2 * n)
    for i, (sx, sy) in enumerate(fuente):
        A[2 * i] = [sx, -sy, 1.0, 0.0]
        y[2 * i] = mapa[i, 0]
        A[2 * i + 1] = [sy, sx, 0.0, 1.0]
        y[2 * i + 1] = mapa[i, 1]
    a, b, tx, ty = np.linalg.lstsq(A, y, rcond=None)[0]
    return np.array([[a, -b, tx], [b, a, ty], [0, 0, 1.0]])


def ajustar_afin(fuente, mapa):
    """Afín: 6 parámetros (escalas distintas por eje + corte + traslación)."""
    Tf, Tm = _normalizador(fuente), _normalizador(mapa)
    f = _aplicar_matriz(Tf, fuente)
    m = _aplicar_matriz(Tm, mapa)
    n = len(f)
    A = np.zeros((2 * n, 6))
    y = np.zeros(2 * n)
    for i, (sx, sy) in enumerate(f):
        A[2 * i] = [sx, sy, 1, 0, 0, 0]
        y[2 * i] = m[i, 0]
        A[2 * i + 1] = [0, 0, 0, sx, sy, 1]
        y[2 * i + 1] = m[i, 1]
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    Hn = np.array([[c[0], c[1], c[2]], [c[3], c[4], c[5]], [0, 0, 1.0]])
    return np.linalg.inv(Tm) @ Hn @ Tf


def ajustar_proyectiva(fuente, mapa):
    """Proyectiva (homografía): 8 parámetros. Modelo de una toma inclinada
    sobre terreno plano. Se resuelve por DLT con normalización."""
    Tf, Tm = _normalizador(fuente), _normalizador(mapa)
    f = _aplicar_matriz(Tf, fuente)
    m = _aplicar_matriz(Tm, mapa)
    A = []
    for (x, y), (X, Y) in zip(f, m):
        A.append([x, y, 1, 0, 0, 0, -X * x, -X * y, -X])
        A.append([0, 0, 0, x, y, 1, -Y * x, -Y * y, -Y])
    _, _, Vt = np.linalg.svd(np.array(A))
    Hn = Vt[-1].reshape(3, 3)
    H = np.linalg.inv(Tm) @ Hn @ Tf
    return H / H[2, 2]


MODELOS = {
    "similitud (Helmert)": dict(ajustar=ajustar_similitud, p=4, min_pts=2),
    "afín":                dict(ajustar=ajustar_afin,      p=6, min_pts=3),
    "proyectiva":          dict(ajustar=ajustar_proyectiva, p=8, min_pts=4),
}


# ---------------------------------------------------------------------------
# Estimadores de error
# ---------------------------------------------------------------------------

def residuos(H, fuente, mapa):
    d = mapa - _aplicar_matriz(H, fuente)
    return np.hypot(d[:, 0], d[:, 1])


def rms(valores):
    return math.sqrt(float(np.mean(np.asarray(valores) ** 2)))


def error_estandar_corregido(res, p):
    """√(S/(n − p/2)). Devuelve None si no quedan grados de libertad."""
    n = len(res)
    gl = n - p / 2
    if gl <= 0:
        return None
    return math.sqrt(float(np.sum(res ** 2)) / gl)


def prediccion_loo(ajustar, p, min_pts, fuente, mapa):
    """Error de predicción leave-one-out real.

    Para cada punto i: se ajusta con los n−1 restantes y se predice i. El
    residuo de esa predicción es el error de predicción de ese punto.
    Devuelve (lista de errores, nota) o (None, motivo) si no es calculable.
    """
    n = len(mapa)
    if n - 1 < min_pts:
        return None, (f"n−1 = {n - 1} puntos no alcanzan el mínimo de "
                      f"{min_pts} que el modelo necesita")
    redundancia = 2 * (n - 1) - p
    errores = []
    for i in range(n):
        k = [j for j in range(n) if j != i]
        H = ajustar(fuente[k], mapa[k])
        errores.append(float(residuos(H, fuente[i:i + 1], mapa[i:i + 1])[0]))
    nota = None
    if redundancia == 0:
        nota = ("el ajuste de cada iteración queda exactamente determinado "
                "(2(n−1) = p): pasa por sus puntos y la predicción es inestable")
    elif redundancia <= 2:
        nota = (f"redundancia mínima en cada iteración (2(n−1) − p = "
                f"{redundancia}): la predicción es inestable")
    return np.array(errores), nota


# ---------------------------------------------------------------------------
# Informes
# ---------------------------------------------------------------------------

def a_geo(x, y):
    if _A_GEO is None:
        return None
    lon, lat = _A_GEO.transform(x, y)
    return lat, lon


def parametros_similitud(H):
    a, b = H[0, 0], H[1, 0]
    return math.hypot(a, b), math.degrees(math.atan2(b, a))


def informe_modelos(nombre, fuente, mapa, residuo_qgis=None, etiquetas=None):
    n = len(mapa)
    etiquetas = etiquetas or list(range(1, n + 1))
    print(f"\n{'=' * 88}\n{nombre}  —  {n} puntos de control\n{'=' * 88}")

    resultados = {}
    for mod, cfg in MODELOS.items():
        p, min_pts = cfg["p"], cfg["min_pts"]
        if n < min_pts:
            print(f"\n  {mod}: no calculable, requiere al menos {min_pts} puntos")
            continue

        H = cfg["ajustar"](fuente, mapa)
        res = residuos(H, fuente, mapa)
        redundancia = 2 * n - p
        ese = error_estandar_corregido(res, p)
        loo, nota_loo = prediccion_loo(cfg["ajustar"], p, min_pts, fuente, mapa)

        print(f"\n  ── {mod}  (p = {p}; ecuaciones 2n = {2*n}; "
              f"redundancia 2n − p = {redundancia})")
        if redundancia == 0:
            print("     ajuste exactamente determinado: el error de ajuste es "
                  "cero por construcción y no mide nada")
        print(f"     RMS de residuos            √(S/n)       {rms(res):>7.2f} m")
        if ese is None:
            print(f"     error estándar corregido   √(S/(n−{p//2}))   "
                  f"no calculable (sin grados de libertad)")
        else:
            print(f"     error estándar corregido   √(S/(n−{p//2}))   {ese:>7.2f} m")
        if loo is None:
            print(f"     error de predicción LOO                 no calculable: {nota_loo}")
        else:
            print(f"     error de predicción LOO                 {rms(loo):>7.2f} m")
            if nota_loo:
                print(f"       nota: {nota_loo}")

        resultados[mod] = dict(H=H, res=res, loo=loo, rms=rms(res), ese=ese,
                               loo_rms=rms(loo) if loo is not None else None,
                               nota_loo=nota_loo, redundancia=redundancia)

    # Detalle por punto
    print(f"\n  Detalle por punto de control")
    cab = f"  {'GCP':>4}  {'latitud':>9} {'longitud':>11}"
    for mod in resultados:
        cab += f" │ {mod[:11]:>11} {'LOO':>7}"
    print(cab)
    print("  " + "─" * (len(cab) - 2))
    for i in range(n):
        g = a_geo(mapa[i, 0], mapa[i, 1])
        coord = (f"{g[0]:>9.5f} {g[1]:>11.5f}" if g
                 else f"{mapa[i,0]:>10.1f} {mapa[i,1]:>10.1f}")
        fila = f"  {etiquetas[i]:>4}  {coord}"
        for mod, r in resultados.items():
            loo_i = (f"{r['loo'][i]:>7.2f}" if r["loo"] is not None else f"{'—':>7}")
            fila += f" │ {r['res'][i]:>11.2f} {loo_i}"
        print(fila)

    # Verificación contra el archivo: los residuos del ajuste de similitud
    # deben coincidir con los que QGIS guardó en el .points. Solo tiene
    # sentido sobre el conjunto completo del archivo: si se retira un punto,
    # el ajuste cambia y los residuos guardados ya no son comparables.
    if residuo_qgis is not None and "similitud (Helmert)" in resultados:
        esc, _ = parametros_similitud(resultados["similitud (Helmert)"]["H"])
        dif = np.abs(resultados["similitud (Helmert)"]["res"] / esc - residuo_qgis)
        print(f"\n  Residuos recalculados frente a los guardados en el .points: "
              f"desviación máxima {dif.max():.4f} px")

    return resultados


def geometria_puntos(clave, fuente, mapa, H_similitud, etiquetas=None):
    """Posición de los puntos respecto del centro de toma y del cuadro."""
    f = FICHA[clave]
    n = len(mapa)
    etiquetas = etiquetas or list(range(1, n + 1))
    esc, rot = parametros_similitud(H_similitud)

    print(f"\n  Geometría de los puntos de control — {clave}")
    print(f"    escala ajustada {esc:.4f} m/px · rotación {rot:+.2f}°")

    # Centro de toma declarado, convertido a UTM con la misma proyección.
    if _A_GEO is not None:
        inv = Transformer.from_crs("EPSG:4326", "EPSG:32613", always_xy=True)
        cx, cy = inv.transform(f["centro_lon"], f["centro_lat"])
        print(f"    centro de toma declarado en el índice de vuelos: "
              f"{f['centro_lat']}, {f['centro_lon']}  →  E {cx:.0f}  N {cy:.0f}")
        print(f"\n    {'GCP':>4}  {'ΔE(m)':>8} {'ΔN(m)':>8}  {'dist(m)':>8}  posición")
        print("    " + "─" * 52)
        for i in range(n):
            de, dn = mapa[i, 0] - cx, mapa[i, 1] - cy
            cuad = ("N" if dn > 0 else "S") + ("E" if de > 0 else "O")
            print(f"    {etiquetas[i]:>4}  {de:>+8.0f} {dn:>+8.0f}  "
                  f"{math.hypot(de, dn):>8.0f}  {cuad}")
        al_norte = int(np.sum(mapa[:, 1] > cy))
        print(f"\n    al norte del centro de toma: {al_norte} de {n}")

    # Envolvente de los puntos dentro del cuadro, en píxeles y en metros.
    ancho_px, alto_px = np.ptp(fuente[:, 0]), np.ptp(fuente[:, 1])
    area_gcp = ancho_px * alto_px
    area_cuadro = f["ancho_px"] * f["alto_px"]
    print(f"\n    envolvente de los puntos: {ancho_px:.0f} × {alto_px:.0f} px"
          f"  =  {ancho_px*esc/1000:.2f} × {alto_px*esc/1000:.2f} km")
    print(f"    cuadro completo:         {f['ancho_px']} × {f['alto_px']} px"
          f"  =  {f['ancho_px']*esc/1000:.2f} × {f['alto_px']*esc/1000:.2f} km")
    print(f"    la envolvente cubre el {100*area_gcp/area_cuadro:.1f} % "
          f"del área del cuadro")


if __name__ == "__main__":
    # ── 1971 ───────────────────────────────────────────────────────────────
    # El archivo trae 6 puntos. El 5.º se marcó sobre una curva de río y se
    # retiró por principio (no se puede sostener la estabilidad de un cauce
    # entre 1971 y la imagen de referencia). Se reportan los dos conjuntos.
    m1, f1, q1 = leer_points(RUTA_1971)

    print("\n" + "#" * 88)
    print("# 1971 — 0237BZ_0056_041 (nominal 1:30,000) — conjunto de 6 puntos")
    print("# incluye el punto sobre la curva del río, retirado del ajuste definitivo")
    print("#" * 88)
    r6 = informe_modelos("1971, 6 puntos (de referencia)", f1, m1, q1)

    keep = [0, 1, 2, 3, 5]
    print("\n" + "#" * 88)
    print("# 1971 — conjunto definitivo de 5 puntos (sin el punto del río)")
    print("#" * 88)
    r5 = informe_modelos("1971, 5 puntos (definitivo)", f1[keep], m1[keep],
                         residuo_qgis=None, etiquetas=[i + 1 for i in keep])
    geometria_puntos("1971", f1[keep], m1[keep],
                     r5["similitud (Helmert)"]["H"],
                     etiquetas=[i + 1 for i in keep])

    # Error de predicción del punto del río bajo el ajuste de los otros cinco:
    # es el leave-one-out de ese punto dentro del conjunto de 6.
    if r6["similitud (Helmert)"]["loo"] is not None:
        print(f"\n  Punto sobre la curva del río (GCP 5), error de predicción "
              f"LOO en el conjunto de 6: "
              f"{r6['similitud (Helmert)']['loo'][4]:.2f} m")

    # ── 1993 ───────────────────────────────────────────────────────────────
    m2, f2, q2 = leer_points(RUTA_1993)
    print("\n" + "#" * 88)
    print("# 1993 — 0713VE_0001_029 (nominal 1:10,000) — 7 puntos")
    print("#" * 88)
    r7 = informe_modelos("1993, 7 puntos", f2, m2, q2)
    geometria_puntos("1993", f2, m2, r7["similitud (Helmert)"]["H"])

    # ── Resumen comparativo ────────────────────────────────────────────────
    print(f"\n{'=' * 88}\nResumen: error de predicción LOO por modelo\n{'=' * 88}")
    print(f"  {'conjunto':<28} {'modelo':<22} {'RMS resid':>10} "
          f"{'est. corr.':>11} {'LOO':>9}")
    print("  " + "─" * 84)
    for etiqueta, r in [("1971, 6 puntos", r6), ("1971, 5 puntos (definitivo)", r5),
                        ("1993, 7 puntos", r7)]:
        for mod, v in r.items():
            ese = f"{v['ese']:>11.2f}" if v["ese"] is not None else f"{'n/d':>11}"
            loo = f"{v['loo_rms']:>9.2f}" if v["loo_rms"] is not None else f"{'n/d':>9}"
            print(f"  {etiqueta:<28} {mod:<22} {v['rms']:>10.2f} {ese} {loo}")
    print("\n  Los tres estimadores son error relativo a la base de referencia,")
    print("  no exactitud. La elección de modelo se resuelve con los puntos de")
    print("  control independientes, que no entran al ajuste.")
