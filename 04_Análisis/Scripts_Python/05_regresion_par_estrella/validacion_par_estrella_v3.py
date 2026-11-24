#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validación del registro cartográfico del par comparativo 1971 / 1993.
Versión 3 (ronda 4 del 20-10).

Qué corrige y qué añade respecto de la versión 2
------------------------------------------------
1. CORRECCIÓN DE IDENTIFICACIÓN. En el archivo de PE1 (6 filas, ids 0–5 en el
   orden del georreferenciador), la fila 2 (id 2; 677024.64 E, 2272104.51 N;
   20.5395 N, −103.3020 O) es el GCP sobre la curva del río; la fila 4
   (id 4; 678689.65 E, 2269622.43 N; 20.5169 N, −103.2863 O) está sobre una
   carretera. La versión 2 retiraba la fila 4 creyendo que era el río (lista
   `keep = [0, 1, 2, 3, 5]`). Verificado visualmente el 20-10 (capturas v6).
   Aquí el conjunto de 5 puntos es el de las filas 0, 1, 3, 4, 5 y los puntos
   se identifican por id y por coordenadas, no por posición en una lista.
2. PE1 se ajusta en dos versiones: (i) los 5 GCP originales sin la fila 2 y
   (ii) los 7 GCP activos del `_v2.points` (los 5 anteriores + C y D).
3. PUNTOS DE CONTROL independientes (no entran nunca al ajuste): error de
   predicción por punto y en RMS, por modelo y por fecha.
4. SENSIBILIDAD: PE1 con y sin el control F; PE2 con y sin P2C0 (que está en
   el techo del edificio Aeroméxico).
5. Regla de elección (fijada antes de ver las cifras): se elige el modelo con
   menor error de predicción en los puntos de control. Si la elección cambia
   entre escenarios, el script NO elige: lo declara.
6. Geometría y huella del cuadro de 1971 sobre el AGEM 14 097 oficial
   (Marco Geoestadístico 2025), con el ajuste correcto.
7. Cotejo del conjunto ERRÓNEO del 28-09 (filas 0, 1, 2, 3, 5), solo para
   confirmar que se reproducen 6.77 / 8.74 / 11.28 y señalar qué cambia.

Nomenclatura (siempre el nombre completo; ver v2)
-------------------------------------------------
- RMS de residuos           √(S/n)
- Error estándar corregido  √(S/(n − p/2))
- Error de predicción LOO   leave-one-out real, sobre los GCP
- Error de predicción en puntos de control: RMS de la distancia entre la
  posición que predice el ajuste (hecho solo con los GCP) y la posición
  marcada en la base, para puntos que NO entraron al ajuste.
Todos son error relativo a la base de referencia (Esri Wayback 2026-08-05),
no exactitud. Con 3 controles (o 2, en las sensibilidades) la estimación del
error de predicción es débil.

Entrada : .points de `_gcp_qgis/`, CRS EPSG:32613
Salida  : texto por pantalla y CSV en
          04_Análisis/Visualizaciones/03_regresion_par_estrella/_validacion_v3/
Requiere: numpy, pyproj, shapely, pyshp (los dos últimos para el AGEM)
Uso     : python validacion_par_estrella_v3.py   (desde esta carpeta o no)
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[2]
sys.path.insert(0, str(AQUI))
sys.path.insert(0, str(AQUI.parent / "06_linea_base"))

from validacion_par_estrella_v2 import (          # noqa: E402
    MODELOS, residuos, rms, error_estandar_corregido, prediccion_loo,
    parametros_similitud, FICHA, _aplicar_matriz,
)
from poligonos_v12 import cargar_oficial          # noqa: E402
from pyproj import Transformer                    # noqa: E402
from shapely.geometry import Polygon              # noqa: E402

GCP = AQUI / "_gcp_qgis"
R_PE1_ORIG = GCP / "PE1_0237BZ_0056_041_0200_F1312-Z15B_1971_02_0030K_270.points"
R_PE1_V2 = GCP / "PE1_0237BZ_0056_041_0200_F1312-Z15B_1971_02_0030K_270_v2.points"
R_PE1_CTRL = GCP / "PE1_checkpoints.points"
R_PE2 = GCP / "PE2_0713VE_0001_029_0116_GUADAL-JAL_1993_02_0010K_123.points"
R_PE2_CTRL_V2 = GCP / "PE2_checkpoints_v2.points"     # preferente si existe
R_PE2_CTRL = GCP / "PE2_checkpoints.points"
SALIDA = RAIZ / "04_Análisis/Visualizaciones/03_regresion_par_estrella/_validacion_v3"

# Área de imagen medida en los rásteres originales (v5b, 14-10), en píxeles.
# Es una medición del fotograma, independiente del ajuste de puntos.
LADO_IMAGEN_PX = {"1971": (12710, 12686), "1993": (14306, 14310)}

A_GEO = Transformer.from_crs("EPSG:32613", "EPSG:4326", always_xy=True)
A_UTM = Transformer.from_crs("EPSG:4326", "EPSG:32613", always_xy=True)
NOMBRES_MODELO = list(MODELOS)


# ---------------------------------------------------------------------------
# Lectura: todas las filas, con id (orden del archivo) y bandera enable
# ---------------------------------------------------------------------------

def leer_filas(ruta):
    with open(ruta, encoding="utf-8") as f:
        lineas = [l for l in f if not l.startswith("#")]
    filas = []
    for i, r in enumerate(csv.DictReader(lineas)):
        filas.append(dict(
            id=i,
            mapa=np.array([float(r["mapX"]), float(r["mapY"])]),
            fuente=np.array([float(r["sourceX"]), float(r["sourceY"])]),
            enable=int(float(r["enable"])),
        ))
    return filas


def a_arrays(filas):
    return (np.array([f["mapa"] for f in filas]),
            np.array([f["fuente"] for f in filas]))


def latlon(p):
    lon, lat = A_GEO.transform(p[0], p[1])
    return lat, lon


# ---------------------------------------------------------------------------
# Definición de conjuntos (por id y con verificación por coordenadas)
# ---------------------------------------------------------------------------

def definir_conjuntos():
    orig = leer_filas(R_PE1_ORIG)
    v2 = leer_filas(R_PE1_V2)
    assert len(orig) == 6 and len(v2) == 8, "se esperaban 6 y 8 filas"
    # Las 6 primeras filas del v2 deben ser idénticas a las del original.
    for a, b in zip(orig, v2[:6]):
        assert np.allclose(a["mapa"], b["mapa"]) and np.allclose(a["fuente"], b["fuente"])
    # Verificación de identidad: la fila 2 es el río (el GCP más al norte) y
    # la fila 4 es la de la carretera. Se comprueba por coordenadas.
    rio = orig[2]["mapa"]
    carretera = orig[4]["mapa"]
    assert abs(rio[1] - 2272104.51) < 0.5 and abs(rio[0] - 677024.64) < 0.5
    assert abs(carretera[1] - 2269622.43) < 0.5 and abs(carretera[0] - 678689.65) < 0.5
    assert v2[2]["enable"] == 0 and all(f["enable"] for i, f in enumerate(v2) if i != 2)

    pe2 = leer_filas(R_PE2)
    assert len(pe2) == 7 and all(f["enable"] for f in pe2)

    ruta_c2 = R_PE2_CTRL_V2 if R_PE2_CTRL_V2.exists() else R_PE2_CTRL
    c1 = leer_filas(R_PE1_CTRL)
    c2 = leer_filas(ruta_c2)
    assert len(c1) == 3 and len(c2) == 3

    return dict(
        pe1_i=[orig[i] for i in (0, 1, 3, 4, 5)],
        pe1_ii=[f for f in v2 if f["enable"]],
        pe1_err=[orig[i] for i in (0, 1, 2, 3, 5)],      # conjunto erróneo 28-09
        pe1_seis=orig,
        pe2=pe2,
        c1=dict(zip(["B", "E", "F"], c1)),
        c2=dict(zip(["P2C0 (techo)", "P2C1 (circular)", "P2C2 (cruce N3)"], c2)),
        ruta_c2=ruta_c2.name,
    )


# ---------------------------------------------------------------------------
# Ajuste, LOO y control
# ---------------------------------------------------------------------------

def evaluar(filas, controles):
    """Ajusta cada modelo con los GCP y evalúa en los controles.

    controles: dict etiqueta -> fila. Devuelve dict modelo -> resultados.
    """
    mapa, fuente = a_arrays(filas)
    n = len(mapa)
    salida = {}
    for mod, cfg in MODELOS.items():
        p, min_pts = cfg["p"], cfg["min_pts"]
        r = dict(p=p, n=n, redundancia=2 * n - p)
        if n < min_pts:
            r["calculable"] = False
            salida[mod] = r
            continue
        H = cfg["ajustar"](fuente, mapa)
        res = residuos(H, fuente, mapa)
        loo, nota = prediccion_loo(cfg["ajustar"], p, min_pts, fuente, mapa)
        r.update(calculable=True, H=H, res=res, rms=rms(res),
                 ese=error_estandar_corregido(res, p),
                 loo=loo, loo_rms=(rms(loo) if loo is not None else None),
                 loo_redundancia=2 * (n - 1) - p, nota_loo=nota)
        # Control
        if controles:
            cm = np.array([c["mapa"] for c in controles.values()])
            cf = np.array([c["fuente"] for c in controles.values()])
            e = residuos(H, cf, cm)
            r["ctrl"] = dict(zip(controles.keys(), e))
        salida[mod] = r
    return salida


def rms_ctrl(r, etiquetas):
    return rms([r["ctrl"][k] for k in etiquetas])


def imprimir_ajuste(titulo, filas, res, etiquetas_gcp=None):
    n = len(filas)
    print(f"\n{'=' * 92}\n{titulo} — {n} GCP\n{'=' * 92}")
    for mod, r in res.items():
        print(f"\n  ── {mod}  (p = {r['p']}; 2n = {2 * n}; redundancia 2n − p = {r['redundancia']})")
        if r["redundancia"] <= 0:
            print("     AJUSTE EXACTAMENTE DETERMINADO (o sobredeterminado de menos): "
                  "el error de ajuste no mide nada")
        if not r["calculable"]:
            print("     no calculable")
            continue
        print(f"     RMS de residuos            √(S/n)         {r['rms']:>8.2f} m")
        ese = r["ese"]
        print(f"     error estándar corregido   √(S/(n−{r['p'] // 2}))     "
              + (f"{ese:>8.2f} m" if ese is not None else "   n/c (sin grados de libertad)"))
        if r["loo"] is None:
            print(f"     error de predicción LOO                   no calculable")
        else:
            txt = f"{r['loo_rms']:>8.2f} m   (redundancia por iteración: {r['loo_redundancia']})"
            print(f"     error de predicción LOO                   {txt}")
            if r["loo_redundancia"] <= 0:
                print("       NO EVALUABLE: cada iteración queda exactamente determinada; "
                      "el valor no es un resultado del modelo")
            elif r["nota_loo"]:
                print(f"       nota: {r['nota_loo']}")


def imprimir_loo_por_punto(filas, res):
    print("\n  Error de predicción LOO por GCP (m)   [id de fila del archivo]")
    cab = f"    {'id':>3}  {'lat':>9} {'lon':>10}"
    for mod in NOMBRES_MODELO:
        cab += f" │ {mod[:11]:>11}"
    print(cab)
    for i, f in enumerate(filas):
        lat, lon = latlon(f["mapa"])
        linea = f"    {f['id']:>3}  {lat:>9.5f} {lon:>10.5f}"
        for mod in NOMBRES_MODELO:
            r = res[mod]
            if r.get("loo") is None:
                linea += f" │ {'—':>11}"
            else:
                m = "*" if r["loo_redundancia"] <= 0 else " "
                linea += f" │ {r['loo'][i]:>10.2f}{m}"
        print(linea)
    print("    (* = iteración exactamente determinada: no evaluable)")


def imprimir_control(titulo, res, etiquetas, filas_gcp, controles):
    print(f"\n  Error de predicción en puntos de control — {titulo}  (n = {len(etiquetas)})")
    cab = f"    {'control':<18} {'lat':>9} {'lon':>10}  {'dmin GCP':>8}"
    for mod in NOMBRES_MODELO:
        cab += f" │ {mod[:11]:>11}"
    print(cab)
    gm, _ = a_arrays(filas_gcp)
    for k in etiquetas:
        c = controles[k]
        lat, lon = latlon(c["mapa"])
        dmin = float(np.min(np.hypot(*(gm - c["mapa"]).T)))
        linea = f"    {k:<18} {lat:>9.5f} {lon:>10.5f}  {dmin:>7.0f}m"
        for mod in NOMBRES_MODELO:
            linea += f" │ {res[mod]['ctrl'][k]:>11.2f}"
        print(linea)
    linea = f"    {'RMS (titular)':<18} {'':>9} {'':>10}  {'':>8}"
    for mod in NOMBRES_MODELO:
        linea += f" │ {rms_ctrl(res[mod], etiquetas):>11.2f}"
    print(linea)


def elegir(res, etiquetas):
    """Modelo con menor RMS de error de predicción en los controles."""
    cifras = {m: rms_ctrl(res[m], etiquetas) for m in NOMBRES_MODELO
              if res[m]["calculable"]}
    orden = sorted(cifras, key=cifras.get)
    return orden[0], cifras


# ---------------------------------------------------------------------------
# Geometría y huella
# ---------------------------------------------------------------------------

def geometria(clave, filas, H_sim):
    mapa, fuente = a_arrays(filas)
    f = FICHA[clave]
    esc, rot = parametros_similitud(H_sim)
    cx, cy = A_UTM.transform(f["centro_lon"], f["centro_lat"])
    al_norte = int(np.sum(mapa[:, 1] > cy))
    ancho_px, alto_px = np.ptp(fuente[:, 0]), np.ptp(fuente[:, 1])
    env = 100 * ancho_px * alto_px / (f["ancho_px"] * f["alto_px"])
    lx, ly = LADO_IMAGEN_PX[clave]
    return dict(escala=esc, rot=rot, al_norte=al_norte, n=len(mapa), env=env,
                envol_km=(ancho_px * esc / 1000, alto_px * esc / 1000),
                media_huella=(lx + ly) / 4 * esc)


def huella(clave, H, poli):
    f = FICHA[clave]
    W, Hh = f["ancho_px"], f["alto_px"]
    px = np.array([[0, 0], [W, 0], [W, -Hh], [0, -Hh]], dtype=float)
    esq = _aplicar_matriz(H, px)
    h = Polygon(esq)
    dentro = 100 * h.intersection(poli).area / h.area
    ys = esq[:, 1]
    y0, y1 = ys.min(), ys.max()
    bandas = []
    for i in range(5):
        a = y1 - (i + 1) * (y1 - y0) / 5
        b = y1 - i * (y1 - y0) / 5
        x0, x1 = esq[:, 0].min() - 1e4, esq[:, 0].max() + 1e4
        banda = h.intersection(Polygon([(x0, a), (x1, a), (x1, b), (x0, b)]))
        bandas.append(100 * banda.intersection(poli).area / banda.area
                      if not banda.is_empty and banda.area > 0 else float("nan"))
    return h.area / 1e6, dentro, bandas


# ---------------------------------------------------------------------------
# Principal
# ---------------------------------------------------------------------------

def main():
    SALIDA.mkdir(parents=True, exist_ok=True)
    C = definir_conjuntos()
    filas_csv_ajuste, filas_csv_ctrl, matrices = [], [], {}

    print("#" * 92)
    print("# validacion_par_estrella_v3.py — ronda 4 del 20-10")
    print("# base de referencia: Esri Wayback 2026-08-05 (release 26334); CRS EPSG:32613")
    print(f"# controles de PE2 leídos de: {C['ruta_c2']}")
    print("#" * 92)

    # --- 0. Cotejo con el conjunto ERRÓNEO del 28-09 --------------------------
    r_err = evaluar(C["pe1_err"], None)
    r_seis = evaluar(C["pe1_seis"], None)
    print("\n" + "=" * 92)
    print("0. COTEJO (solo para reproducir el 28-09; el conjunto de 5 puntos era el equivocado)")
    print("=" * 92)
    s = r_err["similitud (Helmert)"]
    print(f"  Conjunto del 28-09 (filas 0,1,2,3,5: incluye el río, excluye la carretera):")
    print(f"    Helmert: RMS de residuos {s['rms']:.2f} · error estándar corregido {s['ese']:.2f} · LOO {s['loo_rms']:.2f}")
    s6 = r_seis["similitud (Helmert)"]
    print(f"  6 puntos del archivo original (Helmert): RMS de residuos {s6['rms']:.2f} · "
          f"error estándar corregido {s6['ese']:.2f} · LOO {s6['loo_rms']:.2f}")
    print(f"    LOO del punto del río (fila 2) en el conjunto de 6:        {s6['loo'][2]:.2f} m")
    print(f"    LOO del punto de la carretera (fila 4) en el conjunto de 6: {s6['loo'][4]:.2f} m")
    print("    (el 22.40 m del 28-09 era el de la carretera, no el del río)")

    # --- 1. PE1 (i) y (ii); PE2 ---------------------------------------------
    esc = {}
    for clave_esc, titulo, filas, ctrl in [
        ("pe1_i", "1971 (PE1) — (i) 5 GCP sin la fila 2 (río)", C["pe1_i"], C["c1"]),
        ("pe1_ii", "1971 (PE1) — (ii) 7 GCP activos del v2 (los 5 + C + D)", C["pe1_ii"], C["c1"]),
        ("pe2", "1993 (PE2) — 7 GCP originales", C["pe2"], C["c2"]),
    ]:
        res = evaluar(filas, ctrl)
        esc[clave_esc] = (filas, ctrl, res)
        imprimir_ajuste(titulo, filas, res)
        imprimir_loo_por_punto(filas, res)
        for mod, r in res.items():
            if not r["calculable"]:
                continue
            matrices[f"{clave_esc}|{mod}"] = r["H"].tolist()
            filas_csv_ajuste.append(dict(
                escenario=clave_esc, modelo=mod, n_gcp=r["n"], p=r["p"],
                redundancia=r["redundancia"], rms_residuos_m=round(r["rms"], 4),
                error_estandar_corregido_m=(round(r["ese"], 4) if r["ese"] is not None else ""),
                loo_rms_m=(round(r["loo_rms"], 4) if r["loo_rms"] is not None else ""),
                loo_redundancia=r["loo_redundancia"],
                loo_evaluable=("no" if r["loo_redundancia"] <= 0 else "sí")))

    # --- 2. Error de predicción en controles ---------------------------------
    print("\n" + "=" * 92)
    print("2. ERROR DE PREDICCIÓN EN PUNTOS DE CONTROL (no entraron al ajuste)")
    print("=" * 92)
    escenarios = []     # (fecha, etiqueta, clave_esc, lista de controles)
    for clave_esc in ("pe1_i", "pe1_ii"):
        for et, lista in [("con B,E,F", ["B", "E", "F"]), ("sin F (B,E)", ["B", "E"])]:
            escenarios.append(("1971", clave_esc, et, lista))
    for et, lista in [("con 3 controles", ["P2C0 (techo)", "P2C1 (circular)", "P2C2 (cruce N3)"]),
                      ("sin P2C0 (P2C1,P2C2)", ["P2C1 (circular)", "P2C2 (cruce N3)"])]:
        escenarios.append(("1993", "pe2", et, lista))

    elecciones = {}
    for fecha, clave_esc, et, lista in escenarios:
        filas, ctrl, res = esc[clave_esc]
        imprimir_control(f"{clave_esc} {et}", res, lista, filas, ctrl)
        ganador, cifras = elegir(res, lista)
        elegido = ganador
        orden = sorted(cifras, key=cifras.get)
        margen = cifras[orden[1]] / cifras[orden[0]]
        print(f"    → menor error de predicción en controles: {ganador} "
              f"({cifras[ganador]:.2f} m; siguiente: {orden[1]} {cifras[orden[1]]:.2f} m, "
              f"{margen:.2f}×)")
        elecciones[(fecha, clave_esc, et)] = (ganador, cifras)
        for k in lista:
            for mod in NOMBRES_MODELO:
                filas_csv_ctrl.append(dict(
                    fecha=fecha, escenario=clave_esc, controles=et, control=k,
                    modelo=mod, error_prediccion_m=round(res[mod]["ctrl"][k], 4)))
        for mod in NOMBRES_MODELO:
            filas_csv_ctrl.append(dict(
                fecha=fecha, escenario=clave_esc, controles=et, control="RMS",
                modelo=mod, error_prediccion_m=round(rms_ctrl(res[mod], lista), 4)))

    # --- 3. ¿Cambia la elección? ---------------------------------------------
    print("\n" + "=" * 92)
    print("3. SENSIBILIDAD: ¿CAMBIA LA ELECCIÓN DE MODELO?")
    print("=" * 92)
    for fecha in ("1971", "1993"):
        sub = {k: v[0] for k, v in elecciones.items() if k[0] == fecha}
        unicos = set(sub.values())
        print(f"\n  {fecha}:")
        for k, g in sub.items():
            print(f"    {k[1]:<7} {k[2]:<24} → {g}")
        if len(unicos) == 1:
            print(f"    LA ELECCIÓN NO CAMBIA en ningún escenario: {unicos.pop()}")
        else:
            print(f"    LA ELECCIÓN CAMBIA ENTRE ESCENARIOS ({', '.join(sorted(unicos))}): "
                  "no se elige; decide [1.2] con la revisión visual")

    # --- 4. Geometría y huella con el ajuste correcto --------------------------
    print("\n" + "=" * 92)
    print("4. GEOMETRÍA Y HUELLA SOBRE EL AGEM 14 097 (Marco Geoestadístico 2025, oficial)")
    print("=" * 92)
    poli = cargar_oficial()["14097"]
    print(f"  AGEM 14 097: {poli.area / 1e6:,.1f} km² (EPSG:32613)")
    comparables = [("1971", "pe1_err", "ERRÓNEO 28-09 (cotejo)", C["pe1_err"], r_err),
                   ("1971", "pe1_i", "(i) 5 GCP", C["pe1_i"], esc["pe1_i"][2]),
                   ("1971", "pe1_ii", "(ii) 7 GCP", C["pe1_ii"], esc["pe1_ii"][2]),
                   ("1993", "pe2", "1993, 7 GCP", C["pe2"], esc["pe2"][2])]
    for fecha, clave_esc, etiqueta, filas, res in comparables:
        g = geometria(fecha, filas, res["similitud (Helmert)"]["H"])
        print(f"\n  {fecha} — {etiqueta}")
        print(f"    Helmert: resolución {g['escala']:.4f} m/px · rotación {g['rot']:+.2f}° · "
              f"media huella (área de imagen medida) {g['media_huella']:,.0f} m")
        print(f"    GCP al norte del centro de toma: {g['al_norte']} de {g['n']} · "
              f"envolvente {g['env']:.1f} % del cuadro "
              f"({g['envol_km'][0]:.2f} × {g['envol_km'][1]:.2f} km)")
        for mod in NOMBRES_MODELO:
            if not res[mod]["calculable"] or res[mod]["redundancia"] <= 0:
                continue
            area, dentro, bandas = huella(fecha, res[mod]["H"], poli)
            bt = " / ".join(f"{b:.1f}" for b in bandas)
            print(f"    huella con {mod:<20}: {area:6.2f} km² · dentro del AGEM {dentro:5.1f} % "
                  f"· bandas N→S (%): {bt}")

    # --- 5. Controles frente a la envolvente convexa de los GCP ---------------
    # Un modelo con más parámetros se comporta peor fuera de la zona que
    # sostienen los GCP. Se declara antes de interpretar nada qué controles
    # caen dentro y cuáles fuera de la envolvente, y qué tan alargada es la
    # nube de GCP (cociente de valores singulares de las coordenadas centradas).
    from shapely.geometry import MultiPoint, Point
    print("\n" + "=" * 92)
    print("5. CONTROLES FRENTE A LA ENVOLVENTE CONVEXA DE LOS GCP")
    print("=" * 92)
    for clave_esc, nombre in [("pe1_i", "1971 (i) 5 GCP"), ("pe1_ii", "1971 (ii) 7 GCP"),
                              ("pe2", "1993 7 GCP")]:
        filas, ctrl, _ = esc[clave_esc]
        gm, _f = a_arrays(filas)
        casco = MultiPoint([tuple(p) for p in gm]).convex_hull
        sv = np.linalg.svd(gm - gm.mean(axis=0), compute_uv=False)
        print(f"\n  {nombre}: envolvente {casco.area / 1e6:.2f} km² · "
              f"cociente de ejes de la nube de GCP {sv[0] / sv[1]:.1f}")
        for k, c in ctrl.items():
            p = Point(*c["mapa"])
            if casco.contains(p):
                est = f"dentro (a {casco.exterior.distance(p):.0f} m del borde)"
            else:
                est = f"FUERA, a {casco.distance(p):.0f} m de la envolvente"
            print(f"    {k:<18} {est}")

    # --- Archivos ---------------------------------------------------------------
    import json
    for nombre, filas in [("ajuste_por_escenario_y_modelo.csv", filas_csv_ajuste),
                          ("error_prediccion_controles.csv", filas_csv_ctrl)]:
        with open(SALIDA / nombre, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(filas[0]))
            w.writeheader()
            w.writerows(filas)
    with open(SALIDA / "matrices_pixel_a_UTM32613.json", "w", encoding="utf-8") as fh:
        json.dump(matrices, fh, indent=1)
    print(f"\nSalidas en {SALIDA}")
    print("Los estimadores son error relativo a la base, no exactitud. Con 3 controles "
          "(2 en las sensibilidades) la estimación del error de predicción es débil.")


if __name__ == "__main__":
    main()
