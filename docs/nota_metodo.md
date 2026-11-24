# Nota de método

Complemento del capítulo «Los límites de la mirada. Archivo fotogramétrico y transformación socioambiental» (Montes Ruiz y Bermeo Donoso). Reúne el detalle técnico que no cabía en el capítulo. Todas las cifras de esta nota están en `04_Análisis/cifras_validadas_v2.csv`, con su nombre exacto, su tamaño de muestra, el script que la produce y su salvedad.

## 1. Corpus

- **Fuente.** Índice de vuelos del Sistema de Información de Índices de Vuelos (SIIV) del INEGI, consulta del 7 de julio de 2026 sobre los seis municipios centrales de la Zona Metropolitana de Guadalajara (`02_Archivo_INEGI/Metadatos_SINFA/`). El índice registra 2,074 fotogramas en 45 botes y 12 años entre 1969 y 2003.
- **Unidad territorial.** Área geoestadística municipal (AGEM) 14 097 del Marco Geoestadístico 2025 del INEGI, de 671.33 km². Es un área estadística, no un límite político-administrativo oficial. La capa provisional que se usó antes difiere de la oficial en 0.008 km² (distancia de Hausdorff de 12.4 m) y ninguna cifra municipal cambia al sustituirla.
- **Pertenencia.** Un fotograma se cuenta en un municipio cuando su centro de toma declarado en el índice cae dentro del AGEM. Con ese criterio, 226 fotogramas tienen centro en Tlajomulco de Zúñiga: 3 de 1970, 155 de 1971, 43 de 1978, 23 de 1981 y 2 de 1993.

## 2. Conjuntos solicitados al INEGI

- **A.** 12 fotogramas de 1993 que registran territorio del AGEM (umbral de 1,200 m entre el centro de toma y el AGEM; el umbral medido sobre la media huella de 1,209 m da el mismo conjunto).
- **B.** 19 contrapartes previas a la reforma del artículo 27, con centro a menos de 2 km de alguna fotografía de A. La solicitud midió la distancia con la fórmula de haversine sobre una esfera; medida sobre el elipsoide WGS84, el mismo criterio da 20, porque el fotograma 104AZ/4/44/1 queda a 1,998 m. Se analizaron los 19 solicitados.
- **C.** 67 fotogramas de 1971 con centro dentro del AGEM al sur del paralelo 20.48° N.

Las 98 fotografías llegaron como escaneos con su ficha XML. Sus metadatos declaran como restricción de uso «Ninguna». No se incluyen en este repositorio.

## 3. Cobertura del archivo: centro de toma y huella nominal

- **Densidad.** Entre 1969 y 1983 Tlajomulco tiene 0.334 fotogramas por km², en una posición intermedia entre los seis municipios (de 0.248 en El Salto a 0.689 en Guadalajara). Entre 1992 y 1994 tiene 0.003 por km² (2 fotogramas), frente a 1.166 en Guadalajara.
- **Líneas de vuelo de 1992 a 1994.** Son 48. Medidas por la traza que une sus centros de toma, 1 cruza el AGEM, 2 pasan a 2 km o menos y 45 quedan más lejos. La mediana de distancia de sus extremos al AGEM es de 14,760 m. La línea 713VE/1 aporta los dos fotogramas de 1993 con centro en el municipio, en las posiciones 28 y 29 de sus 29 fotos.
- **Huella nominal.** Cuadrado de 0.23 m por el denominador de la escala declarada, centrado en el centro de toma; no corrige inclinación ni relieve. Medido así, en 1993 alcanzan el municipio 18 de 356 fotogramas (14 a 1:10,000 y 4 a 1:75,000, del bote 720VE), y su unión cubre 28.70 km², el 4.27 % del AGEM. En 1992, 1994 y 2003 ninguna huella lo toca. Por huella, las 48 líneas se reparten en 3 que tocan el AGEM, 2 a 2 km o menos y 43 a más de 2 km. El conteo cambia y el resultado no: después de 1992 el archivo casi no registra el municipio.
- **Nombres de vuelo.** 849 de los 2,074 fotogramas (40.9 %) pertenecen a vuelos llamados «Guadalajara»; ninguno lleva el nombre de Tlajomulco.

## 4. Regresión cartográfica

### 4.1 Procedimiento

Los puntos de control terrestre (GCP) se marcaron a ojo en QGIS sobre rasgos visibles en la fotografía y reconocibles en la base de referencia, la versión del 5 de agosto de 2026 (release 26334) de Esri World Imagery Wayback. Los puntos de verificación se marcaron aparte y nunca entraron al ajuste (`_gcp_qgis/*_checkpoints.points`). Python calcula todas las cifras a partir de los archivos `.points`; QGIS solo sirve para el juicio visual. Sistema de referencia: EPSG:32613.

El modelo de cada fotografía se eligió por el menor error de predicción en los puntos de verificación, regla fijada antes de ver las cifras, comparando similitud (Helmert), afín y proyectiva.

### 4.2 Nombres del error

Cada error lleva su nombre exacto. Todos miden la distancia a la base de referencia, cuya exactitud propia no se evalúa; ninguno es una exactitud absoluta.

- **RMS de residuos**, √(S/n): ajuste a los propios GCP.
- **Error estándar corregido**, √(S/(n − p/2)): el mismo, corregido por grados de libertad (Helmert n−2, afín n−3, proyectiva n−4).
- **Error de predicción leave-one-out**: cada GCP se predice con el ajuste de los demás; mide interpolación dentro de la nube de puntos.
- **Error de predicción en puntos de verificación**: RMS en puntos que no entraron al ajuste. Es la cifra titular cuando existe.

### 4.3 Par del aeropuerto

| | 1971 (237BZ L56 F41) | 1993 (713VE L1 F29) |
|---|---|---|
| Modelo elegido | similitud (Helmert) | afín |
| Error de predicción en puntos de verificación | 17.98 m (n = 3) | 3.52 m (n = 3) |
| Error de predicción leave-one-out | 11.48 m (n = 7) | 3.87 m (n = 7) |
| RMS de residuos | 8.21 m | 2.05 m |
| Error estándar corregido | 9.72 m | 2.71 m |
| Resolución en el terreno del ajuste | 0.4830 m/px | 0.1690 m/px |
| Área del cuadro proyectado | 44.59 km² | 6.57 km² |

Salvedades. Con tres puntos de verificación, Helmert y afín no se distinguen con claridad. En las dos fechas los siete GCP quedan al sur del centro de toma y su caja envolvente cubre el 15.6 % (1971) y el 13.1 % (1993) del cuadro, de modo que la superposición es más confiable en esa zona. La escala de 1971 es una cota (1:26,666), porque sus marcas fiduciales no pasaron la prueba de medición; la de 1993 es de 1:10,283 medida entre fiduciales. Las fotografías no se ortorrectificaron.

### 4.4 Presa La Playa

Fotografía 249BZ L37 F5 (febrero de 1971, 1:35,000), con cinco GCP, tres de ellos a menos de 130 m entre sí en el mismo caserío, de modo que el ajuste se apoya en tres ubicaciones independientes, y dos puntos de verificación.

- El modelo afín da un error de predicción de 28.5 m en los dos puntos de verificación (Helmert, 69.7 m). Con dos puntos la estimación es muy débil, y la elección se invierte si se retira el punto de la cerca.
- El error de predicción leave-one-out es de 255 m, porque al retirar cualquiera de los dos puntos aislados el ajuste queda sin apoyo en su zona. El RMS de residuos, de 4.59 m, refleja la geometría de los puntos y no sirve como validación.
- El contorno del vaso, trazado por su borde tonal, encierra unas 416 ha (perímetro de 11,353 m) e incluye márgenes con parcelas y una vereda. Al cambiar de modelo el contorno se desplaza una mediana de 37 m (máximo 87 m) y su superficie cambia 0.7 %.
- La planta registrada en el DENUE queda fuera del contorno, al oeste-noroeste, a 1,182 m de su borde. La distancia de 1,730 m que circuló antes se midió contra una caja envolvente y está retirada.

## 5. Otros sensores sobre el vaso

- **Cobertura del archivo aéreo.** Por huella nominal, 31 fotogramas contienen el centroide del vaso: 22 de 1971, 5 de 1978 y 4 de 1981. Ninguno posterior a 1981. Los de 1978 y 1981 no se solicitaron.
- **JRC Global Surface Water** (agua abierta a 30 m), sobre una caja de 588 ha. De 454 meses, 275 tienen observada al menos el 90 % de la caja. Ninguno de ellos registra agua antes de 2021, tampoco los 28 febreros bien observados. El único episodio con agua es de julio a noviembre de 2021, entre 1.3 y 15.1 ha (máximo 2.6 % de la caja). El producto detecta agua abierta; un humedal con vegetación emergente, agua somera o canales estrechos puede no quedar registrado.
- **Landsat MSS** (60 m, 1972–1985). El inventario lista 132 escenas, pero el material disponible son 12 mosaicos anuales (sin cobertura en 1974 y 1975). Una primera clasificación visual trató las 132 filas como observaciones independientes; eran 12 imágenes vistas varias veces, y esa clasificación se retiró. Medido sobre los 12 mosaicos, el contorno del vaso es más oscuro en el infrarrojo cercano que su anillo de 500 m en los 12 años (probabilidad de 0.720 a 0.969, mediana 0.938), y el área oscura supera 1 km² en 7 de ellos. Un área oscura en el infrarrojo cercano es compatible con agua, suelo húmedo o vegetación densa. En la ventana común, la serie JRC observa solo 2 de 22 meses, ambos sin agua abierta.

Ninguna de estas fuentes permite fechar la desaparición de la lámina de agua que registra la fotografía de 1971 ni atribuirle una causa.

## 6. Límites declarados

- La pertenencia al municipio se define por el AGEM y por el centro de toma declarado; la huella nominal es una aproximación sin corrección de inclinación ni relieve.
- Los errores de georreferenciación son relativos a una base comercial y se calculan con muy pocos puntos.
- La consulta al SIIV no incluye los productos digitales posteriores del INEGI; la ausencia de fotogramas después de 1993 describe lo que registra ese índice.
- No hubo ortorrectificación ni fotointerpretación codificada.
