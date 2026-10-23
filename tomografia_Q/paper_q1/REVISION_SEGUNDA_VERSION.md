# Revisión de la segunda versión

Sugerencias contrastadas con campañas archivadas; no se modificaron los
resultados originales ni se publicó contenido externo.

1. **Tabla 1:** dos paneles legibles, cuatro decimales en los diagnósticos,
   correlación entre componentes en las cuatro filas y Corr(predicción, dato).
   Se declara el mismo vector retenido para RMS y momentos dentro de cada
   campaña. Se explica por qué la regularización impide asumir ortogonalidad
   exacta residuo–predicción. La identidad de covarianza sí se verifica.
2. **Comparación convencional:** RMS inicial de referencia 0,12281 s y
   expandida 0,20672 s; datos 0,07018 y 0,06941 s. Se retira “opposite
   impression”: son mejoras respecto de predicciones iniciales diferentes.
3. **Recorte:** E_RMS = −67,71 % para Q máximo 300 y −56,77 % para 1500.
   Ambas predicciones recortadas ajustan peor que cero, con Vp y términos
   de estación fijos. Son recortes de una única solución sin restricciones;
   la diferencia no valida ni selecciona un techo. Reajustar los términos
   de estación podría reducir los residuos y no se ha hecho en esta prueba.
4. **Bootstrap por eventos:** 10 000 réplicas, semilla 20261005, percentiles
   2,5/97,5. Se mantienen juntos todos los pares de cada evento. La referencia
   tiene 679 eventos con pares retenidos; las tres campañas expandidas, 975.
   Las diferencias expandidas usan los mismos sorteos y pares.

   | Configuración | E_RMS (%) | Intervalo condicional 95 % |
   |---|---:|---:|
   | Referencia | 12,1128 | 10,22–14,15 |
   | Expandida heredada | 12,5855 | 10,90–14,23 |
   | Dominio de rayos | 15,5116 | 14,10–16,99 |
   | Adoptada | 15,0472 | 13,65–16,51 |

   Ganancia pareada del dominio: 2,9261 puntos; intervalo 2,0577–4,1425.
   Cambio con damping adoptado: −0,4644 puntos; intervalo −0,6144 a −0,3296.
   Es variabilidad condicional de predicciones fijas: no se reinvierte cada
   réplica, no estima incertidumbre de Q, ni dependencia entre eventos por
   estaciones comunes o errores sistemáticos.
5. **Experimento controlado 2×2:** mismos 52 914 pares, rejilla, Q inicial,
   cuadratura y regularización original (damping de estación 0,15).

   | Vp del núcleo | Máscara | E_RMS (%) | RMS residual (s) |
   |---|---|---:|---:|
   | Heredado | Heredada | 12,5855 | 0,06068 |
   | Heredado | Rayos | 15,5060 | 0,05865 |
   | Actualizado | Heredada | 12,9083 | 0,06045 |
   | Actualizado | Rayos | 15,5116 | 0,05865 |

   A Vp fijo heredado, la máscara aporta 2,9205 puntos. Dentro de la máscara
   nueva, actualizar Vp añade 0,0056. Hay interacción entre factores: a
   máscara heredada fija, el cambio de Vp aporta 0,3228 puntos. Las diagonales
   reproducen los RMS archivados con diferencias menores de 4,1×10⁻⁹ s;
   pequeñas diferencias de suma y parada iterativa producen máximas
   diferencias puntuales de 5,31×10⁻⁵ s. Los criterios explícitos son
   10⁻⁷ s en RMS y 10⁻⁴ s por predicción. Estas corridas son lineales sin
   restricciones y no validan los mapas recortados.
6. **Rejilla vertical:** nodos 0, 20, …, 220 km listados completos. El soporte
   almacenado excede los dominios y conserva valores iniciales en nodos
   inactivos. Profundidad activa máxima de ambas ramas: 160 km en campañas
   heredadas y 140 km en las de dominio de rayos.
7. **Iluminación:** 1410 nodos tienen ≥1 impacto entre 8000 pares muestreados.
   Ese umbral procede del escalamiento y redondeo; no es el criterio ≥5
   aplicado al conjunto completo para activar nodos en inversión, ni prueba
   resolución suficiente. Se corrige también el resumen y las conclusiones.
8. **Límite inferior:** 0/1169 nodos activos (0 %) alcanzan Q=20, tanto con
   techo 300 como 1500.
9. **Archivo:** se prepara `archivo_reproducible_borrador.zip`, con manifiesto
   SHA-256, entorno, datos procesados y scripts de auditoría. Es un borrador
   local, sin DOI. No se inventa ni reserva un identificador. Para depósito
   público faltan cuenta del autor, licencia y aprobación de redistribución.
   El README delimita la reproducibilidad: auditoría desde tablas/modelos
   archivados; no se incluyen formas de onda crudas ni todos los insumos de
   los experimentos separados de checkerboard y losa.
10. **Novedad:** formulada como auditoría operativa de covarianza, distinción
    entre invariancia del RMS archivado y no-identificabilidad del límite,
    y evaluación forward del modelo acotado presentado. No se reivindica
    descubrir el trade-off conocido ni una interpretación geológica validada.

Puntos menores: typo “this the” corregido; encabezado abreviado; λ_τ y
Q_max verificados; ruido observado identificado como RMS adoptado previo
al recorte. Las leyendas 7/8 estaban completas en la fuente y se comprueban
en el PDF. Zhao 2005 ya incluía 603–627 en el BibTeX y aparece completo en
el PDF; Audemard M. y Fonseca P. mantienen apellidos/partículas según la
fuente editorial verificada. La bibliografía conserva 47 obras verificadas,
sin añadir citas durante esta revisión. Se mantiene la semilla original de
iluminación: es un entero reproducible, no una afirmación cronológica.

Además se corrige el convertidor GJI: sus entornos de revisión no consumen
el argumento opcional `[t]`, que se estaba imprimiendo encima de figuras y
tablas. Se elimina al convertir, sin alterar el formato editorial original.

La declaración solicitada menciona Codex para revisión de lenguaje, sin
atribuir autoría a la herramienta. El nuevo código de diagnóstico está
disponible para revisión y validación de los autores.

## Reproducción

`segunda_revision_diagnosticos.py` genera resultados y sorteos en
`segunda_revision/`. `genera_segunda_revision.py` genera
`cifras_segunda_revision.tex`. Las correcciones editoriales se conservan en
`manuscrito.tex`; ejecutar `convierte_gji.py`, compilar y después
`cuenta_palabras_gji.py`. Los respaldos de la versión anterior están en
`segunda_revision/antes_*`. No ejecutar los restauradores/intermedios de
revisiones anteriores sobre el manuscrito final.
