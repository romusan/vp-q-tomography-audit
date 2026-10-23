# Revisión de las sugerencias recibidas

Fecha: 2026-10-05. Fuente de trabajo: `manuscrito.tex`; entrega GJI:
`manuscrito_gji.tex` y `manuscrito_gji.pdf`. Se conservaron el manuscrito y PDF
anteriores en `revision_cientifica/antes_*`. No se modificaron ni se volvieron
a invertir las campañas archivadas en `output_V2`.

## Hallazgos que cambian la interpretación

1. **La métrica no era varianza explicada.** `descompone_ajuste.py` calcula
   `100*(1 - RMS(residual)/RMS(dato))`. Se conserva el valor de las campañas,
   pero se renombra como reducción de RMS normalizada y se define mediante
   una ecuación. No es R², reducción de error cuadrático ni
   Var(predicción)/Var(datos). Por ejemplo, 15,0 % de esta métrica corresponde
   aproximadamente a 27,8 % de reducción de error cuadrático respecto a cero,
   sin implicar una partición independiente de las componentes.
2. **La anticorrelación persiste.** Se reconstruyeron las componentes sobre
   todas las observaciones retenidas, usando el signo del operador (+τᵢ−τⱼ).
   Correlaciones: referencia −0,919; ampliada con dominio heredado −0,989;
   dominio por rayos −0,976; amortiguamiento adoptado −0,868. La identidad de
   varianza y covarianza se verifica al redondeo de máquina. Se elimina la
   equivalencia entre razones cercanas a uno y buen condicionamiento.
   La componente estructural incluye la predicción inicial; la ecuación
   previa también omitía ese término y el peso de velocidad del kernel.
3. **El RMS archivado precede al recorte.** `run_q_inversion` resuelve con
   LSMR, recorta Q⁻¹ con `np.clip`, pero almacena `pred0 + A @ x` usando `x`
   sin recortar. Cambiar únicamente el techo deja esa predicción invariable
   por construcción. Evaluando los modelos recortados con la misma
   cuadratura: RMS 0,11641244 s para Qmax=300 y 0,10882063 s para Qmax=1500,
   frente a 0,05896906 s de la predicción lineal sin recorte. Se verificó
   que las velocidades y términos de estación de ambas campañas coinciden.
   La reconstrucción de la predicción inicial sobre los 52 914 pares difiere
   del archivo en menos de 1,2×10⁻¹⁵ s. Se retira la conclusión de que el RMS
   idéntico demuestra que los datos no acotan Q por arriba.
4. **La comparación de referencia no aísla receptores.** La rejilla de
   referencia tiene 13×9×12 nodos, frente a 31×32×12 en las campañas ampliadas.
   Ambas tienen espaciamiento nominal de 20 km, pero distinto soporte
   computacional. La referencia parte de un campo de Vp tridimensional con
   variación lateral por capa de hasta 1,05 km/s; las ampliadas parten de un
   perfil unidimensional distinto. Se documentan estas diferencias y se
   retira la atribución causal exclusiva al número de estaciones.

## Tratamiento de los comentarios mayores

| Sugerencia | Cambio aplicado / criterio |
|---|---|
| 1. Diagnóstico de descomposición | Definición exacta de la métrica, ecuación con modelo inicial, kernel ponderado por Vp, identidad de covarianza y correlaciones reales. Se moderan título, resumen, resultados y conclusiones. |
| 2. Referencia no publicada | Se elimina la mención al manuscrito en preparación. La configuración de referencia se describe aquí: dominio, rejilla, campo inicial, observaciones, umbrales y pesos. Se requieren los arrays archivados para reproducir el campo inicial 3D completo. |
| 3. Uso de “joint” | Título cambiado a *Auditing parameterisation and bound effects in sequential P-wave velocity and spectral-ratio attenuation tomography…*. Resumen y Métodos aclaran la secuencia Vp→Q, sin objetivo simultáneo ni retroalimentación Q→Vp. El ajuste de estructura y términos de estación dentro de la rama Q sí es conjunto. |
| 4. Métodos | Se añaden rejillas, dominios, pesos, perfil inicial ampliado, Q inicial, convergencia, umbrales espectrales, todos los pares disponibles por evento, muestreo de rayos, suavizado de la superficie y trilineal para atributos. Se explica que cuatro pasos de rejilla por una longitud de onda de 80 km no equivalen a resolución general de estructuras arbitrarias. Se informa el 4,3 % de nodos Vp en el techo numérico. |
| 5. Calidad espectral | Se declara ausencia de ponderación por incertidumbre o R², SNR mínimo 3, banda candidata 1–12 Hz, al menos diez bins, rechazo de diferencias no finitas o mayores que 0,75 s en valor absoluto. El 47,1 % de diferencias retenidas son negativas; el signo depende del orden de estaciones y no implica Q negativo. Se reconoce dependencia entre pares que comparten registros. |
| 6. Términos reales de estación | Nueva figura con mapas del amortiguamiento original y adoptado en escala común y respecto a la media de la red. Se reconoce que la amortiguación puede suprimir respuesta real de sitio y que no hay calibración independiente para separarla de compensación estructural. |
| 7. Cortes de Q | Nueva figura vectorial con cortes a 0, 20 y 40 km, techo en gris e inactivos en blanco. Se presentan como salidas recortadas provisionales, no como mediciones validadas. |

## Comentarios menores

- **4221 frente a 4216:** son pares producidos frente a pares aceptados en la
  referencia. Se explican ambos; el resumen usa aceptados para comparar con
  la tabla. No se altera el dato para forzar coincidencia.
- **1114 frente a 1112:** el primero es inventario de eventos, el segundo
  eventos con picks P retenidos y dibujados. Se distingue en texto y leyenda.
- **Q 123–157 frente a 158:** 123 y 157 son medianas de los nodos por debajo
  del techo en las capas de 20 y 40 km. 158 es la mediana global de ese
  subconjunto, con percentiles 10–90 de 64–260. Ninguno demuestra por sí solo
  soporte de los datos tras el hallazgo del recorte.
- **Slab2:** se retira “within the stated uncertainty”, pues no se evaluó un
  campo regional de incertidumbre. Se conserva la diferencia descriptiva y
  se evita afirmar independencia estadística de catálogos no comprobada.
- **Rayos:** trayectorias rectas, 8000 pares aceptados muestreados sin
  reemplazo, semilla archivada 20260924. Vp afecta pesos, no curva las rutas.
- **Generalización:** se añaden escenarios donde receptores nuevos sí
  pueden ayudar; el resultado se presenta como caso particular.
- **“Measurements”:** se retira para los valores de Q del modelo recortado.
  Estar por debajo del techo tampoco garantiza estar determinado por los datos.
- **Autores bibliográficos:** se mantienen las formas publicadas y ya
  verificadas. “Audemard M” y “Fonseca P” son las familias en los metadatos
  editoriales; “Castelblanco Ossa” en la tesis y “Castelblanco” en el artículo
  reflejan sus respectivas publicaciones. No se homogeneizan inventando nombres.
- **Figuras:** ocho PDF generados con Matplotlib. La inspección de objetos
  muestra gráficas, texto y polígonos vectoriales; las imágenes pequeñas
  incrustadas corresponden a gradientes de barras de color, no a capturas
  completas de los paneles. Se revisó la composición de todas las páginas.

## Evidencia, regeneración y límites

`revision_sugerencias.py` obtiene las estadísticas de las tablas CSV y arrays
archivados; produce `revision_cientifica/diagnosticos.json`, las predicciones
recortadas, las dos nuevas figuras y `cifras_revision.tex`. No sobreescribe
campañas. `cifras.tex` mantiene los números originales; las nuevas cifras se
generan separadamente. `figura_curva_l.py` conserva los datos del barrido y
corrige las etiquetas de la métrica y de la línea de amplitud unitaria.

`aplica_sugerencias.py` documenta una transformación intermedia desde el respaldo;
la fuente final contiene ajustes adicionales. No usarlo para regenerar o
restaurar la entrega final. Para futuras ediciones, modificar `manuscrito.tex`
y ejecutar `convierte_gji.py`, compilar y ejecutar `cuenta_palabras_gji.py`.
Los informes de palabras finales están en `CONTEO_PALABRAS_GJI.md` y JSON.

Esta revisión resuelve la presentación honesta de los resultados existentes;
no convierte la campaña en una inversión físicamente acotada validada.
Antes de atribuir significado cuantitativo al modelo de Q hace falta resolver
un objetivo con límites dentro de la optimización, evaluar predicciones desde
el modelo devuelto y repetir las pruebas de sensibilidad pertinentes.
El paquete de reproducción completo también necesita un identificador público
persistente antes de la presentación editorial. No se afirma que ya exista.

Los parámetros por defecto de LSMR y su código de parada se contrastaron con
la instalación local y la [documentación oficial SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsmr.html).
