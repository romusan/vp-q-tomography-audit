"""Cuanto del ajuste de Q lo hace la ESTRUCTURA y cuanto los TERMINOS DE
ESTACION.

Por que esta medida y no el RMS
-------------------------------
Las metricas del guion dan `rms_reduction_percent`, que se calcula contra
el modelo de partida. Ese modelo predice PEOR que no predecir nada (RMS
0,207 s frente a 0,069 s del dato en la campana ampliada), de modo que el
porcentaje mide cuanto se recupera de un mal punto de partida, no cuan
bien se ajusta el dato. Comparar campanas por ese numero induce a error:
la ampliada marca 70,7 % frente al 49,8 % de la publicada, y sin embargo
la varianza del dato que explican es casi la misma.

Aqui se mide contra el DATO, que es la referencia que no depende del
punto de partida, y ademas se separan las dos componentes del modelo:

    prediccion = (estructura Q en las celdas) + (t_i - t_j)

Si cada componente por separado tiene una amplitud mucho mayor que el
dato y solo su SUMA lo explica, la inversion esta en un compromiso del
espacio nulo: ha encontrado dos cosas grandes y opuestas de las que solo
la diferencia esta determinada. Esa es la patologia que se vigila.

El signo de la parte de estaciones se deduce del propio fichero, no se
supone: se prueban las dos convenciones y se toma la que reproduce
`pred_final_s`.

Uso: py tomografia_Q/descompone_ajuste.py [campana ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = Path(__file__).resolve().parent
SALIDAS = AQUI / "output_V2"
POR_DEFECTO = ["sgc_final_q_extra_station_tomography", "sgc_q_ampliado_300km"]


def rms(v) -> float:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.sqrt(np.mean(np.square(v)))) if v.size else float("nan")


def analiza(nombre: str) -> dict | None:
    b = SALIDAS / nombre
    fp, ft = (b / "data" / "sgc_1000_q_predictions.csv",
              b / "data" / "sgc_1000_q_station_terms.csv")
    if not (fp.exists() and ft.exists()):
        print(f"  (falta {nombre})")
        return None
    p, t = pd.read_csv(fp), pd.read_csv(ft)
    term = dict(zip(t.station.astype(str),
                    pd.to_numeric(t.station_term_s, errors="coerce")))
    ti = p.station_i.astype(str).map(term)
    tj = p.station_j.astype(str).map(term)
    d = pd.to_numeric(p.delta_tstar_s, errors="coerce")
    pf = pd.to_numeric(p.pred_final_s, errors="coerce")
    pi = pd.to_numeric(p.pred_initial_s, errors="coerce")
    a, c = ti - tj, tj - ti
    est = a if (pf - a).std() <= (pf - c).std() else c
    signo = "t_i − t_j" if (pf - a).std() <= (pf - c).std() else "t_j − t_i"
    cel = pf - est
    ok = d.notna() & pf.notna() & est.notna()
    d, pf, pi, est, cel = d[ok], pf[ok], pi[ok], est[ok], cel[ok]
    r0 = rms(d)
    m = json.loads((b / "sgc_1000_tomography_metrics.json")
                   .read_text(encoding="utf-8"))
    return {
        "campana": nombre, "n": int(len(d)), "signo": signo,
        "celdas_activas_q": m.get("q_pilot", {}).get("active_cells"),
        "terminos_de_estacion": m.get("q_pilot", {}).get("station_terms"),
        "rms_dato": r0,
        "rms_modelo_de_partida": rms(d - pi),
        "rms_final": rms(d - pf),
        "explica_pct": 100 * (1 - rms(d - pf) / r0),
        "explica_solo_estaciones_pct": 100 * (1 - rms(d - est) / r0),
        "explica_solo_estructura_pct": 100 * (1 - rms(d - cel) / r0),
        "desv_parte_estaciones": float(est.std()),
        "desv_parte_estructura": float(cel.std()),
        "desv_dato": float(d.std()),
        "razon_estaciones_sobre_dato": float(est.std() / d.std()),
        "razon_estructura_sobre_dato": float(cel.std() / d.std()),
        "rms_reduction_percent_del_guion":
            m.get("q_pilot", {}).get("rms_reduction_percent"),
    }


def main() -> None:
    nombres = sys.argv[1:] or POR_DEFECTO
    r = [x for x in (analiza(n) for n in nombres) if x]
    if not r:
        sys.exit("  sin campañas que analizar")

    def fila(etq, clave, fmt="{:.4f}"):
        print(f"  {etq:<40s}" + "".join(
            f"{fmt.format(x[clave]) if x[clave] is not None else '—':>16s}"
            for x in r))

    print(f"\n  {'':<40s}" + "".join(f"{x['campana'][:15]:>16s}" for x in r))
    fila("observaciones", "n", "{:,}")
    fila("celdas activas Q", "celdas_activas_q", "{}")
    fila("términos de estación", "terminos_de_estacion", "{}")
    print()
    fila("RMS del dato (s)", "rms_dato")
    fila("RMS contra el modelo de partida (s)", "rms_modelo_de_partida")
    fila("RMS final (s)", "rms_final")
    print()
    print("  --- varianza del DATO que explica cada cosa ---")
    fila("modelo completo (%)", "explica_pct", "{:.1f}")
    fila("solo términos de estación (%)", "explica_solo_estaciones_pct", "{:.1f}")
    fila("solo estructura Q (%)", "explica_solo_estructura_pct", "{:.1f}")
    print()
    print("  --- amplitud de cada componente frente al dato ---")
    fila("desviación del dato (s)", "desv_dato")
    fila("desviación parte de estaciones (s)", "desv_parte_estaciones")
    fila("desviación parte de estructura (s)", "desv_parte_estructura")
    fila("razón estaciones / dato", "razon_estaciones_sobre_dato", "{:.2f}")
    fila("razón estructura / dato", "razon_estructura_sobre_dato", "{:.2f}")
    print()
    fila("«rms_reduction_percent» del guion",
         "rms_reduction_percent_del_guion", "{:.1f}")
    print("\n  Lectura: una razón cercana a 1 es sana. Muy por encima de 1 en")
    print("  AMBAS componentes, con el modelo completo explicando poco, indica")
    print("  que se cancelan entre sí: compromiso del espacio nulo.")

    salida = AQUI / "descompone_ajuste.json"
    salida.write_text(json.dumps(r, indent=2, ensure_ascii=False),
                      encoding="utf-8")
    print(f"\n  {salida.name}")


if __name__ == "__main__":
    main()
