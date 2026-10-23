"""Emite `cifras.tex`: TODAS las cifras del manuscrito, leidas del pipeline.

Por que existe
--------------
Regla del proyecto: ningun resultado numerico se escribe a mano en LaTeX.
Cada cifra del manuscrito tiene que salir de datos archivados o del
pipeline de analisis. Este guion lee los JSON de las campanas y escribe
macros; el manuscrito solo usa macros. Si una campana se recalcula, el
texto se actualiza solo y no queda ninguna cifra huerfana.

Ademas evita el error contrario, que ya ocurrio: citar de memoria una
linea base equivocada (476/822 celdas) cuando el archivo decia otra cosa
(365/442).

Que NO se toma de aqui
----------------------
Nada de `repositorio_q1_tomografia_fmm/v2`. Ese trabajo es la referencia
de FORMA del manuscrito -- damero, marco teorico, geologia del VMM,
densidad de citas -- no una fuente de datos. Los numeros de este papel
salen solo de `tomografia_Q`.

Uso: py tomografia_Q/paper_q1/genera_cifras.py [--adoptada <campana>]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
SALIDAS = RAIZ / "output_V2"

# La campana publicada, contra la que se compara todo.
PREVIA = "sgc_final_q_extra_station_tomography"
# Campanas intermedias que documentan el argumento del articulo.
MASCARA_VIEJA = "sgc_q_ampliado_300km"     # red ampliada, mascara antigua
MASCARA_RAYOS = "sgc_q_mascara_rayos"      # mascara por rayos, damp original


def lee(p: Path):
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def met(campana: str):
    return lee(SALIDAS / campana / "sgc_1000_tomography_metrics.json")


def macro(nombre: str, valor) -> str:
    return f"\\newcommand{{\\{nombre}}}{{{valor}}}\n"


def mil(n) -> str:
    """Separador de millar de LaTeX-friendly: 52\\,914."""
    return f"{int(n):,}".replace(",", "\\,")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--adoptada", default="sgc_q_mascara_rayos_damp400")
    a = ap.parse_args()

    m_ad = met(a.adoptada)
    m_pr = met(PREVIA)
    m_mv = met(MASCARA_VIEJA)
    m_mr = met(MASCARA_RAYOS)
    if m_ad is None:
        sys.exit(f"falta la campaña adoptada: {a.adoptada}")
    if m_pr is None:
        sys.exit(f"falta la campaña previa: {PREVIA}")

    desc = lee(RAIZ / "descompone_ajuste.json") or []
    por_campana = {d["campana"]: d for d in desc}
    rayos = lee(RAIZ / "mascara_por_rayos.json")
    inv = lee(RAIZ / "inventario_cobertura_q.json")
    amp = lee(RAIZ / "amplia_estaciones_q.json")

    L = ["% Generado por paper_q1/genera_cifras.py -- NO editar a mano.\n",
         f"% Campaña adoptada: {a.adoptada}\n\n"]

    # ------------------------------------------------ red y datos de entrada
    L.append("% --- red, eventos y observaciones ---\n")
    L.append(macro("estPrevia", m_pr["stations_used"]))
    L.append(macro("estAdoptada", m_ad["stations_used"]))
    # Las corridas derivadas reutilizan los pares y no reescriben el
    # recuento de eventos: se toma de la campana de origen, que es la misma
    # entrada de datos.
    ev = (m_ad.get("xml_events_available_for_run")
          or m_ad.get("vp_events_used")
          or (m_mv or {}).get("xml_events_available_for_run"))
    L.append(macro("eventosAdoptada", mil(ev) if ev else "??"))
    evp = m_pr.get("xml_events_available_for_run") or m_pr.get("vp_events_used")
    L.append(macro("eventosPrevia", mil(evp) if evp else "??"))
    L.append(macro("picksPrevia", mil(m_pr["vp_picks"])))
    L.append(macro("picksAdoptada", mil(m_ad["vp_picks"])))
    if m_mv:
        L.append(macro("trazasInventario", mil(m_mv["trazas_en_inventario"])))
        L.append(macro("trazasConEspectro", mil(m_mv["trace_spectra"])))
    L.append(macro("paresPrevia", mil(m_pr["q_pair_observations_total"])))
    L.append(macro("paresAdoptada", mil(m_ad["q_pair_observations_total"])))
    if inv and inv.get("hueco_azimutal"):
        h = inv["hueco_azimutal"]
        L.append(macro("huecoPrevio", f"{h['mediana_actual']:.1f}"))
        L.append(macro("huecoAmpliado", f"{h['mediana_ampliado']:.1f}"))
        L.append(macro("fracHuecoPrevio",
                       f"{100 * h['frac_mayor_180_actual']:.0f}"))
        L.append(macro("fracHuecoAmpliado",
                       f"{100 * h['frac_mayor_180_ampliado']:.0f}"))
    if inv:
        tot = inv["lecturas_dentro"] + inv["lecturas_fuera"]
        L.append(macro("estConLecturas", len(inv["estaciones_fuera"])
                       + len(inv["estaciones_usadas"])))
        L.append(macro("fracLecturasDescartadas",
                       f"{100 * inv['lecturas_fuera'] / tot:.0f}"))
    if amp:
        L.append(macro("radioAmpliacion", f"{amp['radio_km']:.0f}"))

    # ----------------------------------------------- parametrizacion y malla
    L.append("\n% --- dominio y parametrización ---\n")
    for etq, m in (("Previa", m_pr), ("Adoptada", m_ad)):
        d = m["domain"]
        L.append(macro(f"lonMin{etq}", f"{d['lon_min']:.2f}"))
        L.append(macro(f"lonMax{etq}", f"{d['lon_max']:.2f}"))
        L.append(macro(f"latMin{etq}", f"{d['lat_min']:.2f}"))
        L.append(macro(f"latMax{etq}", f"{d['lat_max']:.2f}"))
        if d.get("max_model_depth_km"):
            L.append(macro(f"profMax{etq}", f"{d['max_model_depth_km']:.0f}"))
    if m_ad["domain"].get("celdas_en_mascara"):
        L.append(macro("celdasEnMascara",
                       mil(m_ad["domain"]["celdas_en_mascara"])))
    if rayos:
        L.append(macro("celdasIluminadas", mil(rayos["celdas_con_cobertura"])))
        L.append(macro("celdasFueraMascaraVieja",
                       mil(rayos["celdas_con_cobertura_fuera_de_la_actual"])))
        L.append(macro("fracCeldasFuera",
                       f"{100 * rayos['celdas_con_cobertura_fuera_de_la_actual'] / rayos['celdas_con_cobertura']:.0f}"))
        L.append(macro("fracImpactosFuera",
                       f"{100 * rayos['fraccion_impactos_fuera_de_la_actual']:.0f}"))
        L.append(macro("paresTrazados", mil(rayos["pares_trazados"])))

    # -------------------------------------------------------- regularizacion
    reg = m_ad.get("regularizacion")
    if reg:
        L.append("\n% --- regularización ---\n")
        L.append(macro("qLaplace", f"{reg['Q_LAPLACE']:g}"))
        L.append(macro("qCellDamp", f"{reg['Q_CELL_DAMP']:g}"))
        L.append(macro("qStationDamp", f"{reg['Q_STATION_DAMP']:g}"))
        L.append(macro("qStationDampPrevio", f"{reg['Q_STATION_DAMP_previo']:g}"))

    # --------------------------------------------------- resultados Vp y Q
    L.append("\n% --- inversión Vp ---\n")
    for etq, m in (("Previa", m_pr), ("Adoptada", m_ad)):
        v = m["vp"]
        L.append(macro(f"vpCeldas{etq}", mil(v["active_cells"])))
        L.append(macro(f"vpRms{etq}", f"{v['final_rms_s']:.3f}"))
        L.append(macro(f"vpRazon{etq}", f"{v['data_parameter_ratio']:.1f}"))
        L.append(macro(f"vpIter{etq}", v["lsmr_iterations"]))

    L.append("\n% --- inversión Q ---\n")
    # La campana intermedia entra tambien: el articulo compara las tres.
    for etq, m in (("Previa", m_pr), ("MascaraVieja", m_mv),
                   ("MascaraRayos", m_mr), ("Adoptada", m_ad)):
        if m is None:
            continue
        q = m["q_pilot"]
        L.append(macro(f"qCeldas{etq}", mil(q["active_cells"])))
        L.append(macro(f"qRms{etq}", f"{q['final_rms_s']:.4f}"))
        L.append(macro(f"qTerminos{etq}", q["station_terms"]))
        L.append(macro(f"qRazonDatos{etq}",
                       f"{q['data_parameter_ratio_cells_plus_stations']:.1f}"))
        L.append(macro(f"qMedianaRdos{etq}", f"{q['median_r2']:.3f}"))
        L.append(macro(f"qObs{etq}", mil(q["observations_used"])))

    # --------------------------------------- descomposición estructura/estación
    L.append("\n% --- descomposición del ajuste (la medida que vale) ---\n")
    etiquetas = {PREVIA: "Previa", MASCARA_VIEJA: "MascaraVieja",
                 MASCARA_RAYOS: "MascaraRayos", a.adoptada: "Adoptada"}
    for camp, etq in etiquetas.items():
        d = por_campana.get(camp)
        if not d:
            continue
        L.append(macro(f"expl{etq}", f"{d['explica_pct']:.1f}"))
        L.append(macro(f"razonEst{etq}",
                       f"{d['razon_estaciones_sobre_dato']:.2f}"))
        L.append(macro(f"razonEstr{etq}",
                       f"{d['razon_estructura_sobre_dato']:.2f}"))
        L.append(macro(f"rmsDato{etq}", f"{d['rms_dato']:.4f}"))
        L.append(macro(f"rmsPartida{etq}", f"{d['rms_modelo_de_partida']:.4f}"))
        L.append(macro(f"rmsFinal{etq}", f"{d['rms_final']:.4f}"))
        L.append(macro(f"rmsRedGuion{etq}",
                       f"{d['rms_reduction_percent_del_guion']:.1f}"))

    # ------------------------------------------------ damero (resolución)
    dam = lee(RAIZ / f"damero_{a.adoptada}.json")
    dam015 = lee(RAIZ / f"damero_{a.adoptada}_damp015.json")
    if dam:
        L.append("\n% --- damero ---\n")
        L.append(macro("damPerturbVp", f"{100*dam['perturbacion_vp_frac']:.0f}"))
        L.append(macro("damPerturbQ", f"{100*dam['perturbacion_qinv_frac']:.0f}"))
        lams = sorted({r["lambda_km"] for r in dam["resultados"]})
        L.append(macro("damLambdas",
                       ", ".join(f"{x:.0f}" for x in lams) + "~km"))
        romano = {40: "A", 80: "B", 120: "C", 160: "D"}
        for r in dam["resultados"]:
            suf = romano.get(int(r["lambda_km"]))
            if suf is None:
                continue
            etq = r["rama"].capitalize() + ("Bajo" if r["ruido"] == "bajo"
                                            else "Obs") + suf
            L.append(macro(f"dam{etq}", f"{r['correlacion']:.3f}"))
            L.append(macro(f"damAmp{etq}", f"{r['razon_amplitud']:.2f}"))
            if r.get("desv_terminos_estacion") is not None:
                L.append(macro(f"damTerm{etq}",
                               f"{r['desv_terminos_estacion']:.4f}"))
        if dam015:
            for r in dam015["resultados"]:
                suf = romano.get(int(r["lambda_km"]))
                if suf is None or r["ruido"] != "observado":
                    continue
                L.append(macro(f"damViejoQObs{suf}",
                               f"{r['correlacion']:.3f}"))
                L.append(macro(f"damViejoTermQObs{suf}",
                               f"{r['desv_terminos_estacion']:.4f}"))

    # ------------------------------------------- techo de Q: prior vs dato
    techo = lee(RAIZ / "techo_q.json")
    if techo:
        L.append("\n% --- techo de Q ---\n")
        L.append(macro("qTecho", f"{techo['base']['techo']:.0f}"))
        L.append(macro("qTechoPrueba", f"{techo['prueba']['techo']:.0f}"))
        L.append(macro("fracEnTecho",
                       f"{100*techo['base']['frac_en_techo']:.0f}"))
        L.append(macro("fracEnTechoPrueba",
                       f"{100*techo['prueba']['frac_en_techo']:.0f}"))
        L.append(macro("cambioRmsTecho", f"{techo['cambio_rms_s']:.5f}"))
        lb = techo["celdas_libres_del_techo"]
        L.append(macro("celdasLibresTecho", lb["n"]))
        L.append(macro("fracCeldasLibres",
                       f"{100*lb['n']/techo['base']['celdas_activas']:.0f}"))
        L.append(macro("qLibreMediana", f"{lb['q_mediana_base']:.0f}"))
        L.append(macro("qLibreDiez", f"{lb['q_p10_base']:.0f}"))
        L.append(macro("qLibreNoventa", f"{lb['q_p90_base']:.0f}"))
        L.append(macro("cambioQLibre",
                       f"{lb['cambio_relativo_mediano_pct']:.1f}"))
        L.append(macro("fracVuelveATopar",
                       f"{100*techo['celdas_que_topaban']['frac_que_vuelve_a_topar']:.0f}"))
        # Los nombres de macro de LaTeX no admiten digitos: `\qLibreZ0` se
        # leeria como `\qLibreZ` seguido de un `0` suelto. Van en palabras.
        NOMBRE_Z = {0: "Sup", 20: "Veinte", 40: "Cuarenta", 100: "Cien"}
        for r in techo.get("perfil_profundidad", []):
            z = NOMBRE_Z.get(int(r["z_km"]))
            if z and r["q_mediana_libre"]:
                L.append(macro(f"qLibreZ{z}", f"{r['q_mediana_libre']:.0f}"))
                L.append(macro(f"fracTechoZ{z}",
                               f"{100*r['frac_en_techo']:.0f}"))

    # ----------------------------------------------------- losa de Nazca
    losa = lee(RAIZ / f"losa_nazca_{a.adoptada}.json")
    if losa:
        L.append("\n% --- losa de Nazca ---\n")
        L.append(macro("losaControles", mil(losa["controles"])))
        L.append(macro("losaNodos", mil(losa["nodos_superficie"])))
        pr = losa["profundidad_km"]
        L.append(macro("losaProfMin", f"{pr['min']:.0f}"))
        L.append(macro("losaProfMax", f"{pr['max']:.0f}"))
        bz = losa["buzamiento"]
        L.append(macro("losaBuzamiento", f"{bz['local_mediano_grados']:.0f}"))
        L.append(macro("losaBuzamientoNoventa",
                       f"{bz['local_p90_grados']:.0f}"))
        L.append(macro("losaTendencia",
                       f"{bz['tendencia_global_grados']:.0f}"))
        dc = losa["distancia_a_control_km"]
        L.append(macro("losaDistControl", f"{dc['mediana']:.1f}"))
        s2 = losa.get("contra_slab2")
        if s2:
            L.append(macro("slabDifMediana",
                           f"{s2['diferencia_mediana_km']:+.1f}"))
            L.append(macro("slabDifAbs",
                           f"{s2['diferencia_abs_mediana_km']:.1f}"))
            L.append(macro("slabDifNoventa",
                           f"{s2['diferencia_p90_abs_km']:.0f}"))
            L.append(macro("slabNodos", mil(s2["nodos_comparados"])))

    fuera = AQUI / "cifras.tex"
    fuera.write_text("".join(L), encoding="utf-8")
    n = sum(1 for x in L if x.startswith("\\newcommand"))
    print(f"  {fuera.name}: {n} macros")
    print(f"  campaña adoptada : {a.adoptada}")
    print(f"  campañas en la descomposición: "
          f"{', '.join(e for c, e in etiquetas.items() if c in por_campana)}")
    faltan = [c for c in etiquetas if c not in por_campana]
    if faltan:
        print(f"  *** sin descomposición (correr descompone_ajuste.py): "
              f"{', '.join(faltan)}")


if __name__ == "__main__":
    main()
