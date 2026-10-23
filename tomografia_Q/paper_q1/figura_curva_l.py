"""Curva del barrido: inflación de componentes frente a reducción de RMS normalizada. Las razones unitarias no prueban independencia ni buen condicionamiento. Se conservan los datos archivados y se corrigen sus etiquetas tras la revisión científica. La covarianza se informa en el manuscrito y en revision_cientifica/diagnosticos.json."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parent
FIGS = AQUI / "figuras"

# Paleta validada contra la superficie del articulo (fondo blanco).
AZUL, NARANJA, GRIS, ROJO = "#2a78d6", "#eb6834", "#52514e", "#d03b3b"

# Amortiguacion de cada campana del barrido. La clave es el nombre de la
# campana; el valor, el lambda_tau con que se corrio.
DAMP = {
    "sgc_q_mascara_rayos": 0.15,
    "sgc_q_mascara_rayos_damp050": 0.50,
    "sgc_q_mascara_rayos_damp100": 1.00,
    "sgc_q_mascara_rayos_damp200": 2.00,
    "sgc_q_mascara_rayos_damp400": 4.00,
    "sgc_q_mascara_rayos_damp800": 8.00,
}
ADOPTADA = "sgc_q_mascara_rayos_damp400"
HEREDADA = "sgc_q_ampliado_300km"
PREVIA = "sgc_final_q_extra_station_tomography"


def main() -> None:
    f = RAIZ / "descompone_ajuste.json"
    if not f.exists():
        sys.exit("falta descompone_ajuste.json: correr descompone_ajuste.py")
    d = {x["campana"]: x for x in json.loads(f.read_text(encoding="utf-8"))}
    falta = [c for c in DAMP if c not in d]
    if falta:
        sys.exit(f"faltan campañas en el JSON: {', '.join(falta)}")

    orden = sorted(DAMP, key=lambda c: DAMP[c])
    x = [d[c]["razon_estructura_sobre_dato"] for c in orden]
    y = [d[c]["explica_pct"] for c in orden]
    lam = [DAMP[c] for c in orden]

    FIGS.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)

    ax.plot(x, y, "-", color=AZUL, lw=2, zorder=2)
    ax.plot(x, y, "o", color=AZUL, ms=7, zorder=3)
    # El articulo va en ingles: punto decimal, no coma. Y 0.15 y 0.50 caen
    # casi en el mismo punto, asi que sus etiquetas se separan a mano.
    desvios = {0.15: (7, 5), 0.50: (7, -13)}
    for xi, yi, li in zip(x, y, lam):
        ax.annotate(f"{li:g}", (xi, yi), textcoords="offset points",
                    xytext=desvios.get(li, (7, -13)),
                    fontsize=8.5, color=GRIS)

    # La esquina adoptada.
    ia = orden.index(ADOPTADA)
    ax.plot(x[ia], y[ia], "o", ms=13, mfc="none", mec=ROJO, mew=2, zorder=4)
    ax.annotate("adopted", (x[ia], y[ia]), textcoords="offset points",
                xytext=(10, 10), fontsize=9.5, color=ROJO, weight="bold")

    # Las dos referencias externas al barrido.
    # La campana heredada esta en el borde derecho: su etiqueta va hacia
    # dentro, o se sale del lienzo.
    for camp, etq, dx in ((HEREDADA, "expanded network,\ninherited domain", -8),
                          (PREVIA, "reference\nconfiguration", 8)):
        if camp not in d:
            continue
        ax.plot(d[camp]["razon_estructura_sobre_dato"], d[camp]["explica_pct"],
                "s", color=NARANJA, ms=8, zorder=3)
        ax.annotate(etq, (d[camp]["razon_estructura_sobre_dato"],
                          d[camp]["explica_pct"]),
                    textcoords="offset points", xytext=(dx, -4),
                    fontsize=8.5, color=NARANJA, va="top",
                    ha="right" if dx < 0 else "left")

    ax.axvline(1.0, color=GRIS, ls=":", lw=1, zorder=1)
    ax.annotate("unit component scale", (1.0, ax.get_ylim()[0]),
                textcoords="offset points", xytext=(4, 6),
                fontsize=8, color=GRIS, rotation=90, va="bottom")

    ax.set_xlabel(r"$\sigma_{\mathrm{structure}} / \sigma_{\mathrm{data}}$")
    ax.set_ylabel("normalised RMS reduction (%)")
    ax.grid(alpha=.22, lw=.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    salida = FIGS / "curva_l.pdf"
    fig.savefig(salida)
    fig.savefig(salida.with_suffix(".png"), dpi=200)
    plt.close(fig)
    print(f"  {salida.name}")
    print(f"  barrido: " + ", ".join(
        f"λ={l:g} → razón {xi:.2f}, expl {yi:.1f} %"
        for l, xi, yi in zip(lam, x, y)))


if __name__ == "__main__":
    main()
