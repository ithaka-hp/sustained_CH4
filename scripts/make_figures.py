"""
Regenerate Figures 1-4 of the paper into ../figures/ (PNG and PDF).

    python scripts/make_figures.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import model as m
import offsetting as off

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
os.makedirs(OUT, exist_ok=True)
C_SINK = "#3b6ea5"; C_NEW = "#e08214"; C_WARM = "#1a1a1a"; C_DEF = "#c0392b"; C_SUR = "#2e8b7f"

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight", facecolor="white")
    plt.close(fig)

def figure1():
    Y = 100; n = m.years(Y)
    mrem, _ = off.program(1, Y); W = off.flux_warming(Y)
    new_cool = mrem * float(m.irf_year(1)); earlier = W - new_cool; cum = np.cumsum(mrem)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.0, 6.2), dpi=170, sharex=True,
                                 gridspec_kw=dict(height_ratios=[2, 1.1], hspace=0.14))
    a1.fill_between(n, 0, earlier, color=C_SINK, alpha=.85, label="cooling from sinks of earlier years")
    a1.fill_between(n, earlier, W, color=C_NEW, alpha=.95, label="cooling from removal added this year")
    a1.plot(n, W, color=C_WARM, lw=1.6, label="annual warming of the flux")
    a1.set_ylim(0, 1500); a1.set_ylabel("Annual climate effect (t CO$_2$e)")
    a1.legend(loc="center right", fontsize=8.5, framealpha=.95)
    a1.set_title("(a) each year's warming is cancelled by the cooling of all sinks in place", loc="left", fontsize=10)
    a2.fill_between(n, 0, mrem, color=C_NEW, alpha=.9, step="mid")
    a2.set_ylim(0, 125); a2.set_ylabel("Removal added\n(t CO$_2$e per year)", fontsize=9); a2.set_xlabel("Year")
    a2b = a2.twinx(); a2b.plot(n, cum, color="#555", lw=1.4); a2b.set_ylim(0, 3300)
    a2b.set_ylabel("Cumulative removal\n(t CO$_2$e)", color="#555", fontsize=9); a2b.tick_params(axis="y", colors="#555")
    a2.set_title("(b) new removal added each year falls as earlier sinks accumulate", loc="left", fontsize=10)
    a2.set_xlim(1, Y)
    save(fig, "figure1_principle")

def figure2():
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.4, 8.2), dpi=170, sharex=True,
                                 gridspec_kw=dict(height_ratios=[1, 1], hspace=0.22))
    # (a) single emission
    Y = 100; n, D, A, V = off.single_register(1, Y); warm = m.ace_ch4(n)
    a1.fill_between(n, warm, D, color=C_NEW, alpha=.85, label="available cooling (for later emissions)")
    a1.plot(n, D, color=C_SINK, lw=1.9, label="delivered cooling of the sink")
    a1.plot(n, warm, color=C_WARM, lw=1.9, label="emission's warming = assigned cooling")
    a1.set_ylim(0, 120); a1.set_ylabel("Annual climate effect (t CO$_2$e)")
    a1.legend(loc="upper right", fontsize=8.4, framealpha=.95)
    a1.set_title("(a) single emission", loc="left", fontsize=10.5, fontweight="bold")
    # (b) sustained flux register, 50 yr, bars coloured by recency within each year
    Yb = 50; mrem, _ = off.program(1, Yb); W = off.flux_warming(Yb); maxage = Yb - 1
    col_age = lambda age: plt.cm.Blues(0.95 - (age - 1) / maxage * (0.95 - 0.30))
    for yr in range(1, Yb + 1):
        bottom = 0.0
        for v in range(1, yr + 1):
            h = mrem[v - 1] * float(m.irf_year(yr - v + 1))
            col = C_NEW if v == yr else col_age(yr - v)
            a2.bar(yr, h, bottom=bottom, width=0.86, color=col, linewidth=0, zorder=2); bottom += h
    a2.step(np.arange(1, Yb + 1) + 0.5, W, where="post", color="k", lw=1.1, zorder=4)
    handles = [Patch(color=C_NEW, label="C-sink added this year"),
               Patch(color=col_age(1), label="cooling from former C-sinks (darker = more recently added)"),
               plt.Line2D([0], [0], color="k", lw=1.1, label="warming of the methane impulses (matched)")]
    a2.legend(handles=handles, loc="lower right", fontsize=8.4, framealpha=.96)
    a2.set_ylim(0, 1420); a2.set_ylabel("Annual climate effect (t CO$_2$e)"); a2.set_xlabel("Year"); a2.set_xlim(0.4, Yb + 0.6)
    a2.set_title("(b) sustained flux, as a register", loc="left", fontsize=10.5, fontweight="bold")
    save(fig, "figure2_register")

def figure3():
    Y = 300; n = m.years(Y); W = off.flux_warming(Y); Cl = off.gwp_program(m.GWP100, Y)
    fig, ax = plt.subplots(figsize=(8.0, 4.5), dpi=170)
    ax.fill_between(n, Cl, W, where=(W >= Cl), color=C_DEF, alpha=.25, label="net warming (deficit)")
    ax.fill_between(n, W, Cl, where=(Cl > W), color=C_SUR, alpha=.25, label="net cooling (surplus)")
    ax.plot(n, W, color=C_WARM, lw=1.8, label="warming of the methane flux")
    ax.plot(n, Cl, color="#2c6fbb", lw=1.8, label="cooling of the GWP100 removal (27 t yr$^{-1}$)")
    ax.axvline(100, color="#777", ls=":", lw=1.2)
    ax.annotate("balance only\nat year 100", xy=(100, W[99]), xytext=(120, 760), fontsize=8.6,
                arrowprops=dict(arrowstyle="->", color="#555"))
    ax.set_xlim(1, Y); ax.set_ylim(0, 3000); ax.set_xlabel("Year"); ax.set_ylabel("Annual climate effect (t CO$_2$e)")
    ax.legend(loc="upper left", fontsize=8.5, framealpha=.95)
    save(fig, "figure3_gwp100")

def figure4():
    Y = 100; n = m.years(Y); mrem, _ = off.program(1, Y)
    g100 = m.GWP100; g20 = m.GWP20_AR6; gstar = off.gwpstar_steady()
    fig, ax = plt.subplots(figsize=(8.0, 4.5), dpi=170)
    ax.plot(n, mrem, color=C_NEW, lw=2.2, label="global warming-matched offsetting (annual matching)")
    ax.axhline(g20, color="#555", ls="--", lw=1.4, label=f"GWP20 ({g20:.1f} t)")
    ax.axhline(g100, color="#555", ls="-.", lw=1.4, label=f"GWP100 ({g100:.0f} t)")
    ax.axhline(gstar, color="#999", ls=":", lw=1.8, label=f"GWP* residual, steady flux ({gstar:.1f} t)")
    ax.set_xlim(1, Y); ax.set_ylim(-4, 122); ax.set_xlabel("Year")
    ax.set_ylabel("Permanent removal added (t CO$_2$e per t CH$_4$)")
    ax.legend(loc="upper right", fontsize=8.5, framealpha=.95)
    save(fig, "figure4_metrics")

if __name__ == "__main__":
    figure1(); figure2(); figure3(); figure4()
    print("figures written to", OUT)
