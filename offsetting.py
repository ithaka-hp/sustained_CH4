"""
Global warming-matched offsetting: size permanent CO2 removal so that the cooling
it delivers each year equals the warming a methane emission or flux exerts that
year. Builds on model.py.

Functions
---------
removal_single(H)       permanent removal that offsets one tonne of methane over an H-year window
program(H, T)           per-year removal for a sustained 1 t CH4/yr flux, with its delivered cooling
flux_warming(T)         annual warming of a sustained 1 t CH4/yr flux
gwp_program(rate, T)    cooling delivered by a constant rate of permanent removal (e.g. GWP100 at 27 t/yr)
gwpstar_steady()        the GWP* residual for a steady flux (stock term, a quarter of GWP100)
single_register(H, T)   delivered / assigned / available cooling of one sink for one emission (Figure 2a)
"""
import numpy as np
import model as m

def removal_single(H=1):
    """Permanent removal (t CO2e) that offsets 1 t CH4 over an H-year balance window: TCE_CH4(H)/IRFint(H)."""
    return m.tce_ch4(H) / float(m.irfint(H))

def flux_warming(T=1000):
    """Annual warming of a sustained flux of 1 t CH4/yr in years 1..T (t CO2e). Equals cumulative ACE_CH4."""
    return np.cumsum(m.ace_ch4(m.years(T)))

def _unit_cooling(T):
    """Cooling delivered in years 1..T by one tonne of permanent removal placed at the start of year 1 (positive)."""
    return m.irf_year(m.years(T))

def program(H=1, T=1000):
    """
    Global warming-matched offsetting of a sustained 1 t CH4/yr flux.

    Each year's emission (vintage) is offset under the same rule: within its first H
    years the cooling assigned to it must sum to its H-year total climate effect; from
    year H+1 the cooling assigned must equal its warming that year. The surplus cooling
    of earlier sinks is assigned before any new removal is bought, so the removal added
    each year is only the shortfall. H = 1 is annual matching (the paper's default).

    Returns
    -------
    mrem : array, removal added in each year 1..T (t CO2e)
    cool : array, total cooling delivered in each year 1..T (t CO2e)
    """
    HH = T + 200
    u = _unit_cooling(HH)                                   # u[a] = cooling at age a+1 of a unit sink
    f = np.array([float(m.ch4int(k + 1) - m.ch4int(k)) for k in range(HH)])  # methane warming at age a+1
    uL = u[:H].sum()
    pool = np.zeros(HH); mrem = np.zeros(T); cool = np.zeros(T)
    for v in range(T):
        own = np.zeros(HH)
        y0, y1 = v, min(v + H, HH)
        need = f[:y1 - v].sum()
        short = need; taken = 0.0
        for t in range(y0, y1):
            take = min(pool[t], max(0.0, short - taken)); pool[t] -= take; taken += take
        short -= taken
        if short > 1e-9:
            mb = short / uL; mrem[y0] += mb if y0 < T else 0.0; own[y0:] += mb * u[:HH - y0]
        for t in range(v + H, HH):
            a = t - v; net = own[t] - f[a]
            if net >= 0:
                pool[t] += net
            else:
                take = min(pool[t], -net); pool[t] -= take; net += take
                if net < -1e-9:
                    mt = -net / u[0]; own[t:] += mt * u[:HH - t]
                    if t < T: mrem[t] += mt
        cool += own[:T]
    return mrem, cool

def gwp_program(rate, T=1000):
    """Cooling delivered in years 1..T by a constant `rate` of permanent removal bought every year."""
    u = _unit_cooling(T); cool = np.zeros(T)
    for j in range(T):
        cool[j:] += rate * u[:T - j]
    return cool

def gwpstar_steady(gwp100=m.GWP100):
    """GWP* removal for a steady flux: the stock term only, s * GWP100 with s = 0.25 (Lynch et al., 2020)."""
    return 0.25 * gwp100

def single_register(H=1, T=1000):
    """
    Register of one permanent sink that offsets one tonne of methane (Figure 2a).

    Returns year array and delivered, assigned (= methane warming while covered) and
    available (= delivered - assigned) cooling, all in t CO2e.
    """
    n = m.years(T); M = removal_single(H)
    delivered = M * m.irf_year(n)
    warming = m.ace_ch4(n)
    assigned = np.minimum(delivered, warming)
    available = delivered - assigned
    return n, delivered, assigned, available

if __name__ == "__main__":
    for H in (1, 5, 20, 100):
        print(f"H={H:3d}: single-emission removal = {removal_single(H):6.1f} t CO2e")
    mrem, cool = program(1)
    W = flux_warming(len(mrem))
    print("sustained flux, H=1: removal yr 1/20/100 =",
          [round(float(mrem[k - 1]), 1) for k in (1, 20, 100)],
          "| cumulative 20/100/1000 =", [round(float(np.cumsum(mrem)[k - 1])) for k in (20, 100, 1000)])
