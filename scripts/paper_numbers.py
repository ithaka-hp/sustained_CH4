"""
Reproduce every number quoted in the paper, grouped by section.

    python scripts/paper_numbers.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import model as m
import offsetting as off

def r(x, d=1): return round(float(x), d)

print("== 2.1 Measures ==")
print("  IRF_1, IRF_20, IRF_100 :", r(m.irf_year(1), 3), r(m.irf_year(20), 3), r(m.irf_year(100), 3))
print("  IRFint(20), IRFint(100) yr :", r(m.irfint(20), 2), r(m.irfint(100), 2))
print("  ACE_CH4 yr1 / yr20 / yr61 :", r(m.ace_ch4(1)), r(m.ace_ch4(20)), r(m.ace_ch4(61)))
print("  TCE_CH4(100) :", r(m.tce_ch4(100), 0))

print("== 2.2 Global warming-matched offsetting ==")
print("  single-emission removal, H=1 / H=20 / H=100 :",
      r(off.removal_single(1)), r(off.removal_single(20)), r(off.removal_single(100)))
print("  (H=100 equals GWP100 =", r(m.gwp(100)), ")")

print("== 3.1 Sustained flux under GWP100 ==")
Y = 500; n = m.years(Y); W = off.flux_warming(Y); Cl = off.gwp_program(m.GWP100, Y)
net = W - Cl
print("  annual cooling/warming at yr 100 :", r(Cl[99] / W[99], 3))
print("  max annual deficit :", r(net.max(), 0), "t CO2e at yr", int(n[np.argmax(net)]))
print("  annual surplus at yr 200 / 500 :", r(-net[199], 0), r(-net[499], 0))
cumW = np.cumsum(W); cumC = np.cumsum(Cl)
print("  cumulative cooling/warming over 20 / 50 yr :",
      f"{cumC[19] / cumW[19] * 100:.0f}%", f"{cumC[49] / cumW[49] * 100:.0f}%")

print("== 3.2 Rule vs potentials (per t CH4) ==")
mrem, _ = off.program(1, 1000); cum = np.cumsum(mrem)
print("  removal yr 1 / 20 / 100 :", [r(mrem[k - 1]) for k in (1, 20, 100)])
print("  cumulative 20 / 100 / 1000 :", [r(cum[k - 1], 0) for k in (20, 100, 1000)])
print("  GWP100 / GWP20(AR6) / GWP* steady :", m.GWP100, m.GWP20_AR6, r(off.gwpstar_steady()))

print("== 3.3 Global livestock (112 Tg CH4/yr) ==")
F = 112e6
print("  removal yr 1 / 20 / 100 (Gt) :", [r(mrem[k - 1] * F / 1e9) for k in (1, 20, 100)])
print("  cumulative over century (Gt) :", r(cum[99] * F / 1e9, 0))
print("  GWP100 constant (Gt/yr) :", r(m.GWP100 * F / 1e9), " over century (Gt) :", r(m.GWP100 * F / 1e9 * 100, 0))

print("== abstract check ==")
print("  per-tonne removal 114/50/11 :", [r(mrem[k - 1]) for k in (1, 20, 100)])
print("  offset over first 20 / 50 yr :", f"{cumC[19] / cumW[19] * 100:.0f}%", f"{cumC[49] / cumW[49] * 100:.0f}%")
