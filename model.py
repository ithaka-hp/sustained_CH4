"""
Core forcing model for "Offsetting sustained methane emissions: replacing GWP100
accounting by global warming-matched offsetting" (Schmidt and Hagemann).

All quantities are per tonne of the gas named. Masses are in t CO2e; the annual
climate effect (ACE) is in t CO2e; the total climate effect (TCE) is in t CO2e.yr.

Conventions
-----------
Year n runs from time n-1 to n (n = 1, 2, ...). An emission or removal is placed
at the start of year 1. ACE(n) is the mean airborne climate effect during year n;
TCE(H) is the sum of ACE over years 1..H.

The CO2 impulse-response function is the five-term form of Jeltsch-Thommes and
Joos (2020), as applied by Schmidt and Hagemann (2026). Biogenic methane is scaled
so that its 100-year total climate effect equals its 100-year global warming
potential (27.0) times the 100-year integrated response; its perturbation lifetime
is 11.8 yr.
"""
import numpy as np

# --- CO2 impulse-response function (Jeltsch-Thommes and Joos, 2020) ---
A   = [0.008, 0.044, 0.112, 0.224, 0.310, 0.297]      # amplitudes (a0 is the persistent term)
TAU = [np.inf, 68521., 5312., 362., 47., 6.]          # time constants (yr)

def irf(t):
    """Airborne fraction of a unit CO2 pulse a time t (yr) after release."""
    t = np.asarray(t, float)
    out = np.full(np.shape(t), A[0])
    for a, tau in zip(A[1:], TAU[1:]):
        out = out + a * np.exp(-t / tau)
    return np.where(t < 0, 0.0, out)

def irfint(t):
    """Time-integrated IRF from 0 to t (yr); units: yr."""
    t = np.asarray(t, float)
    out = A[0] * t
    for a, tau in zip(A[1:], TAU[1:]):
        out = out + a * tau * (1 - np.exp(-t / tau))
    return np.where(t < 0, 0.0, out)

# --- methane ---
GWP100  = 27.0        # biogenic methane, 100-year global warming potential (Forster et al., 2021)
TAU_CH4 = 11.8        # methane perturbation lifetime (yr)
TCE100  = GWP100 * float(irfint(100.0))                        # 1363 t CO2e.yr per t CH4
F0      = TCE100 / (TAU_CH4 * (1 - np.exp(-100 / TAU_CH4)))    # 115.5, initial forcing in CO2e burden units

def ch4int(t):
    """Time-integrated forcing of a unit methane pulse from 0 to t (yr)."""
    t = np.asarray(t, float)
    return np.where(t < 0, 0.0, F0 * TAU_CH4 * (1 - np.exp(-t / TAU_CH4)))

# --- annual measures (year n = 1, 2, ...) ---
def years(N):
    return np.arange(1, N + 1)

def irf_year(n):
    """IRF_n: mean airborne fraction of a CO2 pulse during its n-th year (0 for n < 1)."""
    n = np.asarray(n, float)
    out = irfint(np.maximum(n, 0)) - irfint(np.maximum(n - 1, 0))
    return np.where(n >= 1, out, 0.0)

def ace_co2(M, n):
    """Annual climate effect in year n of a CO2 mass M (positive emission, negative removal)."""
    return M * irf_year(n)

def tce_co2(M, H):
    """Total climate effect to horizon H of a CO2 mass M; equals M * irfint(H)."""
    return M * float(irfint(H))

def ace_ch4(n):
    """ACE of 1 t biogenic methane in its n-th year. Equals 110.8 * exp(-(n-1)/11.8)."""
    n = np.asarray(n, float)
    out = ch4int(np.maximum(n, 0)) - ch4int(np.maximum(n - 1, 0))
    return np.where(n >= 1, out, 0.0)

def tce_ch4(H):
    """Total climate effect of 1 t biogenic methane over years 1..H (t CO2e.yr)."""
    return float(ace_ch4(years(int(H))).sum())

def gwp(H):
    """H-year global warming potential of biogenic methane in this model: TCE_CH4(H) / IRFint(H)."""
    return tce_ch4(H) / float(irfint(H))

# IPCC AR6 (Forster et al., 2021) 20-year GWP of non-fossil methane, used in Figure 4.
GWP20_AR6 = 79.7

if __name__ == "__main__":
    print(f"IRF_1={float(irf_year(1)):.3f}  IRF_20={float(irf_year(20)):.3f}  IRF_100={float(irf_year(100)):.3f}")
    print(f"IRFint(20)={float(irfint(20)):.2f} yr  IRFint(100)={float(irfint(100)):.2f} yr")
    print(f"ACE_CH4 yr1={float(ace_ch4(1)):.1f}  yr20={float(ace_ch4(20)):.1f}  TCE_CH4(100)={tce_ch4(100):.0f}")
    print(f"GWP5={gwp(5):.1f}  GWP20(model)={gwp(20):.1f}  GWP100(model)={gwp(100):.1f}  GWP20(AR6)={GWP20_AR6}")
