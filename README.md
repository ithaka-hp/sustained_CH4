# sustained_CH4

Code accompanying the manuscript **"Offsetting sustained methane emissions: replacing GWP100 accounting by global warming-matched offsetting"** (Schmidt and Hagemann, in preparation).

It reproduces every figure and every number in the paper from a single forcing model.

## Method

A sustained biogenic methane flux is offset today by converting it to CO2 equivalents with the 100-year global warming potential (GWP100 = 27 t CO2e per t CH4) and cancelling it with the same mass of permanent CO2 removal, counted once. For a sustained flux this is mistimed: the removal delivers cooling equal to the flux's warming only in year 100, under-cooling before and over-cooling after.

**Global warming-matched offsetting** sizes permanent removal so that the cooling it delivers each year equals the warming the flux exerts that year, both measured on the CO2 impulse-response function as an **annual climate effect (ACE)**. A register records, for each sink and year, the cooling delivered, the cooling assigned to an emission and the cooling left available, so that one permanent sink is credited to a succession of emissions, each in the year its cooling meets the warming.

## Install

```
pip install -r requirements.txt
```

## Reproduce the paper

```
python scripts/paper_numbers.py     # prints every number quoted in the paper
python scripts/make_figures.py      # writes Figures 1-4 to figures/ (PNG and PDF)
```

## Key parameters

- CO2 impulse-response function: five-term form of Jeltsch-Thommes and Joos (2020), as applied by Schmidt and Hagemann (2026).
- Biogenic methane: 100-year global warming potential 27.0 (Forster et al., 2021), perturbation lifetime 11.8 yr.
- 20-year GWP of non-fossil methane used in Figure 4: 79.7 (Forster et al., 2021).
- Balance window H: annual matching (H = 1) throughout the paper; H is a free parameter in the code (H = 100 reproduces the GWP100 mass of 27 t).

## Files

- `model.py` - the CO2 impulse-response function and the annual (ACE) and total (TCE) climate-effect measures for CO2 and methane.
- `offsetting.py` - global warming-matched offsetting: single-emission sizing, the sustained-flux program, the single-sink register, and the GWP100, GWP20 and GWP* comparisons.
- `scripts/make_figures.py` - Figures 1-4.
- `scripts/paper_numbers.py` - all paper numbers.
- `figures/` - generated figures.

## Related

The temporary-sink rule this work builds on is in [ithaka-hp/CH4-offset](https://github.com/ithaka-hp/CH4-offset) (Schmidt and Hagemann, 2026, GCB Bioenergy).

## Development

The code was developed with the assistance of Claude (Anthropic) and reviewed and validated by the authors; all numbers and figures are reproduced by the scripts above.

## License

MIT (see `LICENSE`). If you use this code, please cite the paper and this repository (`citation.cff`).
