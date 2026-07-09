# SRG Toolbox matlab

**Scaled Relative Graph Toolbox for MATLAB.**

Authors: Eder Baron-Prada · Adolfo Anta · Thomas Chaffey · Alberto Padoan · Version 1.0 (2026)

![MATLAB R2020b+](https://img.shields.io/badge/MATLAB-R2020b%2B-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green)

---

## Requirements

- MATLAB R2020b+
- Control System Toolbox
- `matlab2tikz` — optional, for `.tikz` export

---

## Installation

```matlab
cd srg-toolbox
run setup       % adds toolbox root and bg/ to path
help Contents   % list all functions
```

---

## Quick start

```matlab
G = tf([1 2], [1 3 2]);
[W, A, gplus, gminus, rawdata] = srg_compute(G, -3, 3, 200, 32);
figure; srg_plot_gauss(G, gplus, gminus, 1, 200);
```


## Citation

```bibtex
@software{Baron_SRG_Toolbox_2026,
  title  = {{SRG Toolbox}: Scaled Relative Graph Analysis for {MATLAB}},
  author = {Baron-Prada, Eder and Anta, Adolfo and Chaffey, Thomas and Padoan, Alberto},
  year   = {2026},
  url    = {https://github.com/eder-baron/srg-toolbox}
}
```
---

## Function reference

| Function | Description |
|---|---|
| `srg_compute` | Frequency-wise SRG — primary entry point |
| `srg_hard` | Hard SRG (unstable / NMP systems) |
| `srg_scaled_graph` | Soft SRG via hyperbolic convex hull |
| `srg_homotopy` | SRG of τ·G across a sweep of τ values |
| `srg_stability_margin` | Frequency-wise min distance between two SRGs |
| `srg_bounding_circle` | Minimum enclosing circle of an SRG set |
| `srg_qsr_multiplier` | Solver-free (Q,S,R) multiplier fit |
| `srg_plot_gauss` | 2-D Gauss plane plot |
| `srg_plot_beltrami` | 2-D Beltrami–Klein disk plot |
| `srg_plot_compare` | Multi-SRG overlay |
| `srg_plot_3d` | 3-D wire/filled plot |
| `srg_plot_3d_surface` | 3-D smooth surface |
| `srg_plot_3d_compare` | Two SRG surfaces with intersection |
| `srg_plot_3d_beltrami` | 3-D BK-disk plot |
| `srg_plot_3d_beltrami_compare` | Two SRG surfaces in BK disk |
| `srg_plot_witness_3d` | Min-distance witness segments (3-D) |
| `srg_style` / `srg_apply_style` | Style configuration (default / dark / publication) |
| `srg_export` | Export: PNG, TikZ, PDF/EPS |

---

## Examples

| Script | Demonstrates |
|---|---|
| `example_01_freqwise_srg.m` | Frequency-wise SRG, Gauss/BK/3-D plots, stability margin |
| `example_02_soft_srg.m` | Soft SRG for SISO and MIMO pairs |
| `example_03_hard_srg.m` | Hard SRG with unstable poles |
| `example_04_qsr_export.m` | QSR multiplier + publication export |

---



---

## License & attribution

MIT License. See [LICENSE](LICENSE).

| File | Author | License |
|---|---|---|
| `bg/field_of_values.m` | N. J. Higham | BSD-style |
| Algorithm in `srg_hard.m` | J. P. J. Krebbekx (`SrgTools.jl`) | BSD 3-Clause |
