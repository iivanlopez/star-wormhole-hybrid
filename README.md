# star-wormhole-hybrid: monopole instability of wormhole throats

Code, data and reports accompanying the paper

> Iván López Martín, *The $d=N$ deficit and the throat ghost: a common mechanism for the monopole instability of wormholes in second-order scalar- and vector-tensor theories* (2026).

The paper is in [`fase2n/articulo/`](fase2n/articulo/):
- English version: `wormhole_monopole.tex` / `.pdf`.
- Spanish version: `gargantas_monopolo.tex` / `.pdf`.

The reports of each phase (`faseXX/faseXX_resultados.md`) are written in Spanish.

## Structure

Each phase folder contains:
- the scripts;
- their output logs (`*_log.txt`);
- numerical results (`resultados/`) and figures (`figuras/`);
- a report `faseXX_resultados.md` whose *Archivos* section describes each script; later phases also have a *Reproducibilidad* section listing the commands in order.

| phase | topic | section of the paper |
|---|---|---|
| `fase1` | theoretical framework | — |
| `fase2`–`fase2e` | star–throat hybrid in Einstein–scalar–Gauss–Bonnet (closed) | introduction |
| `fase2f` | quartic GLPV: symbolic machinery, monopole structure | §2 |
| `fase2g` | $d=N$ trial, charge, perfect fluid | §3, §5.1–5.2 |
| `fase2h` | fluid with a free surface | §5.2 |
| `fase2i` | elastic solid | §5.2 |
| `fase2j` | quintic Horndeski, BCK throats | §3, §5.3 |
| `fase2k` | scalar linear in time $\phi=qt+\psi$, cubic galileon with flux | §4.1, §6 |
| `fase2l` | slow rotation | §5.1 |
| `fase2m` | generalized Proca | §3, §4.2 |
| `fase2n` | common sign rule, lemmas of §6, paper | §4.3, §6 |

## Reproducing

Requires Python 3.11 and the packages in `requirements.txt`:

```bash
pip install -r requirements.txt
```

Run the scripts from inside each phase folder, in the order given in its report (*Reproducibilidad* or *Archivos*). Some symbolic generators take from minutes to about an hour. They write intermediate caches (`resultados/*.pkl`) that the later scripts read. These caches are included, but they can be regenerated: as with any pickle, regenerate them rather than load them if you do not trust their origin.

## License

The code is released under the MIT license (see `LICENSE`). The text of the paper and the reports is © the author.

Symbolic and numerical computations were carried out with the assistance of Claude (Anthropic).
