# Fase 2h (en curso): análisis suficiente del híbrido denso

*Nota de estado, 2026-10-26. Todavía no es un informe.*

## Hecho y validado

- **`fondo.py`: zona lejana canónica.**
  - Para $\ell>\ell_m$ se prescribe $\nu=\Phi'r'-r''$, mezclando con un escalón $C^\infty$ el valor interior y el de un escalar canónico, $\kappa r\phi'^2/2$.
  - $N$ y $r$ salen de T1 escrita como sistema de segundo orden.
  - Un disparo en la amplitud $\beta_\nu$ de una joroba en la zona de mezcla elimina el modo del lapso que crece linealmente, $N\simeq A+B\ell$.
  - Semilla 2 de la Fase 2g: $\beta_\nu=-951$, $M_{\rm ADM}=m_{MS}=5.357$, $N_\infty=0.9988$, y $y=\nu\ell^3\to23$–$26>0$ a partir de $\ell\approx400$ (cola canónica, $K_{\rm nat}>0$). Las condiciones necesarias de la Fase 2g siguen cumpliéndose, con margen 0.09.
- **`teoria.py`: reconstrucción completa con fluido.**
  - $K$ de diseño: ecuación de Riccati dentro de la estrella y ventana de escalón exacto, porque $K_{\rm nat}$ tiene polos en $\ell\approx89$.
  - $G_{2XX}$ a partir de $P$ objetivo.
  - `validar_teoria.py`: ecuaciones del fondo $\lesssim10^{-13}$; $K$ (T2) frente a $K$ de diseño, $10^{-10}$; continuidad $C^3$ en $\ell_m$; $K>0$ en todas partes y condición sin fantasma con fluido cumplida (mín. $4\times10^{-7}$).
- **`espectro2.py`: elementos finitos de Hermite.** Dos campos $(h,c)$ dentro del fluido y uno fuera, con $W=0$.
  - Validación (`validar_espectro.py`): en la garganta de referencia de la Fase 2f reproduce los autovalores pares de Sturm–Liouville, $-0.1174968$ y $-0.0048524$, frente a $-0.11750$ y $-0.00485$.

## Resultado provisional: el híbrido de la semilla 2 no es estable con este diseño

- El espectro completo tiene taquiones grandes ($w\approx-6.3,\,-4.9,\,-4.4,\dots$), iguales en las dos paridades. Viven fuera de la garganta.
- **Diagnóstico** (forma de Sturm–Liouville de la región sin fluido):
  - justo fuera de la superficie de la estrella ($\ell\approx29$–32) $U<0$, con $U/K\approx-6$;
  - también en $\ell\approx100$–150 ($U/K\approx-4$).
- **Causa:** la estrella consume casi todo el presupuesto de $K$ ($B\approx0.05$). Fuera de ella $K\sim10^{-4}$–$10^{-7}$, y cualquier $U<0$ da un taquión rápido.
- Además $T$ queda casi singular (mín. $\sim10^{-8}$, con signo inestable al refinar la malla).

## Próximos pasos

1. Repartir el presupuesto: más $K$ fuera de la estrella (perfil $q$), o más masa ADM relativa a la materia.
2. Añadir a la optimización la condición $U\ge0$ (o $U/K$ acotado) en la región entre la estrella y la zona lejana. Antes hay que comprobar si la $U$ de esa región depende de funciones libres ($F_{4XX}$, $G_3$) o solo del fondo.
3. Repetir el espectro completo sobre las nuevas candidatas, y después las velocidades de propagación y el sector impar.
