# Híbrido estrella–agujero de gusano — Fase 2n: el teorema común y el artículo

*Documento de trabajo v0.1 · 2026-10-27 · continúa [fase2k/fase2k_resultados.md](../fase2k/fase2k_resultados.md) (§15) y [fase2m/fase2m_resultados.md](../fase2m/fase2m_resultados.md) (§4.5)*

## Resumen

**Parte 1 (teorema común, plazo acotado).** El objetivo era demostrar que toda garganta estática, simétrica y asintóticamente plana de una teoría escalar o vectorial–tensorial de segundo orden con simetría de desplazamiento tiene un taquión o un fantasma en $l=0$. **No cierra en el plazo**: el galileón cúbico con flujo ($G_3\neq0$, $\psi'\neq0$) sigue sin teorema. Lo que sí queda establecido:

1. **Enunciado común de lo ya demostrado** (§1), con hipótesis precisas, en dos ramas:
   - *campo radial* (escalar estático en GLPV cuártica u Horndeski quíntico, Proca con $A_0=0$): $E[N]=-2\int N^3\mathcal G\sigma^2<0$ (2g, 2j, 2m);
   - *campo temporal sin flujo* ($\phi=qt$ o $G_3=0$; Proca con $A_1=0$, $G_3=0$): fantasma en la garganta (2k, 2m).
2. **La regla de signos es una sola fórmula** (§2, nuevo, exacto, `regla_comun.py`). Con $A_4=-\mathcal G$, las dos acciones cinéticas se escriben
   $$T=\int\Big[-\frac{2\mathcal G}N\dot d^2+k\Big(\Big(\frac{\mathcal G\dot d}N\Big)'\Big)^2\Big],\qquad k(0)=-\frac{\kappa_0}{\mathcal G(0)\,r''(0)},$$
   con $\kappa_0=2N_0r_0$ (escalar, $k=4r/(Qr')'$) y $\kappa_0=4N_0^3Xr_0/U_0^2$ (vector, $k=k_2$). El abocinamiento $r''(0)>0$ da a los dos términos el signo de $-\mathcal G(0)$.
3. **Lemas nuevos para el caso cúbico con flujo** (§3, exactos, `lemas.py`):
   - **(L4) la garganta es un mínimo de $\beta=\mathcal NG_4$:** si el gradiente del monopolo es sano en la garganta con $A_4<0$, entonces $\beta''(0)>R_2\beta_0/(2R_0)>0$;
   - **(L1–L3) la cola:** $r^3L[r']=(rK)'$ para todo $r$, $(rK)'=-4\pi r'r^3\,d\varepsilon/dr$ con $\varepsilon=m_r/4\pi r^2$, y donde el flujo se desvanece $\mathfrak B\to-2\beta'(rK)'$. Un coeficiente cinético sano exige **densidad efectiva decreciente** donde el flujo se apaga. Eso excluye todas las colas de tipo Ellis, con cualesquiera $A_4$ y lapso;
   - **(L5) obstrucción a la vía (i)** propuesta: $\mathfrak B$ es afín en $Q'$ con coeficiente $r^4w'(r'w'-r''w)$ ($w=v/r$), y ninguna identidad $\int\mu\mathfrak B$ con peso independiente de $Q$ tiene integrando de signo definido.
4. **Numérica nueva** (§4). La clase que dejan abierta los lemas (núcleo de Ellis y cola de energía positiva) **sí** admite $\beta>0$: en la criba, 34 de 600 fondos tienen $eta>0$ en todo el eje, 19 de ellos con el coeficiente cinético sano en toda la cola ($
ho>20$). Esto corrige la Fase 2k §13 ("con colas de energía positiva, $\beta$ cambia siempre de signo"), que era una propiedad de la familia con $\varepsilon>0$ en la garganta. Pero ningún fondo tiene el monopolo localmente sano: $s_2<0$ junto a la garganta y en la zona donde el flujo decae. La evolución diferencial sobre la clase (9 parámetros) encuentra un fondo **hiperbólico en todo el eje**, con $eta>0$ y cola sana, pero con una **capa fantasma** en $0.56\lesssim
ho\lesssim1.7$: allí $s_2<0$ y $s_0>0$ a la vez, el símbolo entero tiene el signo cambiado frente a los gravitones.
5. **Caso mixto de Proca** ($A_0,A_1\neq0$): no tratado; los lemas usan la ecuación de $\beta$ de la GLPV y no se extienden sin la reducción correspondiente.

**Veredicto de la Parte 1.** Queda demostrado el teorema en las dos ramas sin flujo y sin componentes mixtas. El caso cúbico con flujo queda como **conjetura** con: dos lemas que restringen la garganta y la cola, una obstrucción a la estrategia de pesos fijos y unos 3600 fondos (≈1300 de 2k, 600 de la criba y ≈2200 evaluaciones de la optimización de 2n) sin ningún monopolo localmente sano.

**Parte 2 (artículo).** [articulo/gargantas_monopolo.tex](articulo/gargantas_monopolo.tex): las Fases 2f–2m con la tesis del mecanismo común; el caso cúbico con flujo va como conjetura con la evidencia numérica. Bibliografía verificada con la API de arXiv (`arxiv_bib.py`). No hay `pdflatex` en esta máquina: no se ha compilado (§6).

---

## 1. Lo que ya estaba demostrado, en un enunciado

**Hipótesis comunes.** Fondo $ds^2=-N^2dt^2+d\ell^2+r^2d\Omega^2$, $\ell\in\mathbb R$, garganta en $\ell=0$ ($r'(0)=0$, $r''(0)>0$), dos extremos asintóticamente planos ($r\to\infty$, $N\to1$, $G_4\to G_{4\infty}>0$). Teorías de segundo orden con simetría de desplazamiento: GLPV cuártica, Horndeski quíntico, Proca generalizada cuártica. $\mathcal G$ es el coeficiente cinético del sector impar o tensorial: $G_4-2XG_{4X}-4X^2F_4$ (cuártica), el de Kobayashi–Motohashi–Suyama (quíntica), $G_4-2XG_{4X}$ (Proca).

| rama | casos | resultado | fase |
|---|---|---|---|
| radial | escalar estático, cualquier $G_i,F_4$ (y $G_5$); Proca con $A_0=0$ | $E[N]=-2\int N^3\mathcal G\sigma^2<0$; taquión si $\mathcal G>0$, $K,P>0$ | 2g, 2j, 2m |
| temporal sin flujo | $\phi=qt$ ($\psi'=0$) o $A_3'\equiv0$ | $G_4$ cambia de signo (no-go de existencia) | 2k §5 |
| temporal sin flujo | $\phi=qt$; Proca con $A_1=0$, $G_3=0$ | fantasma monopolar si $\mathcal G(0)>0$, tensorial si $\mathcal G(0)<0$ | 2k §6, 2m §4.5 |
| ampliaciones | carga, rotación lenta, fluido, sólido, BCK | el déficit par no cambia de signo (o acoplo fuerte) | 2g, 2h, 2i, 2j, 2l |
| **abierto** | $G_3\neq0$ con $\psi'\neq0$; Proca con $A_0,A_1\neq0$ | sin teorema | 2k II, 2m §6 |

## 2. La regla de signos común (`regla_comun.py`)

En la Fase 2k, $A_4=-G_4+2XG_{4X}+4X^2F_4$. Es exactamente $-\mathcal G$: el coeficiente cinético tensorial del corte $\phi=$ cte es $-A_4=\mathcal G$, el mismo que el del sector impar. Con $Q=2A_4/N=-2\mathcal G/N$:
$$T_{2k}=\int\Big[Q\dot d^2+\frac{r[(Q\dot d)']^2}{(Qr')'}\Big]=\int\Big[-\frac{2\mathcal G}N\dot d^2+\frac{4r}{(Qr')'}\Big(\Big(\frac{\mathcal G\dot d}N\Big)'\Big)^2\Big]\qquad\text{(identidad, simbólica).}$$
Es la forma de la Fase 2m, $T_{2m}=\int[-(2\mathcal G/N)\dot d^2+k_2((\mathcal G\dot d/N)')^2]$. En la garganta:

| | $k(0)$ | $\kappa_0=-k(0)\,\mathcal G(0)\,r''(0)$ |
|---|---|---|
| escalar ($\phi=qt$) | $4r_0/(Q_0r''_0)=-2N_0r_0/(\mathcal G_0r''_0)$ | $2N_0r_0>0$ |
| vector ($A_0$) | $-4r_0^2N_0^4/(U_0^2\tilde c(0))$, $\tilde c(0)=N_0r_0\mathcal G_0r''_0/X$ | $4N_0^3Xr_0/U_0^2>0$ |

**Proposición.** En las dos teorías temporales sin flujo, cerca de una garganta los dos coeficientes cinéticos del monopolo tienen el signo de $-\mathcal G(0)$. Es decir, el abocinamiento convierte la salud tensorial en un fantasma monopolar.

El mecanismo de fondo es el mismo en las dos ramas: la topología de dos extremos. Con campo radial se paga con $\sigma=(\ln r/N)'\ne0$; con campo temporal, con el potencial integrado positivo $\int(1+r'^2)/r^3$ de la ecuación de $\beta$, o con el signo de $r''(0)$ en la energía cinética.

## 3. El caso cúbico con flujo: lemas exactos (`lemas.py`)

Notación de la Fase 2k, Parte II: $\beta=\mathcal NG_4$, $Q=2A_4/\mathcal N<0$, $u=v^2/r^2$, $w=v/r$, $K=1+rr''-r'^2$,
$$L\beta\equiv\Big(\frac{\beta'}r\Big)'+\frac{K\beta}{r^3}=-\frac{(Qr^4u')'}{4r^3},\qquad s_2=-\frac{2\mathfrak B}{rv^2r''^2},\qquad \mathfrak B=-2\beta'(rK)'-Qwr''(r^4w')'+Q'r^4w'(r'w'-r''w).$$
(La forma de $\mathfrak B$ con $w$ se comprueba contra el símbolo exacto `fase2k/resultados/simbolo_perfiles.pkl`.) El coeficiente cinético en el corte $\tau$ es sano si $s_2>0$, es decir $\mathfrak B<0$.

**(L1) Identidad.** Para todo $r(\rho)$, $r^3L[r']=(rK)'$. Por tanto $r'$ es solución homogénea de la ecuación de $\beta$ si y solo si $rK$ es constante. Ejemplo: Schwarzschild en longitud propia, con $rK=3M$ y $\beta=N=r'$.

**(L2) Densidad efectiva.** Con la masa de Misner–Sharp del corte, $r'^2=1-2m(r)/r$:
$$K=\frac{3m}r-m_r,\qquad (rK)'=-r'r^3\frac{d}{dr}\Big(\frac{m_r}{r^2}\Big)=-4\pi r'r^3\frac{d\varepsilon}{dr},\qquad \varepsilon=\frac{m_r}{4\pi r^2}.$$

**(L3) Límite de flujo débil.** Con $(v,v',v'')=\epsilon(a,b,c)$ y $\epsilon\to0$, $\mathfrak B\to-2\beta'(rK)'$ (exacto en el símbolo). Por tanto:
- donde el flujo se apaga (borde de un soporte compacto o cola con flujo de decaimiento rápido), la salud cinética exige $\beta'(rK)'>0$, es decir $\beta'\,d\varepsilon/dr<0$;
- más allá del flujo, $\beta'/r=\int_\rho^\infty K\beta/r^3>0$ si $K,\beta>0$, así que **la densidad efectiva debe decrecer** donde el flujo se apaga;
- **corolario:** las colas de tipo Ellis ($\varepsilon<0$ y creciente hacia 0; en Ellis $(rK)'=-2r'/r^2$) tienen $s_2<0$ allí, con cualesquiera $A_4$ y $\mathcal N$. Es el fallo observado en los 3 fondos hiperbólicos de la Fase 2k §13, ahora como resultado.
- *Precisión a la Fase 2k §13.* La condición del régimen sin flujo, $(Qr')'>0$, sí se cumple en una cola de Ellis si $A_4$ tiene cola propia: por ejemplo $Q=-2-2/\rho^2$ da $(Qr')'\simeq2/\rho^3>0$. Lo que excluye la cola de Ellis es el régimen de número de onda alto (L3), no el sin flujo.
- Además, $\mathfrak B$ es afín en $Q'$ con coeficiente $r^4w'(r'w'-r''w)$ (comprobado).

**(L4) La garganta es un mínimo de $\beta$.** Desarrollando la EDO en $\rho=0$ ($R_k=r^{(k)}(0)$, $V_k=v^{(k)}(0)$):
$$\beta''(0)=-\frac{Q_0V_1^2}2-\frac{R_2\beta_0}{R_0}-\frac{\beta_0}{R_0^2}=-\frac{\beta_0}{R_0^2}\Big[1+R_0R_2+\frac x2\Big],\qquad x=\frac{Q_0R_0^2V_1^2}{\beta_0}<0 .$$
El corchete de $\hat\Lambda$ (Fase 2k §12) es $\beta_0^2[x^2+2R_0R_2x-4(1+R_0R_2)^2]$. El gradiente sano con $A_4<0$ ($\hat\Lambda<0$) exige
$$-x>R_0R_2+\sqrt{R_0^2R_2^2+4(1+R_0R_2)^2}>R_0R_2+2(1+R_0R_2),$$
y entonces $\beta''(0)>R_2\beta_0/(2R_0)>0$. **La garganta es un mínimo local estricto de $\beta=\mathcal NG_4$**: el flujo debe sobrecompensar el empuje de $K(0)\ge1$, que en el caso sin flujo hacía caer a $\beta$ desde la garganta (2k §5).

**(L5) Por qué la vía (i) no cierra con un peso fijo.** Sea $\mu(\rho)\ge0$ un peso que no dependa de $Q$. En $\int\mu\mathfrak B$ se puede eliminar $\beta$:
- escribiendo $2\mu(rK)'=\varphi'/r+\int K\varphi/r^3$, se tiene $\int\mu(-2\beta'(rK)')=\int\beta L\varphi=\int\varphi S$;
- (I1) es el caso $\varphi=1$.

Lo que queda es un funcional lineal en $Q$ y $Q'$. Al quitar $Q'$ por partes queda $\int Q\,C$, y $C$ contiene el término $-2\mu r^4r'w'w''$ (comprobado). Es un término cruzado en $w''$ sin término $w''^2$, así que $C$ no tiene signo definido. Como $A_4$ es libre (vía $F_4$, sujeto solo a $Q<0$), **ninguna desigualdad puntual que use solo $Q<0$ y $\mathfrak B<0$ cierra el caso**. Un teorema necesita más condiciones: hiperbolicidad, $s_0<0$, el régimen sin flujo o $\mathcal N>|v|$.

**La vía (ii)** (la condición de la garganta junto con $\beta>0$) da (L4), pero no una contradicción. Localmente, una garganta sana existe (2k §12: $V_3\lesssim-3V_1$ en Ellis), y $\beta>0$ es una condición global.

## 4. Numérica: la clase que dejan los lemas (`criba_mixta.py`, `optimiza_mixta.py`, `diagnostico_mixta.py`)

(L3) obliga a una cola de densidad decreciente. La familia de la Fase 2k con cola de energía positiva, $m=M-(M-\tfrac12)r^{-n}$, tiene además $\varepsilon>0$ en la garganta, y no consiguió $\beta>0$. La clase nueva combina un núcleo de Ellis con una cola de energía positiva:
$$m(r)=\frac1{2r}+M\big(1-r^{-n}\big)^2,\qquad m(1)=\tfrac12,\ m_r(1)=-\tfrac12\ (r''(0)=1),\qquad \varepsilon\simeq\frac{2Mn}{4\pi}r^{-n-3}\ \text{lejos}\ (n<1).$$
El flujo (dos gaussianas), $A_4$ con pico, el lapso con cola $1/\rho$ y los criterios son los de `fase2k/criba_masa.py`, reutilizado tal cual.

**Criba** (semilla 31; $M\in[0.05,2]$, $n\in[0.2,0.95]$, $a\in[0.3,300]$, $w_q\in[0.1,3]$, $\mathcal N(0)\in[0.1,0.9]$):

| fondos | $eta>0$ | lapso real ($\mathcal N>|v|$) | $(Qr')'>0$ lejos | $s_2>0$ lejos ($
ho>20$) | $rK$ crece lejos | hiperbólico en todo $
ho$ | $s_2>0$ en todo $
ho$ |
|---|---|---|---|---|---|---|---|
| 600 | 34 | 34 | 28 | 19 | 19 | **0** | **0** |

Entre los 34 con $eta>0$, la fracción de puntos con $s_2>0$ tiene mediana 0.76 y máximo 0.92; la de puntos con $s_0<0$ llega a 1. Las 19 colas sanas coinciden exactamente con las 19 en que $rK$ crece lejos, como predice (L3).

- **$\beta>0$ es frecuente** con cola de energía positiva. La Fase 2k §13 decía que no; era una propiedad de su familia ($\varepsilon>0$ en la garganta, $K(0)<3/2$, sin núcleo exótico).
- **Dónde falla** (`diagnostico_mixta.py`, los 9 primeros fondos con $\beta>0$):
  - $s_2<0$ en una capa junto a la garganta ($\rho\approx0.2$–$1.8$);
  - $s_2<0$ también en la zona donde el flujo decae, desde poco antes de su máximo hasta $\rho\approx8$–$60$;
  - capas no hiperbólicas dentro de esas zonas; $s_0>0$ cerca de la garganta en 8 de 9.
  - En la zona de decaimiento los tres términos de $\mathfrak B$ tienden a tener el signo malo: $Q'>0$ (el pico de $|A_4|$ se deshace), $w'<0$, y $(rK)'<0$ mientras $\varepsilon$ sube desde su valor negativo del núcleo.

**Evolución diferencial** (9 parámetros, población 108; objetivo $=$ fracción de puntos con $s_2\le0$ $+$ fracción con $s_0\ge0$ $+$ 0.5 si no es hiperbólico; $>2$ si $\beta$ tiene ceros):

| generación | 1 | 2 | 4 | 15 | 20 |
|---|---|---|---|---|---|
| mejor objetivo | 0.677 | 0.620 | 0.069 | 0.064 | 0.064 |

El óptimo ($M=1.89$, $n=0.71$, $\mathcal N(0)=0.41$, $w=4.86$, $a=156$, $w_q=1.76$, $w_2=0.30$, $\lambda=0.115$; $eta(0)=0.014$, $v_{\max}=0.13$):
- $eta>0$ en todo el eje, lapso real, $(Qr')'>0$ y $s_2>0$ en la cola, $rK$ creciente lejos;
- **hiperbólico en todos los puntos** (el primer fondo del programa con esta propiedad y cola de energía positiva);
- sano en la garganta ($
ho<0.56$) y fuera de $
hopprox1.9$;
- en $0.56\lesssim
ho\lesssim1.7$, $s_2\sim-10^4$ y $s_0\sim+10$: **capa fantasma**. Los conos están inclinados ($|s_1|/\sqrt{|s_2s_0|}$ entre 0.2 y 12), pero las dos velocidades características tienen signos opuestos (su producto es $s_0/s_2<0$; p. ej. $c=-0.041$ y $+0.019$). Por tanto el corte $	au$ es espacial para el monopolo y la energía negativa es un fantasma también con el criterio de Babichev et al.;
- el objetivo se estanca en 0.064 desde la generación 15 (20 generaciones, 2160 evaluaciones).

Es el modo de fallo que faltaba en la lista de la Fase 2k: además de la capa no hiperbólica y de la cola exótica, una **capa fantasma bien planteada** entre la garganta y el máximo del flujo, donde el pico de $|A_4|$ se deshace ($Q'>0$).

## 5. Estado

| | demostrado | evidencia | abierto |
|---|---|---|---|
| campo radial (2g, 2j, 2m R) | teorema $E[N]<0$ | 32 + 9 + 6 fondos | — |
| temporal sin flujo (2k I, 2m E) | no-go de $G_4$; fantasma en la garganta, regla común | reconstrucciones con nodos | — |
| cúbico con flujo (2k II) | (I1), (I2), símbolo cerrado, (L1)–(L5) | ≈1300 (2k) + 600 y ≈2200 (2n), ninguno sano; el mejor tiene una capa fantasma hiperbólica | teorema |
| Proca mixto ($A_0,A_1\neq0$) | — | — | todo |
| quíntico con $\phi=qt+\psi$ | maquinaria | — | reducción ADM |

Un teorema para el caso cúbico debería usar la hiperbolicidad (no solo el signo de $\mathfrak B$) y la ligadura entre $Q'$ y el valor de $Q$ en la garganta. La numérica sugiere dónde buscar: la zona en la que el flujo pasa de "sostener la garganta" a apagarse.

## 6. Parte 2: el artículo

`articulo/gargantas_monopolo.tex` (clase `article`, `babel` español, bibliografía `thebibliography`, sin dependencias externas):

1. introducción y convenciones;
2. teorías y estructura del monopolo;
3. campo radial: Teorema 1 ($d=N$), invariancia de marco;
4. campo temporal: Teorema 2 (cambio de signo de $G_4$, con demostración completa), Teorema 3 (regla común de la garganta, con la forma unificada de §2);
5. salidas que no funcionan: carga y rotación, fluido y sólido, BCK;
6. el caso cúbico con flujo: (I1), (I2), símbolo, Lemas 4–5 de §3, tabla numérica, Conjetura;
7. discusión (mecanismo, literatura, lo que no cubre, consecuencia para el híbrido) y conclusiones;
8. apéndice de reproducibilidad.

- **Bibliografía:** 33 entradas. Las de arXiv se han comprobado con su API (`arxiv_bib.py` → `resultados/arxiv_bib.json`): título, autores y referencia de revista. Las clásicas sin arXiv (Morris–Thorne, Ellis, Bronnikov) son citas estándar.
- **Compilación:** no hay `pdflatex`, `xelatex` ni `tectonic` en esta máquina, así que no se ha compilado. `comprueba_tex.py` revisa llaves, entornos, `\cite`/`\bibitem` y `\ref`/`\label`. Para compilar: `pdflatex gargantas_monopolo.tex` dos veces.
- **Autoría:** la línea de autor es "Iván"; hay que completarla antes de enviar.

## Archivos

| archivo | contenido |
|---|---|
| `lemas.py` | (L1)–(L5), simbólico; `resultados/lemas.json` |
| `regla_comun.py` | forma común de $T$ y $k(0)$ en la garganta (2k y 2m) |
| `criba_mixta.py` | criba de la clase núcleo de Ellis + cola de energía positiva (reutiliza `fase2k/criba_masa.py`) |
| `optimiza_mixta.py` | evolución diferencial sobre esa clase |
| `diagnostico_mixta.py` | intervalos con $s_2\le0$, $s_0\ge0$ y no hiperbólicos |
| `arxiv_bib.py` | verificación de la bibliografía con la API de arXiv |
| `articulo/gargantas_monopolo.tex` | el artículo |
| `articulo/comprueba_tex.py` | comprobación estructural del LaTeX |

## Reproducibilidad

```
cd fase2n
python lemas.py                 # segundos
python regla_comun.py
python criba_mixta.py 600       # ~1-2 h (necesita fase2k/resultados/simbolo_perfiles.pkl)
python optimiza_mixta.py 20 5   # ~2 h con 6 procesos
python diagnostico_mixta.py
python arxiv_bib.py             # red
cd articulo; python comprueba_tex.py; pdflatex gargantas_monopolo.tex; pdflatex gargantas_monopolo.tex
```
