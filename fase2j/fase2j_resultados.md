# Híbrido estrella–agujero de gusano — Fase 2j: términos quínticos y las gargantas de Bakopoulos–Charmousis–Kanti

*Documento de trabajo v0.1 · 2026-10-27 · continúa [fase2i/fase2i_resultados.md](../fase2i/fase2i_resultados.md) (§7, "salir de la GLPV cuártica")*

## Resumen

1. **La pregunta.** Las Fases 2g–2i comparten un obstáculo: en la GLPV cuártica, la prueba $d=N$ (el radio areal reescalado con el lapso) tiene energía $E[N]=-2\int N^3\mathcal G\sigma^2<0$ en toda garganta estática con el sector impar sano, y la materia no lo arregla. ¿Sobrevive ese déficit con términos quínticos? El caso de prueba son las gargantas de vacío de Bakopoulos–Charmousis–Kanti (BCK, [arXiv:2111.09857](https://arxiv.org/abs/2111.09857), JCAP 05 (2022) 022), que usan todos los galileones extendidos.

2. **Maquinaria quíntica** (§2). Ampliamos la lagrangiana covariante de la Fase 2f con
   $$L_5=G_5G_{\mu\nu}\phi^{\mu\nu}-\tfrac16G_{5X}\big[(\Box\phi)^3-3\Box\phi\,\phi_{\mu\nu}^2+2\phi_{\mu\nu}^3\big]$$
   (y opcionalmente $F_5$). Comprobaciones exactas:
   - identidad de Noether del fondo;
   - límite cuártico ($G_5=0$): ecuaciones del fondo y acción cuadrática **idénticas** a las de la Fase 2f;
   - el agujero negro de Lü–Pang (la semilla de BCK) cumple las ecuaciones quínticas con un residuo de $3\times10^{-41}$ (40 cifras).

3. **El monopolo en Horndeski quíntico** (§3). Se mantiene la estructura de la Fase 2f: $L_2=-FW+L_{\rm rest}(c,d)$, $c$ auxiliar, ecuación reducida de segundo orden y $\alpha_0=-\Phi'\alpha_1$ **exacta**. Por tanto, la prueba $d=Nf$ con $c=0$ sigue cumpliendo la ligadura.

4. **El déficit sobrevive** (§3.2). La energía de la familia $c=0$, reducida a su forma canónica (única), es **exactamente**
   $$E[Nf]\;\cong\;\int-2N^3\,\mathcal G\,(f'-\sigma f)^2\,d\ell,\qquad \mathcal G=G_4-2XG_{4X}+X\big(\Phi'\phi'G_{5X}+G_{5\phi}\big),$$
   el coeficiente del sector impar de Horndeski (Kobayashi–Motohashi–Suyama, [arXiv:1202.4893](https://arxiv.org/abs/1202.4893)), salvo derivadas totales. Con $G_5=0$ es la identidad de la Fase 2g; la parte con $G_5$ sin derivar se cancela, como exige la simetría de desplazamiento. **El teorema del monopolo se extiende a Horndeski quíntico:** toda garganta estática con $\mathcal G>0$ tiene $E[N]<0$.

5. **El cambio de marco** (§4). Con $\phi$ estático, una transformación disforme $\bar g=g+D(X)\partial\phi\partial\phi$ (o conforme $C(X)$) deja invariante la prueba $d=N$, $c=0$, $\delta\phi=0$: el lapso no cambia y $\delta X=0$. La energía de la prueba es la misma en los dos marcos. Por eso el resultado de la semilla (Horndeski quíntico) decide el de BCK (más allá de Horndeski).

6. **Las gargantas BCK** (§5). En el marco de la semilla (el exterior $r>r_0$ del agujero negro de Lü–Pang, dos copias), con la densidad sin integrar por partes (en $r_0$ el cambio de marco es singular):
   - $E[N]<0$ en los 9 juegos de parámetros probados ($-0.17$ a $-0.53$ en unidades de $M$), con $\mathcal G\ge1$;
   - **el coeficiente cinético del monopolo se anula idénticamente**, $K\equiv0$, sobre el agujero negro de Lü–Pang (y en fondos quínticos genéricos no: es una propiedad de esta solución);
   - $P<0$ en alguna zona en todos los casos.

   El monopolo de las gargantas BCK no tiene término cinético: está fuertemente acoplado, y el análisis lineal de estabilidad no se puede aplicar. La configuración $d=N$ tiene energía negativa y energía cinética nula.

---

## 1. Planteamiento

**Las gargantas de BCK** (de su artículo, verificado en arXiv):
- **Semilla:** la teoría de Horndeski de Lü–Pang (la "4D Einstein–Gauss–Bonnet" regularizada, de origen Lovelock), con $X=-\tfrac12(\partial\phi)^2$,
  $$G_2=8\alpha X^2,\quad G_3=-8\alpha X,\quad G_4=1+4\alpha X,\quad G_5=-4\alpha\ln|X| ,$$
  y su agujero negro $\bar h=1+\frac{r^2}{2\alpha}\big(1-\sqrt{1+8\alpha M/r^3}\big)$, $\bar\phi_r=(\sqrt{\bar h}-1)/(r\sqrt{\bar h})$. En longitud propia, $\phi'=(\sqrt{\bar h}-1)/r$.
- **Transformación disforme:** $g_{\mu\nu}=\bar g_{\mu\nu}-D(\bar X)\partial_\mu\phi\partial_\nu\phi$, con $D=-b_1/(2b_1\bar X+\sqrt{-2\bar X})$ y $b_1=r_0/\lambda$.
- **El agujero de gusano:** $ds^2=-h\,dt^2+dr^2/(hW^{-1})+r^2d\Omega^2$, con $W^{-1}=1-(r_0/\lambda)\sqrt{-2\bar X}$. La garganta está donde $W^{-1}=0$, es decir, $h(r_0)=(1-\lambda)^2$, y su radio es $r_0=\big[M+\sqrt{M^2-\alpha\lambda^3(2-\lambda)^3}\big]/[\lambda(2-\lambda)]$.
- **Teoría transformada:** GLPV quíntica, con $G_2,\dots,F_5$ y la condición de degeneración $XG_{5X}F_4=3F_5(G_4-2XG_{4X})$. Sus funciones son singulares en la garganta: $G_4\propto(1-\sqrt2b_1y)^{1/2}\to0$ y $F_5\propto(1-\sqrt2b_1y)^{-7/2}\to\infty$, con $y=\sqrt{-\bar X}$.
- El artículo no estudia perturbaciones; no hemos encontrado ningún estudio posterior de su monopolo.

**Estrategia.**
1. Derivar el monopolo de Horndeski quíntico (la semilla) con la maquinaria de la Fase 2f ampliada.
2. Comprobar si la identidad de la prueba $d=N$ sobrevive.
3. Llevar el resultado a BCK por el cambio de marco (§4), que deja invariante esa prueba.

## 2. La maquinaria quíntica (`q_simbolico.py`, `generar_q.py`)

Lagrangiana covariante sobre la métrica esférica general $(A,b,B,C)$ con $\phi(t,\ell)$, como en `fase2f/bH_simbolico.py`, con el tensor de Ricci y
$$L_5=G_5G_{\mu\nu}\phi^{\mu\nu}-\tfrac16G_{5X}\big[(\Box\phi)^3-3\Box\phi\,\phi_{\mu\nu}\phi^{\mu\nu}+2\phi_\mu{}^\nu\phi_\nu{}^\rho\phi_\rho{}^\mu\big]+F_5\,\epsilon\epsilon\,\phi\phi\,\phi_{..}\phi_{..}\phi_{..}$$
(el término $F_5$ queda como opción; para la semilla, Horndeski, $F_4=F_5=0$).

| Comprobación | Resultado |
|---|---|
| Independencia de $\theta$ de la lagrangiana reducida | $10^{-165}$ |
| Identidad de Noether del fondo | idénticamente nula |
| Límite cuártico ($G_5=F_5=0$): ecuaciones del fondo | idénticas a las de la Fase 2f |
| Límite cuártico: acción cuadrática | idéntica a la de la Fase 2f |
| Agujero negro de Lü–Pang en la teoría de Lü–Pang (3 puntos, 40 cifras) | residuo $3\times10^{-41}$ |

La última comprueba de una vez el término $G_5$ (signos, convenciones y la parte con $G_{5X}$): $G_5=-4\alpha\ln|X|$ equivale al acoplo $\alpha\phi\,\mathcal G_{\rm GB}$.

## 3. El monopolo de Horndeski quíntico

### 3.1 Estructura (`reduccion_q.py`)

En el gauge unitario ($\delta\phi=0$), con aritmética exacta en chorros de la capa de masas:
- **La ligadura:** $E_\beta=i\omega W$, con $W=\gamma c+\alpha_1d'+\alpha_0d$, y $\alpha_0+\Phi'\alpha_1=0$ **idénticamente**.
  $$\alpha_1=-\frac{4r}N\,\tilde{\mathcal G},\qquad \tilde{\mathcal G}=G_4+G_{4X}\phi'^2-\frac{r'}{2r}\phi'^3G_{5X}-\tfrac12\phi'^2G_{5\phi},$$
  que es el $\mathcal H/2$ del sector impar de Kobayashi–Motohashi–Suyama.
- **Estructura** $L_2=-FW+L_{\rm rest}(c,d)$, con $F=N^2a'+\dot\beta$, y **$c$ auxiliar**.
- **Problema reducido:** $\int[K\dot d^2-Pd'^2-Ud^2]$ de segundo orden (33, 58 y 521 términos).

La prueba $d=Nf$, $c=0$ cumple la ligadura igual que en la teoría cuártica.

### 3.2 La identidad de la prueba $d=N$ (`identidad_q.py`)

Puntualmente, $-L_{\rm rest}(0,Nf)$ ya no es $-2N^3\mathcal G(f'-\sigma f)^2$: $G_5$ añade términos con $f''$, y en parte con $G_5$ sin derivar. Reducida a su **forma canónica** $Af'^2+Cf^2$ (única, integrando por partes):
$$A=-2N^3\mathcal G,\qquad C=-2N^3\mathcal G\,\sigma^2-\big(2N^3\mathcal G\sigma\big)',\qquad \mathcal G=G_4-2XG_{4X}+X\big(\Phi'\phi'G_{5X}+G_{5\phi}\big),$$
**exactamente** (identidad simbólica fuera de la capa de masas).
- $\mathcal G$ es el $\mathcal G/2$ de Kobayashi–Motohashi–Suyama: el coeficiente del sector impar con $\Phi'$ (y no con $r'/r$, que es el de $\alpha_1$).
- $\partial A/\partial G_5=\partial C/\partial G_5=0$: la parte con $G_5$ sin derivar es una derivada total, como exige la simetría $G_5\to G_5+$ cte.
- Con $G_5=0$ se recupera exactamente la identidad de la Fase 2g.

**Teorema (Horndeski quíntico).** En toda garganta estática, simétrica y asintóticamente plana de una teoría de Horndeski cualquiera ($G_2,\dots,G_5$) con $\mathcal G>0$ y el monopolo bien planteado ($K,P>0$):
$$E[N]=-2\int N^3\mathcal G\,\sigma^2\,d\ell<0 ,$$
y hay un taquión en el monopolo. Las derivadas totales no contribuyen, porque la garganta es lisa y la prueba decae lejos con el corte de la Fase 2g. $\mathcal G>0$ es una de las condiciones de Kobayashi–Motohashi–Suyama para el sector impar, así que el teorema es: **sector impar sano ⇒ taquión en el monopolo**, como en la teoría cuártica. (Horndeski ya tenía un teorema de no existencia de gargantas sanas; este da el taquión explícito.)

## 4. El cambio de marco

**Proposición.** Sea $\bar g_{\mu\nu}=C(X)\,g_{\mu\nu}+D(X)\,\partial_\mu\phi\partial_\nu\phi$ con $\phi$ estático. La prueba $d=Nf$, $c=0$, $\delta\phi=0$ (con $\beta=0$) de un marco es la prueba $d=\bar Nf$, $\bar c=0$, $\delta\phi=0$ del otro, y su energía es la misma.
- *Por qué.* Con $\phi=\phi(\ell)$, $X=-\tfrac12g^{\ell\ell}\phi'^2$, y $\delta g^{\ell\ell}=-2c$ (los términos en $a$ se cancelan), así que $\delta X=0$ con $c=0$, también a segundo orden si $\beta=0$.
- *El mapa de la prueba.* $\delta\bar g_{\ell\ell}=C\,\delta g_{\ell\ell}=0$; $\bar g_{tt}=Cg_{tt}$, de modo que $\bar N=\sqrt C\,N$; $\delta\bar r=\sqrt C\,\delta r=\bar Nf$.
- *La energía.* La acción es la misma escrita en otras variables. Al segundo orden solo añade términos proporcionales a las ecuaciones del fondo, que se anulan, y a la parte de segundo orden del mapa, que es nula para esta prueba.

**Consecuencias.**
- El teorema de la Fase 2g (GLPV cuártica) vale para sus imágenes por $C(X)$, $D(X)$ regulares: la clase cuártica de DHOST que se obtiene así.
- El de §3 vale para las imágenes de Horndeski quíntico, siempre que el mapa sea regular en la garganta.
- **BCK no es de ese tipo:** el mapa es singular justo en la garganta ($D\to\infty$ donde $W^{-1}=0$). Allí hay que calcular la energía en el marco de la semilla sin integrar por partes (§5).

## 5. Las gargantas BCK en el marco de la semilla (`bck_semilla.py`)

**La geometría en la semilla.** El agujero de gusano son dos copias del exterior $r>r_0$ del agujero negro de Lü–Pang, pegadas en $r_0$. En la semilla, $r_0$ es un punto regular ($dr/d\bar\ell=\sqrt{\bar h}\neq0$). Las funciones pares de la longitud propia $\ell$ del agujero de gusano son funciones regulares de $r-r_0\propto\ell^2$.

**Cómo se calcula la energía.** Una derivada total en $\bar\ell$ no lo es en $\ell$: en $r_0$ las dos copias contribuyen con el mismo signo. Por eso la energía se calcula con la densidad tal cual,
$$E[N]=2\int_{r_0}^\infty m_{00}\,d\bar\ell,\qquad m_{00}=-L_{\rm rest}(0,N)\big|_{f=1}=C+b'/2 ,$$
y se separa en volumen, $-4\int N^3\mathcal G\bar\sigma^2$, y borde, $2\big[2N^3\mathcal G\bar\sigma-b/2\big]_{r_0}$.

| $\alpha$ | $\lambda$ | $r_0$ (= BCK) | $\bar\sigma(r_0)$ | $E[N]$ | volumen | borde | $\mathcal G$ | $\min P$ |
|---|---|---|---|---|---|---|---|---|
| 0.01 | 0.5 | 2.664 | $-0.093$ | $-0.4455$ | $-0.4454$ | $-6.6\times10^{-5}$ | 1.000–1.003 | $-0.007$ |
| 0.1 | 0.5 | 2.638 | $-0.086$ | $-0.4472$ | $-0.4466$ | $-5.6\times10^{-4}$ | 1.000–1.028 | $-0.070$ |
| 0.5 | 0.5 | 2.518 | $-0.052$ | $-0.4537$ | $-0.4526$ | $-1.1\times10^{-3}$ | 1.000–1.139 | $-0.357$ |
| 0.9 | 0.5 | 2.383 | $-0.014$ | $-0.4604$ | $-0.4602$ | $-1.5\times10^{-4}$ | 1.000–1.248 | $-0.655$ |
| 0.5 | 0.2 | 5.523 | $+0.105$ | $-0.4092$ | $-0.4033$ | $-6.0\times10^{-3}$ | 1.000–1.004 | $-0.003$ |
| 0.5 | 0.8 | 1.820 | $-0.764$ | $-0.5067$ | $-0.4875$ | $-0.019$ | 1.000–1.962 | $-9.2$ |
| 0.5 | 0.95 | 1.714 | $-3.58$ | $-0.5257$ | $-0.5191$ | $-0.0066$ | 1.000–5.31 | $-175$ |
| 0.2 | 0.05 | 20.51 | $+0.044$ | $-0.1680$ | $-0.1679$ | $-8.2\times10^{-5}$ | 1.000 | $-10^{-6}$ |
| 0.99 | 0.3 | 3.788 | $+0.098$ | $-0.4617$ | $-0.4498$ | $-0.012$ | 1.000–1.040 | $-0.045$ |

($M=1$. La suma volumen + borde reproduce $E[N]$ a $10^{-12}$.)

- **El déficit sobrevive en todas:** $E[N]<0$, también cuando la garganta queda fuera de la esfera de fotones de la semilla ($\bar\sigma(r_0)>0$), donde el término de borde podría haber cambiado el signo.
- **$K\equiv0$.** El coeficiente cinético del problema reducido se anula en todo $r>r_0$, en todos los casos ($|K|\lesssim10^{-12}$ con 15 cifras; ver la comprobación con 40 cifras abajo). En fondos quínticos genéricos (chorros exactos de la capa de masas) $K$ vale 0.10, $-2.05$ y 0.22: **es una propiedad del agujero negro de Lü–Pang**, no de la reducción.
- **$P<0$** en alguna zona en todos los casos.

**$K\equiv0$ es una propiedad del agujero negro, no de la teoría** (`k_lu_pang.py`, `degeneracion_lp.py`):

| fondo | $K$ | $P$ |
|---|---|---|
| agujero negro de Lü–Pang, 4 puntos ($\alpha=0.1$–0.9, $M=1,2$), 40 cifras | $\le5\times10^{-41}$ | $-0.058$, $-0.015$, $-7.8\times10^{-4}$, $-2.3\times10^{-3}$ |
| chorros genéricos de la capa de masas de la teoría quíntica (4) | $0.10$, $-2.05$, $0.22$, $0.49$ | |
| soluciones genéricas de la teoría de Lü–Pang (3 datos iniciales al azar, $E_B=0$) | $+0.004$, $+0.18$, $-0.44$ | $-0.017$, $-1.4$, $-4.1$ |

- El sistema del fondo de la teoría de Lü–Pang no es degenerado: la matriz de $(\Phi'',r'',\phi'')$ tiene rango 3.
- La anulación es propia de la solución. En ella el escalar está atado a la geometría (en longitud propia, $r'=1+r\phi'$).

**Lectura.** El monopolo de las gargantas BCK no tiene término cinético. Por la invariancia de §4, $K_{\rm gusano}\,d\ell=K\,d\bar\ell=0$. La perturbación esférica no es dinámica a orden lineal: está **fuertemente acoplada**, y su evolución la deciden los términos no lineales. Además:
- la configuración $d=N$ tiene energía negativa y energía cinética nula;
- el término de gradiente es negativo en parte del exterior.

Las gargantas de BCK no son, por tanto, agujeros de gusano linealmente estables. El déficit del vacío sobrevive, y el precio de la construcción disforme singular es un monopolo sin dinámica lineal.

## 6. Conclusiones y siguiente paso

**Lo que queda establecido.**
1. **La maquinaria quíntica** (Horndeski con $G_5$; el término $F_5$ como opción), validada con el límite cuártico exacto y con el agujero negro de Lü–Pang a 40 cifras.
2. **El teorema del monopolo se extiende a Horndeski quíntico.** La energía de la familia $c=0$ es, salvo derivadas totales, exactamente $-2N^3\mathcal G(f'-\sigma f)^2$, con el $\mathcal G$ del sector impar. Toda garganta estática con el sector impar sano tiene $E[N]<0$.
3. **El cambio de marco** con $\phi$ estático ($C(X)$, $D(X)$) deja invariante la prueba $d=N$. El teorema vale también para las imágenes regulares de la GLPV cuártica y de Horndeski quíntico.
4. **Las gargantas de Bakopoulos–Charmousis–Kanti:**
   - $E[N]<0$ en todos los parámetros probados;
   - su monopolo tiene $K\equiv0$ (acoplo fuerte) y $P<0$ en parte del exterior.
   No son una vía hacia una garganta estable.

**El mecanismo común (Fases 2g–2j).** La prueba $d=N$ (reescalar el radio areal con el lapso) tiene energía $-2\int N^3\mathcal G\sigma^2$ en toda la familia de teorías escalares-tensoriales degeneradas que hemos podido tratar. Es una propiedad de la geometría de una garganta estática con un escalar estático: $\sigma=(\ln r/N)'$ no puede anularse en un fondo asintóticamente plano. La materia no la compensa (Fases 2g–2i), y los términos quínticos tampoco.

**Lo que no queda establecido.**
- El término $F_5$ (GLPV quíntica con la condición de degeneración) no se ha reducido de forma general. Su caso relevante, BCK, se ha tratado por el cambio de marco.
- El análisis lineal de BCK no puede ir más allá de $K\equiv0$. Haría falta teoría de perturbaciones no lineal.

**Siguiente paso: lo que el teorema no cubre.** La prueba necesita tres cosas: un escalar estático (gauge unitario sin tiempo), un fondo estático y que $c$ sea auxiliar. Las salidas son:
1. **Un escalar con dependencia lineal en el tiempo**, $\phi=qt+\psi(r)$, con métrica estática (posible en teorías con simetría de desplazamiento, como las soluciones "stealth" de Babichev–Charmousis). El gauge unitario fija entonces el tiempo, no la coordenada radial, y la estructura del monopolo es la del fluido de la Fase 2g, no la del vacío. La familia $c=0$ no existe en la misma forma. Riesgo conocido: el acoplo fuerte de las perturbaciones de soluciones "stealth".
2. **Rotación** (gargantas estacionarias), donde el sector $l=0$ se mezcla con el impar.
3. **Campos vectoriales** (Proca generalizada), con más grados de libertad en $l=0$.

Recomendamos el primero: mantiene la maquinaria de las Fases 2f–2j y ataca directamente la hipótesis de la que depende la prueba $d=N$.

## Archivos

| archivo | contenido |
|---|---|
| `q_simbolico.py` | lagrangiana quíntica covariante en simetría esférica ($G_2,\dots,G_5$, $F_4$, $F_5$ opcional) |
| `generar_q.py` | ecuaciones del fondo y acción cuadrática; Noether, límite cuártico, agujero negro de Lü–Pang |
| `chorros_q.py` | chorros exactos de la capa de masas con $G_5$, $F_5$ |
| `reduccion_q.py` | estructura del monopolo en el gauge unitario, $c$ auxiliar, $\alpha_0=-\Phi'\alpha_1$, $K$, $P$, $U$ |
| `identidad_q.py` | forma canónica de la energía de la familia $c=0$: $A=-2N^3\mathcal G$, $C=-2N^3\mathcal G\sigma^2-(2N^3\mathcal G\sigma)'$ |
| `bck_semilla.py` | gargantas BCK en el marco de la semilla: $E[N]$ (volumen y borde), $\mathcal G$, $K$, $P$ |
| `k_lu_pang.py`, `degeneracion_lp.py` | $K$ con 40 cifras en el agujero negro de Lü–Pang y en soluciones genéricas |
| `espectro_semilla.py` | ángulo de Prüfer en $r>r_0$ (preparado; no se usa porque $K\equiv0$) |

## Reproducibilidad

```
cd fase2j
python generar_q.py               # ~26 min: teoria quintica, comprobaciones -> resultados/q_simbolico.pkl
python reduccion_q.py             # ~1 min: monopolo en el gauge unitario -> reduccion_q.pkl
python identidad_q.py             # forma canonica de la familia c = 0
python bck_semilla.py             # ~10 min: gargantas BCK en el marco de la semilla
python k_lu_pang.py               # K en el agujero negro de Lu-Pang (40 cifras) y en chorros genericos
python degeneracion_lp.py         # soluciones genericas de la teoria de Lu-Pang
```
