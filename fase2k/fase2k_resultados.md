# Híbrido estrella–agujero de gusano — Fase 2k: gargantas con un escalar lineal en el tiempo, $\phi=qt+\psi(\ell)$

*Documento de trabajo v0.1 · 2026-10-27 · continúa [fase2j/fase2j_resultados.md](../fase2j/fase2j_resultados.md) (§6, "un escalar con dependencia lineal en el tiempo")*

## Resumen

1. **La pregunta.** El teorema del monopolo de las Fases 2g–2j ($E[N]=-2\int N^3\mathcal G\sigma^2<0$) usa un escalar estático: el gauge unitario $\delta\phi=0$ fija la coordenada radial y el lapso entra como multiplicador. En teorías con simetría de desplazamiento la métrica puede ser estática con $\phi=qt+\psi(\ell)$ (Babichev–Charmousis). Entonces $\delta\phi=0$ fija el **tiempo**. ¿Escapan estas gargantas al teorema?

2. **Maquinaria** (§2, `generar_t.py`). Ampliamos la de las Fases 2f/2j a un fondo con $\phi_t=q$, para la GLPV cuártica con $G_i(X)$ y, como opción, $G_5(X)$. Comprobaciones exactas:
   - identidad de Noether radial;
   - $E_\phi=-J'$ y $e^{2\Phi}E_{t\ell}=q\,J$ (fuera de la capa de masas): $J=0$ resuelve las dos;
   - límite $q=0$ idéntico a las Fases 2f y 2j (fondo y acción cuadrática);
   - el agujero negro "stealth" de Babichev–Charmousis cumple las ecuaciones a $10^{-51}$;
   - los dos modos gauge son exactos en chorros de la capa de masas (cuártico y quíntico).

3. **El monopolo en variables ADM del corte $\phi=$ cte** (§3). Con $\tau=t+\psi/q$, el fondo es estacionario con desplazamiento $v$. En las variables ADM de ese corte, **el lapso y el desplazamiento son auxiliares algebraicos** (comprobado en chorros exactos), como en la cosmología en gauge unitario. Es la estructura del fluido de la Fase 2g, no la del vacío: la ligadura $W=0$ que salvaba la prueba $d=N$ ya no existe.

4. **El sector estático en forma cerrada** (§4, exacta). En la longitud propia $\rho$ del corte, con el lapso $\mathcal N$ del corte y $Q=2A_4/\mathcal N$:
   $$\mathcal L_{\rm est}=\mathcal N r^2A_2-A_3(r^2v)'+\frac{2A_4}{\mathcal N}\,r'(v^2r)'+\mathcal NB_4\big(2-4rr''-2r'^2\big),$$
   $$A_2=G_2,\quad A_3=yG_3-\!\int_0^y\!G_3\,dy\ \ (y=\sqrt{2X}),\quad A_4=-G_4+2XG_{4X}+4X^2F_4,\quad B_4=G_4 .$$
   El operador de Euler–Lagrange de la diferencia con la lagrangiana covariante se anula idénticamente. La ligadura de momento es $r\,A_3'=2v\,(Qr')'$.
   - **Lema.** Sin término cúbico ($A_3'\equiv0$ en el fondo, p. ej. $G_3=0$), una garganta ($r'=0$) obliga a $\psi'\equiv0$ o a $A_4\equiv0$.

5. **No-go para $\phi=qt$** (§5). Si $v=0$, la función $\beta=\mathcal NG_4$ cumple, sin ninguna otra función libre,
   $$\Big(\frac{\beta'}r\Big)'+\frac{K\beta}{r^3}=0,\qquad K=1+rr''-r'^2,\qquad \int\frac K{r^3}\,d\rho=\int\frac{1+r'^2}{r^3}\,d\rho>0 .$$
   Dividiendo por $\beta$ e integrando sale una contradicción si $\beta$ no cambia de signo.
   - **En toda garganta estática (dos extremos asintóticamente planos) con $\phi=qt$, de cualquier GLPV cuártica con simetría de desplazamiento, $G_4$ cambia de signo.**
   - Con el lema, el resultado vale para **cualquier $\psi$ si $G_3=0$**.
   - $G_4=B_4$ es el coeficiente de $R^{(3)}$, el término de gradiente de las ondas gravitatorias en el corte. Donde se anula o cambia de signo hay acoplo fuerte o inestabilidad de gradiente tensorial.
   - Comprobaciones: la ecuación de $\beta$ se cumple exactamente en chorros covariantes. Las reconstrucciones numéricas (3 familias de geometrías) encuentran modos cero, todos con nodos.

6. **El monopolo con $\phi=qt$** (§6, exacto). La acción cuadrática dependiente del tiempo es $L_2=L_{\rm est,2}(\nu,d)-r(Qr')'w^2+2r(Q\dot d)'w+Q\dot d^2$. Eliminando el desplazamiento,
   $$T[\dot d]=\int\Big[Q\dot d^2+\frac{r\,[(Q\dot d)']^2}{(Qr')'}\Big]d\rho ,$$
   sin términos giroscópicos. En la garganta, $Q$ y $(Qr')'=Qr''$ tienen el signo de $A_4(0)$:
   - $A_4(0)<0$ (ondas gravitatorias sanas): **fantasma en el monopolo**, localizado en la garganta;
   - $A_4(0)>0$: fantasma tensorial.

   Es la versión general de "fluido fantasma = fantasma" (con $G_4=1$, $(Qr')'\propto\rho+p$). La **familia $\delta X=0$** (la análoga de $d=N$) tiene energía $-2\int\beta r^2[(d/r)']^2$, la misma forma que la de las Fases 2g–2j, con la ligadura $(rd')'+(1-r'^2-rr'')d/r=0$ en lugar de $d\propto N$.

7. **Con término cúbico** (§7, $G_3\neq0$, $\psi'\neq0$). La ecuación de $\beta$ gana una fuente, $-[Qr^4(v^2/r^2)']'/4r^3$, y el no-go deja de ser automático. La Parte II (§10–§15) estudia este sector.

8. **Parte II: galileón cúbico con flujo** (§10–§15).
   - **(a) No hay teorema.** Las dos identidades integrales exactas de la ecuación de $\beta$ con fuente no fuerzan un cambio de signo:
     $$\int\!\frac{K\beta}{r^3}=\tfrac34\!\int\!\Big(\frac vr\Big)^2(Qr')',\qquad \int\!\frac1r\Big(\frac{\beta'}\beta+\frac{Qr^2u'}{8\beta}\Big)^2+\int\!\frac{1+r'^2}{r^3}=-\tfrac34\!\int\!\frac{Qu'r'}\beta+\frac1{64}\!\int\!\frac{Q^2r^3u'^2}{\beta^2}.$$
     Hay fondos explícitos con $G_4>0$, $A_4<0$ y flujo subluminal: 54 de 500 en la familia de Ellis.
   - **(b) Análisis local del monopolo**, exacto: la lagrangiana ADM dependiente del tiempo con desplazamiento coincide con la covariante, y el símbolo principal se calcula en cualquier punto. Con flujo, el coeficiente cinético es puramente geométrico:
     $$P_2=\frac{8Q^2r}{\mathcal N'^2}\,\mathfrak B,\qquad \mathfrak B=-2\beta'(rK)'-\frac{Qvr''(r^2W)'}r+Q'W(r'W-rvr''),\quad W=r^2(v/r)' .$$
     En la garganta, el símbolo solo depende de datos de Taylor de la geometría y del flujo. Una garganta sana es posible: hace falta un flujo fuerte que se curve cerca de la garganta.
   - **Criba numérica** (≈1300 fondos en tres familias): ninguno pasa a la vez las condiciones necesarias. Fallan por una de dos razones:
     - una capa no hiperbólica junto a la garganta, donde $\mathfrak B$ cambia de signo (velocidades imaginarias, problema mal planteado);
     - una cola fantasma lejos. La cola de Ellis es exótica: $(Qr')'<0$ y el monopolo tiene energía cinética negativa.

     Con colas de energía positiva, $\beta$ cambia de signo cerca de la garganta en los 400 fondos de la criba. La optimización se acerca a cero ($\min\beta/\max|\beta|=-0.029$) en el borde de la caja, así que no es concluyente.
   - Sin un fondo que pase las condiciones locales, el espectro no tiene sentido: con velocidades imaginarias el crecimiento es $\propto k$, sin cota.

**Veredicto.** El escalar lineal en el tiempo no da gargantas sanas sin término cúbico. Con $G_3=0$ (o $A_3'\equiv0$ en el fondo):
- $G_4$ cambia de signo;
- con $\psi'=0$, el monopolo en la garganta es un fantasma o tiene velocidades imaginarias, o el sector tensorial es un fantasma.

Con término cúbico y flujo no hay obstrucción algebraica, pero ninguno de los ≈1300 fondos construidos tiene un monopolo localmente bien planteado y sin fantasma.

---

## 1. Planteamiento

**Fondo.**
$$ds^2=-e^{2\Phi(\ell)}dt^2+d\ell^2+r(\ell)^2d\Omega^2,\qquad \phi=qt+\psi(\ell),\qquad X=\frac{q^2}{2N^2}-\frac{\psi'^2}2,\quad N=e^\Phi .$$
La teoría es la GLPV cuártica con simetría de desplazamiento, $G_2,G_3,G_4,F_4$ funciones de $X$ (convenciones de la Fase 2f), y como opción $G_5(X)$. La lagrangiana no depende de $\phi$, así que las ecuaciones del fondo son independientes de $t$. La ecuación $(t\ell)$ es proporcional a la corriente de Noether $J^\ell$, y $J=0$ es necesario.

**Literatura** (verificada con la API de arXiv):
- Babichev–Charmousis, *Dressing a black hole with a time-dependent Galileon*, [arXiv:1312.3204](https://arxiv.org/abs/1312.3204), JHEP 08 (2014) 106: la condición $J^r=0$, el agujero negro "stealth" y la familia de Schwarzschild–de Sitter autoajustada.
- Ogawa–Kobayashi–Suyama, [arXiv:1510.07400](https://arxiv.org/abs/1510.07400), PRD 93 (2016) 064078: inestabilidad de agujeros negros con pelo en Horndeski con simetría de desplazamiento.
- Babichev–Charmousis–Esposito-Farèse–Lehébel, [arXiv:1712.04398](https://arxiv.org/abs/1712.04398), PRL 120 (2018) 241101, y [arXiv:1803.11444](https://arxiv.org/abs/1803.11444), PRD 98 (2018) 104050: con conos causales inclinados, una energía no acotada en un corte no implica inestabilidad, y el criterio depende del corte.
- de Rham–Zhang, [arXiv:1907.00699](https://arxiv.org/abs/1907.00699), PRD 100 (2019) 124023: acoplo fuerte de las perturbaciones de agujeros negros "stealth" en DHOST.
- Kobayashi–Motohashi–Suyama II, [arXiv:1402.6740](https://arxiv.org/abs/1402.6740), PRD 89 (2014) 084042: sector par con escalar estático.

No hemos encontrado (búsqueda en arXiv: *wormhole* + *time-dependent* + *Horndeski*) estudios de gargantas con $\phi=qt+\psi$.

**Lo que cambia respecto de las Fases 2g–2j.**
- El gauge unitario $\delta\phi=0$ fija ahora el corte $t+\psi/q=$ cte, y queda libre un difeomorfismo radial: $\xi^\ell$ con $\xi^t=-(\psi'/q)\xi^\ell$.
- El lapso deja de ser un multiplicador: entra con término cuadrático (la "densidad + presión" del escalar).
- La prueba $d=N$ ya no cumple las ligaduras automáticamente.

## 2. Maquinaria y comprobaciones (`generar_t.py`, `chorros_t.py`, `verificar_t.py`)

La lagrangiana covariante es la de `fase2j/q_simbolico.py`. Cambios respecto de `fase2f/generar_bH.py`:
- funciones de la teoría solo de $X$ (símbolos $G_{0j}$);
- $p_{10}=q$ en el fondo;
- $D_\ell X_b=-q^2\Phi'/N^2-\psi'\psi''$;
- la lagrangiana no depende de $p_{00}$ (se comprueba).

Las derivadas de la lagrangiana en el fondo no dependen de $t$, así que el fondo sigue siendo $E_u=\sum_J(-D_\ell)^{J_\ell}\partial L/\partial u_J$ y la acción cuadrática tiene coeficientes independientes de $t$. Tamaños: cuártica, 158 entradas no nulas del hessiano y 1608 términos (17 min); quíntica, 44 min.

| Comprobación | Cuártica | Quíntica ($G_5$) |
|---|---|---|
| Identidad de Noether radial $E_AA'+E_CC'+E_pp'-2D_\ell E_B=0$ | idéntica | idéntica |
| $E_p=-D_\ell J$, $J=\partial L/\partial p_{01}-D_\ell\,\partial L/\partial p_{02}$ | idéntica | idéntica |
| $qE_p+D_\ell(e^{2\Phi}E_b)=0$ y $e^{2\Phi}E_b=q\,J$ | idénticas | idénticas |
| Límite $q=0$: fondo y acción cuadrática frente a 2f / 2j | idénticos | idénticos |
| Stealth Schwarzschild ($G_4=\zeta+\beta X$, $h=1-\mu/r$, $\psi_r=\pm q\sqrt{\mu/r}/h$, $X=q^2/2$), 2 signos × 2 puntos, 50 cifras | $\le5\times10^{-51}$ | $\le5\times10^{-51}$ |
| Control: el mismo fondo con $\psi'\times1.1$ | residuo 0.19 | 0.19 |
| Modos gauge $\xi^t$ ($a=i\omega\zeta$, $\beta=-N^2\zeta'$, $\delta\phi=q\zeta$) y $\xi^\ell$ ($a=\Phi'\xi$, $\beta=i\omega\xi$, $c=\xi'$, $d=r'\xi$, $\delta\phi=\psi'\xi$): 5 ecuaciones lineales, 2 chorros exactos | 0 exacto | 0 exacto |

- $G^{\mu\nu}\phi_\mu\phi_\nu\cong XR+(\Box\phi)^2-\phi_{\mu\nu}^2$ (salvo derivada total): la teoría de Babichev–Charmousis $\zeta R+\beta G^{\mu\nu}\phi_\mu\phi_\nu$ es $G_4=\zeta+\beta X$.
- *Chorros de la capa de masas.* El sistema $(E_A,E_B,J)$ es degenerado en sus derivadas más altas: una combinación es una ligadura de orden menor. Se toma $G_2$ como incógnita en cada orden (reconstrucción local, como $G_{2,00}$ en las fases anteriores): nivel $k$, $(r_k,\psi_k,G_{2,0(k-1)})$ de $(D^{k-2}J,D^{k-2}E_A,D^{k-2}E_B)$, con $\Phi_k$ libre.

## 3. Gauge unitario temporal: variables ADM del corte $\phi=$ cte (`adm_t.py`)

Con $\tau=t+\psi/q$ el escalar es $\phi=q\tau$ y la métrica es estacionaria con desplazamiento:
$$g_{\tau\tau}=-N^2,\quad N_\ell=g_{\tau\ell}=N^2s,\quad h_{\ell\ell}=1-N^2s^2,\quad \mathcal N^2=\frac{N^2}{1-N^2s^2},\qquad s=\psi'/q .$$
La métrica inducida en los cortes es riemanniana si $X>0$ (gradiente temporal). Las perturbaciones del corte, en las variables estáticas ($\delta A=2N^2a$, $\delta g_{t\ell}=\beta$, $\delta g_{\ell\ell}=2c$, $\delta r=d$), son:
$$\tilde\beta=\delta N_\ell=\beta+2N^2sa,\qquad h=\tfrac12\delta h_{\ell\ell}=c-s\beta-s^2N^2a,\qquad n=\tfrac12\delta(\mathcal N^2)=N^2a+\frac{N^2s\tilde\beta}{1-N^2s^2}-\frac{N^4s^2h}{(1-N^2s^2)^2}.$$
Ecuaciones de Euler–Lagrange en las variables nuevas ($E_v=M^TE_u$, cambio puntual), en 2 chorros exactos de la capa de masas con $\omega\neq0$:

| ecuación | amplitudes presentes |
|---|---|
| $E_n$ | solo $n_0$ (ni $n'$ ni $i\omega n$); $\tilde\beta_0,\tilde\beta_1,i\omega\tilde\beta_0$; $h_0,h_1,i\omega h_0$; $d$ hasta $d''$, $\omega^2d$ |
| $E_{\tilde\beta}$ | solo $\tilde\beta_0$; $n_0,n_1,i\omega n_0$; $h_0,h_1,i\omega h_0$; $d$ |
| $E_h$ | solo $h_0$; $n,\tilde\beta$ con una derivada; $d$ |
| $E_d$ | todo, hasta segundo orden |

- **El lapso y el desplazamiento son auxiliares algebraicos.** Sus acoplamientos mutuos llevan una derivada: $E_n\ni\tilde\beta'$ y $E_{\tilde\beta}\ni n'$.
- **$h$ tampoco tiene derivadas propias.** Con el gauge $h=0$ (longitud propia del corte), el campo dinámico es $d$, como en el gauge unitario de la cosmología.
- Comparación con el escalar estático. El gauge unitario fijaba $\ell$, el lapso y el desplazamiento entraban solo vía $F=N^2a'+\dot\beta$ (multiplicador de $W$), y $c$ era auxiliar. Ahora $n$ y $\tilde\beta$ tienen términos cuadráticos propios. **Eliminarlos añade a la energía términos que dependen de las funciones libres** (el análogo del $+W^2/4\kappa$ del fluido de 2g y del vector temporal de 2m).

**Del covariante al ADM.** Con $\phi=q\tau$, $n_\mu=-\phi_\mu/y$, $y=\sqrt{2X}=q/\mathcal N$ y la aceleración $a_\mu=D_\mu\ln\mathcal N$:
$$\phi_{\mu\nu}=-yK_{\mu\nu}-n_\mu D_\nu y-n_\nu D_\mu y+n_\mu n_\nu\dot y,\qquad (\Box\phi)^2-\phi_{\mu\nu}^2=y^2(K^2-K_{ij}^2)+2yK\dot y+2(Dy)^2,\qquad E_4=y^4(K^2-K_{ij}^2).$$
Con $R=R^{(3)}+K_{ij}^2-K^2+2\nabla_\mu(n^\mu K-a^\mu)$ (comprobado en FRW), los términos $yK\dot y$ y $(Dy)^2$ se cancelan. Queda
$$L=A_2+A_3K+A_4(K^2-K_{ij}^2)+B_4R^{(3)},$$
con $A_2=G_2$, $A_3=yG_3-\int_0^yG_3\,dy$ ($-G_3\Box\phi=G_3(yK+\dot y)$ y $G_3\dot y\cong-\big(\int G_3dy\big)K$), $A_4=-G_4+2XG_{4X}+4X^2F_4$ y $B_4=G_4$. Es la forma conocida de Gleyzes–Langlois–Piazza–Vernizzi; aquí se comprueba directamente (§4).

## 4. El sector estático en forma cerrada (`adm_estatico.py`, `energia_estatica.py`)

**Lagrangiana estática.** En la longitud propia $\rho$ del corte, con lapso $\mathcal N$ y desplazamiento $v=N^\rho$:
$$ds^2=-(\mathcal N^2-v^2)d\tau^2+2v\,d\tau d\rho+d\rho^2+r^2d\Omega^2 .$$
Para configuraciones estacionarias, $K_\rho{}^\rho=-v'/\mathcal N$ y $K_\theta{}^\theta=-r'v/(r\mathcal N)$, de modo que
$$K^2-K_{ij}^2=4K_\rho K_\theta+2K_\theta^2=\frac{2r'(v^2r)'}{r^2\mathcal N^2},\qquad K=-\frac{(r^2v)'}{r^2\mathcal N},\qquad R^{(3)}=\frac{2-4rr''-2r'^2}{r^2},$$
$$\boxed{\ \mathcal L_{\rm est}=\mathcal Nr^2A_2-A_3(r^2v)'+\frac{2A_4}{\mathcal N}r'(v^2r)'+\mathcal NB_4\big(2-4rr''-2r'^2\big)\ }$$

| Comprobación | Resultado |
|---|---|
| E-L$[\mathcal L_{\rm cov}-\mathcal L_{\rm est}]$ para $\mathcal N,v,r$ (teoría polinómica, 10 coeficientes simbólicos: $G_2$ cúbica, $G_3$ y $G_4$ cuadráticas, $F_4$ lineal) | $0,0,0$ idénticamente |
| RG ($A_4=-1$, $B_4=1$): Schwarzschild en Painlevé–Gullstrand ($\mathcal N=1$, $v=\sqrt{2M/r}$, $r=\rho$) | $0,0,0$ |

Las perturbaciones estáticas (independientes de $t$, que es lo mismo que de $\tau$) son exactamente las variaciones de este funcional. Su energía es $V=-\int L_2^{\rm est}$ sobre configuraciones que cumplen las ecuaciones de los auxiliares. No depende del corte: es la misma en el estático y en el $\tau$, que es el punto de Babichev et al. (§1).

**Ecuaciones del fondo** ($Q=2A_4/\mathcal N$, $p=\mathcal NA_2$, $\beta=\mathcal NB_4$):
- momento: $E_v=r^2A_3'-2rv\,(Qr')'=0$;
- primera integral ($\rho$ no aparece; es la ligadura de la componente $h$ fijada por el gauge):
  $$H=2\beta(r'^2-1)+4\beta'rr'+Qr'(v^2r)'-pr^2=0 ;$$
- $E_r$; y $E_{\mathcal N}$ se sigue de $dH/d\rho=-(E_rr'+E_{\mathcal N}\mathcal N'+E_vv')$.

**Lema (sin término cúbico).** Si $A_3'\equiv0$ en el fondo (p. ej. $G_3=0$), la ligadura de momento es $v(Qr')'=0$.
- Si $v\not\equiv0$, $Qr'$ es constante, y la constante es cero en la garganta, de modo que $A_4\equiv0$.
- Por tanto: **o $\psi'\equiv0$ ($\phi=qt$) o $A_4\equiv0$** (coeficiente de $K^2-K_{ij}^2$ nulo, el caso degenerado).
- El escalar con $\psi'\neq0$ que "fluye" por la garganta necesita el término cúbico.

**La prueba $d=N$ y la ligadura de momento.** Linealizando $E_v$ con $\nu=0$ ($\delta X=0$) y usando la ligadura del fondo (simbólico, `energia_estatica.py`):
$$\delta E_v=-2v\big[Q\,r^2(d/r)'\big]'-2r(Qr')'\,w .$$
La prueba de las Fases 2g–2j ($w=0$, $d=Nf$) la viola salvo que $v=0$ o $N/r$ sea constante. Si se corrige con $w\neq0$ (cuando $(Qr')'\neq0$), aparece la contribución $+S_w^2/2|k_w|$ de los auxiliares, que depende de las funciones libres.

## 5. No-go para $\phi=qt$ y para $G_3=0$ (`nogo_beta.py`, `caso_qt.py`, `reconstruccion_qt.py`, `modo_cero_beta*.py`)

**La ecuación de $\beta$.** Eliminando $p$ entre $H=0$ y $E_r$, y usando la ligadura de momento (simbólico, exacto):
$$\Big(\frac{\beta'}r\Big)'+\frac{K\beta}{r^3}=-\frac{\big[Q\,r^4(v^2/r^2)'\big]'}{4r^3},\qquad K=1+rr''-r'^2,\qquad \beta=\mathcal NG_4 .$$
- Con $v=0$, el lado derecho se anula y **ni $\mathcal N$ ni ninguna función libre aparecen**.
- Con $G_3=0$, el lema da $v\equiv0$ o $Q\equiv0$, y el lado derecho también se anula.
- En RG ($\beta=N$) es la condición de presión isótropa $G^\rho{}_\rho=G^\theta{}_\theta$.

**Teorema.** Sea un fondo estático con $\phi=qt+\psi$ de una GLPV cuártica con simetría de desplazamiento, con $v\equiv0$ o $A_3'\equiv0$, $r>0$ en todo $\rho\in\mathbb R$ y $r'/r^2\to0$, $\beta'/(r\beta)\to0$ en los dos extremos (asintóticamente plano, $G_4\to G_{4\infty}\neq0$). Entonces $\beta=\mathcal NG_4$ cambia de signo. *Demostración:* si $\beta$ no se anula,
$$0\le\int_{\mathbb R}\frac{\beta'^2}{r\beta^2}\,d\rho=\Big[\frac{\beta'}{r\beta}\Big]-\int\frac{K}{r^3}d\rho=-\Big[\frac{r'}{r^2}\Big]-\int\frac{1+r'^2}{r^3}\,d\rho<0,$$
usando $K/r^3=(r'/r^2)'+(1+r'^2)/r^3$. ∎
- No hace falta simetría de reflexión, ni que la garganta sea única.
- Lo esencial es la topología de dos extremos: en una estrella con centro regular, $r(0)=0$ y la integral no está definida.
- Como $\mathcal N>0$, **$G_4$ cambia de signo en el fondo**. $G_4=B_4$ es el coeficiente de $R^{(3)}$ en el gauge unitario, es decir el término de gradiente de los modos tensoriales en el corte $\phi=$ cte. Donde $B_4=0$, la velocidad tensorial $c_T^2=B_4/(-A_4)$ se anula (acoplo fuerte). A un lado del cero, $c_T^2<0$, salvo que $A_4$ cambie de signo en el mismo punto, y entonces es el término cinético el que se anula.

**Comprobaciones.**

| | Resultado |
|---|---|
| Ecuación de $\beta$ en chorros covariantes exactos con $\psi'=0$, $G_3=0$ (`caso_qt.py`, 3 chorros) | 0 exacto |
| Reconstrucción $r=\sqrt{\rho^2+r_0^2}\,(1+\epsilon e^{-\rho^2/w^2})$, $\epsilon$ afinado para $\beta$ acotada (3 modos cero) | 1–2 nodos en todos |
| Escalón $r=\sqrt{\rho^2+1}+A[\ldots\tanh\ldots]$ con zona $r'>1$ ($K<0$), 9 formas × 66 amplitudes | $\beta$ siempre con $\ge1$ cero |
| Garganta afilada $r'=\tanh(\rho/a)(1+Be^{-\rho^2/b^2})$, $r_0r''(0)$ hasta 50, 9 × 41 formas | $\beta$ siempre con un cero |

En la reconstrucción, el lapso $\mathcal N(\rho)$ es libre y $G_4=\beta/\mathcal N$, $G_2=p/\mathcal N$ como funciones de $X=q^2/2\mathcal N^2$. Residuo de $E_r$ $\sim10^{-15}$; el de $E_{\mathcal N}$ es el error de la derivada numérica ($\le10^{-3}$).

*Por qué falla siempre en la garganta.* En la coordenada $dz=d\rho/r$, $u=\beta/r$ cumple $u_{zz}+(1+2rr''-2r'^2)u=0$. La riccati $y=u_z/u$ diverge (un cero de $\beta$) antes de $z\simeq\pi/2\sqrt{K(0)}$, con $K(0)=1+r_0r''(0)\ge1$, salvo que $K$ se haga negativo antes. Y el teorema dice que no basta.

**Relación con la Fase 2m.** Allí, un vector temporal con $G_3=0$ encuentra la obstrucción local $r''(0)=-X\Theta/(2r_0\mathcal G)$ en la garganta. Aquí la obstrucción es global y más fuerte: $G_4$ debe cambiar de signo.

## 6. El monopolo con $\phi=qt$ (`adm_dinamico.py`, `caso_qt.py`)

**Acción cuadrática dependiente del tiempo** (fondo $v=0$, gauge $h_{\rho\rho}=1$, $\mathcal N=N+\nu$, $v=w$, $r\to r+d$). Exacta salvo derivadas totales: E-L de la diferencia $=0$ para $\nu,w,d$.
$$L_2=L_{\rm est,2}(\nu,d)+A_{3,\mathcal N}\big[\nu(2r\dot d-(r^2w)')+2rN'dw\big]-r(Qr')'w^2+2r(Q\dot d)'w+Q\dot d^2 .$$
La ligadura de momento del fondo es $r^2A_3'=r^2A_{3,\mathcal N}N'=0$. Si $N$ no es constante, $A_{3,\mathcal N}=0$ en el fondo. (Con $N$ constante, $r''=pr/2\beta$ con coeficientes constantes y $r$ crece exponencialmente: no es asintóticamente plano.) Queda:
- **El lapso solo está en la parte estática**, $\tfrac12\kappa\nu^2+\nu S$, con
  $$\kappa=r^2\Big(p_{\mathcal N\mathcal N}-\frac{\beta_{\mathcal N\mathcal N}p_{\mathcal N}}{\beta_{\mathcal N}}\Big),\qquad S=-4\beta_{\mathcal N}\mathcal D[d],\qquad \mathcal D[d]=(rd')'+\frac{1-r'^2-rr''}rd .$$
- **El desplazamiento solo está en la parte cinética.** $w=(Q\dot d)'/(Qr')'$ y
  $$T[\dot d]=\int\Big[Q\dot d^2+\frac{r\,[(Q\dot d)']^2}{(Qr')'}\Big],\qquad V[d]=\int\Big[\frac{S^2}{2\kappa}-2\beta r^2\Big[\Big(\frac dr\Big)'\Big]^2\Big],$$
  y $L_{\rm red}=T-V$, sin términos giroscópicos (el corte $\phi=$ cte es el estático).

**Comprobación covariante** (`caso_qt.py`, chorros exactos con $\psi'=0$, $G_3=G_5=0$). En variables $(n,\tilde\beta,h=0,d)$:
- $E_n$ no tiene $i\omega$ y su parte en $d$ es $-4\beta_{\mathcal N}\mathcal D[d]/N$;
- $E_{\tilde\beta}=-2r(Qr')'\tilde\beta+2i\omega r(Qd'+Q'd)$;
- en $E_d$, los coeficientes de $\omega^2d$, $i\omega\tilde\beta'$ e $i\omega\tilde\beta$ son los de la predicción.

Las tres afirmaciones se cumplen exactamente (aritmética racional) en 3 chorros, junto con la ecuación de $\beta$ de §5. Hay que restringir las ecuaciones del fondo a $\psi=0$ antes de derivar. Si no, los chorros de orden 6 tardan más de 50 min; así tardan 19 s.

**Fantasma en la garganta.** En $\rho=0$: $r'=0$, $r''>0$ y $(Qr')'=Qr''$.
- Los dos coeficientes de $T$ tienen el signo de $Q(0)$, es decir de $A_4(0)$. Para perturbaciones localizadas cerca de la garganta, $\mathrm{signo}\,T=\mathrm{signo}\,A_4(0)$.
- $A_4(0)<0$ (término cinético tensorial $-A_4>0$, el de RG): **el monopolo es un fantasma**.
- $A_4(0)>0$: los modos tensoriales lo son.

*Precisión (Parte II, §12).* El símbolo principal en la garganta da las velocidades radiales del monopolo, $c^2=4\mathcal N_0r_0r''_0\beta_{\mathcal N}^2/(A_4\kappa)$. Con $A_4(0)<0$ el término cinético es negativo, y además:
- si $\kappa(0)>0$ (el gradiente tiene el signo sano), las velocidades son imaginarias: el problema está mal planteado en la garganta;
- si $\kappa(0)<0$, el sector entero tiene el signo cambiado y es un fantasma frente a los gravitones.

Con $G_4=1$, $(Qr')'=-2(r'/N)'=(r/N)\,8\pi(\rho+p)$ es la combinación del NEC. $T$ coincide con el $T[h]$ del fluido de la Fase 2g, $\int[8N\mathcal G^2\dot h'^2/(\tilde\rho+\tilde p)-2G_4N\dot h^2]$: "fluido fantasma ⇒ fantasma". En la GLPV, $A_4$ es libre (vía $F_4$, que no entra en el fondo con $v=0$), pero su signo en la garganta decide qué sector es el fantasma.

**La familia $\delta X=0$** ($\nu=0$, $w=0$; la análoga de la prueba $d=N$):
- Cumple las dos ligaduras si $\mathcal D[d]=0$. Con $d=rf$, la ligadura es $(r^3f')'+rf=0$, una EDO sin funciones libres.
- Su energía es exactamente (simbólico, `energia_estatica.py`)
  $$V=-2\int\beta\,r^2\Big[\Big(\frac dr\Big)'\Big]^2d\rho=-2\int\beta r^2f'^2\,d\rho .$$
- Es la misma forma que la de las Fases 2g–2j: $-2\int N^3\mathcal G(f'-\sigma f)^2=-2\int N\mathcal Gr^2[(d/r)']^2$ con $d=Nf$.
  - Lo que cambia es el coeficiente ($\mathcal G$ allí, $B_4=G_4$ aquí) y la ligadura ($d\propto N$ allí, $\mathcal D[d]=0$ aquí).
- Lejos, $d\simeq c_1+c_2\ln\rho$. El corte de la Fase 2g cuesta $O((\ln R)^2/R)\to0$.
- Si $\beta>0$, la familia tendría energía negativa (taquión si $T>0$). Por §5 eso no puede ocurrir en una garganta. Queda la alternativa fantasma/taquión de arriba.

## 7. Con término cúbico: $G_3\neq0$, $\psi'\neq0$ (`beta_cubico.py`)

Con $v\neq0$ y $A_3'\neq0$, la ligadura de momento $rA_3'=2v(Qr')'$ fija $A_3$ y la ecuación de $\beta$ tiene fuente. Para que la configuración sea simétrica y la teoría no tenga la simetría $\phi\to-\phi$ (la rompe $G_3$), $v$ debe ser impar, con $v(0)=0$.
- La ecuación de $\beta$ es lineal e inhomogénea. Dada $(r,v,Q)$, la solución par acotada es única, sin afinar, y proporcional a $QV_0^2$.
- Exploración con $r=\sqrt{\rho^2+1}$, $v=V_0(\rho/w)e^{-\rho^2/2w^2}$ y $Q=Q_0$ constante (40 casos): $\beta>0$ en todo el eje solo con $Q_0>0$, es decir $A_4>0$.
- *Superado por la Parte II (§11).* Con $Q$ constante, eso no dice nada: con un perfil de $A_4$ concentrado en la garganta hay $\beta>0$ con $A_4<0$ en todo el eje. El análisis del monopolo con flujo está en §12–§14.

## 8. Quíntico

La maquinaria con $G_5(X)$ está generada y comprobada (§2: Noether, $J$, límite $q=0$ frente a la Fase 2j, stealth, modos gauge exactos). No se ha hecho la reducción ADM:
- Con $v=0$, el término $B_5K^{ij}G^{(3)}_{ij}$ entra en la ligadura de momento del fondo ($\propto(B_5/\mathcal N)'G^{(3)\rho\rho}$) y acopla $w$ con $d$ en la parte estática.
- La familia $\delta X=0$ deja de cumplir la ligadura del desplazamiento, como en la Fase 2l la rotación con $G_5$ deforma la prueba $d=N$.

## 9. Conclusiones y siguiente paso

**Lo que queda establecido.**
1. **Maquinaria** para $\phi=qt+\psi$ con métrica estática, cuártica y quíntica, validada con la solución stealth de Babichev–Charmousis y con modos gauge exactos.
2. **El monopolo en el gauge unitario temporal** tiene el lapso y el desplazamiento auxiliares algebraicos. La prueba $d=N$ de las Fases 2g–2j no cumple sus ligaduras: el teorema no se aplica tal cual.
3. **Lagrangiana estática cerrada** en las funciones de la EFT ($A_2,A_3,A_4,B_4$), exacta frente a la covariante.
4. **Lema:** sin término cúbico, una garganta con $\psi'\neq0$ exige $A_4\equiv0$.
5. **No-go:** sin término cúbico (o con $\phi=qt$), en toda garganta estática $G_4$ cambia de signo, para cualesquiera $G_2,G_4,F_4$. Es un resultado de existencia, más fuerte que una inestabilidad.
6. **Monopolo con $\phi=qt$:** $T$ tiene el signo de $A_4(0)$ cerca de la garganta, así que hay fantasma en el monopolo o en el sector tensorial. La familia $\delta X=0$ reproduce la forma $-2\int(\text{lapso})(\text{coef.})\,r^2[(d/r)']^2$ de las fases anteriores.

**El mecanismo común.** Con escalar estático, la garganta pagaba con $\sigma=(\ln r/N)'\neq0$. Con escalar temporal, la ecuación de $\beta=\mathcal NG_4$ es un problema de Sturm–Liouville cuyo potencial integrado es positivo ($\int(1+r'^2)/r^3$). Tampoco se anula: en los dos casos es la geometría de dos extremos.

**Lo que no queda establecido** en la Parte I: el sector con $G_3\neq0$ y $\psi'\neq0$ (tratado en la Parte II) y la reducción ADM quíntica (sigue pendiente).

# Parte II — galileón cúbico con flujo ($G_3\neq0$, $\psi'\neq0$)

## 10. Planteamiento

Es el único hueco que deja la Parte I en la GLPV cuártica con simetría de desplazamiento y métrica estática:
- con $A_3'\neq0$, la ligadura de momento $rA_3'=2v(Qr')'$ permite un desplazamiento $v\neq0$ con $A_4\neq0$;
- la ecuación de $\beta=\mathcal NG_4$ gana una fuente:
  $$\Big(\frac{\beta'}r\Big)'+\frac{K\beta}{r^3}=-\frac{(Qr^4u')'}{4r^3},\qquad u=\frac{v^2}{r^2},\quad K=1+rr''-r'^2,\quad Q=\frac{2A_4}{\mathcal N}.$$

Una configuración simétrica tiene $v$ impar, con $v(0)=0$: $G_3$ rompe $\phi\to-\phi$. Las condiciones sanas que se exigen al fondo son:
- $G_4>0$ y $A_4<0$, es decir ondas gravitatorias en el corte $\phi=$ cte sin fantasma ni inestabilidad de gradiente;
- $\mathcal N>|v|$ (lapso estático real);
- asintóticamente plano, con $G_4\to1$ y $A_4\to-1$.

El programa:
- **(a)** buscar una identidad que fuerce un cambio de signo de $\beta$, como en §5;
- **(b)** si no la hay, reconstruir fondos y estudiar su monopolo.

## 11. (a) Identidades integrales: no hay no-go (`identidad_cubica.py`, `reconstruccion_cubica.py`, `busqueda_cubica.py`)

Dos identidades exactas sobre la ecuación (comprobadas simbólicamente como derivadas totales):
$$\text{(I1)}\qquad \int\frac{K\beta}{r^3}\,d\rho=\frac34\int\Big(\frac vr\Big)^2(Qr')'\,d\rho ,$$
$$\text{(I2)}\qquad \int\frac1r\Big(\frac{\beta'}\beta+\frac{Qr^2u'}{8\beta}\Big)^2+\int\frac{1+r'^2}{r^3}=-\frac34\int\frac{Qu'r'}\beta+\frac1{64}\int\frac{Q^2r^3u'^2}{\beta^2}\qquad(\beta>0).$$
- **(I2)** generaliza la de §5. El último término es positivo sea cual sea el signo de $A_4$: un flujo bastante fuerte paga la geometría de dos extremos. No hay no-go algebraico.
- **(I1)** dice dónde tiene que estar el flujo. En geometrías con $K>0$ (todas las de energía positiva), $\beta>0$ exige $\int(v/r)^2(Qr')'>0$: el flujo debe vivir donde $(Qr')'>0$.
  - $(Qr')'>0$ es justo la condición cinética del monopolo sin flujo (§6).
  - En la garganta, $(Qr')'=Qr''<0$ si $A_4<0$.
  - Por tanto el flujo debe estar en el flanco de un perfil de $A_4$ concentrado en la garganta, donde $Q'r'>|Q|r''$.

**Fondos explícitos** con $G_4>0$ y $A_4<0$ en todo el eje. Geometría de Ellis, $r=\sqrt{\rho^2+1}$; flujo gaussiano; $Q=-2[1+a\,e^{-\rho^2/w_q^2}]$; $V_0$ fijado por $G_{4\infty}=1$.

| familia | fondos | $\beta>0$ | ejemplo |
|---|---|---|---|
| un flujo, $Q$ con pico | 72 | 1 | $w=2$, $a=20$, $w_q=1$: $G_4(0)\simeq0.036$, $v_{\max}=0.47$ |
| dos flujos (ancho + estrecho) | 300 | 1 | — |
| dos flujos + lapso algebraico (`criba_cubica.py`) | 500 | 54 | $w=2$, $w_2=0.2$, $\lambda=0.05$, $a=50$: $\beta(0)=0.083$, $v_{\max}=0.31$ |

(a) no cierra. La parte I §7 (indicio de que $\beta>0$ exige $A_4>0$) queda superada.

## 12. (b) El monopolo con flujo: acción dependiente del tiempo y símbolo principal

**Lagrangiana ADM dependiente del tiempo** (`adm_tiempo_cov.py`). Con $\mathcal N,v,r$ funciones de $(\tau,\rho)$, gauge $h_{\rho\rho}=1$ y $\phi=q\tau$, la forma
$$L=\mathcal Nr^2\big[A_2+A_3K+A_4(4K_\rho K_\theta+2K_\theta^2)+B_4R^{(3)}\big],\qquad K_\rho=-\frac{v'}{\mathcal N},\quad K_\theta=\frac{\dot r-vr'}{r\mathcal N},$$
difiere de la covariante en una derivada total: el E-L de la diferencia es $\equiv0$ para $\mathcal N$, $v$ y $r$.

**Estructura** (`monopolo_cubico.py`: operador lineal exacto con perturbaciones $\nu,w,d$).
- El lapso es algebraico.
- El desplazamiento entra con $w'$ acoplado a $\nu$.
- El bloque auxiliar es $\det M_{\rm aux}\simeq-r^2Z^2k^2/\mathcal N^4$, con
  $$Z=\mathcal N^2\big[rA_{3,\mathcal N}-2r'vQ_{\mathcal N}\big] .$$
- Con flujo, $Z\neq0$ y la estructura principal ($k\to\infty$) no es la del caso sin flujo: $\det M\simeq k^4P(c)$, con $c=\omega/k$ y $P=P_2c^2+P_1c+P_0$.
- El símbolo efectivo de $d$ es $\propto P/D$, normalizado como $2A\omega^2-2Bk^2$. La energía en el corte $\tau$ es positiva si $s_2=P_2/D>0$ y $s_0=P_0/D<0$, y el problema es hiperbólico si $P_1^2\ge4P_2P_0$.

**En la garganta** (`simbolo_garganta.py`, `garganta_cubica.py`).
- $P_1=0$: no hay inclinación.
- Al imponer las ecuaciones del fondo ($H$, ecuación de $\beta$ y ligadura de momento hasta el orden necesario), **todas las segundas derivadas libres de la teoría se cancelan**. Con $R_k=r^{(k)}(0)$, $V_k=v^{(k)}(0)$ y $Q_0=Q(0)$:
  $$\tilde\Gamma=\frac{2\mathcal N_0^2}{R_0^3}\Big[Q_0R_0^2V_1\big(R_0^2R_2V_3-V_1(R_0^2R_4+R_2-R_0R_2^2)\big)-2\beta_0(1+R_0R_2)(R_0^2R_4+R_2)\Big],$$
  $$\hat\Lambda=-\frac{\mathcal N_0}{4R_0^3}\Big[Q_0^2R_0^4V_1^4+2Q_0R_0^3R_2V_1^2\beta_0-4\beta_0^2(1+R_0R_2)^2\Big].$$
- Monopolo sano en la garganta $\iff \tilde\Gamma>0$ (cinética) y $\hat\Lambda/A_4>0$ (gradiente).
- Con $A_4<0$, el gradiente exige un flujo fuerte: $|Q_0|V_1^2\gtrsim2\beta_0(1+R_2)$.
- La cinética exige que el flujo se curve muy cerca de la garganta. Para Ellis, $V_3\lesssim-3V_1$, es decir el máximo del flujo a $\rho\lesssim0.8\,r_0$.
- El límite sin flujo reproduce §6: $c^2=4\mathcal N_0R_0R_2\beta_{\mathcal N}^2/(A_4\kappa)$ (ver la precisión allí).
- Comprobación: en el fondo $w=2$, $a=20$, la fórmula da $c^2=-0.0673$, y la evaluación numérica completa del operador da $c=\pm0.2594i$.

**En un punto cualquiera** (`simbolo_general.py`, `simbolo_perfiles.py`; reconstrucción sustituida y $\beta''$ de su EDO):
$$P_2=\frac{8Q^2r}{\mathcal N'^2}\,\mathfrak B,\qquad \mathfrak B=-2\beta'(rK)'-\frac{Qvr''(r^2W)'}{r}+Q'W\big(r'W-rvr''\big),\qquad W=r^2\Big(\frac vr\Big)',\qquad D=-\frac{4Q^2r^2v^2r''^2}{\mathcal N'^2}.$$
- **El coeficiente cinético $s_2=-2\mathfrak B/(rv^2r''^2)$ no depende del lapso**: es geometría del fondo, del flujo y de $A_4$.
- Las rutas "funciones de la teoría" y "perfiles" coinciden a $2\times10^{-12}$ (`verificar_perfiles.py`).
- Donde el flujo es despreciable ($Z\to0$), la estructura sin flujo de §6 gobierna los números de onda moderados: cinética $\propto r[(Q\dot d)']^2/(Qr')'$, que pide $(Qr')'>0$.

## 13. Reconstrucción y criba (`hiperbolicidad_cubica.py`, `criba_cubica.py`, `criba_masa.py`, `optimiza_masa.py`)

**Dos ejemplos en detalle** (operador completo congelado en 80 puntos, $k=2\times10^3$ y $4\times10^3$):

| fondo | garganta | lejos de la garganta |
|---|---|---|
| Ellis, un flujo ($w=2$, $a=20$) | $c=\pm0.26i$: **mal planteado** en $|\rho|\lesssim0.75$ | hiperbólico |
| Ellis, dos flujos ($w_2=0.2$, $\lambda=0.05$, $a=50$) | $c=\pm0.08$, sano ($\tilde\Gamma=347$, $\hat\Lambda=-234$) | **capa mal planteada** $0.36\lesssim|\rho|\lesssim0.80$: $P_2$ cambia de signo dos veces |

**Criba** (símbolo principal en 420–750 puntos por fondo):

| familia | fondos | $\beta>0$ | hiperbólico en todo $\rho$ | sin fantasma |
|---|---|---|---|---|
| Ellis, dos flujos, $A_4$ con pico, lapso algebraico | 500 | 54 | 3 | 0 |
| cola de energía positiva, $m(r)=M-(M-\tfrac12)r^{-n}$, $A_4\to-1$, lapso con cola $1/\rho$ | 400 | 0 | — | — |

- **Los 3 fondos hiperbólicos de la familia de Ellis** tienen energía positiva en el corte $\tau$ solo en el 34–39 % de los puntos. Lejos, $v\to0$ y $Q\to-2$, de modo que $(Qr')'=-2r''=-2/r^3<0$: la cola de Ellis es exótica y el monopolo tiene allí energía cinética negativa en los dos regímenes ($s_2<0$ y $(Qr')'<0$).
  - Ahí los conos no están inclinados, así que el corte $\tau$ es espacial para el escalar y para los gravitones, y las energías tienen signos opuestos. Es un fantasma también según el criterio de Babichev et al.
- **Con cola de energía positiva** (masa de Misner–Sharp creciente hacia fuera y $(r'/\mathcal N)'<0$ lejos), $\beta$ cambia siempre de signo:
  - o es negativa en la garganta ($\rho\lesssim1$) con $\beta_\infty>0$;
  - o al revés.

  Es el mecanismo de §5: $K(0)\ge1$ hace que $\beta$ caiga antes de que el flujo la sostenga.
- Evolución diferencial sobre esa familia (9 parámetros; objetivo $\min\beta/\max|\beta|$):
  - el mejor valor pasa de $-0.071$ a $-0.043$ en 12 generaciones y a $-0.029$ en 40;
  - el óptimo se va hacia los bordes de la caja: $A_4$ con un pico cada vez más estrecho ($w_q\to0.24$, límite 0.2) y débil ($a\to2.2$, límite 2), y $\mathcal N(0)\to0.85$ (límite 0.9).

  - con la caja ampliada ($w_q\ge0.05$, $a\ge0.3$, $\mathcal N(0)\le0.98$, $M\le2.5$; `optimiza_masa.py 35 ancha`), el mejor valor es $-0.026$. El pico de $A_4$, el flujo y el lapso quedan ya en el interior ($w_q=0.27$, $a=0.63$, $\mathcal N(0)=0.64$); solo la masa sigue en el borde ($M=2.48$).

  Parece una meseta ligeramente negativa que mejora poco al crecer $M$. **No es concluyente**: no se descarta $\beta>0$ con cola de energía positiva, quizá en el límite de masa grande. Si lo hubiera, quedarían por comprobar la hiperbolicidad y el signo cinético.

## 14. Por qué no hay espectro

El criterio de corte de Babichev et al. y el espectro de modos tienen sentido para un problema bien planteado. En los fondos construidos ocurre una de dos cosas:
- **velocidades imaginarias** en una región: entonces $\mathrm{Im}\,\omega\propto k$ sin cota. Es una inestabilidad de gradiente (mal planteamiento), no un modo con tasa finita, y el espectro no está definido;
- **energía cinética negativa con conos sin inclinar** en la cola: un fantasma, que ningún cambio de corte arregla.

El operador está listo para un espectro (`monopolo_cubico.py`: coeficientes exactos en $(\nu,w,d)$ con la frecuencia) si aparece un fondo que pase las condiciones locales.

## 15. Conclusiones de la Parte II y siguiente paso

**Establecido.**
1. Sin término cúbico, la Parte I es concluyente. Con término cúbico y flujo, **no hay no-go algebraico**: (I2) deja que un flujo fuerte pague la garganta, y hay fondos con $G_4>0$ y $A_4<0$ en todo el eje.
2. La lagrangiana ADM dependiente del tiempo con desplazamiento es exacta, y el símbolo principal del monopolo tiene forma cerrada.
   - En la garganta depende solo de datos de Taylor ($\tilde\Gamma$, $\hat\Lambda$).
   - En general, el coeficiente cinético $\mathfrak B$ no depende del lapso.
   - (I1) liga la positividad de $G_4$ con la condición cinética sin flujo, $(Qr')'>0$.
3. **Ninguno de los ≈1300 fondos construidos tiene un monopolo localmente sano.** Los mecanismos son dos:
   - una capa no hiperbólica junto a la garganta, donde $\mathfrak B$ cambia de signo;
   - una cola exótica con monopolo fantasma; con cola de energía positiva no se ha conseguido $G_4>0$ (no descartado: el óptimo está en el borde de la caja).

**Coherencia con la Fase 2m** (v0.2, vector temporal: el límite de Stückelberg del caso sin flujo). La acción reducida tiene la misma forma que la $T$ de §6, $-\tfrac{2\mathcal G}N\dot d^2+k_2((\mathcal G\dot d/N)')^2$, y la misma regla de signos en la garganta: fantasma monopolar o tensorial.

**No establecido.** No es un teorema: la clase de perfiles es infinito-dimensional. Las condiciones locales (hiperbolicidad, signo cinético en los dos regímenes, $G_4>0$, $(Qr')'>0$ lejos) no se han demostrado incompatibles en general.

**Siguiente paso.**
- Un teorema combinando (I1) con el signo de $\mathfrak B$. Por ejemplo, $\int\mathfrak B\,\mu$ con un peso que haga aparecer (I1), o la condición de §12 en la garganta junto con $\beta>0$.
- Si no, dar por cerrada la línea escalar-tensorial de segundo orden con simetría de desplazamiento, con este estado: no-go sin término cúbico, y con él ningún candidato sano. Pasar a redactar el artículo de las Fases 2f–2m.
- La reducción ADM quíntica sigue pendiente.

## Archivos

| archivo | contenido |
|---|---|
| `generar_t.py` | fondo $\phi=qt+\psi$ y acción cuadrática l = 0 (cuártica; `G5` quíntica) → `resultados/t_simbolico.pkl`, `t5_simbolico.pkl` |
| `chorros_t.py` | chorros exactos de la capa de masas con $\phi_t=q$ ($G_2$ reconstruido orden a orden) |
| `verificar_t.py` | stealth de Babichev–Charmousis y modos gauge exactos (`G5`: quíntica) |
| `estructura_t.py` | pares de chorros de la acción cuadrática en el gauge unitario (variables estáticas) |
| `adm_t.py` | ecuaciones en las variables ADM del corte: lapso y desplazamiento algebraicos |
| `adm_estatico.py` | lagrangiana estática cerrada; E-L de la diferencia con la covariante $\equiv0$; Painlevé–Gullstrand |
| `energia_estatica.py` | familia $\delta X=0$: $-L_2=-2\beta r^2[(d/r)']^2+$ derivada total; $\delta E_{\mathcal N}=-4\beta_{\mathcal N}\mathcal D[d]$; ligadura de momento linealizada |
| `adm_dinamico.py` | acción cuadrática dependiente del tiempo con $v=0$; E-L de la diferencia $\equiv0$ |
| `caso_qt.py` | comprobación covariante (chorros exactos, $\psi'=0$) de la estructura ADM, de $T$ y de la ecuación de $\beta$ |
| `nogo_beta.py` | ecuación de $\beta$ con la fuente del desplazamiento |
| `reconstruccion_qt.py`, `modo_cero_beta.py`, `modo_cero_beta2.py` | modos cero de $\beta$ en tres familias: siempre con nodos |
| `beta_cubico.py` | $\beta$ con término cúbico y desplazamiento ($Q$ constante; superado por la Parte II) |
| *Parte II* | |
| `identidad_cubica.py` | identidades (I1) e (I2) de la ecuación de $\beta$ con fuente |
| `reconstruccion_cubica.py`, `busqueda_cubica.py` | fondos de Ellis con $G_4>0$, $A_4<0$; condición de garganta sana |
| `adm_tiempo_cov.py` | lagrangiana ADM dependiente del tiempo con desplazamiento frente a la covariante (E-L de la diferencia $\equiv0$) |
| `monopolo_cubico.py` | operador lineal exacto del monopolo con flujo y reconstrucción del fondo |
| `simbolo_garganta.py`, `garganta_cubica.py` | símbolo principal en la garganta; $\tilde\Gamma$, $\hat\Lambda$ con el fondo impuesto |
| `simbolo_general.py`, `simbolo_perfiles.py`, `verificar_perfiles.py` | símbolo en un punto cualquiera; forma de perfiles ($\mathfrak B$, $D$); comprobación entre rutas |
| `hiperbolicidad_cubica.py` | velocidades características a lo largo del fondo (operador completo congelado) |
| `criba_cubica.py`, `criba_masa.py`, `optimiza_masa.py` | cribas (Ellis; cola de energía positiva) y evolución diferencial |

## Reproducibilidad

```
cd fase2k
python generar_t.py          # ~17 min -> resultados/t_simbolico.pkl
python generar_t.py G5       # ~44 min -> resultados/t5_simbolico.pkl
python verificar_t.py        # stealth y modos gauge (cuartica); con G5, quintica
python adm_estatico.py       # lagrangiana estatica ADM frente a la covariante (segundos)
python energia_estatica.py   # familia dX = 0
python adm_dinamico.py       # accion cuadratica dependiente del tiempo (v = 0)
python nogo_beta.py          # ecuacion de beta
python caso_qt.py            # comprobacion covariante del caso phi = qt y de la ecuacion de beta (~20 s)
python adm_t.py              # estructura auxiliar en las variables ADM (chorros, lento)
python reconstruccion_qt.py; python modo_cero_beta.py; python modo_cero_beta2.py; python beta_cubico.py
# Parte II
python identidad_cubica.py; python adm_tiempo_cov.py                 # identidades; ADM con tiempo (segundos)
python simbolo_garganta.py; python garganta_cubica.py                 # simbolo en la garganta (~10 min)
python simbolo_general.py; python simbolo_perfiles.py                 # simbolo general (necesita monopolo_cubico.operador)
python verificar_perfiles.py
python hiperbolicidad_cubica.py; python hiperbolicidad_cubica.py 2.0 50 1.0 1.0 2.0 0.2 0.05
python reconstruccion_cubica.py; python busqueda_cubica.py 300
python criba_cubica.py 500   # ~1 h
python criba_masa.py 400     # ~1 h
python optimiza_masa.py 40   # ~40 min; con "35 ancha", caja ampliada
```
