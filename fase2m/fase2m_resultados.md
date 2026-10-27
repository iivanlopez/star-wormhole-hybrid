# Híbrido estrella–agujero de gusano — Fase 2m: gargantas sostenidas por un campo vectorial (Proca generalizada)

*Documento de trabajo v0.1 · 2026-10-27 · continúa [fase2j/fase2j_resultados.md](../fase2j/fase2j_resultados.md) (§6, salida 3: "campos vectoriales")*

## Resumen

1. **La pregunta.** En la GLPV cuártica (2g) y en Horndeski quíntico (2j), toda garganta estática con el sector impar sano tiene un taquión en el monopolo. La prueba es $d=N$: reescalar el radio areal con el lapso, sin tocar $g_{\ell\ell}$ ni el escalar, con energía $E[N]=-2\int N^3\mathcal G\sigma^2$. ¿Sobrevive con un campo vectorial? En Proca generalizada el sector $l=0$ tiene, además de la métrica, las dos componentes $A_0,A_1$ del vector.

2. **Maquinaria** (§2). Usamos Proca generalizada cuártica,
   $$L=G_2(X,F,Y)+G_3(X)\nabla_\mu A^\mu+G_4(X)R+G_{4,X}\big[(\nabla_\mu A^\mu)^2-\nabla_\rho A_\sigma\nabla^\sigma A^\rho\big],$$
   en simetría esférica con $A_\mu=(A_0,A_1,0,0)(t,\ell)$. En $l=0$ se cumple $Y=4XF$ idénticamente, así que basta $G_2(X,F)$. Comprobaciones exactas:
   - identidades de Noether radial y temporal;
   - Reissner–Nordström, de forma simbólica;
   - con $A_0=0$, el fondo y la acción cuadrática son **idénticos** a los de Horndeski de la Fase 2f;
   - paridad temporal;
   - los modos gauge resuelven las seis ecuaciones lineales en chorros exactos de la capa de masas.

3. **Vector radial ($A_0=0$): el déficit sobrevive exactamente** (§3). El fondo es el de Horndeski con $\phi'=A_1$. En el gauge $\delta A_1=0$ ocurren tres cosas:
   - la prueba $d=N$ con $\delta A=0$ cumple todas las ligaduras;
   - su energía es $-2\int N^3\mathcal G\sigma^2$ con $\mathcal G=G_4-2XG_{4X}$, sin depender de $G_2(X,F)$ ni de $G_3$;
   - $\delta A_0$ solo se acopla a **velocidades**. No cambia la energía potencial: solo rebaja la cinética.

   **El teorema del monopolo vale para todo vector radial.** (El argumento se extiende a $G_5,G_6,g_5$, porque con $\delta A=0$ el vector sigue siendo exacto.)

4. **Vector temporal ($A_1=0$, $G_3=0$): el vector compensa el déficit** (§4). La estructura exacta es
   $$L_2=\kappa F^2-FW+L_{\rm rest},\qquad F=N^2a'+\dot\beta,\qquad \kappa=\frac{U^2r^2(G_{2F}+2FG_{2FF})}{2N^5}>0 .$$
   - La prueba $d=N$ **no** cumple la ligadura de momento: $W_N=\frac{2Ur}{N^2}\big[2G_{4X}N(U/N)'-G_{2F}U'\big]\neq0$.
   - Eliminando $F$ y la perturbación algebraica $\tilde a=\delta A_0-Ua$ queda, exactamente,
     $$E[N]=\int\Big[-2N^3G_4\sigma^2+\frac{2N}{g}\big(2G_{4X}N(U/N)'-G_{2F}U'\big)^2+\frac{J_N^2}{4\tilde m}\Big]d\ell ,\qquad g=G_{2F}+2FG_{2FF}.$$
     El segundo término es la energía del campo eléctrico perturbado y es positivo.
   - En las seis gargantas regulares que construimos, ese término supera al déficit del vacío por un factor de 1.2 a 4.6.

   **El déficit $d=N$ no sobrevive en general.**

5. **Pero aparece una obstrucción nueva, exacta y sin funciones libres** (§4.4). En una garganta simétrica:
   $$r''(0)=-\frac{X\,\Theta}{2r_0\mathcal G},\qquad \tilde c(0)=-\frac{N_0\Theta}{2}=\frac{N_0r_0\mathcal G}{X}\,r''(0),\qquad \Theta\equiv G_{2X}r_0^2+2G_{4X}.$$
   Aquí $\tilde c$ es el coeficiente de $\tilde w^2$, con $\tilde w=\delta A_1+U\beta/N^2$ la perturbación auxiliar impar: localmente, la masa de la componente radial del vector con el signo cambiado.
   - **Que haya garganta ($r''(0)>0$, $\mathcal G>0$) obliga a $\tilde c(0)>0$:** la componente radial tiene masa local negativa en la garganta.
   - Al eliminar $\tilde w$, la energía cinética del monopolo recibe $-J_v^2/4\tilde c<0$. Si el vector es masivo lejos, $\tilde c<0$ allí, así que $\tilde c$ se anula en alguna esfera, donde el problema reducido es singular.
   - Además, la masa estática $\tilde m$ de $\tilde a$ cambia de signo en las 6 gargantas regulares, con $\tilde m(0)<0$.
   - Con $G_4$ constante se reproduce la dicotomía de la k-esencia temporal: la garganta necesita $G_{2X}<0$.

6. **Literatura** (§1, verificada con la API de arXiv). No encontramos gargantas en Proca generalizada ni estudios de su monopolo. Las referencias relevantes:
   - agujeros negros con pelo de Proca generalizada y sus perturbaciones impares y pares: [1705.09662](https://arxiv.org/abs/1705.09662), [1801.01787](https://arxiv.org/abs/1801.01787), [2409.15606](https://arxiv.org/abs/2409.15606);
   - desestabilización por acoplos no mínimos: [2104.08049](https://arxiv.org/abs/2104.08049);
   - gargantas con un vector masivo mínimamente acoplado, donde los autores señalan que hace falta un fantasma: [2412.19408](https://arxiv.org/abs/2412.19408).

---

## 1. Planteamiento y literatura

**Por qué un vector.** La prueba $d=N$ usa tres hipótesis (2j §6): un escalar estático (gauge unitario sin tiempo), un fondo estático y que $c=\delta g_{\ell\ell}/2$ sea auxiliar. Un vector $A_\mu$ tiene dos componentes en $l=0$. La componente radial $A_1$ es la análoga de $\phi'$. La temporal $A_0$ no tiene análogo escalar estático: en el lenguaje de Stückelberg, $A=\tilde A+d\pi$, es a la vez un potencial eléctrico y un $\dot\pi$ (la salida 1 de 2j, el "escalar lineal en el tiempo", en versión vectorial). En $l=0$ sigue habiendo un único grado de libertad, el modo longitudinal: seis campos, menos dos gauge, menos dos multiplicadores, menos $A_0$ auxiliar.

**Paridad temporal.** Bajo $t\to-t$, $A_0\to-A_0$.
- Con $A_0=0$ (caso R) el fondo es invariante bajo T. Los campos impares son $\beta$ y $\delta A_0$, y se desacoplan de los pares en el sector estático.
- Con $A_1=0$ y $G_3=0$ (caso E) la teoría es invariante bajo $A\to-A$ (C), y el fondo lo es bajo TC. Los impares son $\beta$ y $\delta A_1$.
- Con $A_0,A_1\neq0$ no hay ninguna de las dos simetrías: la acción cuadrática tiene 463 términos que las rompen, y el monopolo puede tener términos giroscópicos. Este caso no lo tratamos (§6).

**Literatura** (todas las citas verificadas con la API de arXiv, `arxiv.py`):
- Heisenberg, Kase, Minamitsuji, Tsujikawa, *Hairy black-hole solutions in generalized Proca theories*, [arXiv:1705.09662](https://arxiv.org/abs/1705.09662), PRD 96 (2017) 084049.
- Kase, Minamitsuji, Tsujikawa, Zhang, *Black hole perturbations in vector-tensor theories: the odd-mode analysis*, [arXiv:1801.01787](https://arxiv.org/abs/1801.01787), JCAP 02 (2018) 048. Las soluciones con $A_1\neq0$ procedentes de $G_4\supset\beta_4X^n$ son inestables cerca del horizonte.
- De Felice, Kase, Tsujikawa, *Scrutinizing black hole stability in cubic vector Galileon theories*, [arXiv:2409.15606](https://arxiv.org/abs/2409.15606): inestabilidades laplacianas angulares de los agujeros negros con $A_0,A_1\neq0$.
- García-Saenz, Held, Zhang, *Destabilization of black holes and stars by generalized Proca fields*, [arXiv:2104.08049](https://arxiv.org/abs/2104.08049), PRL 127 (2021) 131104.
- Ganiyeva, Barros, de la Cruz-Dombriz, Lobo, *Vector wormholes as conduits for matter interaction*, [arXiv:2412.19408](https://arxiv.org/abs/2412.19408), EPJC 85 (2025) 914: gargantas con un vector masivo mínimamente acoplado. Los autores señalan que hace falta al menos un grado de libertad fantasma.
- Su, Hao, Ren, Wang, [arXiv:2311.17557](https://arxiv.org/abs/2311.17557): estrellas de Proca con una garganta sostenida por un escalar fantasma (no por el vector).

Ninguno estudia el monopolo de una garganta sostenida por el vector.

## 2. La maquinaria (`proca_simbolico.py`, `generar_proca.py`, `verificar_proca.py`)

La métrica esférica general $(A,b,B,C)$ es la de las Fases 2f–2j, con $A_\mu=(u,w,0,0)(t,\ell)$. El fondo es $A=N^2$, $B=1$, $C=r^2$, $u=U(\ell)$, $w=V(\ell)$, con
$$X_b=\frac{U^2}{2N^2}-\frac{V^2}2,\qquad F_b=\frac{U'^2}{2N^2}.$$
Las perturbaciones son $\delta A=2N^2a$, $\delta b=\beta$, $\delta B=2c$, $\delta C=2rd$, $\delta A_0=\alpha$, $\delta A_1=w$. Las derivadas usan la regla de la cadena con dos invariantes ($X$, $F$). Las ecuaciones del fondo tienen entre 12 y 30 términos; la acción cuadrática, 1242.

| Comprobación | Resultado |
|---|---|
| Independencia de $\theta$ | $2\times10^{-165}$ |
| $Y=4XF$ en $l=0$ | identidad simbólica |
| Noether radial: $E_AA'+E_CC'+E_uU'+E_wV'-2E_B'-(VE_w)'=0$ | idénticamente nula |
| Noether temporal ($t\to t+T(\ell)$): $N^2E_b=UE_w$ | idénticamente nula |
| Reissner–Nordström ($G_4=1$, $G_2=F$, $A_0=2Q/r$) | las seis ecuaciones se anulan (simbólico en $r,M,Q$) |
| $A_0=0$: $E_A,E_B,E_C$ y $E_p=-D_\ell E_w$ frente a la Fase 2f ($\phi'=A_1$, $G_3^H=-G_3$) | idénticas (aritmética exacta) |
| $A_0=0$: acción cuadrática con $\delta A=0$ frente a la de la Fase 2f en el gauge unitario | idénticas |
| $A_0=0$: $E_u\equiv0$; paridad T en $L_2$ | 0 términos la violan |
| $A_1=0$, $G_3=0$: paridad TC en $L_2$ | 0 términos la violan |
| Modos gauge $\mathcal L_\xi$, $\xi=(T,Z)e^{i\omega t}$, en las seis ecuaciones lineales | exactos (2 chorros de la capa de masas por caso: general, R y E) |

**Chorros de la capa de masas** (`chorros_proca.py`). $E_w$ y $E_B$ fijan $G_{2,10}$ y $G_{2,00}$. En cada orden se resuelven $(N_k,r_k,U_k,G_{2,k0})$, con $A_1$ como dato libre: es la reconstrucción de la Fase 2f.

*Nota técnica.* `sp.nsimplify` aproxima los racionales grandes, y por eso un modo gauge daba un residuo de $3\times10^{-18}$ en lugar de 0. Aquí se usa el racional exacto. Los chorros de 2f y 2j usaban `nsimplify`, pero sus comprobaciones eran de orden más bajo.

## 3. Vector radial ($A_0=0$): el teorema sobrevive (`caso_R.py`)

**Fondo.** Es exactamente el de Horndeski con $\phi'=A_1$, $G_2^H(X)=G_2(X,0)$, $G_3^H=-G_3$ y $G_4^H=G_4$. Las gargantas existen en la misma medida que en Horndeski.

**Gauge $\delta A_1=0$.** Es el análogo del unitario: la reparametrización radial actúa como $\delta A_1\to\delta A_1+(ZA_1)'$. Quedan $a,\beta,c,d$ y $\alpha=\delta A_0$, y valen tres resultados exactos:

1. **Identidad de la familia $c=0$** (fuera de la capa de masas, con $G_2(X,F)$ y $G_3$ cualesquiera):
   $$-L_{\rm rest}(c=0,d=Nf)\big|_{\rm est}=-2N^3\mathcal G\,(f'-\sigma f)^2,\qquad \mathcal G=G_4+A_1^2G_{4X}=G_4-2XG_{4X},\qquad\sigma=\Big(\ln\frac rN\Big)' .$$
2. **Ligaduras.** La prueba $d=N$, $c=a=\beta=\alpha=0$ cumple $E_a=E_\beta=E_\alpha=0$ (3 chorros, exacto). Las ecuaciones de $\alpha$ y $\beta$ se cumplen por la paridad T: en el sector estático, los campos impares solo tienen fuentes impares.
3. **$\alpha$ no toca la energía potencial.** Sus acoplos con $a$ son $\kappa_a(a\dot\alpha+\dot a\alpha)=\kappa_a\partial_t(a\alpha)$. Con $\beta$ son $\mu\beta\alpha+\nu(\beta\alpha)'$, con $\nu'=\mu$ en la capa de masas. El resto es
   $$L_\alpha=k_F\alpha'^2+m_\alpha\alpha^2+n\,\alpha\alpha'+\alpha\,J[\dot c,\dot d,\dot c',\dot d'],\qquad k_F=\frac{G_{2F}r^2}{2N}.$$
   Así que
   $$L_2=-FW+L_{\rm rest}(c,d)+L_\alpha .$$
   Eliminar $\alpha$ da $-\tfrac14\langle J,M^{-1}J\rangle$, con $M=-\partial(k_F\partial)+m_\alpha$: es una contribución **solo cinética**, no positiva si $M>0$.

**Teorema (vector radial).** Sea una garganta estática, simétrica y asintóticamente plana de Proca generalizada cuártica ($G_2(X,F,Y)$, $G_3(X)$, $G_4(X)$) con $A_0=0$ y $\mathcal G>0$. Si el monopolo reducido no tiene fantasma ($K_{\rm eff}>0$), tiene un taquión:
$$\omega_0^2\le\frac{-2\int N^3\mathcal G\sigma^2}{\langle N,K_{\rm eff}N\rangle}\le\frac{-2\int N^3\mathcal G\sigma^2}{\int K_HN^2}<0 \qquad(\text{la segunda cota, si }M>0).$$
El vector no puede compensar el déficit, porque solo se acopla a velocidades.
- **Términos de orden superior.** Con $G_5$, $G_6$ y $g_5$, la configuración $\delta A=0$ deja $A$ exacto. $L_6$ y $g_5$ dependen de $F$ al menos cuadráticamente y no contribuyen. $L_5$ se reduce al de Horndeski, y vale la identidad quíntica de 2j con su $\mathcal G$. La paridad T no depende de la teoría.
- **Lo que solo argumentamos.** No hemos calculado el sector impar de Proca. Restringido a $\delta A_{\rm impar}=0$ coincide con el de Horndeski, cuyo coeficiente cinético es $\mathcal G$, pero la condición completa de Proca (Kase et al. 2018) mezcla además el modo vectorial impar.

## 4. Vector temporal ($A_1=0$, $G_3=0$): `reduccion_E.py`, `energia_E.py`, `longitudinal_E.py`

### 4.1 Estructura

**Variables invariantes** bajo las reparametrizaciones temporales: $\tilde a=\alpha-Ua$ y $\tilde w=w+U\beta/N^2$. En ellas (exacto en chorros de la capa de masas; el resto es una derivada total, y sus ecuaciones de Euler–Lagrange se anulan):
$$L_2=\kappa F^2-FW+L_{\rm rest}(c,d,\tilde a,\tilde w),\qquad \kappa=\frac{U^2r^2\,g}{2N^5},\quad g=G_{2F}+2FG_{2FF},$$
$$W=\frac r N\Big[\Big(\frac{gUU'r}{N^2}+4G_4r'\Big)c-4G_4d'\Big]+\frac{2r}{N^2}\Big[2G_4N'-\frac{G_{2F}UU'}{N}+2G_{4X}U\Big(\frac UN\Big)'\Big]d-\frac{Ur^2g}{N^3}\big(\tilde a'-\dot{\tilde w}\big)-\frac{G_{2XF}U^2U'r^2}{N^5}\tilde a$$
(`reduccion_E_log.txt`).
- A diferencia del vacío, $F$ entra **cuadráticamente**, con $\kappa>0$ si el campo eléctrico no es fantasma ($g>0$). La estructura es la del fluido de la Fase 2g, con $\kappa F^2$ en el papel de $\kappa_aa^2$.
- Las ecuaciones de $a$ y $\beta$ dan $\Pi=2\kappa F-W=0$. Entonces
  $$L_{\rm eff}=-\frac{W^2}{4\kappa}+L_{\rm rest}.$$
  La energía eléctrica de la perturbación entra con signo **positivo**, $+W^2/4\kappa$.
- En $L_{\rm eff}$ se cancelan los términos $\tilde a'^2$ y $\dot{\tilde w}^2$. **$\tilde a$ (par) y $\tilde w$ (impar) son auxiliares algebraicos**, con coeficientes $\tilde m\,\tilde a^2$ y $\tilde c\,\tilde w^2$.

### 4.2 La prueba $d=N$

Con $c=\tilde a=\tilde w=0$ y $d=N$ (exacto, fuera de la capa de masas):
$$W_N=\frac{2Ur}{N^2}\Big[2G_{4X}N\Big(\frac UN\Big)'-G_{2F}U'\Big]\neq0,\qquad -L_d(N)=-2N^3G_4\sigma^2+\frac{W_N^2}{4\kappa}=-2N^3G_4\sigma^2+\frac{2N}{g}\Big[2G_{4X}N\Big(\frac UN\Big)'-G_{2F}U'\Big]^2 .$$
- El término de vacío lleva $G_4$, no $\mathcal G$: la diferencia $-2G_{4X}U^2N\sigma^2$ la reabsorbe el vector.
- La ecuación estática de $\tilde a$ es algebraica, $2\tilde m\tilde a+J=0$, con fuente $J_N\neq0$ en la capa de masas. Eliminando $\tilde a$ (en el gauge $c=0$):
  $$E[N]=\int\Big[-2N^3G_4\sigma^2+\frac{2N}{g}\Big(2G_{4X}N\Big(\frac UN\Big)'-G_{2F}U'\Big)^2+\frac{J_N^2}{4\tilde m}\Big]d\ell .$$
- El término eléctrico es positivo. El de $\tilde a$ tiene el signo de $\tilde m$ y es singular donde $\tilde m=0$. **La prueba ya no da un signo universal.**
- *Cola lejana.* En el gauge $c=0$ la prueba no se puede cortar barato: donde $U\to0$, la penalización $W^2/4\kappa$ obliga a $(d/N)'\to0$. Pero $d=N\to r'$ lejos es una traslación casi gauge (como en 2i), y la energía converge: $-2N^3G_4\sigma^2\sim-2/\ell^2$.

### 4.3 Gargantas eléctricas (`garganta_E.py`, `analiza_E.py`)

**Familia y método.**
- Teorías: $G_4=1+bX$, $G_2=F+h(X)$ con $h=m^2X+\lambda_2X^2+\lambda_3X^3$ y $m^2=1$ (Proca masivo lejos).
- En la garganta ($r'=N'=U'=0$), $E_B$ da $h(X_0)=-2G_4(X_0)/r_0^2$. Con $N_0=1$, el parámetro de tiro es $X_0=U_0^2/2$: se ajusta para eliminar el modo de Yukawa creciente.
- Barrido de $b\in\{0,-0.5,-1,-2\}$, $\lambda_2\in\{-1,-3,-10\}$, $\lambda_3\in\{0,1,5\}$ (36 teorías), con refinado por bisección.
- Resultado: 25 soluciones, 16 de ellas gargantas de verdad ($r''(0)>0$, $r$ creciente, $U\to0$). De las 16, 10 tienen $N_0/N_\infty\sim10^5$–$10^{14}$, es decir, están al borde de un horizonte, y no las usamos. Quedan 6 regulares:

| $(b,\lambda_2,\lambda_3)$ | $X_0$ | $r_0$ | $N_0/N_\infty$ | $r''(0)$ | $\Theta$ | $\tilde c(0)$ | $G_{2X}+2XG_{2XX}$ | $\tilde m(0)$ | $E_{\rm vol}$ | $E_{\rm ele}$ | $E_{\rm vol}+E_{\rm ele}$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| $(0,-1,0)$ | 1.469 | 1.70 | 5.06 | 2.43 | $-5.6$ | 2.8 | $-7.8$ | $<0$ | $-175.2$ | 205.8 | $+30.6$ |
| $(0,-3,1)$ | 0.419 | 7.64 | 0.82 | 1.58 | $-57.7$ | 28.9 | $-3.9$ | $<0$ | $-0.157$ | 0.295 | $+0.138$ |
| $(-0.5,-3,0)$ | 0.367 | 6.60 | 1.03 | 1.26 | $-53.4$ | 26.7 | $-5.6$ | $<0$ | $-0.299$ | 1.387 | $+1.09$ |
| $(-0.5,-3,1)$ | 0.439 | 5.35 | 1.19 | 1.05 | $-31.2$ | 15.6 | $-4.0$ | $<0$ | $-0.505$ | 2.121 | $+1.62$ |
| $(-1,-3,0)$ | 0.387 | 4.45 | 1.72 | 0.88 | $-28.1$ | 14.1 | $-6.0$ | $<0$ | $-1.52$ | 6.86 | $+5.34$ |
| $(-1,-3,1)$ | 0.479 | 3.23 | 3.62 | 0.72 | $-14.4$ | 7.2 | $-4.2$ | $<0$ | $-15.1$ | 35.2 | $+20.1$ |

(Las energías están normalizadas con $N_\infty=1$; $E_{\rm vol}$ incluye la cola $-2N_\infty^3/r(L)$. En todas las filas $\mathcal G\ge1$ y $G_4\ge1$.)

- **El término eléctrico supera al déficit en todas:** $E_{\rm ele}/|E_{\rm vol}|=1.18$–$4.6$. Sin la corrección de $\tilde a$, la prueba $d=N$ tiene energía positiva.
- **$\tilde m$ cambia de signo en las 6, con $\tilde m(0)<0$.** La corrección $J_N^2/4\tilde m$ no es integrable, y la eliminación estática de $\tilde a$ degenera en la esfera $\tilde m=0$. Por eso no damos $E[N]$ completo: el problema reducido es singular.

### 4.4 La obstrucción en la garganta (exacta)

En una garganta simétrica, con las ecuaciones del fondo y $E_B$ (`longitudinal_E.py`, simbólico):
$$r''(0)=-\frac{X\Theta}{2r_0\mathcal G},\qquad \tilde c(0)=-\frac{N_0\Theta}{2},\qquad \Theta=G_{2X}r_0^2+2G_{4X},$$
$$\tilde m(0)=\frac{G_4\big[(G_{2X}+2XG_{2XX})r_0^2+2(G_{4X}+2XG_{4XX})\big]+4X^2r_0^2\big(G_{2X}G_{4XX}-G_{2XX}G_{4X}\big)}{2N_0\mathcal G}.$$
(Vale para $G_2(X,F)$ y $G_4(X)$ cualesquiera. $\tilde c=c_{\tilde w\tilde w}-\tfrac12e_{\tilde w\tilde w}'$ incluye la parte $\tilde w\tilde w'$.)

**Lectura.**
- **La garganta obliga a $\tilde c(0)>0$.** Con $\mathcal G>0$, $r''(0)>0$ equivale a $\Theta<0$, y eso a $\tilde c(0)=N_0r_0\mathcal Gr''(0)/X>0$.
- **Qué significa.** En Proca plano $\tilde c=-m^2r^2N/2<0$: $\tilde c(0)>0$ es una masa local negativa de la componente radial del vector. En el límite local (coeficientes congelados, $\tilde a\sim\dot\pi$, $\tilde w\sim\pi'$), el par auxiliar da
  $$\omega^2=-\frac{\tilde c}{\tilde m}k^2-\frac{2\tilde c}{g'}\qquad(g'>0\text{ el coeficiente de Maxwell}).$$
  Con $\tilde c>0$ el modo es taquiónico a longitudes largas y, si además $\tilde m>0$, también a todas las demás. Si $\tilde m<0$, el hamiltoniano local es indefinido.
- **En el monopolo con gravedad.** Eliminar $\tilde w$ suma $-J_v^2/4\tilde c$ a la energía cinética, negativa donde $\tilde c>0$. Si el vector es masivo lejos, $\tilde c$ cambia de signo, y el problema reducido es singular en esa esfera (salvo que $J_v$ se anule allí). Es la versión cinética de la obstrucción de la carga de Maxwell (2g §7), donde $e_0(0)\neq0$ fijaba $\delta r(0)=0$.
- **Con $G_4$ constante.** $\Theta<0$ exige $G_{2X}<0$, y $\tilde m(0)\propto G_{2X}+2XG_{2XX}$. Es la dicotomía de la k-esencia con gradiente temporal ($c_s^2=P_X/(P_X+2XP_{XX})$): fantasma o inestabilidad de gradiente. En las 6 gargantas regulares, $G_{2X}+2XG_{2XX}<0$.

**Lo que no hemos demostrado.** Que $\tilde c(0)>0$ implique un fantasma o un taquión del monopolo completo. Haría falta la reducción completa: el término cinético $K_{\rm eff}$ con $\tilde w$ y $\tilde a$ eliminados, y el espectro. Sí demostramos que la garganta fuerza el signo malo del coeficiente local y que, con un vector masivo, el problema reducido degenera en una esfera.

## 5. Conclusiones

**Lo que queda establecido.**
1. **Maquinaria de Proca generalizada cuártica** en $l=0$, validada con Reissner–Nordström, con la correspondencia exacta con Horndeski de la Fase 2f y con los modos gauge exactos.
2. **Vector radial ($A_0=0$).** El teorema del monopolo se extiende sin cambios: $E[N]=-2\int N^3\mathcal G\sigma^2$, independiente de $G_2(X,F,Y)$ y de $G_3$. $\delta A_0$ solo se acopla a velocidades.
3. **Vector temporal ($A_1=0$).**
   - Estructura exacta $\kappa F^2-FW+L_{\rm rest}$, con la energía eléctrica $+W^2/4\kappa$.
   - La prueba $d=N$ no cumple la ligadura, y su energía corregida contiene $+\int(2N/g)[2G_{4X}N(U/N)'-G_{2F}U']^2$, que en las gargantas construidas supera el déficit del vacío.
   - **El déficit $d=N$ no es universal con un vector temporal.**
4. **Obstrucción nueva, exacta y sin funciones libres.** Una garganta eléctrica simétrica tiene $\tilde c(0)=N_0r_0\mathcal Gr''(0)/X>0$: masa local negativa de la componente radial del vector.
   - Con un vector masivo, el problema reducido degenera en la esfera $\tilde c=0$.
   - En todas las gargantas regulares construidas, también $\tilde m$ cambia de signo.

**El mecanismo común (2g–2m).**
- Mientras la garganta la sostenga un campo "radial" (un escalar estático, o $A_1$), el déficit $d=N$ es inevitable.
- Si la sostiene una componente temporal ($A_0$), la energía eléctrica puede compensarlo. Pero violar la condición nula en la garganta ($\Theta<0$) obliga a que el propio sector longitudinal tenga allí el signo malo.
- Es la misma tensión que en la k-esencia temporal. La salida 1 de 2j (el escalar lineal en el tiempo, Fase 2k) probablemente se encuentre con la versión escalar de esta obstrucción.

**Lo que no queda establecido.**
- El caso mixto ($A_0,A_1\neq0$, el de casi todos los agujeros negros con pelo de Proca): sin paridad temporal y con posibles términos giroscópicos. No lo tratamos.
- La reducción completa del caso E ($K_{\rm eff}$ con $\tilde a$, $\tilde w$ y el espectro), necesaria para convertir la obstrucción local en un enunciado sobre el monopolo.
- Los términos $G_5$, $G_6$, $g_5$, más allá del argumento del caso R.
- El barrido de §4.3 es pequeño (36 teorías). La obstrucción de §4.4 no depende de él.

**Siguiente paso.**
1. Reducir el caso E por completo: eliminar $\tilde w$ con sus fuentes cinéticas, obtener $K_{\rm eff}$ y $P_{\rm eff}$, y estudiar la esfera $\tilde c=0$.
2. Buscar teorías con $\tilde m>0$ en toda la solución (por ejemplo, $G_{4X}>0$ con $G_{2X}r_0^2<-2G_{4X}$ y $G_{2XX}$ grande). Si existieran, el único obstáculo sería $\tilde c(0)>0$.
3. Comparar con el resultado de la Fase 2k (escalar lineal en el tiempo), que es el límite de Stückelberg del caso E.

## Archivos

| archivo | contenido |
|---|---|
| `arxiv.py` | consulta a la API de arXiv (verificación de citas) |
| `proca_simbolico.py` | lagrangiana covariante de Proca generalizada cuártica en simetría esférica ($X$, $F$, $Y$, $\nabla A$) |
| `generar_proca.py` | ecuaciones del fondo y acción cuadrática; $\theta$, $Y=4XF$, Noether radial y temporal → `resultados/proca_simbolico.pkl` |
| `chorros_proca.py` | chorros exactos de la capa de masas (casos general, R y E) |
| `verificar_proca.py` | Reissner–Nordström, correspondencia con Horndeski (2f), paridad, modos gauge |
| `caso_R.py` | vector radial: identidad $c=0$, ligaduras de la prueba, acoplos de $\delta A_0$ solo con velocidades |
| `reduccion_E.py` | vector temporal: $L_2=\kappa F^2-FW+L_{\rm rest}$, $W$, términos cinéticos |
| `energia_E.py` | eliminación de $F$ y $\tilde a$; energía de la prueba $d=N$; $\Delta=-2G_{4X}U^2N\sigma^2$ |
| `longitudinal_E.py` | $r''(0)$, $\tilde m(0)$, $\tilde c(0)$ exactos en la garganta; cancelación de $\tilde a'^2$ y $\dot{\tilde w}^2$ |
| `garganta_E.py`, `analiza_E.py` | gargantas eléctricas (barrido y refinado), condiciones y energías |

## Reproducibilidad

```
cd fase2m
python generar_proca.py      # ~1.5 min -> resultados/proca_simbolico.pkl
python verificar_proca.py    # ~2 min: RN, Horndeski, paridad, modos gauge
python caso_R.py             # vector radial
python reduccion_E.py        # vector temporal: estructura
python energia_E.py          # energia de la prueba d = N
python longitudinal_E.py     # coeficientes en la garganta
python garganta_E.py         # ~40 min: barrido de gargantas electricas -> resultados/garganta_E_barrido.json
python analiza_E.py          # ~15 min: refinado y energias -> resultados/analiza_E.json
```
