# Híbrido estrella–agujero de gusano (caso G2V) — Fase 1: marco teórico

*Documento de trabajo v0.1 · 2026-10-19 · notación LaTeX (renderiza en VS Code/GitHub/Obsidian)*

---

## 0. Estado del arte verificado y "hueco" que vamos a explotar

Las referencias recibidas se verificaron contra la API de arXiv. **Tres identificadores eran incorrectos** (apuntaban a artículos de puntos cuánticos, finanzas y grafos) y un título no correspondía; abajo, la versión corregida, ampliada con los trabajos de estabilidad que resultan decisivos.

| # | Referencia verificada | arXiv | Lo relevante para nosotros |
|---|---|---|---|
| R1 | Dzhunushaliev, Folomeev, Kleihaus, Kunz (DFKK), *A star harbouring a wormhole at its core*, JCAP **04** (2011) 031 | 1102.4454 ✔ | Agujero de gusano relleno de fluido perfecto politrópico; soluciones estáticas regulares, estables frente a perturbaciones lineales *en la región externa*. |
| R2 | DFKK, *Mixed neutron-star-plus-wormhole systems: Equilibrium configurations*, PRD **85**, 124028 (2012) | **1203.3615** (no 1201.3725; no PRD 85, 064028) | Materia neutrónica + escalar fantasma; estabilidad sólo por método energético. |
| R3 | DFKK, *…: Linear stability analysis*, PRD **87**, 104036 (2013) | 1302.5217 | Fantasma cuártico + fluido politrópico: **inestables**. |
| R4 | DFKK, *Hiding a neutron star inside a wormhole*, PRD **89**, 084018 (2014) | 1401.7093 (título real; no "Observational signatures…") | Fantasma sin masa, fluido anisótropo, gargantas dobles, propagación de luz: **todas inestables**. |
| R5 | Aringazin, DFKK, *Magnetic fields in mixed neutron-star-plus-wormhole systems*, JCAP **04** (2015) 005 | 1412.3194 | Con campo magnético: **inestables**. |
| R6 | DFKK, *Can mixed star-plus-wormhole systems mimic black holes?*, JCAP **08** (2016) 030 | 1601.04124 | Discos de acreción: diferencias observables frente a agujeros negros. |
| R7 | DFKK, *…: Rotating configurations*, PRD **107**, 044060 (2023) | **2210.04425** (no 2205.09575; no PRD 106, 024020) | Rotación rápida; gargantas de hasta ~3 km; estabilidad apenas discutida. |
| R8 | Nojiri, Odintsov, Folomeev, *Wormholes inside stars and black holes* (2024) | **2401.15868** (no 2404.02029; arXiv no lista revista) | Dos escalares; fantasmas eliminados con **restricciones vía multiplicadores de Lagrange** (generalización mimética): sector exótico **no dinámico**. |
| R9 | DFKK, *Wormhole solutions with a complex ghost scalar field and their instability*, PRD **97**, 024002 (2018) | 1710.01884 | Inestables. |
| R10 | DFKK + Radu, *Rotating wormholes in five dimensions*, PRD **88**, 124028 (2013) | 1309.2448 | El modo inestable desaparece por encima de un momento angular crítico (5D). |
| R11 | Bellinger *et al.*, *Solar evolution models with a central black hole* (2023) | 2312.06782 | Sol con agujero negro primordial central; ~$10^{-6}M_\odot$ compatible con observaciones actuales. |
| R12 | Dzhunushaliev, Folomeev, Urazalina, *Star-plus-wormhole systems with two interacting scalar fields*, IJMPD **24**, 1550097 (2015) | 1506.03897 | Máximo de densidad desplazado del centro; sistemas de garganta doble. |

**Lectura crítica.**

1. Toda la línea DFKK trata **estrellas de neutrones** ($2GM/Rc^2\sim0.2$–$0.4$) sostenidas por un fantasma **dinámico**, y sus análisis lineales dan **inestabilidad sistemática** (R3, R4, R5, R9).
2. R8 elimina el grado de libertad fantasma, pero (según su resumen) sólo establece existencia, en objetos compactos y sin estudiar la estabilidad.
3. **Nadie ha estudiado el régimen no compacto** (estrella de secuencia principal, $2GM_\odot/R_\odot c^2\simeq4.2\times10^{-6}$). Ahí aparece una separación de escalas $r_0/R_\odot\lesssim4\times10^{-6}$ que transforma el problema (§5).

---

## Convenciones

Signatura $(-,+,+,+)$; $G_{\mu\nu}=\kappa T_{\mu\nu}$ con $\kappa=8\pi G/c^4$. En las derivaciones usamos $G=c=1$ ($\kappa=8\pi$) y restauramos SI en los resultados enmarcados. $\rho$ es densidad de energía en $G=c=1$ y densidad de masa en SI (energía $=\rho c^2$). La prima $'$ es $d/dr$ salvo en las fórmulas en gauge de longitud propia, donde se indica $d/d\ell$ explícitamente o con subíndice $\ell$.

---

## 1. Métrica del espacio-tiempo

### 1.1 Elemento de línea

Estático y esféricamente simétrico (forma de Morris–Thorne, MT):

$$
ds^2=-e^{2\Phi(r)}c^2dt^2+\frac{dr^2}{1-b(r)/r}+r^2\left(d\theta^2+\sin^2\theta\,d\varphi^2\right),\qquad r\ge r_0 .
\tag{1.1}
$$

$b(r)$ es la **función de forma** y $\Phi(r)$ la **función de corrimiento al rojo** (el "potencial" gravitatorio: en el límite débil, $\Phi_N=c^2\Phi$). La coordenada $r$ es areal y es doble-valuada a través de la garganta, así que la carta global usa la **longitud propia radial**

$$
\ell(r)=\pm\int_{r_0}^{r}\frac{d\tilde r}{\sqrt{1-b(\tilde r)/\tilde r}},\qquad
ds^2=-e^{2\Phi(\ell)}c^2dt^2+d\ell^2+r^2(\ell)\,d\Omega^2,\quad \ell\in(-\ell_R^-,\ \ell_R^+).
\tag{1.2}
$$

**Arquitectura del híbrido** (lado $\ell>0$; el lado $\ell<0$ es su copia bajo simetría $\mathbb Z_2$ o, en general, otra estrella u otra región):

| Región | Rango | Contenido |
|---|---|---|
| $\mathcal M_{\rm WH}$ | $r_0\le r\le r_1$ | garganta + capa exótica (soporte de la topología) |
| $\mathcal M_\star$ | $r_0\le r\le R$ | plasma estelar (puede llegar hasta la garganta) |
| $\mathcal M_{\rm ext}$ | $r\ge R$ | Schwarzschild de masa $M$ (Birkhoff) |

### 1.2 Requisitos geométricos de la garganta

De $\left(dr/d\ell\right)^2=1-b/r$, derivando respecto de $\ell$:

$$
2\,r_\ell\,r_{\ell\ell}=\frac{b-b'r}{r^2}\,r_\ell\quad\Longrightarrow\quad r_{\ell\ell}=\frac{b-b'r}{2r^2}
\quad\Longrightarrow\quad \left.r_{\ell\ell}\right|_{\ell=0}=\frac{1-b'(r_0)}{2r_0}.
\tag{1.3}
$$

Las condiciones son:

- **(G1) Garganta:** $b(r_0)=r_0$, es decir, $r_\ell(0)=0$.
- **(G2) Ensanchamiento ("flare-out"):** $b'(r_0)<1$, es decir, $r_{\ell\ell}(0)>0$ (mínimo de área).
- **(G3) Unicidad:** $b(r)<r$ para $r>r_0$ (no hay otra garganta ni horizonte).
- **(G4) Atravesabilidad:** $\Phi$ es finita en todo punto, de modo que $e^{2\Phi}>0$ y no hay horizonte.
- **(G5) Empalme exterior:** en $r=R$ se pega a Schwarzschild (§3).

### 1.3 Función de forma: descomposición física

La componente $tt$ de las ecuaciones de campo (§2) da exactamente $b'=\frac{8\pi G}{c^2}r^2\rho_{\rm tot}$. Integrando desde la garganta:

$$
b(r)=r_0+\frac{2G}{c^2}\int_{r_0}^{r}4\pi\tilde r^2\rho_{\rm tot}\,d\tilde r
\quad\Longleftrightarrow\quad
m(r)\equiv\frac{c^2b(r)}{2G}=\underbrace{\frac{c^2r_0}{2G}}_{M_0}+m_{\rm ex}(r)+m_\star(r).
\tag{1.4}
$$

La garganta aporta **siempre** una "masa de garganta" $M_0=c^2r_0/2G$, porque $b(r_0)=r_0$ es la misma relación que define un radio de Schwarzschild. Definimos la **masa del núcleo** y el **parámetro de apantallamiento** $\epsilon$:

$$
M_{\rm core}\equiv M_0+M_{\rm ex}=\epsilon\,M_0,\qquad M=M_{\rm core}+M_\star .
\tag{1.5}
$$

> **Resultado 1 (presupuesto de masa).** Si la densidad de energía total es no negativa ($\rho_{\rm tot}\ge0$), entonces $m(r)$ es creciente, $M\ge M_0$, y
> $$\boxed{\,r_0\le\frac{2GM}{c^2}\,}\qquad(\text{Sol: } r_0\le2.95\ \text{km}).$$
> Numéricamente, $M_0=6.73\times10^{26}\,\text{kg}\,(r_0/1\,\text{m})=3.39\times10^{-4}M_\odot\,(r_0/1\,\text{m})$. Una garganta **no apantallada** de 1 m equivale a una masa puntual central de $\approx113\,M_\oplus$; una de 1 km, a $0.34\,M_\odot$. Sólo con energía negativa ($M_{\rm ex}<0$, $\epsilon<1$) se puede tener una garganta grande con un núcleo ligero.

### 1.4 Familia paramétrica para el núcleo exótico

Para la Fase 2 proponemos

$$
b(r)=r_0\left[\epsilon+(1-\epsilon)\,\chi(r)\right]+\frac{2G}{c^2}m_\star(r),\qquad \chi(r_0)=1,\ \ \chi'\le0,\ \ \chi(r\ge r_1)\to0,
\tag{1.6}
$$

con $m_\star(r_0)=0$, de modo que (G1) se cumple automáticamente. Por ejemplo, $\chi=r_0/r$ (cola de Ellis–Bronnikov) o un escalón suave de soporte compacto en $[r_0,r_1]$. Entonces

$$
\rho_{\rm ex}=\frac{c^2 r_0(1-\epsilon)}{8\pi G r^2}\,\chi'(r),\qquad
b'(r_0)=r_0(1-\epsilon)\chi'(r_0)+\frac{8\pi G}{c^2}r_0^2\rho_\star(r_0).
\tag{1.7}
$$

Casos límite:

- $\epsilon=0$: tipo Ellis–Bronnikov; núcleo sin masa, $M_{\rm ex}=-M_0$.
- $\epsilon=1$: garganta "tipo Schwarzschild"; $\rho_{\rm ex}=0$, pero con tensión exótica.
- $\epsilon>1$: energía fantasma con $\rho>0$.

### 1.5 La función de corrimiento ("potencial del agujero")

$\Phi$ **no es libre** una vez fijado el modelo de materia: la determina la ecuación $rr$ (§2). En la configuración $\mathbb Z_2$, $\Phi(\ell)$ es par:

$$
\Phi(\ell)=\Phi_0+\tfrac12\,\Phi_{\ell\ell,0}\,\ell^2+O(\ell^4),\qquad
a_{\rm prop}=c^2\,\Phi_\ell\;\Rightarrow\;a_{\rm prop}(0)=0 .
\tag{1.8}
$$

La garganta es el **"centro"** del objeto, con gravedad propia nula. Una garganta **atractiva** (máximo de presión del plasma en ella) exige $\Phi_{\ell\ell,0}>0$. La constante $\Phi_0$ se fija al final por el empalme $e^{2\Phi(R)}=1-2GM/Rc^2$. Esto da el corrimiento al rojo garganta–infinito, $1+z_0=e^{-\Phi_0}$, que entra en los observables de neutrinos y ondas gravitacionales (Fase 3).

---

## 2. Ecuaciones TOV modificadas

### 2.1 Ecuaciones de campo

Fuente total anisótropa en base ortonormal, $T^{\hat\mu}{}_{\hat\nu}=\mathrm{diag}(-\rho,\,p_r,\,p_t,\,p_t)$, que suma el plasma isótropo ($\rho_\star,P$) y el sector exótico ($\rho_{\rm ex},p_r^{\rm ex},p_t^{\rm ex}$). Para (1.1), con $G=c=1$:

$$
\begin{aligned}
G_{\hat t\hat t}&=\frac{b'}{r^2}=8\pi\rho,\\[2pt]
G_{\hat r\hat r}&=-\frac{b}{r^3}+2\left(1-\frac br\right)\frac{\Phi'}{r}=8\pi p_r,\\[2pt]
G_{\hat\theta\hat\theta}&=\left(1-\frac br\right)\left[\Phi''+\Phi'^2-\frac{b'r-b}{2r(r-b)}\Phi'-\frac{b'r-b}{2r^2(r-b)}+\frac{\Phi'}{r}\right]=8\pi p_t .
\end{aligned}
\tag{2.1}
$$

La identidad de Bianchi hace que la tercera ecuación sea equivalente a $\nabla_\mu T^{\mu}{}_{r}=0$. Con $\Gamma^t_{tr}=\Phi'$, $\Gamma^\theta_{\theta r}=\Gamma^\varphi_{\varphi r}=1/r$:

$$
\nabla_\mu T^{\mu}{}_{r}=p_r'+\Big(\Phi'+\Gamma^r_{rr}+\tfrac2r\Big)p_r-\Big(-\Phi'\rho+\Gamma^r_{rr}p_r+\tfrac2r p_t\Big)
=p_r'+(\rho+p_r)\Phi'-\tfrac2r(p_t-p_r)=0 .
\tag{2.2}
$$

### 2.2 Sistema TOV generalizado (fuente total)

Despejando $\Phi'$ de $G_{\hat r\hat r}$ con $b=2Gm/c^2$, y restaurando unidades:

$$
\boxed{
\begin{aligned}
\frac{dm}{dr}&=4\pi r^2\rho_{\rm tot},\qquad m(r_0)=M_0=\frac{c^2r_0}{2G},\\[2pt]
\frac{d\Phi}{dr}&=\frac{G\left(m+4\pi r^3p_r/c^2\right)}{c^2r^2\left(1-\dfrac{2Gm}{rc^2}\right)},\\[2pt]
\frac{dp_r}{dr}&=-\frac{G\left(\rho_{\rm tot}+p_r/c^2\right)\left(m+4\pi r^3p_r/c^2\right)}{r^2\left(1-\dfrac{2Gm}{rc^2}\right)}+\frac{2}{r}\left(p_t-p_r\right).
\end{aligned}}
\tag{2.3}
$$

Respecto del TOV estándar hay tres cambios: (a) la integración empieza en $r_0>0$ con $m(r_0)=M_0\neq0$; (b) aparece el término de anisotropía; (c) **$r_0$ es un punto singular** de la ecuación, porque $1-2Gm/rc^2\to0$.

### 2.3 Condición de regularidad en la garganta y violación de la NEC

**Proposición (tensión de garganta).** Si $\Phi$ es finita (G4), entonces

$$
\boxed{\,p_r^{\rm tot}(r_0)=-\tau_0,\qquad \tau_0=\frac{c^4}{8\pi G r_0^2}=4.82\times10^{42}\ \text{Pa}\left(\frac{1\ \text{m}}{r_0}\right)^2 .}
\tag{2.4}
$$

*Demostración.* Cerca de la garganta, $r-b\simeq(1-b'_0)(r-r_0)$. Sea $N_0\equiv r_0+8\pi r_0^3p_r(r_0)$. Si $N_0\neq0$, (2.1) da $\Phi'\simeq N_0/[2r_0^2(1-b'_0)(r-r_0)]$, y $\Phi$ diverge logarítmicamente. Si $N_0>0$, entonces $\Phi\to-\infty$ y hay un horizonte; si $N_0<0$, entonces $\Phi\to+\infty$ y la configuración es singular. Ambos casos contradicen (G4), luego $N_0=0$. $\square$

Para la escala: la presión central del Sol es $P_{c,\odot}\simeq2.5\times10^{16}$ Pa. Una garganta de 1 km exige $\tau_0/P_{c,\odot}\simeq2\times10^{20}$; ambas se igualan sólo si $r_0\simeq1.4\times10^{13}$ m $\approx93$ UA.

**Corolario (NEC).** De $G_{\hat t\hat t}$, $\rho_0=b'_0c^2/(8\pi G r_0^2)$ (densidad de energía). Por tanto

$$
\boxed{\,\left(\rho c^2+p_r\right)_{r_0}=\frac{c^4}{8\pi G r_0^2}\,\big(b'_0-1\big)<0\,}
\tag{2.5}
$$

por (G2). En RG, la NEC radial se viola **siempre** en la garganta (Morris–Thorne 1988; para gargantas genéricas, Hochberg–Visser 1997). En §4 se examina quién paga esta violación.

### 2.4 Descomposición en dos componentes: el TOV del plasma solar

Si el plasma y el sector exótico sólo interactúan gravitatoriamente, cada uno se conserva por separado. Para el plasma isótropo, $P'=-(\rho_\star c^2+P)\Phi'$. Sustituyendo $\Phi'$ de (2.3) y $m=M_0+m_{\rm ex}+m_\star$:

$$
\boxed{\;
\frac{dP}{dr}=-\,\frac{G\left(\rho_\star+\dfrac{P}{c^2}\right)\left[M_0+m_{\rm ex}(r)+m_\star(r)+\dfrac{4\pi r^3}{c^2}\big(P+p_r^{\rm ex}\big)\right]}{r^2\left[1-\dfrac{2G}{rc^2}\big(M_0+m_{\rm ex}+m_\star\big)\right]}
\;}
\tag{2.6}
$$

$$
\frac{dp_r^{\rm ex}}{dr}=-\left(\rho_{\rm ex}c^2+p_r^{\rm ex}\right)\frac{d\Phi}{dr}+\frac2r\left(p_t^{\rm ex}-p_r^{\rm ex}\right),\qquad
\frac{dm_\star}{dr}=4\pi r^2\rho_\star,\quad \frac{dm_{\rm ex}}{dr}=4\pi r^2\rho_{\rm ex}.
\tag{2.7}
$$

La ecuación (2.6) se lee término a término:

1. $M_0+m_{\rm ex}$: el núcleo actúa como masa central, que puede ser nula si $\epsilon=0$.
2. $4\pi r^3p_r^{\rm ex}/c^2<0$: la tensión exótica **antigravita**. En RG la presión pesa, y una presión negativa resta peso.
3. $1-2Gm/rc^2\to0$: es la región de campo fuerte localizada en torno a la garganta.

**Cierre del sistema (conteo).** Las incógnitas son $\{b,\Phi,P,\rho_\star,T,\rho_{\rm ex},p_r^{\rm ex},p_t^{\rm ex}\}$. Las ecuaciones son $\{tt,\ rr,\ \text{conservación total},\ \text{conservación del plasma},\ \text{EoS del plasma}\ P(\rho_\star,T,X_i),\ \text{transporte de energía}\}$. Faltan **dos relaciones de cierre del sector exótico**: una EoS anisótropa o, mejor, un lagrangiano (escalar fantasma; escalares con ligadura, R8; términos de curvatura, §4).

### 2.5 Límite post-newtoniano fuera de la capa exótica

Para $r>r_1$ se tiene $p_r^{\rm ex}=0$ y $m_{\rm ex}=M_{\rm ex}$. Desarrollando (2.6) en $1/c^2$:

$$
\frac{dP}{dr}=-\frac{G\left[M_{\rm core}+m_\star(r)\right]\rho_\star}{r^2}\left[1+\frac{P}{\rho_\star c^2}+\frac{4\pi r^3P}{(M_{\rm core}+m_\star)c^2}+\frac{2G(M_{\rm core}+m_\star)}{rc^2}+O(c^{-4})\right].
\tag{2.8}
$$

En el núcleo solar las correcciones del corchete son $\sim10^{-6}$.

> **Resultado 3 (equivalencia efectiva).** A la precisión de los modelos solares, y para $r\gg r_1$, el híbrido es una estrella newtoniana con **masa puntual central $M_{\rm core}=\epsilon M_0$**, como las "estrellas de Hawking" de R11. Hay tres diferencias cruciales: (i) no hay acreción ni luminosidad de acreción; (ii) $M_{\rm core}$ y $r_0$ son independientes, mientras que para un agujero negro $r_s=2GM/c^2$; (iii) la garganta puede transportar energía o materia entre las dos bocas (§3.3).

### 2.6 Formulación regular en longitud propia (para la Fase 2)

En el gauge (1.2), con $G=c=1$ y subíndice $\ell$ para $d/d\ell$:

$$
8\pi\rho=\frac{1-r_\ell^2}{r^2}-\frac{2r_{\ell\ell}}{r},\qquad
8\pi p_r=-\frac{1-r_\ell^2}{r^2}+\frac{2r_\ell\Phi_\ell}{r},\qquad
8\pi p_t=\frac{r_{\ell\ell}}{r}+\Phi_{\ell\ell}+\Phi_\ell^2+\frac{r_\ell\Phi_\ell}{r}.
\tag{2.9}
$$

La primera y la tercera ecuación se integran como sistema de evolución, que es **regular en $\ell=0$**. La segunda es una ligadura: se impone en $\ell=0$ y la conservación la preserva.

$$
r_{\ell\ell}=\frac{1-r_\ell^2}{2r}-4\pi r\rho_{\rm tot},\qquad
\Phi_{\ell\ell}=8\pi p_t^{\rm tot}-\frac{r_{\ell\ell}}{r}-\Phi_\ell^2-\frac{r_\ell\Phi_\ell}{r},\qquad
P_\ell=-(\rho_\star+P)\,\Phi_\ell .
\tag{2.10}
$$

> **Resultado 4 (gravedad en la garganta).** Con $\mathbb Z_2$ ($r_\ell=\Phi_\ell=0$ en $\ell=0$), combinando (2.9) con la ligadura $8\pi p_{r,0}=-1/r_0^2$:
> $$\boxed{\;\Phi_{\ell\ell,0}=4\pi\left(\rho+p_r+2p_t\right)_0\;}\qquad\Big(\text{SI: } c^2\Phi_{\ell\ell,0}=\tfrac{4\pi G}{c^2}\big(\rho c^2+p_r+2p_t\big)_0\Big).$$
> Es la **densidad activa de Tolman–Whittaker** en la garganta. La garganta es atractiva (presión del plasma máxima en ella) **si y sólo si** $(\rho+p_r+2p_t)_0>0$, aunque la NEC esté violada. Si la densidad activa es negativa, el máximo de densidad se desplaza fuera del centro, lo que explica las configuraciones de "densidad desplazada" y garganta doble de R4, R5 y R12. Cerca de la garganta, $a(\ell)\simeq c^2\Phi_{\ell\ell,0}\,\ell$. Esto es **tres veces** el valor de un centro regular con la misma densidad, porque la geometría es localmente un cilindro $S^2\times\mathbb R$: el flujo gravitatorio no se diluye y se comporta como el de una lámina plana.

**Ejemplo de cierre (DFKK, escalar fantasma sin masa).** Para un fantasma con $V=0$: $\rho_\phi=p_{r,\phi}=-\tfrac12\phi_\ell^2$ y $p_{t,\phi}=+\tfrac12\phi_\ell^2$. La ecuación $\Box\phi=0$ da $\phi_\ell=D/(r^2e^{\Phi})$, con carga escalar $D$. La ligadura en la garganta fija

$$
\frac{D^2}{2r_0^4e^{2\Phi_0}}=P_0+\frac{1}{8\pi r_0^2},\qquad
r_{\ell\ell,0}=\frac1{r_0}-4\pi r_0\left(\rho_{\star0}-P_0\right)>0,\qquad
\Phi_{\ell\ell,0}=4\pi\left(\rho_{\star0}+3P_0\right).
\tag{2.11}
$$

La densidad activa del fantasma sin masa es **exactamente nula** ($-X-X+2X=0$). En la garganta, el plasma ve sólo su propio $\rho+3P$.

---

## 3. Condiciones de contorno y de empalme

### 3.1 Garganta suave ($\mathbb Z_2$): la garganta como "centro" del objeto

Si la solución $\ell>0$ se refleja en $\ell\mapsto-\ell$, la ausencia de capa delgada (condiciones de Darmois) exige que $g_{\mu\nu}$ sea $C^1$ en $\ell=0$. Por tanto, **todas las derivadas impares se anulan**.

| | Centro regular ($r=0$) | Garganta $\mathbb Z_2$ ($\ell=0$) |
|---|---|---|
| Área | $0$ | $4\pi r_0^2$ |
| Masa | $m(0)=0$ | $m(r_0)=c^2r_0/2G$ |
| Geometría | $r_\ell=1$ | $r_\ell=0,\ \ r_{\ell\ell}=\frac{1-b'_0}{2r_0}>0$ |
| Potencial | $\Phi'(0)=0$ | $\Phi_\ell(0)=0,\ \ \Phi_{\ell\ell}(0)=4\pi(\rho+p_r+2p_t)_0$ |
| Esfuerzo radial | $p_r=P_c$ (libre) | $p_r^{\rm tot}=-c^4/8\pi Gr_0^2$ (fijo) |
| Plasma | $P(0)=P_c,\ P'(0)=0$ | $P(0)=P_0$ (parámetro de disparo), $P_\ell(0)=0$ |
| Energía | $L(0)=0$ | $L(0)=0$ (sin flujo neto por simetría), $T_\ell(0)=0$ |
| NEC | $\rho+p_r>0$ | $\rho+p_r<0$ |

En coordenada $r$, la regla de L'Hôpital da $\Phi'(r_0)=\Phi_{\ell\ell,0}/r_{\ell\ell,0}$, que es finito. Además, $r-r_0=\tfrac12r_{\ell\ell,0}\ell^2+O(\ell^4)$, de modo que una integración en $r$ debe arrancar en $r_0+\delta$ con valores tomados de la serie. Si $r_\ell(0^+)\neq0$ o $\Phi_\ell(0^+)\neq0$, la reflexión produce una **esquina**, y con ella una capa delgada (§3.2).

### 3.2 Garganta de capa delgada (cirugía de Visser) y su estabilidad

Se toman dos copias de la región $\{r\ge a\}$ del interior estelar, con $a>b(a)$, y se pegan en $r=a(\tau)$. Las ecuaciones de Lanczos son $S^i{}_j=-\frac1{8\pi}\left([K^i{}_j]-\delta^i{}_j[K]\right)$ con $[K]=2K^+$. Calculando la curvatura extrínseca con normal $n^r=\sqrt{1-b/r}$:

$$
K^{\theta}{}_{\theta}{}^{\pm}=\pm\frac1a\sqrt{1-\frac{b(a)}{a}+\dot a^2},\qquad
K^{\tau}{}_{\tau}{}^{\pm}\big|_{\rm est}=\pm\,\Phi'(a)\sqrt{1-\frac{b(a)}{a}} .
\tag{3.1}
$$

De ahí, con $G=c=1$:

$$
\boxed{\;\sigma=-\frac{1}{2\pi a}\sqrt{1-\frac{b(a)}{a}+\dot a^2}<0,\qquad
\mathcal P_{\rm est}=\frac{1}{4\pi}\sqrt{1-\frac{b(a)}{a}}\left(\Phi'(a)+\frac1a\right).\;}
\tag{3.2}
$$

*Comprobación.* Con $b=2M$ y $e^{2\Phi}=1-2M/r$ se recupera la capa de Schwarzschild de Visser (1989): $\mathcal P=\frac{1-M/a}{4\pi a\sqrt{1-2M/a}}$.

Hay que destacar que aquí **$M_{\rm core}=m(a^+)$ es libre**, con $0\le2Gm(a)/c^2<a$. La capa, de energía propia negativa ($\approx-2c^2a/G$ en equivalente de masa si $b\ll a$), permite cualquier $\epsilon\in[0,1)$, lo que es consistente con el Resultado 1.

**Dinámica.** La componente $\theta\theta$ da una ecuación de energía:

$$
\dot a^2+V(a)=0,\qquad V(a)=1-\frac{b(a)}{a}-\big[2\pi a\,\sigma(a)\big]^2 .
\tag{3.3}
$$

$\sigma(a)$ sale de la ley de conservación de la capa. Esta incluye **flujo desde el volumen**, porque el plasma toca la capa: $T_{\mu\nu}u^\mu n^\nu=\dot a\,\frac{\sqrt{1-b/a+\dot a^2}}{1-b/a}(\rho+p_r)$ en cada cara. Resulta

$$
\boxed{\;\frac{d\sigma}{da}=-\frac2a\left(\sigma+\mathcal P\right)+\Xi,\qquad
\Xi=\frac{2\sqrt{1-b/a+\dot a^2}}{1-b/a}\,\big(\rho+p_r\big)_{\rm vol}\Big|_{r=a}\;}
\tag{3.4}
$$

(comprobado contra el límite estático de (3.2) en §6). Con una EoS de capa $\mathcal P(\sigma)$ y $\beta^2\equiv d\mathcal P/d\sigma|_{a_0}$:

$$
\text{estable}\iff V(a_0)=V'(a_0)=0\ \ \text{y}\ \ V''(a_0)=-\Big(\frac ba\Big)''-8\pi^2\Big[(u')^2+u\,u''\Big]_{a_0}>0,\qquad u\equiv a\sigma,
\tag{3.5}
$$

donde $\sigma''=-\frac2a(1+\beta^2)\sigma'+\frac2{a^2}(\sigma+\mathcal P)+\Xi'$. La condición es **lineal en $\beta^2$** y da una cota $\beta^2\gtrless\beta^2_{\rm crit}(a_0;\,b,\Phi)$, que la Fase 2 calculará con el perfil solar. En el caso Schwarzschild (Poisson–Visser 1995), la región estable exige $\beta^2$ fuera de $[0,1]$, lo cual es ya una señal de EoS exótica.

*Limitación:* este es un análisis de "capa de prueba en un volumen congelado". El problema completo acopla la capa a los modos radiales del plasma.

### 3.3 Fotosfera ($r=R$) y empalme con Schwarzschild

Darmois–Israel sin capa en $r=R$, con exterior $ds^2=-(1-\frac{2GM}{rc^2})c^2dt^2+\dots$:

$$
\underbrace{e^{2\Phi(R)}=1-\frac{2GM}{Rc^2}}_{[h_{ij}]=0},\qquad
\underbrace{m(R)=M}_{[K^\theta{}_\theta]=0},\qquad
\underbrace{p_r^{\rm tot}(R)=0\ \Rightarrow\ P(R)=0,\ p_r^{\rm ex}(R)=0}_{[K^\tau{}_\tau]=0\ \Leftrightarrow\ [\Phi']=0}.
\tag{3.6}
$$

En RG, $\rho$ puede saltar en $R$; en $f(R)$ no (§4.3). Para el Sol real se sustituye $P(R)=0$ por la atmósfera gris de Eddington en $\tau=2/3$:

$$
P_{\rm ph}\simeq\frac23\frac{g_s}{\bar\kappa},\qquad g_s=\frac{GM}{R^2\sqrt{1-2GM/Rc^2}},\qquad T(R)=T_{\rm eff},\qquad L_\infty=4\pi R^2\sigma_{\rm SB}T_{\rm eff}^4\left(1-\frac{2GM}{Rc^2}\right).
\tag{3.7}
$$

Por el **teorema de Birkhoff**, el exterior es idéntico al de un Sol ordinario de igual $M$. Los tests clásicos del sistema solar son ciegos a la garganta; sólo la **sismología**, los **neutrinos** y la **dinámica** (ondas gravitacionales) la ven (Fase 3).

### 3.4 Ecuaciones de energía y condición térmica en la garganta (anticipo de la Fase 2)

En la generalización relativista (Thorne 1977), con luminosidad desplazada al infinito $L_\infty=Le^{2\Phi}$:

$$
\frac{dL_\infty}{dr}=\frac{4\pi r^2\rho_\star\,\varepsilon_{\rm nuc}\,e^{2\Phi}}{\sqrt{1-2Gm/rc^2}},\qquad
\frac{d\left(Te^{\Phi}\right)}{dr}=-\frac{3\bar\kappa\rho_\star}{64\pi\sigma_{\rm SB}T^3r^2}\,\frac{L_\infty e^{-\Phi}}{\sqrt{1-2Gm/rc^2}},
\tag{3.8}
$$

más el criterio de Schwarzschild/Ledoux para la convección. El equilibrio de Tolman ($Te^\Phi=$ cte. sin flujo) hace que la región de garganta esté **más caliente por un factor $e^{-\Phi_0}$**.

- **Caso $\mathbb Z_2$:** $L_\infty(r_0)=0$.
- **Caso asimétrico:** la continuidad del flujo a través de la garganta da $L_\infty(0^+)=-L_\infty(0^-)\equiv L_{\rm th}$. La garganta es un **conducto térmico entre las dos bocas**. En equilibrio, $L_\gamma=L_{\rm nuc}+L_{\rm th}$, mientras que los neutrinos sólo miden $L_{\rm nuc}$. La comparación entre la luminosidad solar inferida de neutrinos y la fotónica (Bergström *et al.* 2016, JHEP 03, 132; consistencia de orden 10 %, cifra exacta a extraer en la Fase 3) acota directamente $|L_{\rm th}|/L_\odot$.

---

## 4. Materia exótica frente a gravedad modificada

### 4.1 Qué hay que pagar, sea cual sea el mecanismo

Toda teoría que pueda escribirse como $G_{\mu\nu}=\kappa\,T^{\rm eff}_{\mu\nu}$ hereda la geometría de §1–2. Para el vector nulo radial $k^{\hat\mu}=(1,1,0,0)$, (2.9) da en la garganta

$$
G_{\hat\mu\hat\nu}k^{\hat\mu}k^{\hat\nu}\big|_0=-\frac{2r_{\ell\ell,0}}{r_0}=-\frac{1-b'_0}{r_0^2}<0
\quad\Longrightarrow\quad T^{\rm eff}_{\mu\nu}k^\mu k^\nu\big|_0<0 .
\tag{4.1}
$$

Además, la censura topológica (Friedman–Schleich–Witt 1993) exige violar la ANEC a lo largo de las geodésicas que atraviesan la garganta. La única pregunta es **quién** aporta la contribución negativa.

> **Resultado 5 (escala de nueva física).** El mecanismo debe generar un esfuerzo efectivo $\sim\tau_0=c^4/8\pi Gr_0^2$ en la garganta. En gravedad modificada esto significa correcciones de **orden uno** a curvaturas $\mathcal R\sim r_0^{-2}$. Por tanto, la longitud característica de la teoría debe cumplir $\ell_{\rm MG}\gtrsim r_0$. Las cotas experimentales sobre $\ell_{\rm MG}$ se convierten así en **cotas superiores sobre $r_0$**.

### 4.2 Vía A: Relatividad General con sector exótico

**(a) Fluido fantasma anisótropo.** $p_r=w\rho c^2$ con $w<-1$ (Lobo 2005; Sushkov 2005). Es fenomenológico y no tiene microfísica.

**(b) Escalar fantasma dinámico** (DFKK, R1–R9). $\mathcal L_\phi=+\tfrac12(\partial\phi)^2-V$, con

$$
T^{(\phi)}_{\mu\nu}=-\partial_\mu\phi\,\partial_\nu\phi+g_{\mu\nu}\Big[\tfrac12(\partial\phi)^2-V\Big],\qquad
\rho_\phi=-X+V,\quad p_{r,\phi}=-X-V,\quad p_{t,\phi}=X-V,\quad X=\tfrac12g^{rr}\phi'^2 ,
\tag{4.2}
$$

de modo que $\rho_\phi+p_{r,\phi}=-2X<0$. Tiene dos problemas: (i) el vacío cuántico es inestable, porque el hamiltoniano no está acotado inferiormente; (ii) **inestabilidad lineal clásica**, tanto en Ellis–Bronnikov (Gonzalez–Guzmán–Sarbach 2009) como en todos los sistemas mixtos con estrella de neutrones (R3, R4, R5, R9).

**(c) Semiclásico** (Casimir, $\langle T_{\mu\nu}\rangle_{\rm ren}$). Las desigualdades cuánticas (Ford–Roman 1996) confinan la energía negativa en bandas de espesor $\ll r_0$, lo que es inviable para gargantas macroscópicas.

**(d) Capa delgada** (§3.2). Concentra la violación en $\sigma<0$ y exige una EoS con $\beta^2\notin[0,1]$ en el caso Schwarzschild.

**(e) Campos no fantasma en RG.** Espinores de Dirac con campo de Maxwell (Blázquez-Salcedo–Knoll–Radu 2021; Dzhunushaliev–Folomeev 2025, arXiv:2510.03656, 2511.02431) dan soluciones asimétricas sin fantasma. Su escalado a tamaños astrofísicos está abierto.

**(f) Fantasma con ligadura** (R8, Nojiri–Odintsov–Folomeev 2024). Dos escalares con restricciones impuestas por multiplicadores de Lagrange, que generalizan la restricción mimética: el sector exótico existe en el fondo pero **no propaga grados de libertad**. Su estabilidad clásica **no está analizada**, y esa es nuestra oportunidad (§5).

### 4.3 Vía B: gravedad modificada (el plasma cumple la NEC y la geometría paga)

**(i) $f(R)$ métrica.** Con $F\equiv f_R$:

$$
G_{\mu\nu}=\frac1F\Big[\kappa T_{\mu\nu}+\nabla_\mu\nabla_\nu F-g_{\mu\nu}\Box F-\tfrac12g_{\mu\nu}(RF-f)\Big]
\;\Longrightarrow\;
F\,G_{\mu\nu}k^\mu k^\nu=\kappa T_{\mu\nu}k^\mu k^\nu+k^\mu k^\nu\nabla_\mu\nabla_\nu F .
\tag{4.3}
$$

En el gauge $\ell$, $\nabla_{\hat t}\nabla_{\hat t}F=-\Phi_\ell F_\ell$ y $\nabla_{\hat\ell}\nabla_{\hat\ell}F=F_{\ell\ell}$. En la garganta $\mathbb Z_2$:

$$
-\frac{1-b'_0}{r_0^2}\,F_0=8\pi\left(\rho+p_r\right)_{m,0}+F_{\ell\ell,0}.
\tag{4.4}
$$

> **Resultado 6.** Si el plasma cumple la NEC en la garganta y $F_0>0$ (graviton no fantasma), entonces $\;F_{\ell\ell,0}\le-\dfrac{1-b'_0}{r_0^2}F_0<0$: $f_R$ debe tener un **máximo pronunciado** en la garganta.
>
> **Estimación de orden de magnitud** para $f=R+\alpha R^2$: $F=1+2\alpha R$, con $R\sim r_0^{-2}$ y $F_{\ell\ell}\sim2\alpha/r_0^4$, lo que da $r_0\lesssim\sqrt{2|\alpha|}$. Las cotas de laboratorio sobre el coeficiente de $R^2$ son de orden $10^{-9}$ m² (p. ej. Näf–Jetzer 2010; cifra a precisar). Resulta **$r_0\lesssim10^{-4}$ m**: sólo micro-gargantas.
>
> Los modelos cosmológicos (tipo Hu–Sawicki) se reducen a RG a alta curvatura y son inútiles aquí. Además, el apantallamiento camaleónico en el interior denso suprime el escalarón.

Las condiciones de empalme en $f(R)$ son más estrictas (Deruelle–Sasaki–Sendouda 2008; Senovilla 2013): $R$ debe ser continuo, lo que prohíbe saltos de $\rho$ en $R_\odot$, y las capas delgadas llevan términos de doble capa. Eso modificaría §3.2 y §3.3.

**(ii) Escalar–Gauss–Bonnet (4D).** La acción es

$$
S=\frac{c^4}{16\pi G}\int d^4x\sqrt{-g}\Big[R-\tfrac12(\partial\phi)^2+\alpha\,\xi(\phi)\,\mathcal G\Big]+S_m,\qquad
\mathcal G=R^2-4R_{\mu\nu}R^{\mu\nu}+R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}.
\tag{4.5}
$$

El término de Gauss–Bonnet contribuye a $G_{\mu\nu}$ mediante el doble dual del Riemann contraído con $\nabla\nabla\xi$. En la garganta $\mathbb Z_2$ se anulan $\nabla_t\nabla_t\xi\propto\Phi_\ell$ y $\nabla_\theta\nabla_\theta\xi\propto r_\ell$. Sólo sobrevive $\xi_{\ell\ell}$, multiplicado por la curvatura seccional de la 2-esfera, $R_{\hat\theta\hat\varphi\hat\theta\hat\varphi}=1/r_0^2$. Así, la contribución a la NEC es $\sim\alpha\,\xi_{\ell\ell,0}/r_0^2$, y violarla exige $|\alpha|\,\Delta\xi\gtrsim r_0^2$.

Kanti–Kleihaus–Kunz (2011, 2012) construyeron gargantas sin materia exótica, con un subconjunto que declararon linealmente estable frente a perturbaciones radiales. *Corrección posterior:* Cuyubamba–Konoplya–Zhidenko (2018, arXiv:1804.11170) demostraron que son inestables para todos los parámetros, y la Fase 2b lo confirma con la garganta libre y la cáscara de la garganta dinámica ([fase2b/fase2b_resultados.md](../fase2b/fase2b_resultados.md)). Las cotas actuales de ondas gravitacionales y estrellas de neutrones sitúan $\sqrt{|\alpha|}$ en $\lesssim\mathcal O(\text{km})$ (depende del acoplamiento; a precisar). **Gargantas de hasta ~km quedan en el borde de lo permitido**, justo la escala máxima admitida por el Resultado 1. Es una coincidencia notable.

**(iii) Acoplamiento no mínimo curvatura–materia** (Bertolami–Böhmer–Harko–Lobo 2007). Con $\mathcal L=\tfrac12f_1(R)+[1+\lambda f_2(R)]\mathcal L_m$, la materia no se conserva:

$$
\nabla^\mu T_{\mu\nu}=\frac{\lambda F_2}{1+\lambda f_2}\left(g_{\mu\nu}\mathcal L_m-T_{\mu\nu}\right)\nabla^\mu R
\;\;\xrightarrow{\ \nu=r,\ \mathcal L_m=-\rho c^2\ }\;\;
\boxed{\,\frac{dP}{dr}=-\left(\rho_\star c^2+P\right)\frac{d}{dr}\Big[\Phi+\ln\big(1+\lambda f_2(R)\big)\Big].}
\tag{4.6}
$$

> **Resultado 7.** El plasma siente un **potencial efectivo** $\Phi_{\rm eff}=\Phi+\ln(1+\lambda f_2(R))$. Es una modificación *directa* del TOV, relevante sólo donde $R$ varía mucho (la garganta) y despreciable en la envolvente solar, donde $R\simeq\kappa(\rho c^2-3P)$. *Advertencia:* depende de la elección de $\mathcal L_m$; con $\mathcal L_m=P$ la fuerza extra desaparece. Gargantas sostenidas así: Garcia–Lobo 2010.

**(iv) Escalar–tensor sin fantasma.** No hay agujeros de gusano en el marco de Einstein sin fantasma (Bronnikov–Starobinsky 2007), salvo continuaciones conformes.

### 4.4 Balance comparativo

| Mecanismo | Quién viola la NEC | Estabilidad | Restricción observacional | Viabilidad G2V |
|---|---|---|---|---|
| Fantasma dinámico (DFKK) | escalar $\rho+p_r=-2X$ | **inestable** (R3–R5, R9) | Birkhoff: invisible desde fuera | Baja; sirve como *benchmark* de validación |
| Fantasma con ligadura (R8) | sector no dinámico | **no estudiada** | ídem | **Alta: candidato principal** |
| Capa delgada | $\sigma<0$ | condicional en $\beta^2$ | ídem | Media; modelo analítico de juguete |
| $f(R)$ | $F_{\ell\ell}<0$ | según el modelo | laboratorio ⇒ $r_0\lesssim10^{-4}$ m | Baja para gargantas macroscópicas |
| Escalar–Gauss–Bonnet | $\alpha\,\xi_{\ell\ell}/r_0^2$ | inestable (CKZ 2018; Fase 2b), contra lo afirmado por KKK 2012 | ondas gravitacionales ⇒ $\sqrt\alpha\lesssim$ km | Media-alta para $r_0\lesssim$ km |
| No mínimo $f_2(R)\mathcal L_m$ | acoplamiento | abierta | ambigüedad de $\mathcal L_m$ | Media; modifica el TOV directamente |
| Dirac + Maxwell (RG) | campos no fantasma | abierta | — | Exploratoria |

---

## 5. Propuesta: un punto de vista nuevo (síntesis de R1–R11 y de los Resultados 1–7)

La literatura (R1–R10) ataca el problema **globalmente y en régimen compacto**: integra el sistema completo desde la garganta hasta la superficie de una estrella de neutrones con un fantasma dinámico, y encuentra inestabilidad. Proponemos invertir la estrategia para una estrella G2V, apoyándonos en tres ideas.

### 5.1 Idea 1: jerarquía de **tres** escalas y la garganta como condición de contorno interna

Para el Sol, el Resultado 1 impone $r_0\le2GM_\odot/c^2$, y por tanto $r_0/R_\odot\le4.2\times10^{-6}$. Entre ambos extremos aparece una tercera escala que ningún trabajo previo usa: el **radio de Bondi hidrostático** del núcleo, donde la energía gravitatoria de $M_{\rm core}$ iguala la térmica del plasma central:

$$
r_B\equiv\frac{GM_{\rm core}\,\mu m_u}{k_BT_c}
=\frac{\epsilon\,r_0}{2}\,\frac{\mu m_uc^2}{k_BT_c}
\simeq2.9\times10^{5}\,\epsilon\,r_0
\qquad(\mu\simeq0.85,\ T_c\simeq1.57\times10^7\,\text{K}).
\tag{5.1}
$$

La jerarquía es $r_0\ll r_B\ll R_\odot$, salvo que $\epsilon\lesssim3\times10^{-6}$. Una garganta no apantallada de 1 m tiene $r_B\approx290$ km. Esto sugiere un **desarrollo asintótico acoplado** en lugar de una integración global con nueve órdenes de magnitud de rigidez:

| Zona | Rango | Física | Herramienta |
|---|---|---|---|
| Interior relativista | $r_0\le r\lesssim$ pocos $r_0$ | ecuaciones completas (2.10) + sector exótico; plasma ≈ uniforme | integrador propio en $\ell$ (Python) |
| Cúspide de Bondi | $r_0\ll r\lesssim r_B$ | atmósfera hidrostática alrededor de $M_{\rm core}$; $k_BT_{\rm vir}\sim\mu m_uGM_{\rm core}/r$ | post-newtoniano (2.8) + transporte (3.8) |
| Estrella | $r_B\ll r\le R_\odot$ | modelo solar estándar con masa puntual | código solar (p. ej. MESA, con condiciones de contorno centrales $M_c,R_c,L_c$) |

La garganta entra en el modelo solar **sólo a través de unas pocas "cargas"**: $M_{\rm core}=\epsilon M_0$, $L_{\rm th}$, un flujo de masa $\dot M_{\rm th}$ (nulo en el caso estático) y el salto de presión a través de la zona relativista. Esto hace el problema **computacionalmente factible** con herramientas estándar.

**Consecuencia física (hipótesis a verificar en la Fase 2).** Con $\epsilon\sim1$, la temperatura virial cerca de la garganta llega a $k_BT\sim\epsilon\,\mu m_uc^2/5$, en el rango de los cientos de MeV. Eso es incompatible con el Sol. Por tanto, **un híbrido solar viable exige gargantas fuertemente apantalladas ($\epsilon\ll1$) o diminutas**: la viabilidad no la decide el tamaño de la garganta, sino su masa efectiva.

### 5.2 Idea 2: el híbrido como "estrella de Hawking sin acreción"

El Resultado 3 permite **reutilizar el programa de las estrellas de Hawking** (R11), que con masas centrales de ~$10^{-6}M_\odot$ sigue siendo compatible con las observaciones, para un problema distinto. Hay tres discriminantes entre un agujero negro y un agujero de gusano en el centro:

| | Agujero negro central (R11) | Garganta (este trabajo) |
|---|---|---|
| Relación tamaño–masa | $r_s=2GM/c^2$ | $r_0$ y $M_{\rm core}=\epsilon M_0$ independientes |
| Evolución | acreción de Bondi, crecimiento, consumo de la estrella | estática ($\mathbb Z_2$); sin sumidero |
| Energía | luminosidad de acreción | $L_{\rm th}$ de cualquier signo (conducto entre bocas) |
| Centro | singularidad tras un horizonte | "centro" regular con $a=0$ y presión $P_0$ finita (tabla §3.1) |

La **cúspide hidrostática estática** (calentamiento local sin acreción) es la firma distintiva. Por la fuerte dependencia en $T$ de las tasas CNO y ⁸B, debería aparecer antes en los **neutrinos** que en la heliosismología.

### 5.3 Idea 3: la clave de la estabilidad está en combinar R3–R5/R9 con R8

- **Hecho 1:** con un fantasma *dinámico*, todos los sistemas mixtos son linealmente inestables (R3, R4, R5, R9).
- **Hecho 2:** R8 hace el fantasma *no dinámico* mediante ligaduras, pero no estudia la estabilidad.

> **Hipótesis H1.** El modo inestable de DFKK lo transporta el grado de libertad dinámico del fantasma. Con el sector exótico ligado, en simetría esférica el único grado de libertad radial que propaga es el plasma. La estabilidad se reduce entonces a un **problema de Sturm–Liouville para el desplazamiento radial del plasma** sobre el fondo con garganta, con las condiciones de §3.1: un análogo del criterio variacional de Chandrasekhar con un término de garganta. El Resultado 4 sugiere que, si $(\rho+p_r+2p_t)_0>0$, la garganta actúa como un centro atractivo normal y el modo fundamental podría ser estable con $\Gamma_1\simeq5/3$.

**Riesgos conocidos, que convierten la hipótesis en una pregunta con respuesta.** Los multiplicadores de Lagrange de tipo mimético pueden propagar perturbaciones de tipo polvo, con cáusticas e inestabilidades de gradiente, como se discute en la literatura de gravedad mimética. Hay que comprobar que el sector de ligadura no reintroduce la inestabilidad.

**La vía rotacional queda descartada para el Sol.** R10 y R7 sugieren que la rotación puede estabilizar, pero $\Omega_\odot r_0/c\sim10^{-11}$ para $r_0\sim$ km. Sigue siendo relevante sólo para estrellas de neutrones.

### 5.4 Firmas para la Fase 3 que se derivan de este marco

1. **Neutrinos:** el exceso de ⁸B/CNO de la cúspide acota $M_{\rm core}$; el cociente $L_\nu$/$L_\gamma$ acota $L_{\rm th}$ (§3.4).
2. **Modos g y separaciones pequeñas $\delta\nu_{02}$:** sensibles al gradiente de la velocidad del sonido y a la frecuencia de Brunt–Väisälä en el núcleo, que la masa puntual modifica.
3. **Ondas gravitacionales:** la región $2m/r\to1$ en torno a $r_0$ es la única con potencial efectivo fuerte. Sus modos tienen $f\sim c/2\pi r_0\approx5\times10^4$ Hz $(1\,\text{km}/r_0)$, muy por encima de la banda de LIGO. Los ecos "a lo agujero negro" no son esperables para compacidad $4\times10^{-6}$. Esto es una predicción **negativa** útil.

### 5.5 Programa concreto para la Fase 2

1. Integrador interior en $\ell$ (2.10), **validado** reproduciendo DFKK R2/R3 (fantasma + politropo de estrella de neutrones).
2. Cierres del sector exótico: (A) fantasma sin masa (*benchmark*), (B) ligadura tipo R8, (C) capa delgada (3.2)–(3.5), (D) escalar–Gauss–Bonnet.
3. Zona de cúspide más la estrella: modelo solar con masa puntual $M_{\rm core}$ y $L_{\rm th}$, con código propio o MESA.
4. Estabilidad radial: problema de Sturm–Liouville con condiciones de garganta (test de H1).
5. Barrido en el plano $(r_0,\epsilon,L_{\rm th})$ dentro de la región permitida ($r_0\le2.95$ km; §4: $r_0\lesssim10^{-4}$ m en $f(R)$, $\lesssim$ km en escalar–Gauss–Bonnet).

---

## 6. Verificación y limitaciones

**Verificación simbólica** ([verificacion_simbolica.py](verificacion_simbolica.py), SymPy). El script calcula el tensor de Einstein desde la métrica y comprueba:

- (2.1), (2.2), (2.3) y (2.5);
- el gauge $\ell$ (2.9), la ligadura de garganta y el **Resultado 4**;
- el modelo fantasma (2.11);
- la contracción $f(R)$ (4.3);
- la ley de conservación con flujo de la capa delgada (3.4), en el límite estático, y el límite de Schwarzschild de (3.2).

También imprime todas las escalas numéricas usadas en el texto. Ejecutar: `python fase1/verificacion_simbolica.py`.

**Limitaciones y preguntas abiertas.**

1. La hipótesis H1 (estabilidad con sector exótico ligado) no está demostrada: es el objetivo central de la Fase 2.
2. La estructura térmica de la cúspide de Bondi (§5.1) es una estimación de escala; requiere el transporte completo (3.8).
3. Las cotas de §4.3 son de orden de magnitud. Hay que fijar los valores exactos del coeficiente de $R^2$ y de $\sqrt\alpha$.
4. Todo es estático y esférico; el caso asimétrico ($L_{\rm th}\neq0$) sólo está planteado.
5. Falta extraer de Bergström *et al.* (2016) la cifra exacta de $L_\nu/L_\gamma$.

---

## Referencias

*Verificadas en esta sesión vía la API de arXiv:* R1–R12 (§0); Bergström *et al.*, JHEP **03** (2016) 132 [arXiv:1601.00972]; Dzhunushaliev–Folomeev, arXiv:2510.03656 y 2511.02431. *Citadas de memoria* (referencias estándar; **comprobar volumen y página antes de publicar**):

- M. S. Morris, K. S. Thorne, Am. J. Phys. **56**, 395 (1988).
- M. Visser, PRD **39**, 3182 (1989); *Lorentzian Wormholes* (AIP, 1995).
- E. Poisson, M. Visser, PRD **52**, 7318 (1995).
- W. Israel, Nuovo Cim. B **44**, 1 (1966).
- D. Hochberg, M. Visser, PRD **56**, 4745 (1997).
- J. L. Friedman, K. Schleich, D. M. Witt, PRL **71**, 1486 (1993).
- K. S. Thorne, ApJ **212**, 825 (1977).
- H. G. Ellis, J. Math. Phys. **14**, 104 (1973).
- K. A. Bronnikov, Acta Phys. Pol. B **4**, 251 (1973).
- F. S. N. Lobo, PRD **71**, 084011 (2005).
- S. V. Sushkov, PRD **71**, 043520 (2005).
- J. A. González, F. S. Guzmán, O. Sarbach, CQG **26**, 015010 (2009).
- L. H. Ford, T. A. Roman, PRD **53**, 5496 (1996).
- J. L. Blázquez-Salcedo, C. Knoll, E. Radu, PRL **126**, 101102 (2021).
- P. Kanti, B. Kleihaus, J. Kunz, PRL **107**, 271101 (2011); PRD **85**, 044007 (2012).
- N. Deruelle, M. Sasaki, Y. Sendouda, PTP **119**, 237 (2008).
- J. M. M. Senovilla, PRD **88**, 064015 (2013).
- C. Näf, P. Jetzer, PRD **81**, 104003 (2010).
- O. Bertolami, C. G. Böhmer, T. Harko, F. S. N. Lobo, PRD **75**, 104016 (2007).
- N. M. Garcia, F. S. N. Lobo, PRD **82**, 104018 (2010).
- K. A. Bronnikov, A. A. Starobinsky, JETP Lett. **85**, 1 (2007).
