# Cómo publicar el artículo

## Archivos

| archivo | uso |
|---|---|
| `wormhole_monopole.tex` / `.pdf` | versión en inglés: la que se envía |
| `arxiv/wormhole_monopole_arxiv.tar.gz` | paquete para subir a arXiv (solo el `.tex`; la bibliografía va incluida) |
| `gargantas_monopolo.tex` / `.pdf` | versión original en español (para el repositorio) |

Los dos PDF se compilaron con Tectonic 0.17 y salen sin errores, sin referencias rotas y sin desbordes de margen (11 páginas). Para recompilar:

```bash
tectonic -X compile wormhole_monopole.tex
```

## 1. arXiv (preprint)

1. Crea una cuenta en https://arxiv.org/user/register. Usa tu nombre completo real, porque es el que figurará como autor.
2. **Endorsement.** Es tu primer envío a `gr-qc`, así que arXiv te pedirá el aval de alguien que ya publique en esa categoría. Al iniciar el envío, arXiv te da un código de endorsement; envíaselo a un investigador del área. Por ejemplo, alguno de los autores citados que trabaje en estabilidad de gargantas o en Horndeski: con el artículo adjunto y una nota breve.
3. Inicia un envío nuevo (*Start new submission*):
   - Licencia: **CC BY 4.0** (recomendada).
   - Categoría principal: **gr-qc**. Cross-list: **hep-th**.
   - Sube `arxiv/wormhole_monopole_arxiv.tar.gz`. arXiv lo compila con pdflatex; revisa el PDF que genera antes de confirmar.
4. Metadatos:
   - **Title:** The $d=N$ deficit and the throat ghost: a common mechanism for the monopole instability of wormholes in second-order scalar- and vector-tensor theories
   - **Authors:** tu nombre completo
   - **Comments:** 11 pages, 1 table. Code and data at https://github.com/iivanlopez/gusanos
   - **Abstract:** copia el del `.tex`, quitando los `\cite{...}`. arXiv acepta `$...$` en el resumen.
5. Revisa la vista previa y envía. El artículo aparece el siguiente día hábil, tras la moderación.

## 2. Revista con revisión por pares

Cuando esté en arXiv, envíalo a una de estas revistas. Todas aceptan el envío directamente desde el identificador de arXiv:
- **Physical Review D** (APS). Es la más natural por la bibliografía, que cita sobre todo PRD y JCAP. Se envía en https://authors.aps.org/Submissions; puede pedir pasar el texto a la plantilla REVTeX.
- **Classical and Quantum Gravity** (IOP).
- **JCAP**, si se enfoca hacia teorías modificadas de gravedad o cosmología.

## Antes de enviar: responsabilidad del autor

- Las políticas de arXiv y de las revistas exigen que **el autor humano responda de todo el contenido**. Una IA no puede figurar como autora. El uso de Claude ya se declara en la nota al pie del autor; algunas revistas piden además declararlo en la carta de presentación.
- Lee el artículo completo y comprueba que entiendes y suscribes cada afirmación. En especial:
  - el Teorema 1 (déficit $d=N$);
  - el Teorema 3 (regla de signos);
  - la Conjetura 1 y los números de las tablas.

  Todos los cálculos se pueden reproducir con los scripts del repositorio (apéndice A).
- Si tienes una afiliación o un ORCID, añádelos (https://orcid.org).
