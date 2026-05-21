# Guía del Manuscrito — División de Trabajo por Integrante
## Análisis Integrativo del Transportador SLC1A2 (EAAT2) en DEE41

**Límite:** 5–10 páginas | Fuente: 14 pt (títulos), 12 pt (texto), 11 pt (pies de figura)  
**Formato:** Manuscrito científico PDF

---

## Título propuesto

> **Análisis computacional de las mutaciones G82R y L85P en el dominio de trimerización del transportador EAAT2 (SLC1A2) y su impacto sobre la excitabilidad neuronal en la Encefalopatía Epiléptica del Desarrollo (DEE41): un abordaje bioinformático, estructural y de simulación Hodgkin-Huxley**

---

## Distribución por Integrante (9 personas → ~1 sección/subsección cada una)

---

### 🔵 Integrante 1 — Abstract + Título final
**Extensión sugerida:** 200–250 palabras (1 párrafo exacto)

**Debe incluir los 4 elementos obligatorios:**
1. Problema biológico: falla de EAAT2 en DEE41 → excitotoxicidad
2. Herramienta: bioinformática (Python/Biopython) + modelado estructural (PyMOL/FoldX) + simulación HH
3. Hallazgo principal: las mutaciones G82R/L85P aumentan ΔΔG (desestabilización) y en la simulación HH se observa X potenciales de acción más en el mutante
4. Implicación clínica: hiperexcitabilidad neuronal → crisis epilépticas DEE41

**Tip:** Escribir este apartado al final, cuando ya todos los resultados estén listos.

---

### 🔵 Integrante 2 — Introducción: Biofísica del EAAT2
**Extensión sugerida:** ~1 página (400–600 palabras)

**Contenido mínimo:**
- Qué hace EAAT2 en condiciones normales (simportador electrogénico)
- Estequiometría: 3Na⁺ + 1H⁺ + Glu⁻ → interior; 1K⁺ → exterior
- Corriente neta generada (electrogénica → despolarización astrocítica leve)
- Parámetros electrofisiológicos que lo gobiernan:
  - Potencial de inversión estimado (~-30 mV)
  - Dependencia del gradiente electroquímico de Na⁺ (ecuación de Nernst/GHK)
  - Conductancia aniónica intrínseca
- Por qué su falla (DEE41) es crítica: acumulación de Glu → NMDA → hiperexcitabilidad

**Ecuaciones a incluir:**
- Ecuación de Nernst para Na⁺: `E_Na = (RT/F)·ln([Na⁺]_e/[Na⁺]_i)` 
- Descripción cualitativa de GHK para el potencial de membrana

---

### 🔵 Integrante 3 — Metodología: Fase I (Selección y Caracterización)
**Extensión sugerida:** ~0.5 página

**Contenido:**
- Identificación de la proteína: SLC1A2, UniProt P43004, cromosoma 11p13, 574 aa
- Estructura PDB de referencia: 7XR4 (cryo-EM, Zhang et al. 2022)
- Parámetros base obtenidos de la literatura:
  - pI = 6.09, PM = 62.1 kDa
  - Estequiometría de transporte
  - Localización: astrocitos perisinápticos
- Justificación clínica: DEE41 (OMIM #617105), mutaciones G82R y L85P

---

### 🔵 Integrante 4 — Metodología: Fase II (Herramientas Computacionales)
**Extensión sugerida:** ~0.5–0.75 página

**Contenido:**
Describir las 3 herramientas usadas en el proyecto (tabla formato):

| Categoría | Herramienta | Uso específico |
|-----------|-------------|----------------|
| Bioinformática | Python / Biopython | Extracción secuencia P43004, cálculo pI, perfil Kyte-Doolittle, búsqueda de motivos |
| Modelado estructural | PyMOL + FoldX | Carga PDB 7XR4, visualización mutaciones G82R/L85P, cálculo ΔΔG |
| Simulación | Modelo Hodgkin-Huxley (Python/SciPy) | Comparación excitabilidad neuronal WT vs. mutante DEE41 |

- Justificar por qué se eligieron estas herramientas
- Mencionar que el análisis supera la limitación 1D identificada en el trabajo previo

---

### 🔵 Integrante 5 — Metodología: Fase III (Protocolo de Simulación HH)
**Extensión sugerida:** ~0.75 página

**Contenido:**
- Descripción del modelo HH clásico + término I_EAAT2
- Condiciones iniciales: V₀ = -65 mV, estado estacionario m, h, n
- Parámetros del estímulo: pulso de 10 μA/cm², t = 10–90 ms
- Qué variable se modificó: g_EAAT2 (0.08 → 0 mS/cm²) y ΔI_NMDA (+1.5 μA/cm²)
- Criterio de validación: comparar frecuencia de disparo WT con literatura (Hodgkin & Huxley, 1952)
- Para el protocolo de PyMOL/FoldX: resumir los pasos en 3–5 oraciones con referencia al protocolo completo en el repositorio

---

### 🔵 Integrante 6 — Resultados: Análisis Bioinformático (Tarea 2)
**Extensión sugerida:** ~0.75 página

**Contenido (integrar resultados ya calculados en Tarea 2):**
- Tabla de propiedades fisicoquímicas (pI, PM, cargas +/-)
- Resultado de búsqueda de motivos: 0 ocurrencias S4, confirmación HP1 (ASS 362-364), sitios Na⁺ (D399, D486)
- Gráfica de perfil Kyte-Doolittle (imagen del notebook)
- **Sin interpretación aún** — solo describir los números y lo que muestra cada figura

---

### 🔵 Integrante 7 — Resultados: Modelado Estructural (PyMOL/FoldX)
**Extensión sugerida:** ~0.75–1 página

**Contenido:**
- Capturas de PyMOL: WT overview, sitio G82/L85, superposición WT vs. G82R
- Tabla de ΔΔG calculado por FoldX para G82R y L85P
- RMSD de la superposición (valor numérico del comando align)
- Descripción objetiva de qué muestra cada imagen — SIN interpretación todavía
  - Ej: "La Figura X muestra que el residuo 82 en la conformación silvestre es glicina (sin cadena lateral); en la mutante G82R se observa la cadena lateral de arginina ocupando el espacio hidrofóbico de TM2."

---

### 🔵 Integrante 8 — Resultados: Simulación HH
**Extensión sugerida:** ~0.75–1 página

**Contenido:**
- Panel A de la figura (potencial de membrana WT vs. mutante)
- Tabla de resultados cuantitativos (spikes, frecuencia, Vrest, Vpico)
- Panel F (diagrama de fase)
- Descripción objetiva: "La simulación muestra que bajo las condiciones del mutante DEE41, la neurona genera X potenciales de acción vs. Y en el silvestre, para el mismo estímulo de entrada."
- **Sin interpretación** — eso va en Discusión

---

### 🔵 Integrante 9 — Discusión + Conclusión + Bibliografía
**Extensión sugerida:** ~1.5 páginas

**Discusión (la sección de mayor peso — 15% de la rúbrica):**

Conectar los resultados con los fundamentos del curso:

1. **ΔΔG y deformación estructural → pérdida de I_EAAT2**
   - Glicina (G) en posición 82 es el aminoácido más flexible; su sustitución por Arg (+carga) en un ambiente hidrofóbico genera choque estérico → ΔΔG > 2 kcal/mol → dominio de trimerización desestabilizado
   - L85P: la prolina no puede adoptar los ángulos φ/ψ de una α-hélice → ruptura conformacional de TM2

2. **Pérdida de I_EAAT2 → circuito RC de membrana**
   - La corriente electrogénica de EAAT2 actúa como componente estabilizador en el circuito RC
   - Su abolición equivale a aumentar la conductancia de fuga: ΔG_L → τ = Rm·Cm disminuye → membrana más "reactiva"
   - Relacionar con la ecuación: `Cm·dV/dt = I_ext - I_iónicas - I_EAAT2`

3. **Acumulación de Glu → sobreactivación NMDA → potencial de Goldman**
   - El glutamato extracelular activa receptores NMDA (permeables a Na⁺ y Ca²⁺)
   - Esto desplaza el potencial de membrana hacia valores menos negativos (GHK)
   - Resultado de la simulación: X PAs más en el mutante para el mismo estímulo

4. **Hiperexcitabilidad → DEE41**
   - El diagrama de fase (Panel F) muestra mayor área del lazo en el mutante → más energía disipada por ciclo → tren de PAs sostenido → correlato electrofisiológico de las crisis epilépticas

**Conclusión:**
- Resumen del impacto: las mutaciones G82R/L85P no solo cambian la estructura, sino que traducen esa deformación 3D en un cambio funcional cuantificable: X PAs adicionales por estímulo
- Posibles soluciones terapéuticas:
  1. Activadores de EAAT2: ceftriaxona (antibiótico betalactámico que aumenta expresión de SLC1A2)
  2. Antagonistas de NMDA (memantina) para reducir la sobreactivación postsináptica
  3. Terapia génica dirigida al promotor de SLC1A2

**Bibliografía (mínimo 5, formato APA):**
1. Zhang, Z., et al. (2022). Structural basis of ligand binding modes of human EAAT2. *Nature Communications*, 13(1), 3329.
2. Hodgkin, A. L., & Huxley, A. F. (1952). A quantitative description of membrane current. *The Journal of Physiology*, 117(4), 500–544.
3. Tanaka, K., et al. (1997). Epilepsy and exacerbation of brain injury in mice lacking the glutamate transporter GLT-1. *Science*, 276(5319), 1699–1702.
4. Schymkowitz, J., et al. (2005). The FoldX web server: an online force field. *Nucleic Acids Research*, 33(Suppl 2), W382–W388.
5. UniProt Consortium (2025). UniProtKB — P43004 (EAA2_HUMAN). Recuperado de https://www.uniprot.org/uniprotkb/P43004
6. OMIM (2023). #617105 — Developmental and Epileptic Encephalopathy 41; DEE41. Johns Hopkins University.
7. Rothstein, J. D., et al. (2005). β-Lactam antibiotics offer neuroprotection by increasing glutamate transporter expression. *Nature*, 433(7021), 73–77.

---

## Checklist Final (antes de entregar)

- [ ] Fuente: 14 pt (títulos/apartados), 12 pt (texto), 11 pt (pies de figura/tabla)
- [ ] Límite: 5–10 páginas en total
- [ ] Abstract: exactamente 1 párrafo, máximo 250 palabras
- [ ] Mínimo 5 figuras/tablas originales del equipo (no de internet)
- [ ] Bibliografía: mínimo 5 fuentes indexadas (PubMed/Scopus/WoS), formato APA
- [ ] No hay wikis, blogs ni fuentes sin revisión por pares
- [ ] Todas las secciones obligatorias presentes
- [ ] Sin plagio (parafraseado y citado correctamente)
