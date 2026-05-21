# Análisis Integrativo del Transportador SLC1A2 (EAAT2) y su Impacto sobre la Excitabilidad Neuronal en Encefalopatía Epiléptica DEE41

**Universidad de Guadalajara — Centro Universitario de Tlajomulco**  
**Ingeniería Biomédica, 6.° Semestre | Electrofisiología Molecular I — 2026A**  
**Profesora:** Elsa Patricia Magaña Cuevas  

---

## 📋 Descripción del Proyecto

Este repositorio contiene el pipeline computacional completo para el análisis integrativo del transportador de aminoácidos excitatorios EAAT2 (gen *SLC1A2*) en el contexto de la **Encefalopatía Epiléptica del Desarrollo (DEE41)**, causada por mutaciones puntuales G82R y L85P en el dominio de trimerización.

El análisis complementa la caracterización en 1D (bioinformática de secuencia) con una dimensión estructural 3D (modelado en PyMOL/FoldX) y una simulación funcional de excitabilidad neuronal (modelo de Hodgkin-Huxley modificado), respondiendo directamente a la limitación computacional identificada en las tareas previas.

---

## Estructura del Repositorio

```
proyecto_eaat2/
│
├── README.md                          ← Este archivo
│
├── notebooks/
│   ├── 01_bioinformatica_secuencia.ipynb   ← Análisis 1D (pI, hidropatía, motivos) — Tarea 2 integrada
│   └── 02_hodgkin_huxley_eaat2.ipynb       ← Simulación HH: silvestre vs. mutante DEE41
│
├── scripts/
│   ├── hh_eaat2_simulation.py         ← Script standalone de la simulación HH
│   └── analyze_sequence.py            ← Script standalone del análisis de secuencia
│
├── pymol_protocol/
│   ├── PROTOCOLO_PYMOL_FOLDX.md       ← Protocolo paso a paso para PyMOL y FoldX
│   └── pymol_commands.pml             ← Script de comandos PyMOL reproducible
│
├── manuscript/
│   └── GUIA_MANUSCRITO.md             ← Estructura del manuscrito con asignación por integrante
│
└── figures/                           ← Directorio para guardar figuras generadas
    └── .gitkeep
```

---

## Pregunta Científica Central

> **¿Cómo alteran las mutaciones G82R y L85P en el dominio de trimerización de EAAT2 la conductancia del transportador y, en consecuencia, la excitabilidad neuronal medible mediante el modelo de Hodgkin-Huxley?**

---

## Sistema Bioeléctrico Estudiado

| Parámetro | Valor |
|-----------|-------|
| Proteína | EAAT2 / SLC1A2 (GLT-1 en roedor) |
| UniProt | P43004 (EAA2_HUMAN) |
| Estructura PDB referencia | 7XR4 (cryo-EM, Zhang et al. 2022) |
| Mutaciones analizadas | G82R (TM2) y L85P (TM5) |
| Patología | DEE41 — Developmental and Epileptic Encephalopathy 41 |
| Mecanismo de falla | Pérdida de conductancia selectiva → excitotoxicidad → hiperexcitabilidad |

---

##  Fundamentos del Curso Aplicados

1. **Ecuación de Nernst / Goldman-Hodgkin-Katz** — Potencial de equilibrio del Na⁺ y K⁺ bajo condiciones de gradiente alterado por acumulación de glutamato extracelular.
2. **Circuito RC equivalente de membrana** — Modificación de la conductancia de membrana (g_leak) cuando EAAT2 falla y acumula cargas extracelulares.
3. **Cinética de compuertas de Hodgkin-Huxley** — Simulación del umbral de disparo y tren de potenciales de acción bajo condiciones de excitotoxicidad.
4. **Simportador electrogénico** — Corriente neta generada por el cotransporte 3Na⁺/1H⁺/Glu⁻ → 1K⁺.

---

##  Cómo Reproducir el Análisis

### Requisitos
```bash
pip install biopython numpy matplotlib scipy jupyter
```

### Ejecución rápida
```bash
# Simulación HH standalone
python scripts/hh_eaat2_simulation.py

# O abrir los notebooks en Google Colab / Jupyter
jupyter notebook notebooks/02_hodgkin_huxley_eaat2.ipynb
```

### Para el análisis estructural (PyMOL)
Ver instrucciones detalladas en [`pymol_protocol/PROTOCOLO_PYMOL_FOLDX.md`](pymol_protocol/PROTOCOLO_PYMOL_FOLDX.md)

---

##  División del Trabajo (9 integrantes)

Ver [`manuscript/GUIA_MANUSCRITO.md`](manuscript/GUIA_MANUSCRITO.md) para la asignación propuesta por sección.

---

## Referencias Clave

- Zhang, Z., et al. (2022). Structural basis of ligand binding modes of human EAAT2. *Nature Communications*, 13(1), 3329. DOI: 10.1038/s41467-022-31031-x
- Hodgkin, A. L., & Huxley, A. F. (1952). A quantitative description of membrane current. *Journal of Physiology*, 117(4), 500–544.
- Tanaka, K., et al. (1997). Epilepsy and exacerbation of brain injury in mice lacking the glutamate transporter GLT-1. *Science*, 276(5319), 1699–1702.
- UniProt Consortium (2025). P43004 · EAA2_HUMAN.
- OMIM #617105 — Developmental and Epileptic Encephalopathy 41.
