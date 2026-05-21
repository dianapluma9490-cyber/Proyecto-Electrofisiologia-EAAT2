# Protocolo de Análisis Estructural: PyMOL y FoldX
## Comparación Silvestre vs. Mutante EAAT2 (SLC1A2) — DEE41 (G82R / L85P)

**Universidad de Guadalajara — CUTlajomulco**  
**Electrofisiología Molecular I — 2026A**

---

## Contexto

Las mutaciones G82R (glicina→arginina en TM2) y L85P (leucina→prolina en TM5) causan DEE41.
Este protocolo permite **visualizar y cuantificar** cómo esas sustituciones deforman el dominio
de trimerización, alterando la conductancia del transportador.

La estructura de referencia es el modelo **AlphaFold** (AF-P43004-F1) o la estructura cryo-EM
**PDB: 7XR4** (Zhang et al., 2022, *Nature Communications*).

---

## PARTE 1 — PyMOL: Visualización y Comparación

### Requisitos
- PyMOL (versión 2.5+ o PyMOL Open-Source)
- Conexión a internet (para fetch automático desde PDB)
- O descarga manual desde: https://www.rcsb.org/structure/7XR4

---

### Paso 1: Cargar la estructura silvestre

```python
# En la consola de PyMOL:

# Opción A — Fetch directo desde PDB (requiere internet)
fetch 7XR4, name=EAAT2_WT, async=0

# Opción B — Cargar archivo local descargado
load /ruta/a/7XR4.pdb, EAAT2_WT
```

---

### Paso 2: Limpiar y preparar la estructura

```python
# Remover moléculas de agua y heteroátomos no esenciales
remove solvent
remove (EAAT2_WT and resn HOH)

# Solo trabajar con la cadena A (monómero funcional)
create EAAT2_WT_chainA, EAAT2_WT and chain A

# Remover cargas negativas, lipidos y otros (opcional)
# remove (EAAT2_WT and not polymer)

# Mostrar el monómero limpio
hide everything
show cartoon, EAAT2_WT_chainA
color marine, EAAT2_WT_chainA
```

---

### Paso 3: Localizar los residuos G82 y L85

```python
# Seleccionar residuos clave en el dominio TM2
select G82_site, (EAAT2_WT_chainA and resi 82)
select L85_site, (EAAT2_WT_chainA and resi 85)

# Mostrar como sticks con colores diferenciados
show sticks, G82_site
show sticks, L85_site
color yellow,    G82_site
color orange,    L85_site

# Mostrar etiquetas
label G82_site, "%s%s" % (resn, resi)
label L85_site, "%s%s" % (resn, resi)

# Centrar la vista en el dominio TM2
zoom G82_site, 12
```

---

### Paso 4: Generar la mutación in silico (PyMOL Mutagenesis Wizard)

```
Método A (Interfaz gráfica — recomendado para captura de pantalla):
1. Menú superior → Wizard → Mutagenesis
2. Clic sobre el residuo G82 en la estructura 3D
3. En el panel derecho: seleccionar "ARG" (arginina) en el menú desplegable
4. Clic en "Apply" → PyMOL generará el rotámero de menor energía
5. Guardar estructura: File → Export Molecule → EAAT2_G82R.pdb
6. Repetir para L85 → PRO (prolina)

Método B (Línea de comandos):
```

```python
# Duplicar estructura silvestre
copy EAAT2_G82R, EAAT2_WT_chainA
copy EAAT2_L85P, EAAT2_WT_chainA

# Aplicar mutación (PyMOL instalará los átomos del nuevo residuo)
# NOTA: este método es ilustrativo; FoldX (Parte 2) hace la optimización energética
alter (EAAT2_G82R and resi 82), resn="ARG"
alter (EAAT2_L85P and resi 85), resn="PRO"
```

---

### Paso 5: Superposición silvestre vs. mutante

```python
# Superponer la mutante sobre la silvestre (alineamiento por Cα)
align EAAT2_G82R, EAAT2_WT_chainA
# El RMSD reportado indica el grado de deformación estructural

# Colorear por diferencia
color marine,    EAAT2_WT_chainA
color firebrick, EAAT2_G82R

# Mostrar las dos estructuras superpuestas
show cartoon, EAAT2_WT_chainA
show cartoon, EAAT2_G82R

# Resaltar zona de cambio
select region_cambio, (EAAT2_G82R and resi 75-95)
show sticks, region_cambio
```

---

### Paso 6: Visualización del dominio de trimerización

```python
# El dominio de trimerización involucra TM2 (resi ~65-90) y TM5 (resi ~230-260)
# Verificar los rangos exactos en la anotación UniProt P43004

select TM2, (EAAT2_WT_chainA and resi 65-90)
select TM5, (EAAT2_WT_chainA and resi 230-260)

color cyan,   TM2
color salmon, TM5

# Interfaz de trimerización (contactos entre protómeros)
# Si trabajas con el homotrímero completo (3 cadenas A,B,C):
select interfaz_AB, (chain A within 4 of chain B)
color magenta, interfaz_AB
```

---

### Paso 7: Generar imagen para el manuscrito

```python
# Configurar fondo blanco y calidad
bg_color white
set ray_opaque_background, on
set antialias, 2
set ray_shadows, 1

# Orientar la vista
set_view (\
   1.0, 0.0, 0.0,\
   0.0, 1.0, 0.0,\
   0.0, 0.0, 1.0,\
   0.0, 0.0, -100.0,\
   0.0, 0.0, 0.0, 1.0 )

# Renderizar y guardar (resolución alta para figura)
ray 2400, 1800
png figures/pymol_EAAT2_WT_vs_G82R.png, dpi=300
```

---

## PARTE 2 — FoldX: Cuantificación Energética (ΔΔG)

FoldX calcula el cambio en energía libre de estabilidad (ΔΔG) por una mutación.
ΔΔG > 0 kcal/mol → la mutación desestabiliza la proteína.

### Instalación
Descargar FoldX 5 desde: https://foldxsuite.crg.eu/  
(Registro gratuito para uso académico)

---

### Paso 1: Reparar la estructura PDB

Antes de calcular mutaciones, la estructura debe ser reparada/optimizada por FoldX:

```bash
# En la terminal, desde el directorio donde está el PDB:
./foldx --command=RepairPDB --pdb=7XR4.pdb

# Genera: 7XR4_Repair.pdb
```

---

### Paso 2: Crear el archivo de mutaciones

Crear un archivo de texto llamado `individual_list.txt`:

```
# Formato: ResidueWT Chain ResidueNumber MutantAA;
# Para múltiples mutaciones: separar con coma, terminar con ;

GA82R;
LA85P;
```

---

### Paso 3: Ejecutar BuildModel (genera estructuras mutantes)

```bash
./foldx --command=BuildModel \
        --pdb=7XR4_Repair.pdb \
        --mutant-file=individual_list.txt \
        --numberOfRuns=5 \
        --out-pdb=true

# Genera:
# - 7XR4_Repair_1.pdb      (estructura mutante)
# - Average_7XR4_Repair.fxout  (energías promedio)
# - Dif_7XR4_Repair.fxout      (ΔΔG de cada mutación)
```

---

### Paso 4: Leer los resultados de ΔΔG

```bash
# Ver resultados de estabilidad:
cat Dif_7XR4_Repair.fxout
```

Interpretar la columna **"total energy"** en la tabla de salida.
El valor representa ΔΔG en **kcal/mol**:

| Mutación | ΔΔG (kcal/mol) | Interpretación |
|----------|----------------|----------------|
| G82R     | > 2.0          | Muy desestabilizante (valor esperado) |
| L85P     | > 3.0          | Muy desestabilizante — la prolina interrumpe la α-hélice |

---

### Paso 5: Visualizar la mutante de FoldX en PyMOL

```python
# Cargar la estructura mutante generada por FoldX
load 7XR4_Repair_1.pdb, EAAT2_FoldX_mut

# Superponer sobre silvestre
align EAAT2_FoldX_mut, EAAT2_WT_chainA

# Mostrar y comparar
show cartoon, all
color marine,    EAAT2_WT_chainA
color firebrick, EAAT2_FoldX_mut
```

---

## PARTE 3 — Conexión con Fundamentos del Curso

### ¿Cómo se traduce el ΔΔG al circuito RC de membrana?

```
ΔΔG elevado (G82R, L85P)
        ↓
Dominio de trimerización deformado
        ↓
Pérdida de conformación "inward-facing" óptima
        ↓
Reducción o abolición de I_EAAT2
        ↓
[Glu]_extra ↑ → NMDA sobreactivado
        ↓
Incremento de conductancia catiónica (g_NMDA ↑)
        ↓
En el circuito RC: equivale a disminuir Rm (resistencia de membrana)
        ↓
τ = Rm·Cm disminuye → la membrana responde más rápido
        ↓
Umbral de disparo más fácil de alcanzar → hiperexcitabilidad → DEE41
```

### Ecuación de Goldman aplicada

Cuando EAAT2 falla, la homeostasis iónica del espacio sináptico se altera.
La ecuación de Goldman-Hodgkin-Katz predice el potencial de membrana resultante:

```
Vm = (RT/F) · ln[ (P_Na·[Na+]_e + P_K·[K+]_i) / (P_Na·[Na+]_i + P_K·[K+]_e) ]

- La falla de EAAT2 no afecta directamente [Na+] o [K+],
  sino que aumenta [Glu]_e, lo que activa receptores NMDA
  (permeables a Na+ y Ca2+), desplazando el equilibrio hacia
  potenciales de membrana menos negativos.
```

---

## PARTE 4 — Script PyMOL Completo (reproducible)

Ver archivo: `pymol_protocol/pymol_commands.pml`

---

## Referencias

- Zhang, Z., et al. (2022). Structural basis of ligand binding modes of human EAAT2. *Nature Communications*, 13(1), 3329.
- Schymkowitz, J., et al. (2005). The FoldX web server. *Nucleic Acids Research*, 33, W382–W388.
- UniProt P43004 — https://www.uniprot.org/uniprot/P43004
- RCSB PDB 7XR4 — https://www.rcsb.org/structure/7XR4
- OMIM #617105 — DEE41
