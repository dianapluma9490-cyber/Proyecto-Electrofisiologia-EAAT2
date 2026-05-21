# pymol_commands.pml
# Script PyMOL reproducible — EAAT2 (SLC1A2) WT vs. Mutante DEE41
# Electrofisiología Molecular I — 2026A, CUTlajomulco
#
# Uso: En PyMOL → File → Run Script → seleccionar este archivo
#      O en la consola: @pymol_commands.pml

# ─────────────────────────────────────────────────────────────────────────────
# 1. CARGAR ESTRUCTURA
# ─────────────────────────────────────────────────────────────────────────────

# Fetch desde PDB (estructura cryo-EM de EAAT2 humano, Zhang et al. 2022)
fetch 7XR4, name=EAAT2_WT, async=0

# Limpiar: remover solvente
remove solvent

# Crear objeto solo con cadena A (monómero)
create EAAT2_WT_A, EAAT2_WT and chain A

# ─────────────────────────────────────────────────────────────────────────────
# 2. ESTILO BASE
# ─────────────────────────────────────────────────────────────────────────────

bg_color white
hide everything
show cartoon, EAAT2_WT_A
color marine, EAAT2_WT_A
set cartoon_fancy_helices, 1
set ray_opaque_background, on

# ─────────────────────────────────────────────────────────────────────────────
# 3. RESALTAR RESIDUOS CRÍTICOS
# ─────────────────────────────────────────────────────────────────────────────

# Sitios de mutación DEE41
select G82, (EAAT2_WT_A and resi 82)
select L85, (EAAT2_WT_A and resi 85)
show sticks, G82
show sticks, L85
color yellow,  G82
color orange,  L85
label G82, "%s%s" % (resn, resi)
label L85, "%s%s" % (resn, resi)

# Dominio HP1 (transporte de glutamato)
select HP1, (EAAT2_WT_A and resi 362-364)
show sticks, HP1
color green, HP1

# Sitios de coordinación de sodio
select D399, (EAAT2_WT_A and resi 399)
select D486, (EAAT2_WT_A and resi 486)
show sticks, D399
show sticks, D486
color red, D399
color red, D486

# Dominios transmembrana TM2 y TM5
select TM2, (EAAT2_WT_A and resi 65-90)
select TM5, (EAAT2_WT_A and resi 230-260)
color cyan,   TM2
color salmon, TM5

# ─────────────────────────────────────────────────────────────────────────────
# 4. GUARDAR IMAGEN SILVESTRE
# ─────────────────────────────────────────────────────────────────────────────

zoom EAAT2_WT_A
ray 2400, 1800
png figures/pymol_01_EAAT2_WT_overview.png, dpi=300

zoom G82, 10
ray 1200, 900
png figures/pymol_02_EAAT2_G82_site.png, dpi=300

# ─────────────────────────────────────────────────────────────────────────────
# 5. GENERAR MUTANTE G82R CON WIZARD
# (Ejecutar manualmente: Wizard > Mutagenesis > G82 > ARG > Apply)
# O cargar estructura mutante generada por FoldX:
# load 7XR4_Repair_1.pdb, EAAT2_G82R
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# 6. SUPERPOSICIÓN (ejecutar DESPUÉS de cargar mutante)
# ─────────────────────────────────────────────────────────────────────────────

# align EAAT2_G82R, EAAT2_WT_A
# color firebrick, EAAT2_G82R
# zoom all
# ray 2400, 1800
# png figures/pymol_03_superposicion_WT_vs_G82R.png, dpi=300

# ─────────────────────────────────────────────────────────────────────────────
# 7. NOTAS PARA EL MANUSCRITO
# ─────────────────────────────────────────────────────────────────────────────
# - Capturar RMSD reportado por el comando align (valor en consola)
# - El RMSD global es bajo, pero la región TM2 (resi 65-90) mostrará
#   mayor desviación — se puede calcular como:
#   rms_cur (EAAT2_G82R and resi 65-90), (EAAT2_WT_A and resi 65-90)
# - La mutación G82R introduce una cadena lateral cargada (+) de arginina
#   en una región hidrofóbica → desestabiliza empaquetamiento de hélices
# - L85P rompe la α-hélice TM2 (la prolina no puede adoptar ángulos φ/ψ
#   de una hélice regular) → colapsa la interfaz de trimerización

print("Script PyMOL cargado. Ver instrucciones en PROTOCOLO_PYMOL_FOLDX.md")
