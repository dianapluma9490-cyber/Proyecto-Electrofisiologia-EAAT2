"""
=============================================================================
hh_eaat2_simulation.py
Simulación de Hodgkin-Huxley: EAAT2 Silvestre vs. Mutante DEE41 (G82R/L85P)

Universidad de Guadalajara — CUTlajomulco
Electrofisiología Molecular I — 2026A
=============================================================================

CONTEXTO BIOFÍSICO
------------------
El transportador EAAT2 (SLC1A2) elimina glutamato de la hendidura sináptica
mediante simporte electrogénico (3Na⁺ + 1H⁺ + Glu⁻ → intracelular; 1K⁺ → extracelular).
Su falla (mutaciones G82R/L85P en DEE41) impide la recaptura de glutamato,
elevando [Glu]_extra y activando receptores NMDA de forma sostenida.

En términos de circuito RC:
- La corriente EAAT2 funcional actúa como un "sumidero" que reduce la
  conductancia excitatoria efectiva de la membrana postsináptica.
- Su falla equivale a aumentar la conductancia de fuga (g_leak) y el
  potencial de equilibrio efectivo, desplazando el umbral de disparo.

MODELO IMPLEMENTADO
-------------------
Se usa el modelo clásico de Hodgkin-Huxley con un término adicional I_EAAT2
que modela la corriente electrogénica del transportador:

  Cm * dV/dt = I_ext - I_Na - I_K - I_L - I_EAAT2

Donde:
  I_Na  = g_Na * m³ * h * (V - E_Na)
  I_K   = g_K  * n⁴  * (V - E_K)
  I_L   = g_L  *      (V - E_L)
  I_EAAT2 = corriente electrogénica del transportador (positiva = hacia interior)

Condición silvestre:  EAAT2 funcional → I_EAAT2 activa, [Glu]_extra baja
Condición mutante:    EAAT2 no funcional → I_EAAT2 = 0, [Glu]_extra alta
                      → g_NMDA aumenta → mayor I_ext efectiva → hiperexcitabilidad
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import odeint

# ─────────────────────────────────────────────────────────────────────────────
# 1. PARÁMETROS BIOFÍSICOS (unidades: mV, ms, μA/cm², mS/cm², μF/cm²)
# ─────────────────────────────────────────────────────────────────────────────

# Parámetros estándar de Hodgkin-Huxley (axón de calamar, 6.3°C)
Cm    = 1.0     # μF/cm²  — capacitancia de membrana
g_Na  = 120.0   # mS/cm²  — conductancia máxima Na⁺
g_K   = 36.0    # mS/cm²  — conductancia máxima K⁺
g_L   = 0.3     # mS/cm²  — conductancia de fuga (leak)
E_Na  = 50.0    # mV      — potencial de equilibrio Na⁺  (Nernst)
E_K   = -77.0   # mV      — potencial de equilibrio K⁺   (Nernst)
E_L   = -54.4   # mV      — potencial de equilibrio leak

# Parámetros del transportador EAAT2
# La corriente electrogénica se modela como proporcional a la
# disponibilidad funcional del transportador y al gradiente de Na⁺
g_EAAT2_wt  = 0.08   # mS/cm²  — conductancia efectiva EAAT2 silvestre
g_EAAT2_mut = 0.0    # mS/cm²  — conductancia efectiva EAAT2 mutante (DEE41)
E_EAAT2     = -30.0  # mV      — potencial de inversión del simportador
                     #           (estimado considerando 3Na⁺ in / 1K⁺ out)

# Efecto de la acumulación de glutamato en el mutante:
# La falla de EAAT2 → [Glu]_extra aumenta → mayor activación NMDA
# Se modela como un aumento en la corriente externa efectiva (I_ext_NMDA)
delta_I_NMDA = 1.5   # μA/cm²  — corriente NMDA adicional por excitotoxicidad

# ─────────────────────────────────────────────────────────────────────────────
# 2. FUNCIONES DE COMPUERTA α y β (cinética de Hodgkin-Huxley)
# ─────────────────────────────────────────────────────────────────────────────

def alpha_m(V):
    """Tasa de apertura compuerta m (activación Na⁺)"""
    dV = V + 40.0
    if abs(dV) < 1e-7:
        return 1.0
    return 0.1 * dV / (1.0 - np.exp(-dV / 10.0))

def beta_m(V):
    """Tasa de cierre compuerta m (activación Na⁺)"""
    return 4.0 * np.exp(-(V + 65.0) / 18.0)

def alpha_h(V):
    """Tasa de apertura compuerta h (inactivación Na⁺)"""
    return 0.07 * np.exp(-(V + 65.0) / 20.0)

def beta_h(V):
    """Tasa de cierre compuerta h (inactivación Na⁺)"""
    return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))

def alpha_n(V):
    """Tasa de apertura compuerta n (activación K⁺)"""
    dV = V + 55.0
    if abs(dV) < 1e-7:
        return 0.1
    return 0.01 * dV / (1.0 - np.exp(-dV / 10.0))

def beta_n(V):
    """Tasa de cierre compuerta n (activación K⁺)"""
    return 0.125 * np.exp(-(V + 65.0) / 80.0)

# Vectorizar funciones de compuerta para uso con arrays de NumPy
alpha_m_v = np.vectorize(alpha_m)
beta_m_v  = np.vectorize(beta_m)
alpha_h_v = np.vectorize(alpha_h)
beta_h_v  = np.vectorize(beta_h)
alpha_n_v = np.vectorize(alpha_n)
beta_n_v  = np.vectorize(beta_n)

# ─────────────────────────────────────────────────────────────────────────────
# 3. SISTEMA DE ECUACIONES DIFERENCIALES
# ─────────────────────────────────────────────────────────────────────────────

def hodgkin_huxley(y, t, I_ext, g_eaat2, include_nmda_effect):
    """
    Sistema de 4 EDOs del modelo HH con término EAAT2.

    Estado: y = [V, m, h, n]
    
    Parámetros:
    -----------
    I_ext            : float — corriente externa aplicada (μA/cm²)
    g_eaat2          : float — conductancia EAAT2 (0 = mutante, 0.08 = silvestre)
    include_nmda_effect : bool — True para simular excitotoxicidad por acumulación Glu
    """
    V, m, h, n = y

    # Corrientes iónicas principales
    I_Na = g_Na * (m**3) * h * (V - E_Na)
    I_K  = g_K  * (n**4)     * (V - E_K)
    I_L  = g_L  *              (V - E_L)

    # Corriente electrogénica del transportador EAAT2
    # Positiva hacia el interior → hiperpolarizante neta
    I_EAAT2 = g_eaat2 * (V - E_EAAT2)

    # Corriente NMDA extra por acumulación de glutamato (solo en mutante)
    I_NMDA_extra = delta_I_NMDA if include_nmda_effect else 0.0

    # Ecuación de membrana: Cm * dV/dt = I_ext - I_iónicas
    dVdt = (I_ext + I_NMDA_extra - I_Na - I_K - I_L - I_EAAT2) / Cm

    # Ecuaciones de compuertas
    dmdt = alpha_m(V) * (1.0 - m) - beta_m(V) * m
    dhdt = alpha_h(V) * (1.0 - h) - beta_h(V) * h
    dndt = alpha_n(V) * (1.0 - n) - beta_n(V) * n

    return [dVdt, dmdt, dhdt, dndt]

# ─────────────────────────────────────────────────────────────────────────────
# 4. CONDICIONES INICIALES Y TIEMPO DE SIMULACIÓN
# ─────────────────────────────────────────────────────────────────────────────

# Estado de reposo (valores de estado estacionario a V = -65 mV)
V0 = -65.0
m0 = alpha_m(V0) / (alpha_m(V0) + beta_m(V0))
h0 = alpha_h(V0) / (alpha_h(V0) + beta_h(V0))
n0 = alpha_n(V0) / (alpha_n(V0) + beta_n(V0))

y0 = [V0, m0, h0, n0]

# Vector de tiempo: 0 a 100 ms, pasos de 0.01 ms
t = np.arange(0, 100, 0.01)

# Corriente externa: pulso de 10 μA/cm² entre t = 10 ms y t = 90 ms
def get_I_ext(t_array, I_amp=10.0, t_start=10.0, t_end=90.0):
    return np.where((t_array >= t_start) & (t_array <= t_end), I_amp, 0.0)

I_ext_array = get_I_ext(t)

# ─────────────────────────────────────────────────────────────────────────────
# 5. INTEGRACIÓN NUMÉRICA (método adaptativo de SciPy)
# ─────────────────────────────────────────────────────────────────────────────

print("=" * 60)
print("SIMULACIÓN HH — EAAT2 (SLC1A2) DEE41")
print("=" * 60)

# Función auxiliar para integrar con corriente variable
def simulate_hh(g_eaat2, nmda_effect, I_amplitude=10.0):
    """Integra el sistema HH paso a paso con corriente externa variable."""
    sol = np.zeros((len(t), 4))
    sol[0] = y0

    for i in range(1, len(t)):
        dt_span = [t[i-1], t[i]]
        I_step = I_ext_array[i-1]

        result = odeint(
            hodgkin_huxley, sol[i-1], dt_span,
            args=(I_step, g_eaat2, nmda_effect),
            rtol=1e-6, atol=1e-9
        )
        sol[i] = result[-1]

    return sol

# --- Condición 1: EAAT2 Silvestre (WT) ---
print("\n[1/2] Simulando EAAT2 silvestre (g_EAAT2 = {:.3f} mS/cm²)...".format(g_EAAT2_wt))
sol_wt  = simulate_hh(g_EAAT2_wt,  nmda_effect=False)
V_wt, m_wt, h_wt, n_wt = sol_wt.T
print("      ✓ Completado")

# --- Condición 2: EAAT2 Mutante DEE41 ---
print("[2/2] Simulando EAAT2 mutante DEE41 (g_EAAT2 = {:.3f} mS/cm², +NMDA)...".format(g_EAAT2_mut))
sol_mut = simulate_hh(g_EAAT2_mut, nmda_effect=True)
V_mut, m_mut, h_mut, n_mut = sol_mut.T
print("      ✓ Completado")

# ─────────────────────────────────────────────────────────────────────────────
# 6. MÉTRICAS CUANTITATIVAS
# ─────────────────────────────────────────────────────────────────────────────

def count_spikes(V_trace, threshold=0.0):
    """Cuenta potenciales de acción usando cruce ascendente por umbral."""
    crossings = np.where(np.diff((V_trace > threshold).astype(int)) == 1)[0]
    return len(crossings)

def get_firing_rate(V_trace, t_array, threshold=0.0):
    """Calcula la frecuencia de disparo en Hz."""
    n_spikes = count_spikes(V_trace, threshold)
    duration_s = (t_array[-1] - t_array[0]) / 1000.0  # ms → s
    return n_spikes / duration_s

def get_resting_potential(V_trace, t_array, t_before_stim=9.0):
    """Potencial de reposo promedio antes del estímulo."""
    idx = t_array < t_before_stim
    return np.mean(V_trace[idx])

n_spikes_wt  = count_spikes(V_wt)
n_spikes_mut = count_spikes(V_mut)
fr_wt  = get_firing_rate(V_wt, t)
fr_mut = get_firing_rate(V_mut, t)
Vrest_wt  = get_resting_potential(V_wt, t)
Vrest_mut = get_resting_potential(V_mut, t)
Vpeak_wt  = np.max(V_wt)
Vpeak_mut = np.max(V_mut)

print("\n" + "─" * 60)
print("RESULTADOS CUANTITATIVOS")
print("─" * 60)
print(f"{'Parámetro':<35} {'Silvestre (WT)':<18} {'Mutante DEE41':<18}")
print("─" * 60)
print(f"{'Potenciales de acción (#)':<35} {n_spikes_wt:<18} {n_spikes_mut:<18}")
print(f"{'Frecuencia de disparo (Hz)':<35} {fr_wt:<18.1f} {fr_mut:<18.1f}")
print(f"{'Potencial de reposo (mV)':<35} {Vrest_wt:<18.2f} {Vrest_mut:<18.2f}")
print(f"{'Amplitud pico (mV)':<35} {Vpeak_wt:<18.1f} {Vpeak_mut:<18.1f}")
print("─" * 60)

# ─────────────────────────────────────────────────────────────────────────────
# 7. VISUALIZACIÓN
# ─────────────────────────────────────────────────────────────────────────────

# Paleta de colores (inspirada en la paleta azul de los documentos de clase)
COLOR_WT   = '#1F4E79'   # Azul oscuro → silvestre
COLOR_MUT  = '#C0392B'   # Rojo         → mutante patológico
COLOR_STIM = '#95A5A6'   # Gris         → estímulo
COLOR_BG   = '#F8F9FA'

fig = plt.figure(figsize=(16, 14))
fig.patch.set_facecolor(COLOR_BG)

gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.45, wspace=0.35)

# ── Panel A: Potencial de membrana (comparación directa) ──────────────────
ax_V = fig.add_subplot(gs[0, :])
ax_V.set_facecolor('white')
ax_V.plot(t, V_wt,  color=COLOR_WT,  linewidth=1.5, label='EAAT2 Silvestre (WT)')
ax_V.plot(t, V_mut, color=COLOR_MUT, linewidth=1.5, label='EAAT2 Mutante DEE41 (G82R/L85P)', alpha=0.9)
ax_V.axvspan(10, 90, alpha=0.07, color='gold', label='Pulso de corriente (10 μA/cm²)')
ax_V.set_xlabel('Tiempo (ms)', fontsize=11)
ax_V.set_ylabel('Potencial de membrana (mV)', fontsize=11)
ax_V.set_title('A — Potencial de Membrana: EAAT2 Silvestre vs. Mutante DEE41',
               fontsize=13, fontweight='bold', pad=10)
ax_V.legend(loc='upper right', fontsize=9)
ax_V.set_xlim([0, 100])
ax_V.grid(True, alpha=0.3, linestyle='--')
ax_V.axhline(y=0, color='k', linewidth=0.5, linestyle=':')

# Anotación del conteo de spikes
ax_V.annotate(f'WT: {n_spikes_wt} PA  |  Mut: {n_spikes_mut} PA',
              xy=(0.02, 0.92), xycoords='axes fraction',
              fontsize=10, color='#2C3E50',
              bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

# ── Panel B: Compuertas m, h (Na⁺) ───────────────────────────────────────
ax_mh = fig.add_subplot(gs[1, 0])
ax_mh.set_facecolor('white')
ax_mh.plot(t, m_wt**3,       color=COLOR_WT,  linewidth=1.2, label='m³ WT (activación)')
ax_mh.plot(t, m_mut**3,      color=COLOR_MUT, linewidth=1.2, label='m³ Mut', linestyle='--')
ax_mh.plot(t, h_wt,          color=COLOR_WT,  linewidth=1.2, label='h WT (inactivación)', alpha=0.6)
ax_mh.plot(t, h_mut,         color=COLOR_MUT, linewidth=1.2, label='h Mut', linestyle='--', alpha=0.6)
ax_mh.set_xlabel('Tiempo (ms)', fontsize=10)
ax_mh.set_ylabel('Variable de compuerta', fontsize=10)
ax_mh.set_title('B — Cinética de Compuertas Na⁺ (m³, h)', fontsize=11, fontweight='bold')
ax_mh.legend(fontsize=8, ncol=2)
ax_mh.set_xlim([0, 100])
ax_mh.grid(True, alpha=0.3, linestyle='--')

# ── Panel C: Compuerta n (K⁺) ─────────────────────────────────────────────
ax_n = fig.add_subplot(gs[1, 1])
ax_n.set_facecolor('white')
ax_n.plot(t, n_wt**4,  color=COLOR_WT,  linewidth=1.2, label='n⁴ WT')
ax_n.plot(t, n_mut**4, color=COLOR_MUT, linewidth=1.2, label='n⁴ Mutante', linestyle='--')
ax_n.set_xlabel('Tiempo (ms)', fontsize=10)
ax_n.set_ylabel('n⁴ (conductancia K⁺ normalizada)', fontsize=10)
ax_n.set_title('C — Cinética de Compuerta K⁺ (n⁴)', fontsize=11, fontweight='bold')
ax_n.legend(fontsize=9)
ax_n.set_xlim([0, 100])
ax_n.grid(True, alpha=0.3, linestyle='--')

# ── Panel D: Corrientes iónicas WT ────────────────────────────────────────
ax_curr_wt = fig.add_subplot(gs[2, 0])
ax_curr_wt.set_facecolor('white')
I_Na_wt = g_Na * (m_wt**3) * h_wt * (V_wt - E_Na)
I_K_wt  = g_K  * (n_wt**4)        * (V_wt - E_K)
I_L_wt  = g_L  *                    (V_wt - E_L)
I_EAAT2_wt = g_EAAT2_wt * (V_wt - E_EAAT2)

ax_curr_wt.plot(t, -I_Na_wt,    color='#2980B9', linewidth=1.2, label='I_Na')
ax_curr_wt.plot(t, -I_K_wt,     color='#27AE60', linewidth=1.2, label='I_K')
ax_curr_wt.plot(t, -I_L_wt,     color='#8E44AD', linewidth=1.0, label='I_leak', alpha=0.7)
ax_curr_wt.plot(t, -I_EAAT2_wt, color='#F39C12', linewidth=1.2, label='I_EAAT2', linestyle='--')
ax_curr_wt.set_xlabel('Tiempo (ms)', fontsize=10)
ax_curr_wt.set_ylabel('Corriente (μA/cm²)', fontsize=10)
ax_curr_wt.set_title('D — Corrientes Iónicas (WT)', fontsize=11, fontweight='bold')
ax_curr_wt.legend(fontsize=8)
ax_curr_wt.set_xlim([0, 100])
ax_curr_wt.axhline(y=0, color='k', linewidth=0.5)
ax_curr_wt.grid(True, alpha=0.3, linestyle='--')

# ── Panel E: Corrientes iónicas Mutante ───────────────────────────────────
ax_curr_mut = fig.add_subplot(gs[2, 1])
ax_curr_mut.set_facecolor('white')
I_Na_mut = g_Na * (m_mut**3) * h_mut * (V_mut - E_Na)
I_K_mut  = g_K  * (n_mut**4)         * (V_mut - E_K)
I_L_mut  = g_L  *                      (V_mut - E_L)

ax_curr_mut.plot(t, -I_Na_mut, color='#E74C3C', linewidth=1.2, label='I_Na (mut)')
ax_curr_mut.plot(t, -I_K_mut,  color='#E67E22', linewidth=1.2, label='I_K (mut)')
ax_curr_mut.plot(t, -I_L_mut,  color='#8E44AD', linewidth=1.0, label='I_leak', alpha=0.7)
ax_curr_mut.plot(t, np.zeros_like(t), color='#F39C12', linewidth=1.2,
                 label='I_EAAT2 = 0 (abolida)', linestyle='--')
ax_curr_mut.set_xlabel('Tiempo (ms)', fontsize=10)
ax_curr_mut.set_ylabel('Corriente (μA/cm²)', fontsize=10)
ax_curr_mut.set_title('E — Corrientes Iónicas (Mutante DEE41)', fontsize=11, fontweight='bold')
ax_curr_mut.legend(fontsize=8)
ax_curr_mut.set_xlim([0, 100])
ax_curr_mut.axhline(y=0, color='k', linewidth=0.5)
ax_curr_mut.grid(True, alpha=0.3, linestyle='--')

# ── Panel F: Diagrama de fase (V vs dV/dt) ────────────────────────────────
ax_phase = fig.add_subplot(gs[3, :])
ax_phase.set_facecolor('white')

dV_wt  = np.gradient(V_wt,  t)
dV_mut = np.gradient(V_mut, t)

ax_phase.plot(V_wt,  dV_wt,  color=COLOR_WT,  linewidth=1.0, label='WT', alpha=0.8)
ax_phase.plot(V_mut, dV_mut, color=COLOR_MUT, linewidth=1.0, label='Mutante DEE41', alpha=0.8)
ax_phase.axvline(x=0,     color='k', linewidth=0.5, linestyle=':')
ax_phase.axhline(y=0,     color='k', linewidth=0.5, linestyle=':')
ax_phase.set_xlabel('V (mV)', fontsize=11)
ax_phase.set_ylabel('dV/dt (mV/ms)', fontsize=11)
ax_phase.set_title('F — Diagrama de Fase: V vs. dV/dt (Retrato de Excitabilidad)',
                   fontsize=11, fontweight='bold')
ax_phase.legend(fontsize=10)
ax_phase.grid(True, alpha=0.3, linestyle='--')

# Anotación sobre el diagrama de fase
ax_phase.annotate('Mayor área del lazo mutante\n→ mayor excitabilidad (hiperexcitabilidad DEE41)',
                  xy=(30, 200), fontsize=9, color='#C0392B',
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='#FADBD8', alpha=0.8))

# Título general de la figura
fig.suptitle(
    'Simulación Hodgkin-Huxley: Impacto de la Pérdida de Función de EAAT2 (SLC1A2)\n'
    'sobre la Excitabilidad Neuronal — Modelo DEE41 (G82R/L85P)',
    fontsize=14, fontweight='bold', y=0.98
)

plt.savefig('figures/HH_EAAT2_simulacion_completa.png', dpi=200,
            bbox_inches='tight', facecolor=COLOR_BG)
plt.savefig('figures/HH_EAAT2_simulacion_completa.pdf',
            bbox_inches='tight', facecolor=COLOR_BG)

print("\n✓ Figuras guardadas en figures/")
print("  → HH_EAAT2_simulacion_completa.png")
print("  → HH_EAAT2_simulacion_completa.pdf")

plt.show()

# ─────────────────────────────────────────────────────────────────────────────
# 8. TABLA RESUMEN PARA EL MANUSCRITO
# ─────────────────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("TABLA DE RESULTADOS — copiar al manuscrito")
print("=" * 60)
print(f"""
| Parámetro                          | WT (silvestre) | Mutante DEE41  |
|------------------------------------|----------------|----------------|
| g_EAAT2 (mS/cm²)                  | {g_EAAT2_wt:.3f}          | {g_EAAT2_mut:.3f}           |
| I_EAAT2 basal (μA/cm²)            | activa         | abolida (=0)   |
| ΔI_NMDA por excitotoxicidad        | 0 μA/cm²       | +{delta_I_NMDA} μA/cm²       |
| Potenciales de acción (#)          | {n_spikes_wt:<15} | {n_spikes_mut:<15} |
| Frecuencia de disparo (Hz)         | {fr_wt:<15.1f} | {fr_mut:<15.1f} |
| Potencial de reposo (mV)           | {Vrest_wt:<15.2f} | {Vrest_mut:<15.2f} |
| Amplitud pico PA (mV)              | {Vpeak_wt:<15.1f} | {Vpeak_mut:<15.1f} |
""")
print("=" * 60)
print("\nConexión con fundamentos del curso:")
print("• La abolición de I_EAAT2 → equivale a aumentar g_L en el circuito RC")
print("• El ΔI_NMDA desplaza el umbral de disparo → más PAs por el mismo estímulo")
print("• La cinética de compuertas m³h (Na⁺) no cambia intrínsecamente;")
print("  el cambio es en la corriente de fondo que modifica el Vm de partida")
