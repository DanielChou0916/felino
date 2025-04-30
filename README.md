felino
=====

## 📦 Customized Objects Overview

### 📂 Materials

#### 1️⃣ `ADComputePFFStress`
A general-purpose tool for calculating **damaged elastic energy** and **damaged stress** in phase-field fracture models.
- Supports **finite strain** and **anisotropic materials**.
- Offers options for **spectral decomposition** or no decomposition.
- Requires the user to define an uncracked stress material as input.

**Example Usage:**
```ini
[./cracked_stress]
  type = ADComputePFFStress
  c = d
  E_name = E_el
  D_name = degradation
  decomposition = spectral
  use_current_history_variable = true
  uncracked_base_name = uncracked
  finite_strain_model = true
[../]
```

---

#### 2️⃣ `ADComputeFatigueEnergy`
Computes **fatigue energy** at the current step $\Psi_t$, along with its accumulated value $\bar{\alpha}_t$.

- **Three approaches for energy calculation**:
  1. `mean_load`: Default mode.  
     $\Psi_t = 2E\epsilon_{\text{max}}^2\left(\frac{1+R}{2}\right)^2\left(\frac{1+R}{2}\right)^n$,  
     where $\epsilon_{\text{max}}$ is the maximum principal strain.  
     Requires loading ratio $R$ and exponent $n$ to be defined in the generic material block.
  2. `elastic_energy`:  
     $\Psi_t = 0.5 \, \boldsymbol{\sigma} : \boldsymbol{\epsilon}$
  3. `spectral_activation`:  
     $\Psi_t = \boldsymbol{\sigma}^+ : \boldsymbol{\epsilon}$,  
     where $\boldsymbol{\sigma}^+$ is the tensile-activated component of stress obtained from spectral decomposition.

- **Four accumulation modes** for fatigue energy $\bar\alpha_t = \bar{\alpha}_{t-1} + \Delta \bar{\alpha}$:
  1. `Monotonic` — No fatigue accumulation.  
     $\Delta \bar{\alpha} = 0$
  2. `Fatigue` — *(default)*  
     Compares $\Psi_t$ with $\Psi_{t-1}$; if $\Psi_t > \Psi_{t-1}$, accumulate difference:  
     $\Delta \bar{\alpha} = \Psi_t - \Psi_{t-1}$, else 0  
     **Note:** Highly sensitive to time step $\Delta t$; must resolve cyclic peaks properly.
  3. `FatigueCLA` — Cycle-count-based linear accumulation.  
     $\bar{\alpha}_t = \Psi_t \times N_t$, where $N_t$ is the cycle count for the current step.
  4. `FatigueICLA` — *Incremental* cycle-count-based accumulation.  
     Uses increment $\Delta N = N_t - N_{t-1}$:  
     $\Delta \bar{\alpha} = \Psi_t \times \Delta N$, then  
     $\bar\alpha_t = \bar{\alpha}_{t-1} + \Delta \bar{\alpha}$

- **Note 1:** This object only calculates fatigue energy. Users must define additional materials (e.g., fatigue degradation functions) to utilize $\bar{\alpha}$.

- **Note 2:** For `FatigueCLA` and `FatigueICLA`, the loading boundary condition can be set to peak amplitude, and only the cycle count $N$ needs to be updated via an auxiliary variable.

**Example Usage: mean_load and Fatigue accumulation**
```ini
[./fatigue_variable]
  type = ADComputeFatigueEnergy
  uncracked_base_name = uncracked
  finite_strain_model = true
  multiply_by_D = false
  accumulation_mode = Fatigue
  acc_bar_psi_name = bar_alpha      # Output: accumulated energy
  bar_psi_name = current_fatigue    # Output: current energy
[]
```
**Example Usage: spectral_activation and FatigueICLA accumulation**
```ini
[./fatigue_variable]
  type = ADComputeFatigueEnergy
  energy_calculation = spectral_activation
  uncracked_base_name = uncracked
  finite_strain_model = true
  #D_name = #no need to set this if multiply_by_D = false
  multiply_by_D = false
  accumulation_mode = FatigueICLA
  N_cyc_variable = n_cycle
  acc_bar_psi_name = bar_alpha
  bar_psi_name = current_fatigue
[]
```

**Supplement Example: Fatigue function**
```
[./fatigue_function]
  type = ADParsedMaterial
  material_property_names = 'bar_alpha alpha_critical'
  property_name = f_alpha
  expression = 'if(bar_alpha > alpha_critical, (2*alpha_critical/(bar_alpha + alpha_critical))^2, 1)'
[]
```

---

### 📂 Actions & Kernels

#### 3️⃣ `Actions/ADNonconserved`
Automates the setup of Kernels for solving the **Allen-Cahn equation**.
- Use `use_grad_kappa = true` to include spatial gradients of material parameters when they are non-constant.
- Gradient terms ($\nabla \kappa$) are handled via `AuxKernels` and `AuxVariables`.
- Automatically detects mesh dimension (e.g., requires `grad_kappa_z` in 3D).

**Example Usage:**
```ini
[Actions/ADNonconserved]
  [./d]
    free_energy = F
    kappa = kappa_op
    mobility = L
    variable_mobility = false
    use_grad_kappa = true
    grad_kappa_x = dkappa_dx
    grad_kappa_y = dkappa_dy
  [../]
[]
```

