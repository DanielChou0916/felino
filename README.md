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
Computes **fatigue energy** using the **mean load method**.
- Supports three accumulation modes:
  1. `monotonic` — No fatigue accumulation.
  2. `Fatigue` — Time-step-sensitive energy accumulation based on comparing $\Psi_t$ and  $\Psi_{t-1}$. Requires careful control of $\Delta t$ to capture cyclic peaks.
  3. `FatigueCLA` — Cycle-count-based accumulation. Less sensitive to $\Delta t$, but also less accurate.
- **Note:** This object only computes fatigue energy. Users must define additional materials to compute the fatigue function.

**Example Usage:**
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

