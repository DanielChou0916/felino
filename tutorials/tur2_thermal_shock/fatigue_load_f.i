########Equilibrium/Temp System parameters############
E = 21.0e4 #MPa
nu = 0.3

alpha_thermal = 8e-6 # Expansion mm/mm-k

T_ref = 300 #K
T_max = 1000 #K
########Equilibrium/Temp System parameters############

########PFF system############
gc0 = 2.7 # MPa-mm
l = 0.01
xi = 0
C0 = 2
L = 1e4
########PFF system############

########Fatigue############
alpha_critical0 = 45 #MPa
R=0.5
n=0.5
########Fatigue############
[Mesh]
  [gen]
    type = FileMeshGenerator
    file = ../mesh_files/single_notch_square.msh
  []
  [upper_left]
      type = BoundingBoxNodeSetGenerator
      input = gen
      bottom_left = '-0.501 0.4999 0'
      top_right = '-0.37499 0.5001 0'
      new_boundary = 'upper_left'
  []
  [create_sideset]
    type = SideSetsFromNodeSetsGenerator
    input = upper_left
  []
  #uniform_refine = 1
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  initial_steps = 1
  stop_time = 0
  max_h_level = 1
  [Markers]
    [marker]
      type = RotatedBoxMarker
      cx = 0.25
      cy = 0.25
      cz = 0
      lx = 0.51
      ly = 0.51
      lz = 0
      angle_z = 0
      angle_y = 0
      angle_x = 0
      inside = REFINE
      outside = DO_NOTHING
    []
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Actions/PFNonconserved]
  [./d]
    free_energy = F
    kappa = kappa_op
    mobility = L
    variable_mobility=false
    use_automatic_differentiation = true
    use_anisotropic_matrix = false
  [../]
[]

[AuxVariables]
  [./bounds_dummy]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./T]
  []
  ### Fatigue Related ###
  [./current_fatigue]
    order = CONSTANT
    family = MONOMIAL
  []
  [./bar_alpha]
    order = CONSTANT
    family = MONOMIAL
  []
  [./f_alpha]
    order = CONSTANT
    family = MONOMIAL
  []
  [./kappa_op]
    order = FIRST
    family = MONOMIAL
  []
  [./n_cycle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./current_fatigue]
    type = ADMaterialRealAux
    variable = current_fatigue
    property = current_fatigue
  [../]
  [./bar_alpha]
    type = ADMaterialRealAux
    variable = bar_alpha
    property = bar_alpha
    execute_on = timestep_end
  [../]
  [./f_alpha]
    type = ADMaterialRealAux
    variable = f_alpha
    property = f_alpha
  [../]
  [./kappa_op]
    type = ADMaterialRealAux
    variable = kappa_op
    property = kappa_op
  [../]
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
    base_name = uncracked
  []
  [./uncracked_strain]
    type = ADComputeFiniteStrain
    eigenstrain_names = thermal_expansion
    base_name = uncracked
  [../]
  [./trial_stress]
    type = ADComputeFiniteStrainElasticStress
    base_name = uncracked
  [../]
#####################Start : Temperature Coupling Relationship###########################
  [./gc] 
    type = ADParsedMaterial
    property_name = gc
    coupled_variables = 'T'
    expression = 'gc0 * (1 - b1 * ((T-T_ref)/T_max) + b2 * ((T-T_ref)/T_max)^2)'
    constant_names       = 'gc0    b1  b2  T_ref   T_max'
    constant_expressions = '${gc0} 1.8 1.1 ${T_ref} ${T_max}'
  [../]
  [./alpha_critical] 
    type = ADParsedMaterial
    property_name = alpha_critical
    coupled_variables = 'T'
    expression = 'alpha_critical0 * (1 - b1 * ((T-T_ref)/T_max) + b2 * ((T-T_ref)/T_max)^2)'
    constant_names       = 'alpha_critical0    b1  b2  T_ref   T_max'
    constant_expressions = '${alpha_critical0} 1.8 1.1 ${T_ref} ${T_max}'
  [../]
  [expansion]
    type = ADComputeThermalExpansionEigenstrain
    temperature = T
    thermal_expansion_coeff = ${alpha_thermal}
    stress_free_temperature = ${T_ref}
    eigenstrain_name = thermal_expansion
    base_name = uncracked
  []
  #####################End : Temperature Coupling Relationship###########################
  [./public_materials_forPF_model]
    type = ADGenericConstantMaterial
    prop_names =  '  l      xi    C0         L  ' # density k0=(2650 ,   3.1) are not used here, so do not included them
    prop_values = '  ${l}   ${xi} ${C0}      ${L}' #'0 2'#for AT2 # Or use '1 2.6666667' for AT1
  [../]
  [./fatigue_mats]
    type = ADGenericConstantMaterial
    prop_names =  "load_ratio  material_constant_n   "
    prop_values = "${R}        ${n}                  "
  [../]
  [./fatigue_variable]
    type = ADComputeFatigueEnergy
    #energy_calculation = mean_load
    uncracked_base_name = uncracked
    finite_strain_model = true
    multiply_by_D = false
    accumulation_mode = FatigueICLA
    N_cyc_variable = n_cycle
    acc_bar_psi_name = bar_alpha
    bar_psi_name = current_fatigue
  []
  [./fatigue_function]
    type = ADParsedMaterial
    material_property_names = 'bar_alpha alpha_critical'
    property_name = f_alpha
    expression = 'if(bar_alpha > alpha_critical, (2*alpha_critical/(bar_alpha + alpha_critical))^2, 1)'
  [] 


  [./degradation] # Define w(d)
    type = ADDerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'd'
    expression = '(1-d)^p*(1-k)+k'
    constant_names       = 'p k'
    constant_expressions = '2 1e-6'
    derivative_order = 2
  [../]
  [./local_fracture_energy] #Define psi_frac and alpha(d)
    type = ADDerivativeParsedMaterial
    property_name = local_fracture_energy
    coupled_variables = 'd'
    material_property_names = 'gc l xi C0 f_alpha'
    expression = '(xi*d+(1-xi)*d^2)* (gc / l)/C0* f_alpha'
    derivative_order = 2
  [../]
  [./define_kappa]
    type = ADParsedMaterial
    material_property_names = 'gc l C0 f_alpha'
    property_name = kappa_op
    expression = '2 * gc * l / C0* f_alpha'
  [../]
  [./cracked_stress]
    type = ADComputePFFStress
    decomposition = spectral
    #type = ADPFFRCEStress
    c = d
    E_name = E_el
    D_name = degradation
    use_current_history_variable = true
    uncracked_base_name = uncracked
    finite_strain_model = true
  [../]
  [./fracture_driving_energy]
    type = ADDerivativeSumMaterial
    coupled_variables = d
    sum_materials = 'E_el local_fracture_energy'
    derivative_order = 2
    property_name = F
  [../]
[]


[Bounds]
  [./d_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1.0
  [../]
  [./d_lower_bound]
    type = VariableOldValueBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  [../]
[]


[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -snes_type'
  petsc_options_value = 'lu vinewtonrsls'
  automatic_scaling = true
  nl_max_its = 40
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7
[]

[Outputs]
  print_linear_residuals = false
[]