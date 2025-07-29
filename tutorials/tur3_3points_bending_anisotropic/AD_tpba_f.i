############Elastic Phase Field Parameters############
E = 21e4 #MPa # 210 GPa
nu = 0.3     #
gc = 2.7     #KJ/m2 = MPa. mm
l = 0.5   #mm

xi = 0#1
C0 = 2#2.666667
L = 1e4
######################################################

############Fatigue Parameters############
alpha_critical = 62.5 #MPa
R=0.5
n=0.5
##########################################
[Mesh]
  [gen]
    type = FileMeshGenerator
    file = ../mesh_files/three_points_metal.inp
  []
  construct_side_list_from_node_list=true
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  initial_steps = 2
  stop_time = 0
  max_h_level = 2
  [Markers]
    [marker]
      type = RotatedBoxMarker
      cx = 68
      cy = 2.5
      cz = 10
      lx = 15
      ly = 20
      lz = 50
      angle_z = 0
      angle_y = 0
      angle_x = 0
      inside = REFINE
      outside = DO_NOTHING
    []
  []
[]


[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Actions/PFNonconserved]
  [./d]
    free_energy = F
    kappa = kappa_op
    mobility = L
    variable_mobility=false
    use_automatic_differentiation = true
    use_anisotropic_matrix = true
    anisotropic_matrix = A_matrix
  [../]
[]

[AuxVariables]
  [./bounds_dummy]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
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
  [./anisotropy]
    type = ADAnisotropicDirector
    normal = "-0.70710678 0 0.70710678"
    coef = 30
    normalize_director = det_norm
    output_name = A_matrix
  []
  [./uncracked_strain]
    type = ADComputeFiniteStrain
    base_name = uncracked
  [../]
  [./trial_stress]
    type = ADComputeFiniteStrainElasticStress
    base_name = uncracked
  [../]
  [./pfbulkmat]
    type = ADGenericConstantMaterial
    prop_names =  'load_ratio  material_constant_n    alpha_critical'
    prop_values = '${R}        ${n}                   ${alpha_critical}' 
  [../]
  [./public_materials_forPF_model]
    type = ADGenericConstantMaterial
    prop_names =  'gc     l    xi    C0      L  ' 
    prop_values = '${gc}  ${l} ${xi} ${C0}   ${L}' 
  [../]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    poissons_ratio = ${nu}
    youngs_modulus = ${E} 
    base_name = uncracked
  [../]
  [./degradation] 
    type = ADDerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'd'
    expression = '(1-d)^p*(1-k)+k'
    constant_names       = 'p k'
    constant_expressions = '2 1e-6'
    derivative_order = 2
  [../]
  [./local_fracture_energy] 
    type = ADDerivativeParsedMaterial
    property_name = local_fracture_energy
    coupled_variables = 'd'
    material_property_names = 'gc l xi C0 f_alpha'
    expression = '(xi*d+(1-xi)*d^2) * (gc / l)/C0 * f_alpha'
    derivative_order = 2
  [../]
  [./fatigue_variable]
    type = ADComputeFatigueEnergy
    uncracked_base_name = uncracked
    finite_strain_model = true
    #D_name = #no need to set this if multiply_by_D = false
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
  [./define_kappa]
    type = ADParsedMaterial
    material_property_names = 'gc l C0 f_alpha'
    property_name = kappa_op
    expression = '2 * gc * l / C0 * f_alpha'
  [../]
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
  [./fracture_driving_energy]
    type = ADDerivativeSumMaterial
    #block = 1
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
  nl_max_its = 45
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7
[]

[Outputs]
  print_linear_residuals = false
[]