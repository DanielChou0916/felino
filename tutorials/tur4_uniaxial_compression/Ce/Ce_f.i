# Aggregate (typical crushed stone / gravel)
E_agg        = 5.0e4        # MPa, Young's modulus (typical 40–60 GPa)
nu_agg       = 0.20         # Poisson's ratio
gc_agg       = 3e-2       # N/mm, fracture energy (typical 0.01–0.03 N/mm)
sigma_cs_agg = 80           # MPa, compressive strength
sigma_ts_agg = 7            # MPa, tensile strength (~0.08–0.1 fc)
delta_elp_agg= 0.16         # 

# Cement paste (ordinary paste for normal-strength concrete)
E_ce         = 2.5e4        # MPa, Young's modulus
nu_ce        = 0.20         # Poisson's ratio
gc_ce        = 1.25e-2         # N/mm, fracture energy
sigma_cs_ce  = 30           # MPa, compressive strength
sigma_ts_ce  = 2.5          # MPa, tensile strength (~0.08–0.1 fc)
delta_elp_ce = 0.64      # 

# Other PFF parameters
xi  = 1
C0  = 2.666666667
l   = 0.2
L   = 1e4




[Mesh]
  [gen]
    type = FileMeshGenerator
    file = ../two_phase.inp
  []
  [fix_node]
      type = BoundingBoxNodeSetGenerator
      input = gen
      bottom_left = '-0.1 -0.1 0'
      top_right = '0.0001 0.0001 0'
      new_boundary = 'fix_node'
  []
  uniform_refine = 2
  skip_partitioning = true
  construct_side_list_from_node_list=true
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
  [../]
[]

[AuxVariables]
  [./bounds_dummy]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]


[Materials]
  [./PFF_AGG]
    type = ADGenericConstantMaterial
    block = 0
    prop_names =  'gc        l     sigma_ts        sigma_cs         delta_elp '
    prop_values = '${gc_agg} ${l}  ${sigma_ts_agg} ${sigma_cs_agg}  ${delta_elp_agg}'
  [../]
  [./PFF_CE]
    type = ADGenericConstantMaterial
    block = 1
    prop_names =  'gc        l     sigma_ts        sigma_cs         delta_elp' 
    prop_values = '${gc_ce} ${l}   ${sigma_ts_ce}   ${sigma_cs_ce}     ${delta_elp_ce}'
  [../]
  [./public_materials_PFF]
    type = ADGenericConstantMaterial
    block = '0 1'
    prop_names =  'xi    C0      L     ' 
    prop_values = '${xi} ${C0}   ${L}  '
  [../]
  [./uncracked_strain]
    type = ADComputeSmallStrain
  [../]
  [./elasticity_tensorAgg]
    type = ADComputeIsotropicElasticityTensor 
    block = 0
    poissons_ratio = ${nu_agg}
    youngs_modulus = ${E_agg} 
  [../]
  [./elasticity_tensorCE]
    type = ADComputeIsotropicElasticityTensor 
    block = 1
    poissons_ratio = ${nu_ce}
    youngs_modulus = ${E_ce} 
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
    material_property_names = 'gc l xi C0'
    expression = '(xi*d+(1-xi)*d^2)* gc / C0 / l'
    derivative_order = 2
  [../]
  [./define_kappa]
    type = ADParsedMaterial
    material_property_names = 'gc l C0'
    property_name = kappa_op
    expression = '2 * gc * l / C0'
  [../]
  [./cracked_stress]
    type = ADLinearElasticPFFractureStress
    c = d
    E_name = E_el
    D_name = degradation 
    decomposition_type = none
    use_current_history_variable = true
    use_snes_vi_solver = true
  [../]
  [./additional_driving_force]
    type = ADComputeExtraDrivingForce
    D_name = degradation 
    Ce_name = Ce
  [] 
  [./fracture_driving_energy]
    type = ADDerivativeSumMaterialWithConstantOn1stOrder
    coupled_variables = d
    sum_materials = 'E_el local_fracture_energy'
    additional_sum_materials = 'Ce'
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