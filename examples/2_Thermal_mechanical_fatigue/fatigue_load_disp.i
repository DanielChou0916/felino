#reference: A thermo-mechanically coupled phase-field fatigue fracture model
# 3.2 SENT fatigue cracking under nonuniform temperature and cyclic mechanical loading
# See BCs block in fatigue_load_T.i
# (t - 3): Start to heating up the top when t=3, about 15 cycles
########Start: Equilibrium/Temp System parameters############
E = 34.0e4 #MPa
nu = 0.22

alpha_bulk = 8e-6 # Expansion m/m-K
T_ref = 300 #K
bc_stress = 150
##########End: Equilibrium/Temp System parameters#############

########Start: Fatigue cycle counting############
period = 0.2
num_cycle = 400
end_time = ${fparse period * num_cycle}
########End: Fatigue cycle counting############

########Start: Visualization Coupling relationship#############
T_max = 1000 #K
k0 = 3.0 #(J/10s)/mmK=0.3 #(J/s)/mmK
gc0 = 2.7 # MPa-mm
alpha_critical0 = 45 #MPa
########End: Visualization Coupling relationship#############


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
  uniform_refine = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[MultiApps]
  [temp]
    type = TransientMultiApp
    input_files = 'fatigue_load_T.i'
  []
  [crack]
    type = TransientMultiApp
    input_files = 'fatigue_load_f.i'
  []
[]

[Transfers]
  [to_d_T]
    type = MultiAppCopyTransfer
    to_multi_app = 'temp'
    source_variable = 'd'
    variable = 'd'
  []
  [from_T]
    type = MultiAppCopyTransfer
    from_multi_app = 'temp'
    source_variable = 'T'
    variable = 'T'
  []
  #loosely couple with damage
  [to_disp_x]
    type = MultiAppCopyTransfer 
    to_multi_app = 'crack'
    source_variable = 'disp_x'
    variable = 'disp_x'
  []
  [to_disp_y]
    type = MultiAppCopyTransfer
    to_multi_app = 'crack'
    source_variable = 'disp_y'
    variable = 'disp_y'
  []
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'd'
    variable = 'd'
  []
  [to_T_d]
    type = MultiAppCopyTransfer
    to_multi_app = 'crack'
    source_variable = 'T'
    variable = 'T'
  []
  [to_CLA]
    type = MultiAppCopyTransfer
    to_multi_app = 'crack'
    source_variable = 'n_cycle'
    variable = 'n_cycle'
  []
  [from_accumulate_fatigue]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'bar_alpha'
    variable = 'bar_alpha'
  []
  [from_fatigue_function]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'f_alpha'
    variable = 'f_alpha'
  []
[]




[AuxVariables]
  [./T]
  []
  [./d]
  []
  [./n_cycle]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./bar_alpha]
    order = CONSTANT
    family = MONOMIAL
  []
  [./f_alpha]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    add_variables = true
    strain = FINITE
    automatic_eigenstrain_names = false
    eigenstrain_names = thermal_expansion
    incremental = true
    additional_generate_output = 'stress_xx stress_yy stress_xy'
    use_automatic_differentiation=true
    strain_base_name = uncracked
    #decomposition_method = EigenSolution
  []
[]

[AuxKernels]
  [./n_cycle_aux]
    type = FunctionAux
    variable = n_cycle
    function = current_cycle
  [../]
[]

[Functions]
  [./current_cycle]
    type = ParsedFunction
    expression = 't / ${period}'
  [../]
[]
# No need to add stress divergence

[BCs]
  [pull_y]
    type = NeumannBC
    value = ${bc_stress}
    #type = DirichletBC
    #value = 0.001
    variable = disp_y
    boundary = upper_left
  []
  [fix_disp_x]
    type = DirichletBC
    variable = disp_x
    boundary = 3
    value = 0
  []
  [fix_disp_y]
    type = DirichletBC
    variable = disp_y
    boundary = 3
    value = 0
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
    base_name = uncracked
  []
  [./trial_stress]
    type = ADComputeFiniteStrainElasticStress
    base_name = uncracked
  [../]
  [expansion]
    type = ADComputeThermalExpansionEigenstrain
    temperature = T
    thermal_expansion_coeff = ${alpha_bulk}
    stress_free_temperature = ${T_ref}
    eigenstrain_name = thermal_expansion
    base_name = uncracked
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

  [./cracked_stress]
    type = ADComputePFFStress
    decomposition = spectral
    c = d
    E_name = E_el
    D_name = degradation
    use_current_history_variable = true
    uncracked_base_name = uncracked
    finite_strain_model = true
  [../]

  # Additional functional relaiton: conductivity and d
  [thermal_conductivity]
    type = ADParsedMaterial
    property_name = 'thermal_conductivity'
    expression = 'degradation*(k0-k0/50)+k0/50'
    coupled_variables = 'd '
    material_property_names = 'degradation'
    constant_names       = 'k0   '
    constant_expressions = '${k0}'
    outputs = exodus
  []
  [./gc] 
    type = ADParsedMaterial
    property_name = gc
    coupled_variables = 'T'
    expression = 'gc0 * (1 - b1 * ((T-T_ref)/T_max) + b2 * ((T-T_ref)/T_max)^2)'
    constant_names       = 'gc0    b1  b2  T_ref   T_max'
    constant_expressions = '${gc0} 1.8 1.1 ${T_ref} ${T_max}'
    outputs = exodus
  [../]
  [./alpha_critical] 
    type = ADParsedMaterial
    property_name = alpha_critical
    coupled_variables = 'T'
    expression = 'alpha_critical0 * (1 - b1 * ((T-T_ref)/T_max) + b2 * ((T-T_ref)/T_max)^2)'
    constant_names       = 'alpha_critical0    b1  b2  T_ref   T_max'
    constant_expressions = '${alpha_critical0} 1.8 1.1 ${T_ref} ${T_max}'
    outputs = exodus
  [../]
[]


[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [./crack_area]
    type = ElementIntegralVariablePostprocessor
    variable = d
  [../]
  [./max_d]
    type = NodalExtremeValue
    variable = d
  [../]
  [./top_stress_yy]
    type = SideAverageValue
    variable = stress_yy
    boundary = upper_left
  [../]
  [./max_accumulate]
    type = ElementExtremeValue
    variable = bar_alpha
  [../]
  [./cycle_counts]
    type = ElementAverageValue
    variable = n_cycle
  [../]
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package '
  petsc_options_value = 'lu       superlu_dist                  '
  #solve_type = PJFNK
  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'
  automatic_scaling = true

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-1
    optimal_iterations = 12
    cutback_factor = 0.3 
    growth_factor = 1.2
  [../]

  end_time = ${end_time}
  dtmin= 1e-7
  num_steps=1
  fixed_point_max_its = 6
  nl_max_its = 18  
  l_max_its = 24 
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-7
[]

[Outputs]
  file_base=themral_loading_dT700
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval = 4
[]
