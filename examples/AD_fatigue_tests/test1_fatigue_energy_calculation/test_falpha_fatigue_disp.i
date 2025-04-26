E = 21e4 #MPa= 210 GPa
nu = 0.3     #
gc = 2.7     #KJ/m2 = MPa.mm
l = 0.016    #mm

[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[MultiApps]
  [crack]
    type = TransientMultiApp
    input_files = 'test_falpha_fatigue_f.i'
  []
[]

[Transfers]
  [to_disp_x]
    type = MultiAppCopyTransfer #MultiAppGeometricInterpolationTransfer# MultiAppCopyTransfer #
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
  [from_fatigue]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'f_alpha'
    variable = 'f_alpha'
  []
  [from_current_fatigue]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'current_fatigue'
    variable = 'current_fatigue'
  []
  [from_accumulate_fatigue]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'bar_alpha'
    variable = 'bar_alpha'
  []
  [from_dfdx]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'dfatigue_dx'
    variable = 'dfatigue_dx'
  []
  [from_dfdy]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'dfatigue_dy'
    variable = 'dfatigue_dy'
  []
[]

[Mesh]
  file = mesh/mesh_in.e
  uniform_refine = 0
  skip_partitioning = true
  construct_side_list_from_node_list=true
[]

[Physics/SolidMechanics/QuasiStatic]
  [./All]
    add_variables = true
    strain = FINITE
    incremental = true
    additional_generate_output = 'stress_xx stress_yy stress_xy'
    use_automatic_differentiation=true
    strain_base_name = uncracked
    decomposition_method = EigenSolution
  [../]
[../]


[AuxVariables]
  [./d]
  []
  [./bounds_dummy]
  [../]
  [rf_x]
  []
  [rf_y]
  []
  [./f_alpha]
  order = CONSTANT
  family = MONOMIAL
  []
  [./current_fatigue]
    order = CONSTANT
    family = MONOMIAL
  []
  [./bar_alpha]
    order = FIRST
    family = MONOMIAL
  []
  [./dfatigue_dx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dfatigue_dy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    save_in = rf_x
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    save_in = rf_y
  []
[]


[BCs]
  [xdisp]
    type = FunctionDirichletBC
    variable = 'disp_x'
    boundary = 2
    function = 't'
    preset = false
  []
  [yfix]
    type = DirichletBC
    variable = 'disp_y'
    boundary = '2 3'
    value = 0
  []
  [xfix]
    type = DirichletBC
    variable = 'disp_x'
    boundary = 3
    value = 0
  []
[]


[Materials]
  [./pfbulkmat]
    type = ADGenericConstantMaterial
    prop_names =  'gc     l     '
    prop_values = '${gc}  ${l}  ' #Gc:MPa mm
  [../]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor #Constitutive law here
    poissons_ratio = ${nu}
    youngs_modulus = ${E} #MPa
    base_name = uncracked
  [../]
  [./trial_stress]
    type = ADComputeFiniteStrainElasticStress
    base_name = uncracked
  [../]
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
    c = d
    E_name = E_el
    D_name = degradation
    decomposition = spectral
    use_current_history_variable = true
    uncracked_base_name = uncracked
    finite_strain_model = true
  [../]
[]

[Postprocessors]
  [Fy]
    type = NodalSum
    variable = rf_y
    boundary = 2
  []
  [Fx]
    type = NodalSum
    variable = rf_x
    boundary = 2
  []
  [./top_stress_yy]
    type = SideAverageValue
    variable = stress_yy
    boundary = 2
  [../]
  [./av_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = 2
  [../]
  [./crack_area]
    type = ElementIntegralVariablePostprocessor
    variable = d
  [../]
  [./max_d]
    type = NodalExtremeValue
    variable = d
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./z_n_nl_its]
    type = NumNonlinearIterations
    accumulate_over_step = true
  [../]
  [./z_n_picard_its]
    type = NumFixedPointIterations
  [../]
  [./run_time]
    type = PerfGraphData
    data_type = TOTAL
    section_name = Root
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
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
    dt = 5e-5
    optimal_iterations = 12
    cutback_factor = 0.3 
    growth_factor = 1.25
  [../]

  end_time = 2e-2
  num_steps=2
  fixed_point_max_its = 30 
  nl_max_its = 16  
  l_max_its = 20  
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-7
[]

[Outputs]
  file_base=f_normal_spectral
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval =1
[]