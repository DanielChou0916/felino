E = 21e4 #MPa= 210 GPa
nu = 0.3     #
gc = 300     #KJ/m2 = MPa.mm
l = 0.2   #mm

#umax = 0.0083
#bc_stress = 50


period = 0.0001
num_cycle = 10000000
end_time = ${fparse period * num_cycle}
deltat = ${fparse 1000 * period} 
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[MultiApps]
  [crack]
    type = TransientMultiApp
    input_files = 'NA_f.i'
  []
[]

[Transfers]
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
  [to_CLA]
    type = MultiAppCopyTransfer
    to_multi_app = 'crack'
    source_variable = 'n_cycle'
    variable = 'n_cycle'
  []
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'd'
    variable = 'd'
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
  [from_fatigue_function]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'f_alpha'
    variable = 'f_alpha'
  []
  [from_kappa]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'kappa_op'
    variable = 'kappa_op'
  []
  [from_dkdx]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'dkappa_dx'
    variable = 'dkappa_dx'
  []
  [from_dkdy]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'dkappa_dy'
    variable = 'dkappa_dy'
  []
[]

[Mesh]
  type = FileMesh
  file = tensile.inp
#  uniform_refine = 1
#  skip_partitioning = true
  construct_side_list_from_node_list=true
[]

[Physics/SolidMechanics/QuasiStatic]
  [./All]
    add_variables = true
    strain = FINITE
    incremental = true
    additional_generate_output = 'stress_xx stress_yy stress_xy'
    use_automatic_differentiation=false
    strain_base_name = uncracked
    decomposition_method = EigenSolution
    save_in = 'resid_x resid_y'
  [../]
[../]


[AuxVariables]
  [./d]
  []
  [./bounds_dummy]
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
  [./dkappa_dx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dkappa_dy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./n_cycle]
    order = CONSTANT
    family = MONOMIAL
  [../]
    [./resid_x]
  [../]
  [./resid_y]
  [../]
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

[BCs] # 2: top,   3:bottom
  [top_cycle]
    type = FunctionDirichletBC
    variable = 'disp_y'
    boundary = top
    ## Instead of using cosine function for top displacement, fix it as constant
    ## The energy accumulation is taken by [./current_cycle] block in [Functions]
    #function = '${umax} * 0.5 * (cos(2 * 3.1415926 * t / ${period}) + 1)'
    #function = '${umax}'
    function ='0.5 / 60 * t'
  []

 # [pull_y2]
  #  type = NeumannBC
  #  value = ${bc_stress}
  #  #type = DirichletBC
  #  #value = 0.001
  #  variable = disp_y
  #  boundary = 'top'
  #[]



  [yfix]
    type = DirichletBC
    variable = 'disp_y'
    boundary = bottom
    value = 0
  []
  [xfix]
    type = DirichletBC
    variable = 'disp_x'
    boundary = bottom
    value = 0
  []
[]


[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names =  'gc     l     '
    prop_values = '${gc}  ${l}  ' #Gc:MPa mm
  [../]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor #Constitutive law here
    poissons_ratio = ${nu}
    youngs_modulus = ${E} #MPa
    base_name = uncracked
  [../]
  [./trial_stress]
    type = ComputeFiniteStrainElasticStress
    base_name = uncracked
  [../]
  [./degradation] # Define w(d)
    type = DerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'd'
    expression = '(1-d)^p*(1-k)+k'
    constant_names       = 'p k'
    constant_expressions = '2 1e-6'
    derivative_order = 2
  [../]
  [./cracked_stress]
    type = ComputePFFStress
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
  [./cycle_current]
    type = ElementAverageValue
    variable = n_cycle
  [../]
  [./max_current]
    type = ElementExtremeValue
    variable = current_fatigue
  [../]
  [./max_accumulate]
    type = ElementExtremeValue
    variable = bar_alpha
  [../]
  [./top_stress_yy]
    type = SideAverageValue
    variable = stress_yy
    boundary = top
  [../]
  [./av_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  [../]
  [./crack_area]
    type = ElementIntegralVariablePostprocessor
    variable = d
  [../]
  [./max_d]
    type = NodalExtremeValue
    variable = d
  [../]
  #[./dt]
  #  type = TimestepSize
  #[../]
  #[./z_n_nl_its]
  #  type = NumNonlinearIterations
  #  accumulate_over_step = true
  #[../]
  #[./z_n_picard_its]
  #  type = NumFixedPointIterations
  #[../]
  [./run_time]
    type = PerfGraphData
    data_type = TOTAL
    section_name = Root
  [../]
    [./resid_x]
    type = NodalSum
    variable = resid_x
    boundary = top
  [../]
  [./resid_y]
    type = NodalSum
    variable = resid_y
    boundary = top
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
  dt = ${deltat}
  end_time = ${end_time}
  #num_steps=2
  fixed_point_max_its = 12
  nl_max_its = 16  
  l_max_its = 20  
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-7
[]

[Outputs]
  file_base=NA_weak_45_coef10_2
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval = 10
[]