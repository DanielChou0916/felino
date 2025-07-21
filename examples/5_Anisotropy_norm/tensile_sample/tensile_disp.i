E = 21e4 #MPa= 210 GPa
nu = 0.3     #
gc = 300     #KJ/m2 = MPa.mm
l = 0.2   #mm

end_time = 1000 
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[MultiApps]
  [crack]
    type = TransientMultiApp
    input_files = 'tensile_f.i'
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
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'd'
    variable = 'd'
  []
[]

[Mesh]
  type = FileMesh
  file = tensile.inp
  #uniform_refine = 1
  #skip_partitioning = true
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
  [./resid_x]
  [../]
  [./resid_y]
  [../]
[]


[BCs] # 2: top,   3:bottom
  [top_cycle]
    type = FunctionDirichletBC
    variable = 'disp_y'
    boundary = top
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
  end_time = ${end_time}
  #num_steps=2
  fixed_point_max_its = 12
  nl_max_its = 25
  l_max_its = 40 
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-7
[]

[Outputs]
  file_base=none_45
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval = 10
[]