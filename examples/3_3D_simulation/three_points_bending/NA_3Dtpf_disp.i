E = 21e4 #MPa= 210 GPa
nu = 0.3     #
#gc = 2.7     #KJ/m2 = MPa.mm

umax = 0.01
period = 0.0001
num_cycle = 900000
end_time = ${fparse period * num_cycle}
#deltat = ${fparse 100 * period} 
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[MultiApps]
  [crack]
    type = TransientMultiApp
    input_files = 'AD_3Dtpf_f.i'
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
  [to_disp_z]
  type = MultiAppCopyTransfer 
  to_multi_app = 'crack'
  source_variable = 'disp_z'
  variable = 'disp_z'                                                                     
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
  #[from_dkdx]
  #  type = MultiAppCopyTransfer
  #  from_multi_app = 'crack'
  #  source_variable = 'dkappa_dx'
  #  variable = 'dkappa_dx'
  #[]
  #[from_dkdy]
  #  type = MultiAppCopyTransfer
  #  from_multi_app = 'crack'
  #  source_variable = 'dkappa_dy'
  #  variable = 'dkappa_dy'
  #[]
  #[from_dkdz]
  #type = MultiAppCopyTransfer
  #from_multi_app = 'crack'
  #source_variable = 'dkappa_dz'
  #variable = 'dkappa_dz'
  #[]
[]

[Mesh]
  [gen]
    type = FileMeshGenerator
    file = ../../mesh_files/three_points_metal.inp
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


[Physics/SolidMechanics/QuasiStatic]
  [./All]
    add_variables = true
    strain = FINITE
    incremental = true
    additional_generate_output = 'stress_xx stress_yy stress_xy stress_zz stress_xz stress_yz'
    use_automatic_differentiation=false
    strain_base_name = uncracked
    decomposition_method = EigenSolution
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
  #[./dkappa_dx]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./dkappa_dy]
  #  order = CONSTANT
  #  family = MONOMIAL
  #[../]
  #[./dkappa_dz]
  #order = CONSTANT
  #family = MONOMIAL
  #[../]
  [./n_cycle]
    order = CONSTANT
    family = MONOMIAL
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

[BCs]
  [./BOTTOMy]
    type = DirichletBC
    boundary = BACK_SUP2
    variable = disp_y
    value = 0.0
  [../]
  [./BOTTOMx]
    type = DirichletBC
    boundary = BACK_SUP1
    variable = disp_x
    value = 0.0
  [../]
  [./BOTTOMz]
    type = DirichletBC
    boundary = 'BACK_SUP1 BACK_SUP2'
    variable = disp_z
    value = 0.0
  [../]
  [./FRONT_SUP_load]
    type = FunctionDirichletBC
    boundary = FRONT_SUP
    variable = disp_z
    function = ${fparse -1*umax}
  [../]
  #[./FRONT_SUPx]
  #  type = DirichletBC
  #  boundary = FRONT_SUP
  #  variable = disp_x
  #  value = 0.0
  #[../]
  #[./FRONT_SUPy]
  #  type = DirichletBC
  #  boundary = FRONT_SUP
  #  variable = disp_y
  #  value = 0.0
  #[../]
[]


[Materials]
  #[./pfbulkmat]
  #  type = GenericConstantMaterial
  #  prop_names =  'gc     l     '
  #  prop_values = '${gc}  ${l}  ' #Gc:MPa mm
  #[../]
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
  [./side_stress_zz]
    type = SideAverageValue
    variable = stress_zz
    boundary = FRONT_SUP
  [../]
  [./av_disp_z]
    type = AverageNodalVariableValue
    variable = disp_z
    boundary = FRONT_SUP#FRONT_SUP
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
[]


[Preconditioning]
  [./smp]
    type = SMP
    full = true
 [../]
[]


[Executioner]
  type = Transient

  #solve_type = NEWTON
  #petsc_options_iname = '-pc_type -pc_factor_mat_solver_package '
  #petsc_options_value = 'lu       superlu_dist                  '
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  automatic_scaling = true

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-7

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    optimal_iterations = 15
    cutback_factor = 0.8 
    growth_factor = 1.3
  [../]
  #dt = ${deltat}
  end_time = ${end_time}
  num_steps=1
  fixed_point_max_its = 8
  nl_max_its = 30
  l_max_its = 25
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-7
[]

[Outputs]
  file_base=NA_tp3D
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval = 5
[]