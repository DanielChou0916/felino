E = 21e4 #MPa= 210 GPa
nu = 0.3     #
#gc = 2.7     #KJ/m2 = MPa.mm

umax = 0.005
period = 0.001
num_cycle = 15000
end_time = ${fparse period * num_cycle}
#deltat = ${fparse 100 * period} 
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[MultiApps]
  [crack]
    type = TransientMultiApp
    input_files = 'AD_3Dcut_f.i'
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
  [to_disp_z]
  type = MultiAppCopyTransfer #MultiAppGeometricInterpolationTransfer# MultiAppCopyTransfer #
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
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10#160
    ny = 10#32
    nz = 5
    xmax = 20#
    ymax = 20 #
    zmax = 10
  []
  [fix_node1]
      type = BoundingBoxNodeSetGenerator
      input = gen
      bottom_left = '-0.1 -0.1 -0.1'
      top_right = '0.0001 0.0001 0.0001'
      new_boundary = 'fix_node1'
  []
  [fix_node2]
      type = BoundingBoxNodeSetGenerator
      input = fix_node1
      bottom_left = '19.999 19.999 9.999'
      top_right = '20.0001 20.0001 20.0001'
      new_boundary = 'fix_node2'
  []
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  initial_steps = 3
  stop_time = 0
  max_h_level = 3
  [Markers]
    [marker]
      type = RotatedBoxMarker
      cx = 10
      cy = 10
      cz = 5
      lx = 4.1
      ly = 22
      lz = 20
      angle_z = 0
      angle_y = -38
      angle_x = 0
      inside = REFINE
      outside = DO_NOTHING
    []
  []
[]

[ICs]
  [init_d_box]
    type = MultiRotBoundingBoxIC
    variable = d
    cx = '10'
    cy = '16'
    cz = '5'
    lx = '0.24'
    ly = '11'
    lz = '20'
    angle_z = '0'
    angle_y = '-45'
    inside = '1'
    outside = '0'
    int_width = '0.001 '
  [../]
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

[BCs] # 2: top,   3:bottom
  [right_cycle]
    type = FunctionDirichletBC
    variable = 'disp_x'
    boundary = 'right'
    ## It is not recommanded to set BC as explicit periodic function
    ## this well  expose the simulation in no fatigue accumulation risk!
    #function = '${umax} * 0.5 * (cos(2 * 3.1415926 * t / ${period}) + 1)'
    function = '${umax}'
  []
  [left_cycle]
    type = FunctionDirichletBC
    variable = 'disp_x'
    boundary = 'left'
    ## It is not recommanded to set BC as explicit periodic function
    ## this well  expose the simulation in no fatigue accumulation risk!
    #function = '-1*${umax} * 0.5 * (cos(2 * 3.1415926 * t / ${period}) + 1)'
    function = '-1*${umax}'
  []
  [yfix12]
    type = DirichletBC
    variable = 'disp_y'
    boundary = 'fix_node1'
    value = 0
  []
  [zfix12]
    type = DirichletBC
    variable = 'disp_z'
    boundary = 'fix_node2'
    value = 0
  []
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
  [./side_stress_xx]
    type = SideAverageValue
    variable = stress_xx
    boundary = right
  [../]
  [./av_disp_x]
    type = SideAverageValue
    variable = disp_x
    boundary = right
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
  fixed_point_max_its = 15
  nl_max_its = 15
  l_max_its = 18
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-7
[]

[Outputs]
  file_base=NA_bar_3D_jump2
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval = 2
[]