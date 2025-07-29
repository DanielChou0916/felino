# Adjustments 0719
# Add linear increment time stage for adding BCs.
# But this is actually not ideal.
# Since if the material break during adding load stage it means either
# 1. BC is too extreme.
# 2. Material is too weak (inappropriate material parameters)

E = 21e4 #MPa= 210 GPa
nu = 0.3     #

bc_stress = 189 #MPa     #6530N/20 mm*2 mm (plane stress)

#umax = 0.05
period = 0.01
num_cycle = 32000
end_time = ${fparse period * num_cycle}
cycle_count_start = 0 # t= cycle_count_start  -> start to count cycle.
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[MultiApps]
  [crack]
    type = TransientMultiApp
    input_files = 'DB_f0.i'
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
  file = dogbone55.inp
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
[]

[AuxKernels]
  [./n_cycle_aux]
    type = FunctionAux
    variable = n_cycle
    function = current_cycle
  [../]
[]

[Functions]
  [./pulling_ramp]
    type = ParsedFunction
    expression = 'if(t < ${cycle_count_start}, ${bc_stress} * t / ${cycle_count_start}, ${bc_stress})'
  [../]

  [./current_cycle]
    type = ParsedFunction
    expression = 'if(t < ${cycle_count_start}, 0, (t - ${cycle_count_start}) / ${period})'
  [../]
[]



[BCs] # 2: top,   3:bottom
 #  [top_cycle]
 #   type = FunctionDirichletBC
 #   variable = 'disp_y'
 #   boundary = TOP_HOLE
 #   ## Instead of using cosine function for top displacement, fix it as constant
 #   ## The energy accumulation is taken by [./current_cycle] block in [Functions]
 #   #function = '${umax} * 0.5 * (cos(2 * 3.1415926 * t / ${period}) + 1)'
 #   function = '${umax}'
 # []

  [pull_y]   #duvergant stress is added into the hearder!
    type = FunctionNeumannBC
    function = pulling_ramp
    variable = disp_y
    boundary = UP
  []
  [yfix]
    type = DirichletBC
    variable = 'disp_y'
    boundary = BOTTOM
    value = 0
  []
  [xfix]
    type = DirichletBC
    variable = 'disp_x'
    boundary = BOTTOM
    value = 0
  []
[]


[Materials]
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
  #[./anisotropy]
  #  type = ADAnisotropicDirector
  #  normal = "-0.70710678 0.70710678 0"
  #  output_name = A_matrix
  #  outputs = exodus
  #[]
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
    boundary = UP
  [../]
  [./av_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = UP
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
    dt = 1e-2
    optimal_iterations = 12
    cutback_factor = 0.4
    growth_factor = 1.3
  [../]
  #dt = ${deltat}
  end_time = ${end_time}
  dtmin = 1e-8
  #num_steps=8
  fixed_point_max_its = 18
  nl_max_its = 16  
  l_max_its = 20  
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-7
[]

[Outputs]
  file_base=DB_0_det
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval = 3
[]