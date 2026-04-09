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
delta_elp_ce = 0.64     # 

# Other PFF parameters
#xi  = 1
#C0  = 2.666666667
l   = 1.5
#L   = 1e4


[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[MultiApps]
  [crack]
    type = TransientMultiApp
    input_files = 'AD_f.i'
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
    type = MultiAppCopyTransfer
    to_multi_app = 'crack'
    source_variable = 'disp_z'
    variable = 'disp_z'
  []
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'crack'
    source_variable = 'd'
    variable = 'd'
  []
[]

[Mesh]
  [gen]
    type = FileMeshGenerator
    file = ../../mesh_files/cylinder_200_2.inp
  []
  construct_side_list_from_node_list=true
[]



[Physics/SolidMechanics/QuasiStatic]
  [./All]
    add_variables = true
    strain = SMALL
    incremental = true
    #planar_formulation = PLANE_STRAIN
    additional_generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz'
  [../]
[../]


[AuxVariables]
  [./d]
  []
  [./bounds_dummy]
  [../]
  #[rf_x]
  #[]
  #[rf_y]
  #[]
  #[rf_z]
  #[]
[]

#[Kernels]
#  [solid_x]
#    type = StressDivergenceTensors
#    variable = disp_x
#    component = 0
#    save_in = rf_x
#  []
#  [solid_y]
#    type = StressDivergenceTensors
#    variable = disp_y
#    component = 1
#    save_in = rf_y
#  []
#  [solid_z]
#    type = StressDivergenceTensors
#    variable = disp_z
#    component = 2
#    save_in = rf_z
#  []
#[]

[Functions]
  [./tension]
    type = ParsedFunction
    expression ='t'
  [../]
  [./compression]
    type = ParsedFunction
    expression ='-t'
  [../]
  [./dts]
    type = PiecewiseLinear
    y = '1e-1 1e-2 5e-3  1e-3 5e-3 1e-2  5e-2 1e-1'
    x = '0    0.4  0.7   1.0  2  2.8   3.2    4.0'
  [../]
[]


[BCs]
  [./BOTTOMy]
    type = DirichletBC
    boundary = BOTTOM
    variable = disp_y
    value = 0.0
  [../]
  [./BOTTOMx]
    type = DirichletBC
    boundary = BOTTOM
    variable = disp_x
    value = 0.0
  [../]
  [./BOTTOMz]
    type = DirichletBC
    boundary = BOTTOM
    variable = disp_z
    value = 0.0
  [../]
  [./TOP_load]
    type = FunctionDirichletBC
    boundary = TOP
    variable = disp_y
    function = compression
  [../]
  [./TOPx]
    type = DirichletBC
    boundary = TOP
    variable = disp_x
    value = 0.0
  [../]
  [./TOPz]
    type = DirichletBC
    boundary = TOP
    variable = disp_z
    value = 0.0
  [../]
  #[./left_compression]
  #  type = NeumannBC
  #  boundary = 2
  #  variable = disp_x
  #  value = 5 #MPa
  #[../]
  #[./right_compression]
  #  type = NeumannBC
  #  boundary = 1
  #  variable = disp_x
  #  value = -5 #MPa
  #[../]
[]


[Materials]
#Apply d in cement only delete below block:
  [./PFF_AGG]
    type = GenericConstantMaterial
    block = 0
    prop_names =  'gc        l     sigma_ts        sigma_cs         delta_elp' 
    prop_values = '${gc_agg} ${l}  ${sigma_ts_agg} ${sigma_cs_agg}  ${delta_elp_agg}'
  [../]
  [./PFF_CE]
    type = GenericConstantMaterial
    block = 1
    prop_names =  'gc        l     sigma_ts        sigma_cs         delta_elp' 
    prop_values = '${gc_ce} ${l}   ${sigma_ts_ce}   ${sigma_cs_ce}     ${delta_elp_ce}'
  [../]
  [./elasticity_tensorAgg]
    type = ComputeIsotropicElasticityTensor #Constitutive law here
    block = 0
    poissons_ratio = ${nu_agg}
    youngs_modulus = ${E_agg} #MPa
  [../]
  [./elasticity_tensorCE]
    type = ComputeIsotropicElasticityTensor #Constitutive law here
    block = 1
    poissons_ratio = ${nu_ce}
    youngs_modulus = ${E_ce} #MPa
  [../]
  [./degradation] 
    type = DerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'd'
    expression = '(1-d)^p*(1-k)+k'
    constant_names       = 'p k'
    constant_expressions = '2 1e-6 '
    derivative_order = 2
  [../]
  [./additional_driving_force]
    type = ComputeExtraDrivingForce
    D_name = degradation 
    Ce_name = Ce
    outputs = exodus
  [] 
  [./cracked_stress] 
    type = LinearElasticPFFractureStress
    c = d
    E_name = E_el
    D_name = degradation 
    decomposition_type = none
    use_current_history_variable = true
    use_snes_vi_solver = true
  [../]
[]

[Postprocessors]
  #[Fy]
  #  type = NodalSum
  #  variable = rf_y
  #  boundary = TOP#TOP
  #[]
  #[Fx]
  #  type = NodalSum
  #  variable = rf_x
  #  boundary = TOP#TOP
  #[]
  [./TOP_stress_yy]
    type = SideAverageValue
    variable = stress_yy
    boundary = TOP#TOP
  [../]
  [./av_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = TOP#TOP
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

  #solve_type = NEWTON
  #petsc_options_iname = '-pc_type -pc_factor_mat_solver_package '
  #petsc_options_value = 'lu       superlu_dist                  '
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  automatic_scaling = true

  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.002
    optimal_iterations = 12
    cutback_factor = 0.3 
    growth_factor = 1.1
  [../]
  dtmax = 0.002
  #[./TimeStepper]
  #  type = FunctionDT
  #  function = dts
  #[../]
  dtmin = 1e-8
  end_time = 0.1
  fixed_point_max_its = 18
  nl_max_its = 35
  l_max_its = 40
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-7
  fixed_point_abs_tol = 1e-8
[]

[Outputs]
  file_base=bz_3d
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval = 1
[]