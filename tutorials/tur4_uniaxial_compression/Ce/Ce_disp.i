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
delta_elp_ce = 1.0      # 

# Other PFF parameters
#xi  = 1
#C0  = 2.666666667
l   = 0.2
#L   = 1e4


[GlobalParams]
  displacements = 'disp_x disp_y'
[]
[MultiApps]
  [crack]
    type = TransientMultiApp
    input_files = 'Ce_f.i'
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
[]

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



[Physics/SolidMechanics/QuasiStatic]
  [./All]
    add_variables = true
    strain = SMALL
    incremental = true
    additional_generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_yy strain_xy'
    use_automatic_differentiation=true
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

[Functions]
  [./tension]
    type = ParsedFunction
    expression ='t'
  [../]
  [./compression]
    type = ParsedFunction
    expression ='-t'
  [../]
[]

[BCs]
  [./bottomy]
    type = DirichletBC
    boundary = 4#=set id from paraview
    variable = disp_y
    value = 0.0
  [../]
  [./bottomx]
    type = DirichletBC
    boundary = fix_node# =set id from paraview
    variable = disp_x
    value = 0.0
  [../]
  [./pull]
    type = FunctionDirichletBC
    boundary = 3#'NODE_YMAX_SET'
    variable = disp_y
    function = compression
  [../]
[]


[Materials]
  [./PFF_AGG]
    type = ADGenericConstantMaterial
    block = 0
    prop_names =  'gc        l     sigma_ts        sigma_cs         delta_elp' 
    prop_values = '${gc_agg} ${l}  ${sigma_ts_agg} ${sigma_cs_agg}  ${delta_elp_agg}'
  [../]
  [./PFF_CE]
    type = ADGenericConstantMaterial
    block = 1
    prop_names =  'gc        l     sigma_ts        sigma_cs         delta_elp' 
    prop_values = '${gc_ce} ${l}   ${sigma_ts_ce}   ${sigma_cs_ce}     ${delta_elp_ce}'
  [../]
  #[./public_materials_PFF]
  #  type = ADGenericConstantMaterial
  #  block = '0 1'
  #  prop_names =  'xi    C0      L     ' 
  #  prop_values = '${xi} ${C0}   ${L}  '
  #[../]
  [./elasticity_tensorAgg]
    type = ADComputeIsotropicElasticityTensor 
    block = 0
    poissons_ratio = ${nu_agg}
    youngs_modulus = ${E_agg} #MPa
  [../]
  [./elasticity_tensorCE]
    type = ADComputeIsotropicElasticityTensor 
    block = 1
    poissons_ratio = ${nu_ce}
    youngs_modulus = ${E_ce} #MPa
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
    outputs = exodus
  [] 
[]

[Postprocessors]
  [Fy]
    type = NodalSum
    variable = rf_y
    boundary = 3
  []
  [Fx]
    type = NodalSum
    variable = rf_x
    boundary = 3
  []
  [./top_stress_yy]
    type = SideAverageValue
    variable = stress_yy
    boundary = 3
  [../]
  [./av_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = 3
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

  #[./TimeStepper]
  #  type = IterationAdaptiveDT
  #  dt = 5e-3
  #  optimal_iterations = 12
  #  cutback_factor = 0.3 
  #  growth_factor = 1.3
  #[../]
  dt = 4e-2
  dtmin=1e-10
  end_time = 0.5
  #num_steps = 4
  fixed_point_max_its = 25
  nl_max_its = 20  
  l_max_its = 24 
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-7
[]

[Outputs]
  file_base=Ce
  exodus = true
  #perf_graph = true
  csv = true
  time_step_interval =1
[]