########Temp System parameters############
rho = 7.85e-6 #kg/mm^3 #7850 kg/m^3
c = 460 #J/(kg·K)
k0 = 0.5 #(J/10s)/mmK=50 (J/s)/mK
dT_bc = 700
T_ref = 300
########Temp System parameters############
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
  #uniform_refine = 1
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  initial_steps = 1
  stop_time = 0
  max_h_level = 1
  [Markers]
    [marker]
      type = RotatedBoxMarker
      cx = 0.25
      cy = 0.25
      cz = 0
      lx = 0.51
      ly = 0.51
      lz = 0
      angle_z = 0
      angle_y = 0
      angle_x = 0
      inside = REFINE
      outside = DO_NOTHING
    []
  []
[]

[Variables]
  [T]
    initial_condition = ${T_ref}
  []
[]

[AuxVariables]
  [./d]
  []
[]

[Kernels]
  [heat_conduction]
    type = ADHeatConduction
    variable = T
  []
  [time_derivative]
    type = ADHeatConductionTimeDerivative
    variable = T
  []
  #[heat_source]
  #  type = ADBodyForce
  #  variable = T
  #  value = 5e4
  #[]
[]

[BCs]
  [t_load]
    ## (t - 5): Start to heating up the top when t=5, about 25 cycles
    type = FunctionDirichletBC
    variable = T
    function = '${T_ref} + ${dT_bc} / (1 + exp(-3 * (t - 3)))'
    boundary = 2
  []
  [fix_T]
    type = DirichletBC
    variable = T
    boundary = 3
    value = ${T_ref}
  []
[]

[Materials]
  [thermal_materials]
    type = ADGenericConstantMaterial
    prop_names =  'density  k0     specific_heat'
    prop_values = '${rho}   ${k0}  ${c}        '
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
  # Additional functional relaiton: conductivity and d
  [thermal_conductivity]
    type = ADParsedMaterial
    property_name = 'thermal_conductivity'
    expression = 'degradation*(k0-k0/50)+k0/50'
    coupled_variables = 'd '
    material_property_names = 'degradation'
    constant_names       = 'k0   '
    constant_expressions = '${k0}'
  []
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
