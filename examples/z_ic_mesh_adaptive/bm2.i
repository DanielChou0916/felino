[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 12
    ny = 24
    ymax = 2.0
    xmax = 1.0
  []
[]

[Adaptivity]
  marker = marker
  initial_marker = marker
  initial_steps = 4
  stop_time = 0
  max_h_level = 4
  [Markers]
    [marker]
      type = RotatedBoxMarker
      cx = 0.5
      cy = 1.0
      cz = 0.0
      lx = 2
      ly = 0.4
      lz = 0.0001
      angle_z = 0
      angle_y = 0
      angle_x = 0
      inside = REFINE
      outside = DO_NOTHING
    []
  []
[]

[ICs]
  [./boxes]
    type = MultiRotBoundingBoxIC
    variable = d

    cx = '0.2 0.8'
    cy = '0.95 1.05'
    cz = '0 0'
    lx = '0.402 0.402'
    ly = '0.01 0.01'
    lz = '0.0 0'
    angle_x = '0 0'
    angle_y = '0 0'
    angle_z = '0 0'
    inside = '1 1'
    outside = '0.0 0.0'
    int_width = '0 0'
  [../]
[]

[Variables]
  [d]
  []
[]
[Kernels]
  [diff]
    type = ADDiffusion
    variable = d
  []
[]

#[Materials]
#  [./phase_normal]
#    type = PhaseNormalTensor
#    phase = d
#    normal_tensor_name = dir_tensor
#    outputs = exodus
#  [../]
#[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '
  automatic_scaling = true

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  dt = 1e-5
  end_time = 0#0.001#6.5e-3

  fixed_point_max_its = 20
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-8
[]

[Outputs]
  file_base=bm2
  exodus = true
  #csv = true
  print_linear_residuals = false
  time_step_interval = 1
[]
