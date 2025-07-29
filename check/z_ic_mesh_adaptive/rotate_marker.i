[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 20#160
    ny = 20#32
    nz = 10
    xmin = 40
    xmax = 60#
    ymax = 20 #
    zmax = 10
  []
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
      cx = 50
      cy = 10
      cz = 5
      lx = 3.8
      ly = 22
      lz = 20
      angle_z = 0
      angle_y = -47.5
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
    cx = '50'
    cy = '16'
    cz = '5'
    lx = '0.25'
    ly = '11'
    lz = '20'
    angle_z = '0'
    angle_y = '-45'
    inside = '1'
    outside = '0'
    int_width = '0.001 '
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
  file_base=bar_cut3
  exodus = true
  #csv = true
  print_linear_residuals = false
  time_step_interval = 1
[]
