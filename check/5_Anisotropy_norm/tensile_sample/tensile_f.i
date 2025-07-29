E = 21e4 #MPa # 210 GPa
nu = 0.3     #
gc = 10     #KJ/m2 = MPa. mm
l = 0.1    #mm

xi = 1#0
C0 = 2.666667#2
L = 1e4

#alpha_critical = 1e10 #MPa
#R=1
#n=0.1
[Mesh]
  type = FileMesh
  file = tensile.inp
  #uniform_refine = 1
  #skip_partitioning = true
  construct_side_list_from_node_list=true
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Actions/PFNonconserved]
  [./d]
    free_energy = F
    kappa = kappa_op
    mobility = L
    variable_mobility=false
    use_automatic_differentiation = false
    use_grad_kappa = false
    use_anisotropic_matrix = true
    anisotropic_matrix = A_matrix
  [../]
[]

[AuxVariables]
  [./bounds_dummy]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]


[Materials]
  [./anisotropy]
    type = AnisotropicDirector
    input_type = xy_angle
    xy_angle_deg = 45
    normalize_director = det_norm
    coef = 100
    output_name = A_matrix
  []
  [./uncracked_strain]
    type = ComputeFiniteStrain
    base_name = uncracked
  [../]
  [./trial_stress]
    type = ComputeFiniteStrainElasticStress
    base_name = uncracked
  [../]
  [./public_materials_forPF_model]
    type = GenericConstantMaterial
    prop_names =  'gc     l    xi    C0      L  ' 
    prop_values = '${gc}  ${l} ${xi} ${C0}   ${L}' #'0 2'#for AT2 # Or use '1 2.6666667' for AT1
  [../]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor #Constitutive law here
    poissons_ratio = ${nu}
    youngs_modulus = ${E} #MPa
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
  [./local_fracture_energy] #Define psi_frac and alpha(d)
    type = DerivativeParsedMaterial
    property_name = local_fracture_energy
    coupled_variables = 'd'
    material_property_names = 'gc l xi C0'
    expression = '(xi*d+(1-xi)*d^2) * (gc / l)/C0'
    derivative_order = 2
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc l C0'
    property_name = kappa_op
    expression = '2 * gc * l / C0 '
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
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    #block = 1
    coupled_variables = d
    sum_materials = 'E_el local_fracture_energy'
    derivative_order = 2
    property_name = F
  [../]
[]


[Bounds]
  [./d_upper_bound]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1.0
  [../]
  [./d_lower_bound]
    type = VariableOldValueBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  [../]
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