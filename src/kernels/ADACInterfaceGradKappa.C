#include "ADACInterfaceGradKappa.h"

registerMooseObject("PhaseFieldApp", ADACInterfaceGradKappa);

InputParameters
ADACInterfaceGradKappa::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Extended ACInterface Kernel including grad(kappa) term.");
  params.addParam<MaterialPropertyName>("mob_name", "L", "Mobility used in the kernel");
  params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "Interfacial parameter kappa");
  params.addCoupledVar("args", "Coupled nonlinear variables affecting mobility");
  params.addParam<bool>("variable_L",
                        true,
                        "Set to false if L is constant over the domain.");
  params.addParam<bool>("use_grad_kappa",
                        false,
                        "Set to false if L is constant over the domain.");
  params.addCoupledVar("grad_kappa_x", "Gradient of kappa in x direction");
  params.addCoupledVar("grad_kappa_y", "Gradient of kappa in y direction");
  params.addCoupledVar("grad_kappa_z", "Gradient of kappa in z direction (only for 3D)");  
  params.addParam<bool>("use_anisotropic_matrix",
                        false,
                        "Set to false for isotropic PFF");  
  params.addParam<MaterialPropertyName>("anisotropic_matrix_name", "A", "The anisotropic matrix used with the kernel");                                    
  return params;
}

ADACInterfaceGradKappa::ADACInterfaceGradKappa(const InputParameters & parameters)
  : ADKernel(parameters),
    _prop_L(getADMaterialProperty<Real>("mob_name")),
    _name_L(getParam<MaterialPropertyName>("mob_name")),
    _kappa(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("kappa_name"))),
    _use_grad_kappa(getParam<bool>("use_grad_kappa")),
    _grad_kappa_x((_use_grad_kappa && parameters.isParamValid("grad_kappa_x"))
                  ? &coupledValue("grad_kappa_x"): nullptr),
    _grad_kappa_y((_use_grad_kappa && parameters.isParamValid("grad_kappa_y"))
                  ? &coupledValue("grad_kappa_y"): nullptr),
    _grad_kappa_z((_use_grad_kappa && parameters.isParamValid("grad_kappa_z"))
                  ? &coupledValue("grad_kappa_z"): nullptr),
    _variable_L(getParam<bool>("variable_L")),
    _dLdop(_variable_L
               ? &getADMaterialProperty<Real>(derivativePropertyNameFirst(_name_L, _var.name()))
               : nullptr),
    _nvar(Coupleable::_coupled_standard_moose_vars.size()),
    _dLdarg(_nvar),
    _gradarg(_nvar),
    _use_anisotropic_matrix(getParam<bool>("use_anisotropic_matrix")),// ✅ 2025/07/06
    _A_ptr((nullptr))// ✅ 2025/07/06
{
  if (_variable_L)
    for (unsigned int i = 0; i < _nvar; ++i)
    {
      MooseVariable * ivar = _coupled_standard_moose_vars[i];
      const VariableName iname = ivar->name();

      if (iname == _var.name())
        paramError("args", "The kernel variable should not be included in coupled 'args'.");

      _dLdarg[i] = &getADMaterialProperty<Real>(derivativePropertyNameFirst(_name_L, iname));
      _gradarg[i] = &(ivar->adGradSln());
    }
// constructor
  if (_use_grad_kappa && (!parameters.isParamValid("grad_kappa_x") || !parameters.isParamValid("grad_kappa_y")))
  {
    if (processor_id() == 0)
      mooseWarning("use_grad_kappa is true, but grad_kappa_x or grad_kappa_y not provided → skipped");
  }
  //anisotropic matrix A
  if (_use_anisotropic_matrix)// ✅ 2025/07/06
  { 
    const std::string & name = getParam<MaterialPropertyName>("anisotropic_matrix_name");
    if (hasADMaterialProperty<RankTwoTensor>(name))
      {
        mooseInfo("Anisotropic matrix loaded and applied in residual");
        _A_ptr = &getADMaterialProperty<RankTwoTensor>(name);
      }
    else if (processor_id() == 0)
      mooseWarning("Anisotropic matrix material '", name, "' not found. Falling back to identity.");
  }
  else
    mooseInfo("Use isotropic PFF!");
}

ADReal
ADACInterfaceGradKappa::computeQpResidual()
{
  ADRankTwoTensor A = (_use_anisotropic_matrix && _A_ptr) ? (*_A_ptr)[_qp] : ADRankTwoTensor::initIdentity;
  ADRealVectorValue nabla_Lpsi = _prop_L[_qp] * _grad_test[_i][_qp];

  if (_variable_L)
  {
    ADRealVectorValue grad_L = _grad_u[_qp] * (*_dLdop)[_qp];
    for (unsigned int i = 0; i < _nvar; ++i)
      grad_L += (*_gradarg[i])[_qp] * (*_dLdarg[i])[_qp];

    nabla_Lpsi += grad_L * _test[_i][_qp];
  }

  ADReal residual = (A * _grad_u[_qp]) * _kappa[_qp] * nabla_Lpsi;

  //ADRealVectorValue grad_kappa_vec;

  if (_use_grad_kappa)
  {
    mooseInfo("use_grad_kappa = true");
  
    if (_grad_kappa_x && _grad_kappa_y)
    {
      ADRealVectorValue grad_kappa_vec;
      grad_kappa_vec(0) = (*_grad_kappa_x)[_qp];
      grad_kappa_vec(1) = (*_grad_kappa_y)[_qp];
  
      if (LIBMESH_DIM == 3 && _grad_kappa_z)
        grad_kappa_vec(2) = (*_grad_kappa_z)[_qp];
  
      mooseInfo("grad_kappa loaded and applied in residual");
      residual -= (grad_kappa_vec * (A * _grad_u[_qp])) * _prop_L[_qp] * _test[_i][_qp];
    }
    //else
    //{
    //  mooseWarning("use_grad_kappa is true, but grad_kappa_x or grad_kappa_y not provided → skipped");
    //}
  }
  return residual;
}
