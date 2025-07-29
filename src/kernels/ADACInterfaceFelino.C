#include "ADACInterfaceFelino.h"

registerMooseObject("PhaseFieldApp", ADACInterfaceFelino);

InputParameters
ADACInterfaceFelino::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Extended ACInterface Kernel including grad(kappa) term.");
  params.addParam<MaterialPropertyName>("mob_name", "L", "Mobility used in the kernel");
  params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "Interfacial parameter kappa");
  params.addCoupledVar("args", "Coupled nonlinear variables affecting mobility");
  params.addParam<bool>("variable_L",
                        true,
                        "Set to false if L is constant over the domain.");
  params.addParam<bool>("use_anisotropic_matrix",
                        false,
                        "Set to false for isotropic PFF");  
  params.addParam<MaterialPropertyName>("anisotropic_matrix_name", "A", "The anisotropic matrix used with the kernel");                                    
  return params;
}

ADACInterfaceFelino::ADACInterfaceFelino(const InputParameters & parameters)
  : ADKernel(parameters),
    _prop_L(getADMaterialProperty<Real>("mob_name")),
    _name_L(getParam<MaterialPropertyName>("mob_name")),
    _kappa(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("kappa_name"))),
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
ADACInterfaceFelino::computeQpResidual()
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

  return residual;
}
