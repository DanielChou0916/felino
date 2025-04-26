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
  params.addCoupledVar("grad_kappa_x", "Gradient of kappa in x direction");
  params.addCoupledVar("grad_kappa_y", "Gradient of kappa in y direction");
  params.addCoupledVar("grad_kappa_z", "Gradient of kappa in z direction (only for 3D)");                      
  return params;
}

ADACInterfaceGradKappa::ADACInterfaceGradKappa(const InputParameters & parameters)
  : ADKernel(parameters),
    _prop_L(getADMaterialProperty<Real>("mob_name")),
    _name_L(getParam<MaterialPropertyName>("mob_name")),
    _kappa(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("kappa_name"))),
    _grad_kappa_x(coupledValue("grad_kappa_x")),
    _grad_kappa_y(coupledValue("grad_kappa_y")),
    _grad_kappa_z(parameters.isParamValid("grad_kappa_z") ? &coupledValue("grad_kappa_z") : nullptr),  
    _variable_L(getParam<bool>("variable_L")),
    _dLdop(_variable_L
               ? &getADMaterialProperty<Real>(derivativePropertyNameFirst(_name_L, _var.name()))
               : nullptr),
    _nvar(Coupleable::_coupled_standard_moose_vars.size()),
    _dLdarg(_nvar),
    _gradarg(_nvar)
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
}

ADReal
ADACInterfaceGradKappa::computeQpResidual()
{
  ADRealVectorValue nabla_Lpsi = _prop_L[_qp] * _grad_test[_i][_qp];

  if (_variable_L)
  {
    ADRealVectorValue grad_L = _grad_u[_qp] * (*_dLdop)[_qp];
    for (unsigned int i = 0; i < _nvar; ++i)
      grad_L += (*_gradarg[i])[_qp] * (*_dLdarg[i])[_qp];

    nabla_Lpsi += grad_L * _test[_i][_qp];
  }

  ADReal residual = _grad_u[_qp] * _kappa[_qp] * nabla_Lpsi;

  ADRealVectorValue grad_kappa_vec;

  grad_kappa_vec(0) = _grad_kappa_x[_qp];
  
  grad_kappa_vec(1) = _grad_kappa_y[_qp];
  if (LIBMESH_DIM == 3 && _grad_kappa_z)
    grad_kappa_vec(2) = (*_grad_kappa_z)[_qp];
  
  // Add grad(kappa) term
  residual -= (grad_kappa_vec * _grad_u[_qp]) * _prop_L[_qp] * _test[_i][_qp];

  return residual;
}
