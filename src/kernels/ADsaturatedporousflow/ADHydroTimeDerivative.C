#include "ADHydroTimeDerivative.h"

registerMooseObject("MooseApp", ADHydroTimeDerivative);

InputParameters
ADHydroTimeDerivative::validParams()
{
  InputParameters params = ADTimeKernelValue::validParams();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  params.addParam<MaterialPropertyName>("biot_modulus", "biot_modulus", "Name of biot modulus as an input.");
  return params;
}

ADHydroTimeDerivative::ADHydroTimeDerivative(const InputParameters & parameters)
  : ADTimeKernelValue(parameters),
  _biot_modulus(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("biot_modulus") ))
{
}

ADReal
ADHydroTimeDerivative::precomputeQpResidual()
{
  return _u_dot[_qp] / _biot_modulus[_qp];
}
