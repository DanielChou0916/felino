#include "ADHydroTimeDerivative.h"

registerMooseObject("MooseApp", ADHydroTimeDerivative);

InputParameters
ADHydroTimeDerivative::validParams()
{
  InputParameters params = ADTimeKernelValue::validParams();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  params.addParam<MaterialPropertyName>("fluid_density", "fluid_density", "Name of fluid density as an input.");
  params.addParam<MaterialPropertyName>("biot_modulus", "biot_modulus", "Name of biot modulus as an input.");
  return params;
}

ADHydroTimeDerivative::ADHydroTimeDerivative(const InputParameters & parameters)
  : ADTimeKernelValue(parameters),
  _fluid_density(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("fluid_density") )),
  _biot_modulus(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("biot_modulus") ))
{
}

ADReal
ADHydroTimeDerivative::precomputeQpResidual()
{
  return _fluid_density[_qp] *_u_dot[_qp] / _biot_modulus[_qp];
}
