#include "ADHydroVolStrainCoupling.h"

registerMooseObject("MooseApp", ADHydroVolStrainCoupling);

InputParameters
ADHydroVolStrainCoupling::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The operator of strain rate coupling term in porous elasticity.");
  params.addParam<MaterialPropertyName>("fluid_density", "fluid_density", "Name of fluid density as an input.");
  params.addParam<MaterialPropertyName>("biot_coefficient", "biot_coefficient", "Name of biot coefficient as an input.");
  params.addParam<MaterialPropertyName>("strain_rate", "strain_rate", "Name of strain rate as an input.");
  return params;
}

ADHydroVolStrainCoupling::ADHydroVolStrainCoupling(const InputParameters & parameters)
  : ADKernel(parameters),
  _fluid_density(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("fluid_density") )),
  _biot_coefficient(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("biot_coefficient") )),
  _strain_rate(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("strain_rate") ))
{
}

ADReal
ADHydroVolStrainCoupling::computeQpResidual()
{
  return _test[_i][_qp] * _fluid_density[_qp]*_biot_coefficient[_qp]*_strain_rate[_qp];
}
