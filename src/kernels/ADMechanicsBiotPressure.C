// ADMechanicsBiotPressure.C
#include "ADMechanicsBiotPressure.h"
#include "MooseMesh.h"

registerMooseObject("felinoApp", ADMechanicsBiotPressure);

InputParameters
ADMechanicsBiotPressure::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("biot_coef_name",
        "Name of material property for Biot coefficient (e.g. 'biot_coefficient').");
  params.addRequiredCoupledVar("p_eff_variable",
        "Aux/coupled variable holding effective fluid pressure (lagged).");
  params.addRequiredParam<unsigned int>("component",
        "Displacement component (0:x, 1:y, 2:z).");
  return params;
}

ADMechanicsBiotPressure::ADMechanicsBiotPressure(const InputParameters & params)
  : ADKernel(params),
    _alpha(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("biot_coef_name"))),
    _p_eff(coupledValue("p_eff_variable")),
    _component(getParam<unsigned int>("component")),
    _rz(getBlockCoordSystem() == Moose::COORD_RZ)
{
  if (_component >= _mesh.dimension())
    paramError("component", "component exceeds mesh dimension");
}

ADReal
ADMechanicsBiotPressure::computeQpResidual()
{
  if (_rz && _component == 0)
    return -_alpha[_qp] * _p_eff[_qp] *
           (_grad_test[_i][_qp](0) + _test[_i][_qp] / _q_point[_qp](0));
  return -_alpha[_qp] * _p_eff[_qp] * _grad_test[_i][_qp](_component);
}
