#include "PorousFlowEffectiveStressCouplingVariableBiot.h"

#include "Function.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

registerMooseObject("PorousFlowApp", PorousFlowEffectiveStressCouplingVariableBiot);

InputParameters
PorousFlowEffectiveStressCouplingVariableBiot::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Implements the weak form of the expression biot_coefficient * "
                             "grad(effective fluid pressure)");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  //params.addRangeCheckedParam<Real>(
  //    "biot_coefficient", 1, "biot_coefficient>=0&biot_coefficient<=1", "Biot coefficient");
  params.addParam<MaterialPropertyName>(
      "biot_coef_name", "biot_coefficient", "Name of material property for biot coefficient.");//✅0402  
  params.addRequiredParam<unsigned int>("component",
                                        "The component (0 for x, 1 for y and 2 for z) of grad(P)");
  return params;
}

PorousFlowEffectiveStressCouplingVariableBiot::PorousFlowEffectiveStressCouplingVariableBiot(
    const InputParameters & parameters)
  : Kernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    //_coefficient(getParam<Real>("biot_coefficient")),
    _biot_coef_name(getParam<MaterialPropertyName>("biot_coef_name")),// ✅
    _coefficient(getMaterialProperty<Real>(_biot_coef_name)),// ✅
    _component(getParam<unsigned int>("component")),
    _pf(getMaterialProperty<Real>("PorousFlow_effective_fluid_pressure_qp")),
    _dpf_dvar(
        getMaterialProperty<std::vector<Real>>("dPorousFlow_effective_fluid_pressure_qp_dvar")),
    _rz(getBlockCoordSystem() == Moose::COORD_RZ)
{
  if (_component >= _mesh.dimension())
    paramError("component", "The component cannot be greater than the mesh dimension");
}

Real
PorousFlowEffectiveStressCouplingVariableBiot::computeQpResidual()
{
  if (_rz && _component == 0)
    return -_coefficient[_qp] * _pf[_qp] * (_grad_test[_i][_qp](0) + _test[_i][_qp] / _q_point[_qp](0));
  return -_coefficient[_qp] * _pf[_qp] * _grad_test[_i][_qp](_component);
}

Real
PorousFlowEffectiveStressCouplingVariableBiot::computeQpJacobian()
{
  if (_dictator.notPorousFlowVariable(_var.number()))
    return 0.0;
  const unsigned int pvar = _dictator.porousFlowVariableNum(_var.number());
  if (_rz && _component == 0)
    return -_coefficient[_qp] * _phi[_j][_qp] * _dpf_dvar[_qp][pvar] *
           (_grad_test[_i][_qp](0) + _test[_i][_qp] / _q_point[_qp](0));
  return -_coefficient[_qp] * _phi[_j][_qp] * _dpf_dvar[_qp][pvar] * _grad_test[_i][_qp](_component);
}

Real
PorousFlowEffectiveStressCouplingVariableBiot::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);
  if (_rz && _component == 0)
    return -_coefficient[_qp] * _phi[_j][_qp] * _dpf_dvar[_qp][pvar] *
           (_grad_test[_i][_qp](0) + _test[_i][_qp] / _q_point[_qp](0));
  return -_coefficient[_qp] * _phi[_j][_qp] * _dpf_dvar[_qp][pvar] * _grad_test[_i][_qp](_component);
}
