#include "ADHydroStrainRate.h"
registerMooseObject("MooseApp", ADHydroStrainRate);

InputParameters ADHydroStrainRate::validParams()
{
  auto p = Material::validParams();
  p.addParam<std::string>(
      "base_name", "",
      "Prefix matching the TensorMechanics strain calculator (e.g. 'tm_'). "
      "Will read '<base_name>total_strain'.");
  p.addParam<MaterialPropertyName>(
      "vol_strain_rate_name", "vol_strain_rate",
      "Name of the AD property to declare for volumetric strain rate.");
  p.addClassDescription("Compute volumetric strain rate: (tr(total_strain) - tr(total_strain_old)) / dt");
  return p;
}

ADHydroStrainRate::ADHydroStrainRate(const InputParameters & params)
  : Material(params),
    _base_name(getParam<std::string>("base_name")),
    _total_strain(getADMaterialProperty<RankTwoTensor>(_base_name + "total_strain")),
    _total_strain_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "total_strain")),
    _vol_strain_rate(declareADProperty<Real>(getParam<MaterialPropertyName>("vol_strain_rate_name")))
{
}

void ADHydroStrainRate::computeQpProperties()
{
  const ADReal eps_v     = _total_strain[_qp].trace();          // AD
  const Real   eps_v_old = _total_strain_old[_qp].trace();       // 常數（相對本步）
  _vol_strain_rate[_qp] = (eps_v - eps_v_old) / _dt;                 // AD 自動帶 off-diag
}
