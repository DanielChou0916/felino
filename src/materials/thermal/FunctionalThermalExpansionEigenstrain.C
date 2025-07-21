

#include "FunctionalThermalExpansionEigenstrain.h"

registerMooseObject("SolidMechanicsApp", FunctionalThermalExpansionEigenstrain);
registerMooseObject("SolidMechanicsApp", ADFunctionalThermalExpansionEigenstrain);

template <bool is_ad>
InputParameters
FunctionalThermalExpansionEigenstrainTempl<is_ad>::validParams()
{
  InputParameters params = ComputeThermalExpansionEigenstrainBaseTempl<is_ad>::validParams();
  params.addClassDescription("Computes eigenstrain due to thermal expansion "
                             "with a constant coefficient");
  params.addRequiredParam<MaterialPropertyName>("thermal_expansion_coeff", "Thermal expansion coefficient");

  return params;
}

template <bool is_ad>
FunctionalThermalExpansionEigenstrainTempl<is_ad>::FunctionalThermalExpansionEigenstrainTempl(
    const InputParameters & parameters)
  : ComputeThermalExpansionEigenstrainBaseTempl<is_ad>(parameters),
    //_thermal_expansion_coeff(this->template getParam<Real>("thermal_expansion_coeff"))
    _thermal_expansion_coeff(this->template getGenericMaterialProperty<Real, is_ad>("thermal_expansion_coeff"))
{
}

template <bool is_ad>
ValueAndDerivative<is_ad>
FunctionalThermalExpansionEigenstrainTempl<is_ad>::computeThermalStrain()
{
  return _thermal_expansion_coeff[_qp] * (_temperature[_qp] - _stress_free_temperature[_qp]);
}

template class FunctionalThermalExpansionEigenstrainTempl<false>;
template class FunctionalThermalExpansionEigenstrainTempl<true>;
