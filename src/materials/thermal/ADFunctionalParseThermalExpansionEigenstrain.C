#include "ADFunctionalParseThermalExpansionEigenstrain.h"

registerMooseObject("SolidMechanicsApp", ADFunctionalParseThermalExpansionEigenstrain);

InputParameters
ADFunctionalParseThermalExpansionEigenstrain::validParams()
{
  InputParameters params = ComputeThermalExpansionEigenstrainBaseTempl<true>::validParams();
  params.addClassDescription("Computes eigenstrain due to thermal expansion with variable coefficient.");
  params.addRequiredParam<MaterialPropertyName>("thermal_expansion_coeff", "Thermal expansion coefficient");
  return params;
}

ADFunctionalParseThermalExpansionEigenstrain::ADFunctionalParseThermalExpansionEigenstrain(
    const InputParameters & parameters)
  : ComputeThermalExpansionEigenstrainBaseTempl<true>(parameters),
    _thermal_expansion_coeff(this->template getADMaterialProperty<Real>("thermal_expansion_coeff"))
{
}

ADReal
ADFunctionalParseThermalExpansionEigenstrain::computeThermalStrain()
{
  return _thermal_expansion_coeff[_qp] * (_temperature[_qp] - _stress_free_temperature[_qp]);
}
