#pragma once

#include "ComputeThermalExpansionEigenstrainBase.h"

/**
 * ADFunctionalParseThermalExpansionEigenstrain computes an eigenstrain for thermal expansion
 * with a spatially or temporally varying coefficient (via material property).
 */
class ADFunctionalParseThermalExpansionEigenstrain
  : public ComputeThermalExpansionEigenstrainBaseTempl<true>
{
public:
  static InputParameters validParams();

  ADFunctionalParseThermalExpansionEigenstrain(const InputParameters & parameters);

protected:
  virtual ADReal computeThermalStrain() override;

  const ADMaterialProperty<Real> & _thermal_expansion_coeff;

  using ComputeThermalExpansionEigenstrainBaseTempl<true>::_qp;
  using ComputeThermalExpansionEigenstrainBaseTempl<true>::_temperature;
  using ComputeThermalExpansionEigenstrainBaseTempl<true>::_stress_free_temperature;
};
