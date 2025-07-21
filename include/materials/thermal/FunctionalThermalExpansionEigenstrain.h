
#pragma once

#include "ComputeThermalExpansionEigenstrainBase.h"

/**
 * ComputeThermalExpansionEigenstrain computes an eigenstrain for thermal expansion
 * with a functional expansion coefficient.
 */
template <bool is_ad>
class FunctionalThermalExpansionEigenstrainTempl
  : public ComputeThermalExpansionEigenstrainBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  FunctionalThermalExpansionEigenstrainTempl(const InputParameters & parameters);

protected:
  virtual ValueAndDerivative<is_ad> computeThermalStrain() override;

  //const Real & _thermal_expansion_coeff;
  const GenericMaterialProperty<Real, is_ad> & _thermal_expansion_coeff;
  using ComputeThermalExpansionEigenstrainBaseTempl<is_ad>::_qp;
  using ComputeThermalExpansionEigenstrainBaseTempl<is_ad>::_temperature;
  using ComputeThermalExpansionEigenstrainBaseTempl<is_ad>::_stress_free_temperature;
};

typedef FunctionalThermalExpansionEigenstrainTempl<false> FunctionalThermalExpansionEigenstrain;
typedef FunctionalThermalExpansionEigenstrainTempl<true> ADFunctionalThermalExpansionEigenstrain;
