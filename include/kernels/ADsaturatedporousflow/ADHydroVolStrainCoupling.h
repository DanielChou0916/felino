

#pragma once

#include "ADKernel.h"

class ADHydroVolStrainCoupling : public ADKernel
{
public:
  static InputParameters validParams();

  ADHydroVolStrainCoupling(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  const ADMaterialProperty<Real> & _fluid_density;
  const ADMaterialProperty<Real> & _biot_coefficient;
  const ADMaterialProperty<Real> &  _vol_strain_rate;
};
