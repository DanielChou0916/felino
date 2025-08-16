

#pragma once

#include "ADTimeKernelValue.h"

class ADHydroTimeDerivative : public ADTimeKernelValue
{
public:
  static InputParameters validParams();

  ADHydroTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  const ADMaterialProperty<Real> & _fluid_density;
  const ADMaterialProperty<Real> & _biot_modulus;
};
