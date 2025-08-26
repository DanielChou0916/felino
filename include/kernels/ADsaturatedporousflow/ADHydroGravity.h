

#pragma once

#include "ADKernel.h"
#include "RankTwoTensor.h"

class ADHydroGravity : public ADKernel
{
public:
  static InputParameters validParams();

  ADHydroGravity(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<Real> & _fluid_density;

  const ADMaterialProperty<Real> & _fluid_viscosity;

  const ADMaterialProperty<RankTwoTensor> & _permeability_tensor;
  
  const RealVectorValue _gravity;

  const bool _mass_flux; // ✅ 0826
};
