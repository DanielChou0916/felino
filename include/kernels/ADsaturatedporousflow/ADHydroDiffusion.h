

#pragma once

#include "ADKernel.h"
#include "RankTwoTensor.h"

/**
 * This ADkernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class ADHydroDiffusion : public ADKernel
{
public:
  static InputParameters validParams();

  ADHydroDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<Real> & _fluid_density;

  const ADMaterialProperty<Real> & _fluid_viscosity;

  const ADMaterialProperty<RankTwoTensor> & _permeability_tensor;
};
