#pragma once

#include "Material.h"
#include "RankTwoTensor.h"

/**
 * Pure-AD permeability material (no PorousFlow).
 * k = K_aniso * perm  (+ optional model term)
 *   - model = 'directional': k = K_aniso * perm + A * (I - P),  P = crack_direction_property
 *   - model = 'exponential': k = K_aniso * perm * exp(A)
 */
class ADPermeabilityFromVar : public Material
{
public:
  static InputParameters validParams();
  ADPermeabilityFromVar(const InputParameters & params);

protected:
  void computeQpProperties() override;

  // inputs
  const ADMaterialProperty<Real> & _perm; // AD coupled variable
  const RealTensorValue _k_aniso;

  const bool _use_directional;
  const bool _use_exponential;

  // optional inputs (pointers may be null if not provided)
  const ADMaterialProperty<RankTwoTensor> * _direction = nullptr; // projection (e.g. n⊗n)
  const ADMaterialProperty<Real>        * _A_prop   = nullptr;    // A as material property
  const Real                              _A_const;               // fallback constant A

  // output
  ADMaterialProperty<RankTwoTensor> & _k_tensor;
};
