#pragma once

#include "Material.h"
#include "RankTwoTensor.h"

/**
 * Generate anisotropic matrix A = I + coef * (I - n ⊗ n)
 * from user-defined weak plane normal vector n
 */
class ADAnisotropicDirector : public Material
{
public:
  static InputParameters validParams();
  ADAnisotropicDirector(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Input: manually defined weak plane normal
  const RealVectorValue _normal;

  /// Input: anisotropy scaling coefficient
  const Real _coef;

  /// Output: directional tensor A
  const MaterialPropertyName _output_name;
  ADMaterialProperty<RankTwoTensor> & _directional_tensor;
};
