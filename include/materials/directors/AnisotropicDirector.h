#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "MooseEnum.h"

/**
 * Generate anisotropic matrix A = I + coef * (I - n ⊗ n)
 * from user-defined weak plane normal vector n
 */
class AnisotropicDirector : public Material
{
public:
  static InputParameters validParams();
  AnisotropicDirector(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Input: manually defined weak plane normal
  RealVectorValue _normal;

  /// Input: anisotropy scaling coefficient
  const Real _coef;

  /// Output name
  const MaterialPropertyName _output_name;

  /// Output: directional tensor A
  MaterialProperty<RankTwoTensor> & _directional_tensor;

  /// Input mode
  enum class input_type
  {
    normal_vector,
    xy_angle
  } _input_type;

  /// Normalization mode
  enum class norm_type
  { 
    factorial_norm,
    trace_norm,
    det_norm,
    none
  } _normalize_director;

  const Real _factor;
};