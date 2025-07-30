#pragma once

#include "PorousFlowPermeabilityBase.h"

/**
 * Material to compute permeability tensor as:
 *   k = k_ijk * k0
 * where k_ijk is an anisotropy tensor and k0 is a scalar variable.
 * Depending on the selected model, an additional term is added to account
 * for fracture plane permeability changes.
 */
template <bool is_ad>
class PorousFlowPermeabilityTensorFromVarPFFTempl : public PorousFlowPermeabilityBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  PorousFlowPermeabilityTensorFromVarPFFTempl(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Coupled scalar permeability
  const VariableValue & _perm;

  /// Tensor multiplier k_ijk
  const RealTensorValue _k_anisotropy;

  /// Parameter names
  const MaterialPropertyName _crack_direction_name;
  const MaterialPropertyName _coef_name;

  /// Selected model ("directional" or "exponential")
  std::string _model;

  /// Material properties (initialized conditionally)
  const MaterialProperty<RankTwoTensor> * _crack_direction = nullptr;
  const MaterialProperty<Real> * _A = nullptr;

  usingPorousFlowPermeabilityBaseMembers;
};

typedef PorousFlowPermeabilityTensorFromVarPFFTempl<false> PorousFlowPermeabilityTensorFromVarPFF;
//typedef PorousFlowPermeabilityTensorFromVarPFFTempl<true> ADPorousFlowPermeabilityTensorFromVarPFF;
