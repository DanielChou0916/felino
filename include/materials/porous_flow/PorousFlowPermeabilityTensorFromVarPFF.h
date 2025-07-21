

#pragma once

#include "PorousFlowPermeabilityBase.h"

/**
 * Material designed to provide the permeability tensor which is calculated
 * from a tensor multiplied by a scalar:
 * k = k_ijk * k0
 * where k_ijk is a tensor providing the anisotropy, and k0 is a scalar
 * variable.
 */
template <bool is_ad>
class PorousFlowPermeabilityTensorFromVarPFFTempl : public PorousFlowPermeabilityBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  PorousFlowPermeabilityTensorFromVarPFFTempl(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  /// Permeability components
  /// Note: these can only be constant (Real constant Monomial auxvariables) so no AD version
  const VariableValue & _perm;

  /// Tensor multiplier k_ijk
  const RealTensorValue _k_anisotropy;
  /// Fracture normal orientation matrix n_d cross n_d
  ///MaterialProperty<RankTwoTensor> & _crack_direction;// ✅
  const MaterialPropertyName _crack_direction_name;// ✅
  const MaterialProperty<RankTwoTensor> & _crack_direction;// ✅
  const MaterialPropertyName _coef_name;// ✅
  const MaterialProperty<Real> & _A;// ✅
  const std::string _model;// ✅
  usingPorousFlowPermeabilityBaseMembers;
};

typedef PorousFlowPermeabilityTensorFromVarPFFTempl<false> PorousFlowPermeabilityTensorFromVarPFF;
//typedef PorousFlowPermeabilityTensorFromVarPFFTempl<true> ADPorousFlowPermeabilityTensorFromVarPFF;
