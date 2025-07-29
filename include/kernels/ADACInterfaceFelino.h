#pragma once

#include "ADKernel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Extended ACInterface Kernel that includes an additional gradient kappa term.
 *
 * Weak form:
 *   ∫ ( kappa * ∇d ) · ∇(L * phi) dΩ
 * + ∫ ( ∇kappa · ∇d ) * L * phi dΩ
 */
class ADACInterfaceFelino : public ADKernel, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();
  ADACInterfaceFelino(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// Mobility L
  const ADMaterialProperty<Real> & _prop_L;
  const MaterialPropertyName & _name_L;

  /// Interfacial parameter kappa
  const ADMaterialProperty<Real> & _kappa;
  
  /// Flag to indicate if L is variable-dependent
  const bool _variable_L;

  /// Derivative of L w.r.t. order parameter
  const ADMaterialProperty<Real> * const _dLdop;

  /// Number of coupled variables
  const unsigned int _nvar;

  /// Derivatives of L w.r.t. coupled variables
  std::vector<const ADMaterialProperty<Real> *> _dLdarg;

  /// Gradients of coupled variables
  std::vector<const ADVariableGradient *> _gradarg;

  const bool _use_anisotropic_matrix;// ✅ 2025/07/06
  const ADMaterialProperty<RankTwoTensor> * _A_ptr;// ✅ 2025/07/06
};
