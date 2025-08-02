

#pragma once

#include "DerivativeFunctionMaterialBase.h"

template <bool is_ad>
class DerivativeSumMaterialWithConstantOn1stOrderTempl : public DerivativeFunctionMaterialBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  DerivativeSumMaterialWithConstantOn1stOrderTempl(const InputParameters & parameters);

  virtual void initialSetup();

protected:
  usingDerivativeFunctionMaterialBaseMembers(is_ad);

  virtual void computeProperties();

  std::vector<std::string> _sum_materials;// ✅
  unsigned int _num_materials;// ✅

  /// arguments to construct a sum of the form \f$ c+\gamma\sum_iF_i \f$
  std::vector<Real> _prefactor;
  Real _constant;

  /// Flag to optionally turn on or off validateCoupling
  const bool _validate_coupling;

  /// Function values of the summands.
  std::vector<const GenericMaterialProperty<Real, is_ad> *> _summand_F;// ✅

  /// Derivatives of the summands with respect to arg[i]
  std::vector<std::vector<const GenericMaterialProperty<Real, is_ad> *>> _summand_dF;

  /// Second derivatives of the summands.
  std::vector<std::vector<std::vector<const GenericMaterialProperty<Real, is_ad> *>>> _summand_d2F;

  /// Third derivatives of the summands.
  std::vector<std::vector<std::vector<std::vector<const GenericMaterialProperty<Real, is_ad> *>>>>
      _summand_d3F;
  /// Additional varlable: for extra driving force
  /// Additional function values (not needing derivatives)
  std::vector<std::string> _additional_summands_materials;// ✅
  unsigned int _num_constant_materials;// ✅
  std::vector<const GenericMaterialProperty<Real, is_ad> *> _add_summand_F;// ✅

  /// For outputting individual summand derivatives
  const bool _output_individual_derivatives; // ✅ 0801
  std::vector<GenericMaterialProperty<Real, is_ad> *> _summand_dF_out;// ✅0801

};

typedef DerivativeSumMaterialWithConstantOn1stOrderTempl<false> DerivativeSumMaterialWithConstantOn1stOrder;
typedef DerivativeSumMaterialWithConstantOn1stOrderTempl<true> ADDerivativeSumMaterialWithConstantOn1stOrder;
