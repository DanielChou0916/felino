#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Computes energy and modifies the stress for phase field fracture. Can be used with any
 * constitutive model or elastic symmetry.
 */
class ComputePFFStress :public Material, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  ComputePFFStress(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();
  /// Base name of the stress after being modified to include cracks
  const std::string _base_name;

  /// Base name of the uncracked stress and strain
  const std::string _uncracked_base_name;

  /// Indicator if finite strain model is used, to determine if mechanical_strain or elastic_strain should be used
  bool _finite_strain_model;

  /// Use current value of history variable
  bool _use_current_hist;

  /// Mechanical_strain if finite_strain_model = false, otherwise elastic_strain
  const MaterialProperty<RankTwoTensor> & _strain;

  /// Uncracked stress calculated by another material
  const MaterialProperty<RankTwoTensor> & _uncracked_stress;

  /// Uncracked Jacobian_mult calculated by another material
  const MaterialProperty<RankFourTensor> & _uncracked_Jacobian_mult;// ❌

  /// Variable defining the phase field damage parameter
  //const VariableValue & _c;// ❌
  const VariableName _c_name;// ✅

  /// Stress being computed by this kernel
  MaterialProperty<RankTwoTensor> & _stress;

  const MaterialPropertyName _E_name;// ✅
  MaterialProperty<Real> & _E;// ✅
  MaterialProperty<Real> & _dEdc;// ✅
  MaterialProperty<Real> & _d2Ed2c;// ✅
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;
  MaterialProperty<RankTwoTensor> & _dstress_dc;

  /// history variable storing the maximum positive deformation energy
  MaterialProperty<Real> & _hist;
  const MaterialProperty<Real> & _hist_old;

  /// Material property defining pressure, declared elsewhere
  const OptionalMaterialProperty<Real> & _pressure;// ✅
  /// Material property for damage indicator function
  const MaterialPropertyName _I_name;// ✅
  const OptionalMaterialProperty<Real> & _I;

  /// Derivative of damage indicator function w.r.t damage variable
  const OptionalMaterialProperty<Real> & _dIdc;

  /// Second-order derivative of damage indicator function w.r.t damage variable
  const OptionalMaterialProperty<Real> & _d2Id2c;
  /// Material property for energetic degradation function
  const MaterialPropertyName _D_name;// ✅
  const MaterialProperty<Real> & _D;
  /// Derivative of degradation function w.r.t damage variable
  const MaterialProperty<Real> & _dDdc;

  /// Second-order derivative of degradation w.r.t damage variable
  const MaterialProperty<Real> & _d2Dd2c;
  const std::string _decomposition;
  /// derivative of stress w.r.t. strain (_dstress_dstrain)
  MaterialProperty<RankFourTensor> & _Jacobian_mult;// ❌

};
