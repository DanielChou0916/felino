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
class ADPFFRCEStress :public Material, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  ADPFFRCEStress(const InputParameters & parameters);

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
  const ADMaterialProperty<RankTwoTensor> & _strain;

  /// Uncracked stress calculated by another material
  const ADMaterialProperty<RankTwoTensor> & _uncracked_stress;

  /// Uncracked Jacobian_mult calculated by another material
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;
  //const MaterialProperty<RankFourTensor> & _uncracked_Jacobian_mult;// ❌

  /// Variable defining the phase field damage parameter
  //const VariableValue & _c;// ❌
  const VariableName _c_name;// ✅
  const VariableValue & _c;//

  /// Stress being computed by this kernel
  ADMaterialProperty<RankTwoTensor> & _stress;

  const MaterialPropertyName _E_name;// ✅
  ADMaterialProperty<Real> & _E;// ✅
  ADMaterialProperty<Real> & _dEdc;// ✅
  ADMaterialProperty<Real> & _d2Ed2c;// ✅
  ADMaterialProperty<RankTwoTensor> & _d2Fdcdstrain;
  ADMaterialProperty<RankTwoTensor> & _dstress_dc;

  /// history variable storing the maximum positive deformation energy
  ADMaterialProperty<Real> & _hist;
  const MaterialProperty<Real> & _hist_old;

  /// Material property defining pressure, declared elsewhere
  const OptionalADMaterialProperty<Real> & _pressure;// ✅
  /// Material property for damage indicator function
  const MaterialPropertyName _I_name;// ✅
  const OptionalADMaterialProperty<Real> & _I;

  /// Derivative of damage indicator function w.r.t damage variable
  const OptionalADMaterialProperty<Real> & _dIdc;

  /// Second-order derivative of damage indicator function w.r.t damage variable
  const OptionalADMaterialProperty<Real> & _d2Id2c;
  /// Material property for energetic degradation function
  const MaterialPropertyName _D_name;// ✅
  const ADMaterialProperty<Real> & _D;
  /// Derivative of degradation function w.r.t damage variable
  const ADMaterialProperty<Real> & _dDdc;

  /// Second-order derivative of degradation w.r.t damage variable
  const ADMaterialProperty<Real> & _d2Dd2c;
  const std::string _normal_vector_mode;
  /// derivative of stress w.r.t. strain (_dstress_dstrain)
  //MaterialProperty<RankFourTensor> & _Jacobian_mult;// ❌

};
