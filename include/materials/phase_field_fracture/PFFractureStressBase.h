//* This file is the customized phase field fracture model.

#pragma once

#include "ComputeStressBase.h"
//#include "DerivativeMaterialPropertyNameInterface.h"
/**
  PFFractureStressBase is the base class for stress in phase field fracture model
 */
class PFFractureStressBase : public ComputeStressBase//,
                                  //public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  PFFractureStressBase(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  /// Extra: Degradated Stiffness
  //MaterialProperty<RankFourTensor> & _Jacobian_mult;// ✅
  /// Coupled order parameter defining the crack
  //old: const VariableValue & _c;
  const VariableName _c_name;// ✅
  /// Material property defining crack width, declared elsewhere
  //const MaterialProperty<Real> & _l;

  /// Material property defining pressure, declared elsewhere
  const OptionalMaterialProperty<Real> & _pressure;// ✅

  /// Use current value of history variable
  bool _use_current_hist;

  /// Use PETSc's VI (Reduced space active set solvers for variational inequalities based on Newton's method) solver
  bool _use_snes_vi_solver;

  /// History variable that prevents crack healing, declared in this material
  MaterialProperty<Real> & _H;

  /// Old value of history variable
  const MaterialProperty<Real> & _H_old;

  /// material property for fracture energy barrier
  const OptionalMaterialProperty<Real> & _barrier;

  /// Material property for elastic energy
  const MaterialPropertyName _E_name;// ✅
  MaterialProperty<Real> & _E;

  /// Derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _dEdc;

  /// Second-order derivative of elastic energy w.r.t damage variable
  MaterialProperty<Real> & _d2Ed2c;

  /// Derivative of stress w.r.t damage variable
  MaterialProperty<RankTwoTensor> & _dstress_dc;

  /// Second-order derivative of elastic energy w.r.t damage variable and strain
  MaterialProperty<RankTwoTensor> & _d2Fdcdstrain;

  /// Material property for energetic degradation function
  const MaterialPropertyName _D_name;// ✅
  const MaterialProperty<Real> & _D;

  /// Derivative of degradation function w.r.t damage variable
  const MaterialProperty<Real> & _dDdc;

  /// Second-order derivative of degradation w.r.t damage variable
  const MaterialProperty<Real> & _d2Dd2c;

  /// Material property for damage indicator function
  const MaterialPropertyName _I_name;// ✅
  const OptionalMaterialProperty<Real> & _I;

  /// Derivative of damage indicator function w.r.t damage variable
  const OptionalMaterialProperty<Real> & _dIdc;

  /// Second-order derivative of damage indicator function w.r.t damage variable
  const OptionalMaterialProperty<Real> & _d2Id2c;
  // Member variables for the new parameters
  const OptionalMaterialProperty<Real> & _B;
//
};
