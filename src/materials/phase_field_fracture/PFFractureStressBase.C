//* This file is the customized AD phase field fracture model.

#include "PFFractureStressBase.h"

InputParameters
PFFractureStressBase::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<bool>("use_snes_vi_solver",
                        false,
                        "Use PETSc's SNES variational inequalities solver to enforce damage "
                        "irreversibility condition and restrict damage value <= 1.");
  params.addParam<MaterialPropertyName>("barrier_energy",
                                        "Name of material property for fracture energy barrier.");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<MaterialPropertyName>(
      "I_name", "indicator", "Name of material property for damage indicator function.");
  // Add new parameters here
  params.addParam<MaterialPropertyName>("B", "Coefficient for Drucker Prager model");
  return params;
}

PFFractureStressBase::PFFractureStressBase(const InputParameters & parameters)
  : ComputeStressBase(parameters), 
    //DerivativeMaterialPropertyNameInterface(),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    //_Jacobian_mult(declareProperty<RankFourTensor>("degraded_stiffness_tensor")),
    //old: _c(coupledValue("c")),
    _c_name(getVar("c", 0)->name()),
    //_l(getMaterialProperty<Real>("l")),
    _pressure(getOptionalMaterialProperty<Real>("fracture_pressure")), // ✅Optional variable
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _use_snes_vi_solver(getParam<bool>("use_snes_vi_solver")),
    _H(declareProperty<Real>("hist")),
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _barrier(getOptionalMaterialProperty<Real>("barrier_energy")), //✅Optional variable

    _E_name(getParam<MaterialPropertyName>("E_name")),  // ✅ 
    _E(declareProperty<Real>(_E_name)),  // ✅ The elastic energy is output calcuated by this material, hence use declareProperty for initialization.
    _dEdc(declareProperty<Real>(derivativePropertyName(_E_name, {_c_name}))),  // ✅ create a derivative by declareProperty
    _d2Ed2c(declareProperty<Real>(derivativePropertyName(_E_name, {_c_name, _c_name}))),  // ✅ 

    _dstress_dc(declareProperty<RankTwoTensor>(derivativePropertyName(_base_name + "stress", {_c_name}))), // ✅
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),

    _D_name(getParam<MaterialPropertyName>("D_name")),// ✅
    _D(getMaterialProperty<Real>(_D_name)),// ✅degradation function is an input, need to be created else where in input file, hence use getMaterialProperty instead of declareProperty
    //_dDdc(declareProperty<Real>(derivativePropertyName(_D_name, {_c_name}))),// ✅
    //_d2Dd2c(declareProperty<Real>(derivativePropertyName(_D_name, {_c_name, _c_name}))),// ✅
    _dDdc(getMaterialProperty<Real>(derivativePropertyName(_D_name, {_c_name}))),// ✅ ✅
    _d2Dd2c(getMaterialProperty<Real>(derivativePropertyName(_D_name, {_c_name, _c_name}))),// ✅ ✅

    _I_name(getParam<MaterialPropertyName>("I_name")), // ✅
    _I(getOptionalMaterialProperty<Real>(_I_name)), // ✅ Optional variable
    _dIdc(getOptionalMaterialProperty<Real>(derivativePropertyName(_I_name, {_c_name}))), // ✅
    _d2Id2c(getOptionalMaterialProperty<Real>(derivativePropertyName(_I_name, {_c_name, _c_name}))), // ✅
    // Initialize new parameters here
    _B(getOptionalMaterialProperty<Real>("B"))
{
}

void
PFFractureStressBase::initQpStatefulProperties()
{
  _H[_qp] = 0.0;
}


