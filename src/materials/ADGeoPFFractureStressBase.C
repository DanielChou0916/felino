//* This file is the customized AD phase field fracture model.

#include "ADGeoPFFractureStressBase.h"

InputParameters
ADGeoPFFractureStressBase::validParams()
{
  InputParameters params = ADComputeStressBase::validParams();
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
  params.addParam<MaterialPropertyName>(
      "F_name",
      "local_fracture_energy",
      "Name of material property for local fracture energy function.");
  // Add new parameters here
  return params;
}

ADGeoPFFractureStressBase::ADGeoPFFractureStressBase(const InputParameters & parameters)
  : ADComputeStressBase(parameters), 
    DerivativeMaterialPropertyNameInterface(),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getADMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    //_Jacobian_mult(declareADProperty<RankFourTensor>("degraded_stiffness_tensor")),
    //old: _c(coupledValue("c")),
    _c_name(getVar("c", 0)->name()),
    _l(getADMaterialProperty<Real>("l")),
    _pressure(getOptionalADMaterialProperty<Real>("fracture_pressure")), // ✅Optional variable
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _use_snes_vi_solver(getParam<bool>("use_snes_vi_solver")),
    _H(declareADProperty<Real>("hist")),
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _barrier(getOptionalADMaterialProperty<Real>("barrier_energy")), //✅Optional variable

    _E_name(getParam<MaterialPropertyName>("E_name")),  // ✅ 
    _E(declareADProperty<Real>(_E_name)),  // ✅ The elastic energy is output calcuated by this material, hence use declareADProperty for initialization.
    _dEdc(declareADProperty<Real>(derivativePropertyName(_E_name, {_c_name}))),  // ✅ create a derivative by declareADProperty
    _d2Ed2c(declareADProperty<Real>(derivativePropertyName(_E_name, {_c_name, _c_name}))),  // ✅ 

    _dstress_dc(declareADProperty<RankTwoTensor>(derivativePropertyName(_base_name + "stress", {_c_name}))), // ✅
    _d2Fdcdstrain(declareADProperty<RankTwoTensor>("d2Fdcdstrain")),

    _D_name(getParam<MaterialPropertyName>("D_name")),// ✅
    _D(getADMaterialProperty<Real>(_D_name)),// ✅degradation function is an input, need to be created else where in input file, hence use getADMaterialProperty instead of declareADProperty
    //_dDdc(declareADProperty<Real>(derivativePropertyName(_D_name, {_c_name}))),// ✅
    //_d2Dd2c(declareADProperty<Real>(derivativePropertyName(_D_name, {_c_name, _c_name}))),// ✅
    _dDdc(getADMaterialProperty<Real>(derivativePropertyName(_D_name, {_c_name}))),// ✅ ✅
    _d2Dd2c(getADMaterialProperty<Real>(derivativePropertyName(_D_name, {_c_name, _c_name}))),// ✅ ✅

    _I_name(getParam<MaterialPropertyName>("I_name")), // ✅
    _I(getOptionalADMaterialProperty<Real>(_I_name)), // ✅ Optional variable
    _dIdc(getOptionalADMaterialProperty<Real>(derivativePropertyName(_I_name, {_c_name}))), // ✅
    _d2Id2c(getOptionalADMaterialProperty<Real>(derivativePropertyName(_I_name, {_c_name, _c_name}))), // ✅
    // Initialize new parameters here
    //Drucker-Prager Based sigmoid approach
    _internal_friction(getOptionalADMaterialProperty<Real>("internal_friction")),
    _cohesion(getOptionalADMaterialProperty<Real>("cohesion")),
    _R_res(getOptionalADMaterialProperty<Real>("R_res")),
    _R_init(getOptionalADMaterialProperty<Real>("R_init"))
{
}

void
ADGeoPFFractureStressBase::initQpStatefulProperties()
{
  _H[_qp] = 0.0;

}


