#include "ADComputePFFStress.h"

registerMooseObject("SolidMechanicsApp", ADComputePFFStress);

InputParameters
ADComputePFFStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Computes energy and modifies the stress for phase field fracture");
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addParam<bool>("finite_strain_model", false, "The model is using finite strain");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<MaterialPropertyName>(
      "I_name", "indicator", "Name of material property for damage indicator function.");  
  params.addParam<std::string>("base_name", "The base name used to save the cracked stress");
  params.addRequiredParam<std::string>("uncracked_base_name",
                                       "The base name used to calculate the original stress");
  params.addParam<std::string>("decomposition", "none", 
      "Choose the energy decomposition: 'none' or 'spectral'");
  return params;
}

ADComputePFFStress::ADComputePFFStress(const InputParameters & parameters)
  : Material(parameters), DerivativeMaterialPropertyNameInterface(),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _uncracked_base_name(getParam<std::string>("uncracked_base_name") + "_"),
    _finite_strain_model(getParam<bool>("finite_strain_model")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _strain(
        _finite_strain_model
            ? getADMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "elastic_strain")
            : getADMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "mechanical_strain")),
    _uncracked_stress(getADMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "stress")),
    //_uncracked_Jacobian_mult(// ❌
    //    getMaterialPropertyByName<RankFourTensor>(_uncracked_base_name + "Jacobian_mult")),// ❌
    _c_name(getVar("c", 0)->name()),
    _stress(declareADProperty<RankTwoTensor>(_base_name + "stress")),
    _E_name(getParam<MaterialPropertyName>("E_name")),  // ✅ 
    _E(declareADProperty<Real>(_E_name)),  // ✅ The elastic energy is output calcuated by this material, hence use declareADProperty for initialization.
    _dEdc(declareADProperty<Real>(derivativePropertyName(_E_name, {_c_name}))),  // ✅ create a derivative by declareADProperty
    _d2Ed2c(declareADProperty<Real>(derivativePropertyName(_E_name, {_c_name, _c_name}))),  // ✅ 
    _d2Fdcdstrain(declareADProperty<RankTwoTensor>("d2Fdcdstrain")),
    _dstress_dc(declareADProperty<RankTwoTensor>(derivativePropertyName(_base_name + "stress", {_c_name}))), // ✅
    _hist(declareADProperty<Real>("hist")),
    _hist_old(getMaterialPropertyOld<Real>("hist")),
    _pressure(getOptionalADMaterialProperty<Real>("fracture_pressure")), // ✅Optional variable
    _I_name(getParam<MaterialPropertyName>("I_name")), // ✅
    _I(getOptionalADMaterialProperty<Real>(_I_name)), // ✅ Optional variable
    _dIdc(getOptionalADMaterialProperty<Real>(derivativePropertyName(_I_name, {_c_name}))), // ✅
    _d2Id2c(getOptionalADMaterialProperty<Real>(derivativePropertyName(_I_name, {_c_name, _c_name}))), // ✅
    _D_name(getParam<MaterialPropertyName>("D_name")),// ✅
    _D(getADMaterialProperty<Real>(_D_name)),// ✅degradation function is an input, need to be created else where in input file, hence use getADMaterialProperty instead of declareADProperty
    _dDdc(getADMaterialProperty<Real>(derivativePropertyName(_D_name, {_c_name}))),// ✅ ✅
    _d2Dd2c(getADMaterialProperty<Real>(derivativePropertyName(_D_name, {_c_name, _c_name}))),// ✅ ✅
    _decomposition(getParam<std::string>("decomposition"))
    //_Jacobian_mult(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult")),// ❌
{
}

void
ADComputePFFStress::initQpStatefulProperties()
{
  _stress[_qp].zero();
  _hist[_qp] = 0.0;
}

void
ADComputePFFStress::computeQpProperties()
{ ADReal G0_pos, G0_neg;
  ADRankTwoTensor stress0pos, stress0neg;
  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  if (_decomposition == "spectral")
  {
    mooseInfo("Using spectral decomposition.");
    // Create the positive and negative projection tensors
    ADRankFourTensor I4sym(ADRankFourTensor::initIdentitySymmetricFour);
    std::vector<ADReal> eigval;
    ADRankTwoTensor eigvec;
    ADRankFourTensor Ppos = _uncracked_stress[_qp].positiveProjectionEigenDecomposition(eigval, eigvec);
    ADRankFourTensor Pneg = I4sym - Ppos;

    // Project the positive and negative stresses
    stress0pos = Ppos * _uncracked_stress[_qp];
    stress0neg = Pneg * _uncracked_stress[_qp];

    // Compute the positive and negative elastic energies
    G0_pos = (stress0pos).doubleContraction(_strain[_qp]) / 2.0;
    G0_neg = (stress0neg).doubleContraction(_strain[_qp]) / 2.0;
  }
  else if (_decomposition == "none")
  {
    mooseInfo("By default, no decomposition.");
    stress0pos = _uncracked_stress[_qp];
     // Compute the positive and negative elastic energies
    G0_pos = (stress0pos).doubleContraction(_strain[_qp]) / 2.0;
    G0_neg = 0.0;   
  }

  // Update the history variable
  if (G0_pos > _hist_old[_qp])
    _hist[_qp] = G0_pos;
  else
    _hist[_qp] = _hist_old[_qp];

  ADReal hist_variable = _hist_old[_qp];
  if (_use_current_hist)
    hist_variable = _hist[_qp];

  // Compute stress and its derivatives
  if (_I_name.empty() || !(_I) || !(_pressure))
  { //Stress
    _stress[_qp] = _D[_qp] * stress0pos + stress0neg;
    _dstress_dc[_qp] = stress0pos * _dDdc[_qp];
    //Energy
    _E[_qp] = hist_variable * _D[_qp] + G0_neg;
    _dEdc[_qp] = hist_variable * _dDdc[_qp];
    _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];
  }
  else   
  { //Stress
    _stress[_qp] = _D[_qp] * stress0pos - _pressure[_qp] * I2 * _I[_qp] + stress0neg;
    _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
    //Energy
    _E[_qp] =
        hist_variable * _D[_qp] + G0_neg - _pressure[_qp] * _strain[_qp].trace() * _I[_qp];
    _dEdc[_qp] =
        hist_variable * _dDdc[_qp] - _pressure[_qp] * _strain[_qp].trace() * _dIdc[_qp];
    _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp] -
                  _pressure[_qp] * _strain[_qp].trace() * _d2Id2c[_qp];    
  }


  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];
}
