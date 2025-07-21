#include "ADPFFRCEStress.h"

registerMooseObject("SolidMechanicsApp", ADPFFRCEStress);

InputParameters
ADPFFRCEStress::validParams()
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
  params.addParam<std::string>("normal_vector_mode", "stress_based", 
      "Choose the crack normal formulation: 'stress_based' or 'strain_based'");
  return params;
}

ADPFFRCEStress::ADPFFRCEStress(const InputParameters & parameters)
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
    _elasticity_tensor(getADMaterialPropertyByName<RankFourTensor>(_uncracked_base_name + "elasticity_tensor")),
    //_uncracked_Jacobian_mult(// ❌
    //    getMaterialPropertyByName<RankFourTensor>(_uncracked_base_name + "Jacobian_mult")),// ❌
    _c_name(getVar("c", 0)->name()),
    _c(coupledValue("c")),
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
    _normal_vector_mode(getParam<std::string>("normal_vector_mode"))
    //_Jacobian_mult(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult")),// ❌
{
}

void
ADPFFRCEStress::initQpStatefulProperties()
{
  _stress[_qp].zero();
  _hist[_qp] = 0.0;
}

void
ADPFFRCEStress::computeQpProperties()
{ ADReal psi_0, psi_res,psi_act;
  ADRankTwoTensor stress_res,stress_act;
  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  mooseInfo("Using RCE formulation.");
  // Create the positive and negative projection tensors
  //ADRankFourTensor I4sym(ADRankFourTensor::initIdentitySymmetricFour);
  std::vector<ADReal> eigvals;
  ADRankTwoTensor eigvecs;
  if (_normal_vector_mode == "stress_based")
  {
    _uncracked_stress[_qp].symmetricEigenvaluesEigenvectors(eigvals, eigvecs);
  }
  else if (_normal_vector_mode == "strain_based")
  {
    _strain[_qp].symmetricEigenvaluesEigenvectors(eigvals, eigvecs);
  }
  else
    mooseError("Invalid input for calculating crack normal vector.");
  //Define projections
  std::vector<ADRankTwoTensor> projections(LIBMESH_DIM);
  for (const auto i : make_range(Moose::dim))
    projections[i] =0.5*
                    (ADRankTwoTensor::outerProduct(eigvecs.column(Moose::dim-1),eigvecs.column(i))
                    +ADRankTwoTensor::outerProduct(eigvecs.column(i),eigvecs.column(Moose::dim-1)));
  // Define Lambda
  // Isotropic case 
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  std::vector<ADReal> Lam(LIBMESH_DIM);
  for (const auto i : make_range(Moose::dim))
    if (i==Moose::dim-1)
    {
      ADReal L1 = lambda/(lambda+2*mu)* _strain[_qp].trace()
                + 2*mu/(lambda+2*mu)* (_strain[_qp]).doubleContraction(projections[i]);
      Lam[i] = std::max(L1,0.0);
    }
    else 
    {
      Lam[i] = 2* (_strain[_qp]).doubleContraction(projections[i]);
    }

  // compute jump strain
  ADRankTwoTensor strain_jump;
  for (const auto i : make_range(Moose::dim))
    strain_jump += Lam[i]*projections[i];

  ADRankTwoTensor delta_strain = _strain[_qp] - strain_jump;

  // Project the positive and negative stresses
  stress_res = _elasticity_tensor[_qp]*delta_strain;
             //_elasticity_tensor[_qp]*_strain[_qp]-_elasticity_tensor[_qp]*strain_jump;
              //_uncracked_stress[_qp]-_elasticity_tensor[_qp]*strain_jump;
  stress_act = _elasticity_tensor[_qp]*strain_jump;//_uncracked_stress[_qp] - stress_res;

  psi_res = stress_res.doubleContraction(delta_strain) / 2.0;
  psi_0 = (_elasticity_tensor[_qp]*_strain[_qp]).doubleContraction(_strain[_qp]) / 2.0;//_uncracked_stress[_qp].doubleContraction(_strain[_qp]) / 2.0;
  psi_act = psi_0 - psi_res;

  // Update the history variable
  if (_c[_qp] < 0.85)
    _hist[_qp] = psi_act;
  else
    _hist[_qp] = std::max(psi_act, _hist_old[_qp]);

  ADReal hist_variable = _hist_old[_qp];
  if (_use_current_hist)
    hist_variable = _hist[_qp];

  // Compute stress and its derivatives
  if (_I_name.empty() || !(_I) || !(_pressure))
  { //Stress
    _stress[_qp] = _D[_qp] * stress_act + stress_res;
    _dstress_dc[_qp] = stress_act * _dDdc[_qp];
    //Energy
    _E[_qp] = hist_variable * _D[_qp] + psi_res;
    _dEdc[_qp] = hist_variable * _dDdc[_qp];
    _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];
  }
  else   
  { //Stress
    _stress[_qp] = _D[_qp] * stress_act - _pressure[_qp] * I2 * _I[_qp] + stress_res;
    _dstress_dc[_qp] = stress_act * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
    //Energy
    _E[_qp] =
        hist_variable * _D[_qp] + psi_res - _pressure[_qp] * _strain[_qp].trace() * _I[_qp];
    _dEdc[_qp] =
        hist_variable * _dDdc[_qp] - _pressure[_qp] * _strain[_qp].trace() * _dIdc[_qp];
    _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp] -
                  _pressure[_qp] * _strain[_qp].trace() * _d2Id2c[_qp];    
  }


  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress_act * _dDdc[_qp];
}
