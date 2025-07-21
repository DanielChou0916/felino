//

#include "ADComputeExtraDrivingForcefromStress.h"
#include "MathUtils.h"

registerMooseObject("SolidMechanicsApp", ADComputeExtraDrivingForcefromStress);

InputParameters
ADComputeExtraDrivingForcefromStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("uncracked_base_name", "",
                                "Optional parameter that allows the user to define "
                                "multiple mechanics material systems on the same "
                                "block, i.e. for multiple phases");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");//✅0318
  params.addParam<std::string>("model", "finite_e", 
      "Choose the energy model: 'finite_e' or 'asymptotic'");
  params.addParam<MaterialPropertyName>(
      "Ce_name", "driving_force", "Name of material property for extra driving force");
  params.addParam<bool>("finite_strain_model", false, "The model is using finite strain");
  // Add new parameters here
  return params;
}

ADComputeExtraDrivingForcefromStress::ADComputeExtraDrivingForcefromStress(const InputParameters & parameters)
  : Material(parameters), 
    GuaranteeConsumer(this),
    _uncracked_base_name(getParam<std::string>("uncracked_base_name").empty() ? "" : getParam<std::string>("uncracked_base_name") + "_"),
    _elasticity_tensor_name(_uncracked_base_name + "elasticity_tensor"),
    _elasticity_tensor(getADMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _stress(getADMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "stress")),
    //_finite_strain_model(getParam<bool>("finite_strain_model")),
    //_strain(
    //  _finite_strain_model
    //      ? getADMaterialPropertyByName<RankTwoTensor>(_base_name + "elastic_strain")
    //      : getADMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _D_name(getParam<MaterialPropertyName>("D_name")),
    _D(getADMaterialProperty<Real>(_D_name)),
    _l(getADMaterialProperty<Real>("l")),
    _gc(getADMaterialProperty<Real>("gc")),
    _sigma_ts(getADMaterialProperty<Real>("sigma_ts")),
    _sigma_cs(getADMaterialProperty<Real>("sigma_cs")),
    _delta_elp(getADMaterialProperty<Real>("delta_elp")),
    _model(getParam<std::string>("model")),
    _Ce_name(getParam<MaterialPropertyName>("Ce_name")),
    _Ce(declareADProperty<Real>(_Ce_name))
{
}

void
ADComputeExtraDrivingForcefromStress::initialSetup()
{
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("This function right now can only be used with "
               "isotropic elasticity tensor materials.");
}
void
ADComputeExtraDrivingForcefromStress::computeQpProperties()
{
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const ADReal k = lambda + 2.0 * mu / LIBMESH_DIM;
  const ADReal E = (mu * (3.0 * lambda + 2.0 * mu)) / (lambda + mu);
  const ADReal nu = lambda / (2.0 * (lambda + mu));
  // Remember v=1-d
  ADReal wts = _sigma_ts[_qp]*_sigma_ts[_qp]/2/E;
  ADReal sigma_hs = (2*_sigma_ts[_qp]*_sigma_cs[_qp])/ 
                    (3*std::max((_sigma_cs[_qp] - _sigma_ts[_qp]), 1e-9) );
  ADReal whs = sigma_hs*sigma_hs/2/k;

  ADReal beta1 = -1/sigma_hs*_delta_elp[_qp]*_gc[_qp]/8/_l[_qp]+2*whs/3/sigma_hs;


  ADReal beta2 = -(std::sqrt(3)*(3*sigma_hs-_sigma_ts[_qp])/sigma_hs/_sigma_ts[_qp] ) *_delta_elp[_qp]*_gc[_qp]/8/_l[_qp]
               - 2*whs/std::sqrt(3)/sigma_hs
               + 2*std::sqrt(3)*whs/_sigma_ts[_qp];

  // Method 1 I1 and J2 from mechanical strain
  //ADReal I1 = (3*lambda + 2*mu)*_D[_qp] * _strain[_qp].trace();
  //ADRankTwoTensor strain_dev = _strain[_qp].deviatoric(); 
  //ADReal J2_strain = strain_dev.doubleContraction(strain_dev);
  //ADReal J2 = 2.0 * mu * mu * _D[_qp]*_D[_qp] * J2_strain + 1e-14;

  // Method 2 I1 and J2 from damaged_stress
  ADReal I1 = _D[_qp] * _stress[_qp].trace();
  ADRankTwoTensor stress_dev = _stress[_qp].deviatoric(); 
  ADReal J2= 1/2*_D[_qp]*_D[_qp]* stress_dev.doubleContraction(stress_dev)+1e-14;
  
  if (_model == "finite_e")
  {
    mooseInfo("Using finite_e model");
    _Ce[_qp] = beta2 * std::sqrt(J2) + beta1 * I1
             + (1 / (_D[_qp] * std::sqrt(_D[_qp])))
             * (I1 < 0 ? 2 : 0) * (J2 * (1 + nu) / E + I1 * I1 * (1 - 2 * nu) / 6 / E);
  }
  else if (_model == "asymptotic")
  {
    mooseInfo("Using asymptotic model");
    _Ce[_qp] = beta2 * std::sqrt(J2) + beta1 * I1;
  }
  else
  {
    mooseError("Model is not specified properly, please use 'finite_e' or 'asymptotic'");
  }
}