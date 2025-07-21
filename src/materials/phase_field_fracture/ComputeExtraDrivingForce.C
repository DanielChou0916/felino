//

#include "ComputeExtraDrivingForce.h"
#include "MathUtils.h"

registerMooseObject("SolidMechanicsApp", ComputeExtraDrivingForce);

InputParameters
ComputeExtraDrivingForce::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("uncracked_base_name", "",
                                "Optional parameter that allows the user to define "
                                "multiple mechanics material systems on the same "
                                "block, i.e. for multiple phases");
  params.addParam<bool>("finite_strain_model", false, "The model is using finite strain");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");//✅0318
  params.addParam<std::string>("model", "finite_e", 
      "Choose the energy model: 'finite_e' or 'asymptotic'");
  params.addParam<MaterialPropertyName>(
      "Ce_name", "driving_force", "Name of material property for extra driving force");
  //params.addParam<bool>("output_Ce_aux",false,"If true, output Ce as an AuxVariable.");
  // Add new parameters here
  return params;
}

ComputeExtraDrivingForce::ComputeExtraDrivingForce(const InputParameters & parameters)
  : Material(parameters), 
    GuaranteeConsumer(this),
    _uncracked_base_name(getParam<std::string>("uncracked_base_name").empty() ? "" : getParam<std::string>("uncracked_base_name") + "_"),
    _finite_strain_model(getParam<bool>("finite_strain_model")),
    _strain(
        _finite_strain_model
            ? getMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "elastic_strain")
            : getMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "mechanical_strain")),
    _elasticity_tensor_name(_uncracked_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _D_name(getParam<MaterialPropertyName>("D_name")),// ✅
    _D(getMaterialProperty<Real>(_D_name)),// ✅degradation function is an input, need to be created else where in input file, hence use getMaterialProperty instead of declareProperty
    _l(getMaterialProperty<Real>("l")),
    _gc(getMaterialProperty<Real>("gc")),
    _sigma_ts(getMaterialProperty<Real>("sigma_ts")),
    _sigma_cs(getMaterialProperty<Real>("sigma_cs")),
    _delta_elp(getMaterialProperty<Real>("delta_elp")),
    _model(getParam<std::string>("model")),  // ✅
    _Ce_name(getParam<MaterialPropertyName>("Ce_name")),  // ✅ 
    _Ce(declareProperty<Real>(_Ce_name))  // ✅ 
    //_output_Ce_aux(getParam<bool>("output_Ce_aux"))
{
}

void
ComputeExtraDrivingForce::initialSetup()
{
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("This function right now can only be used with "
               "isotropic elasticity tensor materials.");
}
void
ComputeExtraDrivingForce::computeQpProperties()
{
  const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const Real k = lambda + 2.0 * mu / LIBMESH_DIM;
  const Real E = (mu * (3.0 * lambda + 2.0 * mu)) / (lambda + mu);
  const Real nu = lambda / (2.0 * (lambda + mu));

  // Remember v=1-d
  Real wts = _sigma_ts[_qp]*_sigma_ts[_qp]/2/E;
  Real sigma_hs = (2*_sigma_ts[_qp]*_sigma_cs[_qp])/ 
                    (3*std::max((_sigma_cs[_qp] - _sigma_ts[_qp]), 1e-9) );
  Real whs = sigma_hs*sigma_hs/2/k;

  Real beta1 = -1/sigma_hs*_delta_elp[_qp]*_gc[_qp]/8/_l[_qp]+2*whs/3/sigma_hs;


  Real beta2 = -(std::sqrt(3)*(3*sigma_hs-_sigma_ts[_qp])/sigma_hs/_sigma_ts[_qp] ) *_delta_elp[_qp]*_gc[_qp]/8/_l[_qp]
               - 2*whs/std::sqrt(3)/sigma_hs
               + 2*std::sqrt(3)*whs/_sigma_ts[_qp];

  // Method 1 I1 and J2 from mechanical strain
  Real I1 = (3*lambda + 2*mu)*_D[_qp] * _strain[_qp].trace();
  RankTwoTensor strain_dev = _strain[_qp].deviatoric(); 
  Real J2_strain = strain_dev.doubleContraction(strain_dev);
  Real J2 = 2.0 * mu * mu * _D[_qp]*_D[_qp] * J2_strain + 1e-14;
  //
  // Method 2 I1 and J2 from damaged_stress
  //RankTwoTensor stress = _elasticity_tensor[_qp]*_strain[_qp];
  //Real I1 = _D[_qp] * stress.trace();
  //RankTwoTensor stress_dev = stress.deviatoric(); 
  //Real J2_tri = 1/2*_D[_qp]*_D[_qp]* stress_dev.doubleContraction(stress_dev);
  //Real J2;
  //if (J2_tri<1e-8)
  //{
  //  J2 = 0 ;
  //}
  //else
  //{
  //  J2 = J2_tri;
  //}
  //
  //Method 3: get Stress property from domain
  // 1. need to include more header file to access cracked stress.
  // 2. need to ensure this function is working AFTER cracked stress computation 
  //RankTwoTensor cracked_stress = stress[_qp]; //Extract cracked stress from domain
  //Real I1 = cracked_stress.trace();
  //RankTwoTensor cracked_stress_dev = cracked_stress.deviatoric(); 
  //Real J2 = 1/2* cracked_stress_dev.doubleContraction(cracked_stress_dev);
  
  
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