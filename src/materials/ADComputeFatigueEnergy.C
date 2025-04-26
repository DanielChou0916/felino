//

#include "ADComputeFatigueEnergy.h"
#include "MathUtils.h"
#include <algorithm>
registerMooseObject("SolidMechanicsApp", ADComputeFatigueEnergy);

InputParameters
ADComputeFatigueEnergy::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>("uncracked_base_name", "",
    "Optional prefix for uncracked stress-strain properties. Leave empty if using default names.");
  params.addParam<bool>("finite_strain_model", false, "The model is using finite strain");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");//✅0318
  params.addParam<std::string>("type", "mean_load", 
      "Choose the energy accumulation: 'mean_load' or 'elastic_energy'");
  params.addParam<MaterialPropertyName>(
      "bar_psi_name", "current_energy", "Name of material property for fatigue energy");
  params.addParam<MaterialPropertyName>(
      "acc_bar_psi_name", "accumulated_energy", "Name of material property for fatigue energy");
  params.addParam<bool>(
      "multiply_by_D",
      false,
      "If true, multiply current energy by degradation function.");
  params.addParam<std::string>("accumulation_mode", "Fatigue", 
    "Choose accumulation mode: 'Monotonic', 'Fatigue', 'FatigueCLA'");
  params.addCoupledVar("N_cyc_variable", "Optional variable representing current N_cycle");
  return params;
}

ADComputeFatigueEnergy::ADComputeFatigueEnergy(const InputParameters & parameters)
  : Material(parameters), 
    GuaranteeConsumer(this),
    _uncracked_base_name(isParamValid("uncracked_base_name") ? getParam<std::string>("uncracked_base_name") + "_" : ""),
    _finite_strain_model(getParam<bool>("finite_strain_model")),
    _strain(
        _finite_strain_model
            ? getADMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "elastic_strain")
            : getADMaterialPropertyByName<RankTwoTensor>(_uncracked_base_name + "mechanical_strain")),
    _elasticity_tensor_name(_uncracked_base_name + "elasticity_tensor"),
    _elasticity_tensor(getADMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _D_name(getParam<MaterialPropertyName>("D_name")),// ✅
    _D(getADMaterialProperty<Real>(_D_name)),// ✅degradation function is an input, need to be created else where in input file, hence use getADMaterialProperty instead of declareADProperty
    _type(getParam<std::string>("type")),  // ✅ We will extend this later!
    _bar_psi_name(getParam<MaterialPropertyName>("bar_psi_name")),  // ✅ 
    _bar_psi(declareADProperty<Real>(_bar_psi_name)),  // ✅ 
    _bar_psi_old(getMaterialPropertyOld<Real>(_bar_psi_name)),

    _acc_bar_psi_name(getParam<MaterialPropertyName>("acc_bar_psi_name")),  // ✅ 
    _acc_bar_psi(declareADProperty<Real>(_acc_bar_psi_name)),  // ✅ 
    _acc_bar_psi_old(getMaterialPropertyOld<Real>(_acc_bar_psi_name)),

    _n(getADMaterialProperty<Real>("material_constant_n")),
    _R(getADMaterialProperty<Real>("load_ratio")),
    _multiply_by_degradation(getParam<bool>("multiply_by_D")),
    _accumulation_mode(getParam<std::string>("accumulation_mode"))
{
  if (parameters.isParamValid("N_cyc_variable"))
  {
    _N_cyc_var = &coupledValue("N_cyc_variable");
    _N_cyc_var_old = &coupledValueOld("N_cyc_variable");
  }
  else
  {
    _N_cyc_var = nullptr;
    _N_cyc_var_old = nullptr;
  }
}
void
ADComputeFatigueEnergy::initQpStatefulProperties()
{
  _acc_bar_psi[_qp] = _acc_bar_psi_old[_qp];
}
void
ADComputeFatigueEnergy::initialSetup()
{ // From Daniel: As the elastic modulus E is adopted in eq 24, I suppose this method is based on isotropy scenario.
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("This function right now can only be used with "
               "isotropic elasticity tensor materials.");
}
void
ADComputeFatigueEnergy::computeQpProperties()
{
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  //const ADReal k = lambda + 2.0 * mu / LIBMESH_DIM;
  const ADReal E = (mu * (3.0 * lambda + 2.0 * mu)) / (lambda + mu);
  //const ADReal nu = lambda / (2.0 * (lambda + mu));

  // Step 1 Spectral Decomposition on mechanical strain
  ADRankTwoTensor eigvec;
  std::vector<ADReal> eigval(LIBMESH_DIM);
  ADRankFourTensor Ppos =
                  _strain[_qp].positiveProjectionEigenDecomposition(eigval, eigvec);
  ADReal e_min = *std::min_element(eigval.begin(), eigval.end());
  ADReal e_max = *std::max_element(eigval.begin(), eigval.end());
      
  //ADReal R = e_min / (e_max + 1e-12); //1e-12 for numerical stability;

  // Step 2 Compute Fatigue Energy
  // Add conditions here to incoporate other fatigue energy formation(in the future!!)
  ADReal current_psi = 2 * E * e_max * e_max 
  * std::pow((1 + _R[_qp]) / 2, 2)
  * std::pow((1 - _R[_qp]) / 2, _n[_qp]);
  // 
  ADReal D_constant = _multiply_by_degradation ? _D[_qp] : 1.0;

  _bar_psi[_qp] = D_constant * current_psi;
  //
  // Accumulation rule
  //Step 3 Accumulate energy
  if (_accumulation_mode == "Monotonic")
  {
    //delta_energy = 0.0;
    _acc_bar_psi[_qp] = 0.0;
  }
  else if (_accumulation_mode == "Fatigue")
  {
    ADReal delta_energy = (_bar_psi[_qp] > _bar_psi_old[_qp]) ? (_bar_psi[_qp] - _bar_psi_old[_qp]) : 0.0;
    _acc_bar_psi[_qp] = _acc_bar_psi_old[_qp] + delta_energy;
  }
  else if (_accumulation_mode == "FatigueCLA")
  {
      if (!_N_cyc_var)
        mooseError("N_cyc_variable must be provided when using FatigueCLA mode!");
      else
      {
        Real delta_N = (*_N_cyc_var)[_qp] - (*_N_cyc_var_old)[_qp];
        if (delta_N < 0.0)
          delta_N = 0.0;
        ADReal delta_energy = _bar_psi[_qp] * delta_N;
        _acc_bar_psi[_qp] = _acc_bar_psi_old[_qp] + delta_energy;
        //_acc_bar_psi[_qp] = _bar_psi[_qp] * (*_N_cyc_var)[_qp];
      }
  }
  else
    mooseError("Unknown accumulation_mode: ", _accumulation_mode);
}