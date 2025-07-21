

#include "PorousFlowPermeabilityTensorFromVarPFF.h"

registerMooseObject("PorousFlowApp", PorousFlowPermeabilityTensorFromVarPFF);
//registerMooseObject("PorousFlowApp", ADPorousFlowPermeabilityTensorFromVarPFF);

template <bool is_ad>
InputParameters
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addRequiredCoupledVar("perm", "The scalar permeability");
  params.addParam<RealTensorValue>("k_anisotropy",
                                   "A tensor to multiply the scalar "
                                   "permeability, in order to obtain anisotropy if "
                                   "required. Defaults to isotropic permeability "
                                   "if not specified.");
  params.addParam<MaterialPropertyName>(
      "crack_direction_name", "direction_matrix", "Name of material property for crack direction matrix.");// ✅
  params.addParam<MaterialPropertyName>(
      "coef_name", "A", "Name of coefficient");// ✅
  params.addParam<std::string>("model", "directional", 
        "Choose the model: 'directional' or 'exponential'");
  params.addClassDescription(
      "This Material calculates the permeability tensor from a coupled variable "
      "multiplied by a tensor, additional term is adding to account change of permeability of fracture plane.");
  return params;
}

template <bool is_ad>
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::PorousFlowPermeabilityTensorFromVarPFFTempl(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBaseTempl<is_ad>(parameters),
    _perm(coupledValue("perm")),
    _k_anisotropy(parameters.isParamValid("k_anisotropy")
                      ? this->template getParam<RealTensorValue>("k_anisotropy")
                      : RealTensorValue(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)),
  _crack_direction_name(this->template getParam<MaterialPropertyName>("crack_direction_name")),// ✅
  _crack_direction(this->template getMaterialProperty<RankTwoTensor>(_crack_direction_name)),// ✅
  _coef_name(this->template getParam<MaterialPropertyName>("coef_name")),// ✅
  _A(this->template getMaterialProperty<Real>(_coef_name)),  // ✅ 
  _model(this->template getParam<std::string>("model"))// ✅ 
{
}

template <bool is_ad>
void
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::computeQpProperties()
{ 
  RankTwoTensor I2(RankTwoTensor::initIdentity);// ✅
  if (_model == "directional")
  {  
    mooseInfo("Using directional permeability function");
    _permeability_qp[_qp] = _k_anisotropy * _perm[_qp] + _A[_qp]*(I2-_crack_direction[_qp]);
  }
  else if (_model == "exponential")
  {
    mooseInfo("Using exponential permeability function");
    _permeability_qp[_qp] = _k_anisotropy * _perm[_qp] * std::exp(_A[_qp]);
  }

  if (!is_ad)
  {
    (*_dpermeability_qp_dvar)[_qp].resize(_num_var, RealTensorValue());
    (*_dpermeability_qp_dgradvar)[_qp].resize(LIBMESH_DIM);

    for (const auto i : make_range(Moose::dim))
      (*_dpermeability_qp_dgradvar)[_qp][i].resize(_num_var, RealTensorValue());
  }
}

template class PorousFlowPermeabilityTensorFromVarPFFTempl<false>;
//template class PorousFlowPermeabilityTensorFromVarPFFTempl<true>;
