#include "PorousFlowPermeabilityTensorFromVarPFF.h"

registerMooseObject("PorousFlowApp", PorousFlowPermeabilityTensorFromVarPFF);
registerMooseObject("PorousFlowApp", ADPorousFlowPermeabilityTensorFromVarPFF);

template <bool is_ad>
InputParameters
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addRequiredCoupledVar("perm", "The scalar permeability variable.");
  params.addParam<RealTensorValue>("k_anisotropy", RealTensorValue(1,0,0, 0,1,0, 0,0,1),
                                   "Anisotropy tensor to multiply the scalar permeability.");
  params.addParam<MaterialPropertyName>("crack_direction_name", "direction_matrix",
                                        "Material property name for the crack direction tensor.");
  params.addParam<MaterialPropertyName>("coef_name", "A",
                                        "Material property name for the coefficient.");
  params.addParam<std::string>("model", "directional", "directional or exponential");
  params.addParam<MaterialPropertyName>(
    "permeability_tensor_name",
    "permeability",
    "Name of the permeability tensor property to declare/expose");
  return params;
}

template <bool is_ad>
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::
PorousFlowPermeabilityTensorFromVarPFFTempl(const InputParameters & parameters)
  : PorousFlowPermeabilityBaseTempl<is_ad>(parameters),
    _perm(this->template coupledGenericValue<is_ad>("perm")),
    _k_anisotropy(parameters.isParamValid("k_anisotropy")
                    ? this->template getParam<RealTensorValue>("k_anisotropy")
                    : RealTensorValue(1,0,0, 0,1,0, 0,0,1)),
    _crack_direction_name(this->template getParam<MaterialPropertyName>("crack_direction_name")),
    _coef_name(this->template getParam<MaterialPropertyName>("coef_name")),
    _model(this->template getParam<std::string>("model")),
    // 這個一定有，用「參考」一次到位（注意模板參數順序 <T, is_ad>）
    _A(this->template getGenericMaterialProperty<Real, is_ad>(_coef_name)),
    // 先設成 nullptr，等會兒看模型要不要抓
    _crack_direction(nullptr),
    _perm_out(this->template declareGenericProperty<RankTwoTensor, is_ad>(
                this->template getParam<MaterialPropertyName>("permeability_tensor_name")))
{
  if (_model == "directional")
  {
    // 只有 directional 才需要方向張量
    _crack_direction =
      &this->template getGenericMaterialProperty<RankTwoTensor, is_ad>(_crack_direction_name);
  }
  else if (_model != "exponential")
  {
    mooseWarning("Unknown permeability model '" + _model + "'. Defaulting to 'exponential'.");
    _model = "exponential";
  }

  mooseInfo("Permeability model = " + _model);
}

template <bool is_ad>
void
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::computeQpProperties()
{
  const RankTwoTensor I2(RankTwoTensor::initIdentity);

  if (_model == "directional")
    this->_permeability_qp[_qp] =
        _k_anisotropy * _perm[_qp] + _A[_qp] * (I2 - (*_crack_direction)[_qp]);
  else // exponential
    this->_permeability_qp[_qp] =
        _k_anisotropy * _perm[_qp] * std::exp(_A[_qp]);

  _perm_out[_qp]              = this->_permeability_qp[_qp];

  
  if (!is_ad)
  {
    (*this->_dpermeability_qp_dvar)[_qp].assign(this->_num_var, RealTensorValue());
    (*this->_dpermeability_qp_dgradvar)[_qp].resize(LIBMESH_DIM);
    for (const auto i : make_range(Moose::dim))
      (*this->_dpermeability_qp_dgradvar)[_qp][i].assign(this->_num_var, RealTensorValue());
  }
}

// 明確實例化
template class PorousFlowPermeabilityTensorFromVarPFFTempl<false>;
template class PorousFlowPermeabilityTensorFromVarPFFTempl<true>;
