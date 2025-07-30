#include "PorousFlowPermeabilityTensorFromVarPFF.h"

registerMooseObject("PorousFlowApp", PorousFlowPermeabilityTensorFromVarPFF);
//registerMooseObject("PorousFlowApp", ADPorousFlowPermeabilityTensorFromVarPFF);

template <bool is_ad>
InputParameters
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::validParams()
{
  InputParameters params = PorousFlowPermeabilityBase::validParams();
  params.addRequiredCoupledVar("perm", "The scalar permeability variable.");
  params.addParam<RealTensorValue>(
      "k_anisotropy",
      "Anisotropy tensor to multiply the scalar permeability. "
      "Defaults to identity tensor for isotropic permeability.");
  params.addParam<MaterialPropertyName>(
      "crack_direction_name", "direction_matrix",
      "Material property name for the crack direction tensor.");
  params.addParam<MaterialPropertyName>(
      "coef_name", "A", "Material property name for the coefficient.");
  params.addParam<std::string>(
      "model", "directional",
      "Choose the model: 'directional' or 'exponential'.");
  params.addClassDescription(
      "Computes permeability tensor from a coupled scalar variable and an anisotropy tensor. "
      "An additional term can be added to account for fracture plane permeability changes.");
  return params;
}

template <bool is_ad>
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::PorousFlowPermeabilityTensorFromVarPFFTempl(
    const InputParameters & parameters)
  : PorousFlowPermeabilityBaseTempl<is_ad>(parameters),
    _perm(coupledValue("perm")),
    _k_anisotropy(parameters.isParamValid("k_anisotropy")
                      ? this->template getParam<RealTensorValue>("k_anisotropy")
                      : RealTensorValue(1.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0)),
    _crack_direction_name(this->template getParam<MaterialPropertyName>("crack_direction_name")),
    _coef_name(this->template getParam<MaterialPropertyName>("coef_name")),
    _model(this->template getParam<std::string>("model"))
{
  if (_model == "directional")
  {
    bool has_direction = parameters.isParamValid("crack_direction_name");
    bool has_coef = parameters.isParamValid("coef_name");

    if (has_direction && has_coef)
    {
      _crack_direction =
          &this->template getMaterialProperty<RankTwoTensor>(_crack_direction_name);
      _A = &this->template getMaterialProperty<Real>(_coef_name);
    }
    else
    {
      mooseWarning("Directional permeability model selected but required parameters are missing. "
                   "Switching to exponential permeability model.");
      _model = "exponential";
      _A = &this->template getMaterialProperty<Real>(_coef_name);
    }
  }
  else if (_model == "exponential")
  {
    _A = &this->template getMaterialProperty<Real>(_coef_name);
  }
  else
  {
    mooseWarning("Unknown permeability model '" + _model + "'. Defaulting to exponential permeability model.");
    _model = "exponential";
    _A = &this->template getMaterialProperty<Real>(_coef_name);
  }

  mooseInfo("Permeability model = " + _model);
}

template <bool is_ad>
void
PorousFlowPermeabilityTensorFromVarPFFTempl<is_ad>::computeQpProperties()
{
  RankTwoTensor I2(RankTwoTensor::initIdentity);

  if (_model == "directional")
    _permeability_qp[_qp] =
        _k_anisotropy * _perm[_qp] + (*_A)[_qp] * (I2 - (*_crack_direction)[_qp]);
  else
    _permeability_qp[_qp] =
        _k_anisotropy * _perm[_qp] * std::exp((*_A)[_qp]);

  if (!is_ad)
  {
    (*_dpermeability_qp_dvar)[_qp].assign(_num_var, RealTensorValue());
    (*_dpermeability_qp_dgradvar)[_qp].resize(LIBMESH_DIM);
    for (const auto i : make_range(Moose::dim))
      (*_dpermeability_qp_dgradvar)[_qp][i].assign(_num_var, RealTensorValue());
  }
}

template class PorousFlowPermeabilityTensorFromVarPFFTempl<false>;
//template class PorousFlowPermeabilityTensorFromVarPFFTempl<true>;
