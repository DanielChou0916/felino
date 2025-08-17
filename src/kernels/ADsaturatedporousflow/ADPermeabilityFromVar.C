#include "ADPermeabilityFromVar.h"

registerMooseObject("MooseApp", ADPermeabilityFromVar);

InputParameters
ADPermeabilityFromVar::validParams()
{
  InputParameters params = Material::validParams();

  // required AD-coupled scalar permeability k0 (dimension of scalar permeability)
  // validParams()
  params.addParam<MaterialPropertyName>("perm_property", "perm",
    "Name of the scalar permeability material property");


  // anisotropy multiplier (default = identity)
  params.addParam<RealTensorValue>("k_anisotropy",
                                   RealTensorValue(1,0,0, 0,1,0, 0,0,1),
                                   "Anisotropy tensor K_aniso multiplying the scalar perm.");

  // model selection
  params.addParam<std::string>("model", "exponential",
                               "Permeability model: 'directional' or 'exponential'.");

  // optional properties/consts used by models
  params.addParam<MaterialPropertyName>("crack_direction_property",
                                        "",
                                        "AD property name for direction projection tensor (P = n⊗n). "
                                        "Required if model='directional'.");

  params.addParam<MaterialPropertyName>("coef_property",
                                        "",
                                        "AD property name for coefficient A (optional). "
                                        "Used by both models if provided.");

  params.addParam<Real>("coef", 0.0, "Constant A (used if 'coef_property' is not provided).");

  // output property name
  params.addParam<MaterialPropertyName>("permeability_tensor_name",
                                        "permeability_tensor",
                                        "Name of the AD permeability tensor property to declare/expose.");

  params.addClassDescription("Pure-AD permeability tensor from a scalar coupled variable.");
  return params;
}

ADPermeabilityFromVar::ADPermeabilityFromVar(const InputParameters & params)
  : Material(params),
    _perm(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("perm_property"))),
    _k_aniso(getParam<RealTensorValue>("k_anisotropy")),
    _use_directional(getParam<std::string>("model") == "directional"),
    _use_exponential(getParam<std::string>("model") == "exponential"),
    _A_const(getParam<Real>("coef")),
    _k_tensor(declareADProperty<RankTwoTensor>(
                getParam<MaterialPropertyName>("permeability_tensor_name")))
{
  // model sanity
  if (!_use_directional && !_use_exponential)
    paramError("model", "Must be 'directional' or 'exponential'.");

  // optional A as property
  if (isParamValid("coef_property"))
  {
    const auto & aname = getParam<MaterialPropertyName>("coef_property");
    if (!aname.empty())
      _A_prop = &getADMaterialProperty<Real>(aname);
  }

  // directional requires projection tensor property
  if (_use_directional)
  {
    if (!isParamValid("crack_direction_property") ||
        getParam<MaterialPropertyName>("crack_direction_property").empty())
      paramError("crack_direction_property",
                 "Required when model='directional' (projection tensor P = n⊗n).");

    const auto & pname = getParam<MaterialPropertyName>("crack_direction_property");
    _direction = &getADMaterialProperty<RankTwoTensor>(pname);
  }
}

void
ADPermeabilityFromVar::computeQpProperties()
{
  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  if (_use_directional)
  {
    const ADReal A = _A_prop ? (*_A_prop)[_qp] : _A_const; // AD if prop, else lifts const
    const ADRankTwoTensor & P = (*_direction)[_qp];          // projection tensor
    _k_tensor[_qp] = _k_aniso * _perm[_qp] + A * (I2 - P);
  }
  else // exponential
  {
    const ADReal A = _A_prop ? (*_A_prop)[_qp] : _A_const;
    _k_tensor[_qp] = _k_aniso * _perm[_qp] * exp(A);
  }
}
