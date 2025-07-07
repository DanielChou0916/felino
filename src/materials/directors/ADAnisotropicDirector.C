#include "ADAnisotropicDirector.h"

registerMooseObject("PhaseFieldApp", ADAnisotropicDirector);

InputParameters
ADAnisotropicDirector::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Generates anisotropic tensor A = I + coef*(I - n ⊗ n)");
  params.addRequiredParam<RealVectorValue>("normal", "Manually specified weak plane normal");
  params.addParam<Real>("coef", 1.0, "Intensity of anisotropy: A = I + coef*(I - n ⊗ n)");
  params.addParam<MaterialPropertyName>("output_name", "A", "Name of output RankTwoTensor A");
  return params;
}

ADAnisotropicDirector::ADAnisotropicDirector(const InputParameters & parameters)
  : Material(parameters),
    _normal(getParam<RealVectorValue>("normal")),
    _coef(getParam<Real>("coef")),
    _output_name(getParam<MaterialPropertyName>("output_name")),  // ✅
    _directional_tensor(declareADProperty<RankTwoTensor>(_output_name))
{
}

void
ADAnisotropicDirector::computeQpProperties()
{
  const ADRankTwoTensor I = ADRankTwoTensor::initIdentity;
  const ADRankTwoTensor n_outer = ADRankTwoTensor::selfOuterProduct(_normal.unit());
  _directional_tensor[_qp] = I + _coef * (I - n_outer);
}
