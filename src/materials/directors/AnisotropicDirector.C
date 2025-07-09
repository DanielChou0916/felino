#include "AnisotropicDirector.h"

registerMooseObject("PhaseFieldApp", AnisotropicDirector);

InputParameters
AnisotropicDirector::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Generates anisotropic tensor A = I + coef*(I - n ⊗ n)");

  MooseEnum input_type("normal_vector xy_angle", "normal_vector");
  params.addParam<MooseEnum>("input_type", input_type, "Input mode: normal_vector or xy_angle");

  params.addParam<RealVectorValue>("normal", RealVectorValue(1, 0, 0),
                                   "Weak plane normal (used if input_type = normal_vector)");
  params.addParam<Real>("xy_angle_deg", 0.0,
                        "Angle from x-axis in degrees (used if input_type = xy_angle)");
  params.addParam<Real>("coef", 1.0, "Intensity of anisotropy: A = I + coef*(I - n ⊗ n)");
  params.addParam<MaterialPropertyName>("output_name", "A", "Name of output RankTwoTensor A");
  params.addParam<bool>("normalize_director", false, "If true, normalize anisotropic tensor A by its trace.");
  return params;
}

AnisotropicDirector::AnisotropicDirector(const InputParameters & parameters)
  : Material(parameters),
    _coef(getParam<Real>("coef")),
    _output_name(getParam<MaterialPropertyName>("output_name")),
    _directional_tensor(declareProperty<RankTwoTensor>(_output_name)),
    _input_type(getParam<MooseEnum>("input_type").getEnum<input_type>()),
    _normalize_director(getParam<bool>("normalize_director"))
{
  if (_input_type == input_type::xy_angle)
  {
    const Real angle_deg = getParam<Real>("xy_angle_deg");
    const Real theta = (angle_deg + 90) * libMesh::pi / 180.0;
    _normal = RealVectorValue(std::cos(theta), std::sin(theta), 0.0);
  }
  else
    _normal = getParam<RealVectorValue>("normal");
}

void
AnisotropicDirector::computeQpProperties()
{
  const RankTwoTensor I = RankTwoTensor::initIdentity;
  const RankTwoTensor n_outer = RankTwoTensor::selfOuterProduct(_normal.unit());
  RankTwoTensor A = I + _coef * (I - n_outer);
  if (_normalize_director)
  {
    const Real tr = A.trace();
    A /= (tr / 3.0); // normalize so that trace(A) == 3
  } 
  _directional_tensor[_qp] = A;
}
