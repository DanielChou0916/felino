#include "ADAnisotropicDirector.h"

registerMooseObject("PhaseFieldApp", ADAnisotropicDirector);

InputParameters
ADAnisotropicDirector::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Generates anisotropic tensor A = I + coef*(I - n ⊗ n)");

  MooseEnum input_type("normal_vector xy_angle", "normal_vector");
  params.addParam<MooseEnum>("input_type", input_type, 
    "Input mode: normal_vector or xy_angle");

  MooseEnum normalize_type("factorial_norm trace_norm det_norm none", "factorial_norm");
  params.addParam<MooseEnum>("normalize_director", normalize_type, 
    "Normalization mode: trace_norm, det_norm, factorial_norm, or none");

  params.addParam<RealVectorValue>("normal", RealVectorValue(1, 0, 0),
                                   "Weak plane normal (used if input_type = normal_vector)");
  params.addParam<Real>("xy_angle_deg", 0.0,
                        "Angle from x-axis in degrees (used if input_type = xy_angle)");
  params.addParam<Real>("coef", 1.0, "Intensity of anisotropy: A = I + coef*(I - n ⊗ n)");
  params.addParam<Real>("factor", 0.01, "Factor for manual scaling A.");
  params.addParam<MaterialPropertyName>("output_name", "A", "Name of output RankTwoTensor A");
  return params;
}

ADAnisotropicDirector::ADAnisotropicDirector(const InputParameters & parameters)
  : Material(parameters),
    _coef(getParam<Real>("coef")),
    _output_name(getParam<MaterialPropertyName>("output_name")),
    _directional_tensor(declareADProperty<RankTwoTensor>(_output_name)),
    _input_type(getParam<MooseEnum>("input_type").getEnum<input_type>()),
    _normalize_director(getParam<MooseEnum>("normalize_director").getEnum<norm_type>()),
    _factor(getParam<Real>("factor"))
{
  if (_input_type == input_type::xy_angle)
  { 
    if (_mesh.dimension() == 2)
    {
      const Real angle_deg = getParam<Real>("xy_angle_deg");
      const Real theta = (angle_deg + 90.0) * libMesh::pi / 180.0;
      _normal = RealVectorValue(std::cos(theta), std::sin(theta), 0.0);
      mooseInfo("xy_angle as input.");
    }
    else mooseError("Angle input of only for 2D simulation!");
  }
  else
    _normal = getParam<RealVectorValue>("normal");
}

void
ADAnisotropicDirector::computeQpProperties()
{
  const ADRankTwoTensor I = ADRankTwoTensor::initIdentity;
  const ADRankTwoTensor n_outer = ADRankTwoTensor::selfOuterProduct(_normal.unit());

  ADRankTwoTensor A = I + _coef * (I - n_outer);
  if (_normalize_director == norm_type::trace_norm)
  {
    mooseInfo("Ensuring trace of A = DIM.");
    const ADReal tr = A.trace();
    A /= (tr / LIBMESH_DIM); // normalize so that trace(A) == LIBMESH_DIM
  } 
  else if (_normalize_director == norm_type::det_norm)
  {
    mooseInfo("Ensuring determinant of A = 1.");
    const ADReal detA = A.det();
    A /= std::pow(detA, 1.0 / LIBMESH_DIM);
  } 
  else if (_normalize_director == norm_type::factorial_norm)
  {
    mooseInfo("Manual scaling is used.");
    A /= _factor;
  } 

  _directional_tensor[_qp] = A;
}
