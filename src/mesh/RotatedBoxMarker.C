#include "RotatedBoxMarker.h"
#include "MooseUtils.h"

registerMooseObject("MooseApp", RotatedBoxMarker);

InputParameters
RotatedBoxMarker::validParams()
{
  InputParameters params = Marker::validParams();

  params.addRequiredParam<Real>("cx", "Center x-coordinate of the box");
  params.addRequiredParam<Real>("cy", "Center y-coordinate of the box");
  params.addRequiredParam<Real>("cz", "Center z-coordinate of the box");

  params.addRequiredParam<Real>("lx", "Length of the box in x-direction");
  params.addRequiredParam<Real>("ly", "Length of the box in y-direction");
  params.addRequiredParam<Real>("lz", "Length of the box in z-direction");

  params.addParam<Real>("angle_x", 0.0, "Rotation angle around x-axis in degrees");
  params.addParam<Real>("angle_y", 0.0, "Rotation angle around y-axis in degrees");
  params.addParam<Real>("angle_z", 0.0, "Rotation angle around z-axis in degrees");

  MooseEnum marker_states = Marker::markerStates();
  params.addRequiredParam<MooseEnum>("inside", marker_states, "How to mark elements inside the box.");
  params.addRequiredParam<MooseEnum>("outside", marker_states, "How to mark elements outside the box.");

  params.addClassDescription("Marks a rotated box region for refinement or coarsening using center/size/angle definition.");
  return params;
}

RotatedBoxMarker::RotatedBoxMarker(const InputParameters & parameters)
  : Marker(parameters),
    _inside(parameters.get<MooseEnum>("inside").getEnum<MarkerValue>()),
    _outside(parameters.get<MooseEnum>("outside").getEnum<MarkerValue>()),
    _cx(parameters.get<Real>("cx")),
    _cy(parameters.get<Real>("cy")),
    _cz(parameters.get<Real>("cz")),
    _lx(parameters.get<Real>("lx")),
    _ly(parameters.get<Real>("ly")),
    _lz(parameters.get<Real>("lz")),
    _angle_x(parameters.get<Real>("angle_x")),
    _angle_y(parameters.get<Real>("angle_y")),
    _angle_z(parameters.get<Real>("angle_z"))
{
}

Point
RotatedBoxMarker::rotatePoint(const Point & p) const
{
  Real rad_x = _angle_x * libMesh::pi / 180.0;
  Real rad_y = _angle_y * libMesh::pi / 180.0;
  Real rad_z = _angle_z * libMesh::pi / 180.0;

  Eigen::Quaterniond qx(Eigen::AngleAxisd(rad_x, Eigen::Vector3d::UnitX()));
  Eigen::Quaterniond qy(Eigen::AngleAxisd(rad_y, Eigen::Vector3d::UnitY()));
  Eigen::Quaterniond qz(Eigen::AngleAxisd(rad_z, Eigen::Vector3d::UnitZ()));

  Eigen::Quaterniond q = qz * qy * qx;
  Eigen::Vector3d vec(p(0) - _cx, p(1) - _cy, p(2) - _cz);
  Eigen::Vector3d rotated = q.inverse() * vec;

  return Point(rotated.x() + _cx, rotated.y() + _cy, rotated.z() + _cz);
}

Marker::MarkerValue
RotatedBoxMarker::computeElementMarker()
{
  const Point centroid = _current_elem->vertex_average();
  const Point local = rotatePoint(centroid);

  const Real dx = std::abs(local(0) - _cx);
  const Real dy = std::abs(local(1) - _cy);
  const Real dz = std::abs(local(2) - _cz);

  if (dx <= _lx / 2.0 && dy <= _ly / 2.0 && dz <= _lz / 2.0)
    return _inside;

  return _outside;
}
