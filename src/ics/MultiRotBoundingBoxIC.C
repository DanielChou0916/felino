#include "MultiRotBoundingBoxIC.h"
#include <Eigen/Geometry> // 用于四元数旋转

registerMooseObject("MooseApp", MultiRotBoundingBoxIC);

InputParameters
MultiRotBoundingBoxIC::validParams()
{
  InputParameters params = InitialCondition::validParams();

  // 使用向量定义参数
  params.addRequiredParam<std::vector<Real>>("cx", "The x coordinates of the centers of the boxes");
  params.addRequiredParam<std::vector<Real>>("cy", "The y coordinates of the centers of the boxes");
  params.addParam<std::vector<Real>>("cz", std::vector<Real>(), "The z coordinates of the centers of the boxes");

  params.addRequiredParam<std::vector<Real>>("lx", "The x lengths of the boxes");
  params.addRequiredParam<std::vector<Real>>("ly", "The y lengths of the boxes");
  params.addParam<std::vector<Real>>("lz", std::vector<Real>(), "The z lengths of the boxes");

  params.addParam<std::vector<Real>>("angle_x", std::vector<Real>(), "The rotation angles around the x-axis");
  params.addParam<std::vector<Real>>("angle_y", std::vector<Real>(), "The rotation angles around the y-axis");
  params.addParam<std::vector<Real>>("angle_z", std::vector<Real>(), "The rotation angles around the z-axis");

  params.addParam<std::vector<Real>>("inside", std::vector<Real>(), "The values inside the boxes");
  params.addParam<std::vector<Real>>("outside", std::vector<Real>(), "The values outside the boxes");

  params.addParam<std::vector<Real>>("int_width", std::vector<Real>(), "The widths of the diffuse interfaces");

  return params;
}

MultiRotBoundingBoxIC::MultiRotBoundingBoxIC(const InputParameters & parameters)
  : InitialCondition(parameters)
{
  const std::vector<Real> & cx = getParam<std::vector<Real>>("cx");
  const std::vector<Real> & cy = getParam<std::vector<Real>>("cy");
  const std::vector<Real> & cz = getParam<std::vector<Real>>("cz");
  const std::vector<Real> & lx = getParam<std::vector<Real>>("lx");
  const std::vector<Real> & ly = getParam<std::vector<Real>>("ly");
  const std::vector<Real> & lz = getParam<std::vector<Real>>("lz");

  const std::vector<Real> & angle_x = getParam<std::vector<Real>>("angle_x");
  const std::vector<Real> & angle_y = getParam<std::vector<Real>>("angle_y");
  const std::vector<Real> & angle_z = getParam<std::vector<Real>>("angle_z");

  const std::vector<Real> & inside = getParam<std::vector<Real>>("inside");
  const std::vector<Real> & outside = getParam<std::vector<Real>>("outside");
  const std::vector<Real> & int_width = getParam<std::vector<Real>>("int_width");

  unsigned int num_boxes = cx.size();  // 使用 cx 的长度确定盒子数量

  for (unsigned int i = 0; i < num_boxes; ++i)
  {
    Box box;
    box.cx = cx[i];
    box.cy = cy[i];
    box.cz = (cz.size() > i) ? cz[i] : 0.0;

    box.lx = lx[i];
    box.ly = ly[i];
    box.lz = (lz.size() > i) ? lz[i] : 0.0;

    box.angle_x = (angle_x.size() > i) ? angle_x[i] : 0.0;
    box.angle_y = (angle_y.size() > i) ? angle_y[i] : 0.0;
    box.angle_z = (angle_z.size() > i) ? angle_z[i] : 0.0;

    box.inside = (inside.size() > i) ? inside[i] : 0.0;
    box.outside = (outside.size() > i) ? outside[i] : 0.0;
    box.int_width = (int_width.size() > i) ? int_width[i] : 0.0;

    _boxes.push_back(box);
  }
}

Point
MultiRotBoundingBoxIC::rotatePoint(const Point & p, const Real angle_x, const Real angle_y, const Real angle_z, const Real cx, const Real cy, const Real cz) const
{
  // Convert angles to radians
  Real rad_x = angle_x * M_PI / 180.0;
  Real rad_y = angle_y * M_PI / 180.0;
  Real rad_z = angle_z * M_PI / 180.0;

  // Define quaternions for each rotation
  Eigen::Quaterniond qx(Eigen::AngleAxisd(rad_x, Eigen::Vector3d::UnitX()));
  Eigen::Quaterniond qy(Eigen::AngleAxisd(rad_y, Eigen::Vector3d::UnitY()));
  Eigen::Quaterniond qz(Eigen::AngleAxisd(rad_z, Eigen::Vector3d::UnitZ()));

  // Combine the rotations into one quaternion
  Eigen::Quaterniond q = qz * qy * qx;

  // Translate the point to the box center (cx, cy, cz)
  Eigen::Vector3d point(p(0) - cx, p(1) - cy, p(2) - cz);

  // Apply the rotation
  Eigen::Vector3d rotated_point = q.inverse() * point;

  // Translate the point back to the original coordinate system
  return Point(rotated_point.x() + cx, rotated_point.y() + cy, rotated_point.z() + cz);
}


Real
MultiRotBoundingBoxIC::value(const Point & p)
{
  Real result = 0.0;

  for (const auto & box : _boxes)
  {
    Real box_value = evaluateBox(box, p);
    result = std::max(result, box_value);  // ✅ 只取最大，不累加，防止 d > 1
  }

  return result;
}

Real
MultiRotBoundingBoxIC::evaluateBox(const Box & box, const Point & p) const
{
  // Now passing cx, cy, cz to rotatePoint for rotation around box center
  Point rotated_p = rotatePoint(p, box.angle_x, box.angle_y, box.angle_z, box.cx, box.cy, box.cz);

  Point bottom_left(box.cx - box.lx / 2, box.cy - box.ly / 2, box.cz - box.lz / 2);
  Point top_right(box.cx + box.lx / 2, box.cy + box.ly / 2, box.cz + box.lz / 2);

  if (box.int_width == 0.0)
  {
    for (unsigned int i = 0; i < Moose::dim; ++i)
      if (rotated_p(i) < bottom_left(i) || rotated_p(i) > top_right(i))
        return box.outside;

    return box.inside;
  }
  else
  {
    Real f_in = 1.0;
    for (unsigned int i = 0; i < Moose::dim; ++i)
      if (bottom_left(i) != top_right(i))
        f_in *= 0.5 * (std::tanh(2.0 * (rotated_p(i) - bottom_left(i)) / box.int_width) -
                       std::tanh(2.0 * (rotated_p(i) - top_right(i)) / box.int_width));

    return box.outside + (box.inside - box.outside) * f_in;
  }
}

