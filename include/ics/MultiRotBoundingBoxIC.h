#pragma once

#include "InitialCondition.h"
#include "libmesh/point.h"

class MultiRotBoundingBoxIC : public InitialCondition
{
public:
  static InputParameters validParams();
  MultiRotBoundingBoxIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:
  struct Box {
    Real cx, cy, cz;
    Real lx, ly, lz;
    Real angle_x, angle_y, angle_z;
    Real inside, outside;
    Real int_width;
  };

  std::vector<Box> _boxes;  // 存储多个盒子

  // 用于旋转点的辅助函数
  Point rotatePoint(const Point & p, const Real angle_x, const Real angle_y, const Real angle_z, const Real cx, const Real cy, const Real cz) const;

  // 评估每个盒子对当前点的贡献
  Real evaluateBox(const Box & box, const Point & p) const;
};
