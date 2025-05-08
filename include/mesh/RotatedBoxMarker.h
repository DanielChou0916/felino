#pragma once

#include "Marker.h"
#include <Eigen/Geometry>

class RotatedBoxMarker : public Marker
{
public:
  static InputParameters validParams();
  RotatedBoxMarker(const InputParameters & parameters);

protected:
  virtual MarkerValue computeElementMarker() override;

  // 幫助函數：將世界座標點轉為旋轉後的 box 區域坐標
  Point rotatePoint(const Point & p) const;

  MarkerValue _inside;
  MarkerValue _outside;

  // box 定義參數
  Real _cx, _cy, _cz;
  Real _lx, _ly, _lz;

  // 旋轉角度（以度為單位）
  Real _angle_x, _angle_y, _angle_z;
};
