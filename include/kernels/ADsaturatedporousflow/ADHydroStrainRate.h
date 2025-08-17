#pragma once
#include "Material.h"
#include "RankTwoTensor.h"

class ADHydroStrainRate : public Material
{
public:
  static InputParameters validParams();
  ADHydroStrainRate(const InputParameters & params);

protected:
  void computeQpProperties() override;

  // 與 TensorMechanics 對齊的 base_name（可留空）
  const std::string _base_name;

  // 直接讀 TensorMechanics 提供的總應變（現在步）：AD，讓導數自動傳遞
  const ADMaterialProperty<RankTwoTensor> & _total_strain;

  // 舊步總應變：用 non-AD Old property 就好
  const MaterialProperty<RankTwoTensor> & _total_strain_old;

  // 我們要輸出的體積應變率 ε_v_dot
  ADMaterialProperty<Real> & _vol_strain_rate;
};
