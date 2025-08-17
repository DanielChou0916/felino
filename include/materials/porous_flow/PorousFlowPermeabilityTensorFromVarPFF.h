#include "PorousFlowPermeabilityBase.h"

// TEMPLATED Material
template <bool is_ad>
class PorousFlowPermeabilityTensorFromVarPFFTempl
  : public PorousFlowPermeabilityBaseTempl<is_ad>
{
public:
  static InputParameters validParams();
  PorousFlowPermeabilityTensorFromVarPFFTempl(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  // 變數（AD/非AD 通吃）
  const GenericVariableValue<is_ad> & _perm;

  // 常數各向異性張量
  const RealTensorValue _k_anisotropy;

  // 參數名稱
  const MaterialPropertyName _crack_direction_name;
  const MaterialPropertyName _coef_name;

  // 模型
  std::string _model;

  // 這個「係數 A」無論哪個模型都會用到 => 用「參考」
  const GenericMaterialProperty<Real,         is_ad> & _A;

  // crack 方向只在 directional 模型需要 => 用「指標」，沒用就 nullptr
  const GenericMaterialProperty<RankTwoTensor, is_ad> * _crack_direction;
  // 對外可讀取的同名/自訂名 property
  GenericMaterialProperty<RankTwoTensor, is_ad> & _perm_out;

  usingPorousFlowPermeabilityBaseMembers;
};

// 型別別名 + 註冊
using PorousFlowPermeabilityTensorFromVarPFF  = PorousFlowPermeabilityTensorFromVarPFFTempl<false>;
using ADPorousFlowPermeabilityTensorFromVarPFF = PorousFlowPermeabilityTensorFromVarPFFTempl<true>;
