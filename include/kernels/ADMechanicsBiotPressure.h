// ADMechanicsBiotPressure.h
#pragma once
#include "ADKernel.h"

class ADMechanicsBiotPressure : public ADKernel
{
public:
  static InputParameters validParams();
  ADMechanicsBiotPressure(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<Real> & _alpha;   // biot_coefficient (AD material)
  const VariableValue & _p_eff;            // lagged p_eff Aux variable (coupled)
  const unsigned int _component;             // 0:x, 1:y, 2:z
  const bool _rz;
};
