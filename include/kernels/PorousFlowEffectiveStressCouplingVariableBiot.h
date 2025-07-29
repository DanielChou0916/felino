#pragma once

#include "Kernel.h"
#include "PorousFlowDictator.h"

/**
 * PorousFlowEffectiveStressCouplingVariableBiot computes
 * -coefficient*effective_porepressure*grad_component(test)
 * where component is the spatial component (not
 * a fluid component!)
 */
class PorousFlowEffectiveStressCouplingVariableBiot : public Kernel
{
public:
  static InputParameters validParams();

  PorousFlowEffectiveStressCouplingVariableBiot(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// The PorousFlow dictator that holds global info about the simulation
  const PorousFlowDictator & _dictator;

  /// Biot coefficient
  //const Real _coefficient;
  const MaterialPropertyName _biot_coef_name;// ✅0402
  const MaterialProperty<Real> & _coefficient;// ✅0402
  /// The spatial component
  const unsigned int _component;

  /// Effective porepressure
  const MaterialProperty<Real> & _pf;

  /// d(effective porepressure)/(d porflow variable)
  const MaterialProperty<std::vector<Real>> & _dpf_dvar;

  /// Whether an RZ coordinate system is being used
  const bool _rz;
};
