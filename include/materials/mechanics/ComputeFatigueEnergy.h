//

#pragma once

#include "Material.h"
#include "GuaranteeConsumer.h"
class ComputeFatigueEnergy : public Material,
                                   public GuaranteeConsumer
{
public:
  static InputParameters validParams();

  ComputeFatigueEnergy(const InputParameters & parameters);
  void initialSetup() override;
protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
/// Additional: expand to finite strain
  /// Base name of the uncracked stress and strain
  const std::string _uncracked_base_name;
  /// Indicator if finite strain model is used, to determine if mechanical_strain or elastic_strain should be used
  bool _finite_strain_model;
  //const MaterialProperty<RankTwoTensor> & _mechanical_strain;
  const MaterialProperty<RankTwoTensor> & _strain;
  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
  /// Degradation Function
  const MaterialPropertyName _D_name;// ✅
  const MaterialProperty<Real> & _D;// ✅
  /// Accumulation energy type ✅
  const std::string _energy_calculation;
  ///II. Outputs from this program:
  /// 1. current energy
  const MaterialPropertyName _bar_psi_name;// ✅
  MaterialProperty<Real> & _bar_psi;
  const MaterialProperty<Real> & _bar_psi_old;
  /// 2. accumulate energy
  const MaterialPropertyName _acc_bar_psi_name;// ✅
  MaterialProperty<Real> & _acc_bar_psi;
  /// Old value of history variable
  const MaterialProperty<Real> & _acc_bar_psi_old;
  /// Other parameters
  const OptionalMaterialProperty<Real> & _n;
  const OptionalMaterialProperty<Real> & _R;  
  /// multiply by degradation or not
  const bool _multiply_by_degradation;
  /// accumulation mode
  const std::string _accumulation_mode;
  const VariableValue * _N_cyc_var;
  const VariableValue * _N_cyc_var_old;
//
};
