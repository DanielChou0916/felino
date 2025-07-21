//

#pragma once

#include "Material.h"
#include "GuaranteeConsumer.h"
class ADComputeExtraDrivingForce : public Material,
                                   public GuaranteeConsumer
{
public:
  static InputParameters validParams();

  ADComputeExtraDrivingForce(const InputParameters & parameters);
  void initialSetup() override;
protected:
  virtual void computeQpProperties() override;
  const std::string _uncracked_base_name;
  bool _finite_strain_model;
  const ADMaterialProperty<RankTwoTensor> & _strain;
  /// Name of the elasticity tensor material property
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;
  /// Degradation Function
  const MaterialPropertyName _D_name;// ✅0318
  const ADMaterialProperty<Real> & _D;// ✅0318
  /// Material property defining crack width, declared elsewhere
  const ADMaterialProperty<Real> & _l;
  /// Energy release rate
  const ADMaterialProperty<Real> & _gc;
  /// Tensile and compressive strength
  const ADMaterialProperty<Real> & _sigma_ts;
  const ADMaterialProperty<Real> & _sigma_cs;
  /// Unitless coefficient
  const ADMaterialProperty<Real> & _delta_elp;
  /// Model type (finite_e or asymptotic) ✅
  const std::string _model;
  ///II. Outputs from this program:
  /// Material property for extra driving force
  const MaterialPropertyName _Ce_name;// ✅
  ADMaterialProperty<Real> & _Ce;
  //const bool _output_Ce_aux;
//
};
