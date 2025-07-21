//

#pragma once

#include "Material.h"
#include "GuaranteeConsumer.h"
class ADComputeExtraDrivingForcefromStress : public Material,
                                   public GuaranteeConsumer
{
public:
  static InputParameters validParams();

  ADComputeExtraDrivingForcefromStress(const InputParameters & parameters);
  void initialSetup() override;
protected:
  virtual void computeQpProperties() override;
    /// Base name of the input elasticity tensor 
  const std::string _uncracked_base_name;
  const std::string _elasticity_tensor_name;
  /// Elasticity tensor material property
  const ADMaterialProperty<RankFourTensor> & _elasticity_tensor;
  /// Indicator if finite strain model is used, to determine if mechanical_strain or elastic_strain should be used
  //bool _finite_strain_model;
  /// Mechanical_strain if finite_strain_model = false, otherwise elastic_strain
  //const ADMaterialProperty<RankTwoTensor> & _strain;
  /// Uncracked stress calculated by another material
  const ADMaterialProperty<RankTwoTensor> & _stress;
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
//
};
