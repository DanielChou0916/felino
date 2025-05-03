//* This file is the customized AD phase field fracture model.

#pragma once

#include "ADGeoPFFractureStressBase.h"
#include "MooseEnum.h"
#include "GuaranteeConsumer.h"

/**
 * Phase-field fracture
 * This class computes the stress and energy contribution for the
 * small strain Linear Elastic formulation of phase field fracture
 */
class ADGeoLinearElasticPFFractureStress : public ADGeoPFFractureStressBase,
                                             public GuaranteeConsumer
{
public:
  static InputParameters validParams();

  ADGeoLinearElasticPFFractureStress(const InputParameters & parameters);

  void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /**
   * Method to split elastic energy based on strain spectral decomposition
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStrainSpectral(ADReal & F_pos, ADReal & F_neg);

  /**
   * Method to split elastic energy based on strain volumetric/deviatoric decomposition
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStrainVolDev(ADReal & F_pos, ADReal & F_neg);

  /**
   * Method to split elastic energy based on stress spectral decomposition
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStressSpectral(ADReal & F_pos, ADReal & F_neg);

  // New method for decomposition

  /**
   * Method to split elastic energy based on Drucker Prager model
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStressDev(ADReal & F_pos, ADReal & F_neg);

  /**
   * Method to split elastic energy based on Drucker Prager model
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStressSpDP(ADReal & F_pos, ADReal & F_neg);   

  /**
   * Method to split elastic energy based on Drucker Prager model
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStressSpVol(ADReal & F_pos, ADReal & F_neg);  
     
  /**
   * Method to split elastic energy based on Drucker Prager model
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStressSpVolR(ADReal & F_pos, ADReal & F_neg);  

  /// Decomposittion type
  enum class Decomposition_type
  {
    strain_spectral,
    strain_vol_dev,
    stress_spectral,
    stress_dev,
    proto1,
    proto2,
    proto3,
    none
  } _decomposition_type;
};
