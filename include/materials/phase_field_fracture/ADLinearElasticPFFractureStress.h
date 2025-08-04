//* This file is the customized AD phase field fracture model.

#pragma once

#include "ADPFFractureStressBase.h"
#include "MooseEnum.h"
#include "GuaranteeConsumer.h"

/**
 * Phase-field fracture
 * This class computes the stress and energy contribution for the
 * small strain Linear Elastic formulation of phase field fracture
 */
class ADLinearElasticPFFractureStress : public ADPFFractureStressBase,
                                             public GuaranteeConsumer
{
public:
  static InputParameters validParams();

  ADLinearElasticPFFractureStress(const InputParameters & parameters);

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
   * Method to split elastic energy based on stress dev decomposition
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStressDev(ADReal & F_pos, ADReal & F_neg);

  /**
   * Method to split elastic energy based on RCE model
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStressRCE(ADReal & F_pos, ADReal & F_neg);
  /**
   * Method to split elastic energy based on Drucker Prager model
   * @param F_pos tensile part of total elastic energy
   * @param F_neg compressive part of total elastic energy
   */
  void computeStrainDruckerPrager(ADReal & F_pos, ADReal & F_neg);  
  /// Decomposittion type
  enum class Decomposition_type
  {
    strain_spectral,
    strain_vol_dev,
    stress_spectral,
    stress_dev,
    rce,
    strain_dp,
    none
  } _decomposition_type;
};
