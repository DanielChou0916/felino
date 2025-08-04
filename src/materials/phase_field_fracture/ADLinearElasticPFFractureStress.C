//* This file is the customized phase field fracture model.

#include "ADLinearElasticPFFractureStress.h"
#include "MathUtils.h"

registerMooseObject("SolidMechanicsApp", ADLinearElasticPFFractureStress);

InputParameters
ADLinearElasticPFFractureStress::validParams()
{
  InputParameters params = ADPFFractureStressBase::validParams();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field "
                             "fracture model, with small strain");
  MooseEnum Decomposition("strain_spectral strain_vol_dev stress_spectral stress_dev rce strain_dp none", "none");
  params.addParam<MooseEnum>("decomposition_type",
                             Decomposition,
                             "Decomposition approaches. Choices are: " +
                                 Decomposition.getRawNames());
  return params;
}

ADLinearElasticPFFractureStress::ADLinearElasticPFFractureStress(
    const InputParameters & parameters)
  : ADPFFractureStressBase(parameters),
    GuaranteeConsumer(this),
    _decomposition_type(getParam<MooseEnum>("decomposition_type").getEnum<Decomposition_type>())
{
}

void
ADLinearElasticPFFractureStress::initialSetup()
{
  if ((_decomposition_type == Decomposition_type::strain_vol_dev  ||
       _decomposition_type == Decomposition_type::strain_spectral ||
       _decomposition_type == Decomposition_type::rce ||
       _decomposition_type == Decomposition_type::strain_dp) &&
      !hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("Decomposition approach of strain_vol_dev and strain_spectral can only be used with "
               "isotropic elasticity tensor materials, use stress_spectral of stress_dev for anistropic "
               "elasticity tensor materials");
}


void
ADLinearElasticPFFractureStress::computeStrainSpectral(ADReal & F_pos, ADReal & F_neg)
{
  // Isotropic elasticity is assumed and should be enforced
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);

  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  // Compute eigenvectors and eigenvalues of mechanical strain and projection tensor
  ADRankTwoTensor eigvec;
  std::vector<ADReal> eigval(LIBMESH_DIM);
  ADRankFourTensor Ppos =
      _mechanical_strain[_qp].positiveProjectionEigenDecomposition(eigval, eigvec);
  ADRankFourTensor I4sym(ADRankFourTensor::initIdentitySymmetricFour);

  // Calculate tensors of outerproduct of eigen vectors
  std::vector<ADRankTwoTensor> etens(LIBMESH_DIM);

  for (const auto i : make_range(Moose::dim))
    etens[i] = ADRankTwoTensor::selfOuterProduct(eigvec.column(i));

  // Separate out positive and negative eigen values
  std::vector<ADReal> epos(LIBMESH_DIM), eneg(LIBMESH_DIM);
  for (const auto i : make_range(Moose::dim))
  {
    epos[i] = (std::abs(eigval[i]) + eigval[i]) / 2.0;
    eneg[i] = -(std::abs(eigval[i]) - eigval[i]) / 2.0;
  }

  // Seprate positive and negative sums of all eigenvalues
  ADReal etr = 0.0;
  for (const auto i : make_range(Moose::dim))
    etr += eigval[i];

  const ADReal etrpos = (std::abs(etr) + etr) / 2.0;
  const ADReal etrneg = -(std::abs(etr) - etr) / 2.0;

  // Calculate the tensile (postive) and compressive (negative) parts of stress
  ADRankTwoTensor stress0pos, stress0neg;
  for (const auto i : make_range(Moose::dim))
  {
    stress0pos += etens[i] * (lambda * etrpos + 2.0 * mu * epos[i]);
    stress0neg += etens[i] * (lambda * etrneg + 2.0 * mu * eneg[i]);
  }

  // sum squares of epos and eneg
  ADReal pval(0.0), nval(0.0);
  for (const auto i : make_range(Moose::dim))
  {
    pval += epos[i] * epos[i];
    nval += eneg[i] * eneg[i];
  }

  if (_I_name.empty() || !(_I) || !(_pressure))
  {
      _stress[_qp] = stress0pos * _D[_qp] + stress0neg;
      _dstress_dc[_qp] = stress0pos * _dDdc[_qp];
  }
  else   
  {
      _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;
      _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }
  // Energy with positive principal strains
  F_pos = lambda * etrpos * etrpos / 2.0 + mu * pval;
  F_neg = lambda * etrneg * etrneg / 2.0 + mu * nval;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];



  //_Jacobian_mult[_qp] = (I4sym - (1 - _D[_qp]) * Ppos) * _elasticity_tensor[_qp];
}

void
ADLinearElasticPFFractureStress::computeStressSpectral(ADReal & F_pos, ADReal & F_neg)
{
  // Compute Uncracked stress
  ADRankTwoTensor stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp];

  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  // Create the positive and negative projection tensors
  ADRankFourTensor I4sym(ADRankFourTensor::initIdentitySymmetricFour);
  std::vector<ADReal> eigval;
  ADRankTwoTensor eigvec;
  ADRankFourTensor Ppos = stress.positiveProjectionEigenDecomposition(eigval, eigvec);

  // Project the positive and negative stresses
  ADRankTwoTensor stress0pos = Ppos * stress;
  ADRankTwoTensor stress0neg = stress - stress0pos;



  if (_I_name.empty() || !(_I) || !(_pressure))
  {
      _stress[_qp] = stress0pos * _D[_qp] + stress0neg;
      _dstress_dc[_qp] = stress0pos * _dDdc[_qp];
  }
  else   
  {
    _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;
    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }

  // Compute the positive and negative elastic energies
  F_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  F_neg = (stress0neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;


  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  //_Jacobian_mult[_qp] = (I4sym - (1 - _D[_qp]) * Ppos) * _elasticity_tensor[_qp];
}

void
ADLinearElasticPFFractureStress::computeStrainVolDev(ADReal & F_pos, ADReal & F_neg)
{
  // Isotropic elasticity is assumed and should be enforced
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const ADReal k = lambda + 2.0 * mu / LIBMESH_DIM;

  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankFourTensor I2I2 = I2.outerProduct(I2);

  //ADRankFourTensor Jacobian_pos, Jacobian_neg;
  ADRankTwoTensor strain0vol, strain0dev;
  ADRankTwoTensor stress0pos, stress0neg;
  ADReal strain0tr, strain0tr_neg, strain0tr_pos;

  strain0dev = _mechanical_strain[_qp].deviatoric();
  strain0vol = _mechanical_strain[_qp] - strain0dev;
  strain0tr = _mechanical_strain[_qp].trace();
  strain0tr_neg = std::min(strain0tr, 0.0);
  strain0tr_pos = strain0tr - strain0tr_neg;
  stress0neg = k * strain0tr_neg * I2;
  stress0pos = _elasticity_tensor[_qp] * _mechanical_strain[_qp] - stress0neg;
  // Energy with positive principal strains
  ADRankTwoTensor strain0dev2 = strain0dev * strain0dev;

  if (_I_name.empty() || !(_I) || !(_pressure))
  {
    _stress[_qp] = stress0pos * _D[_qp] + stress0neg;
    _dstress_dc[_qp] = stress0pos * _dDdc[_qp];
  }
  else   
  {
    _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;
    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }

  F_pos = 0.5 * k * strain0tr_pos * strain0tr_pos + mu * strain0dev2.trace();
  F_neg = 0.5 * k * strain0tr_neg * strain0tr_neg;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  //if (strain0tr < 0)
  //  Jacobian_neg = k * I2I2;
  //Jacobian_pos = _elasticity_tensor[_qp] - Jacobian_neg;
  //_Jacobian_mult[_qp] = _D[_qp] * Jacobian_pos + Jacobian_neg;
}


void
ADLinearElasticPFFractureStress::computeStressDev(ADReal & F_pos, ADReal & F_neg)
{
  //const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  //const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  //const ADReal k = lambda + 2.0 * mu / LIBMESH_DIM;
  // Compute Uncracked stress
  ADRankTwoTensor stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
  
  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankFourTensor I2I2 = I2.outerProduct(I2);
  
  ADRankTwoTensor stress0dev;
  ADRankFourTensor Jacobian_pos, Jacobian_neg;
  stress0dev = stress.deviatoric();
  ADReal stress0hydro = stress.trace()/LIBMESH_DIM;//should be /LIBMESH_DIM or /3??
  ADReal stress0hydro_neg = std::min(stress0hydro, 0.0);
  ADReal stress0hydro_pos = stress0hydro - stress0hydro_neg;
  // Create the positive and negative stress parts
  ADRankTwoTensor stress0pos = stress0hydro_pos * I2 + stress0dev;
  ADRankTwoTensor stress0neg = stress - stress0pos;



  if (_I_name.empty() || !(_I) || !(_pressure))
  {
      _stress[_qp] = stress0pos * _D[_qp] + stress0neg;
      _dstress_dc[_qp] = stress0pos * _dDdc[_qp];
  }
  else   
  {
    _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;
    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }

  // Compute the positive and negative elastic energies
  F_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  F_neg = (stress0neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;


  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  //if (stress0hydro < 0)
  //  Jacobian_neg = (1.0 / LIBMESH_DIM) * k * I2I2;
  //Jacobian_pos = _elasticity_tensor[_qp] - Jacobian_neg;
  //_Jacobian_mult[_qp] = _D[_qp] * Jacobian_pos + Jacobian_neg;
}

void
ADLinearElasticPFFractureStress::computeStressRCE(ADReal & F_pos, ADReal & F_neg)
{
  // Compute Uncracked stress
  ADRankTwoTensor stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp];

  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  // Create the positive and negative projection tensors
  //ADRankFourTensor I4sym(ADRankFourTensor::initIdentitySymmetricFour);
  std::vector<ADReal> eigval;
  ADRankTwoTensor eigvec;
  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);

  //Define projections
  std::vector<ADRankTwoTensor> projections(LIBMESH_DIM);
  for (const auto i : make_range(Moose::dim))
    projections[i] =0.5*
                    (ADRankTwoTensor::outerProduct(eigvec.column(Moose::dim-1),eigvec.column(i))
                    +ADRankTwoTensor::outerProduct(eigvec.column(i),eigvec.column(Moose::dim-1)));
  
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  std::vector<ADReal> Lam(LIBMESH_DIM);
  for (const auto i : make_range(Moose::dim))
    if (i==Moose::dim-1)
    { 
      ADReal L1 = lambda/(lambda+2*mu)* _mechanical_strain[_qp].trace()
                  + 2*mu/(lambda+2*mu)* (_mechanical_strain[_qp]).doubleContraction(projections[i]);
      Lam[i] = std::max(L1,0.0);
      //if Lam[2]<0 
        //Lam[2] = 0.0;  
    }
    else 
    {
      Lam[i] = 2* (_mechanical_strain[_qp]).doubleContraction(projections[i]);
      //if Lam[2]<0 
        //Lam[0 and 1] will account some shear resistance ;  
    }

  // compute jump strain
  ADRankTwoTensor strain_jump;
  for (const auto i : make_range(Moose::dim))
    strain_jump += Lam[i]*projections[i];

  ADRankTwoTensor delta_strain = _mechanical_strain[_qp] - strain_jump; 
  
  // Project the positive and negative stresses
  ADRankTwoTensor stress0neg = _elasticity_tensor[_qp]*delta_strain;
  ADRankTwoTensor stress0pos = stress-stress0neg;

  ADReal F0 = stress.doubleContraction(_mechanical_strain[_qp]) / 2.0;
  F_neg = stress0neg.doubleContraction(delta_strain) / 2.0;
  F_pos = F0 - F_neg;  

  if (_I_name.empty() || !(_I) || !(_pressure))
  {
      _stress[_qp] = stress0pos * _D[_qp] + stress0neg;
      _dstress_dc[_qp] = stress0pos * _dDdc[_qp];
  }
  else   
  {
    _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;
    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }



  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  //_Jacobian_mult[_qp] = (I4sym - (1 - _D[_qp]) * Ppos) * _elasticity_tensor[_qp];
}

void
ADLinearElasticPFFractureStress::computeStrainDruckerPrager(ADReal & F_pos, ADReal & F_neg)
{
  // Isotropic elasticity is assumed and should be enforced
  const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  const ADReal k = lambda + 2.0 * mu / LIBMESH_DIM;

  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  ADRankFourTensor I2I2 = I2.outerProduct(I2);
  //RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  //RankFourTensor Jacobian_pos;
  //RankFourTensor Jacobian_neg(RankFourTensor::initIdentitySymmetricFour);
  ADRankTwoTensor stress0pos, stress0neg;
  
  ADRankTwoTensor strain0dev = _mechanical_strain[_qp].deviatoric();
  ADReal strain0tr = _mechanical_strain[_qp].trace();
  ADReal J2_strain  = 0.5 * strain0dev.doubleContraction(strain0dev);
  
  ADReal J2_sqrt = std::sqrt(std::max(J2_strain, 0.0));
  ADReal J2_sqrt_safe = (J2_strain < 1e-14) ? 1e-7 : J2_sqrt;
  // Stress for intact material
  ADRankTwoTensor stress0 = k * strain0tr * I2 + 2 * mu * strain0dev;
  // Strain energy density
  ADReal F0 = 0.5 * k * strain0tr * strain0tr + mu * strain0dev.doubleContraction(strain0dev);

  // Essential terms for if-else conditions
  //ADReal term1 = -6 * _B[_qp] * J2_sqrt;
  //ADReal term2 = 2 * mu * J2_sqrt;
  //ADReal term3 = 3 * _B[_qp] * k * strain0tr;

  if (-6 * _B[_qp] * J2_sqrt < strain0tr) 
  {
      // Case 1: Strain is primarily volumetric
      F_neg = 0.0;
      stress0neg = 1e-10 * I2;  // Negligible stress in this case

      // Jacobian_neg for Case 1
      //Jacobian_neg = 0.0*Jacobian_neg;  // Zero tensor
  } 
  else if (2 * mu * J2_sqrt < 3 * _B[_qp] * k * strain0tr)
  {
      // Case 3: Strain is primarily deviatoric
      F_neg = 0.5 * k * strain0tr * strain0tr + mu * strain0dev.doubleContraction(strain0dev);

      stress0neg = k * strain0tr * I2 + 2 * mu * strain0dev;

      // Jacobian_neg for Case 3
      //Jacobian_neg = k * I2I2 + 2 * mu * I4sym;
  }
  else
  {
      // Case 2: Combined volumetric and deviatoric strain
      F_neg = (k * mu) / (18 * _B[_qp] * _B[_qp] * k + 2 * mu) *
              (strain0tr + 6 * _B[_qp] * J2_sqrt_safe) *
              (strain0tr + 6 * _B[_qp] * J2_sqrt_safe);

      stress0neg = (k * mu) / (18 * _B[_qp] * _B[_qp] * k + 2 * mu) *
                  (2 * (strain0tr + 6 * _B[_qp] * J2_sqrt_safe) * I2 +
                    6 * _B[_qp] / J2_sqrt_safe *
                    (strain0tr + 6 * _B[_qp] * J2_sqrt_safe) * strain0dev);

      // Jacobian_neg for Case 2
      //Real c1 = (-18 * k * k * _B[_qp] * _B[_qp]) / (18 * _B[_qp] * _B[_qp] * k + 2 * mu);
      //Real c2 = (2 * k * mu * _B[_qp]) / (18 * _B[_qp] * _B[_qp] * k + 2 * mu);
      //Real c3 = (4 * mu * mu) / (18 * _B[_qp] * _B[_qp] * k + 2 * mu);
      //Real c4 = (2 * mu) / (18 * _B[_qp] * _B[_qp] * k + 2 * mu);

      //Jacobian_neg = c1 * I2.outerProduct(I2) +
      //              c2 * I2.outerProduct(strain0dev) +
      //              c3 * strain0dev.outerProduct(strain0dev) +
      //              c4 * I4sym;
  } 

  stress0pos = stress0 - stress0neg;
  F_pos=F0-F_neg;

  if (_I_name.empty() || !(_I) || !(_pressure))
  {
      _stress[_qp] = stress0pos * _D[_qp] + stress0neg;
      _dstress_dc[_qp] = stress0pos * _dDdc[_qp];
  }
  else   
  {
    _stress[_qp] = stress0pos * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg;
    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = stress0pos * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }
  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0pos * _dDdc[_qp];

  //Jacobian_pos = _elasticity_tensor[_qp] - Jacobian_neg;
  //_Jacobian_mult[_qp] = _D[_qp] * Jacobian_pos + Jacobian_neg;
}


void
ADLinearElasticPFFractureStress::computeQpStress()
{
  ADReal F_pos, F_neg;
  ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  //Start:Check is optional variables exist here
  
  //End:Check is optional variables exist here
  switch (_decomposition_type)
  {
    case Decomposition_type::strain_spectral:
      computeStrainSpectral(F_pos, F_neg);
      break;
    case Decomposition_type::strain_vol_dev:
      computeStrainVolDev(F_pos, F_neg);
      break;
    case Decomposition_type::stress_spectral:
      computeStressSpectral(F_pos, F_neg);
      break;
    case Decomposition_type::stress_dev:
      computeStressDev(F_pos, F_neg);
      break;
    case Decomposition_type::rce:
      computeStressRCE(F_pos, F_neg);
      break;
    case Decomposition_type::strain_dp:
      computeStrainDruckerPrager(F_pos, F_neg);
      if (!(_B))
        mooseError("Drucker-Prager material parameter 'B' is missing! You must set 'B' in the input file.");
      break;
    default:
    {
      ADRankTwoTensor stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
      F_pos = stress.doubleContraction(_mechanical_strain[_qp]) / 2.0;
      F_neg = 0.0;
      if (_use_current_hist)
        _d2Fdcdstrain[_qp] = stress * _dDdc[_qp];

      if (_I_name.empty() || !(_I) || !(_pressure))
      {
        _stress[_qp] = _D[_qp] * stress;
        _dstress_dc[_qp] = stress * _dDdc[_qp];
      }
      else   
      {
        _stress[_qp] = _D[_qp] * stress - _pressure[_qp] * I2 * _I[_qp];
        _dstress_dc[_qp] = stress * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
      }

      //_Jacobian_mult[_qp] = _D[_qp] * _elasticity_tensor[_qp];
    }
  }

  // // Assign history variable
  ADReal hist_variable = _H_old[_qp];
  if (_use_snes_vi_solver)
  {
    _H[_qp] = F_pos;

    if (_use_current_hist)
      hist_variable = _H[_qp];
  }
  else
  {
    if (F_pos > _H_old[_qp])
      _H[_qp] = F_pos;
    else
      _H[_qp] = _H_old[_qp];

    if (_use_current_hist)
      hist_variable = _H[_qp];

    if (_barrier)
    {
      if (hist_variable < _barrier[_qp])
        hist_variable = _barrier[_qp];
    }
  }

  // Elastic free energy density

  if (_I_name.empty() || !(_I) || !(_pressure))
  {
    _E[_qp] = hist_variable * _D[_qp] + F_neg;
    _dEdc[_qp] = hist_variable * _dDdc[_qp];
    _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp];
  }
  else   
  {
    _E[_qp] =
        hist_variable * _D[_qp] + F_neg - _pressure[_qp] * _mechanical_strain[_qp].trace() * _I[_qp];
    _dEdc[_qp] =
        hist_variable * _dDdc[_qp] - _pressure[_qp] * _mechanical_strain[_qp].trace() * _dIdc[_qp];
    _d2Ed2c[_qp] = hist_variable * _d2Dd2c[_qp] -
                  _pressure[_qp] * _mechanical_strain[_qp].trace() * _d2Id2c[_qp];
  }

}
