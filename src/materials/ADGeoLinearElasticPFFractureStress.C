//* This file is the customized phase field fracture model.

#include "ADGeoLinearElasticPFFractureStress.h"
#include "MathUtils.h"

registerMooseObject("SolidMechanicsApp", ADGeoLinearElasticPFFractureStress);

InputParameters
ADGeoLinearElasticPFFractureStress::validParams()
{
  InputParameters params = ADGeoPFFractureStressBase::validParams();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field "
                             "fracture model, with small strain");
  MooseEnum Decomposition("strain_spectral strain_vol_dev stress_spectral stress_dev proto1 proto2 proto3 none", "none");
  params.addParam<MooseEnum>("decomposition_type",
                             Decomposition,
                             "Decomposition approaches. Choices are: " +
                                 Decomposition.getRawNames());
  return params;
}

ADGeoLinearElasticPFFractureStress::ADGeoLinearElasticPFFractureStress(
    const InputParameters & parameters)
  : ADGeoPFFractureStressBase(parameters),
    GuaranteeConsumer(this),
    _decomposition_type(getParam<MooseEnum>("decomposition_type").getEnum<Decomposition_type>())
{
}

void
ADGeoLinearElasticPFFractureStress::initialSetup()
{
  if ((_decomposition_type == Decomposition_type::strain_vol_dev ||
       _decomposition_type == Decomposition_type::strain_spectral) &&
      !hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("Decomposition approach of strain_vol_dev and strain_spectral can only be used with "
               "isotropic elasticity tensor materials, use stress_spectral of stress_dev for anistropic "
               "elasticity tensor materials");
}


void
ADGeoLinearElasticPFFractureStress::computeStrainSpectral(ADReal & F_pos, ADReal & F_neg)
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
ADGeoLinearElasticPFFractureStress::computeStressSpectral(ADReal & F_pos, ADReal & F_neg)
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
ADGeoLinearElasticPFFractureStress::computeStrainVolDev(ADReal & F_pos, ADReal & F_neg)
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
ADGeoLinearElasticPFFractureStress::computeStressDev(ADReal & F_pos, ADReal & F_neg)
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
ADGeoLinearElasticPFFractureStress::computeStressSpDP(ADReal & F_pos, ADReal & F_neg)
{ mooseInfo("Now using prototype 1: Drucker-Prager factorial decomposition ");
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
  //1. Initialize variables : stick on negative part!!!!
  ADReal I1 = stress0neg.trace();
  ADReal J2 = stress0neg.secondInvariant();
  ADReal phi = _internal_friction[_qp] * 3.141592654 / 180.0;
  ADReal B = -1.0 / std::sqrt(3.0) * ((std::sin(phi)) / (1.0 - std::sin(phi)));
  ADReal A = _cohesion[_qp] / std::sqrt(3.0) * (2.0 * (std::cos(phi)) / (1.0 - std::sin(phi)));
  ADReal criterion, factor;
  if (J2<1e-8)
  { //mooseInfo("J2 value at qp ", _qp, " : ", J2);
    //mooseWarning("J2 is very small at qp ", _qp, ": ", J2);
    criterion = 0.0;
    factor = 0.0;
  }
  else 
  { 
    criterion = std::sqrt(J2)-A-B*I1;
    factor = std::max(1e-8, (1 - _R_res[_qp]) / (1 + exp(-30 * (criterion - _R_init[_qp]))));
    //mooseWarning("factor", _qp, ": ", factor);
  }
  //mooseInfo("factor value at qp ", _qp, " : ", factor);
  //2. Split stress0neg
  ADRankTwoTensor stress0neg_active, stress0neg_residual, stress0act;
  stress0neg_active = factor * stress0neg;
  stress0neg_residual = stress0neg - stress0neg_active;
  stress0act = stress0pos + stress0neg_active;

  if (_I_name.empty() || !(_I) || !(_pressure))
  {
      _stress[_qp] = stress0act * _D[_qp] + stress0neg_residual;
      _dstress_dc[_qp] = stress0act * _dDdc[_qp];
  }
  else   
  {
    _stress[_qp] = stress0act * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg_residual;
    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = stress0act * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }

  // Compute the positive and negative elastic energies
  F_pos = (stress0act).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  F_neg = (stress0neg_residual).doubleContraction(_mechanical_strain[_qp]) / 2.0;


  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0act * _dDdc[_qp];

}

void
ADGeoLinearElasticPFFractureStress::computeStressSpVol(ADReal & F_pos, ADReal & F_neg)
{ mooseInfo("Now using prototype 2: mixed decomposition ");
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
  ADRankTwoTensor stress0neg_dev= stress0neg.deviatoric();
  ADRankTwoTensor stress0neg_active, stress0neg_residual, stress0act;
  stress0neg_active = stress0neg_dev;
  stress0neg_residual = stress0neg - stress0neg_active;
  stress0act = stress0pos + stress0neg_active;

  if (_I_name.empty() || !(_I) || !(_pressure))
  {
      _stress[_qp] = stress0act * _D[_qp] + stress0neg_residual;
      _dstress_dc[_qp] = stress0act * _dDdc[_qp];
  }
  else   
  {
    _stress[_qp] = stress0act * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg_residual;
    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = stress0act * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }

  // Compute the positive and negative elastic energies
  F_pos = (stress0act).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  F_neg = (stress0neg_residual).doubleContraction(_mechanical_strain[_qp]) / 2.0;


  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0act * _dDdc[_qp];

}

void
ADGeoLinearElasticPFFractureStress::computeStressSpVolR(ADReal & F_pos, ADReal & F_neg)
{ mooseInfo("Now using prototype 3: mixed decomposition with ratio");
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
  ADRankTwoTensor stress0neg_dev= stress0neg.deviatoric();
  ADRankTwoTensor stress0neg_active, stress0neg_residual, stress0act;
  stress0neg_active = stress0neg_dev;
  //stress0neg_residual = stress0neg - stress0neg_active;
  stress0act = _D[_qp]*stress0pos + (1-_D[_qp])*stress0neg_active;// ✅
  stress0neg_residual = stress - stress0act;// ✅
  if (_I_name.empty() || !(_I) || !(_pressure))
  {
      _stress[_qp] = stress0act * _D[_qp] + stress0neg_residual;
      _dstress_dc[_qp] = stress0act * _dDdc[_qp];
  }
  else   
  {
    _stress[_qp] = stress0act * _D[_qp] - _pressure[_qp] * I2 * _I[_qp] + stress0neg_residual;
    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = stress0act * _dDdc[_qp] - _pressure[_qp] * I2 * _dIdc[_qp];
  }

  // Compute the positive and negative elastic energies
  F_pos = (stress0act).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  F_neg = (stress0neg_residual).doubleContraction(_mechanical_strain[_qp]) / 2.0;


  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = stress0act * _dDdc[_qp];

}

void
ADGeoLinearElasticPFFractureStress::computeQpStress()
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
    case Decomposition_type::proto1:
      computeStressSpDP(F_pos, F_neg);
      break;
    case Decomposition_type::proto2:
      computeStressSpVol(F_pos, F_neg);
      break;
    case Decomposition_type::proto3:
      computeStressSpVolR(F_pos, F_neg);
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
