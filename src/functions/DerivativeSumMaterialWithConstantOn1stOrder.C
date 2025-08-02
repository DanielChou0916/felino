#include "DerivativeSumMaterialWithConstantOn1stOrder.h"

#include "libmesh/quadrature.h"

registerMooseObject("felinoApp", DerivativeSumMaterialWithConstantOn1stOrder);
registerMooseObject("felinoApp", ADDerivativeSumMaterialWithConstantOn1stOrder);

template <bool is_ad>
InputParameters
DerivativeSumMaterialWithConstantOn1stOrderTempl<is_ad>::validParams()
{ //mooseInfo("=== DEBUG: Entering validParams for DerivativeSumMaterialWithConstantOn1stOrderTempl ===");
  InputParameters params = DerivativeFunctionMaterialBaseTempl<is_ad>::validParams();
  params.addClassDescription("Meta-material to sum up multiple derivative materials");
  params.addParam<std::vector<std::string>>("sum_materials",
                                            "Base name of the parsed sum material property");
  // Additional constant variable
  params.addParam<std::vector<std::string>>("additional_sum_materials",
      "Base name of the additional material property summed only on the first order");
  // All arguments of the parsed expression (free energy) being summed
  params.addRequiredCoupledVar("args", "Vector of names of variables being summed");
  params.deprecateCoupledVar("args", "coupled_variables", "02/07/2024");

  params.addCoupledVar("displacement_gradients",
                       "Vector of displacement gradient variables (see "
                       "Modules/PhaseField/DisplacementGradients "
                       "action)");

  // Advanced arguments to construct a sum of the form \f$ c+\gamma\sum_iF_i \f$
  params.addParam<std::vector<Real>>("prefactor", {}, "Prefactor to multiply the sum term with.");
  params.addParam<Real>("constant", 0.0, "Constant to be added to the prefactor multiplied sum.");

  params.addParam<bool>("validate_coupling",
                        true,
                        "Check if all variables the specified materials depend on are listed in "
                        "the `coupled_variables` parameter.");
  params.addParamNamesToGroup("prefactor constant", "Advanced");
  params.addParam<bool>(
    "output_individual_derivatives",
    false,
    "Whether to output individual first derivatives for each summand material."
  );

  //mooseInfo("=== DEBUG: Done adding parameters in validParams ===");
  return params;
}

template <bool is_ad>
DerivativeSumMaterialWithConstantOn1stOrderTempl<is_ad>::DerivativeSumMaterialWithConstantOn1stOrderTempl(const InputParameters & parameters)
  : DerivativeFunctionMaterialBaseTempl<is_ad>(parameters),
    _sum_materials(this->template getParam<std::vector<std::string>>("sum_materials")),
    _num_materials(_sum_materials.size()),
    _prefactor(_num_materials, 1.0),
    _constant(this->template getParam<Real>("constant")),
    _validate_coupling(this->template getParam<bool>("validate_coupling")),
    // Additional constant variable
    // ✅
    _additional_summands_materials(this->template getParam<std::vector<std::string>>("additional_sum_materials")),
    _num_constant_materials(_additional_summands_materials.size()),
    _output_individual_derivatives(this->template getParam<bool>("output_individual_derivatives"))
{ // we need at least one constant material in the sum
  if (_num_constant_materials == 0)// ✅
    mooseError("Please supply at least one constant material, otherwise use DerivativeSumMaterial!", name());// ✅
  _add_summand_F.resize(_num_constant_materials);// ✅

  for (unsigned int n = 0; n < _num_constant_materials; ++n)// ✅
  {
    this->template getGenericMaterialProperty<Real, is_ad>(_additional_summands_materials[n]);
    _add_summand_F[n] = &this->template getGenericMaterialProperty<Real, is_ad>(_additional_summands_materials[n]);// ✅
  }
  // we need at least one material in the sum
  if (_num_materials == 0)
    mooseError("Please supply at least one material to sum in DerivativeSumMaterial ", name());

  // get prefactor values if not 1.0
  std::vector<Real> p = this->template getParam<std::vector<Real>>("prefactor");

  // if prefactor is used we need the same number of prefactors as sum materials
  if (_num_materials == p.size())
    _prefactor = p;
  else if (p.size() != 0)
    mooseError("Supply the same number of sum materials and prefactors.");

  // reserve space for summand material properties
  _summand_F.resize(_num_materials);
  _summand_dF.resize(_num_materials);
  _summand_d2F.resize(_num_materials);
  _summand_d3F.resize(_num_materials);

  for (unsigned int n = 0; n < _num_materials; ++n)
  {
    _summand_F[n] = &this->template getGenericMaterialProperty<Real, is_ad>(_sum_materials[n]);
    _summand_dF[n].resize(_nargs);
    _summand_d2F[n].resize(_nargs);
    _summand_d3F[n].resize(_nargs);

    for (unsigned int i = 0; i < _nargs; ++i)
    {
      _summand_dF[n][i] = &this->template getMaterialPropertyDerivative<Real, is_ad>(
          _sum_materials[n], _arg_names[i]);
      _summand_d2F[n][i].resize(_nargs);

      if (_third_derivatives)
        _summand_d3F[n][i].resize(_nargs);

      for (unsigned int j = 0; j < _nargs; ++j)
      {
        _summand_d2F[n][i][j] = &this->template getMaterialPropertyDerivative<Real, is_ad>(
            _sum_materials[n], _arg_names[i], _arg_names[j]);

        if (_third_derivatives)
        {
          _summand_d3F[n][i][j].resize(_nargs);

          for (unsigned int k = 0; k < _nargs; ++k)
            _summand_d3F[n][i][j][k] = &this->template getMaterialPropertyDerivative<Real, is_ad>(
                _sum_materials[n], _arg_names[i], _arg_names[j], _arg_names[k]);
        }
      }
    }
  }
  if (_output_individual_derivatives)
  {
    _summand_dF_out.reserve(_num_materials);
    for (unsigned int n = 0; n < _num_materials; ++n)
      _summand_dF_out.push_back(
          &this->template declareGenericProperty<Real, is_ad>(
              "first_derivative_" + _sum_materials[n]));
  }
}

template <bool is_ad>
void
DerivativeSumMaterialWithConstantOn1stOrderTempl<is_ad>::initialSetup()
{
  if (_validate_coupling)
    for (unsigned int n = 0; n < _num_materials; ++n)
      this->template validateCoupling<Real>(_sum_materials[n]);
}

template <bool is_ad>
void
DerivativeSumMaterialWithConstantOn1stOrderTempl<is_ad>::computeProperties()
{
  unsigned int i, j, k;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    // set function value
    if (_prop_F)
    {
      (*_prop_F)[_qp] = (*_summand_F[0])[_qp] * _prefactor[0];
      for (unsigned int n = 1; n < _num_materials; ++n)
        (*_prop_F)[_qp] += (*_summand_F[n])[_qp] * _prefactor[n];
    }

    for (i = 0; i < _nargs; ++i)
    {
      // set first derivatives
      if (_prop_dF[i])
      {
        (*_prop_dF[i])[_qp] = (*_summand_dF[0][i])[_qp] * _prefactor[0];
        for (unsigned int n = 1; n < _num_materials; ++n)
          (*_prop_dF[i])[_qp] += (*_summand_dF[n][i])[_qp] * _prefactor[n];
        // ✅ `Output Derivatives` 0801
        if (_output_individual_derivatives)
        {
          for (unsigned int n = 0; n < _num_materials; ++n)
            (*_summand_dF_out[n])[_qp] = (*_summand_dF[n][0])[_qp];
        }
        // [0] stands for only using the first derivative
        // ✅ `_additional_constant`
        for (unsigned int n = 0; n < _num_constant_materials; ++n) // ✅ 
        {
          (*_prop_dF[i])[_qp] += (*_add_summand_F[n])[_qp];

          // ✅ 
          //mooseInfo("qp=", _qp, " i=", i, " n=", n, " additional_sum_materials = ", (*_add_summand_F[n])[_qp]);
          //mooseInfo("This is the verification message that additional constant is applied successfully!", (*_add_summand_F[n])[_qp]);
          mooseInfo("This is the verification message that additional constant is applied successfully!");
        }
      }

      // second derivatives
      for (j = i; j < _nargs; ++j)
      {
        if (_prop_d2F[i][j])
        {
          (*_prop_d2F[i][j])[_qp] = (*_summand_d2F[0][i][j])[_qp] * _prefactor[0];
          for (unsigned int n = 1; n < _num_materials; ++n)
            (*_prop_d2F[i][j])[_qp] += (*_summand_d2F[n][i][j])[_qp] * _prefactor[n];
        }

        // third derivatives
        if (_third_derivatives)
        {
          for (k = j; k < _nargs; ++k)
            if (_prop_d3F[i][j][k])
            {
              (*_prop_d3F[i][j][k])[_qp] = (*_summand_d3F[0][i][j][k])[_qp] * _prefactor[0];
              for (unsigned int n = 1; n < _num_materials; ++n)
                (*_prop_d3F[i][j][k])[_qp] += (*_summand_d3F[n][i][j][k])[_qp] * _prefactor[n];
            }
        }
      }
    }
  }
}
