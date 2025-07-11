#include "ACInterfaceGradKappa.h"

registerMooseObject("PhaseFieldApp", ACInterfaceGradKappa);

InputParameters
ACInterfaceGradKappa::validParams()
{
  InputParameters params = JvarMapKernelInterface<Kernel>::validParams();
  params.addClassDescription("Gradient energy Allen-Cahn Kernel");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "The kappa used with the kernel");
  params.addParam<bool>("variable_L",
                        true,
                        "The mobility is a function of any MOOSE variable (if "
                        "this is set to false L must be constant over the "
                        "entire domain!)");
  params.addParam<bool>("use_grad_kappa",
                        false,
                        "Set to false if L is constant over the domain.");
  params.addCoupledVar("grad_kappa_x", "AuxVariable for ∂κ/∂x");
  params.addCoupledVar("grad_kappa_y", "AuxVariable for ∂κ/∂y");
  params.addCoupledVar("grad_kappa_z", "AuxVariable for ∂κ/∂z (only in 3D)");  
  params.addParam<bool>("use_anisotropic_matrix",
                        false,
                        "Set to false for isotropic PFF");  
  params.addParam<MaterialPropertyName>("anisotropic_matrix_name", "A", "The anisotropic matrix used with the kernel");                    
  return params;
}

ACInterfaceGradKappa::ACInterfaceGradKappa(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    _L(getMaterialProperty<Real>("mob_name")),
    _kappa(getMaterialProperty<Real>("kappa_name")),
    _variable_L(getParam<bool>("variable_L")),
    _dLdop(getMaterialPropertyDerivative<Real>("mob_name", _var.name())),
    _d2Ldop2(getMaterialPropertyDerivative<Real>("mob_name", _var.name(), _var.name())),
    _dkappadop(getMaterialPropertyDerivative<Real>("kappa_name", _var.name())),
    _dLdarg(_n_args),
    _d2Ldargdop(_n_args),
    _d2Ldarg2(_n_args),
    _dkappadarg(_n_args),
    _gradarg(_n_args),
    _use_grad_kappa(getParam<bool>("use_grad_kappa")),
    _grad_kappa_x((_use_grad_kappa && parameters.isParamValid("grad_kappa_x"))
                  ? &coupledValue("grad_kappa_x"): nullptr),
    _grad_kappa_y((_use_grad_kappa && parameters.isParamValid("grad_kappa_y"))
                  ? &coupledValue("grad_kappa_y"): nullptr),
    _grad_kappa_z((_use_grad_kappa && parameters.isParamValid("grad_kappa_z"))
                  ? &coupledValue("grad_kappa_z"): nullptr),
    _use_anisotropic_matrix(getParam<bool>("use_anisotropic_matrix")),// ✅ 2025/07/06
    _A_ptr((nullptr))// ✅ 2025/07/06
{
  // Get mobility and kappa derivatives and coupled variable gradients
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    MooseVariable * ivar = _coupled_standard_moose_vars[i];
    const VariableName iname = ivar->name();
    if (iname == _var.name())
    {
      if (isCoupled("args"))
        paramError("args",
                   "The kernel variable should not be specified in the coupled `args` parameter.");
      else
        paramError("coupled_variables",
                   "The kernel variable should not be specified in the coupled `coupled_variables` "
                   "parameter.");
    }

    _dLdarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", i);
    _dkappadarg[i] = &getMaterialPropertyDerivative<Real>("kappa_name", i);
    _d2Ldargdop[i] = &getMaterialPropertyDerivative<Real>("mob_name", iname, _var.name());

    _gradarg[i] = &(ivar->gradSln());

    _d2Ldarg2[i].resize(_n_args);
    for (unsigned int j = 0; j < _n_args; ++j)
      _d2Ldarg2[i][j] = &getMaterialPropertyDerivative<Real>("mob_name", i, j);
  }
// constructor
  if (_use_grad_kappa && (!parameters.isParamValid("grad_kappa_x") || !parameters.isParamValid("grad_kappa_y")))
  {
    if (processor_id() == 0)
      mooseWarning("use_grad_kappa is true, but grad_kappa_x or grad_kappa_y not provided → skipped");
  }
//anisotropic matrix A
  if (_use_anisotropic_matrix)// ✅ 2025/07/06
  { 
    const std::string & name = getParam<MaterialPropertyName>("anisotropic_matrix_name");
    if (hasMaterialProperty<RankTwoTensor>(name))
      {
        mooseInfo("Anisotropic matrix loaded and applied in residual");
        _A_ptr = &getMaterialProperty<RankTwoTensor>(name);
      }
    else if (processor_id() == 0)
      mooseWarning("Anisotropic matrix material '", name, "' not found. Falling back to identity.");
  }
  else
    mooseInfo("Use isotropic PFF!");
}

void
ACInterfaceGradKappa::initialSetup()
{
  validateCoupling<Real>("mob_name");
  validateCoupling<Real>("kappa_name");
}

RealGradient
ACInterfaceGradKappa::gradL()
{
  RealGradient g = _grad_u[_qp] * _dLdop[_qp];
  for (unsigned int i = 0; i < _n_args; ++i)
    g += (*_gradarg[i])[_qp] * (*_dLdarg[i])[_qp];
  return g;
}

RealGradient
ACInterfaceGradKappa::nablaLPsi()
{
  // sum is the product rule gradient \f$ \nabla (L\psi) \f$
  RealGradient sum = _L[_qp] * _grad_test[_i][_qp];

  if (_variable_L)
    sum += gradL() * _test[_i][_qp];

  return sum;
}

RealGradient
ACInterfaceGradKappa::kappaNablaLPsi()
{
  return _kappa[_qp] * nablaLPsi();
}

RealGradient
ACInterfaceGradKappa::gradKappa()
{
  RealGradient g; // default (0,0,0)

  if (_use_grad_kappa)
  {
    mooseInfo("use_grad_kappa = true");
    if (_grad_kappa_x && _grad_kappa_y)
      {
      g(0) = (*_grad_kappa_x)[_qp];
      g(1) = (*_grad_kappa_y)[_qp];

      if (LIBMESH_DIM == 3 && _grad_kappa_z)
        g(2) = (*_grad_kappa_z)[_qp];
      mooseInfo("grad_kappa loaded and applied in residual");
      }
  }
  return g;
}

Real
ACInterfaceGradKappa::computeQpResidual()
{
  RankTwoTensor A = (_use_anisotropic_matrix && _A_ptr) ? (*_A_ptr)[_qp] : RankTwoTensor::initIdentity;// ✅ 2025/07/06
  Real r = (A* _grad_u[_qp]) * kappaNablaLPsi();
  RealGradient gk = gradKappa();
  // -(∇κ · ∇u) * L * test
  r -= (gk * (A* _grad_u[_qp])) * _L[_qp] * _test[_i][_qp];
  return r;
}

Real
ACInterfaceGradKappa::computeQpJacobian()
{
  // --------------------------------------------------------------------------
  // 1) Original Jacobian from the κ ∇(L ψ) · ∇ test term
  // dsum = ∂/∂η [ κ ∇(L ψ) ] · ∇ test
  // --------------------------------------------------------------------------
  RealGradient dsum =
      (_dkappadop[_qp] * _L[_qp] + _kappa[_qp] * _dLdop[_qp])   // dκ/dη * L + κ * dL/dη
      * _phi[_j][_qp]                                           // multiply by φ_j
      * _grad_test[_i][_qp];                                    // and ∇ test

  // if L depends on the solution, add its second-derivative contributions
  if (_variable_L)
  {
    // derivative of ∇L = ∇φ_j * dL/dη + ∇u * φ_j * d²L/dη²
    RealGradient dgradL =
        _grad_phi[_j][_qp] * _dLdop[_qp]                       // ∇φ_j · dL/dη
      + _grad_u[_qp]     * _phi[_j][_qp] * _d2Ldop2[_qp];      // ∇u · φ_j · d²L/dη²

    // include cross-derivatives for any other coupled variables
    for (unsigned int i = 0; i < _n_args; ++i)
      dgradL += (*_gradarg[i])[_qp]                            // ∇ arg_i
              * _phi[_j][_qp]                                  // times φ_j
              * (*_d2Ldargdop[i])[_qp];                        // times ∂²L/(∂arg_i ∂η)

    // add to dsum: κ * d(∇L)/dη + dκ/dη * ∇L
    dsum += (_kappa[_qp] * dgradL
          + _dkappadop[_qp] * _phi[_j][_qp] * gradL())         // gradL() computes ∇(L ψ)
          * _test[_i][_qp];                                    // times test
  }

  // --------------------------------------------------------------------------
  // 2) Off-diagonal block for same κ∇Lψ term (no change needed here if 
  //    only added grad κ to the residual; only the diagonal η–η block gets it)
  // --------------------------------------------------------------------------
  // (already have computeQpOffDiagJacobian for other variables.)

  // --------------------------------------------------------------------------
  // 3) New Jacobian contribution from the grad(κ) term
  //
  // Residual added: – (∇κ · ∇u) * L * test
  // must linearize this:
  //   a) ∇u → ∇φ_j
  //   b) L  → dL/dη
  // --------------------------------------------------------------------------
  RankTwoTensor A = (_use_anisotropic_matrix && _A_ptr) ? (*_A_ptr)[_qp] : RankTwoTensor::initIdentity;// ✅ 2025/07/06
  // build ∇κ vector from coupled AuxVariables
  RealGradient gk = gradKappa();

  // derivative of ∇u term: gk·∇φ_j
  Real d_res_grad_u = gk* (A * _grad_phi[_j][_qp]) * _L[_qp];// ✅ 2025/07/06

  // derivative of L term: gk·∇u  times dL/dη
  Real d_res_L      = gk * (A * _grad_u[_qp])     * _dLdop[_qp];// ✅ 2025/07/06

  // combine both and apply the test function
  Real extra = -(d_res_grad_u + d_res_L) * _test[_i][_qp];

  // --------------------------------------------------------------------------
  // 4) assemble full Jacobian: original + new grad(κ) part
  // --------------------------------------------------------------------------
  
  return (_grad_phi[_j][_qp]) * (A * kappaNablaLPsi())  // from original term// ✅ 2025/07/06
       + _grad_u[_qp]       * dsum              // from original term
       + extra;                                 // from grad(κ) term
}

Real
ACInterfaceGradKappa::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Map the Jacobian index back to the coupled-variable index
  const unsigned int cvar = mapJvarToCvar(jvar);

  // --------------------------------------------------------------------------
  // 1) Linearization of the original κ ∇(L ψ) · ∇ test term w.r.t. coupled var
  // --------------------------------------------------------------------------
  RealGradient dsum =
      ((*_dkappadarg[cvar])[_qp] * _L[_qp]                          // ∂κ/∂arg * L
     + _kappa[_qp] * (*_dLdarg[cvar])[_qp])                        // + κ * ∂L/∂arg
      * _phi[_j][_qp]                                              //   * φ_j
      * _grad_test[_i][_qp];                                       //   · ∇ test

  // If L depends on the solution, include its second‐derivative terms
  if (_variable_L)
  {
    // d(∇L)/d arg = ∇φ_j * ∂L/∂arg + ∇u * φ_j * ∂²L/(∂ arg ∂ η)
    RealGradient dgradL =
        _grad_phi[_j][_qp] * (*_dLdarg[cvar])[_qp]
      + _grad_u[_qp]     * _phi[_j][_qp] * (*_d2Ldargdop[cvar])[_qp];

    // Cross‐derivatives ∂²L/(∂ arg ∂ other_args)
    for (unsigned int i = 0; i < _n_args; ++i)
      dgradL += (*_gradarg[i])[_qp]
              * _phi[_j][_qp]
              * (*_d2Ldarg2[cvar][i])[_qp];

    // Add: κ * d(∇L)/d arg + ∂κ/∂η * φ_j * ∇L
    dsum += (_kappa[_qp] * dgradL
          + _dkappadop[_qp] * _phi[_j][_qp] * gradL())             // gradL() = ∇(L ψ)
          * _test[_i][_qp];
  }

  // Original off‐diagonal block
  Real jac = _grad_u[_qp] * dsum;

  // --------------------------------------------------------------------------
  // 2) New grad(κ) term: residual had –(∇κ·∇u) * L * test
  //    Its derivative w.r.t. a coupled var 'arg' only enters via L
  // --------------------------------------------------------------------------
  RealGradient gk = gradKappa();
  RankTwoTensor A = (_use_anisotropic_matrix && _A_ptr) ? (*_A_ptr)[_qp] : RankTwoTensor::initIdentity;// ✅ 2025/07/06
  
  Real extra = -(
      (gk * (A * _grad_u[_qp]))       // ∇κ·(A·∇u)
    * (*_dLdarg[cvar])[_qp]           // ∂L/∂arg
    * _test[_i][_qp]
  );// ✅ 2025/07/06

  return jac + extra;
}
