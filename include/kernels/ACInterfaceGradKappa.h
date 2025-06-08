#pragma once

#include "Kernel.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"

/**
 * Compute the Allen-Cahn interface term with the weak form residual
 * \f$ \left( \kappa_i \nabla\eta_i, \nabla (L_i \psi) \right) \f$
 */
class ACInterfaceGradKappa : public DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>
{
public:
  static InputParameters validParams();

  ACInterfaceGradKappa(const InputParameters & parameters);
  virtual void initialSetup();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  RealGradient gradL();
  RealGradient gradKappa();

  /// the \f$ \nabla(L\psi) \f$ term
  RealGradient nablaLPsi();

  /// the \f$ \kappa\nabla(L\psi) \f$ term
  RealGradient kappaNablaLPsi();

  /// Mobility
  const MaterialProperty<Real> & _L;
  /// Interfacial parameter
  const MaterialProperty<Real> & _kappa;

  /// flag set if L is a function of non-linear variables in args
  const bool _variable_L;

  /// @{ Mobility derivatives w.r.t. order parameter
  const MaterialProperty<Real> & _dLdop;
  const MaterialProperty<Real> & _d2Ldop2;
  /// @}

  /// kappa derivative w.r.t. order parameter
  const MaterialProperty<Real> & _dkappadop;

  /// @{ Mobility derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dLdarg;
  std::vector<const MaterialProperty<Real> *> _d2Ldargdop;
  std::vector<std::vector<const MaterialProperty<Real> *>> _d2Ldarg2;
  /// @}

  /// kappa derivative w.r.t. other coupled variables
  std::vector<const MaterialProperty<Real> *> _dkappadarg;

  /// Gradients for all coupled variables
  std::vector<const VariableGradient *> _gradarg;
  /// gradient of kappa in each direction
  const VariableValue & _grad_kappa_x;
  const VariableValue & _grad_kappa_y;
  const VariableValue * _grad_kappa_z;
};
