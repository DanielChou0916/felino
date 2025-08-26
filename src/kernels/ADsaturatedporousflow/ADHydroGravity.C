

#include "ADHydroGravity.h"

registerMooseObject("MooseApp", ADHydroGravity);

InputParameters
ADHydroGravity::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Hydro-Gravity Term");
  params.addRequiredParam<RealVectorValue>("gravity",
                                           "Gravitational acceleration vector(m/s^2)");
  params.addParam<MaterialPropertyName>("fluid_density", "fluid_density", "Name of fluid density as an input.");
  params.addParam<MaterialPropertyName>("fluid_viscosity", "fluid_viscosity", "Name of fluid viscosity as an input.");
  params.addParam<MaterialPropertyName>("permeability_tensor", "permeability_tensor", "Name of permeability tensor as an input.");
  params.addParam<bool>(
  "mass_flux",
  false,
  "Whether to multiply density into porousflow governing equation."
  );
  return params;
}

ADHydroGravity::ADHydroGravity(const InputParameters & parameters) 
: ADKernel(parameters), 
  _fluid_density(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("fluid_density") )),
  _fluid_viscosity(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("fluid_viscosity") )),
  _permeability_tensor(getADMaterialProperty<RankTwoTensor>(getParam<MaterialPropertyName>("permeability_tensor") )),
  _gravity(getParam<RealVectorValue>("gravity")),
  _mass_flux(getParam<bool>("mass_flux"))
{ 
}

ADReal
ADHydroGravity::computeQpResidual()
{
  ADRankTwoTensor prop;
  if (_mass_flux)
    prop = _fluid_density[_qp]/_fluid_viscosity[_qp] * _permeability_tensor[_qp];
  else
    prop = _permeability_tensor[_qp]/_fluid_viscosity[_qp];

  return _grad_test[_i][_qp] * (prop * _fluid_density[_qp]*_gravity);
}

