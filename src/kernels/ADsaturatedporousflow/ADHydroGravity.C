

#include "ADHydroDiffusion.h"

registerMooseObject("MooseApp", ADHydroDiffusion);

InputParameters
ADHydroDiffusion::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addParam<MaterialPropertyName>("fluid_density", "fluid_density", "Name of fluid density as an input.");
  params.addParam<MaterialPropertyName>("fluid_viscosity", "fluid_viscosity", "Name of fluid viscosity as an input.");
  params.addParam<MaterialPropertyName>("permeability_tensor", "permeability_tensor", "Name of permeability tensor as an input.");
  return params;
}

ADHydroDiffusion::ADHydroDiffusion(const InputParameters & parameters) 
: ADKernel(parameters), 
  _fluid_density(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("fluid_density") )),
  _fluid_viscosity(getADMaterialProperty<Real>(getParam<MaterialPropertyName>("fluid_viscosity") )),
  _permeability_tensor(getADMaterialProperty<RankTwoTensor>(getParam<MaterialPropertyName>("permeability_tensor") ))
{ 
}

ADReal
ADHydroDiffusion::computeQpResidual()
{
  ADRankTwoTensor prop = _fluid_density[_qp]/_fluid_viscosity[_qp] * _permeability_tensor[_qp];
  return _grad_test[_i][_qp] * (prop * _grad_u[_qp]);
}

