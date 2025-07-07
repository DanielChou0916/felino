#pragma once

// MOOSE includes
#include "Action.h"
#include "AddVariableAction.h"

#include "libmesh/fe_type.h"

class PFNonconservedAction : public Action
{
public:
  static InputParameters validParams();

  PFNonconservedAction(const InputParameters & params);

  virtual void act();

protected:
  /// Name of the variable being created
  const NonlinearVariableName _var_name;
  /// FEType for the variable being created
  const FEType _fe_type;
};
