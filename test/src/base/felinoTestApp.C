//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "felinoTestApp.h"
#include "felinoApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
felinoTestApp::validParams()
{
  InputParameters params = felinoApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

felinoTestApp::felinoTestApp(InputParameters parameters) : MooseApp(parameters)
{
  felinoTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

felinoTestApp::~felinoTestApp() {}

void
felinoTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  felinoApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"felinoTestApp"});
    Registry::registerActionsTo(af, {"felinoTestApp"});
  }
}

void
felinoTestApp::registerApps()
{
  registerApp(felinoApp);
  registerApp(felinoTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
felinoTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  felinoTestApp::registerAll(f, af, s);
}
extern "C" void
felinoTestApp__registerApps()
{
  felinoTestApp::registerApps();
}
