#include "felinoApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
felinoApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

felinoApp::felinoApp(InputParameters parameters) : MooseApp(parameters)
{
  felinoApp::registerAll(_factory, _action_factory, _syntax);
}

felinoApp::~felinoApp() {}

void
felinoApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<felinoApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"felinoApp"});
  Registry::registerActionsTo(af, {"felinoApp"});

  /* register custom execute flags, action syntax, etc. here */
  syntax.registerActionSyntax("ADNonconservedAction", "Actions/ADNonconserved/*"); // ✅
  syntax.registerActionSyntax("NonconservedActionGrad", "Actions/NonconservedG/*"); // ✅
  syntax.registerActionSyntax("PFNonconservedAction", "Actions/PFNonconserved/*"); // ✅
}

void
felinoApp::registerApps()
{
  registerApp(felinoApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
felinoApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  felinoApp::registerAll(f, af, s);
}
extern "C" void
felinoApp__registerApps()
{
  felinoApp::registerApps();
}
