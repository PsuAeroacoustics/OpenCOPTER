module opencopter.liftmodels;

public import opencopter.liftmodels.steadylift;

import std.meta;

alias lift_model_list = aliasSeqOf!(["steady_lift_model"]);
