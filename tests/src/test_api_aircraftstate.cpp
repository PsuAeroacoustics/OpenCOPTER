/* ------------------------------------------------------------------ */
/*  test_api_aircraftstate.cpp                                         */
/*                                                                    */
/*  Tests for AircraftState API:                                       */
/*    - Null-safe destroy                                             */
/*    - Rotor CT/CQ null-safe queries                                 */
/*    - Freestream set/get (full chain if constructible)               */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v; v.x = x; v.y = y; v.z = z;
    return v;
}

/* ================================================================== */
TEST(AircraftState, NullSafeDestroy) {
    oc_aircraft_state_destroy(nullptr); /* should not crash */
}

/* ================================================================== */
TEST(AircraftState, RotorCTNullSafe) {
    double out_val;
    // Passing nullptr state: function may do nothing or return 0
    oc_aircraft_state_get_rotor_C_T(nullptr, 0, &out_val);
}

/* ================================================================== */
TEST(AircraftState, RotorCQNullSafe) {
    double out_val;
    oc_aircraft_state_get_rotor_C_Q(nullptr, 0, &out_val);
}