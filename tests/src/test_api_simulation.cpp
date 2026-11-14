/* ------------------------------------------------------------------ */
/*  test_api_simulation.cpp                                            */
/*                                                                    */
/*  Tests for Simulation / dynamics API:                               */
/*    - basic_single_rotor_dynamics null-safe                          */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>

/* ================================================================== */
TEST(Simulation, BasicSingleRotorDynamicsNullSafe) {
    double result = oc_basic_single_rotor_dynamics(nullptr, 0.01);
    // With null input, expect 0.0 or some safe default
    EXPECT_NEAR(result, 0.0, 1e-6);
}

/* ================================================================== */
TEST(Simulation, BasicSingleRotorDynamicsValid) {
    // Build a minimal RotorInputState via AircraftInputState
    size_t num_blades_arr[1] = {2};
    OC_AircraftInputState* ac_input = oc_aircraft_input_state_create(1, num_blades_arr, 0);
    ASSERT_NE(ac_input, nullptr);

    OC_RotorInputState* rotor_input = oc_aircraft_input_get_rotor_input(ac_input, 0);
    ASSERT_NE(rotor_input, nullptr);

    // Initialize all required fields so dynamics has valid input state
    oc_rotor_input_set_angular_velocity(rotor_input, 100.0);
    oc_rotor_input_set_angular_accel(rotor_input, 0.0);
    oc_rotor_input_set_azimuth(rotor_input, 0.0);

    double dt_result = oc_basic_single_rotor_dynamics(rotor_input, 0.01);
    // The result should not be NaN
    EXPECT_FALSE(std::isnan(dt_result));

    oc_aircraft_input_state_destroy(ac_input);
}

/* ================================================================== */
TEST(Simulation, DirectionHelpers) {
    OC_Direction cw = oc_direction_clockwise();
    OC_Direction ccw = oc_direction_counter_clockwise();
    // Just verify they don't crash; values depend on implementation
    EXPECT_NE(cw, ccw);
}