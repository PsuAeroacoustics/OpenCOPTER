/* ------------------------------------------------------------------ */
/*  test_memory_aircraft.cpp                                           */
/*                                                                    */
/*  Memory safety tests for the Aircraft composite API:                */
/*    - Aircraft with rotors                                          */
/*    - oc_aircraft_set_rotors memory handling                        */
/*    - AircraftInputState setter / getter patterns                    */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cstddef>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v = {x, y, z};
    return v;
}

/* ================================================================== */
/*  TEST: Aircraft with rotors full lifecycle                          */
/* ================================================================== */
TEST(MemoryAircraft, AircraftWithRotors) {
    OC_Aircraft* ac = oc_aircraft_create(1, 0);
    ASSERT_NE(ac, nullptr);

    OC_RotorGeometry* rotor = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 2.0, 0.5);
    ASSERT_NE(rotor, nullptr);

    OC_RotorGeometry* rotors[1] = {rotor};
    oc_aircraft_set_rotors(ac, rotors, 1);

    oc_rotor_geometry_destroy(rotor);
    oc_aircraft_destroy(ac);
}

/* ================================================================== */
/*  TEST: Aircraft with multiple rotor slots                           */
/* ================================================================== */
TEST(MemoryAircraft, AircraftMultipleRotors) {
    const int NUM_ROTORS = 4;
    OC_Aircraft* ac = oc_aircraft_create(NUM_ROTORS, 0);
    ASSERT_NE(ac, nullptr);

    OC_RotorGeometry* rotors[NUM_ROTORS] = {};
    for (int i = 0; i < NUM_ROTORS; ++i) {
        rotors[i] = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 2.0, 0.5);
        ASSERT_NE(rotors[i], nullptr);
    }

    oc_aircraft_set_rotors(ac, rotors, NUM_ROTORS);

    for (int i = 0; i < NUM_ROTORS; ++i) {
        oc_rotor_geometry_destroy(rotors[i]);
    }
    oc_aircraft_destroy(ac);
}

/* ================================================================== */
/*  TEST: Aircraft with wings                                          */
/* ================================================================== */
TEST(MemoryAircraft, AircraftWithWings) {
    OC_Aircraft* ac = oc_aircraft_create(0, 1);
    ASSERT_NE(ac, nullptr);

    OC_WingGeometry* wing = oc_wing_geometry_create(1, make_vec3(0, 0, 0), 4.0);
    ASSERT_NE(wing, nullptr);

    oc_wing_geometry_destroy(wing);
    oc_aircraft_destroy(ac);
}

/* ================================================================== */
/*  TEST: get_root_frame returns non-null                              */
/* ================================================================== */
TEST(MemoryAircraft, GetRootFrame) {
    OC_Aircraft* ac = oc_aircraft_create(0, 0);
    ASSERT_NE(ac, nullptr);

    OC_Frame* root = oc_aircraft_get_root_frame(ac);
    EXPECT_NE(root, nullptr);

    oc_aircraft_destroy(ac);
}

/* ================================================================== */
/*  TEST: AircraftInputState with RotorInputState setters              */
/* ================================================================== */
TEST(MemoryAircraft, InputStateSettersAndGetters) {
    size_t num_blades = 4;
    OC_AircraftInputState* input = oc_aircraft_input_state_create(1, &num_blades, 0);
    ASSERT_NE(input, nullptr);

    OC_RotorInputState* rotor_in = oc_aircraft_input_get_rotor_input(input, 0);
    ASSERT_NE(rotor_in, nullptr);

    oc_rotor_input_set_angular_velocity(rotor_in, 109.12);
    EXPECT_DOUBLE_EQ(oc_rotor_input_get_angular_velocity(rotor_in), 109.12);

    oc_rotor_input_set_angular_accel(rotor_in, 0.0);
    EXPECT_DOUBLE_EQ(oc_rotor_input_get_angular_accel(rotor_in), 0.0);

    oc_rotor_input_set_azimuth(rotor_in, 1.57);
    EXPECT_DOUBLE_EQ(oc_rotor_input_get_azimuth(rotor_in), 1.57);

    double r0_data[4] = {0.1, 0.2, 0.3, 0.4};
    oc_rotor_input_set_r_0(rotor_in, r0_data, 4);
    double r0_out[4] = {};
    oc_rotor_input_get_r_0(rotor_in, r0_out, 4);
    for (int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(r0_out[i], r0_data[i]);
    }

    double flap[4] = {0.0, 0.0, 0.0, 0.0};
    oc_rotor_input_set_blade_flapping(rotor_in, flap, 4);
    double flap_out[4] = {};
    oc_rotor_input_get_blade_flapping(rotor_in, flap_out, 4);

    double fRate[4] = {0.1, 0.2, 0.3, 0.4};
    oc_rotor_input_set_blade_flapping_rate(rotor_in, fRate, 4);
    double fRate_out[4] = {};
    oc_rotor_input_get_blade_flapping_rate(rotor_in, fRate_out, 4);

    oc_aircraft_input_set_blade_pitch(input, 0, 0, 0.5);
    EXPECT_DOUBLE_EQ(oc_aircraft_input_get_blade_pitch(input, 0, 0), 0.5);

    oc_aircraft_input_state_destroy(input);
}

/* ================================================================== */
/*  TEST: Multiple create/destroy cycles                               */
/* ================================================================== */
TEST(MemoryAircraft, MultipleCycles) {
    for (int i = 0; i < 10; ++i) {
        OC_Aircraft* ac = oc_aircraft_create(2, 1);
        ASSERT_NE(ac, nullptr);

        OC_RotorGeometry* r0 = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 2.0, 0.5);
        OC_RotorGeometry* r1 = oc_rotor_geometry_create(4, make_vec3(1, 0, 0), 2.0, 0.5);
        OC_RotorGeometry* rotors[2] = {r0, r1};
        oc_aircraft_set_rotors(ac, rotors, 2);

        OC_WingGeometry* w = oc_wing_geometry_create(1, make_vec3(0, 0, 0), 4.0);
        oc_wing_geometry_destroy(w);

        oc_rotor_geometry_destroy(r0);
        oc_rotor_geometry_destroy(r1);
        oc_aircraft_destroy(ac);
    }
}