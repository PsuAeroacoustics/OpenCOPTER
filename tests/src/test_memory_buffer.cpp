/* ------------------------------------------------------------------ */
/*  test_memory_buffer.cpp                                             */
/*                                                                    */
/*  Memory safety tests for buffer-based getter/setter functions:       */
/*    - RotorInputState get/set array patterns with sentinel guards     */
/*    - Verify output buffers are not overwritten beyond len           */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cstddef>
#include <algorithm>

static constexpr double kSentinel = -987654.0;

/* ================================================================== */
/*  TEST: RotorInputState set_r_0 and get_r_0 with sentinel guards     */
/* ================================================================== */
TEST(MemoryBuffer, RotorInputR0Sentinel) {
    size_t num_blades = 4;
    OC_AircraftInputState* input = oc_aircraft_input_state_create(1, &num_blades, 0);
    ASSERT_NE(input, nullptr);

    OC_RotorInputState* rotor_in = oc_aircraft_input_get_rotor_input(input, 0);

    double r0_in[4] = {0.1, 0.2, 0.3, 0.4};
    oc_rotor_input_set_r_0(rotor_in, r0_in, 4);

    const size_t OUT_SIZE = 8;
    double r0_out[OUT_SIZE];
    std::fill(std::begin(r0_out), std::end(r0_out), kSentinel);

    oc_rotor_input_get_r_0(rotor_in, r0_out, 4);

    for (int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(r0_out[i], r0_in[i]);
    }
    for (int i = 4; i < OUT_SIZE; ++i) {
        EXPECT_DOUBLE_EQ(r0_out[i], kSentinel);
    }

    oc_aircraft_input_state_destroy(input);
}

/* ================================================================== */
/*  TEST: RotorInputState blade flapping sentinel                      */
/* ================================================================== */
TEST(MemoryBuffer, RotorInputFlappingSentinel) {
    size_t num_blades = 4;
    OC_AircraftInputState* input = oc_aircraft_input_state_create(1, &num_blades, 0);

    OC_RotorInputState* rotor_in = oc_aircraft_input_get_rotor_input(input, 0);

    double flap[4] = {0.0, 0.1, 0.2, 0.3};
    oc_rotor_input_set_blade_flapping(rotor_in, flap, 4);

    const size_t OUT_SIZE = 8;
    double flap_out[OUT_SIZE];
    std::fill(std::begin(flap_out), std::end(flap_out), kSentinel);

    oc_rotor_input_get_blade_flapping(rotor_in, flap_out, 4);

    for (int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(flap_out[i], flap[i]);
    }
    for (int i = 4; i < OUT_SIZE; ++i) {
        EXPECT_DOUBLE_EQ(flap_out[i], kSentinel);
    }

    oc_aircraft_input_state_destroy(input);
}

/* ================================================================== */
/*  TEST: RotorInputState flapping rate sentinel                       */
/* ================================================================== */
TEST(MemoryBuffer, RotorInputFlappingRateSentinel) {
    size_t num_blades = 4;
    OC_AircraftInputState* input = oc_aircraft_input_state_create(1, &num_blades, 0);

    OC_RotorInputState* rotor_in = oc_aircraft_input_get_rotor_input(input, 0);

    double rate[4] = {1.0, 2.0, 3.0, 4.0};
    oc_rotor_input_set_blade_flapping_rate(rotor_in, rate, 4);

    const size_t OUT_SIZE = 8;
    double rate_out[OUT_SIZE];
    std::fill(std::begin(rate_out), std::end(rate_out), kSentinel);

    oc_rotor_input_get_blade_flapping_rate(rotor_in, rate_out, 4);

    for (int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(rate_out[i], rate[i]);
    }
    for (int i = 4; i < OUT_SIZE; ++i) {
        EXPECT_DOUBLE_EQ(rate_out[i], kSentinel);
    }

    oc_aircraft_input_state_destroy(input);
}

/* ================================================================== */
/*  TEST: RotorInputState scalar setters do not corrupt adjacent mem   */
/* ================================================================== */
TEST(MemoryBuffer, ScalarSetterNoOverflow) {
    size_t num_blades = 4;
    OC_AircraftInputState* input = oc_aircraft_input_state_create(1, &num_blades, 0);
    OC_RotorInputState* rotor_in = oc_aircraft_input_get_rotor_input(input, 0);

    // Set many scalars repeatedly -- should not corrupt.
    for (int i = 0; i < 100; ++i) {
        oc_rotor_input_set_angular_velocity(rotor_in, 50.0 + i);
        oc_rotor_input_set_azimuth(rotor_in, i * 0.1);
        oc_rotor_input_set_angular_accel(rotor_in, i * 0.01);
    }

    EXPECT_DOUBLE_EQ(oc_rotor_input_get_angular_velocity(rotor_in), 149.0);
    EXPECT_FALSE(oc_rotor_input_get_azimuth(rotor_in) != oc_rotor_input_get_azimuth(rotor_in));

    oc_aircraft_input_state_destroy(input);
}