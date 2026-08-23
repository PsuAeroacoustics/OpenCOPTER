#include "opencopter.hpp"
#include "opencopter.h"
#include <gtest/gtest.h>
#include <cmath>

using namespace opencopter;

namespace {

// Helper: create a 1-rotor, 2-blade, 16-chunk input state
class BladeInputTestFixture : public ::testing::Test {
protected:
    AircraftInputState input;
    RotorInputState rotor;

    void SetUp() override {
        input = AircraftInputState(1, {2}, 0, {16});
        rotor = input.get_rotor_input(0);
    }
};

// ========================================================================
// Per-blade scalar round-trip tests
// ========================================================================

TEST_F(BladeInputTestFixture, BladeInputPitchRoundTrip) {
    double val = 0.5;
    rotor.set_blade_input_pitch(0, val);
    EXPECT_DOUBLE_EQ(rotor.get_blade_input_pitch(0), val);
}

TEST_F(BladeInputTestFixture, BladeInputPitchBlade1RoundTrip) {
    double val = -0.3;
    rotor.set_blade_input_pitch(1, val);
    EXPECT_DOUBLE_EQ(rotor.get_blade_input_pitch(1), val);
}

TEST_F(BladeInputTestFixture, BladeInputFlappingRoundTrip) {
    double val = 0.12;
    rotor.set_blade_input_flapping(0, val);
    EXPECT_DOUBLE_EQ(rotor.get_blade_input_flapping(0), val);
}

TEST_F(BladeInputTestFixture, BladeInputFlappingRateRoundTrip) {
    double val = 1.5;
    rotor.set_blade_input_flapping_rate(0, val);
    EXPECT_DOUBLE_EQ(rotor.get_blade_input_flapping_rate(0), val);
}

TEST_F(BladeInputTestFixture, BladeInputR0RoundTrip) {
    double val = 0.05;
    rotor.set_blade_input_r0(1, val);
    EXPECT_DOUBLE_EQ(rotor.get_blade_input_r0(1), val);
}

TEST_F(BladeInputTestFixture, BladeInputScalarsOOBSafe) {
    // Out-of-bounds blade index should not crash
    rotor.set_blade_input_pitch(99, 0.5);
    // Should return default (infinity sentinel = inf)
    double r = rotor.get_blade_input_pitch(99);
    EXPECT_TRUE(std::isinf(r) || r == 0.0);
}

// ========================================================================
// Per-station array span tests (zero-copy writable spans)
// ========================================================================

TEST_F(BladeInputTestFixture, FlapDeflectionSpanSize) {
    auto span = rotor.blade_flap_deflection(0);
    EXPECT_EQ(span.size(), 128u);  // 16 chunks * 8 doubles
    EXPECT_NE(span.data(), nullptr);
}

TEST_F(BladeInputTestFixture, FlapDeflectionSpanReadWrite) {
    auto span = rotor.blade_flap_deflection(0);
    // Write through span
    span[0] = 0.1;
    span[8] = 0.2;
    span[127] = 0.3;
    // Read back
    EXPECT_DOUBLE_EQ(span[0], 0.1);
    EXPECT_DOUBLE_EQ(span[8], 0.2);
    EXPECT_DOUBLE_EQ(span[127], 0.3);
    // Verify via fresh span (same backing memory)
    auto span2 = rotor.blade_flap_deflection(0);
    EXPECT_DOUBLE_EQ(span2[0], 0.1);
    EXPECT_DOUBLE_EQ(span2[8], 0.2);
    EXPECT_DOUBLE_EQ(span2[127], 0.3);
}

TEST_F(BladeInputTestFixture, LagDeflectionSpanReadWrite) {
    auto span = rotor.blade_lag_deflection(1);
    EXPECT_EQ(span.size(), 128u);
    span[4] = -0.05;
    span[100] = 0.7;
    auto span2 = rotor.blade_lag_deflection(1);
    EXPECT_DOUBLE_EQ(span2[4], -0.05);
    EXPECT_DOUBLE_EQ(span2[100], 0.7);
}

TEST_F(BladeInputTestFixture, TwistDeflectionSpanReadWrite) {
    auto span = rotor.blade_twist_deflection(0);
    EXPECT_EQ(span.size(), 128u);
    span[4] = M_PI / 12.0;
    auto span2 = rotor.blade_twist_deflection(0);
    EXPECT_DOUBLE_EQ(span2[4], M_PI / 12.0);
}

TEST_F(BladeInputTestFixture, FlapVelocitySpanReadWrite) {
    auto span = rotor.blade_flap_velocity(0);
    EXPECT_EQ(span.size(), 128u);
    for (size_t i = 0; i < 128; ++i) span[i] = 0.001 * static_cast<double>(i);
    auto span2 = rotor.blade_flap_velocity(0);
    EXPECT_DOUBLE_EQ(span2[0], 0.0);
    EXPECT_DOUBLE_EQ(span2[127], 0.001 * 127.0);
}

TEST_F(BladeInputTestFixture, LagVelocitySpanReadWrite) {
    auto span = rotor.blade_lag_velocity(1);
    EXPECT_EQ(span.size(), 128u);
    std::fill(span.begin(), span.end(), 0.5);
    auto span2 = rotor.blade_lag_velocity(1);
    for (size_t i = 0; i < 128; ++i) EXPECT_DOUBLE_EQ(span2[i], 0.5);
}

TEST_F(BladeInputTestFixture, SpanZeroCopySameBacking) {
    // Verify that two spans to the same blade share the same backing memory
    auto s1 = rotor.blade_flap_deflection(0);
    auto s2 = rotor.blade_flap_deflection(0);
    EXPECT_EQ(s1.data(), s2.data());
}

TEST_F(BladeInputTestFixture, SpanDifferentBladesDifferentBacking) {
    auto s0 = rotor.blade_flap_deflection(0);
    auto s1 = rotor.blade_flap_deflection(1);
    EXPECT_NE(s0.data(), s1.data());
}

// ========================================================================
// Null-safety tests
// ========================================================================

TEST_F(BladeInputTestFixture, NullRotorSpansThrow) {
    RotorInputState null_rotor;
    EXPECT_THROW(null_rotor.blade_flap_deflection(0), std::runtime_error);
    EXPECT_THROW(null_rotor.blade_lag_deflection(0), std::runtime_error);
    EXPECT_THROW(null_rotor.blade_twist_deflection(0), std::runtime_error);
    EXPECT_THROW(null_rotor.blade_flap_velocity(0), std::runtime_error);
    EXPECT_THROW(null_rotor.blade_lag_velocity(0), std::runtime_error);
}

TEST_F(BladeInputTestFixture, NullScalarSettersThrow) {
    RotorInputState null_rotor;
    EXPECT_THROW(null_rotor.set_blade_input_pitch(0, 1.0), std::runtime_error);
    EXPECT_DOUBLE_EQ(null_rotor.get_blade_input_pitch(0), 0.0);
}

// ========================================================================
// C-level API tests (using opencopter.h directly)
// ========================================================================

TEST_F(BladeInputTestFixture, CLevelCreateWithChunks) {
    size_t num_blades[1] = {3};
    size_t num_chunks[1] = {8};
    auto* raw = oc_aircraft_input_state_create_with_chunks(1, num_blades, 0, num_chunks);
    ASSERT_NE(raw, nullptr);

    auto* rotor = oc_aircraft_input_get_rotor_input(raw, 0);
    ASSERT_NE(rotor, nullptr);

    // Test scalar set/get
    oc_rotor_input_set_blade_input_pitch(rotor, 0, 0.7);
    EXPECT_DOUBLE_EQ(oc_rotor_input_get_blade_input_pitch(rotor, 0), 0.7);

    // Test span access via _ref
    size_t len = 0;
    double* ptr = oc_rotor_input_get_blade_flap_deflection_ref(rotor, 0, &len);
    ASSERT_NE(ptr, nullptr);
    EXPECT_EQ(len, 64u);  // 8 chunks * 8
    for (size_t i = 0; i < 64; ++i) ptr[i] = 0.01;
    // Read back
    double* ptr2 = oc_rotor_input_get_blade_flap_deflection_ref(rotor, 0, &len);
    for (size_t i = 0; i < 64; ++i) {
        EXPECT_DOUBLE_EQ(ptr2[i], 0.01) << "C-level mismatch at index " << i;
    }

    oc_aircraft_input_state_destroy(raw);
}

TEST_F(BladeInputTestFixture, CLevelOOBSafe) {
    size_t num_blades[1] = {1};
    size_t num_chunks[1] = {4};
    auto* raw = oc_aircraft_input_state_create_with_chunks(1, num_blades, 0, num_chunks);
    ASSERT_NE(raw, nullptr);

    auto* rotor = oc_aircraft_input_get_rotor_input(raw, 0);
    ASSERT_NE(rotor, nullptr);

    // OOB blade index should not crash
    oc_rotor_input_set_blade_input_pitch(rotor, 99, 1.0);
    double r = oc_rotor_input_get_blade_input_pitch(rotor, 99);
    EXPECT_TRUE(std::isinf(r) || r == 0.0);

    // OOB _ref should return null
    size_t len = 999;
    double* ptr = oc_rotor_input_get_blade_flap_deflection_ref(rotor, 99, &len);
    EXPECT_EQ(ptr, nullptr);
    EXPECT_EQ(len, 0u);

    oc_aircraft_input_state_destroy(raw);
}

}  // namespace