/* ------------------------------------------------------------------ */
/*  test_api_sentinel.cpp                                              */
/*                                                                    */
/*  Extended sentinel / boundary guard tests for new APIs:             */
/*    - Buffer boundary guards                                         */
/*    - Null-pointer sentinels on extended functions                   */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v; v.x = x; v.y = y; v.z = z;
    return v;
}

/* ================================================================== */
TEST(Sentinel, BladeGeoNullSetters) {
    // All setters with nullptr should not crash
    oc_blade_geometry_set_twist(nullptr, nullptr, 0);
    oc_blade_geometry_set_chord(nullptr, nullptr, 0);
    oc_blade_geometry_set_radius(nullptr, nullptr, 0);
    oc_blade_geometry_set_C_l_alpha(nullptr, nullptr, 0);
    oc_blade_geometry_set_alpha_0(nullptr, nullptr, 0);
    oc_blade_geometry_set_sweep(nullptr, nullptr, 0);
}

/* ================================================================== */
TEST(Sentinel, BladeGeoNullSetters2) {
    oc_blade_geometry_set_xi(nullptr, nullptr, 0);
    oc_blade_geometry_set_thickness(nullptr, nullptr, 0);
    oc_blade_geometry_set_xi_p(nullptr, nullptr, 0);
    oc_blade_geometry_compute_vectors(nullptr);
}

/* ================================================================== */
TEST(Sentinel, RotorGeoNullSetters) {
    oc_rotor_geometry_set_solidity(nullptr, 0.0);
    oc_rotor_geometry_set_frame(nullptr, nullptr);
}

/* ================================================================== */
TEST(Sentinel, AircraftNullSetters) {
    oc_aircraft_set_rotors(nullptr, nullptr, 0);
}

/* ================================================================== */
TEST(Sentinel, FrameNullSafe) {
    oc_frame_destroy(nullptr);
}

/* ================================================================== */
TEST(Sentinel, Mat3Mat4Helpers) {
    OC_Mat3 m3 = oc_mat3_identity();
    // data is flat [row*3+col], diagonal at indices 0, 4, 8
    EXPECT_NEAR(m3.data[0], 1.0, 1e-6);
    EXPECT_NEAR(m3.data[4], 1.0, 1e-6);

    OC_Mat4 m4 = oc_mat4_identity();
    // data is flat [row*4+col], diagonal at indices 0, 5, 10, 15
    EXPECT_NEAR(m4.data[0], 1.0, 1e-6);
    EXPECT_NEAR(m4.data[15], 1.0, 1e-6);
}

/* ================================================================== */
TEST(Sentinel, ChunkSizePositive) {
    size_t cs = oc_chunk_size();
    EXPECT_GT(cs, (size_t)0);
}

/* ================================================================== */
TEST(Sentinel, RotorInputNullSafe) {
    // Null-pointer sentinels on RotorInputState setters/getters
    EXPECT_DOUBLE_EQ(oc_rotor_input_get_angular_velocity(nullptr), 0.0);
    oc_rotor_input_set_angular_velocity(nullptr, 999.0);
}