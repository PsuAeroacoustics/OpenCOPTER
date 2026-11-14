/* ------------------------------------------------------------------ */
/*  test_api_config.cpp                                                */
/*                                                                    */
/*  Tests for configuration and value-type APIs:                       */
/*    - oc_chunk_size()                                                */
/*    - oc_mat3_identity()                                             */
/*    - Struct layout verification                                     */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>

/* ================================================================== */
/*  TEST: Chunk size returns a positive value                         */
/* ================================================================== */
TEST(APIStructLayout, ChunkSize) {
    size_t cs = oc_chunk_size();
    EXPECT_GT(cs, static_cast<size_t>(0));
}

/* ================================================================== */
/*  TEST: Mat3 identity has correct diagonal                          */
/* ================================================================== */
TEST(APIHelpers, Mat3Identity) {
    OC_Mat3 m = oc_mat3_identity();
    // Diagonal elements should be 1.0
    for (int i = 0; i < 3; ++i) {
        EXPECT_DOUBLE_EQ(m.data[i * 3 + i], 1.0);
    }
    // Off-diagonal elements should be 0.0
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            if (r != c) EXPECT_DOUBLE_EQ(m.data[r * 3 + c], 0.0);
        }
    }
}

/* ================================================================== */
/*  TEST: ValueType sizes are reasonable                              */
/* ================================================================== */
TEST(APIStructLayout, ValueTypeSizes) {
    EXPECT_GT(sizeof(OC_Vec3), static_cast<size_t>(0));
    EXPECT_GT(sizeof(OC_Vec4), static_cast<size_t>(0));
    EXPECT_GT(sizeof(OC_Mat3), static_cast<size_t>(0));
    EXPECT_GT(sizeof(OC_Mat4), static_cast<size_t>(0));
    EXPECT_GT(sizeof(OC_Atmosphere), static_cast<size_t>(0));
    EXPECT_GT(sizeof(OC_InducedVelocities), static_cast<size_t>(0));

    // Verify OC_Vec3 has 3 doubles
    EXPECT_EQ(sizeof(OC_Vec3), 3 * sizeof(double));
    // Verify OC_Vec4 has 4 doubles
    EXPECT_EQ(sizeof(OC_Vec4), 4 * sizeof(double));
    // Verify OC_Mat3 has 9 doubles
    EXPECT_EQ(sizeof(OC_Mat3), 9 * sizeof(double));
    // Verify OC_Mat4 has 16 doubles
    EXPECT_EQ(sizeof(OC_Mat4), 16 * sizeof(double));
}

/* ================================================================== */
/*  TEST: Atmosphere and InducedVelocities layout                     */
/* ================================================================== */
TEST(APIStructLayout, AtmosphereAndInducedVelocities) {
    OC_Atmosphere atmo;
    atmo.density = 1.225;
    atmo.dynamic_viscosity = 0.0000181;
    atmo.kinematic_viscosity = 0.0000148;
    atmo.speed_of_sound = 343.0;

    EXPECT_DOUBLE_EQ(atmo.density, 1.225);
    EXPECT_DOUBLE_EQ(atmo.dynamic_viscosity, 0.0000181);
    EXPECT_DOUBLE_EQ(atmo.kinematic_viscosity, 0.0000148);
    EXPECT_DOUBLE_EQ(atmo.speed_of_sound, 343.0);

    OC_InducedVelocities iv;
    for (int i = 0; i < 8; ++i) {
        iv.v_x[i] = static_cast<double>(i);
        iv.v_y[i] = static_cast<double>(i * 2);
        iv.v_z[i] = static_cast<double>(i * 3);
    }
    EXPECT_DOUBLE_EQ(iv.v_x[4], 4.0);
    EXPECT_DOUBLE_EQ(iv.v_y[4], 8.0);
    EXPECT_DOUBLE_EQ(iv.v_z[4], 12.0);
}