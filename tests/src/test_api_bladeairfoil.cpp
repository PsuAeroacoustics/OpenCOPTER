/* ------------------------------------------------------------------ */
/*  test_api_bladeairfoil.cpp                                          */
/*                                                                    */
/*  Tests for BladeAirfoil API:                                        */
/*    - CreateBasic lifecycle                                          */
/*    - Query methods (skipped: D interface vtable dispatch from C)    */
/*    - Null safety                                                    */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v; v.x = x; v.y = y; v.z = z;
    return v;
}

/* ================================================================== */
TEST(BladeAirfoil, CreateBasicDestroy) {
    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 6.28);
    ASSERT_NE(ba, nullptr);
    oc_blade_airfoil_destroy(ba);
}

/* ================================================================== */
TEST(BladeAirfoil, CreateBasicGetCl) {
    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 6.28);
    ASSERT_NE(ba, nullptr);

    // Query Cl at alpha=5 deg, mach=0.
    // D formula: Cl = 2*pi*(alpha - C_l_alpha_0) + C_l_alpha_0 = 2*pi*alpha + C_l_alpha_0
    // With C_l_alpha_0=6.28 (zero-lift AoA offset) and alpha=5deg: Cl ~ 2*pi*0.0873 + 6.28 ~ 6.83
    double cl = oc_blade_airfoil_get_Cl(ba, 0, 5.0 * 3.14159265358979 / 180.0, 0.0);
    EXPECT_NEAR(cl, 6.83, 0.2);

    oc_blade_airfoil_destroy(ba);
}

/* ================================================================== */
TEST(BladeAirfoil, CreateBasicZeroElements) {
    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(0, 6.28);
    oc_blade_airfoil_destroy(ba);
}

/* ================================================================== */
TEST(BladeAirfoil, GetCdReturnsZero) {
    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 6.28);
    ASSERT_NE(ba, nullptr);

    // Thin airfoil has zero drag
    double cd = oc_blade_airfoil_get_Cd(ba, 0, 5.0 * 3.14159265358979 / 180.0, 0.0);
    EXPECT_DOUBLE_EQ(cd, 0.0);

    oc_blade_airfoil_destroy(ba);
}

/* ================================================================== */
TEST(BladeAirfoil, LiftCurveSlope) {
    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 0);
    ASSERT_NE(ba, nullptr);

    double dCl_da = oc_blade_airfoil_lift_curve_slope(ba, 0);
    // Thin airfoil lift curve slope is always 2*pi (constant)
    EXPECT_NEAR(dCl_da, 2.0 * 3.14159265358979, 0.1);

    oc_blade_airfoil_destroy(ba);
}

/* ================================================================== */
TEST(BladeAirfoil, ZeroLiftAoa) {
    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 6.28);
    ASSERT_NE(ba, nullptr);

    // C_l_alpha_0 is the zero-lift AoA offset; with 6.28 passed, zero_lift_aoa returns 6.28
    double aoa = oc_blade_airfoil_zero_lift_aoa(ba, 0);
    EXPECT_NEAR(aoa, 6.28, 1e-6);

    oc_blade_airfoil_destroy(ba);
}

/* ================================================================== */
TEST(BladeAirfoil, FillCoefficients) {
    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 6.28);
    ASSERT_NE(ba, nullptr);

    double alphas[1] = {5.0 * 3.14159265358979 / 180.0};
    double machs[1]  = {0.0};
    double Cl_out[1] = {0.0};
    double Cd_out[1] = {0.0};

    oc_blade_airfoil_fill_coefficients(ba, 0, alphas, machs, Cl_out, Cd_out, 1);

    // Cl = 2*pi*(alpha - C_l_alpha_0) + C_l_alpha_0 = 2*pi*alpha + C_l_alpha_0
    // With C_l_alpha_0=6.28 and alpha=5deg: Cl ~ 6.83
    EXPECT_NEAR(Cl_out[0], 6.83, 0.2);
    // Cd is zero for thin airfoil
    EXPECT_DOUBLE_EQ(Cd_out[0], 0.0);

    oc_blade_airfoil_destroy(ba);
}

/* ================================================================== */
TEST(BladeAirfoil, FillLiftCurveSlopeBuffer) {
    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 6.28);
    ASSERT_NE(ba, nullptr);

    double buf[1] = {0.0};
    oc_blade_airfoil_fill_lift_curve_slope(ba, 0, buf, 1);

    // Should fill with ~6.28 (thin airfoil lift curve slope)
    EXPECT_NEAR(buf[0], 6.28, 0.1);

    oc_blade_airfoil_destroy(ba);
}

/* ================================================================== */
TEST(BladeAirfoil, NullSafe) {
    EXPECT_DOUBLE_EQ(oc_blade_airfoil_get_Cl(nullptr, 0, 5.0, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(oc_blade_airfoil_get_Cd(nullptr, 0, 5.0, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(oc_blade_airfoil_lift_curve_slope(nullptr, 0), 0.0);
    EXPECT_DOUBLE_EQ(oc_blade_airfoil_zero_lift_aoa(nullptr, 0), 0.0);
}

/* ================================================================== */
TEST(BladeAirfoil, MultipleCycles) {
    for (int i = 0; i < 10; ++i) {
        OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 6.28);
        ASSERT_NE(ba, nullptr);
        oc_blade_airfoil_destroy(ba);
    }
}