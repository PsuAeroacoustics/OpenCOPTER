/* ------------------------------------------------------------------ */
/*  test_api_bladegroup.cpp                                            */
/*                                                                    */
/*  Tests for BladeGeometry API:                                       */
/*    - Create/destroy lifecycle (with null airfoil)                   */
/*    - All setter methods                                             */
/*    - compute_vectors, frame get/set                                 */
/*    - Null safety                                                    */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>
#include <cstring>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v; v.x = x; v.y = y; v.z = z;
    return v;
}

/* Helper: create a valid BladeAirfoil for BladeGeometry tests */
static OC_BladeAirfoil* make_test_airfoil() {
    return oc_blade_airfoil_create_basic(8, 6.28);
}

/* ================================================================== */
TEST(BladeGeo, CreateDestroy) {
    OC_BladeAirfoil* af = make_test_airfoil();
    ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);
    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, CreateZeroElements) {
    // D side allows creating BladeGeometry with 0 elements - it succeeds.
    OC_BladeAirfoil* af = make_test_airfoil();
    ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(0, 0.0, 0.3, af, 0.5);
    EXPECT_NE(bg, nullptr);
    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetTwist) {
    OC_BladeAirfoil* af = make_test_airfoil();
    ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double twist[8];
    for (int i = 0; i < 8; ++i) twist[i] = i * 0.1;
    oc_blade_geometry_set_twist(bg, twist, 8);

    oc_blade_geometry_destroy(bg);
}

/* ================================================================== */
TEST(BladeGeo, SetChord) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double chord[8];
    for (int i = 0; i < 8; ++i) chord[i] = 0.1 + i * 0.01;
    oc_blade_geometry_set_chord(bg, chord, 8);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetRadius) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double radius[8];
    for (int i = 0; i < 8; ++i) radius[i] = 0.1 + i * 0.1;
    oc_blade_geometry_set_radius(bg, radius, 8);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetClAlpha) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double cl_alpha[8];
    for (int i = 0; i < 8; ++i) cl_alpha[i] = 6.28;
    oc_blade_geometry_set_C_l_alpha(bg, cl_alpha, 8);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetAlpha0) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double alpha0[8];
    for (int i = 0; i < 8; ++i) alpha0[i] = 0.0;
    oc_blade_geometry_set_alpha_0(bg, alpha0, 8);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetSweep) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double sweep[8];
    for (int i = 0; i < 8; ++i) sweep[i] = 0.0;
    oc_blade_geometry_set_sweep(bg, sweep, 8);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetXi) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double xi[8];
    for (int i = 0; i < 8; ++i) xi[i] = 0.25;
    oc_blade_geometry_set_xi(bg, xi, 8);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetThickness) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double thick[8];
    for (int i = 0; i < 8; ++i) thick[i] = 0.12;
    oc_blade_geometry_set_thickness(bg, thick, 8);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetXiP) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    double xip[8];
    for (int i = 0; i < 8; ++i) xip[i] = 0.5;
    oc_blade_geometry_set_xi_p(bg, xip, 8);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, SetBladeLength) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    oc_blade_geometry_set_blade_length(bg, 1.0);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, ComputeVectors) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    // Should not crash
    oc_blade_geometry_compute_vectors(bg);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, FrameGetSet) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);

    OC_Vec3 zero = make_vec3(0, 0, 0);
    OC_Vec3 zaxis = make_vec3(0, 0, 1);
    OC_Frame* f = oc_frame_create(zaxis, 0.0, zero, nullptr, "blade_frame", 0);
    ASSERT_NE(f, nullptr);

    oc_blade_geometry_set_frame(bg, f);
    OC_Frame* got = oc_blade_geometry_get_frame(bg);
    EXPECT_NE(got, nullptr);

    oc_blade_geometry_destroy(bg);
    oc_frame_destroy(f);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, NullSafeDestroy) {
    oc_blade_geometry_destroy(nullptr); /* should not crash */
}

/* ================================================================== */
TEST(BladeGeo, MultipleCycles) {
    for (int i = 0; i < 10; ++i) {
        OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);
        OC_BladeGeometry* bg = oc_blade_geometry_create(8, 0.0, 0.3, af, 0.5);
        ASSERT_NE(bg, nullptr);
        oc_blade_geometry_destroy(bg);
        oc_blade_airfoil_destroy(af);
    }
}
