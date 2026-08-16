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
TEST(BladeGeo, AzimuthOffsetRoundTrip) {
    OC_BladeAirfoil* af = make_test_airfoil(); ASSERT_NE(af, nullptr);

    // Constructor value is retrievable
    const double init = 0.785;
    OC_BladeGeometry* bg = oc_blade_geometry_create(8, init, 0.3, af, 0.5);
    ASSERT_NE(bg, nullptr);
    EXPECT_DOUBLE_EQ(oc_blade_geometry_get_azimuth_offset(bg), init);

    // Setter updates the value
    const double newval = 1.5708;
    oc_blade_geometry_set_azimuth_offset(bg, newval);
    EXPECT_DOUBLE_EQ(oc_blade_geometry_get_azimuth_offset(bg), newval);

    // Negative value round-trips
    oc_blade_geometry_set_azimuth_offset(bg, -0.5);
    EXPECT_DOUBLE_EQ(oc_blade_geometry_get_azimuth_offset(bg), -0.5);

    oc_blade_geometry_destroy(bg);
    oc_blade_airfoil_destroy(af);
}

/* ================================================================== */
TEST(BladeGeo, AzimuthOffsetNullSafe) {
    // Null getter returns 0.0, null setter does not crash
    EXPECT_DOUBLE_EQ(oc_blade_geometry_get_azimuth_offset(nullptr), 0.0);
    oc_blade_geometry_set_azimuth_offset(nullptr, 1.0); /* should not crash */
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

/* ================================================================== */
/*  Tests for new Step 2.4 C API functions                             */
/* ================================================================== */

/* Helper: create a rotor with 1 blade */
static OC_RotorGeometry* make_test_rotor(size_t num_blades) {
    return oc_rotor_geometry_create(num_blades, make_vec3(0,0,0), 0.5, 0.05);
}

TEST(AircraftAccessors, GetNumRotors) {
    OC_Aircraft* ac = oc_aircraft_create(2, 0);
    ASSERT_NE(ac, nullptr);
    size_t n = oc_aircraft_get_num_rotors(ac);
    EXPECT_EQ(n, 2u);
    oc_aircraft_destroy(ac);
}

TEST(AircraftAccessors, GetNumRotorsZero) {
    OC_Aircraft* ac = oc_aircraft_create(0, 0);
    ASSERT_NE(ac, nullptr);
    size_t n = oc_aircraft_get_num_rotors(ac);
    EXPECT_EQ(n, 0u);
    oc_aircraft_destroy(ac);
}

TEST(AircraftAccessors, GetRotorValidIndex) {
    OC_Aircraft* ac = oc_aircraft_create(2, 0);
    ASSERT_NE(ac, nullptr);
    OC_RotorGeometry* r0 = oc_aircraft_get_rotor(ac, 0);
    OC_RotorGeometry* r1 = oc_aircraft_get_rotor(ac, 1);
    EXPECT_NE(r0, nullptr);
    EXPECT_NE(r1, nullptr);
    EXPECT_NE(r0, r1); /* different rotors */
    oc_aircraft_destroy(ac);
}

TEST(AircraftAccessors, GetRotorOutOfBounds) {
    OC_Aircraft* ac = oc_aircraft_create(1, 0);
    ASSERT_NE(ac, nullptr);
    OC_RotorGeometry* r = oc_aircraft_get_rotor(ac, 5); /* out of bounds */
    EXPECT_EQ(r, nullptr);
    oc_aircraft_destroy(ac);
}

TEST(AircraftAccessors, GetRotorNullAircraft) {
    OC_RotorGeometry* r = oc_aircraft_get_rotor(nullptr, 0);
    EXPECT_EQ(r, nullptr);
}

TEST(RotorGeo, GetFrame) {
    OC_RotorGeometry* rotor = make_test_rotor(1);
    OC_Frame* frame = oc_frame_create(make_vec3(0,0,1), 0.0, make_vec3(0,0,0), nullptr, "test", 2);
    ASSERT_NE(frame, nullptr);
    oc_rotor_geometry_set_frame(rotor, frame);
    OC_Frame* got = oc_rotor_geometry_get_frame(rotor);
    EXPECT_NE(got, nullptr);
    EXPECT_EQ(got, frame); /* same pointer */
    oc_rotor_geometry_destroy(rotor);
    oc_frame_destroy(frame);
}

TEST(RotorGeo, GetFrameNull) {
    OC_Frame* got = oc_rotor_geometry_get_frame(nullptr);
    EXPECT_EQ(got, nullptr);
}
