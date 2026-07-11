/* ------------------------------------------------------------------ */
/*  test_api_rotorgeo_advanced.cpp                                     */
/*                                                                    */
/*  Tests for RotorGeometry advanced API:                              */
/*    - set_solidity, set_frame, set_blades                            */
/*    - Composite lifecycle (aircraft -> rotors -> blades)             */
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
TEST(RotorGeoAdvanced, SetSolidity) {
    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(2, origin, 1.0, 0.3);
    ASSERT_NE(rotor, nullptr);

    oc_rotor_geometry_set_solidity(rotor, 0.5);

    oc_rotor_geometry_destroy(rotor);
}

/* ================================================================== */
TEST(RotorGeoAdvanced, SetFrame) {
    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(2, origin, 1.0, 0.3);
    ASSERT_NE(rotor, nullptr);

    OC_Vec3 zero = make_vec3(0, 0, 0);
    OC_Vec3 zaxis = make_vec3(0, 0, 1);
    OC_Frame* f = oc_frame_create(zaxis, 0.0, zero, nullptr, "rotor_frame", 0);
    ASSERT_NE(f, nullptr);

    oc_rotor_geometry_set_frame(rotor, f);

    oc_rotor_geometry_destroy(rotor);
    oc_frame_destroy(f);
}

/* ================================================================== */
TEST(RotorGeoAdvanced, SetBlades) {
    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(2, origin, 1.0, 0.3);
    ASSERT_NE(rotor, nullptr);

    OC_BladeGeometry* b0 = oc_blade_geometry_create(8, 0.0, 0.3, nullptr, 0.5);
    OC_BladeGeometry* b1 = oc_blade_geometry_create(8, 1.57, 0.3, nullptr, 0.5);
    ASSERT_NE(b0, nullptr);
    ASSERT_NE(b1, nullptr);

    OC_BladeGeometry* blades[2] = {b0, b1};
    oc_rotor_geometry_set_blades(rotor, blades, 2);

    oc_rotor_geometry_destroy(rotor);
    oc_blade_geometry_destroy(b0);
    oc_blade_geometry_destroy(b1);
}

/* ================================================================== */
TEST(RotorGeoAdvanced, CompositeLifecycle) {
    // Aircraft(1 rotor, 0 wings) -> Rotor -> 2 Blades
    OC_Aircraft* ac = oc_aircraft_create(1, 0);
    ASSERT_NE(ac, nullptr);

    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(2, origin, 1.0, 0.3);
    ASSERT_NE(rotor, nullptr);

    OC_BladeGeometry* b0 = oc_blade_geometry_create(8, 0.0, 0.3, nullptr, 0.5);
    OC_BladeGeometry* b1 = oc_blade_geometry_create(8, 1.57, 0.3, nullptr, 0.5);
    ASSERT_NE(b0, nullptr);
    ASSERT_NE(b1, nullptr);

    OC_BladeGeometry* blades[2] = {b0, b1};
    oc_rotor_geometry_set_blades(rotor, blades, 2);

    OC_RotorGeometry* rotors[1] = {rotor};
    oc_aircraft_set_rotors(ac, rotors, 1);

    // Destroy in reverse order
    oc_aircraft_destroy(ac);
    oc_rotor_geometry_destroy(rotor);
    oc_blade_geometry_destroy(b0);
    oc_blade_geometry_destroy(b1);
}

/* ================================================================== */
TEST(RotorGeoAdvanced, MultipleRotorsOnAircraft) {
    OC_Aircraft* ac = oc_aircraft_create(2, 0);
    ASSERT_NE(ac, nullptr);

    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* r0 = oc_rotor_geometry_create(2, origin, 1.0, 0.3);
    OC_RotorGeometry* r1 = oc_rotor_geometry_create(2, origin, 1.0, 0.3);
    ASSERT_NE(r0, nullptr);
    ASSERT_NE(r1, nullptr);

    OC_BladeAirfoil* ba = oc_blade_airfoil_create_basic(8, 0);

    OC_BladeGeometry* b0 = oc_blade_geometry_create(8, 0.0, 0.3, ba, 0.5);
    OC_BladeGeometry* b1 = oc_blade_geometry_create(8, 1.57, 0.3, ba, 0.5);
    OC_BladeGeometry* b2 = oc_blade_geometry_create(8, 0.0, 0.3, ba, 0.5);
    OC_BladeGeometry* b3 = oc_blade_geometry_create(8, 1.57, 0.3, ba, 0.5);

    OC_BladeGeometry* blades_r0[2] = {b0, b1};
    OC_BladeGeometry* blades_r1[2] = {b2, b3};
    oc_rotor_geometry_set_blades(r0, blades_r0, 2);
    oc_rotor_geometry_set_blades(r1, blades_r1, 2);

    OC_RotorGeometry* rotors[2] = {r0, r1};
    oc_aircraft_set_rotors(ac, rotors, 2);

    oc_aircraft_destroy(ac);
    oc_rotor_geometry_destroy(r0);
    oc_rotor_geometry_destroy(r1);
    oc_blade_geometry_destroy(b0);
    oc_blade_geometry_destroy(b1);
    oc_blade_geometry_destroy(b2);
    oc_blade_geometry_destroy(b3);
}

/* ================================================================== */
TEST(RotorGeoAdvanced, NullSafeSetBlades) {
    OC_BladeGeometry* b = oc_blade_geometry_create(8, 0.0, 0.3, nullptr, 0.5);
    OC_BladeGeometry* blades[1] = {b};
    oc_rotor_geometry_set_blades(nullptr, blades, 1); /* should not crash */
    oc_blade_geometry_destroy(b);
}