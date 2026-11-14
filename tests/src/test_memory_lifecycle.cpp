/* ------------------------------------------------------------------ */
/*  test_memory_lifecycle.cpp                                          */
/*                                                                    */
/*  Memory safety tests for basic lifecycle operations:                */
/*    - struct layout / alignment                                     */
/*    - Direction enum helpers                                        */
/*    - Aircraft create / destroy                                     */
/*    - RotorGeometry / WingGeometry create / destroy                  */
/*    - AircraftInputState create / destroy                            */
/*    - RotorInputState getter                                        */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cstddef>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v = {x, y, z};
    return v;
}

/* ================================================================== */
/*  TEST: verify struct layout                                        */
/* ================================================================== */
TEST(MemoryLifecycle, StructLayout) {
    EXPECT_EQ(sizeof(OC_Vec3), sizeof(double) * 3);
    EXPECT_EQ(sizeof(OC_Vec4), sizeof(double) * 4);
    EXPECT_EQ(sizeof(OC_Mat3),  sizeof(double) * 9);
    EXPECT_EQ(sizeof(OC_Mat4),  sizeof(double) * 16);
}

/* ================================================================== */
/*  TEST: direction helpers return valid values                       */
/* ================================================================== */
TEST(MemoryLifecycle, DirectionHelpers) {
    OC_Direction cw = oc_direction_clockwise();
    OC_Direction ccw = oc_direction_counter_clockwise();
    // Both should be small enum-like int values.
    EXPECT_TRUE(cw == 0 || cw == 1);
    EXPECT_TRUE(ccw == 0 || ccw == 1);
    // Clockwise and counter-clockwise should differ.
    EXPECT_NE(cw, ccw);
}

/* ================================================================== */
/*  TEST: Aircraft create / destroy                                   */
/* ================================================================== */
TEST(MemoryLifecycle, AircraftCreateDestroy) {
    OC_Aircraft* ac = oc_aircraft_create(0, 0);
    ASSERT_NE(ac, nullptr);
    oc_aircraft_destroy(ac);
}

/* ================================================================== */
/*  TEST: RotorGeometry create / destroy                               */
/* ================================================================== */
TEST(MemoryLifecycle, RotorGeometryCreateDestroy) {
    OC_RotorGeometry* r = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 2.0, 0.5);
    ASSERT_NE(r, nullptr);
    oc_rotor_geometry_destroy(r);
}

/* ================================================================== */
/*  TEST: WingGeometry create / destroy                                */
/* ================================================================== */
TEST(MemoryLifecycle, WingGeometryCreateDestroy) {
    OC_WingGeometry* w = oc_wing_geometry_create(1, make_vec3(0, 0, 0), 4.0);
    ASSERT_NE(w, nullptr);
    oc_wing_geometry_destroy(w);
}

/* ================================================================== */
/*  TEST: AircraftInputState create / destroy                          */
/* ================================================================== */
TEST(MemoryLifecycle, AircraftInputStateCreateDestroy) {
    size_t num_blades = 4;
    OC_AircraftInputState* s = oc_aircraft_input_state_create(1, &num_blades, 0);
    ASSERT_NE(s, nullptr);
    oc_aircraft_input_state_destroy(s);
}

/* ================================================================== */
/*  TEST: RotorInputState accessor                                     */
/* ================================================================== */
TEST(MemoryLifecycle, RotorInputStateAccessor) {
    size_t num_blades = 4;
    OC_AircraftInputState* s = oc_aircraft_input_state_create(1, &num_blades, 0);
    ASSERT_NE(s, nullptr);

    OC_RotorInputState* r = oc_aircraft_input_get_rotor_input(s, 0);
    ASSERT_NE(r, nullptr);

    oc_rotor_input_set_angular_velocity(r, 100.0);
    EXPECT_DOUBLE_EQ(oc_rotor_input_get_angular_velocity(r), 100.0);

    oc_aircraft_input_state_destroy(s);
}

/* ================================================================== */
/*  TEST: Multiple create / destroy cycles (no leaks)                  */
/* ================================================================== */
TEST(MemoryLifecycle, MultipleAircraftCycles) {
    for (int i = 0; i < 10; ++i) {
        OC_Aircraft* ac = oc_aircraft_create(0, 0);
        ASSERT_NE(ac, nullptr);
        oc_aircraft_destroy(ac);
    }
}

/* ================================================================== */
/*  TEST: RotorGeometry multiple cycles                                */
/* ================================================================== */
TEST(MemoryLifecycle, MultipleRotorCycles) {
    for (int i = 0; i < 10; ++i) {
        OC_RotorGeometry* r = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 2.0, 0.5);
        ASSERT_NE(r, nullptr);
        oc_rotor_geometry_destroy(r);
    }
}

/* ================================================================== */
/*  TEST: WingGeometry multiple cycles                                 */
/* ================================================================== */
TEST(MemoryLifecycle, MultipleWingCycles) {
    for (int i = 0; i < 10; ++i) {
        OC_WingGeometry* w = oc_wing_geometry_create(1, make_vec3(0, 0, 0), 4.0);
        ASSERT_NE(w, nullptr);
        oc_wing_geometry_destroy(w);
    }
}