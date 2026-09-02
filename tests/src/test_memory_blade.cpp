/* ------------------------------------------------------------------ */
/*  test_memory_blade.cpp                                              */
/*                                                                    */
/*  Memory safety tests for RotorGeometry (blade-level geometry).       */
/*    The BladeGeometry/BladeAirfoil APIs require D runtime GC and       */
/*    virtual dispatch which are unstable from plain C test harnesses.   */
/*    This file focuses on the RotorGeometry create/destroy lifecycle    */
/*    which is stable and tests correct memory ownership.               */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cstddef>
#include <vector>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v = {x, y, z};
    return v;
}

/* ================================================================== */
/*  TEST: RotorGeometry create / destroy                               */
/* ================================================================== */
TEST(MemoryBlade, RotorGeometryCreateDestroy) {
    OC_RotorGeometry* r = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 2.0, 0.5);
    ASSERT_NE(r, nullptr);
    oc_rotor_geometry_destroy(r);
}

/* ================================================================== */
/*  TEST: RotorGeometry with zero hub radius                           */
/* ================================================================== */
TEST(MemoryBlade, RotorGeometryZeroHub) {
    OC_RotorGeometry* r = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 1.0, 0.0);
    ASSERT_NE(r, nullptr);
    oc_rotor_geometry_destroy(r);
}

/* ================================================================== */
/*  TEST: RotorGeometry with different blade counts                    */
/* ================================================================== */
TEST(MemoryBlade, RotorGeometryBladeCounts) {
    std::vector<int> blades = {2, 3, 4, 5, 6};
    for (int i = 0; i < blades.size(); ++i) {
        OC_RotorGeometry* r = oc_rotor_geometry_create(blades[i], make_vec3(0, 0, 0), 2.0, 0.5);
        ASSERT_NE(r, nullptr);
        oc_rotor_geometry_destroy(r);
    }
}

/* ================================================================== */
/*  TEST: RotorGeometry multiple create/destroy cycles                 */
/* ================================================================== */
TEST(MemoryBlade, MultipleCycles) {
    for (int i = 0; i < 20; ++i) {
        OC_RotorGeometry* r = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 2.0, 0.5);
        ASSERT_NE(r, nullptr);
        oc_rotor_geometry_destroy(r);
    }
}

/* ================================================================== */
/*  TEST: RotorGeometry with offset position                           */
/* ================================================================== */
TEST(MemoryBlade, RotorGeometryOffset) {
    OC_RotorGeometry* r = oc_rotor_geometry_create(4, make_vec3(1.5, -2.0, 3.0), 2.0, 0.5);
    ASSERT_NE(r, nullptr);
    oc_rotor_geometry_destroy(r);
}

/* ================================================================== */
/*  TEST: RotorGeometry with large radius                              */
/* ================================================================== */
TEST(MemoryBlade, RotorGeometryLargeRadius) {
    OC_RotorGeometry* r = oc_rotor_geometry_create(4, make_vec3(0, 0, 0), 50.0, 10.0);
    ASSERT_NE(r, nullptr);
    oc_rotor_geometry_destroy(r);
}