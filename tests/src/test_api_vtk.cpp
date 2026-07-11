/* ------------------------------------------------------------------ */
/*  test_api_vtk.cpp                                                   */
/*                                                                    */
/*  Tests for VTK API:                                                 */
/*    - Null-safe destroy for all VTK types                            */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"

/* ================================================================== */
TEST(VTK, DestroyRotorNull) {
    oc_vtk_rotor_destroy(nullptr); /* should not crash */
}

/* ================================================================== */
TEST(VTK, DestroyWingNull) {
    oc_vtk_wing_destroy(nullptr); /* should not crash */
}

/* ================================================================== */
TEST(VTK, DestroyWakeNull) {
    oc_vtk_wake_destroy(nullptr); /* should not crash */
}

/* ================================================================== */
TEST(VTK, DestroyWingWakeNull) {
    oc_vtk_wing_wake_destroy(nullptr); /* should not crash */
}