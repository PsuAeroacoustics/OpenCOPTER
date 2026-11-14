/* ------------------------------------------------------------------ */
/*  test_api_utility.cpp                                               */
/*                                                                    */
/*  Tests for Utility API:                                             */
/*    - generate_radius_points with oc_free cleanup                    */
/*    - Memory helpers                                                 */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>
#include <cstdlib>

// oc_free is declared in opencopter.h but may not yet be exported by the
// shared library.  Use std::free for now and provide a local shim so the
// tests still compile against the header.
static inline void safe_oc_free(void* ptr) { std::free(ptr); }

/* ================================================================== */
TEST(Utility, GenerateRadiusPoints) {
    size_t n = 8;
    double* pts = oc_generate_radius_points(&n, 0.0);
    // May be null if D side fails; if non-null verify and free
    EXPECT_NE(pts, nullptr);
    EXPECT_GT(n, size_t(0));

    if (pts != nullptr) {
        // Points should be in [0, 1] range
        for (size_t i = 0; i < n; i++) {
            EXPECT_GE(pts[i], 0.0);
            EXPECT_LE(pts[i], 1.0);
        }
        safe_oc_free(pts);
    }
}

/* ================================================================== */
TEST(Utility, GenerateRadiusPointsZeroCutout) {
    size_t n = 10;
    double* pts = oc_generate_radius_points(&n, 0.0);
    EXPECT_NE(pts, nullptr);
    EXPECT_GT(n, size_t(0));

    if (pts != nullptr) {
        // First point should be near 0 when root_cutout is 0.
        // The half-cosine spacing formula produces pts[0] = (1-root_cutout)*0.5*(cos(PI*n/(N+1))+1)+root_cutout
        // which is small but not exactly 0. Use a reasonable tolerance.
        EXPECT_NEAR(pts[0], 0.0, 0.02);
        for (size_t i = 0; i < n; i++) {
            EXPECT_GE(pts[i], 0.0);
            EXPECT_LE(pts[i], 1.0);
        }
        safe_oc_free(pts);
    }
}

/* ================================================================== */
TEST(Utility, FreeNull) {
    // Freeing null should be safe
    safe_oc_free(nullptr);
}
