/* ------------------------------------------------------------------ */
/*  test_api_wake.cpp                                                  */
/*                                                                    */
/*  Tests for Wake API:                                                */
/*    - Create/destroy lifecycle                                       */
/*    - Rotor wake access                                              */
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
TEST(Wake, CreateDestroy) {
    size_t shed_history = 5;
    size_t shed_release = 8;
    OC_Wake* wake = oc_wake_create(1, 2, 3, 4, &shed_history, &shed_release);
    // May be null if D side fails; if non-null destroy it
    if (wake != nullptr) {
        oc_wake_destroy(wake);
    }
    EXPECT_NE(wake, nullptr);
}


/* ================================================================== */
TEST(Wake, GetRotorWake) {
    size_t shed_history = 5;
    size_t shed_release = 8;
    OC_Wake* wake = oc_wake_create(1, 2, 3, 4, &shed_history, &shed_release);
    EXPECT_NE(wake, nullptr);

    if (wake != nullptr) {
        OC_RotorWake* rw = oc_wake_get_rotor_wake(wake, 0);
        // May be null or valid; the key is no crash
        (void)rw;
        oc_wake_destroy(wake);
    }
}


/* ================================================================== */
TEST(Wake, NullSafeDestroy) {
    oc_wake_destroy(nullptr); /* should not crash */
}

/* ================================================================== */
TEST(Wake, HistoryCreateDestroy) {
    size_t shed_history = 5;
    size_t shed_release = 8;
    OC_WakeHistory* hist = oc_wake_history_create(1, 2, 3, 4, 8, &shed_history, &shed_release, 0.5, 0);
    ASSERT_NE(hist, nullptr);

    oc_wake_history_destroy(hist);
}

/* ================================================================== */
TEST(Wake, HistoryPushBackAndGet) {
    size_t shed_history = 5;
    size_t shed_release = 8;
    OC_WakeHistory* hist = oc_wake_history_create(1, 2, 3, 4, 8, &shed_history, &shed_release, 0.5, 0);
    ASSERT_NE(hist, nullptr);

    oc_wake_history_push_back(hist);
    OC_Wake* w = oc_wake_history_get_wake(hist, 0);
    // May be null or valid; just verify no crash

    oc_wake_history_destroy(hist);
}

/* ================================================================== */
TEST(Wake, MultipleCycles) {
    size_t shed_history = 5;
    size_t shed_release = 8;
    for (int i = 0; i < 5; ++i) {
        OC_Wake* wake = oc_wake_create(1, 2, 3, 4, &shed_history, &shed_release);
        if (wake != nullptr) {
            oc_wake_destroy(wake);
        }
    }
    // At least the factory should be callable multiple times without crashing
}
