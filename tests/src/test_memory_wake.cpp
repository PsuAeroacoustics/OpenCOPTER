/* ------------------------------------------------------------------ */
/*  test_memory_wake.cpp                                               */
/*                                                                    */
/*  Memory safety tests for WakeHistory API.                           */
/*    oc_wake_create throws D exceptions when parameters don't meet     */
/*    internal constraints, so we skip it and test only WakeHistory     */
/*    which was observed to be stable from C.                          */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cstddef>

/* ================================================================== */
/*  TEST: WakeHistory create / destroy                                 */
/* ================================================================== */
TEST(MemoryWake, WakeHistoryCreateDestroy) {
    size_t shed_history = 10;
    size_t shed_release = 1;
    OC_WakeHistory* hist = oc_wake_history_create(
        1, 4, 256, 4, 47,
        &shed_history, &shed_release,
        6.5e-5, 0);
    ASSERT_NE(hist, nullptr);
    oc_wake_history_destroy(hist);
}

/* ================================================================== */
/*  TEST: WakeHistory push_back multiple times                         */
/* ================================================================== */
TEST(MemoryWake, WakeHistoryPushBack) {
    size_t shed_history = 10;
    size_t shed_release = 1;
    OC_WakeHistory* hist = oc_wake_history_create(
        1, 4, 256, 4, 47,
        &shed_history, &shed_release,
        6.5e-5, 0);
    ASSERT_NE(hist, nullptr);

    for (int i = 0; i < 10; ++i) {
        oc_wake_history_push_back(hist);
    }

    for (int i = 0; i < 10; ++i) {
        OC_Wake* w = oc_wake_history_get_wake(hist, i);
        (void)w;
    }

    oc_wake_history_destroy(hist);
}

/* ================================================================== */
/*  TEST: WakeHistory with different parameters                        */
/* ================================================================== */
TEST(MemoryWake, WakeHistoryDifferentParams) {
    size_t shed_history[2] = {50, 50};
    size_t shed_release[2] = {2, 2};
    OC_WakeHistory* hist = oc_wake_history_create(
        2, 8, 1024, 8, 94,
        &shed_history[0], &shed_release[0],
        6.5e-5, 0);
    ASSERT_NE(hist, nullptr);
    oc_wake_history_destroy(hist);
}

/* ================================================================== */
/*  TEST: WakeHistory hybrid mode                                      */
/* ================================================================== */
TEST(MemoryWake, WakeHistoryHybridMode) {
    size_t shed_history = 10;
    size_t shed_release = 1;
    OC_WakeHistory* hist = oc_wake_history_create(
        1, 4, 256, 4, 47,
        &shed_history, &shed_release,
        6.5e-5, 1);
    ASSERT_NE(hist, nullptr);

    oc_wake_history_push_back(hist);
    oc_wake_history_destroy(hist);
}

/* ================================================================== */
/*  TEST: Multiple WakeHistory create/destroy cycles                   */
/* ================================================================== */
TEST(MemoryWake, MultipleCycles) {
    for (int i = 0; i < 10; ++i) {
        size_t shed_history = 10;
        size_t shed_release = 1;
        OC_WakeHistory* hist = oc_wake_history_create(
            1, 4, 256, 4, 47,
            &shed_history, &shed_release,
            6.5e-5, 0);
        ASSERT_NE(hist, nullptr);
        oc_wake_history_destroy(hist);
    }
}