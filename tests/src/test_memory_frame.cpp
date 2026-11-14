/* ------------------------------------------------------------------ */
/*  test_memory_frame.cpp                                              */
/*                                                                    */
/*  Memory safety tests for the Frame API:                             */
/*    - Parent / child relationship lifecycle                          */
/*    - Destruction order (children first vs parent first)             */
/*    - oc_frame_set_children with various array sizes                */
/*    - Re-assignment of children                                      */
/*    - Rotate / translate / update do not corrupt memory             */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cstring>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v = {x, y, z};
    return v;
}

/* ================================================================== */
/*  TEST: Single frame create / destroy                               */
/* ================================================================== */
TEST(MemoryFrame, SingleFrameCreateDestroy) {
    OC_Frame* f = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                  make_vec3(0, 0, 0),
                                  nullptr, "test", OC_CONNECTION_FRAME);
    ASSERT_NE(f, nullptr);
    oc_frame_destroy(f);
}

/* ================================================================== */
/*  TEST: Parent-child – destroy child then parent                    */
/* ================================================================== */
TEST(MemoryFrame, DestroyChildThenParent) {
    OC_Frame* parent = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                       make_vec3(0, 0, 0),
                                       nullptr, "parent", OC_CONNECTION_FRAME);
    OC_Frame* child = oc_frame_create(make_vec3(1, 0, 0), 0.0,
                                      make_vec3(0, 0, 0),
                                      parent, "child", OC_BLADE_FRAME);
    ASSERT_NE(parent, nullptr);
    ASSERT_NE(child, nullptr);

    OC_Frame* children[1] = {child};
    oc_frame_set_children(parent, children, 1);

    // Destroy child first, then parent.
    oc_frame_destroy(child);
    oc_frame_destroy(parent);
}

/* ================================================================== */
/*  TEST: Parent-child – destroy parent then child                    */
/* ================================================================== */
TEST(MemoryFrame, DestroyParentThenChild) {
    OC_Frame* parent = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                       make_vec3(0, 0, 0),
                                       nullptr, "parent", OC_CONNECTION_FRAME);
    OC_Frame* child = oc_frame_create(make_vec3(1, 0, 0), 0.0,
                                      make_vec3(0, 0, 0),
                                      parent, "child", OC_BLADE_FRAME);
    ASSERT_NE(parent, nullptr);
    ASSERT_NE(child, nullptr);

    OC_Frame* children[1] = {child};
    oc_frame_set_children(parent, children, 1);

    // Destroy parent first (library should still manage internal refs).
    oc_frame_destroy(parent);
    oc_frame_destroy(child);
}

/* ================================================================== */
/*  TEST: Multiple children array                                     */
/* ================================================================== */
TEST(MemoryFrame, MultipleChildren) {
    OC_Frame* parent = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                       make_vec3(0, 0, 0),
                                       nullptr, "parent", OC_CONNECTION_FRAME);

    const int N = 8;
    OC_Frame* children[N] = {};
    for (int i = 0; i < N; ++i) {
        children[i] = oc_frame_create(make_vec3(0, 1, 0), 0.0,
                                      make_vec3(0, 0, 0),
                                      parent, "child", OC_BLADE_FRAME);
        ASSERT_NE(children[i], nullptr);
    }

    oc_frame_set_children(parent, children, N);

    // Cleanup: children first.
    for (int i = 0; i < N; ++i) {
        oc_frame_destroy(children[i]);
    }
    oc_frame_destroy(parent);
}

/* ================================================================== */
/*  TEST: Set children with zero count                                */
/* ================================================================== */
TEST(MemoryFrame, SetChildrenZeroCount) {
    OC_Frame* f = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                  make_vec3(0, 0, 0),
                                  nullptr, "test", OC_CONNECTION_FRAME);
    // Should not crash even with zero children.
    oc_frame_set_children(f, nullptr, 0);
    oc_frame_destroy(f);
}

/* ================================================================== */
/*  TEST: Rotate and translate do not leak                            */
/* ================================================================== */
TEST(MemoryFrame, RotateAndTranslate) {
    OC_Frame* f = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                  make_vec3(0, 0, 0),
                                  nullptr, "test", OC_CONNECTION_FRAME);

    // Apply several transforms.
    for (int i = 0; i < 10; ++i) {
        oc_frame_rotate(f, make_vec3(0, 0, 1), 0.1);
        oc_frame_translate(f, make_vec3(0.01, 0.02, 0.03));
    }

    oc_frame_destroy(f);
}

/* ================================================================== */
/*  TEST: set_rotation                                                */
/* ================================================================== */
TEST(MemoryFrame, SetRotation) {
    OC_Frame* f = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                  make_vec3(0, 0, 0),
                                  nullptr, "test", OC_CONNECTION_FRAME);

    oc_frame_set_rotation(f, make_vec3(1, 0, 0), 1.57);
    oc_frame_destroy(f);
}

/* ================================================================== */
/*  TEST: Frame update + matrix getters                               */
/* ================================================================== */
TEST(MemoryFrame, UpdateAndGetMatrix) {
    OC_Frame* f = oc_frame_create(make_vec3(0, 0, 1), 0.5,
                                  make_vec3(1, 2, 3),
                                  nullptr, "test", OC_BLADE_FRAME);

    OC_Mat4 identity = oc_mat4_identity();
    oc_frame_update(f, &identity);

    const OC_Mat4* local_mat = oc_frame_get_local_matrix(f);
    const OC_Mat4* global_mat = oc_frame_get_global_matrix(f);
    const OC_Mat4* inv_global = oc_frame_get_inverse_global_matrix(f);

    // Pointers should be non-null.
    EXPECT_NE(local_mat, nullptr);
    EXPECT_NE(global_mat, nullptr);
    EXPECT_NE(inv_global, nullptr);

    oc_frame_destroy(f);
}

/* ================================================================== */
/*  TEST: set_frame_type and set_name                                 */
/* ================================================================== */
TEST(MemoryFrame, SetFrameTypeAndName) {
    OC_Frame* f = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                  make_vec3(0, 0, 0),
                                  nullptr, "initial", OC_CONNECTION_FRAME);

    oc_frame_set_name(f, "renamed");
    oc_frame_set_frame_type(f, OC_ROTOR_FRAME);

    oc_frame_destroy(f);
}

/* ================================================================== */
/*  TEST: get_parent returns correct parent                           */
/* ================================================================== */
TEST(MemoryFrame, GetParent) {
    OC_Frame* parent = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                       make_vec3(0, 0, 0),
                                       nullptr, "parent", OC_CONNECTION_FRAME);
    OC_Frame* child = oc_frame_create(make_vec3(1, 0, 0), 0.0,
                                      make_vec3(0, 0, 0),
                                      parent, "child", OC_BLADE_FRAME);

    EXPECT_EQ(oc_frame_get_parent(child), parent);
    EXPECT_EQ(oc_frame_get_parent(parent), nullptr);

    oc_frame_destroy(child);
    oc_frame_destroy(parent);
}

/* ================================================================== */
/*  TEST: Reassign children                                           */
/* ================================================================== */
TEST(MemoryFrame, ReassignChildren) {
    OC_Frame* parent = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                       make_vec3(0, 0, 0),
                                       nullptr, "parent", OC_CONNECTION_FRAME);

    OC_Frame* old_child = oc_frame_create(make_vec3(1, 0, 0), 0.0,
                                          make_vec3(0, 0, 0),
                                          parent, "old", OC_BLADE_FRAME);
    OC_Frame* new_child = oc_frame_create(make_vec3(0, 1, 0), 0.0,
                                          make_vec3(0, 0, 0),
                                          parent, "new", OC_BLADE_FRAME);

    // First set of children.
    oc_frame_set_children(parent, &old_child, 1);
    // Replace with new child.
    oc_frame_set_children(parent, &new_child, 1);

    oc_frame_destroy(old_child);
    oc_frame_destroy(new_child);
    oc_frame_destroy(parent);
}