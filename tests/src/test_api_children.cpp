/* ------------------------------------------------------------------ */
/*  test_api_children.cpp                                              */
/*                                                                    */
/*  Tests for the new children retrieval API:                          */
/*    - oc_frame_get_children / oc_frame_get_children_count (C API)    */
/*    - Frame::children() (C++ wrapper)                                */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include "opencopter.hpp"
#include <vector>
#include <functional>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v = {x, y, z};
    return v;
}

/* ================================================================== */
/*  C API Tests                                                        */
/* ================================================================== */

/* TEST: Get children on frame with no children */
TEST(ChildrenC, GetChildrenEmpty) {
    OC_Frame* f = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                  make_vec3(0, 0, 0),
                                  nullptr, "empty", OC_CONNECTION_FRAME);
    ASSERT_NE(f, nullptr);

    EXPECT_EQ(oc_frame_get_children_count(f), 0u);
    EXPECT_EQ(oc_frame_get_children(f), nullptr);

    oc_frame_destroy(f);
}

/* TEST: Get single child */
TEST(ChildrenC, GetChildrenSingle) {
    OC_Frame* parent = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                       make_vec3(0, 0, 0),
                                       nullptr, "parent", OC_AIRCRAFT_FRAME);
    OC_Frame* child = oc_frame_create(make_vec3(1, 0, 0), 0.0,
                                      make_vec3(0, 0, 0),
                                      parent, "child", OC_ROTOR_FRAME);
    ASSERT_NE(parent, nullptr);
    ASSERT_NE(child, nullptr);

    OC_Frame* children[1] = {child};
    oc_frame_set_children(parent, children, 1);

    EXPECT_EQ(oc_frame_get_children_count(parent), 1u);
    OC_Frame** got = oc_frame_get_children(parent);
    ASSERT_NE(got, nullptr);
    EXPECT_EQ(got[0], child);

    oc_frame_destroy(child);
    oc_frame_destroy(parent);
}

/* TEST: Get multiple children */
TEST(ChildrenC, GetChildrenMultiple) {
    OC_Frame* parent = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                       make_vec3(0, 0, 0),
                                       nullptr, "parent", OC_AIRCRAFT_FRAME);

    const int N = 5;
    OC_Frame* children[N] = {};
    for (int i = 0; i < N; ++i) {
        children[i] = oc_frame_create(make_vec3(0, 1, 0), 0.0,
                                      make_vec3(0, 0, 0),
                                      parent, "child", OC_CONNECTION_FRAME);
        ASSERT_NE(children[i], nullptr);
    }

    oc_frame_set_children(parent, children, N);

    EXPECT_EQ(oc_frame_get_children_count(parent), static_cast<size_t>(N));
    OC_Frame** got = oc_frame_get_children(parent);
    ASSERT_NE(got, nullptr);
    for (int i = 0; i < N; ++i) {
        EXPECT_EQ(got[i], children[i]);
    }

    for (int i = 0; i < N; ++i) {
        oc_frame_destroy(children[i]);
    }
    oc_frame_destroy(parent);
}

/* TEST: Get children after reassignment */
TEST(ChildrenC, GetChildrenAfterReassign) {
    OC_Frame* parent = oc_frame_create(make_vec3(0, 0, 1), 0.0,
                                       make_vec3(0, 0, 0),
                                       nullptr, "parent", OC_AIRCRAFT_FRAME);

    OC_Frame* old_child = oc_frame_create(make_vec3(1, 0, 0), 0.0,
                                          make_vec3(0, 0, 0),
                                          parent, "old", OC_BLADE_FRAME);
    OC_Frame* new_child = oc_frame_create(make_vec3(0, 1, 0), 0.0,
                                          make_vec3(0, 0, 0),
                                          parent, "new", OC_BLADE_FRAME);

    // First set.
    oc_frame_set_children(parent, &old_child, 1);
    EXPECT_EQ(oc_frame_get_children_count(parent), 1u);

    // Replace.
    oc_frame_set_children(parent, &new_child, 1);
    EXPECT_EQ(oc_frame_get_children_count(parent), 1u);
    OC_Frame** got = oc_frame_get_children(parent);
    ASSERT_NE(got, nullptr);
    EXPECT_EQ(got[0], new_child);

    oc_frame_destroy(old_child);
    oc_frame_destroy(new_child);
    oc_frame_destroy(parent);
}

/* TEST: Null parent safety */
TEST(ChildrenC, GetChildrenNullParent) {
    EXPECT_EQ(oc_frame_get_children_count(nullptr), 0u);
    EXPECT_EQ(oc_frame_get_children(nullptr), nullptr);
}

/* ================================================================== */
/*  C++ API Tests                                                      */
/* ================================================================== */

/* TEST: Cpp children empty */
TEST(ChildrenCpp, ChildrenEmpty) {
    opencopter::Frame frame(opencopter::Vec3{0, 0, 1}, 0.0,
                            opencopter::Vec3{0, 0, 0},
                            nullptr, "empty", opencopter::FrameType::Connection);

    auto children = frame.children();
    EXPECT_TRUE(children.empty());
}

/* TEST: Cpp single child */
TEST(ChildrenCpp, ChildrenSingle) {
    opencopter::Frame parent(opencopter::Vec3{0, 0, 1}, 0.0,
                             opencopter::Vec3{0, 0, 0},
                             nullptr, "parent", opencopter::FrameType::Aircraft);

    opencopter::Frame child(opencopter::Vec3{1, 0, 0}, 0.0,
                            opencopter::Vec3{0, 0, 0},
                            &parent, "child", opencopter::FrameType::Rotor);

    const opencopter::Frame* children_arr[1] = {&child};
    parent.set_children({children_arr, 1});

    auto got = parent.children();
    EXPECT_EQ(got.size(), 1u);
    EXPECT_TRUE(got[0]); // non-null wrapper
}

/* TEST: Cpp multiple children */
TEST(ChildrenCpp, ChildrenMultiple) {
    opencopter::Frame parent(opencopter::Vec3{0, 0, 1}, 0.0,
                             opencopter::Vec3{0, 0, 0},
                             nullptr, "parent", opencopter::FrameType::Aircraft);

    // Use stack-allocated frames to avoid vector reallocation issues
    opencopter::Frame c0(opencopter::Vec3{0, 1, 0}, 0.0,
                         opencopter::Vec3{0, 0, 0},
                         &parent, "c0", opencopter::FrameType::Connection);
    opencopter::Frame c1(opencopter::Vec3{0, 1, 0}, 0.0,
                         opencopter::Vec3{0, 0, 0},
                         &parent, "c1", opencopter::FrameType::Connection);
    opencopter::Frame c2(opencopter::Vec3{0, 1, 0}, 0.0,
                         opencopter::Vec3{0, 0, 0},
                         &parent, "c2", opencopter::FrameType::Connection);
    opencopter::Frame c3(opencopter::Vec3{0, 1, 0}, 0.0,
                         opencopter::Vec3{0, 0, 0},
                         &parent, "c3", opencopter::FrameType::Connection);
    opencopter::Frame c4(opencopter::Vec3{0, 1, 0}, 0.0,
                         opencopter::Vec3{0, 0, 0},
                         &parent, "c4", opencopter::FrameType::Connection);

    const int N = 5;
    const opencopter::Frame* children_arr[N] = {&c0, &c1, &c2, &c3, &c4};
    parent.set_children({children_arr, N});

    auto got = parent.children();
    EXPECT_EQ(got.size(), static_cast<size_t>(N));
}

/* TEST: Cpp recursive tree walk */
TEST(ChildrenCpp, RecursiveWalk) {
    // Build a small tree:
    //   root
    //   ├── childA
    //   │   ├── grandA1
    //   │   └── grandA2
    //   └── childB
    //       ├── grandB1
    //       └── grandB2

    opencopter::Frame root(opencopter::Vec3{0, 0, 1}, 0.0,
                           opencopter::Vec3{0, 0, 0},
                           nullptr, "root", opencopter::FrameType::Aircraft);

    opencopter::Frame childA(opencopter::Vec3{1, 0, 0}, 0.0,
                             opencopter::Vec3{0, 0, 0},
                             &root, "childA", opencopter::FrameType::Connection);
    opencopter::Frame childB(opencopter::Vec3{0, 1, 0}, 0.0,
                             opencopter::Vec3{0, 0, 0},
                             &root, "childB", opencopter::FrameType::Connection);

    const opencopter::Frame* rootChildren[2] = {&childA, &childB};
    root.set_children({rootChildren, 2});

    EXPECT_EQ(root.children().size(), 2u);

    // Grandchildren for A
    opencopter::Frame grandA1(opencopter::Vec3{0, 0, 1}, 0.0,
                              opencopter::Vec3{0, 0, 0},
                              &childA, "grandA1", opencopter::FrameType::Rotor);
    opencopter::Frame grandA2(opencopter::Vec3{0, 0, 1}, 0.0,
                              opencopter::Vec3{0, 0, 0},
                              &childA, "grandA2", opencopter::FrameType::Rotor);

    const opencopter::Frame* aChildren[2] = {&grandA1, &grandA2};
    childA.set_children({aChildren, 2});

    // Grandchildren for B
    opencopter::Frame grandB1(opencopter::Vec3{0, 0, 1}, 0.0,
                              opencopter::Vec3{0, 0, 0},
                              &childB, "grandB1", opencopter::FrameType::Rotor);
    opencopter::Frame grandB2(opencopter::Vec3{0, 0, 1}, 0.0,
                              opencopter::Vec3{0, 0, 0},
                              &childB, "grandB2", opencopter::FrameType::Rotor);

    const opencopter::Frame* bChildren[2] = {&grandB1, &grandB2};
    childB.set_children({bChildren, 2});

    // Verify tree structure via recursive walk
    size_t total_nodes_walked = 0;
    std::function<void(const opencopter::Frame&)> walk =
        [&](const opencopter::Frame& frame) {
            total_nodes_walked++;
            for (const auto& child : frame.children()) {
                if (child) walk(child);
            }
        };

    walk(root);
    // root + 2 children + 4 grandchildren = 7
    EXPECT_EQ(total_nodes_walked, 7u);

    // Verify each level
    EXPECT_EQ(root.children().size(), 2u);
    auto root_children = root.children();
    EXPECT_EQ(root_children[0].children().size(), 2u);
    EXPECT_EQ(root_children[1].children().size(), 2u);
}

/* TEST: Cpp children are non-owning wrappers */
TEST(ChildrenCpp, ChildrenNonOwning) {
    opencopter::Frame parent(opencopter::Vec3{0, 0, 1}, 0.0,
                             opencopter::Vec3{0, 0, 0},
                             nullptr, "parent", opencopter::FrameType::Aircraft);

    opencopter::Frame child(opencopter::Vec3{1, 0, 0}, 0.0,
                            opencopter::Vec3{0, 0, 0},
                            &parent, "child", opencopter::FrameType::Rotor);

    const opencopter::Frame* children_arr[1] = {&child};
    parent.set_children({children_arr, 1});

    // Get children - these are non-owning wrappers
    auto got = parent.children();
    EXPECT_EQ(got.size(), 1u);
    EXPECT_TRUE(got[0]); // wrapper wraps the raw pointer, not owned

    // Calling children() again should return same count (idempotent)
    auto got2 = parent.children();
    EXPECT_EQ(got2.size(), 1u);
}
