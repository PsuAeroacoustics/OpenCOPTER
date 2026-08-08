/* ------------------------------------------------------------------ */
/*  test_api_inflow.cpp                                                */
/*                                                                    */
/*  Tests for Inflow API:                                              */
/*    - NullInflow lifecycle (with valid rotor + rotor_input)          */
/*    - Null input handling                                            */
/*    - Frame access, wake_skew query                                  */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v; v.x = x; v.y = y; v.z = z;
    return v;
}

/* ================================================================== */
TEST(Inflow, NullInflowCreateDestroy) {
    size_t num_rotors = 1;
    OC_Aircraft* aircraft = oc_aircraft_create(num_rotors, 0);
    ASSERT_NE(aircraft, nullptr);

    OC_Frame* root_frame = oc_aircraft_get_root_frame(aircraft);
    ASSERT_NE(root_frame, nullptr);

    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_fixed_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        root_frame,
        "rotor_0_fixed",
        OC_CONNECTION_FRAME
    );

    ASSERT_NE(rotor_fixed_frame, nullptr);

    OC_Frame* root_fixed_children[] = { rotor_fixed_frame };
    oc_frame_set_children(root_frame, root_fixed_children, 1);

    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        rotor_fixed_frame,
        "rotor_0",
        OC_ROTOR_FRAME
    );

    ASSERT_NE(rotor_frame, nullptr);

    OC_Frame* rotor_fixed_children[] = { rotor_frame };
    oc_frame_set_children(rotor_fixed_frame, rotor_fixed_children, 1);

    // Build a minimal RotorGeometry + RotorInputState so we can create a NullInflow.
    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(2, origin, 1.0, 0.05);
    ASSERT_NE(rotor, nullptr);

    oc_rotor_geometry_set_frame(rotor, rotor_frame);

    // Create an AircraftInputState with 1 rotor, 2 blades to get a RotorInputState pointer
    size_t num_blades_arr[1] = {2};
    OC_AircraftInputState* ac_input = oc_aircraft_input_state_create(1, num_blades_arr, 0);
    ASSERT_NE(ac_input, nullptr);

    OC_RotorInputState* rotor_input = oc_aircraft_input_get_rotor_input(ac_input, 0);
    ASSERT_NE(rotor_input, nullptr);

    // Now create the NullInflow with valid rotor + rotor_input
    OC_Inflow* inflow = oc_null_inflow_create(rotor, rotor_input);
    // May be null if D side fails; if non-null destroy it
    if (inflow != nullptr) {
        oc_inflow_destroy(inflow);
    }
    EXPECT_NE(inflow, nullptr);

    oc_aircraft_input_state_destroy(ac_input);
    oc_rotor_geometry_destroy(rotor);
}


/* ================================================================== */
TEST(Inflow, NullInflowReturnsNullOnBadInput) {
    OC_Inflow* inflow = oc_null_inflow_create(nullptr, nullptr);
    EXPECT_EQ(inflow, nullptr);
}

/* ================================================================== */
TEST(Inflow, HuangPetersReturnsNullOnBadInput) {
    OC_Inflow* inflow = oc_huang_peters_create(10, 20, nullptr, nullptr, 0.01);
    EXPECT_EQ(inflow, nullptr);
}

/* ================================================================== */
TEST(Inflow, DestroyIdempotent) {
    size_t num_rotors = 1;
    OC_Aircraft* aircraft = oc_aircraft_create(num_rotors, 0);
    ASSERT_NE(aircraft, nullptr);

    OC_Frame* root_frame = oc_aircraft_get_root_frame(aircraft);
    ASSERT_NE(root_frame, nullptr);
    
    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_fixed_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        root_frame,
        "rotor_0_fixed",
        OC_CONNECTION_FRAME
    );

    ASSERT_NE(rotor_fixed_frame, nullptr);

    OC_Frame* root_fixed_children[] = { rotor_fixed_frame };
    oc_frame_set_children(root_frame, root_fixed_children, 1);

    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        rotor_fixed_frame,
        "rotor_0",
        OC_ROTOR_FRAME
    );

    ASSERT_NE(rotor_frame, nullptr);

    OC_Frame* rotor_fixed_children[] = { rotor_frame };
    oc_frame_set_children(rotor_fixed_frame, rotor_fixed_children, 1);

    // Build minimal objects for inflow creation
    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(2, origin, 1.0, 0.05);
    ASSERT_NE(rotor, nullptr);

    oc_rotor_geometry_set_frame(rotor, rotor_frame);
    
    size_t num_blades_arr[1] = {2};
    OC_AircraftInputState* ac_input = oc_aircraft_input_state_create(1, num_blades_arr, 0);
    ASSERT_NE(ac_input, nullptr);

    OC_RotorInputState* rotor_input = oc_aircraft_input_get_rotor_input(ac_input, 0);
    ASSERT_NE(rotor_input, nullptr);

    OC_Inflow* inflow = oc_null_inflow_create(rotor, rotor_input);
    if (inflow != nullptr) {
        // Calling destroy twice should be safe
        oc_inflow_destroy(inflow);
        oc_inflow_destroy(inflow);
    }

    oc_aircraft_input_state_destroy(ac_input);
    oc_rotor_geometry_destroy(rotor);
}


/* ================================================================== */
TEST(Inflow, NullInflowGetFrame) {
    size_t num_rotors = 1;
    OC_Aircraft* aircraft = oc_aircraft_create(num_rotors, 0);
    ASSERT_NE(aircraft, nullptr);

    OC_Frame* root_frame = oc_aircraft_get_root_frame(aircraft);
    ASSERT_NE(root_frame, nullptr);
    
    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_fixed_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        root_frame,
        "rotor_0_fixed",
        OC_CONNECTION_FRAME
    );

    ASSERT_NE(rotor_fixed_frame, nullptr);

    OC_Frame* root_fixed_children[] = { rotor_fixed_frame };
    oc_frame_set_children(root_frame, root_fixed_children, 1);

    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        rotor_fixed_frame,
        "rotor_0",
        OC_ROTOR_FRAME
    );

    ASSERT_NE(rotor_frame, nullptr);

    OC_Frame* rotor_fixed_children[] = { rotor_frame };
    oc_frame_set_children(rotor_fixed_frame, rotor_fixed_children, 1);

    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(2, origin, 1.0, 0.05);
    ASSERT_NE(rotor, nullptr);

    oc_rotor_geometry_set_frame(rotor, rotor_frame);

    size_t num_blades_arr[1] = {2};
    OC_AircraftInputState* ac_input = oc_aircraft_input_state_create(1, num_blades_arr, 0);
    ASSERT_NE(ac_input, nullptr);

    OC_RotorInputState* rotor_input = oc_aircraft_input_get_rotor_input(ac_input, 0);
    ASSERT_NE(rotor_input, nullptr);

    OC_Inflow* inflow = oc_null_inflow_create(rotor, rotor_input);
    if (inflow != nullptr) {
        OC_Frame* frame = oc_inflow_get_frame(inflow);
        // Frame may be null or valid; the key is no crash occurs
        (void)frame;
        oc_inflow_destroy(inflow);
    }
    EXPECT_NE(inflow, nullptr);

    oc_aircraft_input_state_destroy(ac_input);
    oc_rotor_geometry_destroy(rotor);
}


/* ================================================================== */
TEST(Inflow, NullInflowWakeSkew) {
    size_t num_rotors = 1;
    OC_Aircraft* aircraft = oc_aircraft_create(num_rotors, 0);
    ASSERT_NE(aircraft, nullptr);

    OC_Frame* root_frame = oc_aircraft_get_root_frame(aircraft);
    ASSERT_NE(root_frame, nullptr);
    
    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_fixed_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        root_frame,
        "rotor_0_fixed",
        OC_CONNECTION_FRAME
    );

    ASSERT_NE(rotor_fixed_frame, nullptr);

    OC_Frame* root_fixed_children[] = { rotor_fixed_frame };
    oc_frame_set_children(root_frame, root_fixed_children, 1);

    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        rotor_fixed_frame,
        "rotor_0",
        OC_ROTOR_FRAME
    );

    ASSERT_NE(rotor_frame, nullptr);

    OC_Frame* rotor_fixed_children[] = { rotor_frame };
    oc_frame_set_children(rotor_fixed_frame, rotor_fixed_children, 1);

    OC_Vec3 origin = make_vec3(0, 0, 0);
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(2, origin, 1.0, 0.05);
    ASSERT_NE(rotor, nullptr);

    oc_rotor_geometry_set_frame(rotor, rotor_frame);

    size_t num_blades_arr[1] = {2};
    OC_AircraftInputState* ac_input = oc_aircraft_input_state_create(1, num_blades_arr, 0);
    ASSERT_NE(ac_input, nullptr);

    OC_RotorInputState* rotor_input = oc_aircraft_input_get_rotor_input(ac_input, 0);
    ASSERT_NE(rotor_input, nullptr);

    OC_Inflow* inflow = oc_null_inflow_create(rotor, rotor_input);
    if (inflow != nullptr) {
        // NullInflow may return NaN for wake_skew since it has no inflow model
        double skew = oc_inflow_wake_skew(inflow);
        EXPECT_FALSE(std::isinf(skew));
        oc_inflow_destroy(inflow);
    }
    EXPECT_NE(inflow, nullptr);

    oc_aircraft_input_state_destroy(ac_input);
    oc_rotor_geometry_destroy(rotor);
}

/* ================================================================== */
TEST(Inflow, WingInflowNullInput) {
    OC_Inflow* inflow = oc_wing_inflow_create(nullptr, nullptr, nullptr);
    EXPECT_EQ(inflow, nullptr);
}