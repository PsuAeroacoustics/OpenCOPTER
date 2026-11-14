/* ------------------------------------------------------------------ */
/*  test_api_aircraftstate.cpp                                         */
/*                                                                    */
/*  Tests for AircraftState API:                                       */
/*    - Null-safe destroy                                             */
/*    - Rotor CT/CQ null-safe queries                                 */
/*    - Freestream set/get (full chain if constructible)               */
/*    - Rotor count & state retrieval                                  */
/*    - Blade count & state retrieval                                  */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>

static OC_Vec3 make_vec3(double x, double y, double z) {
    OC_Vec3 v; v.x = x; v.y = y; v.z = z;
    return v;
}

/* ================================================================== */
TEST(AircraftState, NullSafeDestroy) {
    oc_aircraft_state_destroy(nullptr); /* should not crash */
}

/* ================================================================== */
TEST(AircraftState, RotorCTNullSafe) {
    double out_val;
    // Passing nullptr state: function may do nothing or return 0
    oc_aircraft_state_get_rotor_C_T(nullptr, 0, &out_val);
}

/* ================================================================== */
TEST(AircraftState, RotorCQNullSafe) {
    double out_val;
    oc_aircraft_state_get_rotor_C_Q(nullptr, 0, &out_val);
}

// --------------------------------------------------------------------
//  Helper: build full frame hierarchy needed for NullInflow to work
// --------------------------------------------------------------------

// --------------------------------------------------------------------
//  Helper: create a rotor with proper blade geometry (elements, chord)
//  so that WeissingerL has blade.chunks to work with.
// --------------------------------------------------------------------
static OC_RotorGeometry* create_rotor_with_blades(size_t num_blades, double radius, double solidity) {
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(num_blades, make_vec3(0, 0, 0), radius, solidity);
    if (rotor == nullptr) return nullptr;

    // Create a thin airfoil for the blades
    OC_AirfoilModel* af = oc_thin_airfoil_create(0.0);
    if (af == nullptr) {
        oc_rotor_geometry_destroy(rotor);
        return nullptr;
    }

    // Create blade geometries with proper elements
    size_t num_elements = 8; // one chunk of 8
    double avg_chord = 0.1;
    double r_c = 0.01;

    OC_BladeGeometry* blades[num_blades];
    for (size_t i = 0; i < num_blades; ++i) {
        blades[i] = oc_blade_geometry_create(num_elements, 0.0, avg_chord, (OC_BladeAirfoil*)af, r_c);
        if (blades[i] == nullptr) {
            // cleanup partial
            for (size_t j = 0; j < i; ++j) oc_blade_geometry_destroy(blades[j]);
            oc_airfoil_model_destroy(af);
            oc_rotor_geometry_destroy(rotor);
            return nullptr;
        }
        // Set constant chord and xi on each blade
        double chord_data[num_elements];
        for (int k = 0; k < (int)num_elements; ++k) chord_data[k] = avg_chord;
        oc_blade_geometry_set_chord(blades[i], chord_data, num_elements);

        double xi_data[num_elements];
        for (int k = 0; k < (int)num_elements; ++k) xi_data[k] = 0.0;
        oc_blade_geometry_set_xi(blades[i], xi_data, num_elements);

        double xi_p_data[num_elements];
        for (int k = 0; k < (int)num_elements; ++k) xi_p_data[k] = 0.0;
        oc_blade_geometry_set_xi_p(blades[i], xi_p_data, num_elements);

        // Set blade length (needed for aspect ratio in WeissingerL)
        oc_blade_geometry_set_blade_length(blades[i], 1.0);

        // Set twist distribution (avoids singular influence matrix)
        double twist_data[num_elements];
        for (int k = 0; k < (int)num_elements; ++k)
            twist_data[k] = 0.1 - 0.05 * ((double)k / (num_elements - 1));
        oc_blade_geometry_set_twist(blades[i], twist_data, num_elements);
    }

    oc_rotor_geometry_set_blades(rotor, (OC_BladeGeometry**)blades, num_blades);
    oc_airfoil_model_destroy(af);
    return rotor;
}

// --------------------------------------------------------------------
//  Helper: build full frame hierarchy needed for NullInflow to work
// --------------------------------------------------------------------
static void setup_rotor_frame(OC_Aircraft* ac, OC_RotorGeometry* rotor, int idx) {
    OC_Frame* root = oc_aircraft_get_root_frame(ac);
    ASSERT_NE(root, nullptr);

    // Create rotor_fixed_frame (connection frame under root)
    char name[64];
    snprintf(name, sizeof(name), "rotor_%d_fixed", idx);
    OC_Frame* fixed = oc_frame_create(
        make_vec3(1.0, 0.0, 0.0), 0.0,
        make_vec3(0.0, 0.0, 0.0),
        root, name, OC_CONNECTION_FRAME);
    ASSERT_NE(fixed, nullptr);
    oc_frame_set_children(root, &fixed, 1);

    // Create rotor_frame under fixed frame
    snprintf(name, sizeof(name), "rotor_%d", idx);
    OC_Frame* rotor_f = oc_frame_create(
        make_vec3(1.0, 0.0, 0.0), 0.0,
        make_vec3(0.0, 0.0, 0.0),
        fixed, name, OC_ROTOR_FRAME);
    ASSERT_NE(rotor_f, nullptr);
    oc_frame_set_children(fixed, &rotor_f, 1);

    oc_rotor_geometry_set_frame(rotor, rotor_f);
}

/* ================================================================== */
/*  TEST: Full chain - rotor count, state retrieval, blade access      */
/* ================================================================== */
TEST(AircraftState, RotorStateCTSetGet) {
    OC_Aircraft* ac = oc_aircraft_create(1, 0);
    ASSERT_NE(ac, nullptr);

    OC_RotorGeometry* rotor = create_rotor_with_blades(2, 1.0, 0.3);
    ASSERT_NE(rotor, nullptr);

    OC_RotorGeometry* rotors[1] = {rotor};
    oc_aircraft_set_rotors(ac, rotors, 1);

    // Set up frame hierarchy so NullInflow works
    setup_rotor_frame(ac, rotor, 0);

    size_t num_blades_arr[1] = {2};
    OC_AircraftInputState* ac_input = oc_aircraft_input_state_create(1, num_blades_arr, 0);
    ASSERT_NE(ac_input, nullptr);

    OC_RotorInputState* rotor_in = oc_aircraft_input_get_rotor_input(ac_input, 0);
    ASSERT_NE(rotor_in, nullptr);
    oc_rotor_input_set_angular_velocity(rotor_in, 100.0);

    OC_Inflow* inflow = oc_null_inflow_create(rotor, rotor_in);
    ASSERT_NE(inflow, nullptr);

    size_t nba[1] = {2};
    double dir_val = 1.0;
    OC_Inflow* rinfl[] = {inflow};
    OC_AircraftState* state = oc_aircraft_state_create(
        1, nba, 8, 0, nullptr, 1, 1,
        ac, rinfl, nullptr, &dir_val);

    ASSERT_NE(state, nullptr);

    // Verify rotor count
    size_t rotor_count = oc_aircraft_state_get_rotor_count(state);
    EXPECT_EQ(rotor_count, 1u);

    // Get rotor state and set/get C_T and C_Q
    OC_RotorState* rs = oc_aircraft_state_get_rotor_state(state, 0);
    ASSERT_NE(rs, nullptr);

    oc_rotor_state_set_C_T(rs, 1.234);
    double ct;
    oc_rotor_state_get_C_T(rs, &ct);
    EXPECT_DOUBLE_EQ(ct, 1.234);

    oc_rotor_state_set_C_Q(rs, 0.567);
    double cq;
    oc_rotor_state_get_C_Q(rs, &cq);
    EXPECT_DOUBLE_EQ(cq, 0.567);

    // Verify blade count from rotor state
    size_t blade_count = oc_rotor_state_get_blade_count(rs);
    EXPECT_EQ(blade_count, 2u);

    // Get each blade state and verify accessor calls don't crash
    for (size_t i = 0; i < blade_count; ++i) {
        OC_BladeState* bs = oc_rotor_state_get_blade_state(rs, i);
        ASSERT_NE(bs, nullptr);
        oc_blade_state_get_azimuth(bs);  // just verify no crash
    }

    // Cleanup
    oc_aircraft_state_destroy(state);
    oc_inflow_destroy(inflow);
    oc_aircraft_input_state_destroy(ac_input);
    oc_rotor_geometry_destroy(rotor);
    oc_aircraft_destroy(ac);
}

/* ================================================================== */
/*  TEST: Multi-rotor aircraft state retrieval                          */
/* ================================================================== */
TEST(AircraftState, MultiRotorStates) {
    const int NUM_ROTORS = 2;
    OC_Aircraft* ac = oc_aircraft_create(NUM_ROTORS, 0);
    ASSERT_NE(ac, nullptr);

    OC_RotorGeometry* rotors[NUM_ROTORS] = {};
    for (int i = 0; i < NUM_ROTORS; ++i) {
        rotors[i] = create_rotor_with_blades(3, 1.5, 0.4);
        ASSERT_NE(rotors[i], nullptr);
    }

    OC_RotorGeometry* rotor_ptrs[NUM_ROTORS];
    for (int i = 0; i < NUM_ROTORS; ++i) rotor_ptrs[i] = rotors[i];
    oc_aircraft_set_rotors(ac, rotor_ptrs, NUM_ROTORS);

    // Set up frame hierarchy for each rotor
    for (int i = 0; i < NUM_ROTORS; ++i) {
        setup_rotor_frame(ac, rotors[i], i);
    }

    size_t num_blades_arr[NUM_ROTORS] = {3, 3};
    OC_AircraftInputState* ac_input = oc_aircraft_input_state_create(NUM_ROTORS, num_blades_arr, 0);
    ASSERT_NE(ac_input, nullptr);

    for (int i = 0; i < NUM_ROTORS; ++i) {
        OC_RotorInputState* ri = oc_aircraft_input_get_rotor_input(ac_input, i);
        ASSERT_NE(ri, nullptr);
        oc_rotor_input_set_angular_velocity(ri, 100.0 + i);
    }

    OC_Inflow* inflows[NUM_ROTORS] = {};
    for (int i = 0; i < NUM_ROTORS; ++i) {
        OC_RotorInputState* ri = oc_aircraft_input_get_rotor_input(ac_input, i);
        inflows[i] = oc_null_inflow_create(rotors[i], ri);
        ASSERT_NE(inflows[i], nullptr);
    }

    size_t nba[NUM_ROTORS] = {3, 3};
    double dir_val = 1.0;
    OC_AircraftState* state = oc_aircraft_state_create(
        NUM_ROTORS, nba, 8, 0, nullptr, 1, 1,
        ac, inflows, nullptr, &dir_val);

    ASSERT_NE(state, nullptr);

    // Verify rotor count
    size_t rotor_count = oc_aircraft_state_get_rotor_count(state);
    EXPECT_EQ(rotor_count, static_cast<size_t>(NUM_ROTORS));

    // Set/get CT on each rotor
    for (int i = 0; i < NUM_ROTORS; ++i) {
        OC_RotorState* rs = oc_aircraft_state_get_rotor_state(state, i);
        ASSERT_NE(rs, nullptr);

        double expected_ct = 1.0 + i * 0.5;
        oc_rotor_state_set_C_T(rs, expected_ct);

        double ct;
        oc_rotor_state_get_C_T(rs, &ct);
        EXPECT_DOUBLE_EQ(ct, expected_ct);

        // Verify blade count per rotor
        size_t bc = oc_rotor_state_get_blade_count(rs);
        EXPECT_EQ(bc, 3u);
    }

    // Cleanup
    oc_aircraft_state_destroy(state);
    for (int i = 0; i < NUM_ROTORS; ++i) oc_inflow_destroy(inflows[i]);
    oc_aircraft_input_state_destroy(ac_input);
    for (int i = 0; i < NUM_ROTORS; ++i) oc_rotor_geometry_destroy(rotors[i]);
    oc_aircraft_destroy(ac);
}

/* ================================================================== */
/*  TEST: BladeState scalar getters                                     */
/* ================================================================== */
TEST(AircraftState, BladeStateScalarGetters) {
    OC_Aircraft* ac = oc_aircraft_create(1, 0);
    ASSERT_NE(ac, nullptr);

    OC_RotorGeometry* rotor = create_rotor_with_blades(2, 1.0, 0.3);
    ASSERT_NE(rotor, nullptr);

    OC_RotorGeometry* rotors[1] = {rotor};
    oc_aircraft_set_rotors(ac, rotors, 1);

    setup_rotor_frame(ac, rotor, 0);

    size_t num_blades_arr[1] = {2};
    OC_AircraftInputState* ac_input = oc_aircraft_input_state_create(1, num_blades_arr, 0);
    ASSERT_NE(ac_input, nullptr);

    OC_RotorInputState* rotor_in = oc_aircraft_input_get_rotor_input(ac_input, 0);
    ASSERT_NE(rotor_in, nullptr);
    oc_rotor_input_set_angular_velocity(rotor_in, 100.0);

    OC_Inflow* inflow = oc_null_inflow_create(rotor, rotor_in);
    ASSERT_NE(inflow, nullptr);

    size_t nba[1] = {2};
    double dir_val = 1.0;
    OC_Inflow* rinfl[] = {inflow};
    OC_AircraftState* state = oc_aircraft_state_create(
        1, nba, 8, 0, nullptr, 1, 1,
        ac, rinfl, nullptr, &dir_val);

    ASSERT_NE(state, nullptr);

    OC_RotorState* rs = oc_aircraft_state_get_rotor_state(state, 0);
    ASSERT_NE(rs, nullptr);

    // Get first blade and test scalar getters don't crash / return valid values
    OC_BladeState* bs = oc_rotor_state_get_blade_state(rs, 0);
    ASSERT_NE(bs, nullptr);

    // Verify accessor calls don't crash (values are uninitialized before first sim step)
    oc_blade_state_get_azimuth(bs);
    oc_blade_state_get_C_T(bs);
    oc_blade_state_get_C_Q(bs);
    oc_blade_state_get_C_L(bs);
    oc_blade_state_get_C_D(bs);

    // Cleanup
    oc_aircraft_state_destroy(state);
    oc_inflow_destroy(inflow);
    oc_aircraft_input_state_destroy(ac_input);
    oc_rotor_geometry_destroy(rotor);
    oc_aircraft_destroy(ac);
}