module opencopter.cbindings;

/**
 * C API for OpenCOPTER.
 *
 * This module provides a stable C interface to the OpenCOPTER simulation
 * library. The interface uses opaque pointers for all complex types,
 * allowing C++ (or any other language) consumers to interact with the
 * D library without needing D runtime dependencies or tooling.
 *
 * A C++ wrapper layer (in opencopter.h + cpp_wrapper.cpp) provides RAII
 * smart pointers and idiomatic C++ classes on top of this C API.
 *
 * Memory management:
 *   Since D is garbage-collected, objects allocated with `new` are pinned
 *   in the GC using GC.addRoot so they persist across FFI boundaries.
 *   Destroy functions use GC.removeRoot to allow reclamation.
 *
 * Author: OpenCOPTER Team
 * License: MIT
 */

import opencopter.aircraft;
import opencopter.aircraft : get_wing_state_array, get_wing_state_matrix, WingPartGeometryChunk, WingPartCtrlPointChunk, WingPartStateChunk, WingPartCtrlPointStateChunk;
import opencopter.airfoilmodels;
import opencopter.atmosphere;
import opencopter.bladeelement;
import opencopter.bwi;
import opencopter.config;
import opencopter.inflow;
import opencopter.math;
import opencopter.memory;
import opencopter.vortexlattice;
import opencopter.wake;

static import opencopter.aircraft.geometry;
static import opencopter.wake;

import core.memory : GC;
import core.stdc.stdlib : malloc;

import std.algorithm;
import std.array;
import std.conv : to;
import std.exception : enforce;
import std.math : abs, fmod, PI;
import std.string : fromStringz;
import std.traits : isBasicType;
import std.string : toStringz;

import numd.linearalgebra.matrix;

// ========================================================================
//  OC_ Opaque Type Definitions — D-side FFI declarations
//  These struct types match the C header definitions so they can be
//  cast to/from the actual D objects across the FFI boundary.
// ========================================================================

struct OC_Vec3 {
    double x;
    double y;
    double z;
}

struct OC_Vec4 {
    double x;
    double y;
    double z;
    double w;
}

struct OC_Mat3 {
    double[9] data;
}

struct OC_Mat4 {
    double[16] data;
}

enum OC_Direction : int {
    clockwise = 0,
    counter_clockwise = 1
}

enum OC_FrameType : int {
    aircraft = 0,
    connection = 1,
    rotor = 2,
    blade = 3,
    wing = 4
}

struct OC_Atmosphere {
    double density;
    double dynamic_viscosity;
    double kinematic_viscosity;
    double speed_of_sound;
}

struct OC_InducedVelocities {
    double[8] v_x;
    double[8] v_y;
    double[8] v_z;
}

// Forward declarations for opaque pointer types (empty structs).
// These allow the D compiler to recognize the type names.
struct OC_Aircraft               {}
struct OC_AircraftInputState     {}
struct OC_AircraftState          {}
struct OC_BladeGeometry          {}
struct OC_BladeState             {}
struct OC_BladeAirfoil           {}
struct OC_AirfoilModel           {}
struct OC_Frame                  {}
struct OC_Inflow                 {}
struct OC_HuangPeters            {}
struct OC_NullInflow             {}
struct OC_WingInflow             {}
struct OC_RotorGeometry          {}
struct OC_RotorInputState        {}
struct OC_RotorState             {}
struct OC_RotorWake              {}
struct OC_VortexFilament         {}
struct OC_Wake                   {}
struct OC_WakeHistory            {}
struct OC_WingGeometry           {}
struct OC_WingInputState         {}
struct OC_WingLiftSurf           {}
struct OC_WingPartGeometry       {}
struct OC_WingState              {}
struct OC_VtkRotor               {}
struct OC_VtkWing                {}
struct OC_VtkWake                {}
struct OC_VtkWingWake            {}

// ========================================================================
//  Helper functions: convert between OC_ value types and D numd types
// ========================================================================

/** Convert OC_Vec3 to Vec3 (Matrix!(3,1)) */
Vec3 oc_vec3_to_vec3(const OC_Vec3 v) {
    Vec3 _out;
    _out[0] = v.x;
    _out[1] = v.y;
    _out[2] = v.z;
    return _out;
}

/** Convert Vec3 to OC_Vec3 */
OC_Vec3 vec3_to_oc_vec3(const Vec3 v) {
    OC_Vec3 _out;
    _out.x = v[0];
    _out.y = v[1];
    _out.z = v[2];
    return _out;
}

/** Convert OC_Vec4 to Vec4 (Matrix!(4,1)) */
Vec4 oc_vec4_to_vec4(const OC_Vec4 v) {
    Vec4 _out;
    _out[0] = v.x;
    _out[1] = v.y;
    _out[2] = v.z;
    _out[3] = v.w;
    return _out;
}

/** Convert Vec4 to OC_Vec4 */
OC_Vec4 vec4_to_oc_vec4(const Vec4 v) {
    OC_Vec4 _out;
    _out.x = v[0];
    _out.y = v[1];
    _out.z = v[2];
    _out.w = v[3];
    return _out;
}

/** Convert OC_Mat4 to Mat4 (Matrix!(4,4)) */
Mat4 oc_mat4_to_mat4(const OC_Mat4 m) {
    Mat4 _out;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            _out[r, c] = m.data[r * 4 + c];
        }
    }
    return _out;
}

/** Convert Mat4 to OC_Mat4 */
OC_Mat4 mat4_to_oc_mat4(const Mat4 m) {
    OC_Mat4 _out;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            _out.data[r * 4 + c] = m[r, c];
        }
    }
    return _out;
}

/** Convert OC_Mat3 to Mat3-like (embedded in Mat4 for D side) */
Mat4 oc_mat3_to_mat3(const OC_Mat3 m) {
    Mat4 _out;
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            _out[r, c] = m.data[r * 3 + c];
        }
    }
    return _out;
}

// ========================================================================
//  Configuration
// ========================================================================

extern(C) size_t oc_chunk_size() {
    return chunk_size;
}

// ========================================================================
//  Direction helpers
// ========================================================================

extern(C) OC_Direction oc_direction_clockwise() {
    return OC_Direction.clockwise;
}

extern(C) OC_Direction oc_direction_counter_clockwise() {
    return OC_Direction.counter_clockwise;
}

// ========================================================================
//  Matrix helpers
// ========================================================================

extern(C) OC_Mat3 oc_mat3_identity() {
    Mat4 ident4 = Mat4.identity();
    OC_Mat3 mat3;
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            mat3.data[r * 3 + c] = ident4[r, c];
        }
    }
    return mat3;
}

extern(C) OC_Mat4 oc_mat4_identity() {
    auto ident = Mat4.identity();
    OC_Mat4 mat4;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            mat4.data[r * 4 + c] = ident[r, c];
        }
    }
    return mat4;
}

// ========================================================================
//  Frame API — all Vec3/Mat4 parameters converted to OC_ types
// ========================================================================

extern(C) OC_Frame* oc_frame_create(OC_Vec3 axis, double angle, OC_Vec3 translation,
                                      OC_Frame* parent, const(char)* name, int frame_type) {
    string str_name = (name !is null) ? fromStringz(name).idup : "";
    auto ftype = cast(FrameType)frame_type;
    auto d_axis = oc_vec3_to_vec3(axis);
    auto d_trans = oc_vec3_to_vec3(translation);
    auto d_parent = cast(Frame*)parent;
    auto frame = new Frame(d_axis, angle, d_trans, d_parent, str_name, ftype);
    GC.addRoot(frame);
    return cast(OC_Frame*)frame;
}

extern(C) void oc_frame_destroy(OC_Frame* frame) {
    if (frame !is null) {
        GC.removeRoot(cast(Frame*)frame);
    }
}

extern(C) void oc_frame_rotate(OC_Frame* frame, OC_Vec3 axis, double angle) {
    auto f = cast(Frame*)frame;
    if (f !is null) f.rotate(oc_vec3_to_vec3(axis), angle);
}

extern(C) void oc_frame_set_rotation(OC_Frame* frame, OC_Vec3 axis, double angle) {
    auto f = cast(Frame*)frame;
    if (f !is null) f.set_rotation(oc_vec3_to_vec3(axis), angle);
}

extern(C) void oc_frame_translate(OC_Frame* frame, OC_Vec3 translation) {
    auto f = cast(Frame*)frame;
    if (f !is null) f.translate(oc_vec3_to_vec3(translation));
}

extern(C) void oc_frame_update(OC_Frame* frame, const OC_Mat4* parent_global_mat) {
    auto f = cast(Frame*)frame;
    if (f !is null && parent_global_mat !is null) {
        Mat4 parent_mat = oc_mat4_to_mat4(*parent_global_mat);
        f.update(parent_mat);
    }
}

extern(C) const(OC_Mat4*) oc_frame_get_local_matrix(const OC_Frame* frame) {
    auto f = cast(const(Frame)*)frame;
    if (f !is null) {
        return cast(const OC_Mat4*)&(*f).local_matrix;
    }
    return null;
}

extern(C) const(OC_Mat4*) oc_frame_get_global_matrix(const OC_Frame* frame) {
    auto f = cast(const(Frame)*)frame;
    if (f !is null) {
        return cast(const OC_Mat4*)&(*f).global_matrix;
    }
    return null;
}

extern(C) const(OC_Mat4*) oc_frame_get_inverse_global_matrix(const OC_Frame* frame) {
    auto f = cast(const(Frame)*)frame;
    if (f !is null) {
        return cast(const OC_Mat4*)&(*f).inverse_global_matrix;
    }
    return null;
}

// ========================================================================
//  Aircraft API
// ========================================================================

extern(C) OC_Aircraft* oc_aircraft_create(size_t num_rotors, size_t num_wings) {
    auto ac = new Aircraft(num_rotors, num_wings);
    GC.addRoot(ac);
    return cast(OC_Aircraft*)ac;
}

extern(C) void oc_aircraft_destroy(OC_Aircraft* ac) {
    if (ac !is null) {
        GC.removeRoot(cast(Aircraft*)ac);
    }
}

// ========================================================================
//  RotorGeometry API
// ========================================================================

extern(C) OC_RotorGeometry* oc_rotor_geometry_create(size_t num_blades, OC_Vec3 origin,
                                                      double radius, double solidity) {
    auto rotor = new RotorGeometry(num_blades, oc_vec3_to_vec3(origin), radius, solidity);
    GC.addRoot(rotor);
    return cast(OC_RotorGeometry*)rotor;
}

extern(C) void oc_rotor_geometry_destroy(OC_RotorGeometry* geom) {
    if (geom !is null) {
        GC.removeRoot(cast(RotorGeometry*)geom);
    }
}

extern(C) void oc_rotor_geometry_set_solidity(OC_RotorGeometry* rotor, double solidity) {
    auto r = cast(RotorGeometry*)rotor;
    if (r !is null) r.solidity = solidity;
}

// ========================================================================
//  BladeGeometry API
// ========================================================================

extern(C) OC_BladeGeometry* oc_blade_geometry_create(size_t num_elements, double azimuth_offset,
                                                      double average_chord, OC_BladeAirfoil* airfoil,
                                                      double r_c) {
    // BladeAirfoil is a D class (reference type). Cast OC_BladeAirfoil* directly to BladeAirfoil.
    auto blade_af = cast(BladeAirfoil)airfoil;
    if(blade_af !is null) {
        // BladeGeometryT is an extern(C++) struct allocated on the heap via new.
        auto blade = new BladeGeometry(num_elements, azimuth_offset, average_chord, blade_af, r_c);
        GC.addRoot(blade);
        return cast(OC_BladeGeometry*)blade;
    }
    else {
        return null;
    }
}

extern(C) void oc_blade_geometry_destroy(OC_BladeGeometry* geom) {
    if (geom !is null) {
        GC.removeRoot(cast(BladeGeometry*)geom);
    }
}

extern(C) void oc_blade_geometry_set_twist(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"twist"(data[0..len]);
}

extern(C) void oc_blade_geometry_set_chord(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"chord"(data[0..len]);
}

extern(C) void oc_blade_geometry_set_radius(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"r"(data[0..len]);
}

extern(C) void oc_blade_geometry_set_C_l_alpha(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"C_l_alpha"(data[0..len]);
}

extern(C) void oc_blade_geometry_set_alpha_0(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"alpha_0"(data[0..len]);
}

extern(C) void oc_blade_geometry_set_sweep(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"sweep"(data[0..len]);
}

extern(C) void oc_blade_geometry_set_xi(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"xi"(data[0..len]);
}

extern(C) void oc_blade_geometry_set_thickness(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"thickness"(data[0..len]);
}

extern(C) void oc_blade_geometry_set_xi_p(OC_BladeGeometry* geom, double* data, size_t len) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null && data !is null) g.set_geometry_array!"xi_p"(data[0..len]);
}

extern(C) void oc_blade_geometry_compute_vectors(OC_BladeGeometry* geom) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null) compute_blade_vectors(*g);
}

// ========================================================================
//  WingGeometry API
// ========================================================================

extern(C) OC_WingGeometry* oc_wing_geometry_create(size_t num_parts, OC_Vec3 origin, double wing_span) {
    auto wing = new WingGeometry(num_parts, oc_vec3_to_vec3(origin), wing_span);
    GC.addRoot(wing);
    return cast(OC_WingGeometry*)wing;
}

extern(C) void oc_wing_geometry_destroy(OC_WingGeometry* geom) {
    if (geom !is null) GC.removeRoot(cast(WingGeometry*)geom);
}

extern(C) void oc_wing_part_geometry_set_chord(OC_WingPartGeometry* wg, double* data, size_t len) {
    auto g = cast(WingPartGeometry*)wg;
    if (g !is null && data !is null) g.set_geometry_array!"chord"(data[0..len]);
}

extern(C) void oc_wing_part_geometry_set_twist(OC_WingPartGeometry* wg, double* data, size_t len) {
    auto g = cast(WingPartGeometry*)wg;
    if (g !is null && data !is null) g.set_geometry_array!"twist"(data[0..len]);
}

extern(C) void oc_wing_part_geometry_set_sweep(OC_WingPartGeometry* wg, double* data, size_t len) {
    auto g = cast(WingPartGeometry*)wg;
    if (g !is null && data !is null) g.set_geometry_array!"sweep"(data[0..len]);
}

extern(C) void oc_wing_part_geometry_set_y_span(OC_WingPartGeometry* wg, double* data, size_t len) {
    auto g = cast(WingPartGeometry*)wg;
    if (g !is null && data !is null) g.set_geometry_array!"y_span"(data[0..len]);
}

extern(C) void oc_wing_geometry_set_ctrl_points(OC_WingGeometry* wing, size_t spanwise_nodes,
                                                  size_t chordwise_nodes, double camber) {
    auto w = cast(WingGeometry*)wing;
    if (w !is null) set_wing_ctrl_pt_geometry(w, spanwise_nodes, chordwise_nodes, camber);
}

// ========================================================================
//  AircraftInputState API
// ========================================================================

extern(C) OC_AircraftInputState* oc_aircraft_input_state_create(size_t num_rotors, size_t* num_blades, size_t num_wings) {
    auto input = new AircraftInputState(num_rotors, num_blades[0..num_rotors], num_wings);
    GC.addRoot(input);
    return cast(OC_AircraftInputState*)input;
}

extern(C) void oc_aircraft_input_state_destroy(OC_AircraftInputState* input) {
    if (input !is null) GC.removeRoot(cast(AircraftInputState*)input);
}

// ========================================================================
//  AircraftState API
// ========================================================================

extern(C) OC_AircraftState* oc_aircraft_state_create(
    size_t num_rotors,
    const size_t* num_blades,
    size_t num_elements,
    size_t num_wings,
    const size_t* num_wing_parts,
    size_t num_span_nodes,
    size_t num_chord_nodes,
    OC_Aircraft* aircraft,
    OC_Inflow** rotor_inflows,
    OC_Inflow** wing_inflows,
    const double* direction
)
{
    auto oc_ac = cast(Aircraft*)aircraft;

    auto oc_rotor_inflows = new InflowT!(ArrayContainer.none)[num_rotors];
    auto oc_wing_inflows  = new InflowT!(ArrayContainer.none)[num_wings];

    if (rotor_inflows !is null) {
        for (size_t i = 0; i < num_rotors; i++) {
            auto inflow_val = cast(Inflow)rotor_inflows[i];
            if (inflow_val !is null) oc_rotor_inflows[i] = inflow_val;
        }
    }
    if (wing_inflows !is null) {
        for (size_t i = 0; i < num_wings; i++) {
            auto inflow_val = cast(Inflow)wing_inflows[i];
            if (inflow_val !is null) oc_wing_inflows[i] = inflow_val;
        }
    }

    auto dir = (direction !is null) ? direction[0..num_rotors] : new double[num_rotors];

    auto state = new AircraftState(
        num_rotors,
        (num_blades !is null) ? num_blades[0..num_rotors] : null,
        num_elements,
        num_wings,
        (num_wing_parts !is null) ? num_wing_parts[0..num_wings] : null,
        num_span_nodes,
        num_chord_nodes,
        *oc_ac,
        oc_rotor_inflows,
        oc_wing_inflows,
        dir
    );

    GC.addRoot(state);
    return cast(OC_AircraftState*)state;
}

extern(C) void oc_aircraft_state_destroy(OC_AircraftState* state) {
    if (state !is null) GC.removeRoot(cast(AircraftState*)state);
}

// ========================================================================
//  State accessors — RotorState / AircraftState
// ========================================================================

extern(C) void oc_rotor_state_get_C_T(const OC_RotorState* state, double* result_out) {
    auto s = cast(const(RotorState)*)state;
    if (s !is null && result_out !is null) *result_out = (*s).C_T;
}
extern(C) void oc_rotor_state_set_C_T(OC_RotorState* state, double C_T) {
    auto s = cast(RotorState*)state;
    if (s !is null) (*s).C_T = C_T;
}
extern(C) void oc_rotor_state_get_C_Q(const OC_RotorState* state, double* result_out) {
    auto s = cast(const(RotorState)*)state;
    if (s !is null && result_out !is null) *result_out = (*s).C_Q;
}
extern(C) void oc_rotor_state_set_C_Q(OC_RotorState* state, double C_Q) {
    auto s = cast(RotorState*)state;
    if (s !is null) (*s).C_Q = C_Q;
}

extern(C) void oc_aircraft_state_get_rotor_C_T(const OC_AircraftState* state, size_t rotor_idx, double* result_out) {
    auto s = cast(const(AircraftState)*)state;
    if (s !is null && result_out !is null && rotor_idx < (*s).rotor_states.length) *result_out = (*s).rotor_states[rotor_idx].C_T;
}
extern(C) void oc_aircraft_state_get_rotor_C_Q(const OC_AircraftState* state, size_t rotor_idx, double* result_out) {
    auto s = cast(const(AircraftState)*)state;
    if (s !is null && result_out !is null && rotor_idx < (*s).rotor_states.length) *result_out = (*s).rotor_states[rotor_idx].C_Q;
}

// ========================================================================
//  BladeState API — fill data FROM D into caller's buffer
// ========================================================================

extern(C) void oc_blade_state_fill_dC_T(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_T"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_Db(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_Db"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_Db_profile(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dynamic_dC_Db_profile"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_Db_induced(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dynamic_dC_Db_induced"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dynamic_dC_Db_profile(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dynamic_dC_Db_profile"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dynamic_dC_Db_induced(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dynamic_dC_Db_induced"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_N(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_N"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_c(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_c"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_D(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_D"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_T_dot(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_T_dot"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_Q(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_Q"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_L(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_L"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_l(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_l"(data[0..len]);
}
// extern(C) void oc_blade_state_fill_dC_D_dot(const OC_BladeState* bs, double* data, size_t len) {
//     auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_D_dot"(data[0..len]);
// }
extern(C) void oc_blade_state_fill_dC_Mz(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_Mz"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dC_My(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_My"(data[0..len]);
}
extern(C) void oc_blade_state_fill_u_p(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"u_p"(data[0..len]);
}
extern(C) void oc_blade_state_fill_dynamic_u_p(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dynamic_u_p"(data[0..len]);
}
extern(C) void oc_blade_state_fill_u_t(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"u_t"(data[0..len]);
}
extern(C) void oc_blade_state_fill_aoa(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"aoa"(data[0..len]);
}
extern(C) void oc_blade_state_fill_aoa_eff(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"aoa_eff"(data[0..len]);
}
extern(C) void oc_blade_state_fill_gamma(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"gamma"(data[0..len]);
}
extern(C) void oc_blade_state_fill_r_c(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"r_c"(data[0..len]);
}
extern(C) void oc_blade_state_fill_x(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"x"(data[0..len]);
}
extern(C) void oc_blade_state_fill_y(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"y"(data[0..len]);
}
extern(C) void oc_blade_state_fill_z(const OC_BladeState* bs, double* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"z"(data[0..len]);
}

// Float variants
extern(C) void oc_blade_state_fill_dC_Df(const OC_BladeState* bs, float* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_D"(cast(float[])(data[0..len]));
}
extern(C) void oc_blade_state_fill_dC_Nf(const OC_BladeState* bs, float* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_N"(cast(float[])(data[0..len]));
}
extern(C) void oc_blade_state_fill_dC_cf(const OC_BladeState* bs, float* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_c"(cast(float[])(data[0..len]));
}
extern(C) void oc_blade_state_fill_dC_Tf(const OC_BladeState* bs, float* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_T"(cast(float[])(data[0..len]));
}
extern(C) void oc_blade_state_fill_dC_Qf(const OC_BladeState* bs, float* data, size_t len) {
    auto b = cast(const(BladeState)*)bs; if (b !is null && data !is null) (*b).get_state_array!"dC_Q"(cast(float[])(data[0..len]));
}

// Blade state scalar getters
extern(C) double oc_blade_state_get_azimuth(const OC_BladeState* bs) {
    auto b = cast(const(BladeState)*)bs; if (b !is null) return (*b).azimuth; return 0.0;
}
extern(C) double oc_blade_state_get_C_T(const OC_BladeState* bs) {
    auto b = cast(const(BladeState)*)bs; if (b !is null) return (*b).C_T; return 0.0;
}
extern(C) double oc_blade_state_get_C_Q(const OC_BladeState* bs) {
    auto b = cast(const(BladeState)*)bs; if (b !is null) return (*b).C_Q; return 0.0;
}
extern(C) double oc_blade_state_get_C_L(const OC_BladeState* bs) {
    auto b = cast(const(BladeState)*)bs; if (b !is null) return (*b).C_L; return 0.0;
}
extern(C) double oc_blade_state_get_C_D(const OC_BladeState* bs) {
    auto b = cast(const(BladeState)*)bs; if (b !is null) return (*b).C_D; return 0.0;
}

// ========================================================================
//  Wake API
// ========================================================================

extern(C) OC_Wake* oc_wake_create(size_t num_rotors, size_t num_blades,
                                   size_t wake_history, size_t radial_elements,
                                   const size_t* shed_history, const size_t* shed_release)
{
    auto wake = new Wake(
        num_rotors, num_blades, wake_history, radial_elements,
        (shed_history !is null) ? shed_history[0..num_rotors] : null,
        (shed_release !is null) ? shed_release[0..num_rotors] : null
    );
    GC.addRoot(wake);
    return cast(OC_Wake*)wake;
}
extern(C) void oc_wake_destroy(OC_Wake* wake) {
    if (wake !is null) GC.removeRoot(cast(Wake*)wake);
}

extern(C) void oc_vortex_filament_fill_x(const OC_VortexFilament* fil, double* data, size_t len) {
    auto f = cast(const(VortexFilament)*)fil;
    if (f !is null && data !is null) { auto arr = get_wake_component!"x"(*f); foreach(i, val; arr[0..len]) data[i] = val; }
}
extern(C) void oc_vortex_filament_fill_y(const OC_VortexFilament* fil, double* data, size_t len) {
    auto f = cast(const(VortexFilament)*)fil;
    if (f !is null && data !is null) { auto arr = get_wake_component!"y"(*f); foreach(i, val; arr[0..len]) data[i] = val; }
}
extern(C) void oc_vortex_filament_fill_z(const OC_VortexFilament* fil, double* data, size_t len) {
    auto f = cast(const(VortexFilament)*)fil;
    if (f !is null && data !is null) { auto arr = get_wake_component!"z"(*f); foreach(i, val; arr[0..len]) data[i] = val; }
}
extern(C) void oc_vortex_filament_fill_gamma(const OC_VortexFilament* fil, double* data, size_t len) {
    auto f = cast(const(VortexFilament)*)fil;
    if (f !is null && data !is null) { auto arr = get_wake_component!"gamma"(*f); foreach(i, val; arr[0..len]) data[i] = val; }
}
extern(C) void oc_vortex_filament_fill_r_c(const OC_VortexFilament* fil, double* data, size_t len) {
    auto f = cast(const(VortexFilament)*)fil;
    if (f !is null && data !is null) { auto arr = get_wake_component!"r_c"(*f); foreach(i, val; arr[0..len]) data[i] = val; }
}
extern(C) void oc_vortex_filament_fill_v_z(const OC_VortexFilament* fil, double* data, size_t len) {
    auto f = cast(const(VortexFilament)*)fil;
    if (f !is null && data !is null) { auto arr = get_wake_component!"v_z"(*f); foreach(i, val; arr[0..len]) data[i] = val; }
}

// ========================================================================
//  WakeHistory API
// ========================================================================

extern(C) OC_WakeHistory* oc_wake_history_create(
    size_t num_rotors, size_t num_blades,
    size_t wake_history, size_t time_history, size_t radial_elements,
    const size_t* shed_history, const(size_t)* shed_release,
    double a1, int hybrid)
{
    // Pass num_blades as an array (one value per rotor) to use the constructor
    // overload that avoids an uninit-array bug in the scalar-num_blades path.
    auto nb = new size_t[num_rotors];
    foreach(i; 0..num_rotors) nb[i] = num_blades;

    auto history = new WakeHistory(
        num_rotors, nb, wake_history, time_history, radial_elements,
        (shed_history !is null) ? shed_history[0..num_rotors] : null,
        (shed_release !is null) ? shed_release[0..num_rotors] : null,
        a1, hybrid != 0
    );
    GC.addRoot(history);
    return cast(OC_WakeHistory*)history;
}
extern(C) void oc_wake_history_destroy(OC_WakeHistory* history) {
    if (history !is null) GC.removeRoot(cast(WakeHistory*)history);
}
extern(C) void oc_wake_history_push_back(OC_WakeHistory* history) {
    auto h = cast(WakeHistory*)history;
    if (h !is null) h.push_back_wake();
}

// ========================================================================
//  Inflow Model API
// ========================================================================

// extern(C) void* oc_inflow_get_wrapped(OC_Inflow* inflow) {
//     auto i = cast(Inflow*)inflow;
//     if (i !is null) return (*i).get_wrapped_inflow();
//     return null;
// }
extern(C) double oc_inflow_wake_skew(OC_Inflow* inflow) {
    auto i = cast(Inflow)inflow;
    if (i !is null) return i.wake_skew();
    return 0.0;
}
extern(C) OC_Frame* oc_inflow_get_frame(OC_Inflow* inflow) {
    auto i = cast(Inflow)inflow;
    if (i !is null) return cast(OC_Frame*)i.frame();
    return null;
}
extern(C) const(OC_Mat4*) oc_inflow_get_inverse_global_frame(OC_Inflow* inflow) {
    auto i = cast(Inflow)inflow;
    if (i !is null) {
        auto f = i.frame();
        if (f !is null) return cast(const OC_Mat4*)&f.inverse_global_matrix;
    }
    return null;
}
extern(C) void oc_inflow_update(OC_Inflow* inflow, OC_AircraftState* ac_state, OC_Wake* wake, double dt) {
    auto i = cast(Inflow)inflow; auto s = cast(AircraftState*)ac_state; auto w = cast(Wake*)wake;
    if (i !is null && s !is null && w !is null) i.update(*s, *w, dt);
}

extern(C) void oc_inflow_at(OC_Inflow* inflow, const double* x, const double* y, const double* z, double* result_out, size_t len) {
    auto i = cast(Inflow)inflow;
    if (i !is null && x !is null && y !is null && z !is null && result_out !is null) {
        // Build Vector!(4, Chunk) from x/y/z arrays in chunks
        import opencopter.config : chunk_size;
        size_t processed = 0;
        while (processed < len) {
            size_t remaining = len - processed;
            size_t chunk_len = (remaining >= chunk_size) ? chunk_size : remaining;
            // Pad to chunk_size if needed
            Chunk _x, _y, _z;
            _x[] = 0; _y[] = 0; _z[] = 0;
            foreach(j; 0..chunk_len) {
                _x[j] = x[processed + j];
                _y[j] = y[processed + j];
                _z[j] = z[processed + j];
            }
            Vector!(4, Chunk) xyz_vec;
            xyz_vec[0][] = _x[];
            xyz_vec[1][] = _y[];
            xyz_vec[2][] = _z[];
            xyz_vec[3][] = 0;
            auto result_chunk = i.inflow_at(xyz_vec);
            foreach(j; 0..chunk_len) {
                result_out[(processed + j) * 3] = result_chunk[j];
            }
            processed += chunk_len;
        }
    }
}

// ========================================================================
//  New C bindings for feature parity with Python API
// ========================================================================

/**
 * Generate radial distribution points using the OpenCOPTER half-cosine method.
 * Writes results into the caller-provided buffer and returns the actual number
 * of points written (may be larger than requested due to chunk alignment).
 */
extern(C) double* oc_generate_radius_points(size_t* n_sections, double root_cutout) {
    if (n_sections !is null && *n_sections > 0) {
        auto pts = opencopter.aircraft.geometry.generate_radius_points(*n_sections, root_cutout);
        double* buf = cast(double*)malloc(pts.length*double.sizeof);
        *n_sections = pts.length;
        foreach(i; 0..pts.length) {
            buf[i] = pts[i];
        }
        return buf;
    }
    return null;
}

/**
 * Set children array for a Frame.
 */
extern(C) void oc_frame_set_children(OC_Frame* frame, OC_Frame** children, size_t num_children) {
    auto f = cast(Frame*)frame;
    if (f !is null && children !is null) {
        Frame*[] newChildren;
        for (size_t i = 0; i < num_children; i++) {
            auto child = cast(Frame*)children[i];
            if (child !is null) {
                newChildren ~= child;
                child.parent = f;
            }
        }
        f.children = newChildren;
    }
}

/**
 * Set rotors on an Aircraft.
 */
extern(C) void oc_aircraft_set_rotors(OC_Aircraft* ac, OC_RotorGeometry** rotors, size_t num_rotors) {
    auto a = cast(Aircraft*)ac;
    if (a !is null && rotors !is null) {
        RotorGeometry[] rotorArray;
        for (size_t i = 0; i < num_rotors; i++) {
            auto r = cast(RotorGeometry*)rotors[i];
            if (r !is null) {
                rotorArray ~= *r;
            }
        }
        a.rotors = rotorArray;
    }
}

/**
 * Set blade_length on a BladeGeometry.
 */
extern(C) void oc_blade_geometry_set_blade_length(OC_BladeGeometry* geom, double length) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null) g.blade_length = length;
}

/**
 * Set frame on a RotorGeometry.
 */
extern(C) void oc_rotor_geometry_set_frame(OC_RotorGeometry* rotor, OC_Frame* frame) {
    auto r = cast(RotorGeometry*)rotor;
    auto f = cast(Frame*)frame;
    if (r !is null && f !is null) r.frame = f;
}

/**
 * Get parent frame.
 */
extern(C) OC_Frame* oc_frame_get_parent(OC_Frame* frame) {
    auto f = cast(Frame*)frame;
    if (f !is null && f.parent !is null) return cast(OC_Frame*)f.parent;
    return null;
}

/**
 * Get children array for a Frame.
 * Returns a pointer to the first element of the children array.
 */
extern(C) OC_Frame** oc_frame_get_children(OC_Frame* frame) {
    auto f = cast(Frame*)frame;
    if (f !is null && f.children.length > 0) {
        // Cast the D array pointer to OC_Frame**
        return cast(OC_Frame**)f.children.ptr;
    }
    return null;
}

/**
 * Get the number of children for a Frame.
 */
extern(C) size_t oc_frame_get_children_count(OC_Frame* frame) {
    auto f = cast(Frame*)frame;
    if (f !is null) return f.children.length;
    return 0;
}

/**
 * Set frame type for a frame (as OC_FrameType enum int).
 */
extern(C) void oc_frame_set_frame_type(OC_Frame* frame, int frame_type) {
    auto f = cast(Frame*)frame;
    if (f !is null) f.frame_type = cast(FrameType)frame_type;
}

/**
 * Set name for a frame.
 */
extern(C) void oc_frame_set_name(OC_Frame* frame, const(char)* name) {
    auto f = cast(Frame*)frame;
    if (f !is null && name !is null) {
        f.name = fromStringz(name).idup;
    }
}
extern(C) void oc_inflow_update_wing_circulation(OC_Inflow* inflow, OC_WingState* wing_state) {
    auto i = cast(Inflow*)inflow; auto ws = cast(WingState*)wing_state;
    if (i !is null && ws !is null) (*i).update_wing_circulation(*ws);
}
extern(C) void oc_inflow_update_wing_dC_L(OC_Inflow* inflow, OC_WingState* wing_state) {
    auto i = cast(Inflow*)inflow; auto ws = cast(WingState*)wing_state;
    if (i !is null && ws !is null) (*i).update_wing_dC_L(*ws);
}
extern(C) OC_InducedVelocities oc_inflow_compute_wing_induced_vel_on_blade(OC_Inflow* inflow, const double* x, const double* y, const double* z) {
    auto i = cast(Inflow*)inflow; OC_InducedVelocities result;
    if (i !is null && x !is null && y !is null && z !is null) {
        auto iv = (*i).compute_wing_induced_vel_on_blade(x[0..8], y[0..8], z[0..8]);
        foreach(j; 0..8) { result.v_x[j] = iv.v_x[j]; result.v_y[j] = iv.v_y[j]; result.v_z[j] = iv.v_z[j]; }
    } else { foreach(j; 0..8) { result.v_x[j] = 0.0; result.v_y[j] = 0.0; result.v_z[j] = 0.0; } }
    return result;
}

extern(C) OC_Inflow* oc_huang_peters_create(long _Mo, long _Me, OC_RotorGeometry* rotor, OC_RotorInputState* rotor_input, double dt) {
    auto r = cast(RotorGeometry*)rotor; auto inp = cast(RotorInputState*)rotor_input;
    if (r !is null && inp !is null) {
        Inflow inflow = new HuangPetersInflow(_Mo, _Me, r, inp, dt);
        GC.addRoot(cast(void*)inflow);
        return cast(OC_Inflow*)inflow;
    }
    return null;
}
extern(C) OC_Inflow* oc_null_inflow_create(OC_RotorGeometry* rotor, OC_RotorInputState* rotor_input) {
    auto r = cast(RotorGeometry*)rotor;
    auto inp = cast(RotorInputState*)rotor_input;

    if (r !is null && inp !is null) {
        Inflow inflow = new NullInflow!(ArrayContainer.none)(r, inp);
        GC.addRoot(cast(void*)inflow);
        return cast(OC_Inflow*)inflow;
    }
    return null;
}
extern(C) OC_Inflow* oc_wing_inflow_create(OC_WingGeometry* wing, OC_WingInputState* wing_inputs, OC_WingLiftSurf* wing_lift_surf) {
    auto w = cast(WingGeometry*)wing; auto wi = cast(WingInputState*)wing_inputs; auto wl = cast(WingLiftSurf*)wing_lift_surf;
    if (w !is null && wi !is null && wl !is null) {
        Inflow inflow = new WingInflow(w, wi, wl);
        GC.addRoot(cast(void*)inflow);
        return cast(OC_Inflow*)inflow;
    }
    return null;
}

extern(C) void oc_inflow_destroy(OC_Inflow* inflow) {
    GC.removeRoot(inflow);
}

// ========================================================================
//  WingLiftSurf API
// ========================================================================

extern(C) void oc_wing_set_vortex_geometry(OC_WingLiftSurf* lift_surf, OC_WingGeometry* wing, size_t spanwise_chunks, size_t chordwise_nodes) {
    auto ls = cast(WingLiftSurf*)lift_surf; auto w = cast(WingGeometry*)wing;
    if (ls !is null && w !is null) set_wing_vortex_geometry(ls, w, spanwise_chunks, chordwise_nodes);
}

// ========================================================================
//  Simulation API
// ========================================================================

extern(C) void oc_simulation_step(OC_AircraftState* ac_state, OC_Aircraft* aircraft,
                                    OC_AircraftInputState* ac_input_state, OC_WakeHistory* wake_history,
                                    const OC_Atmosphere* atmo, size_t iteration, double dt,
                                    int track_bwi_events, int converged) {
    auto s = cast(AircraftState*)ac_state; auto a = cast(Aircraft*)aircraft;
    auto i = cast(AircraftInputState*)ac_input_state; auto h = cast(WakeHistory*)wake_history;
    if (s !is null && a !is null && i !is null && h !is null && atmo !is null) {
        immutable Atmosphere atmo_val = Atmosphere(atmo.density, atmo.dynamic_viscosity, atmo.speed_of_sound);
        step(*s, *a, *i, *h, atmo_val, iteration, dt, track_bwi_events != 0, converged != 0);
    }
}

extern(C) void oc_basic_aircraft_rotor_dynamics(OC_AircraftInputState* input, double dt) {
    auto i = cast(AircraftInputState*)input; if (i !is null) basic_aircraft_rotor_dynamics(i, dt);
}
extern(C) double oc_basic_single_rotor_dynamics(OC_RotorInputState* input, double dt) {
    auto i = cast(RotorInputState*)input; if (i !is null) return basic_single_rotor_dynamics(i, dt);
    return 0.0;
}

// ========================================================================
//  AircraftInputState helpers
// ========================================================================

extern(C) double get_blade_pitch(OC_AircraftInputState* input, size_t rotor_idx, size_t blade_idx) {
    auto i = cast(AircraftInputState*)input;
    if (i !is null && rotor_idx < (*i).rotor_inputs.length) return i.rotor_inputs[rotor_idx].blade_pitches[blade_idx];
    assert(0);
}

extern(C) void oc_aircraft_input_set_blade_pitch(OC_AircraftInputState* input, size_t rotor_idx, size_t blade_idx, double pitch) {
    auto i = cast(AircraftInputState*)input;
    if (i !is null && rotor_idx < (*i).rotor_inputs.length) i.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = pitch;
}

// ========================================================================
//  Utility Functions
// ========================================================================

extern(C) void oc_set_wing_ctrl_pt_geometry(OC_WingGeometry* wing, size_t spanwise_nodes, size_t chordwise_nodes, double camber) {
    auto w = cast(WingGeometry*)wing; if (w !is null) set_wing_ctrl_pt_geometry(w, spanwise_nodes, chordwise_nodes, camber);
}

// ========================================================================
//  VTK Output API
// ========================================================================

import opencopter.vtk;

extern(C) OC_VtkRotor* oc_build_vtu_rotor(OC_RotorGeometry* rotor) {
    if (rotor !is null) {
        auto r = cast(RotorGeometry*)rotor;
        auto vtk_rotor = build_base_vtu_rotor(*r);
        if (vtk_rotor !is null) GC.addRoot(cast(void*)vtk_rotor);
        return cast(OC_VtkRotor*)vtk_rotor;
    }
    return null;
}

/**
 * oc_write_rotor_vtu - Write a single rotor's VTU file.
 * Signature matches C header:
 *   void oc_write_rotor_vtu(const char*, size_t, size_t, OC_VtkRotor*, OC_AircraftState*, OC_RotorGeometry*);
 *
 * NOTE: VtkRotor is a D class (reference type). The OC_VtkRotor pointer from C
 * is the same bit pattern as the original class reference. Cast directly to VtkRotor,
 * NOT VtkRotor* (which would be a double-indirection).
 */
extern(C) void oc_write_rotor_vtu(const(char)* filename, size_t iteration, size_t rotor_idx,
                                    OC_VtkRotor* vtk_rotor, OC_AircraftState* ac_state,
                                    OC_RotorGeometry* rotor_geom) {
    auto fname = (filename !is null) ? fromStringz(filename).idup : "";
    // VtkRotor is a D class - cast directly to get the reference back
    auto v = cast(VtkRotor)vtk_rotor;
    auto s = cast(AircraftState*)ac_state;
    RotorGeometry* rg = cast(RotorGeometry*)rotor_geom;

    if (v is null || s is null || rg is null) return;
    if (rotor_idx >= (*s).rotor_states.length) return;

    auto rs = (*s).rotor_states[rotor_idx];
    // VtkRotor is a D class (reference type) - pass directly, no dereference needed
    write_rotor_vtu(fname, iteration, rotor_idx, v, rs, *rg);
}
/**
 * oc_write_rotors_vtu - Write all rotors' VTU files in a batch.
 * Iterates over the provided VtkRotor array and delegates to the
 * already-working singular write_rotor_vtu for each rotor.
 */
extern(C) void oc_write_rotors_vtu(const(char)* filename, size_t iteration,
                                    OC_VtkRotor** vtk_rotors, size_t num_rotors,
                                    OC_AircraftState* ac_state, OC_Aircraft* aircraft) {
    auto fname = (filename !is null) ? fromStringz(filename).idup : "";
    auto s = cast(AircraftState*)ac_state;
    auto a = cast(Aircraft*)aircraft;
    if (s is null || a is null || vtk_rotors is null) return;

    foreach (r_idx; 0 .. num_rotors) {
        if (vtk_rotors[r_idx] is null) continue;
        if (r_idx >= (*s).rotor_states.length) continue;
        if (r_idx >= a.rotors.length) continue;

        auto v = cast(VtkRotor)vtk_rotors[r_idx];
        auto rs = (*s).rotor_states[r_idx];
        auto rg = &a.rotors[r_idx];
        write_rotor_vtu(fname, iteration, r_idx, v, rs, *rg);
    }
}
extern(C) void oc_vtk_rotor_destroy(OC_VtkRotor* vtk) { if (vtk !is null) GC.removeRoot(cast(typeof(vtk)*)vtk); }
extern(C) OC_VtkWing* oc_build_vtu_wing(OC_WingGeometry* wing) { return null; }
extern(C) void oc_write_wing_vtu(const(char)*, size_t, size_t, OC_VtkWing*, OC_WingState*, OC_WingGeometry*) {}
extern(C) void oc_vtk_wing_destroy(OC_VtkWing* vtk) { if (vtk !is null) GC.removeRoot(cast(typeof(vtk)*)vtk); }
/**
 * oc_build_vtu_wake - Build a VtkWake from a Wake struct.
 * VtkWake is a D class (reference type), so the result is GC-pinned
 * and returned as an opaque C pointer.
 */
extern(C) OC_VtkWake* oc_build_vtu_wake(OC_Wake* wake) {
    if (wake !is null) {
        auto w = cast(Wake*)wake;
        auto vtk_wake = build_base_vtu_wake(*w);
        if (vtk_wake !is null) GC.addRoot(cast(void*)vtk_wake);
        return cast(OC_VtkWake*)vtk_wake;
    }
    return null;
}

/**
 * oc_write_wake_vtu - Write a wake's VTU file.
 * VtkWake is a D class - cast directly to get the reference back, no dereference.
 */
extern(C) void oc_write_wake_vtu(const(char)* filename, size_t iteration, OC_VtkWake* vtk_wake, OC_Wake* wake) {
    auto fname = (filename !is null) ? fromStringz(filename).idup : "";
    // VtkWake is a D class - cast directly to get the reference back
    auto v = cast(VtkWake)vtk_wake;
    auto w = cast(Wake*)wake;
    if (v is null || w is null) return;
    write_wake_vtu(fname, iteration, v, *w);
}
extern(C) void oc_vtk_wake_destroy(OC_VtkWake* vtk) { if (vtk !is null) GC.removeRoot(cast(typeof(vtk)*)vtk); }
extern(C) OC_VtkWingWake* oc_build_vtu_wing_wake(OC_WingGeometry*, OC_WingLiftSurf*) { return null; }
extern(C) void oc_write_wing_wake_vtu(const(char)*, size_t, size_t, OC_VtkWingWake*, OC_WingGeometry*, OC_WingLiftSurf*, OC_WingInputState*) {}
extern(C) void oc_vtk_wing_wake_destroy(OC_VtkWingWake* vtk) { if (vtk !is null) GC.removeRoot(cast(typeof(vtk)*)vtk); }
extern(C) void oc_write_wake_field_vtu(const(char)*, OC_AircraftState*, OC_Wake*, double, double, double, double, double, double, size_t, size_t, size_t) {}

// ========================================================================
//  Aircraft accessor functions
// ========================================================================

extern(C) OC_Frame* oc_aircraft_get_root_frame(OC_Aircraft* ac) {
    auto a = cast(Aircraft*)ac;
    if (a !is null) return cast(OC_Frame*)(*a).root_frame;
    return null;
}

// ========================================================================
//  RotorGeometry accessor and setter functions
// ========================================================================

extern(C) void oc_rotor_geometry_set_blades(OC_RotorGeometry* rotor, OC_BladeGeometry** blades, size_t num_blades) {
    auto r = cast(RotorGeometry*)rotor;
    if (r !is null && blades !is null) {
        // BladeGeometry is a struct allocated on the heap via `new`.
        // oc_blade_geometry_create returns cast(OC_BladeGeometry*)blade where blade = new BladeGeometry(...).
        // So we need to dereference the pointer to get the struct value.
        auto blade_array = new BladeGeometry[num_blades];
        for (size_t i = 0; i < num_blades; i++) {
            blade_array[i] = *cast(BladeGeometry*)(blades[i]);
        }
        r.blades = blade_array;
    }
}

// ========================================================================
//  BladeGeometry accessor and setter functions
// ========================================================================

extern(C) OC_Frame* oc_blade_geometry_get_frame(const OC_BladeGeometry* geom) {
    auto g = cast(BladeGeometry*)geom;
    if (g !is null) return cast(OC_Frame*)(*g).frame;
    return null;
}

extern(C) void oc_blade_geometry_set_frame(OC_BladeGeometry* geom, OC_Frame* frame) {
    auto g = cast(BladeGeometry*)geom;
    auto f = cast(Frame*)frame;
    if (g !is null && f !is null) (*g).frame = f;
}

// ========================================================================
//  AircraftInputState accessor and setter functions
// ========================================================================

extern(C) OC_RotorInputState* oc_aircraft_input_get_rotor_input(OC_AircraftInputState* input, size_t rotor_idx) {
    auto i = cast(AircraftInputState*)input;
    if (i !is null && rotor_idx < (*i).rotor_inputs.length)
        return cast(OC_RotorInputState*)(&(*i).rotor_inputs[rotor_idx]);
    return null;
}

// extern(C) void oc_aircraft_input_set_blade_pitch(OC_AircraftInputState* input, size_t rotor_idx, size_t blade_idx, double pitch) {
//     auto i = cast(AircraftInputState*)input;
//     if (i !is null && rotor_idx < (*i).rotor_inputs.length) 
//         i.rotor_inputs[rotor_idx].blade_pitches[blade_idx] = pitch;
// }

extern(C) double oc_aircraft_input_get_blade_pitch(OC_AircraftInputState* input, size_t rotor_idx, size_t blade_idx) {
    auto i = cast(AircraftInputState*)input;
    if (i !is null && rotor_idx < (*i).rotor_inputs.length) 
        return i.rotor_inputs[rotor_idx].blade_pitches[blade_idx];
    return 0.0;
}

extern(C) void oc_rotor_input_set_angular_velocity(OC_RotorInputState* input, double omega) {
    auto i = cast(RotorInputState*)input;
    if (i !is null) i.angular_velocity = omega;
}

extern(C) double oc_rotor_input_get_angular_velocity(OC_RotorInputState* input) {
    auto i = cast(RotorInputState*)input;
    if (i !is null) return i.angular_velocity;
    return 0.0;
}

extern(C) void oc_rotor_input_set_angular_accel(OC_RotorInputState* input, double alpha) {
    auto i = cast(RotorInputState*)input;
    if (i !is null) i.angular_accel = alpha;
}

extern(C) double oc_rotor_input_get_angular_accel(OC_RotorInputState* input) {
    auto i = cast(RotorInputState*)input;
    if (i !is null) return i.angular_accel;
    return 0.0;
}

extern(C) void oc_rotor_input_set_azimuth(OC_RotorInputState* input, double azimuth) {
    auto i = cast(RotorInputState*)input;
    if (i !is null) i.azimuth = azimuth;
}

extern(C) double oc_rotor_input_get_azimuth(OC_RotorInputState* input) {
    auto i = cast(RotorInputState*)input;
    if (i !is null) return i.azimuth;
    return 0.0;
}

extern(C) void oc_rotor_input_set_r_0(OC_RotorInputState* input, double* r_0, size_t len) {
    auto i = cast(RotorInputState*)input;
    if (i !is null && r_0 !is null) {
        foreach(idx; 0..min(len, i.r_0.length)) {
            i.r_0[idx] = r_0[idx];
        }
    }
}

extern(C) void oc_rotor_input_get_r_0(const OC_RotorInputState* input, double* result_out, size_t len) {
    auto i = cast(const(RotorInputState)*)input;
    if (i !is null && result_out !is null) {
        foreach(idx; 0..min(len, (*i).r_0.length)) {
            result_out[idx] = (*i).r_0[idx];
        }
    }
}

extern(C) void oc_rotor_input_set_blade_flapping(OC_RotorInputState* input, double* flapping, size_t len) {
    auto i = cast(RotorInputState*)input;
    if (i !is null && flapping !is null) {
        foreach(idx; 0..min(len, i.blade_flapping.length)) {
            i.blade_flapping[idx] = flapping[idx];
        }
    }
}

extern(C) void oc_rotor_input_get_blade_flapping(const OC_RotorInputState* input, double* result_out, size_t len) {
    auto i = cast(const(RotorInputState)*)input;
    if (i !is null && result_out !is null) {
        foreach(idx; 0..min(len, (*i).blade_flapping.length)) {
            result_out[idx] = (*i).blade_flapping[idx];
        }
    }
}

extern(C) void oc_rotor_input_set_blade_flapping_rate(OC_RotorInputState* input, double* flapping_rate, size_t len) {
    auto i = cast(RotorInputState*)input;
    if (i !is null && flapping_rate !is null) {
        foreach(idx; 0..min(len, i.blade_flapping_rate.length)) {
            i.blade_flapping_rate[idx] = flapping_rate[idx];
        }
    }
}

extern(C) void oc_rotor_input_get_blade_flapping_rate(const OC_RotorInputState* input, double* result_out, size_t len) {
    auto i = cast(const(RotorInputState)*)input;
    if (i !is null && result_out !is null) {
        foreach(idx; 0..min(len, (*i).blade_flapping_rate.length)) {
            result_out[idx] = (*i).blade_flapping_rate[idx];
        }
    }
}

// ========================================================================
//  AircraftState accessor and setter functions
// ========================================================================

extern(C) void oc_aircraft_state_set_freestream(OC_AircraftState* state, const OC_Vec4* freestream) {
    auto s = cast(AircraftState*)state;
    if (s !is null && freestream !is null) (*s).freestream = oc_vec4_to_vec4(*freestream);
}

extern(C) void oc_aircraft_state_get_freestream(const OC_AircraftState* state, OC_Vec4* result_out) {
    auto s = cast(const(AircraftState)*)state;
    if (s !is null && result_out !is null) *result_out = vec4_to_oc_vec4((*s).freestream);
}

// ========================================================================
//  WakeHistory accessor functions
// ========================================================================

extern(C) OC_Wake* oc_wake_history_get_wake(OC_WakeHistory* history, size_t idx) {
    auto h = cast(WakeHistory*)history;
    if (h !is null && idx < (*h).history.length) return cast(OC_Wake*)&(*h).history[idx];
    return null;
}

extern(C) OC_RotorWake* oc_wake_get_rotor_wake(OC_Wake* wake, size_t rotor_idx) {
    auto w = cast(Wake*)wake;
    if (w !is null && rotor_idx < (*w).rotor_wakes.length) return cast(OC_RotorWake*)&(*w).rotor_wakes[rotor_idx];
    return null;
}

extern(C) OC_VortexFilament* oc_rotor_wake_get_tip_vortex(OC_RotorWake* rotor_wake, size_t blade_idx) {
    auto rw = cast(RotorWake*)rotor_wake;
    if (rw !is null && blade_idx < (*rw).tip_vortices.length) return cast(OC_VortexFilament*)&(*rw).tip_vortices[blade_idx];
    return null;
}

// ========================================================================
//  Airfoil Model API
// ========================================================================

// --- ThinAirfoil ---

extern(C) OC_AirfoilModel* oc_thin_airfoil_create(double C_l_alpha_0) {
    auto af = new ThinAirfoil(C_l_alpha_0);
    GC.addRoot(cast(void*)af);
    return cast(OC_AirfoilModel*)af;
}

// --- AeroDAS ---

extern(C) OC_AirfoilModel* oc_aero_das_create(double* alpha, size_t alpha_len,
                                               double* CL, size_t cl_len,
                                               double* CD, size_t cd_len,
                                               double tbyc, double AR) {
    import std.stdio : writeln;

    if (alpha !is null && CL !is null && CD !is null && alpha_len == cl_len && alpha_len == cd_len) {
        debug writeln("Creating aerodas");
        try {
            auto af = new AeroDAS(alpha[0..alpha_len], CL[0..cl_len], CD[0..cd_len], tbyc, AR);
            debug writeln("Created aerodas");
            GC.addRoot(cast(void*)af);
            debug writeln("Rooted aerodas");
            return cast(OC_AirfoilModel*)af;
        } catch(Exception ex) {
            
            debug writeln("Caught exception creating aerodas airfoil: ", ex.msg);
            return null;
        }
    }
    return null;
}

extern(C) OC_AirfoilModel* oc_aero_das_from_xfoil_polar(const(char)* filename, double tbyc) {
    if (filename !is null) {
        try {
            string fname = fromStringz(filename).idup;
            auto af = create_aerodas_from_xfoil_polar(fname, tbyc);
            GC.addRoot(cast(void*)af);
            return cast(OC_AirfoilModel*)af;
        } catch(Exception ex) {
            return null;
        }
    }
    return null;
}

// --- C81 ---

extern(C) OC_AirfoilModel* oc_c81_from_file(const(char)* filename) {
    if (filename !is null) {
        try {
            string fname = fromStringz(filename).idup;
            auto af = load_c81_file(fname);
            GC.addRoot(cast(void*)af);
            return cast(OC_AirfoilModel*)af;
        } catch(Exception ex) {
            return null;
        }
    }
    return null;
}

// --- AirfoilModel destroy ---

extern(C) void oc_airfoil_model_destroy(OC_AirfoilModel* af) {
    if (af !is null) {
        auto p = cast(AirfoilModel)af;
        GC.removeRoot(cast(void*)p);
    }
}

// --- AirfoilModel query methods (scalar) ---

extern(C) double oc_airfoil_get_Cl(OC_AirfoilModel* af, double alpha, double mach) {
    auto a = cast(AirfoilModel)af;
    if (a !is null) return a.get_Cl(alpha, mach);
    return 0.0;
}

extern(C) double oc_airfoil_get_Cd(OC_AirfoilModel* af, double alpha, double mach) {
    auto a = cast(AirfoilModel)af;
    if (a !is null) return a.get_Cd(alpha, mach);
    return 0.0;
}

extern(C) double oc_airfoil_lift_curve_slope(OC_AirfoilModel* af) {
    auto a = cast(AirfoilModel)af;
    if (a !is null) return a.lift_curve_slope();
    return 0.0;
}

extern(C) double oc_airfoil_zero_lift_aoa(OC_AirfoilModel* af) {
    auto a = cast(AirfoilModel)af;
    if (a !is null) return a.zero_lift_aoa();
    return 0.0;
}

// --- BladeAirfoil ---

/**
 * Convenience function to create a BladeAirfoil with a single thin airfoil
 * model spanning all elements [0..num_elements). Uses the given lift-curve slope.
 */
extern(C) OC_BladeAirfoil* oc_blade_airfoil_create_basic(size_t num_elements, double C_l_alpha_0) {
    if (num_elements > 0) {
        auto af = new ThinAirfoil(C_l_alpha_0);
        GC.addRoot(cast(void*)af);

        AirfoilModel[] af_models;
        af_models ~= af;

        size_t[2][] ext;
        ext ~= [size_t(0), num_elements - 1];

        auto blade_af = new BladeAirfoil(af_models, ext);
        GC.addRoot(cast(void*)blade_af);
        return cast(OC_BladeAirfoil*)blade_af;
    }
    return null;
}

extern(C) OC_BladeAirfoil* oc_blade_airfoil_create(OC_AirfoilModel** models, const size_t* extents, size_t num_af) {
    if (models !is null && extents !is null && num_af > 0) {
        auto af_models = new AirfoilModel[num_af];
        for (size_t i = 0; i < num_af; i++) {
            af_models[i] = cast(AirfoilModel)models[i];
        }
        auto ext = new size_t[2][num_af];
        for (size_t i = 0; i < num_af; i++) {
            ext[i][0] = extents[i * 2];
            ext[i][1] = extents[i * 2 + 1];
        }
        auto blade_af = new BladeAirfoil(af_models, ext);
        GC.addRoot(cast(void*)blade_af);
        return cast(OC_BladeAirfoil*)blade_af;
    }
    return null;
}

extern(C) void oc_blade_airfoil_destroy(OC_BladeAirfoil* blade_af) {
    if (blade_af !is null) {
        auto p = cast(BladeAirfoil*)blade_af;
        GC.removeRoot(cast(void*)p);
    }
}

extern(C) double oc_blade_airfoil_get_Cl(OC_BladeAirfoil* blade_af, size_t chunk_idx, double alpha, double mach) {
    auto b = cast(BladeAirfoil)blade_af;
    if (b !is null) {
        Chunk _alpha, _mach;
        _alpha[] = alpha;
        _mach[] = mach;
        auto state = b.compute_coeffiecients(chunk_idx, _alpha, _mach);
        return state.C_l[0];
    }
    return 0.0;
}

extern(C) double oc_blade_airfoil_get_Cd(OC_BladeAirfoil* blade_af, size_t chunk_idx, double alpha, double mach) {
    auto b = cast(BladeAirfoil)blade_af;
    if (b !is null) {
        Chunk _alpha, _mach;
        _alpha[] = alpha;
        _mach[] = mach;
        auto state = b.compute_coeffiecients(chunk_idx, _alpha, _mach);
        return state.C_d[0];
    }
    return 0.0;
}

extern(C) double oc_blade_airfoil_lift_curve_slope(OC_BladeAirfoil* blade_af, size_t chunk_idx) {
    auto b = cast(BladeAirfoil)blade_af;
    if (b !is null) {
        auto val = b.lift_curve_slope(chunk_idx);
        return val[0];
    }
    return 0.0;
}

extern(C) double oc_blade_airfoil_zero_lift_aoa(OC_BladeAirfoil* blade_af, size_t chunk_idx) {
    auto b = cast(BladeAirfoil)blade_af;
    if (b !is null) {
        auto val = b.zero_lift_aoa(chunk_idx);
        return val[0];
    }
    return 0.0;
}

extern(C) void oc_blade_airfoil_fill_lift_curve_slope(OC_BladeAirfoil* blade_af, size_t chunk_idx, double* result_out, size_t len) {
    auto b = cast(BladeAirfoil)blade_af;
    if (b !is null && result_out !is null) {
        auto val = b.lift_curve_slope(chunk_idx);
        foreach(i; 0..min(len, val.length)) {
            result_out[i] = val[i];
        }
    }
}

extern(C) void oc_blade_airfoil_fill_zero_lift_aoa(OC_BladeAirfoil* blade_af, size_t chunk_idx, double* result_out, size_t len) {
    auto b = cast(BladeAirfoil)blade_af;
    if (b !is null && result_out !is null) {
        auto val = b.zero_lift_aoa(chunk_idx);
        foreach(i; 0..min(len, val.length)) {
            result_out[i] = val[i];
        }
    }
}

extern(C) void oc_blade_airfoil_fill_coefficients(OC_BladeAirfoil* blade_af, size_t chunk_idx,
                                                   const double* alphas, const double* machs,
                                                   double* Cl_out, double* Cd_out, size_t len) {
    auto b = cast(BladeAirfoil)blade_af;
    if (b !is null && alphas !is null && machs !is null && Cl_out !is null && Cd_out !is null) {
        import opencopter.config : chunk_size;
        size_t processed = 0;
        while (processed < len) {
            size_t remaining = len - processed;
            size_t chunk_len = (remaining >= chunk_size) ? chunk_size : remaining;
            Chunk _alpha, _mach;
            _alpha[] = 0;
            _mach[] = 0;
            foreach(j; 0..chunk_len) {
                _alpha[j] = alphas[processed + j];
                _mach[j] = machs[processed + j];
            }
            auto state = b.compute_coeffiecients(chunk_idx, _alpha, _mach);
            foreach(j; 0..chunk_len) {
                Cl_out[processed + j] = state.C_l[j];
                Cd_out[processed + j] = state.C_d[j];
            }
            processed += chunk_len;
        }
    }
}
