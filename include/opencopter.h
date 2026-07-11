#ifndef OPENCOPTER_H
#define OPENCOPTER_H

/**
 * C API for OpenCOPTER.
 *
 * This header provides a stable C interface to the OpenCOPTER simulation
 * library. The interface uses opaque pointers for all complex types,
 * allowing C++ (or any other language) consumers to interact with the
 * D library without needing D runtime dependencies or tooling.
 *
 * Author: OpenCOPTER Team
 * License: MIT
 */

#include <stddef.h>
#include <stdint.h>

// ========================================================================
//  Value Type Definitions
// ========================================================================

/**
 * 3D vector with double precision components
 */
typedef struct {
    double x;
    double y;
    double z;
} OC_Vec3;

/**
 * 4D vector with double precision components
 */
typedef struct {
    double x;
    double y;
    double z;
    double w;
} OC_Vec4;

/**
 * 3x3 matrix stored as a flat array of 9 doubles
 */
typedef struct {
    double data[9];
} OC_Mat3;

/**
 * 4x4 matrix stored as a flat array of 16 doubles
 */
typedef struct {
    double data[16];
} OC_Mat4;

/**
 * Atmospheric properties
 */
typedef struct {
    double density;
    double dynamic_viscosity;
    double kinematic_viscosity;
    double speed_of_sound;
} OC_Atmosphere;

/**
 * Induced velocities for 8 points
 */
typedef struct {
    double v_x[8];
    double v_y[8];
    double v_z[8];
} OC_InducedVelocities;

// ========================================================================
//  Enum Definitions
// ========================================================================

/**
 * Direction enumeration
 */
typedef enum {
    OC_CLOCKWISE = 0,
    OC_COUNTER_CLOCKWISE = 1
} OC_Direction;

/**
 * Frame type enumeration
 */
typedef enum {
    OC_AIRCRAFT_FRAME = 0,
    OC_CONNECTION_FRAME = 1,
    OC_ROTOR_FRAME = 2,
    OC_BLADE_FRAME = 3,
    OC_WING_FRAME = 4
} OC_FrameType;

// ========================================================================
//  Opaque Type Definitions
// ========================================================================

typedef struct OC_Aircraft OC_Aircraft;
typedef struct OC_AircraftInputState OC_AircraftInputState;
typedef struct OC_AircraftState OC_AircraftState;
typedef struct OC_BladeGeometry OC_BladeGeometry;
typedef struct OC_BladeState OC_BladeState;
typedef struct OC_BladeAirfoil OC_BladeAirfoil;
typedef struct OC_AirfoilModel OC_AirfoilModel;
typedef struct OC_Frame OC_Frame;
typedef struct OC_Inflow OC_Inflow;
typedef struct OC_HuangPeters OC_HuangPeters;
typedef struct OC_NullInflow OC_NullInflow;
typedef struct OC_WingInflow OC_WingInflow;
typedef struct OC_RotorGeometry OC_RotorGeometry;
typedef struct OC_RotorInputState OC_RotorInputState;
typedef struct OC_RotorState OC_RotorState;
typedef struct OC_RotorWake OC_RotorWake;
typedef struct OC_VortexFilament OC_VortexFilament;
typedef struct OC_Wake OC_Wake;
typedef struct OC_WakeHistory OC_WakeHistory;
typedef struct OC_WingGeometry OC_WingGeometry;
typedef struct OC_WingInputState OC_WingInputState;
typedef struct OC_WingLiftSurf OC_WingLiftSurf;
typedef struct OC_WingPartGeometry OC_WingPartGeometry;
typedef struct OC_WingState OC_WingState;
typedef struct OC_VtkRotor OC_VtkRotor;
typedef struct OC_VtkWing OC_VtkWing;
typedef struct OC_VtkWake OC_VtkWake;
typedef struct OC_VtkWingWake OC_VtkWingWake;

#ifdef __cplusplus
extern "C" {
#endif

// Helper function to create a vector
inline OC_Vec3 vec3(double x, double y, double z) {
    OC_Vec3 v;
    v.x = x;
    v.y = y;
    v.z = z;
    return v;
}

// ========================================================================
//  Configuration
// ========================================================================

size_t oc_chunk_size();

// ========================================================================
//  Direction helpers
// ========================================================================

OC_Direction oc_direction_clockwise();
OC_Direction oc_direction_counter_clockwise();

// ========================================================================
//  Matrix helpers
// ========================================================================

OC_Mat3 oc_mat3_identity();
OC_Mat4 oc_mat4_identity();

// ========================================================================
//  Frame API
// ========================================================================

OC_Frame* oc_frame_create(OC_Vec3 axis, double angle, OC_Vec3 translation,
                          OC_Frame* parent, const char* name, int frame_type);
void oc_frame_destroy(OC_Frame* frame);
void oc_frame_set_rotation(OC_Frame* frame, OC_Vec3 axis, double angle);
void oc_frame_rotate(OC_Frame* frame, OC_Vec3 axis, double angle);
void oc_frame_translate(OC_Frame* frame, OC_Vec3 translation);
void oc_frame_update(OC_Frame* frame, const OC_Mat4* parent_global_mat);
const OC_Mat4* oc_frame_get_local_matrix(const OC_Frame* frame);
const OC_Mat4* oc_frame_get_global_matrix(const OC_Frame* frame);
const OC_Mat4* oc_frame_get_inverse_global_matrix(const OC_Frame* frame);

/** Set children array for a Frame */
void oc_frame_set_children(OC_Frame* frame, OC_Frame** children, size_t num_children);

/** Get parent frame */
OC_Frame* oc_frame_get_parent(OC_Frame* frame);

/** Set frame type for a frame */
void oc_frame_set_frame_type(OC_Frame* frame, int frame_type);

/** Set name for a frame */
void oc_frame_set_name(OC_Frame* frame, const char* name);

// ========================================================================
//  Aircraft API
// ========================================================================

OC_Aircraft* oc_aircraft_create(size_t num_rotors, size_t num_wings);
void oc_aircraft_destroy(OC_Aircraft* ac);
OC_Frame* oc_aircraft_get_root_frame(OC_Aircraft* ac);

/** Set rotors on an Aircraft */
void oc_aircraft_set_rotors(OC_Aircraft* ac, OC_RotorGeometry** rotors, size_t num_rotors);

// ========================================================================
//  RotorGeometry API
// ========================================================================

OC_RotorGeometry* oc_rotor_geometry_create(size_t num_blades, OC_Vec3 origin,
                                           double radius, double solidity);
void oc_rotor_geometry_destroy(OC_RotorGeometry* geom);
void oc_rotor_geometry_set_solidity(OC_RotorGeometry* rotor, double solidity);
void oc_rotor_geometry_set_blades(OC_RotorGeometry* rotor, OC_BladeGeometry** blades, size_t num_blades);

/** Set frame on a RotorGeometry */
void oc_rotor_geometry_set_frame(OC_RotorGeometry* rotor, OC_Frame* frame);

// ========================================================================
//  BladeGeometry API
// ========================================================================

OC_BladeGeometry* oc_blade_geometry_create(size_t num_elements, double azimuth_offset,
                                           double average_chord, OC_BladeAirfoil* airfoil,
                                           double r_c);
void oc_blade_geometry_destroy(OC_BladeGeometry* geom);
void oc_blade_geometry_set_twist(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_set_chord(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_set_radius(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_set_C_l_alpha(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_set_alpha_0(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_set_sweep(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_set_xi(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_set_thickness(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_set_xi_p(OC_BladeGeometry* geom, double* data, size_t len);
void oc_blade_geometry_compute_vectors(OC_BladeGeometry* geom);
OC_Frame* oc_blade_geometry_get_frame(const OC_BladeGeometry* geom);
void oc_blade_geometry_set_frame(OC_BladeGeometry* geom, OC_Frame* frame);

/** Set blade_length on a BladeGeometry */
void oc_blade_geometry_set_blade_length(OC_BladeGeometry* geom, double length);

// ========================================================================
//  WingGeometry API
// ========================================================================

OC_WingGeometry* oc_wing_geometry_create(size_t num_parts, OC_Vec3 origin, double wing_span);
void oc_wing_geometry_destroy(OC_WingGeometry* geom);
void oc_wing_part_geometry_set_chord(OC_WingPartGeometry* wg, double* data, size_t len);
void oc_wing_part_geometry_set_twist(OC_WingPartGeometry* wg, double* data, size_t len);
void oc_wing_part_geometry_set_sweep(OC_WingPartGeometry* wg, double* data, size_t len);
void oc_wing_part_geometry_set_y_span(OC_WingPartGeometry* wg, double* data, size_t len);
void oc_wing_geometry_set_ctrl_points(OC_WingGeometry* wing, size_t spanwise_nodes,
                                      size_t chordwise_nodes, double camber);

// ========================================================================
//  AircraftInputState API
// ========================================================================

OC_AircraftInputState* oc_aircraft_input_state_create(size_t num_rotors, size_t* num_blades, size_t num_wings);
void oc_aircraft_input_state_destroy(OC_AircraftInputState* input);
OC_RotorInputState* oc_aircraft_input_get_rotor_input(OC_AircraftInputState* input, size_t rotor_idx);
void oc_aircraft_input_set_blade_pitch(OC_AircraftInputState* input, size_t rotor_idx, size_t blade_idx, double pitch);
double oc_aircraft_input_get_blade_pitch(OC_AircraftInputState* input, size_t rotor_idx, size_t blade_idx);

// ========================================================================
//  RotorInputState API
// ========================================================================

void oc_rotor_input_set_angular_velocity(OC_RotorInputState* input, double omega);
double oc_rotor_input_get_angular_velocity(OC_RotorInputState* input);
void oc_rotor_input_set_angular_accel(OC_RotorInputState* input, double alpha);
double oc_rotor_input_get_angular_accel(OC_RotorInputState* input);
void oc_rotor_input_set_azimuth(OC_RotorInputState* input, double azimuth);
double oc_rotor_input_get_azimuth(OC_RotorInputState* input);
void oc_rotor_input_set_r_0(OC_RotorInputState* input, double* r_0, size_t len);
void oc_rotor_input_get_r_0(const OC_RotorInputState* input, double* result_out, size_t len);
void oc_rotor_input_set_blade_flapping(OC_RotorInputState* input, double* flapping, size_t len);
void oc_rotor_input_get_blade_flapping(const OC_RotorInputState* input, double* result_out, size_t len);
void oc_rotor_input_set_blade_flapping_rate(OC_RotorInputState* input, double* flapping_rate, size_t len);
void oc_rotor_input_get_blade_flapping_rate(const OC_RotorInputState* input, double* result_out, size_t len);

// ========================================================================
//  AircraftState API
// ========================================================================

OC_AircraftState* oc_aircraft_state_create(
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
);
void oc_aircraft_state_destroy(OC_AircraftState* state);
void oc_aircraft_state_set_freestream(OC_AircraftState* state, const OC_Vec4* freestream);
void oc_aircraft_state_get_freestream(const OC_AircraftState* state, OC_Vec4* result_out);

// ========================================================================
//  State accessors
// ========================================================================

void oc_rotor_state_get_C_T(const OC_RotorState* state, double* result_out);
void oc_rotor_state_set_C_T(OC_RotorState* state, double C_T);
void oc_rotor_state_get_C_Q(const OC_RotorState* state, double* result_out);
void oc_rotor_state_set_C_Q(OC_RotorState* state, double C_Q);
void oc_aircraft_state_get_rotor_C_T(const OC_AircraftState* state, size_t rotor_idx, double* result_out);
void oc_aircraft_state_get_rotor_C_Q(const OC_AircraftState* state, size_t rotor_idx, double* result_out);

// ========================================================================
//  BladeState API
// ========================================================================

void oc_blade_state_fill_dC_T(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_Db(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_Db_profile(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_Db_induced(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dynamic_dC_Db_profile(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dynamic_dC_Db_induced(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_N(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_c(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_D(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_T_dot(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_Q(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_L(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_l(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_Mz(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dC_My(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_u_p(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_dynamic_u_p(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_u_t(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_aoa(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_aoa_eff(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_gamma(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_r_c(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_x(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_y(const OC_BladeState* bs, double* data, size_t len);
void oc_blade_state_fill_z(const OC_BladeState* bs, double* data, size_t len);

// Float variants
void oc_blade_state_fill_dC_Df(const OC_BladeState* bs, float* data, size_t len);
void oc_blade_state_fill_dC_Nf(const OC_BladeState* bs, float* data, size_t len);
void oc_blade_state_fill_dC_cf(const OC_BladeState* bs, float* data, size_t len);
void oc_blade_state_fill_dC_Tf(const OC_BladeState* bs, float* data, size_t len);
void oc_blade_state_fill_dC_Qf(const OC_BladeState* bs, float* data, size_t len);

// Blade state scalar getters
double oc_blade_state_get_azimuth(const OC_BladeState* bs);
double oc_blade_state_get_C_T(const OC_BladeState* bs);
double oc_blade_state_get_C_Q(const OC_BladeState* bs);
double oc_blade_state_get_C_L(const OC_BladeState* bs);
double oc_blade_state_get_C_D(const OC_BladeState* bs);

// ========================================================================
//  Wake API
// ========================================================================

OC_Wake* oc_wake_create(size_t num_rotors, size_t num_blades,
                        size_t wake_history, size_t radial_elements,
                        const size_t* shed_history, const size_t* shed_release);
void oc_wake_destroy(OC_Wake* wake);
void oc_vortex_filament_fill_x(const OC_VortexFilament* fil, double* data, size_t len);
void oc_vortex_filament_fill_y(const OC_VortexFilament* fil, double* data, size_t len);
void oc_vortex_filament_fill_z(const OC_VortexFilament* fil, double* data, size_t len);
void oc_vortex_filament_fill_gamma(const OC_VortexFilament* fil, double* data, size_t len);
void oc_vortex_filament_fill_r_c(const OC_VortexFilament* fil, double* data, size_t len);
void oc_vortex_filament_fill_v_z(const OC_VortexFilament* fil, double* data, size_t len);

// ========================================================================
//  WakeHistory API
// ========================================================================

OC_WakeHistory* oc_wake_history_create(
    size_t num_rotors, size_t num_blades,
    size_t wake_history, size_t time_history, size_t radial_elements,
    const size_t* shed_history, const size_t* shed_release,
    double a1, int hybrid);
void oc_wake_history_destroy(OC_WakeHistory* history);
void oc_wake_history_push_back(OC_WakeHistory* history);
OC_Wake* oc_wake_history_get_wake(OC_WakeHistory* history, size_t idx);

// ========================================================================
//  RotorWake / VortexFilament accessors
// ========================================================================

OC_RotorWake* oc_wake_get_rotor_wake(OC_Wake* wake, size_t rotor_idx);
OC_VortexFilament* oc_rotor_wake_get_tip_vortex(OC_RotorWake* rotor_wake, size_t blade_idx);

// ========================================================================
//  Inflow Model API
// ========================================================================

double oc_inflow_wake_skew(OC_Inflow* inflow);
OC_Frame* oc_inflow_get_frame(OC_Inflow* inflow);
const OC_Mat4* oc_inflow_get_inverse_global_frame(OC_Inflow* inflow);
void oc_inflow_update(OC_Inflow* inflow, OC_AircraftState* ac_state, OC_Wake* wake, double dt);
void oc_inflow_at(OC_Inflow* inflow, const double* x, const double* y, const double* z, double* result_out, size_t len);
void oc_inflow_update_wing_circulation(OC_Inflow* inflow, OC_WingState* wing_state);
void oc_inflow_update_wing_dC_L(OC_Inflow* inflow, OC_WingState* wing_state);
OC_InducedVelocities oc_inflow_compute_wing_induced_vel_on_blade(OC_Inflow* inflow, const double* x, const double* y, const double* z);
void oc_inflow_destroy(OC_Inflow* inflow);

// Inflow factory functions
OC_Inflow* oc_huang_peters_create(long _Mo, long _Me, OC_RotorGeometry* rotor, OC_RotorInputState* rotor_input, double dt);
OC_Inflow* oc_null_inflow_create(OC_RotorGeometry* rotor, OC_RotorInputState* rotor_input);
OC_Inflow* oc_wing_inflow_create(OC_WingGeometry* wing, OC_WingInputState* wing_inputs, OC_WingLiftSurf* wing_lift_surf);

// ========================================================================
//  WingLiftSurf API
// ========================================================================

void oc_wing_set_vortex_geometry(OC_WingLiftSurf* lift_surf, OC_WingGeometry* wing, size_t spanwise_chunks, size_t chordwise_nodes);

// ========================================================================
//  Simulation API
// ========================================================================

void oc_simulation_step(OC_AircraftState* ac_state, OC_Aircraft* aircraft,
                        OC_AircraftInputState* ac_input_state, OC_WakeHistory* wake_history,
                        const OC_Atmosphere* atmo, size_t iteration, double dt,
                        int track_bwi_events, int converged);
void oc_basic_aircraft_rotor_dynamics(OC_AircraftInputState* input, double dt);
double oc_basic_single_rotor_dynamics(OC_RotorInputState* input, double dt);

// ========================================================================
//  Utility Functions
// ========================================================================

double get_blade_pitch(OC_AircraftInputState* input, size_t rotor_idx, size_t blade_idx);
void oc_set_wing_ctrl_pt_geometry(OC_WingGeometry* wing, size_t spanwise_nodes, size_t chordwise_nodes, double camber);

/**
 * Generate radial distribution points using the OpenCOPTER half-cosine method.
 * Returns actual number of points written (may be padded for chunk alignment).
 */
double* oc_generate_radius_points(size_t* n_sections, double root_cutout);

// ========================================================================
//  VTK Output API
// ========================================================================

OC_VtkRotor* oc_build_vtu_rotor(OC_RotorGeometry* rotor);
void oc_write_rotor_vtu(const char*, size_t, size_t, OC_VtkRotor*, OC_RotorState*, OC_RotorGeometry*);
void oc_write_rotors_vtu(const char*, size_t, OC_VtkRotor**, size_t, OC_AircraftState*, OC_Aircraft*);
void oc_vtk_rotor_destroy(OC_VtkRotor* vtk);
OC_VtkWing* oc_build_vtu_wing(OC_WingGeometry* wing);
void oc_write_wing_vtu(const char*, size_t, size_t, OC_VtkWing*, OC_WingState*, OC_WingGeometry*);
void oc_vtk_wing_destroy(OC_VtkWing* vtk);
OC_VtkWake* oc_build_vtu_wake(OC_Wake* wake);
void oc_write_wake_vtu(const char*, size_t, OC_VtkWake*, OC_Wake*);
void oc_vtk_wake_destroy(OC_VtkWake* vtk);
OC_VtkWingWake* oc_build_vtu_wing_wake(OC_WingGeometry*, OC_WingLiftSurf*);
void oc_write_wing_wake_vtu(const char*, size_t, size_t, OC_VtkWingWake*, OC_WingGeometry*, OC_WingLiftSurf*, OC_WingInputState*);
void oc_vtk_wing_wake_destroy(OC_VtkWingWake* vtk);
void oc_write_wake_field_vtu(const char*, OC_AircraftState*, OC_Wake*, double, double, double, double, double, double, size_t, size_t, size_t);

// ========================================================================
//  Airfoil Model API
// ========================================================================

OC_AirfoilModel* oc_thin_airfoil_create(double C_l_alpha_0);
OC_AirfoilModel* oc_aero_das_create(double* alpha, size_t alpha_len,
                                    double* CL, size_t cl_len,
                                    double* CD, size_t cd_len,
                                    double tbyc, double AR);
OC_AirfoilModel* oc_aero_das_from_xfoil_polar(const char* filename, double tbyc);
OC_AirfoilModel* oc_c81_from_file(const char* filename);
void oc_airfoil_model_destroy(OC_AirfoilModel* af);
double oc_airfoil_get_Cl(const OC_AirfoilModel* af, double alpha, double mach);
double oc_airfoil_get_Cd(const OC_AirfoilModel* af, double alpha, double mach);
double oc_airfoil_lift_curve_slope(const OC_AirfoilModel* af);
double oc_airfoil_zero_lift_aoa(const OC_AirfoilModel* af);

// ========================================================================
//  BladeAirfoil API
// ========================================================================

OC_BladeAirfoil* oc_blade_airfoil_create_basic(size_t num_elements, double C_l_alpha_0);
OC_BladeAirfoil* oc_blade_airfoil_create(OC_AirfoilModel** models, const size_t* extents, size_t num_af);
void oc_blade_airfoil_destroy(OC_BladeAirfoil* blade_af);
double oc_blade_airfoil_get_Cl(const OC_BladeAirfoil* blade_af, size_t chunk_idx, double alpha, double mach);
double oc_blade_airfoil_get_Cd(const OC_BladeAirfoil* blade_af, size_t chunk_idx, double alpha, double mach);
double oc_blade_airfoil_lift_curve_slope(const OC_BladeAirfoil* blade_af, size_t chunk_idx);
double oc_blade_airfoil_zero_lift_aoa(const OC_BladeAirfoil* blade_af, size_t chunk_idx);
void oc_blade_airfoil_fill_lift_curve_slope(const OC_BladeAirfoil* blade_af, size_t chunk_idx, double* result_out, size_t len);
void oc_blade_airfoil_fill_zero_lift_aoa(const OC_BladeAirfoil* blade_af, size_t chunk_idx, double* result_out, size_t len);
void oc_blade_airfoil_fill_coefficients(const OC_BladeAirfoil* blade_af, size_t chunk_idx,
                                        const double* alphas, const double* machs,
                                        double* Cl_out, double* Cd_out, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* OPENCOPTER_H */