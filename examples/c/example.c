#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "opencopter.h"

// Mimic Python's build_blade function: create blade geometry + frame hierarchy
// Returns the blade geometry pointer; caller must free it.
OC_BladeGeometry* build_blade(size_t b_idx, double d_azimuth, double R, double r_c,
                              size_t elements, OC_BladeAirfoil* airfoil,
                              double avg_chord,
                              double* r, double* c, double* alpha_0,
                              double* xi, double* xi_p, double* twist,
                              double* C_l_alpha, double* sweep
    ) {

    // Create blade geometry
    OC_BladeGeometry* blade = oc_blade_geometry_create(elements, 0.0, avg_chord, airfoil, r_c);

    // Set blade properties from arrays
    oc_blade_geometry_set_radius(blade, r, elements);
    oc_blade_geometry_set_twist(blade, twist, elements);
    oc_blade_geometry_set_alpha_0(blade, alpha_0, elements);
    oc_blade_geometry_set_C_l_alpha(blade, C_l_alpha, elements);
    oc_blade_geometry_set_chord(blade, c, elements);
    oc_blade_geometry_set_sweep(blade, sweep, elements);
    oc_blade_geometry_set_xi(blade, xi, elements);
    oc_blade_geometry_set_xi_p(blade, xi_p, elements);

    // Set blade_length (matching Python: blade.blade_length = R*(1.0 - r_c))
    oc_blade_geometry_set_blade_length(blade, R * (1.0 - r_c));

    // Compute blade vectors
    oc_blade_geometry_compute_vectors(blade);

    // Build frame hierarchy mirroring Python:
    //   azimuth_offset_frame -> root_cutout_frame -> blade_frame

    // azimuth_offset_frame = Frame(Vec3([0,0,1]), b_idx*d_azimuth, Vec3([0,0,0]), null, ...)
    // Parent (rotor_frame) is not yet known here, so we pass NULL and link later via set_children.
    OC_Frame* azimuth_offset_frame = oc_frame_create(
        vec3(0.0, 0.0, 1.0),
        (double)b_idx * d_azimuth,
        vec3(0.0, 0.0, 0.0),
        NULL,
        "blade_azimuth",
        OC_CONNECTION_FRAME
    );

    // root_cutout_frame = Frame(Vec3([1,0,0]), 0, Vec3([R*r_c, 0, 0]), azimuth_offset_frame, ...)
    OC_Frame* root_cutout_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(R * r_c, 0.0, 0.0),
        azimuth_offset_frame,
        "blade_cutout",
        OC_CONNECTION_FRAME
    );

    // blade_frame = Frame(Vec3([1,0,0]), 0, Vec3([0,0,0]), root_cutout_frame, ...)
    OC_Frame* blade_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        root_cutout_frame,
        "blade",
        OC_BLADE_FRAME
    );

    // Set parent-child relationships:
    //   root_cutout_frame.children = [blade_frame]
    // (already linked via parent param, but set_children ensures bidirectional refs)
    OC_Frame* cutout_children[1] = { blade_frame };
    oc_frame_set_children(root_cutout_frame, cutout_children, 1);

    //   azimuth_offset_frame.children = [root_cutout_frame]
    OC_Frame* azimuth_children[1] = { root_cutout_frame };
    oc_frame_set_children(azimuth_offset_frame, azimuth_children, 1);

    // blade.frame = blade_frame
    oc_blade_geometry_set_frame(blade, blade_frame);

    return blade;
}

int main() {
    // Simulation parameters (matching Python example)
    size_t iterations = 5400;
    size_t wake_history_length = 1*1024;

    size_t requested_elements = 45;
    size_t num_rotors = 1;
    size_t num_blades = 4;
    double R = 2.0;

    double theta_75 = 3.2*(M_PI/180.0);

    double density = 1.125;
    double omega = 109.12;
    double AR = 16.5;
    double sos = 343.0;

    double V_inf = 32.9;
    double aoa = 6.5*(M_PI/180.0);
    double theta_tw_1 = -6.5*(M_PI/180.0);

    double d_psi = 1.0;
    double dt = d_psi*(M_PI/180.0)/fabs(omega);

    double shed_history_angle = 45.0;
    size_t shed_history_val = (size_t)round(shed_history_angle/d_psi);

    double d_azimuth = 2.0*M_PI/(double)num_blades;

    // Root cutout
    double r_c = 0.22;

    // Generate spanwise distributions using OpenCOPTER library function
    // (matching Python: r = generate_radius_points(requested_elements, r_c))
    // Allocate a large enough buffer (chunk-aligned size could be larger)
    size_t elements = requested_elements;
    double* r_buf = oc_generate_radius_points(&elements, r_c);
    printf("requested_elements: %zu, actual elements: %zu\n", requested_elements, elements);

    // Allocate remaining arrays using the actual element count
    double* c = (double*)malloc(elements * sizeof(double));
    double* alpha_0 = (double*)malloc(elements * sizeof(double));
    double* xi = (double*)malloc(elements * sizeof(double));
    double* xi_p = (double*)malloc(elements * sizeof(double));
    double* twist = (double*)malloc(elements * sizeof(double));
    double* C_l_alpha = (double*)malloc(elements * sizeof(double));
    double* sweep = (double*)malloc(elements * sizeof(double));

    // Chord distribution: c = (1.0/AR)*ones(elements)
    double avg_chord = 1.0 / AR;
    for (size_t i = 0; i < elements; i++) {
        c[i] = avg_chord;
    }

    for (size_t i = 0; i < elements; i++) {
        alpha_0[i] = 0.0;
    }

    for (size_t i = 0; i < elements; i++) {
        xi[i] = 0.0;
    }

    for (size_t i = 0; i < elements; i++) {
        xi_p[i] = 0.0;
    }

    // Twist: twist = [(_r - 0.75)*theta_tw_1 for _r in r]
    for (size_t i = 0; i < elements; i++) {
        twist[i] = (r_buf[i] - 0.75)*theta_tw_1;
    }

    for (size_t i = 0; i < elements; i++) {
        C_l_alpha[i] = 2.0*M_PI;
    }

    for (size_t i = 0; i < elements; i++) {
        sweep[i] = 0.0;
    }

    // Create atmosphere
    OC_Atmosphere atmo;
    atmo.density = density;
    atmo.dynamic_viscosity = 18.03e-6;
    atmo.kinematic_viscosity = 0.0;
    atmo.speed_of_sound = sos;

    // Create aircraft (matching Python: Aircraft(num_rotors, 0))
    OC_Aircraft* aircraft = oc_aircraft_create(num_rotors, 0);

    // Origins of rotors (matching Python)
    OC_Vec3 origin_0 = vec3(0.0, 0.0, 0.0);

    // Create blade airfoil using the aerodas model
    const char* polar_path = "../polars/NACA23012mod_1000000_polar.dat";
    OC_AirfoilModel* af = oc_aero_das_from_xfoil_polar(polar_path, 0.12);
    OC_AirfoilModel* models[] = { af };
    size_t extents[2] = {0, 47};
    OC_BladeAirfoil* airfoil = oc_blade_airfoil_create(models, extents, 1);


    // Build rotor (matching Python's build_rotor function)
    OC_RotorGeometry* rotor = oc_rotor_geometry_create(num_blades, origin_0, R, 0.0);

    // Create rotor frame: Frame(Vec3([1,0,0]), 0, Vec3([0,0,0]), aircraft.root_frame, 'rotor_0', FrameType.rotor())
    OC_Frame* root_frame = oc_aircraft_get_root_frame(aircraft);

    // rotor_frame = Frame with aircraft root_frame as parent
    OC_Frame* rotor_fixed_frame = oc_frame_create(
        vec3(1.0, 0.0, 0.0),
        0.0,
        vec3(0.0, 0.0, 0.0),
        root_frame,
        "rotor_0_fixed",
        OC_CONNECTION_FRAME
    );

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

    OC_Frame* rotor_fixed_children[] = { rotor_frame };
    oc_frame_set_children(rotor_fixed_frame, rotor_fixed_children, 1);

    // Build blades with proper frame hierarchy
    OC_BladeGeometry** blades = (OC_BladeGeometry**)malloc(num_blades * sizeof(OC_BladeGeometry*));
    OC_Frame** azimuth_frames = (OC_Frame**)malloc(num_blades * sizeof(OC_Frame*));

    for (size_t b_idx = 0; b_idx < num_blades; b_idx++) {
        blades[b_idx] =
            build_blade(
                b_idx,
                d_azimuth,
                R,
                r_c,
                elements,
                airfoil,
                R * avg_chord,
                r_buf,
                c,
                alpha_0,
                xi,
                xi_p,
                twist,
                C_l_alpha,
                sweep
            );

        // Store the azimuth_offset_frame (parent.parent of blade_frame) for rotor_frame.children
        // In build_blade: blade->frame = blade_frame, blade_frame.parent = root_cutout_frame,
        // root_cutout_frame.parent = azimuth_offset_frame
        // So we need oc_frame_get_parent(oc_frame_get_parent(blade_frame))
        OC_Frame* blade_frame = oc_blade_geometry_get_frame(blades[b_idx]);
        OC_Frame* root_cutout_frame = oc_frame_get_parent(blade_frame);
        OC_Frame* az_frame = oc_frame_get_parent(root_cutout_frame);
        azimuth_frames[b_idx] = az_frame;
    }

    // Set blades for rotor
    oc_rotor_geometry_set_blades(rotor, blades, num_blades);

    // Calculate and set solidity (matching Python: num_blades*blade[0].average_chord/(pi*radius))
    double solidity = num_blades * (R * avg_chord) / (M_PI * R);
    printf("solidity: %f\n", solidity);
    oc_rotor_geometry_set_solidity(rotor, solidity);

    // rotor_frame.children = [b.frame.parent.parent for b in rotor.blades]  (azimuth_frames)
    oc_frame_set_children(rotor_frame, azimuth_frames, num_blades);

    // rotor.frame = rotor_frame
    oc_rotor_geometry_set_frame(rotor, rotor_frame);

    // aircraft.rotors = [rotor]  and  aircraft.root_frame.children = [rotor.frame]
    OC_RotorGeometry* rotors_arr[1] = { rotor };
    oc_aircraft_set_rotors(aircraft, rotors_arr, num_rotors);

    // OC_Frame* root_children[1] = { rotor_frame };
    // oc_frame_set_children(root_frame, root_children, 1);

    // Set aircraft frame properties (matching Python)
    oc_frame_set_name(root_frame, "Aircraft frame");
    oc_frame_set_frame_type(root_frame, OC_AIRCRAFT_FRAME);

    // Apply angle of attack rotation and update transform tree
    oc_frame_rotate(root_frame, vec3(0.0, 1.0, 0.0), -aoa);
    OC_Mat4 identity_mat = oc_mat4_identity();
    oc_frame_update(root_frame, &identity_mat);

    // Create input state (matching Python)
    size_t num_blades_arr = num_blades;
    OC_AircraftInputState* ac_input_state = oc_aircraft_input_state_create(num_rotors, &num_blades_arr, 0);

    // Set rotor inputs
    OC_RotorInputState* rotor_input = oc_aircraft_input_get_rotor_input(ac_input_state, 0);

    oc_rotor_input_set_angular_velocity(rotor_input, omega);
    oc_rotor_input_set_angular_accel(rotor_input, 0.0);
    oc_rotor_input_set_azimuth(rotor_input, 0.0);

    // Set r_0 for each blade (matching Python: 1.0*aircraft.rotors[0].blades[b_idx].average_chord)
    double* r_0_array = (double*)malloc(num_blades * sizeof(double));
    for (size_t b_idx = 0; b_idx < num_blades; b_idx++) {
        r_0_array[b_idx] = 1.0 * R * avg_chord;
        printf("r_0_array[b_idx]: %f\n", r_0_array[b_idx]);
    }
    oc_rotor_input_set_r_0(rotor_input, r_0_array, num_blades);

    // Set blade flapping and pitch
    double* blade_flapping = (double*)malloc(num_blades * sizeof(double));
    double* blade_flapping_rate = (double*)malloc(num_blades * sizeof(double));
    for (size_t b_idx = 0; b_idx < num_blades; b_idx++) {
        blade_flapping[b_idx] = 0.0;
        blade_flapping_rate[b_idx] = 0.0;
        oc_aircraft_input_set_blade_pitch(ac_input_state, 0, b_idx, theta_75);
    }

    oc_rotor_input_set_blade_flapping(rotor_input, blade_flapping, num_blades);
    oc_rotor_input_set_blade_flapping_rate(rotor_input, blade_flapping_rate, num_blades);

    free(r_0_array);
    free(blade_flapping);
    free(blade_flapping_rate);

    // Create inflow models
    OC_Inflow* inflow = oc_huang_peters_create(4, 2, rotor, rotor_input, dt);
    OC_Inflow** inflows = &inflow;

    // Set rotation direction for each rotor (sign of omega)
    double direction[1];
    direction[0] = copysign(1.0, omega);

    // Create aircraft state
    size_t nb_arr[1] = {num_blades};
    OC_AircraftState* ac_state =
        oc_aircraft_state_create(
            num_rotors,
            nb_arr,
            elements,
            0,
            NULL,
            0,
            0,
            aircraft,
            inflows,
            NULL,
            direction
        );

    // Set freestream velocity (matching Python: Vec4([V_inf, 0, 0, 0]))
    OC_Vec4 freestream;
    freestream.x = V_inf;
    freestream.y = 0.0;
    freestream.z = 0.0;
    freestream.w = 0.0;
    oc_aircraft_state_set_freestream(ac_state, &freestream);

    // Setup wake history (matching Python)
    size_t shed_history_arr[1] = {shed_history_val};
    size_t shed_release_arr[1] = {1};

    double a1 = 6.5e-5;
    bool hybrid = false;
    OC_WakeHistory* wake_history =
        oc_wake_history_create(
            num_rotors,
            num_blades_arr,
            wake_history_length,
            2,
            elements,
            shed_history_arr,
            shed_release_arr,
            a1,
            hybrid
        );

    // Perform simulation
    printf("Starting simulation with %zu iterations\n", iterations);
    double C_T;
    for (size_t iteration = 0; iteration < iterations; iteration++) {
        if (iteration % 100 == 0) {
            printf("Iteration: %zu\n", iteration);
        }

        // Advance rotor dynamics
        oc_basic_aircraft_rotor_dynamics(ac_input_state, dt);
        double azimuth = oc_rotor_input_get_azimuth(rotor_input);
        oc_frame_set_rotation(
            rotor_frame,
            vec3(0, 0, 1),
            azimuth
        );

        //printf("azimuth: %f\n", azimuth);

        // Run simulation step
        oc_simulation_step(
            ac_state,
            aircraft,
            ac_input_state,
            wake_history,
            &atmo,
            iteration,
            dt,
            0,
            0
        );

    }

    // Print results
    oc_aircraft_state_get_rotor_C_T(ac_state, 0, &C_T);
    printf("Rotor 0 C_T: %f\n", C_T);

    // Cleanup frames (children first to avoid dangling parent pointers)
    for (size_t b_idx = 0; b_idx < num_blades; b_idx++) {
        OC_Frame* blade_frame = oc_blade_geometry_get_frame(blades[b_idx]);
        OC_Frame* root_cutout_frame = oc_frame_get_parent(blade_frame);
        OC_Frame* az_frame = oc_frame_get_parent(root_cutout_frame);

        oc_frame_destroy(blade_frame);
        oc_frame_destroy(root_cutout_frame);
        oc_frame_destroy(az_frame);
    }

    oc_frame_destroy(rotor_frame);

    oc_inflow_destroy(inflow);

    for (size_t b_idx = 0; b_idx < num_blades; b_idx++) {
        oc_blade_geometry_destroy(blades[b_idx]);
    }
    free(blades);
    free(azimuth_frames);

    oc_rotor_geometry_destroy(rotor);
    oc_blade_airfoil_destroy(airfoil);
    oc_aircraft_state_destroy(ac_state);
    oc_aircraft_input_state_destroy(ac_input_state);
    oc_wake_history_destroy(wake_history);
    oc_aircraft_destroy(aircraft);

    free(r_buf);
    free(c);
    free(alpha_0);
    free(xi);
    free(xi_p);
    free(twist);
    free(C_l_alpha);
    free(sweep);

    printf("Simulation completed successfully!\n");

    return 0;
}