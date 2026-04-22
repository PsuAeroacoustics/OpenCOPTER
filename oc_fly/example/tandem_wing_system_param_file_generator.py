import itertools

# -----------------------------
# Input parameter arrays
# -----------------------------
tilt_angles = [-90]
V_inf_values = [30, 40]
aoa_values = [0]

# -----------------------------
# Template for one flight condition
# -----------------------------
def create_flight_condition(tilt, V_inf, aoa):

    name = f"Tilt_{tilt}_Vinf_{V_inf}_AOA_{aoa}"

    return f"""
        {{
            name: "{name}",
            aoa: {aoa},
            density: 1.25,
            sos: 338.0,
            c_t: [0.02, 0.02],
            collective_force_components: [2, 2],
            V_inf: {V_inf},
            run_for: 25.0,
            observer_file: "externalFile_observers.json5",
            acoustics_file: "acoustics_total_noise.json",
            motion: [
                {{
                    frame: "main_wing",
                    axis_angle_function: "static",
                    angle: 0.0
                }},
                {{
                    frame: "archer_tilter_0",
                    blade_element_func: "rotation",
                    omega: -110.0
                }},
                {{
                    frame: "archer_lifter_0",
                    blade_element_func: "rotation",
                    omega: 100
                }},
                {{
                    frame: "tilter_0_tilt_hinge",
                    axis_angle_function: "static",
                    angle: {tilt}
                }},
                {{
                    frame: "lifter_0_tilt_hinge",
                    axis_angle_function: "static",
                    angle: 0.0
                }}
            ]
        }}
    """


# -----------------------------
# Generate all combinations
# -----------------------------
flight_conditions_list = []

for tilt, V_inf, aoa in itertools.product(tilt_angles, V_inf_values, aoa_values):
    fc_str = create_flight_condition(tilt, V_inf, aoa)
    flight_conditions_list.append(fc_str)

# Join with commas
flight_conditions_block = ",\n".join(flight_conditions_list)

# -----------------------------
# Full JSON5 file content
# -----------------------------
full_content = f"""
{{
    computational_parameters: {{
        d_psi: 1.0,
        spanwise_elements: 48,
        span_elements: 32,
        chord_elements: 6,
        shed_history_angle: [15, 15],
        wake_trail_distance: [5, 10],
        num_trailing_vortices: [7, 7],
        shed_release_angle: 1.0,
        convergence_criteria: 1.0e-2,
        trim_algo: "lympany",
        convergence_type: "run_for",
        r_0: [1.0, 1.0],
        a1: 1.0e-5,
        post_conv_revolutions: 2,
        is_half_wing: true,
        trackBWIevents: false
    }},

    results: {{
        inflow_slices: [
            {{
                resolution: [256, 256, 101],
                slice_start: [-2.0, 0.0, -1.0],
                slice_size: [4.0, 4, 1.0]
            }}
        ],
        wake_slices: [
            {{
                resolution: [456, 152, 80],
                slice_start: [-3.0, -3.0, -1],
                slice_size: [18, 6, 3]
            }}
        ],
        spanwise_time_series: [0.87]
    }},

    flight_conditions: [
{flight_conditions_block}
    ],
}}
"""

# -----------------------------
# Write to file
# -----------------------------
output_file = "generated_flight_conditions.json5"

with open(output_file, "w") as f:
    f.write(full_content)

print(f"File written: {output_file}")