/**
 * C++ Wrapper for OpenCOPTER.
 *
 * This header provides an idiomatic C++ interface to the OpenCOPTER simulation
 * library. All types are native C++ — no opaque OC_* types leak into the public
 * API.  RAII memory management via private void* pointers with explicit destructors.
 * Constructors replace factory functions for a natural object-creation model.
 *
 * Author: OpenCOPTER Team
 * License: MIT
 */

#ifndef OPENCOPTER_HPP_
#define OPENCOPTER_HPP_

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <span>

namespace opencopter {

// ========================================================================
//  Value Types (native C++ — no C typedef aliases)
// ========================================================================

struct Vec3 { double x, y, z; };
struct Vec4 { double x, y, z, w; };
struct Mat3 { double data[9]; };
struct Mat4 { double data[16]; };

struct AtmosphereData {
    double density;
    double dynamic_viscosity;
    double kinematic_viscosity;
    double speed_of_sound;
};

struct InducedVelocities {
    double v_x[8];
    double v_y[8];
    double v_z[8];
};

// ========================================================================
//  Enum Classes (matching D enums)
// ========================================================================

enum class Direction : int { Clockwise = 0, CounterClockwise = 1 };

enum class FrameType : int {
    Aircraft   = 0,
    Connection = 1,
    Rotor      = 2,
    Blade      = 3,
    Wing       = 4
};

// ========================================================================
//  Forward declarations of wrapper classes
// ========================================================================

class Frame;
class Aircraft;
class RotorGeometry;
class BladeGeometry;
class BladeAirfoil;
class AirfoilModel;
class WingGeometry;
class WingPartGeometry;
class AircraftInputState;
class RotorInputState;
class AircraftState;
class RotorState;
class BladeState;
class Wake;
class WakeHistory;
class RotorWake;
class VortexFilament;
class Inflow;
class HuangPetersInflow;
class NullInflow;
class WingInflow;
class WingInputState;
class WingState;
class WingLiftSurf;
class VtkRotor;
class VtkWing;
class VtkWake;
class VtkWingWake;

// ========================================================================
//  Internal: pointer accessor (friend struct for implementation file)
//     Defined at the bottom of this header after all classes are complete.
// ========================================================================

namespace detail {
    struct ptr_accessor;  // forward declared, defined below all classes
}

// ========================================================================
//  Frame
// ========================================================================

class Frame {
protected:
    void* ptr_ = nullptr;
    bool owned_ = false;

    friend struct detail::ptr_accessor;
    friend class Aircraft;
    friend class RotorGeometry;
    friend class BladeGeometry;
    friend class Inflow;
    friend class VortexFilament;  // not needed but safe

    explicit Frame(void*, bool owned = true);

public:
    Frame() = default;

    Frame(Vec3 axis, double angle, Vec3 translation,
          const Frame* parent, std::string_view name, FrameType frame_type);

    ~Frame();

    Frame(const Frame&) = delete;
    Frame& operator=(const Frame&) = delete;
    Frame(Frame&& other) noexcept;
    Frame& operator=(Frame&& other) noexcept;

    void set_rotation(Vec3 axis, double angle);
    void rotate(Vec3 axis, double angle);
    void translate(Vec3 translation);
    void update(const Mat4& parent_global_mat);

    void set_children(std::span<const Frame*> children);
    void set_frame_type(FrameType ft);
    void set_name(std::string_view name);

    const Mat4* local_matrix() const;
    const Mat4* global_matrix() const;
    const Mat4* inverse_global_matrix() const;

    Frame parent();

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  Aircraft
// ========================================================================

class Aircraft {
private:
    friend struct detail::ptr_accessor;
    friend class AircraftState;

    void* ptr_ = nullptr;

    explicit Aircraft(void*);

public:
    Aircraft() = default;
    Aircraft(size_t num_rotors, size_t num_wings);

    ~Aircraft();

    Aircraft(const Aircraft&) = delete;
    Aircraft& operator=(const Aircraft&) = delete;
    Aircraft(Aircraft&&) noexcept;
    Aircraft& operator=(Aircraft&&) noexcept;

    Frame root_frame();
    void set_rotors(std::span<const RotorGeometry*> rotors);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  AirfoilModel
// ========================================================================

class AirfoilModel {
private:
    friend struct detail::ptr_accessor;
    friend class BladeAirfoil;

    void* ptr_ = nullptr;

    explicit AirfoilModel(void*);

public:
    AirfoilModel() = default;

    static AirfoilModel thin_airfoil(double C_l_alpha_0);
    static AirfoilModel aero_das(const std::vector<double>& alpha,
                                 const std::vector<double>& CL,
                                 const std::vector<double>& CD,
                                 double tbyc, double AR);
    static AirfoilModel aero_das_from_xfoil_polar(std::string_view filename, double tbyc);
    static AirfoilModel c81_from_file(std::string_view filename);

    ~AirfoilModel();

    AirfoilModel(const AirfoilModel&) = delete;
    AirfoilModel& operator=(const AirfoilModel&) = delete;
    AirfoilModel(AirfoilModel&&) noexcept;
    AirfoilModel& operator=(AirfoilModel&&) noexcept;

    double get_Cl(double alpha, double mach) const;
    double get_Cd(double alpha, double mach) const;
    double lift_curve_slope() const;
    double zero_lift_aoa() const;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  BladeAirfoil
// ========================================================================

class BladeAirfoil {
private:
    friend struct detail::ptr_accessor;
    friend class BladeGeometry;

    void* ptr_ = nullptr;

    explicit BladeAirfoil(void*);

public:
    BladeAirfoil() = default;

    static BladeAirfoil create_basic(size_t num_elements, double C_l_alpha_0);
    static BladeAirfoil create(const std::vector<AirfoilModel>& models,
                              const std::vector<size_t>& extents);

    ~BladeAirfoil();

    BladeAirfoil(const BladeAirfoil&) = delete;
    BladeAirfoil& operator=(const BladeAirfoil&) = delete;
    BladeAirfoil(BladeAirfoil&&) noexcept;
    BladeAirfoil& operator=(BladeAirfoil&&) noexcept;

    double get_Cl(size_t chunk_idx, double alpha, double mach) const;
    double get_Cd(size_t chunk_idx, double alpha, double mach) const;
    double lift_curve_slope(size_t chunk_idx) const;
    double zero_lift_aoa(size_t chunk_idx) const;

    std::vector<double> fill_lift_curve_slope(size_t chunk_idx) const;
    std::vector<double> fill_zero_lift_aoa(size_t chunk_idx) const;

    void fill_coefficients(size_t chunk_idx,
                           const std::vector<double>& alphas,
                           const std::vector<double>& machs,
                           std::vector<double>& Cl_out,
                           std::vector<double>& Cd_out) const;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  BladeGeometry
// ========================================================================

class BladeGeometry {
private:
    friend struct detail::ptr_accessor;
    friend class RotorGeometry;

    void* ptr_ = nullptr;

    explicit BladeGeometry(void*);

public:
    BladeGeometry() = default;

    BladeGeometry(size_t num_elements, double azimuth_offset,
                  double average_chord, const BladeAirfoil& airfoil,
                  double r_c);

    ~BladeGeometry();

    BladeGeometry(const BladeGeometry&) = delete;
    BladeGeometry& operator=(const BladeGeometry&) = delete;
    BladeGeometry(BladeGeometry&&) noexcept;
    BladeGeometry& operator=(BladeGeometry&&) noexcept;

    void set_twist(const std::vector<double>& data);
    void set_chord(const std::vector<double>& data);
    void set_radius(const std::vector<double>& data);
    void set_C_l_alpha(const std::vector<double>& data);
    void set_alpha_0(const std::vector<double>& data);
    void set_sweep(const std::vector<double>& data);
    void set_xi(const std::vector<double>& data);
    void set_thickness(const std::vector<double>& data);
    void set_xi_p(const std::vector<double>& data);

    void compute_vectors();

    Frame get_frame() const;
    void set_frame(const Frame& frame);
    void set_blade_length(double length);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  RotorGeometry
// ========================================================================

class RotorGeometry {
private:
    friend struct detail::ptr_accessor;
    friend class Aircraft;
    friend class HuangPetersInflow;
    friend class NullInflow;
    friend class VtkRotor;

    void* ptr_ = nullptr;

    explicit RotorGeometry(void*);

public:
    RotorGeometry() = default;

    RotorGeometry(size_t num_blades, Vec3 origin,
                  double radius, double solidity);

    ~RotorGeometry();

    RotorGeometry(const RotorGeometry&) = delete;
    RotorGeometry& operator=(const RotorGeometry&) = delete;
    RotorGeometry(RotorGeometry&&) noexcept;
    RotorGeometry& operator=(RotorGeometry&&) noexcept;

    void set_solidity(double solidity);
    void set_blades(std::span<const BladeGeometry*> blades);
    void set_frame(const Frame& frame);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  WingGeometry
// ========================================================================

class WingGeometry {
private:
    friend struct detail::ptr_accessor;
    friend class WingInflow;
    friend class VtkWing;
    friend class VtkWingWake;
    friend class WingLiftSurf;

    void* ptr_ = nullptr;

    explicit WingGeometry(void*);

public:
    WingGeometry() = default;

    WingGeometry(size_t num_parts, Vec3 origin, double wing_span);

    ~WingGeometry();

    WingGeometry(const WingGeometry&) = delete;
    WingGeometry& operator=(const WingGeometry&) = delete;
    WingGeometry(WingGeometry&&) noexcept;
    WingGeometry& operator=(WingGeometry&&) noexcept;

    void set_ctrl_points(size_t spanwise_nodes, size_t chordwise_nodes, double camber);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  WingPartGeometry (non-owning — owned by WingGeometry)
// ========================================================================

class WingPartGeometry {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit WingPartGeometry(void*);

public:
    WingPartGeometry() = default;

    WingPartGeometry(const WingPartGeometry&) = default;
    WingPartGeometry& operator=(const WingPartGeometry&) = default;

    void set_chord(const std::vector<double>& data);
    void set_twist(const std::vector<double>& data);
    void set_sweep(const std::vector<double>& data);
    void set_y_span(const std::vector<double>& data);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  Inflow (base)
// ========================================================================

class Inflow {
protected:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;
    bool owned_ = false;

    explicit Inflow(void*, bool owned = true);

public:
    Inflow() = default;
    virtual ~Inflow();

    Inflow(const Inflow&) = delete;
    Inflow& operator=(const Inflow&) = delete;
    Inflow(Inflow&&) noexcept;
    Inflow& operator=(Inflow&&) noexcept;

    double wake_skew() const;
    Frame frame();
    const Mat4* inverse_global_frame() const;

    void update(const AircraftState& ac_state, const Wake& wake, double dt);

    std::vector<double> inflow_at(const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& z) const;

    void update_wing_circulation(WingState& wing_state);
    void update_wing_dC_L(WingState& wing_state);

    InducedVelocities compute_wing_induced_vel_on_blade(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z) const;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  HuangPetersInflow
// ========================================================================

class HuangPetersInflow : public Inflow {
public:
    HuangPetersInflow() = default;

    HuangPetersInflow(long _Mo, long _Me,
                      const RotorGeometry& rotor,
                      const RotorInputState& rotor_input,
                      double dt);

    ~HuangPetersInflow() override;

    HuangPetersInflow(const HuangPetersInflow&) = delete;
    HuangPetersInflow& operator=(const HuangPetersInflow&) = delete;
    HuangPetersInflow(HuangPetersInflow&&) noexcept;
    HuangPetersInflow& operator=(HuangPetersInflow&&) noexcept;
};

// ========================================================================
//  NullInflow
// ========================================================================

class NullInflow : public Inflow {
public:
    NullInflow() = default;

    NullInflow(const RotorGeometry& rotor,
               const RotorInputState& rotor_input);

    ~NullInflow() override;

    NullInflow(const NullInflow&) = delete;
    NullInflow& operator=(const NullInflow&) = delete;
    NullInflow(NullInflow&&) noexcept;
    NullInflow& operator=(NullInflow&&) noexcept;
};

// ========================================================================
//  WingInflow
// ========================================================================

class WingInflow : public Inflow {
public:
    WingInflow() = default;

    WingInflow(const WingGeometry& wing,
               const WingInputState& wing_inputs,
               const WingLiftSurf& wing_lift_surf);

    ~WingInflow() override;

    WingInflow(const WingInflow&) = delete;
    WingInflow& operator=(const WingInflow&) = delete;
    WingInflow(WingInflow&&) noexcept;
    WingInflow& operator=(WingInflow&&) noexcept;
};

// ========================================================================
//  RotorInputState (non-owning — owned by AircraftInputState)
// ========================================================================

class RotorInputState {
private:
    friend struct detail::ptr_accessor;
    friend class HuangPetersInflow;
    friend class NullInflow;
    friend class AircraftInputState;  // calls constructor from get_rotor_input

    void* ptr_ = nullptr;

    explicit RotorInputState(void*);

public:
    RotorInputState() = default;

    RotorInputState(const RotorInputState&) = default;
    RotorInputState& operator=(const RotorInputState&) = default;

    void set_angular_velocity(double omega);
    double angular_velocity() const;

    void set_angular_accel(double alpha);
    double angular_accel() const;

    void set_azimuth(double azimuth);
    double azimuth() const;

    void set_r_0(const std::vector<double>& r_0);
    std::vector<double> get_r_0(size_t len) const;

    void set_blade_flapping(const std::vector<double>& flapping);
    std::vector<double> get_blade_flapping(size_t len) const;

    void set_blade_flapping_rate(const std::vector<double>& rate);
    std::vector<double> get_blade_flapping_rate(size_t len) const;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  AircraftInputState
// ========================================================================

class AircraftInputState {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit AircraftInputState(void*);

public:
    AircraftInputState() = default;

    AircraftInputState(size_t num_rotors,
                      const std::vector<size_t>& num_blades,
                      size_t num_wings);

    ~AircraftInputState();

    AircraftInputState(const AircraftInputState&) = delete;
    AircraftInputState& operator=(const AircraftInputState&) = delete;
    AircraftInputState(AircraftInputState&&) noexcept;
    AircraftInputState& operator=(AircraftInputState&&) noexcept;

    RotorInputState get_rotor_input(size_t rotor_idx);
    void set_blade_pitch(size_t rotor_idx, size_t blade_idx, double pitch);
    double get_blade_pitch(size_t rotor_idx, size_t blade_idx) const;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  RotorState (non-owning — owned by AircraftState)
// ========================================================================

class RotorState {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit RotorState(void*);

public:
    RotorState() = default;

    RotorState(const RotorState&) = default;
    RotorState& operator=(const RotorState&) = default;

    double get_C_T() const;
    void set_C_T(double C_T);

    double get_C_Q() const;
    void set_C_Q(double C_Q);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  BladeState (non-owning — owned by AircraftState)
// ========================================================================

class BladeState {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit BladeState(void*);

public:
    BladeState() = default;

    BladeState(const BladeState&) = default;
    BladeState& operator=(const BladeState&) = default;

    double azimuth() const;
    double C_T() const;
    double C_Q() const;
    double C_L() const;
    double C_D() const;

    std::vector<double> dC_T(size_t len) const;
    std::vector<double> dC_Db(size_t len) const;
    std::vector<double> dC_Db_profile(size_t len) const;
    std::vector<double> dC_Db_induced(size_t len) const;
    std::vector<double> dynamic_dC_Db_profile(size_t len) const;
    std::vector<double> dynamic_dC_Db_induced(size_t len) const;
    std::vector<double> dC_N(size_t len) const;
    std::vector<double> dC_c(size_t len) const;
    std::vector<double> dC_D(size_t len) const;
    std::vector<double> dC_T_dot(size_t len) const;
    std::vector<double> dC_Q(size_t len) const;
    std::vector<double> dC_L(size_t len) const;
    std::vector<double> dC_l(size_t len) const;
    std::vector<double> dC_Mz(size_t len) const;
    std::vector<double> dC_My(size_t len) const;
    std::vector<double> u_p(size_t len) const;
    std::vector<double> dynamic_u_p(size_t len) const;
    std::vector<double> u_t(size_t len) const;
    std::vector<double> aoa(size_t len) const;
    std::vector<double> aoa_eff(size_t len) const;
    std::vector<double> gamma(size_t len) const;
    std::vector<double> r_c(size_t len) const;
    std::vector<double> x(size_t len) const;
    std::vector<double> y(size_t len) const;
    std::vector<double> z(size_t len) const;

    std::vector<float> dC_Df(size_t len) const;
    std::vector<float> dC_Nf(size_t len) const;
    std::vector<float> dC_cf(size_t len) const;
    std::vector<float> dC_Tf(size_t len) const;
    std::vector<float> dC_Qf(size_t len) const;

    void fill_dC_T(double* data, size_t len) const;
    void fill_dC_Db(double* data, size_t len) const;
    void fill_dC_N(double* data, size_t len) const;
    void fill_dC_D(double* data, size_t len) const;
    void fill_dC_L(double* data, size_t len) const;
    void fill_dC_Q(double* data, size_t len) const;
    void fill_u_p(double* data, size_t len) const;
    void fill_u_t(double* data, size_t len) const;
    void fill_aoa(double* data, size_t len) const;
    void fill_gamma(double* data, size_t len) const;
    void fill_x(double* data, size_t len) const;
    void fill_y(double* data, size_t len) const;
    void fill_z(double* data, size_t len) const;
    void fill_r_c(double* data, size_t len) const;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  AircraftState
// ========================================================================

class AircraftState {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit AircraftState(void*);

public:
    AircraftState() = default;

    AircraftState(size_t num_rotors,
                  const std::vector<size_t>& num_blades,
                  size_t num_elements,
                  size_t num_wings,
                  const std::vector<size_t>& num_wing_parts,
                  size_t num_span_nodes,
                  size_t num_chord_nodes,
                  const Aircraft& aircraft,
                  std::span<Inflow*> rotor_inflows,
                  std::span<Inflow*> wing_inflows,
                  Direction direction);

    ~AircraftState();

    AircraftState(const AircraftState&) = delete;
    AircraftState& operator=(const AircraftState&) = delete;
    AircraftState(AircraftState&&) noexcept;
    AircraftState& operator=(AircraftState&&) noexcept;

    void set_freestream(const Vec4& freestream);
    Vec4 get_freestream() const;

    double rotor_C_T(size_t rotor_idx);
    double rotor_C_Q(size_t rotor_idx);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  VortexFilament (non-owning)
// ========================================================================

class VortexFilament {
private:
    friend struct detail::ptr_accessor;
    friend class RotorWake;  // calls constructor from get_tip_vortex

    void* ptr_ = nullptr;

    explicit VortexFilament(void*);

public:
    VortexFilament() = default;

    VortexFilament(const VortexFilament&) = default;
    VortexFilament& operator=(const VortexFilament&) = default;

    std::vector<double> x(size_t len) const;
    std::vector<double> y(size_t len) const;
    std::vector<double> z(size_t len) const;
    std::vector<double> gamma(size_t len) const;
    std::vector<double> r_c(size_t len) const;
    std::vector<double> v_z(size_t len) const;

    void fill_x(double* data, size_t len) const;
    void fill_y(double* data, size_t len) const;
    void fill_z(double* data, size_t len) const;
    void fill_gamma(double* data, size_t len) const;
    void fill_r_c(double* data, size_t len) const;
    void fill_v_z(double* data, size_t len) const;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  RotorWake (non-owning)
// ========================================================================

class RotorWake {
private:
    friend struct detail::ptr_accessor;
    friend class Wake;  // calls constructor from get_rotor_wake

    void* ptr_ = nullptr;

    explicit RotorWake(void*);

public:
    RotorWake() = default;

    RotorWake(const RotorWake&) = default;
    RotorWake& operator=(const RotorWake&) = default;

    VortexFilament get_tip_vortex(size_t blade_idx);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  Wake
// ========================================================================

class Wake {
protected:
    friend struct detail::ptr_accessor;
    friend class WakeHistory;  // WakeHistory::get_wake creates Wake objects

    void* ptr_ = nullptr;
    bool owned_ = false;

    explicit Wake(void*, bool owned = true);

public:
    Wake(size_t num_rotors, size_t num_blades,
         size_t wake_history, size_t radial_elements,
         const std::vector<size_t>& shed_history,
         const std::vector<size_t>& shed_release);

    virtual ~Wake();

    Wake(const Wake&) = delete;
    Wake& operator=(const Wake&) = delete;
    Wake(Wake&&) noexcept;
    Wake& operator=(Wake&&) noexcept;

    RotorWake get_rotor_wake(size_t rotor_idx);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  WakeHistory
// ========================================================================

class WakeHistory {
private:
    friend struct detail::ptr_accessor;
    friend class Wake;  // WakeHistory::get_wake creates Wake objects

    void* ptr_ = nullptr;

    explicit WakeHistory(void*);

public:
    WakeHistory(size_t num_rotors, size_t num_blades,
                size_t wake_history, size_t time_history, size_t radial_elements,
                const std::vector<size_t>& shed_history,
                const std::vector<size_t>& shed_release,
                double a1, bool hybrid);

    ~WakeHistory();

    WakeHistory(const WakeHistory&) = delete;
    WakeHistory& operator=(const WakeHistory&) = delete;
    WakeHistory(WakeHistory&&) noexcept;
    WakeHistory& operator=(WakeHistory&&) noexcept;

    void push_back();

    Wake get_wake(size_t idx);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  WingInputState (non-owning — owned by AircraftInputState)
// ========================================================================

class WingInputState {
private:
    friend struct detail::ptr_accessor;
    friend class WingInflow;

    void* ptr_ = nullptr;

    explicit WingInputState(void*);

public:
    WingInputState() = default;

    WingInputState(const WingInputState&) = default;
    WingInputState& operator=(const WingInputState&) = default;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  WingState (non-owning)
// ========================================================================

class WingState {
private:
    friend struct detail::ptr_accessor;
    friend class Inflow;

    void* ptr_ = nullptr;

    explicit WingState(void*);

public:
    WingState() = default;

    WingState(const WingState&) = default;
    WingState& operator=(const WingState&) = default;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  WingLiftSurf (non-owning — lifetime managed externally)
// ========================================================================

class WingLiftSurf {
private:
    friend struct detail::ptr_accessor;
    friend class WingInflow;

    void* ptr_ = nullptr;

    explicit WingLiftSurf(void*);

public:
    WingLiftSurf() = default;

    WingLiftSurf(const WingLiftSurf&) = default;
    WingLiftSurf& operator=(const WingLiftSurf&) = default;

    void set_vortex_geometry(const WingGeometry& wing,
                             size_t spanwise_chunks, size_t chordwise_nodes);

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  VTK Output Types
// ========================================================================

class VtkRotor {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit VtkRotor(void*);

public:
    VtkRotor() = default;

    static VtkRotor build(const RotorGeometry& rotor);

    ~VtkRotor();

    VtkRotor(const VtkRotor&) = delete;
    VtkRotor& operator=(const VtkRotor&) = delete;
    VtkRotor(VtkRotor&&) noexcept;
    VtkRotor& operator=(VtkRotor&&) noexcept;

    operator bool() const { return ptr_ != nullptr; }
};

class VtkWing {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit VtkWing(void*);

public:
    VtkWing() = default;

    static VtkWing build(const WingGeometry& wing);

    ~VtkWing();

    VtkWing(const VtkWing&) = delete;
    VtkWing& operator=(const VtkWing&) = delete;
    VtkWing(VtkWing&&) noexcept;
    VtkWing& operator=(VtkWing&&) noexcept;

    operator bool() const { return ptr_ != nullptr; }
};

class VtkWake {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit VtkWake(void*);

public:
    VtkWake() = default;

    static VtkWake build(const Wake& wake);

    ~VtkWake();

    VtkWake(const VtkWake&) = delete;
    VtkWake& operator=(const VtkWake&) = delete;
    VtkWake(VtkWake&&) noexcept;
    VtkWake& operator=(VtkWake&&) noexcept;

    operator bool() const { return ptr_ != nullptr; }
};

class VtkWingWake {
private:
    friend struct detail::ptr_accessor;

    void* ptr_ = nullptr;

    explicit VtkWingWake(void*);

public:
    VtkWingWake() = default;

    static VtkWingWake build(const WingGeometry& wing, const WingLiftSurf& lift_surf);

    ~VtkWingWake();

    VtkWingWake(const VtkWingWake&) = delete;
    VtkWingWake& operator=(const VtkWingWake&) = delete;
    VtkWingWake(VtkWingWake&&) noexcept;
    VtkWingWake& operator=(VtkWingWake&&) noexcept;

    operator bool() const { return ptr_ != nullptr; }
};

// ========================================================================
//  VTK Write Functions
// ========================================================================

void write_rotor_vtu(std::string_view filename, size_t step, size_t iteration,
                     const VtkRotor& vtk, const RotorState& state,
                     const RotorGeometry& geom);

void write_rotors_vtu(std::string_view filename, size_t step,
                      std::span<const VtkRotor*> vtks,
                      const AircraftState& ac_state, const Aircraft& aircraft);

void write_wing_vtu(std::string_view filename, size_t step, size_t iteration,
                    const VtkWing& vtk, const WingState& state,
                    const WingGeometry& geom);

void write_wake_vtu(std::string_view filename, size_t step,
                    const VtkWake& vtk, const Wake& wake);

void write_wing_wake_vtu(std::string_view filename, size_t step, size_t iteration,
                         const VtkWingWake& vtk, const WingGeometry& wing,
                         const WingLiftSurf& lift_surf, const WingInputState& input);

void write_wake_field_vtu(std::string_view filename,
                          const AircraftState& ac_state, const Wake& wake,
                          double x_min, double x_max,
                          double y_min, double y_max,
                          double z_min, double z_max,
                          size_t nx, size_t ny, size_t nz);

// ========================================================================
//  Utility / Free Functions
// ========================================================================

size_t chunk_size();
Direction direction_clockwise();
Direction direction_counter_clockwise();
Mat3 mat3_identity();
Mat4 mat4_identity();
std::vector<double> generate_radius_points(size_t n_sections, double root_cutout);

void simulation_step(const AircraftState& ac_state, const Aircraft& aircraft,
                     const AircraftInputState& ac_input_state,
                     const WakeHistory& wake_history,
                     const AtmosphereData& atmo, size_t iteration, double dt,
                     bool track_bwi_events, bool converged);

void basic_aircraft_rotor_dynamics(AircraftInputState& input, double dt);
double basic_single_rotor_dynamics(RotorInputState& input, double dt);

inline void set_wing_ctrl_pt_geometry(WingGeometry& wing,
                                     size_t spanwise_nodes,
                                     size_t chordwise_nodes,
                                     double camber) {
    wing.set_ctrl_points(spanwise_nodes, chordwise_nodes, camber);
}

void test_func();

// ========================================================================
//  ptr_accessor definitions (after all classes are complete)
// ========================================================================

namespace detail {
    struct ptr_accessor {
        static void* get(Frame& x) { return x.ptr_; }
        static void* get(Aircraft& x) { return x.ptr_; }
        static void* get(RotorGeometry& x) { return x.ptr_; }
        static void* get(BladeGeometry& x) { return x.ptr_; }
        static void* get(BladeAirfoil& x) { return x.ptr_; }
        static void* get(AirfoilModel& x) { return x.ptr_; }
        static void* get(WingGeometry& x) { return x.ptr_; }
        static void* get(Inflow& x) { return x.ptr_; }
        static void* get(AircraftInputState& x) { return x.ptr_; }
        static void* get(RotorInputState& x) { return x.ptr_; }
        static void* get(AircraftState& x) { return x.ptr_; }
        static void* get(RotorState& x) { return x.ptr_; }
        static void* get(BladeState& x) { return x.ptr_; }
        static void* get(Wake& x) { return x.ptr_; }
        static void* get(WakeHistory& x) { return x.ptr_; }
        static void* get(RotorWake& x) { return x.ptr_; }
        static void* get(VortexFilament& x) { return x.ptr_; }
        static void* get(WingInputState& x) { return x.ptr_; }
        static void* get(WingState& x) { return x.ptr_; }
        static void* get(WingLiftSurf& x) { return x.ptr_; }
        static void* get(VtkRotor& x) { return x.ptr_; }
        static void* get(VtkWing& x) { return x.ptr_; }
        static void* get(VtkWake& x) { return x.ptr_; }
        static void* get(VtkWingWake& x) { return x.ptr_; }

        static bool owned(const Frame& x) { return x.owned_; }
        static bool owned(const Inflow& x) { return x.owned_; }
        static bool owned(const Wake& x) { return x.owned_; }
    };

    // Allow all cppbindings.cpp code to construct objects via ptr_accessor friendship
    namespace bindings_impl { struct access; }
} // namespace detail

} // namespace opencopter

#endif /* OPENCOPTER_HPP_ */
