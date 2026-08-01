/**
 * Implementation of the idiomatic C++ wrapper for OpenCOPTER.
 *
 * This file includes opencopter.h (C API) and implements all C++ wrapper
 * methods. The C types (OC_*) are entirely contained here.
 */

#include "opencopter.hpp"
#include "opencopter.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>

namespace opencopter {

// ========================================================================
//  Internal helpers: convert between native C++ value types and OC_* types
// ========================================================================

static OC_Vec3 to_oc(Vec3 v) { return OC_Vec3{v.x, v.y, v.z}; }
static Vec3 from_oc(OC_Vec3 v) { return Vec3{v.x, v.y, v.z}; }
static OC_Vec4 to_oc(Vec4 v) { return OC_Vec4{v.x, v.y, v.z, v.w}; }
static Vec4 from_oc(OC_Vec4 v) { return Vec4{v.x, v.y, v.z, v.w}; }

static Mat3 from_oc(OC_Mat3 m) { Mat3 r{}; std::memcpy(r.data, m.data, sizeof(r.data)); return r; }
static Mat4 from_oc(OC_Mat4 m) { Mat4 r{}; std::memcpy(r.data, m.data, sizeof(r.data)); return r; }

namespace detail2 {
// cast helpers: void* -> OC_* (each gets a unique function name)
#define CAST(type) static type* cast_##type(void* p){ return static_cast<type*>(p); }
CAST(OC_Frame)
CAST(OC_Aircraft)
CAST(OC_RotorGeometry)
CAST(OC_BladeGeometry)
CAST(OC_BladeAirfoil)
CAST(OC_AirfoilModel)
CAST(OC_WingGeometry)
CAST(OC_Inflow)
CAST(OC_AircraftInputState)
CAST(OC_RotorInputState)
CAST(OC_AircraftState)
CAST(OC_RotorState)
CAST(OC_BladeState)
CAST(OC_Wake)
CAST(OC_WakeHistory)
CAST(OC_RotorWake)
CAST(OC_VortexFilament)
CAST(OC_WingInputState)
CAST(OC_WingState)
CAST(OC_WingLiftSurf)
CAST(OC_WingPartGeometry)
CAST(OC_VtkRotor)
CAST(OC_VtkWing)
CAST(OC_VtkWake)
CAST(OC_VtkWingWake)
#undef CAST

    // short aliases using ptr_accessor::get(obj) -> OC_*
    static inline OC_Frame* fp(Frame& x){ return cast_OC_Frame(detail::ptr_accessor::get(x)); }
    static inline OC_Frame* fp(const Frame& x){ return cast_OC_Frame(detail::ptr_accessor::get(const_cast<Frame&>(x))); }
    static inline OC_Aircraft* ap(Aircraft& x){ return cast_OC_Aircraft(detail::ptr_accessor::get(x)); }
    static inline OC_Aircraft* ap(const Aircraft& x){ return cast_OC_Aircraft(detail::ptr_accessor::get(const_cast<Aircraft&>(x))); }
    static inline OC_RotorGeometry* rgp(RotorGeometry& x){ return cast_OC_RotorGeometry(detail::ptr_accessor::get(x)); }
    static inline OC_RotorGeometry* rgp(const RotorGeometry& x){ return cast_OC_RotorGeometry(detail::ptr_accessor::get(const_cast<RotorGeometry&>(x))); }
    static inline OC_BladeGeometry* bgp(BladeGeometry& x){ return cast_OC_BladeGeometry(detail::ptr_accessor::get(x)); }
    static inline OC_BladeGeometry* bgp(const BladeGeometry& x){ return cast_OC_BladeGeometry(detail::ptr_accessor::get(const_cast<BladeGeometry&>(x))); }
    static inline OC_BladeAirfoil* bap(BladeAirfoil& x){ return cast_OC_BladeAirfoil(detail::ptr_accessor::get(x)); }
    static inline OC_BladeAirfoil* bap(const BladeAirfoil& x){ return cast_OC_BladeAirfoil(detail::ptr_accessor::get(const_cast<BladeAirfoil&>(x))); }
    static inline OC_AirfoilModel* amp(AirfoilModel& x){ return cast_OC_AirfoilModel(detail::ptr_accessor::get(x)); }
    static inline OC_AirfoilModel* amp(const AirfoilModel& x){ return cast_OC_AirfoilModel(detail::ptr_accessor::get(const_cast<AirfoilModel&>(x))); }
    static inline OC_WingGeometry* wg(WingGeometry& x){ return cast_OC_WingGeometry(detail::ptr_accessor::get(x)); }
    static inline OC_WingGeometry* wg(const WingGeometry& x){ return cast_OC_WingGeometry(detail::ptr_accessor::get(const_cast<WingGeometry&>(x))); }
    static inline OC_Inflow* inl(Inflow& x){ return cast_OC_Inflow(detail::ptr_accessor::get(x)); }
    static inline OC_Inflow* inl(const Inflow& x){ return cast_OC_Inflow(detail::ptr_accessor::get(const_cast<Inflow&>(x))); }
    static inline OC_AircraftInputState* ais(AircraftInputState& x){ return cast_OC_AircraftInputState(detail::ptr_accessor::get(x)); }
    static inline OC_AircraftInputState* ais(const AircraftInputState& x){ return cast_OC_AircraftInputState(detail::ptr_accessor::get(const_cast<AircraftInputState&>(x))); }
    static inline OC_RotorInputState* ris(RotorInputState& x){ return cast_OC_RotorInputState(detail::ptr_accessor::get(x)); }
    static inline OC_RotorInputState* ris(const RotorInputState& x){ return cast_OC_RotorInputState(detail::ptr_accessor::get(const_cast<RotorInputState&>(x))); }
    static inline OC_AircraftState* ast(AircraftState& x){ return cast_OC_AircraftState(detail::ptr_accessor::get(x)); }
    static inline OC_AircraftState* ast(const AircraftState& x){ return cast_OC_AircraftState(detail::ptr_accessor::get(const_cast<AircraftState&>(x))); }
    static inline OC_RotorState* rst(RotorState& x){ return cast_OC_RotorState(detail::ptr_accessor::get(x)); }
    static inline OC_RotorState* rst(const RotorState& x){ return cast_OC_RotorState(detail::ptr_accessor::get(const_cast<RotorState&>(x))); }
    static inline OC_BladeState* bst(BladeState& x){ return cast_OC_BladeState(detail::ptr_accessor::get(x)); }
    static inline OC_BladeState* bst(const BladeState& x){ return cast_OC_BladeState(detail::ptr_accessor::get(const_cast<BladeState&>(x))); }
    static inline OC_Wake* wk(Wake& x){ return cast_OC_Wake(detail::ptr_accessor::get(x)); }
    static inline OC_Wake* wk(const Wake& x){ return cast_OC_Wake(detail::ptr_accessor::get(const_cast<Wake&>(x))); }
    static inline OC_WakeHistory* wh(WakeHistory& x){ return cast_OC_WakeHistory(detail::ptr_accessor::get(x)); }
    static inline OC_WakeHistory* wh(const WakeHistory& x){ return cast_OC_WakeHistory(detail::ptr_accessor::get(const_cast<WakeHistory&>(x))); }
    static inline OC_RotorWake* rw(RotorWake& x){ return cast_OC_RotorWake(detail::ptr_accessor::get(x)); }
    static inline OC_RotorWake* rw(const RotorWake& x){ return cast_OC_RotorWake(detail::ptr_accessor::get(const_cast<RotorWake&>(x))); }
    static inline OC_VortexFilament* vf(VortexFilament& x){ return cast_OC_VortexFilament(detail::ptr_accessor::get(x)); }
    static inline OC_VortexFilament* vf(const VortexFilament& x){ return cast_OC_VortexFilament(detail::ptr_accessor::get(const_cast<VortexFilament&>(x))); }
    static inline OC_WingInputState* wis(WingInputState& x){ return cast_OC_WingInputState(detail::ptr_accessor::get(x)); }
    static inline OC_WingInputState* wis(const WingInputState& x){ return cast_OC_WingInputState(detail::ptr_accessor::get(const_cast<WingInputState&>(x))); }
    static inline OC_WingState* ws(WingState& x){ return cast_OC_WingState(detail::ptr_accessor::get(x)); }
    static inline OC_WingState* ws(const WingState& x){ return cast_OC_WingState(detail::ptr_accessor::get(const_cast<WingState&>(x))); }
    static inline OC_WingLiftSurf* wls(WingLiftSurf& x){ return cast_OC_WingLiftSurf(detail::ptr_accessor::get(x)); }
    static inline OC_WingLiftSurf* wls(const WingLiftSurf& x){ return cast_OC_WingLiftSurf(detail::ptr_accessor::get(const_cast<WingLiftSurf&>(x))); }
    static inline OC_VtkRotor* vr(VtkRotor& x){ return cast_OC_VtkRotor(detail::ptr_accessor::get(x)); }
    static inline OC_VtkRotor* vr(const VtkRotor& x){ return cast_OC_VtkRotor(detail::ptr_accessor::get(const_cast<VtkRotor&>(x))); }
    static inline OC_VtkWing* vw(VtkWing& x){ return cast_OC_VtkWing(detail::ptr_accessor::get(x)); }
    static inline OC_VtkWing* vw(const VtkWing& x){ return cast_OC_VtkWing(detail::ptr_accessor::get(const_cast<VtkWing&>(x))); }
    static inline OC_VtkWake* vkw(VtkWake& x){ return cast_OC_VtkWake(detail::ptr_accessor::get(x)); }
    static inline OC_VtkWake* vkw(const VtkWake& x){ return cast_OC_VtkWake(detail::ptr_accessor::get(const_cast<VtkWake&>(x))); }
    static inline OC_VtkWingWake* vww(VtkWingWake& x){ return cast_OC_VtkWingWake(detail::ptr_accessor::get(x)); }
    static inline OC_VtkWingWake* vww(const VtkWingWake& x){ return cast_OC_VtkWingWake(detail::ptr_accessor::get(const_cast<VtkWingWake&>(x))); }
} // namespace detail2

using namespace detail2;

// ========================================================================
//  Utility / Free Functions
// ========================================================================

size_t chunk_size() { return oc_chunk_size(); }
Direction direction_clockwise() { return static_cast<Direction>(oc_direction_clockwise()); }
Direction direction_counter_clockwise() { return static_cast<Direction>(oc_direction_counter_clockwise()); }
Mat3 mat3_identity() { return from_oc(oc_mat3_identity()); }
Mat4 mat4_identity() { return from_oc(oc_mat4_identity()); }

std::span<double> generate_radius_points(size_t n_sections, double root_cutout) {

    double* buff = oc_generate_radius_points(&n_sections, root_cutout);

    return std::span<double>(buff, n_sections);
}

void simulation_step(const AircraftState& ac_state, const Aircraft& aircraft,
                     const AircraftInputState& ac_input_state,
                     const WakeHistory& wake_history,
                     const AtmosphereData& atmo, size_t iteration, double dt,
                     bool track_bwi_events, bool converged) {
    OC_Atmosphere oca{atmo.density, atmo.dynamic_viscosity,
                      atmo.kinematic_viscosity, atmo.speed_of_sound};
    oc_simulation_step(ast(ac_state), ap(aircraft), ais(ac_input_state), wh(wake_history),
                       &oca, iteration, dt,
                       track_bwi_events ? 1 : 0, converged ? 1 : 0);
}

void basic_aircraft_rotor_dynamics(AircraftInputState& input, double dt) {
    oc_basic_aircraft_rotor_dynamics(ais(input), dt);
}

double basic_single_rotor_dynamics(RotorInputState& input, double dt) {
    return oc_basic_single_rotor_dynamics(ris(input), dt);
}

// ========================================================================
//  Frame
// ========================================================================

Frame::Frame(void* p, bool owned) : ptr_(p), owned_(owned) {}

Frame::Frame(Vec3 axis, double angle, Vec3 translation,
               const Frame* parent, std::string_view name, FrameType frame_type) {
    OC_Frame* raw = oc_frame_create(to_oc(axis), angle, to_oc(translation),
                                      parent ? fp(*parent) : nullptr,
                                      name.data(), static_cast<int>(frame_type));
    ptr_ = raw;
    // owned_ defaults to false since frames are typically part of a hierarchy
    // where the Aircraft owns the entire tree. Setting owned_=true would cause
    // GC.removeRoot to be called on a frame still referenced by its parent's
    // .children array, corrupting the D GC pinned list. Users who need explicit
    // cleanup can call oc_frame_destroy via the C API directly.
    owned_ = false;
}

Frame::~Frame() {
    if (owned_ && ptr_) { oc_frame_destroy(fp(*this)); ptr_ = nullptr; }
}

Frame::Frame(Frame&& other) noexcept : ptr_(other.ptr_), owned_(other.owned_) {
    other.ptr_ = nullptr; other.owned_ = false;
}

Frame& Frame::operator=(Frame&& other) noexcept {
    if (this != &other) {
        if (owned_ && ptr_) oc_frame_destroy(fp(*this));
        ptr_ = other.ptr_; owned_ = other.owned_;
        other.ptr_ = nullptr; other.owned_ = false;
    }
    return *this;
}

void Frame::set_rotation(Vec3 axis, double angle) { if (ptr_) oc_frame_set_rotation(fp(*this), to_oc(axis), angle); }
void Frame::rotate(Vec3 axis, double angle) { if (ptr_) oc_frame_rotate(fp(*this), to_oc(axis), angle); }
void Frame::translate(Vec3 translation) { if (ptr_) oc_frame_translate(fp(*this), to_oc(translation)); }

void Frame::update(const Mat4& parent_global_mat) {
    if (ptr_) { OC_Mat4 m{}; std::memcpy(m.data, parent_global_mat.data, sizeof(m.data));
        oc_frame_update(fp(*this), &m); }
}

void Frame::set_children(std::span<const Frame*> children) {
    if (!ptr_) return;
    size_t n = children.size();
    std::vector<OC_Frame*> raw(n);
    for (size_t i = 0; i < n; ++i) raw[i] = children[i] ? fp(*children[i]) : nullptr;
    oc_frame_set_children(fp(*this), raw.data(), n);
}

void Frame::set_frame_type(FrameType ft) { if (ptr_) oc_frame_set_frame_type(fp(*this), static_cast<int>(ft)); }
void Frame::set_name(std::string_view name) { if (ptr_) oc_frame_set_name(fp(*this), name.data()); }

const Mat4* Frame::local_matrix() const {
    if (!ptr_) return nullptr;
    OC_Mat4* r = const_cast<OC_Mat4*>(oc_frame_get_local_matrix(fp(*this)));
    return reinterpret_cast<const Mat4*>(r);
}
const Mat4* Frame::global_matrix() const {
    if (!ptr_) return nullptr;
    OC_Mat4* r = const_cast<OC_Mat4*>(oc_frame_get_global_matrix(fp(*this)));
    return reinterpret_cast<const Mat4*>(r);
}
const Mat4* Frame::inverse_global_matrix() const {
    if (!ptr_) return nullptr;
    OC_Mat4* r = const_cast<OC_Mat4*>(oc_frame_get_inverse_global_matrix(fp(*this)));
    return reinterpret_cast<const Mat4*>(r);
}

Frame Frame::parent() {
    OC_Frame* p = ptr_ ? oc_frame_get_parent(fp(*this)) : nullptr;
    return Frame(p, false);
}

// ========================================================================
//  Aircraft
// ========================================================================

Aircraft::Aircraft(void* p) : ptr_(p) {}
Aircraft::Aircraft(size_t num_rotors, size_t num_wings) { ptr_ = oc_aircraft_create(num_rotors, num_wings); }
Aircraft::~Aircraft() { if (ptr_) oc_aircraft_destroy(ap(*this)); }
Aircraft::Aircraft(Aircraft&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
Aircraft& Aircraft::operator=(Aircraft&& o) noexcept {
    if (this != &o) { if (ptr_) oc_aircraft_destroy(ap(*this)); ptr_ = o.ptr_; o.ptr_ = nullptr; }
    return *this;
}

Frame Aircraft::root_frame() {
    OC_Frame* f = ptr_ ? oc_aircraft_get_root_frame(ap(*this)) : nullptr;
    return Frame(f, false);
}

void Aircraft::set_rotors(std::span<const RotorGeometry*> rotors) {
    if (!ptr_) return;
    size_t n = rotors.size();
    std::vector<OC_RotorGeometry*> raw(n);
    for (size_t i = 0; i < n; ++i) raw[i] = rotors[i] ? rgp(*rotors[i]) : nullptr;
    oc_aircraft_set_rotors(ap(*this), raw.data(), n);
}

// ========================================================================
//  RotorGeometry
// ========================================================================

RotorGeometry::RotorGeometry(void* p) : ptr_(p) {}
RotorGeometry::RotorGeometry(size_t nb, Vec3 origin, double radius, double solidity)
    { ptr_ = oc_rotor_geometry_create(nb, to_oc(origin), radius, solidity); }
RotorGeometry::~RotorGeometry() { if (ptr_) oc_rotor_geometry_destroy(rgp(*this)); }
RotorGeometry::RotorGeometry(RotorGeometry&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
RotorGeometry& RotorGeometry::operator=(RotorGeometry&& o) noexcept {
    if (this != &o) { if (ptr_) oc_rotor_geometry_destroy(rgp(*this)); ptr_ = o.ptr_; o.ptr_ = nullptr; }
    return *this;
}

void RotorGeometry::set_solidity(double s) { if (ptr_) oc_rotor_geometry_set_solidity(rgp(*this), s); }

void RotorGeometry::set_blades(std::span<const BladeGeometry*> blades) {
    if (!ptr_) return; size_t n = blades.size();
    std::vector<OC_BladeGeometry*> raw(n);
    for (size_t i = 0; i < n; ++i) raw[i] = blades[i] ? bgp(*blades[i]) : nullptr;
    oc_rotor_geometry_set_blades(rgp(*this), raw.data(), n);
}

void RotorGeometry::set_frame(const Frame& frame) { if (ptr_) oc_rotor_geometry_set_frame(rgp(*this), fp(frame)); }

// ========================================================================
//  BladeGeometry
// ========================================================================

BladeGeometry::BladeGeometry(void* p) : ptr_(p) {}
BladeGeometry::BladeGeometry(size_t ne, double ao, double ac, const BladeAirfoil& af, double rc) {
    ptr_ = oc_blade_geometry_create(ne, ao, ac, bap(af), rc);
}
BladeGeometry::~BladeGeometry() { if (ptr_) oc_blade_geometry_destroy(bgp(*this)); }
BladeGeometry::BladeGeometry(BladeGeometry&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
BladeGeometry& BladeGeometry::operator=(BladeGeometry&& o) noexcept {
    if (this != &o) { if (ptr_) oc_blade_geometry_destroy(bgp(*this)); ptr_ = o.ptr_; o.ptr_ = nullptr; }
    return *this;
}

#define BG_SET(name) void BladeGeometry::set_##name(const std::vector<double>& d){ \
    if(ptr_) oc_blade_geometry_set_##name(bgp(*this), const_cast<double*>(d.data()), d.size()); } \
    void BladeGeometry::set_##name(const std::span<double>& d){ \
    if(ptr_) oc_blade_geometry_set_##name(bgp(*this), const_cast<double*>(d.data()), d.size()); }
    
BG_SET(twist); BG_SET(chord); BG_SET(radius); BG_SET(C_l_alpha); BG_SET(alpha_0);
BG_SET(sweep); BG_SET(xi); BG_SET(thickness); BG_SET(xi_p);
#undef BG_SET

void BladeGeometry::compute_vectors() { if (ptr_) oc_blade_geometry_compute_vectors(bgp(*this)); }

Frame BladeGeometry::get_frame() const {
    OC_Frame* f = ptr_ ? oc_blade_geometry_get_frame(bgp(*this)) : nullptr;
    return Frame(f, false);
}
void BladeGeometry::set_frame(const Frame& frame) { if (ptr_) oc_blade_geometry_set_frame(bgp(*this), fp(frame)); }
void BladeGeometry::set_blade_length(double l) { if (ptr_) oc_blade_geometry_set_blade_length(bgp(*this), l); }

// ========================================================================
//  WingGeometry
// ========================================================================

WingGeometry::WingGeometry(void* p) : ptr_(p) {}
WingGeometry::WingGeometry(size_t np, Vec3 origin, double ws)
    { ptr_ = oc_wing_geometry_create(np, to_oc(origin), ws); }
WingGeometry::~WingGeometry() { if (ptr_) oc_wing_geometry_destroy(wg(*this)); }
WingGeometry::WingGeometry(WingGeometry&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
WingGeometry& WingGeometry::operator=(WingGeometry&& o) noexcept {
    if (this != &o) { if (ptr_) oc_wing_geometry_destroy(wg(*this)); ptr_ = o.ptr_; o.ptr_ = nullptr; }
    return *this;
}

void WingGeometry::set_ctrl_points(size_t sn, size_t cn, double camber) {
    if (ptr_) oc_wing_geometry_set_ctrl_points(wg(*this), sn, cn, camber);
}

// ========================================================================
//  WingPartGeometry
// ========================================================================

WingPartGeometry::WingPartGeometry(void* p) : ptr_(p) {}
void WingPartGeometry::set_chord(const std::vector<double>& d){ if(ptr_) oc_wing_part_geometry_set_chord(cast_OC_WingPartGeometry(ptr_), const_cast<double*>(d.data()), d.size()); }
void WingPartGeometry::set_twist(const std::vector<double>& d){ if(ptr_) oc_wing_part_geometry_set_twist(cast_OC_WingPartGeometry(ptr_), const_cast<double*>(d.data()), d.size()); }
void WingPartGeometry::set_sweep(const std::vector<double>& d){ if(ptr_) oc_wing_part_geometry_set_sweep(cast_OC_WingPartGeometry(ptr_), const_cast<double*>(d.data()), d.size()); }
void WingPartGeometry::set_y_span(const std::vector<double>& d){ if(ptr_) oc_wing_part_geometry_set_y_span(cast_OC_WingPartGeometry(ptr_), const_cast<double*>(d.data()), d.size()); }

// ========================================================================
//  AirfoilModel
// ========================================================================

AirfoilModel::AirfoilModel(void* p) : ptr_(p) {}
AirfoilModel::~AirfoilModel() { if (ptr_) oc_airfoil_model_destroy(amp(*this)); }
AirfoilModel::AirfoilModel(AirfoilModel&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
AirfoilModel& AirfoilModel::operator=(AirfoilModel&& o) noexcept {
    if (this != &o) { if (ptr_) oc_airfoil_model_destroy(amp(*this)); ptr_ = o.ptr_; o.ptr_ = nullptr; }
    return *this;
}

AirfoilModel AirfoilModel::thin_airfoil(double c) { return AirfoilModel(oc_thin_airfoil_create(c)); }

AirfoilModel AirfoilModel::aero_das(const std::vector<double>& a, const std::vector<double>& cl,
    const std::vector<double>& cd, double tbyc, double ar) {
    return AirfoilModel(oc_aero_das_create(const_cast<double*>(a.data()), a.size(),
        const_cast<double*>(cl.data()), cl.size(), const_cast<double*>(cd.data()), cd.size(), tbyc, ar));
}

AirfoilModel AirfoilModel::aero_das_from_xfoil_polar(std::string_view fn, double tbyc) {
    return AirfoilModel(oc_aero_das_from_xfoil_polar(fn.data(), tbyc));
}
AirfoilModel AirfoilModel::c81_from_file(std::string_view fn) {
    return AirfoilModel(oc_c81_from_file(fn.data()));
}

double AirfoilModel::get_Cl(double a, double m) const { return ptr_ ? oc_airfoil_get_Cl(amp(*this), a, m) : 0.0; }
double AirfoilModel::get_Cd(double a, double m) const { return ptr_ ? oc_airfoil_get_Cd(amp(*this), a, m) : 0.0; }
double AirfoilModel::lift_curve_slope() const { return ptr_ ? oc_airfoil_lift_curve_slope(amp(*this)) : 0.0; }
double AirfoilModel::zero_lift_aoa() const { return ptr_ ? oc_airfoil_zero_lift_aoa(amp(*this)) : 0.0; }

// ========================================================================
//  BladeAirfoil
// ========================================================================

BladeAirfoil::BladeAirfoil(void* p) : ptr_(p) {}
BladeAirfoil::~BladeAirfoil() { if (ptr_) oc_blade_airfoil_destroy(bap(*this)); }
BladeAirfoil::BladeAirfoil(BladeAirfoil&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
BladeAirfoil& BladeAirfoil::operator=(BladeAirfoil&& o) noexcept {
    if (this != &o) { if (ptr_) oc_blade_airfoil_destroy(bap(*this)); ptr_ = o.ptr_; o.ptr_ = nullptr; }
    return *this;
}

BladeAirfoil BladeAirfoil::create_basic(size_t ne, double c) {
    return BladeAirfoil(oc_blade_airfoil_create_basic(ne, c));
}

BladeAirfoil BladeAirfoil::create(const std::vector<AirfoilModel>& models, const std::vector<size_t>& extents) {
    size_t n = models.size();
    std::vector<OC_AirfoilModel*> raw(n);
    for (size_t i = 0; i < n; ++i) raw[i] = amp(models[i]);
    return BladeAirfoil(oc_blade_airfoil_create(raw.data(), extents.data(), n));
}

double BladeAirfoil::get_Cl(size_t ci, double a, double m) const { return ptr_ ? oc_blade_airfoil_get_Cl(bap(*this), ci, a, m) : 0.0; }
double BladeAirfoil::get_Cd(size_t ci, double a, double m) const { return ptr_ ? oc_blade_airfoil_get_Cd(bap(*this), ci, a, m) : 0.0; }
double BladeAirfoil::lift_curve_slope(size_t ci) const { return ptr_ ? oc_blade_airfoil_lift_curve_slope(bap(*this), ci) : 0.0; }
double BladeAirfoil::zero_lift_aoa(size_t ci) const { return ptr_ ? oc_blade_airfoil_zero_lift_aoa(bap(*this), ci) : 0.0; }

std::vector<double> BladeAirfoil::fill_lift_curve_slope(size_t ci) const {
    std::vector<double> r(8);
    if (ptr_) oc_blade_airfoil_fill_lift_curve_slope(bap(*this), ci, r.data(), r.size());
    return r;
}

std::vector<double> BladeAirfoil::fill_zero_lift_aoa(size_t ci) const {
    std::vector<double> r(8);
    if (ptr_) oc_blade_airfoil_fill_zero_lift_aoa(bap(*this), ci, r.data(), r.size());
    return r;
}

void BladeAirfoil::fill_coefficients(size_t ci, const std::vector<double>& alphas,
    const std::vector<double>& machs, std::vector<double>& Cl_out, std::vector<double>& Cd_out) const {
    if (!ptr_ || alphas.size() != machs.size()) return;
    size_t len = alphas.size();
    if (Cl_out.size() < len) Cl_out.resize(len);
    if (Cd_out.size() < len) Cd_out.resize(len);
    oc_blade_airfoil_fill_coefficients(bap(*this), ci, alphas.data(), machs.data(), Cl_out.data(), Cd_out.data(), len);
}

// ========================================================================
//  Inflow
// ========================================================================

Inflow::Inflow(void* p, bool owned) : ptr_(p), owned_(owned) {}
Inflow::~Inflow() { if (owned_ && ptr_) { oc_inflow_destroy(inl(*this)); ptr_ = nullptr; } }
Inflow::Inflow(Inflow&& o) noexcept : ptr_(o.ptr_), owned_(o.owned_) { o.ptr_ = nullptr; o.owned_ = false; }
Inflow& Inflow::operator=(Inflow&& o) noexcept {
    if (this != &o) {
        if (owned_ && ptr_) oc_inflow_destroy(inl(*this));
        ptr_ = o.ptr_; owned_ = o.owned_; o.ptr_ = nullptr; o.owned_ = false;
    }
    return *this;
}

double Inflow::wake_skew() const { return ptr_ ? oc_inflow_wake_skew(inl(*this)) : 0.0; }

Frame Inflow::frame() {
    OC_Frame* f = ptr_ ? oc_inflow_get_frame(inl(*this)) : nullptr;
    return Frame(f, false);
}

const Mat4* Inflow::inverse_global_frame() const {
    if (!ptr_) return nullptr;
    OC_Mat4* r = const_cast<OC_Mat4*>(oc_inflow_get_inverse_global_frame(inl(*this)));
    return reinterpret_cast<const Mat4*>(r);
}

void Inflow::update(const AircraftState& ac_state, const Wake& wake, double dt) {
    if (ptr_) oc_inflow_update(inl(*this), ast(ac_state), wk(wake), dt);
}

std::vector<double> Inflow::inflow_at(const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<double>& z) const {
    if (!ptr_ || x.size() != y.size() || y.size() != z.size()) return {};
    size_t len = x.size(); std::vector<double> result(len * 3);
    oc_inflow_at(inl(*this), x.data(), y.data(), z.data(), result.data(), len);
    return result;
}

void Inflow::update_wing_circulation(WingState& wing_state) {
    if (ptr_) oc_inflow_update_wing_circulation(inl(*this), ws(wing_state));
}
void Inflow::update_wing_dC_L(WingState& wing_state) {
    if (ptr_) oc_inflow_update_wing_dC_L(inl(*this), ws(wing_state));
}

InducedVelocities Inflow::compute_wing_induced_vel_on_blade(const std::vector<double>& x,
    const std::vector<double>& y, const std::vector<double>& z) const {
    InducedVelocities result{};
    if (ptr_ && x.size() >= 8 && y.size() >= 8 && z.size() >= 8) {
        OC_InducedVelocities ocr = oc_inflow_compute_wing_induced_vel_on_blade(inl(*this), x.data(), y.data(), z.data());
        std::memcpy(result.v_x, ocr.v_x, sizeof(result.v_x));
        std::memcpy(result.v_y, ocr.v_y, sizeof(result.v_y));
        std::memcpy(result.v_z, ocr.v_z, sizeof(result.v_z));
    }
    return result;
}

// ========================================================================
//  HuangPetersInflow
// ========================================================================

HuangPetersInflow::HuangPetersInflow(long mMo, long mMe, const RotorGeometry& rotor,
    const RotorInputState& ri, double dt) {
    OC_Inflow* raw = oc_huang_peters_create(mMo, mMe, rgp(rotor), ris(ri), dt);
    ptr_ = raw; owned_ = true;
}
HuangPetersInflow::~HuangPetersInflow() = default;
HuangPetersInflow::HuangPetersInflow(HuangPetersInflow&&) noexcept = default;
HuangPetersInflow& HuangPetersInflow::operator=(HuangPetersInflow&&) noexcept = default;

// ========================================================================
//  NullInflow
// ========================================================================

NullInflow::NullInflow(const RotorGeometry& rotor, const RotorInputState& ri) {
    OC_Inflow* raw = oc_null_inflow_create(rgp(rotor), ris(ri));
    ptr_ = raw; owned_ = true;
}
NullInflow::~NullInflow() = default;
NullInflow::NullInflow(NullInflow&&) noexcept = default;
NullInflow& NullInflow::operator=(NullInflow&&) noexcept = default;

// ========================================================================
//  WingInflow
// ========================================================================

WingInflow::WingInflow(const WingGeometry& wing, const WingInputState& wi, const WingLiftSurf& wls_) {
    OC_Inflow* raw = oc_wing_inflow_create(wg(wing), wis(wi), wls(wls_));
    ptr_ = raw; owned_ = true;
}
WingInflow::~WingInflow() = default;
WingInflow::WingInflow(WingInflow&&) noexcept = default;
WingInflow& WingInflow::operator=(WingInflow&&) noexcept = default;

// ========================================================================
//  RotorInputState
// ========================================================================

RotorInputState::RotorInputState(void* p) : ptr_(p) {}

void RotorInputState::set_angular_velocity(double o) { if(ptr_) oc_rotor_input_set_angular_velocity(ris(*this), o); }
double RotorInputState::angular_velocity() const { return ptr_ ? oc_rotor_input_get_angular_velocity(ris(*this)) : 0.0; }
void RotorInputState::set_angular_accel(double a) { if(ptr_) oc_rotor_input_set_angular_accel(ris(*this), a); }
double RotorInputState::angular_accel() const { return ptr_ ? oc_rotor_input_get_angular_accel(ris(*this)) : 0.0; }
void RotorInputState::set_azimuth(double a) { if(ptr_) oc_rotor_input_set_azimuth(ris(*this), a); }
double RotorInputState::azimuth() const { return ptr_ ? oc_rotor_input_get_azimuth(ris(*this)) : 0.0; }

void RotorInputState::set_r_0(const std::vector<double>& d) { if(ptr_) oc_rotor_input_set_r_0(ris(*this), const_cast<double*>(d.data()), d.size()); }
std::vector<double> RotorInputState::get_r_0(size_t len) const { std::vector<double> r(len); if(ptr_) oc_rotor_input_get_r_0(ris(*this), r.data(), len); return r; }

void RotorInputState::set_blade_flapping(const std::vector<double>& d) { if(ptr_) oc_rotor_input_set_blade_flapping(ris(*this), const_cast<double*>(d.data()), d.size()); }
std::vector<double> RotorInputState::get_blade_flapping(size_t len) const { std::vector<double> r(len); if(ptr_) oc_rotor_input_get_blade_flapping(ris(*this), r.data(), len); return r; }

void RotorInputState::set_blade_flapping_rate(const std::vector<double>& d) { if(ptr_) oc_rotor_input_set_blade_flapping_rate(ris(*this), const_cast<double*>(d.data()), d.size()); }
std::vector<double> RotorInputState::get_blade_flapping_rate(size_t len) const { std::vector<double> r(len); if(ptr_) oc_rotor_input_get_blade_flapping_rate(ris(*this), r.data(), len); return r; }

// ========================================================================
//  AircraftInputState
// ========================================================================

AircraftInputState::AircraftInputState(void* p) : ptr_(p) {}
AircraftInputState::AircraftInputState(size_t nr, const std::vector<size_t>& nb, size_t nw) {
    std::vector<size_t> m(nb);
    ptr_ = oc_aircraft_input_state_create(nr, m.data(), nw);
}
AircraftInputState::~AircraftInputState() { if(ptr_) oc_aircraft_input_state_destroy(ais(*this)); }
AircraftInputState::AircraftInputState(AircraftInputState&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
AircraftInputState& AircraftInputState::operator=(AircraftInputState&& o) noexcept {
    if (this != &o) { if(ptr_) oc_aircraft_input_state_destroy(ais(*this)); ptr_ = o.ptr_; o.ptr_ = nullptr; }
    return *this;
}

RotorInputState AircraftInputState::get_rotor_input(size_t i) {
    OC_RotorInputState* raw = ptr_ ? oc_aircraft_input_get_rotor_input(ais(*this), i) : nullptr;
    return RotorInputState(raw);
}
void AircraftInputState::set_blade_pitch(size_t r, size_t b, double p) {
    if(ptr_) oc_aircraft_input_set_blade_pitch(ais(*this), r, b, p);
}
double AircraftInputState::get_blade_pitch(size_t r, size_t b) const {
    return ptr_ ? oc_aircraft_input_get_blade_pitch(ais(*this), r, b) : 0.0;
}

// ========================================================================
//  RotorState
// ========================================================================

RotorState::RotorState(void* p) : ptr_(p) {}
double RotorState::get_C_T() const { double o=0; if(ptr_) oc_rotor_state_get_C_T(rst(*this),&o); return o; }
void RotorState::set_C_T(double c) { if(ptr_) oc_rotor_state_set_C_T(rst(*this),c); }
double RotorState::get_C_Q() const { double o=0; if(ptr_) oc_rotor_state_get_C_Q(rst(*this),&o); return o; }
void RotorState::set_C_Q(double c) { if(ptr_) oc_rotor_state_set_C_Q(rst(*this),c); }

// ========================================================================
//  BladeState
// ========================================================================

BladeState::BladeState(void* p) : ptr_(p) {}
double BladeState::azimuth() const { return ptr_ ? oc_blade_state_get_azimuth(bst(*this)) : 0; }
double BladeState::C_T() const { return ptr_ ? oc_blade_state_get_C_T(bst(*this)) : 0; }
double BladeState::C_Q() const { return ptr_ ? oc_blade_state_get_C_Q(bst(*this)) : 0; }
double BladeState::C_L() const { return ptr_ ? oc_blade_state_get_C_L(bst(*this)) : 0; }
double BladeState::C_D() const { return ptr_ ? oc_blade_state_get_C_D(bst(*this)) : 0; }

#define BS_VEC(name) std::vector<double> BladeState::name(size_t len) const { \
    std::vector<double> d(len); if(ptr_) oc_blade_state_fill_##name(bst(*this),d.data(),len); return d; }
BS_VEC(dC_T); BS_VEC(dC_Db); BS_VEC(dC_Db_profile); BS_VEC(dC_Db_induced);
BS_VEC(dynamic_dC_Db_profile); BS_VEC(dynamic_dC_Db_induced); BS_VEC(dC_N);
BS_VEC(dC_c); BS_VEC(dC_D); BS_VEC(dC_T_dot); BS_VEC(dC_Q); BS_VEC(dC_L);
BS_VEC(dC_l); BS_VEC(dC_Mz); BS_VEC(dC_My); BS_VEC(u_p); BS_VEC(dynamic_u_p);
BS_VEC(u_t); BS_VEC(aoa); BS_VEC(aoa_eff); BS_VEC(gamma); BS_VEC(r_c);
BS_VEC(x); BS_VEC(y); BS_VEC(z);
#undef BS_VEC

#define BS_VECF(name) std::vector<float> BladeState::name(size_t len) const { \
    std::vector<float> d(len); if(ptr_) oc_blade_state_fill_##name(bst(*this),d.data(),len); return d; }
BS_VECF(dC_Df); BS_VECF(dC_Nf); BS_VECF(dC_cf); BS_VECF(dC_Tf); BS_VECF(dC_Qf);
#undef BS_VECF

#define BS_FILL(name) void BladeState::fill_##name(double* data, size_t len) const { \
    if(ptr_) oc_blade_state_fill_##name(bst(*this),data,len); }
BS_FILL(dC_T); BS_FILL(dC_Db); BS_FILL(dC_N); BS_FILL(dC_D); BS_FILL(dC_L);
BS_FILL(dC_Q); BS_FILL(u_p); BS_FILL(u_t); BS_FILL(aoa); BS_FILL(gamma);
BS_FILL(x); BS_FILL(y); BS_FILL(z); BS_FILL(r_c);
#undef BS_FILL

// ========================================================================
//  AircraftState
// ========================================================================

AircraftState::AircraftState(void* p) : ptr_(p) {}

AircraftState::AircraftState(size_t nr, const std::vector<size_t>& nb, size_t ne, size_t nw,
    const std::vector<size_t>& nwp, size_t sns, size_t cns,
    const Aircraft& aircraft, std::span<Inflow*> ri, std::span<Inflow*> wi, Direction dir) {
    std::vector<OC_Inflow*> rri(ri.size());
    for (size_t i = 0; i < ri.size(); ++i) rri[i] = ri[i] ? inl(*ri[i]) : nullptr;
    std::vector<OC_Inflow*> wri(wi.size());
    for (size_t i = 0; i < wi.size(); ++i) wri[i] = wi[i] ? inl(*wi[i]) : nullptr;
    std::vector<size_t> mnb(nb), mw(nwp);
    double dv = static_cast<double>(static_cast<int>(dir));
    ptr_ = oc_aircraft_state_create(nr, mnb.data(), ne, nw, mw.data(), sns, cns,
        ap(aircraft), rri.data(), wri.data(), &dv);
}

AircraftState::~AircraftState() { if(ptr_) oc_aircraft_state_destroy(ast(*this)); }
AircraftState::AircraftState(AircraftState&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
AircraftState& AircraftState::operator=(AircraftState&& o) noexcept {
    if(this!=&o){if(ptr_)oc_aircraft_state_destroy(ast(*this));ptr_=o.ptr_;o.ptr_=nullptr;} return *this;
}

void AircraftState::set_freestream(const Vec4& f) { if(ptr_){OC_Vec4 v=to_oc(f); oc_aircraft_state_set_freestream(ast(*this),&v);} }
Vec4 AircraftState::get_freestream() const { OC_Vec4 o{}; if(ptr_) oc_aircraft_state_get_freestream(ast(*this),&o); return from_oc(o); }

double AircraftState::rotor_C_T(size_t i) { double o=0; if(ptr_) oc_aircraft_state_get_rotor_C_T(ast(*this),i,&o); return o; }
double AircraftState::rotor_C_Q(size_t i) { double o=0; if(ptr_) oc_aircraft_state_get_rotor_C_Q(ast(*this),i,&o); return o; }

// ========================================================================
//  VortexFilament
// ========================================================================

VortexFilament::VortexFilament(void* p) : ptr_(p) {}

#define VF_VEC(name) std::vector<double> VortexFilament::name(size_t len) const { \
    std::vector<double> d(len); if(ptr_) oc_vortex_filament_fill_##name(vf(*this),d.data(),len); return d; }
VF_VEC(x); VF_VEC(y); VF_VEC(z); VF_VEC(gamma); VF_VEC(r_c); VF_VEC(v_z);
#undef VF_VEC

#define VF_FILL(name) void VortexFilament::fill_##name(double* data, size_t len) const { \
    if(ptr_) oc_vortex_filament_fill_##name(vf(*this),data,len); }
VF_FILL(x); VF_FILL(y); VF_FILL(z); VF_FILL(gamma); VF_FILL(r_c); VF_FILL(v_z);
#undef VF_FILL

// ========================================================================
//  RotorWake
// ========================================================================

RotorWake::RotorWake(void* p) : ptr_(p) {}

VortexFilament RotorWake::get_tip_vortex(size_t b) {
    OC_VortexFilament* r = ptr_ ? oc_rotor_wake_get_tip_vortex(rw(*this), b) : nullptr;
    return VortexFilament(r);
}

// ========================================================================
//  Wake
// ========================================================================

Wake::Wake(void* p, bool owned) : ptr_(p), owned_(owned) {}

Wake::Wake(size_t nr, size_t nb, size_t wh_, size_t re,
    const std::vector<size_t>& sh, const std::vector<size_t>& sr) {
    std::vector<size_t> mh(sh), mr(sr);
    ptr_ = oc_wake_create(nr, nb, wh_, re, mh.data(), mr.data());
    owned_ = true;
}

Wake::~Wake() { if(owned_&&ptr_){oc_wake_destroy(wk(*this)); ptr_=nullptr;} }
Wake::Wake(Wake&& o) noexcept : ptr_(o.ptr_), owned_(o.owned_) { o.ptr_=nullptr; o.owned_=false; }
Wake& Wake::operator=(Wake&& o) noexcept {
    if(this!=&o){if(owned_&&ptr_)oc_wake_destroy(wk(*this)); ptr_=o.ptr_; owned_=o.owned_; o.ptr_=nullptr; o.owned_=false;} return *this;
}

RotorWake Wake::get_rotor_wake(size_t i) {
    OC_RotorWake* r = ptr_ ? oc_wake_get_rotor_wake(wk(*this), i) : nullptr;
    return RotorWake(r);
}

// ========================================================================
//  WakeHistory
// ========================================================================

WakeHistory::WakeHistory(void* p) : ptr_(p) {}

WakeHistory::WakeHistory(size_t nr, size_t nb, size_t wh, size_t th, size_t re,
    const std::vector<size_t>& sh, const std::vector<size_t>& sr, double a1, bool hybrid) {
    std::vector<size_t> mh(sh), mr(sr);
    ptr_ = oc_wake_history_create(nr, nb, wh, th, re, mh.data(), mr.data(), a1, hybrid ? 1 : 0);
}

WakeHistory::~WakeHistory() { if(ptr_) oc_wake_history_destroy(wh(*this)); }
WakeHistory::WakeHistory(WakeHistory&& o) noexcept : ptr_(o.ptr_) { o.ptr_ = nullptr; }
WakeHistory& WakeHistory::operator=(WakeHistory&& o) noexcept {
    if(this!=&o){if(ptr_)oc_wake_history_destroy(wh(*this));ptr_=o.ptr_;o.ptr_=nullptr;} return *this;
}

void WakeHistory::push_back() { if(ptr_) oc_wake_history_push_back(wh(*this)); }

Wake WakeHistory::get_wake(size_t i) {
    OC_Wake* r = ptr_ ? oc_wake_history_get_wake(wh(*this), i) : nullptr;
    return Wake(r, false);
}

// ========================================================================
//  WingInputState / WingState / WingLiftSurf
// ========================================================================

WingInputState::WingInputState(void* p) : ptr_(p) {}
WingState::WingState(void* p) : ptr_(p) {}
WingLiftSurf::WingLiftSurf(void* p) : ptr_(p) {}

void WingLiftSurf::set_vortex_geometry(const WingGeometry& wing, size_t sc, size_t cn) {
    if(ptr_) oc_wing_set_vortex_geometry(wls(*this), wg(wing), sc, cn);
}

// ========================================================================
//  VTK Types
// ========================================================================

VtkRotor::VtkRotor(void* p) : ptr_(p) {}
VtkRotor::~VtkRotor() { if(ptr_) oc_vtk_rotor_destroy(vr(*this)); }
VtkRotor::VtkRotor(VtkRotor&& o) noexcept : ptr_(o.ptr_) { o.ptr_=nullptr; }
VtkRotor& VtkRotor::operator=(VtkRotor&& o) noexcept {
    if(this!=&o){if(ptr_)oc_vtk_rotor_destroy(vr(*this));ptr_=o.ptr_;o.ptr_=nullptr;} return *this;
}
VtkRotor VtkRotor::build(const RotorGeometry& r) { return VtkRotor(oc_build_vtu_rotor(rgp(r))); }

VtkWing::VtkWing(void* p) : ptr_(p) {}
VtkWing::~VtkWing() { if(ptr_) oc_vtk_wing_destroy(vw(*this)); }
VtkWing::VtkWing(VtkWing&& o) noexcept : ptr_(o.ptr_) { o.ptr_=nullptr; }
VtkWing& VtkWing::operator=(VtkWing&& o) noexcept {
    if(this!=&o){if(ptr_)oc_vtk_wing_destroy(vw(*this));ptr_=o.ptr_;o.ptr_=nullptr;} return *this;
}
VtkWing VtkWing::build(const WingGeometry& w) { return VtkWing(oc_build_vtu_wing(wg(w))); }

VtkWake::VtkWake(void* p) : ptr_(p) {}
VtkWake::~VtkWake() { if(ptr_) oc_vtk_wake_destroy(vkw(*this)); }
VtkWake::VtkWake(VtkWake&& o) noexcept : ptr_(o.ptr_) { o.ptr_=nullptr; }
VtkWake& VtkWake::operator=(VtkWake&& o) noexcept {
    if(this!=&o){if(ptr_)oc_vtk_wake_destroy(vkw(*this));ptr_=o.ptr_;o.ptr_=nullptr;} return *this;
}
VtkWake VtkWake::build(const Wake& w) { return VtkWake(oc_build_vtu_wake(wk(w))); }

VtkWingWake::VtkWingWake(void* p) : ptr_(p) {}
VtkWingWake::~VtkWingWake() { if(ptr_) oc_vtk_wing_wake_destroy(vww(*this)); }
VtkWingWake::VtkWingWake(VtkWingWake&& o) noexcept : ptr_(o.ptr_) { o.ptr_=nullptr; }
VtkWingWake& VtkWingWake::operator=(VtkWingWake&& o) noexcept {
    if(this!=&o){if(ptr_)oc_vtk_wing_wake_destroy(vww(*this));ptr_=o.ptr_;o.ptr_=nullptr;} return *this;
}
VtkWingWake VtkWingWake::build(const WingGeometry& w, const WingLiftSurf& l) {
    return VtkWingWake(oc_build_vtu_wing_wake(wg(w), wls(l)));
}

// ========================================================================
//  VTK Write Functions
// ========================================================================

void write_rotor_vtu(std::string_view filename, size_t step, size_t iteration,
    const VtkRotor& vtk, const AircraftState& ac_state, const RotorGeometry& geom) {
    oc_write_rotor_vtu(filename.data(), step, iteration, vr(vtk), ast(ac_state), rgp(geom));
}

void write_rotors_vtu(std::string_view filename, size_t step,
    std::span<const VtkRotor*> vtks, const AircraftState& ac_state, const Aircraft& aircraft) {
    size_t n = vtks.size(); std::vector<OC_VtkRotor*> raw(n);
    for(size_t i=0;i<n;++i) raw[i] = vtks[i] ? vr(*vtks[i]) : nullptr;
    oc_write_rotors_vtu(filename.data(), step, raw.data(), n, ast(ac_state), ap(aircraft));
}

void write_wing_vtu(std::string_view filename, size_t step, size_t iteration,
    const VtkWing& vtk, const AircraftState& ac_state, const WingGeometry& geom) {
    oc_write_wing_vtu(filename.data(), step, iteration, vw(vtk), ast(ac_state), wg(geom));
}

void write_wake_vtu(std::string_view filename, size_t step, const VtkWake& vtk, const Wake& wake) {
    oc_write_wake_vtu(filename.data(), step, vkw(vtk), wk(wake));
}

void write_wing_wake_vtu(std::string_view filename, size_t step, size_t iteration,
    const VtkWingWake& vtk, const WingGeometry& wing, const WingLiftSurf& lift_surf, const WingInputState& input) {
    oc_write_wing_wake_vtu(filename.data(), step, iteration, vww(vtk), wg(wing), wls(lift_surf), wis(input));
}

void write_wake_field_vtu(std::string_view filename, const AircraftState& ac_state, const Wake& wake,
    double xmin,double xmax, double ymin,double ymax, double zmin,double zmax, size_t nx,size_t ny,size_t nz) {
    oc_write_wake_field_vtu(filename.data(), ast(ac_state), wk(wake),
        xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
}

// ========================================================================
//  Test function
// ========================================================================

void test_func() { std::cout << "cpp test from OpenCOPTER" << std::endl; }

} // namespace opencopter