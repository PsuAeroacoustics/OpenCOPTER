/* ------------------------------------------------------------------ */
/*  test_cpp_api.cpp                                                   */
/*                                                                    */
/*  C++-level tests for the RAII wrapper API (opencopter.hpp).        */
/*  Exercises:                                                         */
/*    - Aircraft::num_rotors() / get_rotor()                          */
/*    - RotorGeometry::frame()                                        */
/*    - Span overloads (object-span, const span)                       */
/*    - OC_CHECK exception safety (null → std::runtime_error)         */
/* ------------------------------------------------------------------ */

#include <gtest/gtest.h>
#include "opencopter.hpp"
#include <stdexcept>
#include <array>
#include <span>
#include <vector>

using namespace opencopter;

/* ================================================================== */
/*  Step A1: Aircraft::num_rotors()                                   */
/* ================================================================== */

TEST(CPP_Aircraft, NumRotors) {
    Aircraft ac(2, 0);
    ASSERT_TRUE(static_cast<bool>(ac));
    EXPECT_EQ(ac.num_rotors(), 2u);
}

TEST(CPP_Aircraft, NumRotorsZero) {
    Aircraft ac(0, 0);
    ASSERT_TRUE(static_cast<bool>(ac));
    EXPECT_EQ(ac.num_rotors(), 0u);
}

TEST(CPP_Aircraft, NumRotorsNull) {
    Aircraft ac;  // default-constructed → null pointer
    EXPECT_FALSE(static_cast<bool>(ac));
    EXPECT_EQ(ac.num_rotors(), 0u);  // getter returns 0, no throw
}

/* ================================================================== */
/*  Step A2: Aircraft::get_rotor()                                     */
/* ================================================================== */

TEST(CPP_Aircraft, GetRotorValid) {
    Aircraft ac(2, 0);
    ASSERT_TRUE(static_cast<bool>(ac));
    RotorGeometry r0 = ac.get_rotor(0);
    RotorGeometry r1 = ac.get_rotor(1);
    EXPECT_TRUE(static_cast<bool>(r0));
    EXPECT_TRUE(static_cast<bool>(r1));
}

TEST(CPP_Aircraft, GetRotorOOB) {
    Aircraft ac(1, 0);
    ASSERT_TRUE(static_cast<bool>(ac));
    RotorGeometry r = ac.get_rotor(99);  // out of bounds
    EXPECT_FALSE(static_cast<bool>(r));  // null wrapper, no throw
}

TEST(CPP_Aircraft, GetRotorNullAircraft) {
    Aircraft ac;  // default-constructed → null pointer
    RotorGeometry r = ac.get_rotor(0);
    EXPECT_FALSE(static_cast<bool>(r));  // null wrapper, no throw
}

/* ================================================================== */
/*  Step A3: RotorGeometry::frame()                                    */
/* ================================================================== */

TEST(CPP_RotorGeo, FrameSetGet) {
    RotorGeometry rg(1, Vec3{0,0,0}, 0.5, 0.05);
    ASSERT_TRUE(static_cast<bool>(rg));

    Frame f(Vec3{0,0,1}, 0.0, Vec3{0,0,0}, nullptr, "test_rotor_frame", FrameType::Rotor);
    ASSERT_TRUE(static_cast<bool>(f));

    rg.set_frame(f);
    Frame got = rg.frame();
    EXPECT_TRUE(static_cast<bool>(got));
}

TEST(CPP_RotorGeo, FrameNull) {
    RotorGeometry rg;  // default-constructed → null pointer
    Frame f = rg.frame();
    EXPECT_FALSE(static_cast<bool>(f));  // null Frame, no throw
}

/* ================================================================== */
/*  Step A4: Span overloads                                            */
/* ================================================================== */

TEST(CPP_Aircraft, SetRotorsSpan) {
    // Object span: std::span<const RotorGeometry>
    std::array<RotorGeometry, 2> rotors{};
    rotors[0] = RotorGeometry(1, Vec3{0,0,0}, 0.5, 0.05);
    rotors[1] = RotorGeometry(1, Vec3{1,0,0}, 0.5, 0.05);

    Aircraft ac(2, 0);
    ASSERT_TRUE(static_cast<bool>(ac));

    std::span<const RotorGeometry> rotor_span(rotors.data(), 2);
    ac.set_rotors(rotor_span);  // should not throw
}

/* NOTE: CPP_BladeGeo.SetTwistSpan and CPP_RotorGeo.SetBladesSpan are
   skipped because they require BladeAirfoil::create (C++ wrapper over
   oc_blade_airfoil_create), which hits a pre-existing D-side bug where
   the BladeAirfoil constructor throws a D exception across the FFI
   boundary. The underlying C functions (oc_blade_geometry_set_twist,
   oc_rotor_geometry_set_blades) are already validated at C level in
   test_api_bladegroup.cpp and test_api_rotorgeo_advanced.cpp. */

/* ================================================================== */
/*  Step A5: OC_CHECK exception safety                                 */
/* ================================================================== */

TEST(CPP_Errors, SetRotorsNull) {
    Aircraft ac;  // default-constructed → null pointer
    std::array<RotorGeometry, 1> rotors{};
    rotors[0] = RotorGeometry(1, Vec3{0,0,0}, 0.5, 0.05);
    std::span<const RotorGeometry> rotor_span(rotors.data(), 1);
    EXPECT_THROW(ac.set_rotors(rotor_span), std::runtime_error);
}

TEST(CPP_Errors, SetSolidityNull) {
    RotorGeometry rg;  // default-constructed → null pointer
    EXPECT_THROW(rg.set_solidity(0.1), std::runtime_error);
}

TEST(CPP_Errors, SetFrameNull) {
    Frame f(Vec3{0,0,1}, 0.0, Vec3{0,0,0}, nullptr, "f", FrameType::Rotor);
    RotorGeometry rg;  // default-constructed → null pointer
    EXPECT_THROW(rg.set_frame(f), std::runtime_error);
}

/* ================================================================== */
/*  Span overloads                                                    */
/* ================================================================== */

/* Helper: build a minimal valid BladeAirfoil (move-only types need push_back).
 * The extent must span a whole chunk: BladeAirfoil sizes itself as
 * (extent_end - extent_start + 1)/chunk_size, so a sub-chunk extent yields
 * zero chunks. chunk_size is 8, matching the C-API tests. */
static BladeAirfoil make_test_blade_airfoil() {
    AirfoilModel af = AirfoilModel::thin_airfoil(6.28);
    std::vector<AirfoilModel> models;
    models.push_back(std::move(af));
    std::vector<size_t> extents = {0, 7};
    return BladeAirfoil::create(models, extents);
}

TEST(CPP_BladeGeo, SetTwistSpan) {
    BladeAirfoil ba = make_test_blade_airfoil();
    ASSERT_TRUE(static_cast<bool>(ba));

    BladeGeometry bg(8, 0.0, 0.3, ba, 0.5);
    ASSERT_TRUE(static_cast<bool>(bg));

    // Use std::span<const double> overload
    double twist[8] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    std::span<const double> twist_span(twist, 8);
    bg.set_twist(twist_span);  // should not throw
}

TEST(CPP_RotorGeo, SetBladesSpan) {
    // Object span: std::span<const BladeGeometry>
    BladeAirfoil ba = make_test_blade_airfoil();
    ASSERT_TRUE(static_cast<bool>(ba));

    std::array<BladeGeometry, 1> blades{};
    blades[0] = BladeGeometry(8, 0.0, 0.3, ba, 0.5);
    ASSERT_TRUE(static_cast<bool>(blades[0]));

    RotorGeometry rg(1, Vec3{0,0,0}, 0.5, 0.05);
    ASSERT_TRUE(static_cast<bool>(rg));

    std::span<const BladeGeometry> blade_span(blades.data(), 1);
    rg.set_blades(blade_span);  // should not throw
}
