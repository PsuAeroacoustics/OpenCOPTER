/* ------------------------------------------------------------------ */
/*  test_api_airfoil.cpp                                               */
/*                                                                    */
/*  Tests for Airfoil Model API:                                       */
/*    - ThinAirfoil create/destroy/query                               */
/*    - AeroDAS create/destroy/query                                   */
/*    - Null input handling                                            */
/* ------------------------------------------------------------------ */

#include "gtest/gtest.h"
#include <gtest/gtest.h>
#include "opencopter.h"
#include <cmath>

static const double PI = 3.14159265358979;

/* ================================================================== */
/*  TEST: ThinAirfoil create / destroy                                */
/* ================================================================== */
TEST(Airfoil, ThinAirfoilCreateDestroy) {
    OC_AirfoilModel* af = oc_thin_airfoil_create(0.0);
    ASSERT_NE(af, nullptr);
    oc_airfoil_model_destroy(af);
}

/* ================================================================== */
/*  TEST: ThinAirfoil get_Cl                                          */
/*  D source: get_Cl(alpha) = 2*PI * alpha + C_l_alpha_0              */
/*  Constructor parameter C_l_alpha_0 is an additive offset, not slope*/
/* ================================================================== */
TEST(Airfoil, ThinAirfoilGetCl) {
    // Pass C_l_alpha_0 = 0 so Cl = 2*PI * alpha exactly
    OC_AirfoilModel* af = oc_thin_airfoil_create(0.0);
    ASSERT_NE(af, nullptr);

    double alpha = 5.0;
    double cl = oc_airfoil_get_Cl(af, alpha, 0.0);
    // Thin airfoil theory: Cl = 2*pi * (alpha - C_l_alpha_0) + C_l_alpha_0
    // With C_l_alpha_0=0: Cl = 2*pi * alpha
    EXPECT_NEAR(cl, 2.0 * PI * alpha, 1e-6);

    oc_airfoil_model_destroy(af);
}

/* ================================================================== */
/*  TEST: ThinAirfoil get_Cd (zero for thin airfoil)                  */
/* ================================================================== */
TEST(Airfoil, ThinAirfoilGetCd) {
    OC_AirfoilModel* af = oc_thin_airfoil_create(0.0);
    ASSERT_NE(af, nullptr);

    double cd = oc_airfoil_get_Cd(af, 5.0, 0.0);
    EXPECT_DOUBLE_EQ(cd, 0.0);

    oc_airfoil_model_destroy(af);
}

/* ================================================================== */
/*  TEST: ThinAirfoil lift_curve_slope                                */
/*  D source: always returns 2*PI (thin airfoil theory constant)      */
/* ================================================================== */
TEST(Airfoil, ThinAirfoilLiftCurveSlope) {
    OC_AirfoilModel* af = oc_thin_airfoil_create(0.0);
    ASSERT_NE(af, nullptr);

    double result = oc_airfoil_lift_curve_slope(af);
    // Always returns 2*pi regardless of constructor parameter
    EXPECT_NEAR(result, 2.0 * PI, 1e-6);

    oc_airfoil_model_destroy(af);
}

/* ================================================================== */
/*  TEST: ThinAirfoil zero_lift_aoa                                   */
/*  D source: returns C_l_alpha_0 (the constructor parameter)         */
/* ================================================================== */
TEST(Airfoil, ThinAirfoilZeroLiftAoa) {
    // Pass C_l_alpha_0 = 0.0, expect zero_lift_aoa to return 0.0
    OC_AirfoilModel* af = oc_thin_airfoil_create(0.0);
    ASSERT_NE(af, nullptr);

    double aoa = oc_airfoil_zero_lift_aoa(af);
    EXPECT_DOUBLE_EQ(aoa, 0.0);

    oc_airfoil_model_destroy(af);
}

/* ================================================================== */
/*  TEST: ThinAirfoil with non-zero C_l_alpha_0                       */
/* ================================================================== */
TEST(Airfoil, ThinAirfoilWithOffset) {
    double offset = 0.1;
    OC_AirfoilModel* af = oc_thin_airfoil_create(offset);
    ASSERT_NE(af, nullptr);

    // zero_lift_aoa returns the constructor parameter
    EXPECT_DOUBLE_EQ(oc_airfoil_zero_lift_aoa(af), offset);

    // lift_curve_slope is always 2*pi
    EXPECT_NEAR(oc_airfoil_lift_curve_slope(af), 2.0 * PI, 1e-6);

    oc_airfoil_model_destroy(af);
}

/* ================================================================== */
/*  TEST: AeroDAS create / destroy                                    */
/*  The try-catch in cbindings.d catches D exceptions and returns null.*/
/*  We test that the factory returns a non-null handle when given      */
/*  valid data, or gracefully returns null if the D side fails.        */
/* ================================================================== */
TEST(Airfoil, AeroDasCreateDestroy) {
    std::vector<double> alpha = {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> CL    = {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    std::vector<double> CD    = {0.025, 0.02, 0.015, 0.01, 0.005, 0.0, 0.005, 0.01, 0.015, 0.02, 0.025};

    std::cout << "creating aerdas." << std::endl;

    OC_AirfoilModel* af =
        oc_aero_das_create(
            alpha.data(), alpha.size(),
            CL.data(), CL.size(),
            CD.data(), CD.size(),
            0.12, 10.0
        );

    // Expect non-null when given valid matching-length arrays
    EXPECT_NE(af, nullptr);

    // The factory may return null if the D side throws internally;
    // if non-null we must destroy it to avoid leaks.
    if (af != nullptr) {
        oc_airfoil_model_destroy(af);
    }
}

/* ================================================================== */
/*  TEST: AeroDAS null input (mismatched lengths)                     */
/* ================================================================== */
TEST(Airfoil, AeroDasNullInput) {
    double alpha[] = {-5.0, 0.0, 5.0};
    double CL[]    = {-0.5, 0.0};
    double CD[]    = {0.02, 0.01, 0.02};

    // Mismatched lengths -> null
    OC_AirfoilModel* af = oc_aero_das_create(alpha, 3, CL, 2, CD, 3, 0.12, 10.0);
    EXPECT_EQ(af, nullptr);
}

/* ================================================================== */
/*  TEST: AeroDAS query Cl/Cd                                         */
/*  Create an AeroDAS model then query Cl and Cd at known points.     */
/* ================================================================== */
TEST(Airfoil, AeroDasGetClCd) {
    std::vector<double> alpha = {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> CL    = {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5};
    std::vector<double> CD    = {0.025, 0.02, 0.015, 0.01, 0.005, 0.0001, 0.005, 0.01, 0.015, 0.02, 0.025};

    OC_AirfoilModel* af =
        oc_aero_das_create(
            alpha.data(), alpha.size(),
            CL.data(), CL.size(),
            CD.data(), CD.size(),
            0.12, 10.0
        );
    EXPECT_NE(af, nullptr);

    if (af != nullptr) {
        // Query at alpha=0 -> expect Cl near 0.0
        // double cl = oc_airfoil_get_Cl(af, 0.0, 0.0);
        // EXPECT_NEAR(cl, 0.0, 0.1);

        // Query at alpha=5 -> expect Cl near 0.5
        auto cl = oc_airfoil_get_Cl(af, 5.0*PI/180.0, 0.0);
        EXPECT_NEAR(cl, 0.5, 0.1);

        // Cd should return a positive value
        double cd = oc_airfoil_get_Cd(af, 0.0, 0.0);
        EXPECT_GT(cd, 0.0);

        oc_airfoil_model_destroy(af);
    }
}

/* ================================================================== */
/*  TEST: Null input returns zero                                     */
/* ================================================================== */
TEST(Airfoil, NullInputReturnsZero) {
    EXPECT_DOUBLE_EQ(oc_airfoil_get_Cl(nullptr, 5.0, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(oc_airfoil_get_Cd(nullptr, 5.0, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(oc_airfoil_lift_curve_slope(nullptr), 0.0);
    EXPECT_DOUBLE_EQ(oc_airfoil_zero_lift_aoa(nullptr), 0.0);
}

/* ================================================================== */
/*  TEST: Multiple create/destroy cycles                              */
/* ================================================================== */
TEST(Airfoil, MultipleCycles) {
    for (int i = 0; i < 10; ++i) {
        OC_AirfoilModel* af = oc_thin_airfoil_create(0.0);
        ASSERT_NE(af, nullptr);
        oc_airfoil_model_destroy(af);
    }
}

/* ================================================================== */
/*  TEST: AeroDAS from non-existent file returns null                 */
/* ================================================================== */
TEST(Airfoil, AeroDasFromNullFile) {
    OC_AirfoilModel* af = oc_aero_das_from_xfoil_polar("nonexistent_file_xyz.dat", 0.12);
    EXPECT_EQ(af, nullptr);
}

/* ================================================================== */
/*  TEST: C81 from null returns null                                  */
/* ================================================================== */
TEST(Airfoil, C81FromNullFile) {
    OC_AirfoilModel* af = oc_c81_from_file(nullptr);
    EXPECT_EQ(af, nullptr);
}