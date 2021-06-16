// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/**
 * @file tests/mmstests/elasticity/TestIsotropicLinearElasticitySmallStrain.hh
 *
 * @brief C++ TestIsotropicLinearElasticitySmallStrain object
 *
 * C++ unit testing for IsotropicLinearElasticitySmallStrain.
 */

#if !defined(pylith_mmstests_testisotropiclinearelasticitySmallStrain_hh)
#define pylith_mmstests_testisotropiclinearelasticitySmallStrain_hh

#include "TestElasticity.hh" // ISA TestElasticity

#include "pylith/materials/IsotropicLinearElasticitySmallStrain.hh" // HOLDSA IsotropicLinearElasticitySmallStrain

namespace pylith {
    namespace mmstests {
        class TestIsotropicLinearElasticitySmallStrain;
    } // tests/mmstests
} // pylith

/// C++ unit testing for IsotropicLinearElasticitySmallStrain
class pylith::mmstests::TestIsotropicLinearElasticitySmallStrain : public pylith::mmstests::TestElasticity {
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Setup testing data.
    void setUp(void);

    /// Deallocate testing data.
    void tearDown(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::materials::IsotropicLinearElasticitySmallStrain* _rheology; ///< Rheology for testing.

}; // class TestIsotropicLinearElasticitySmallStrain

#endif // pylith_mmstests_testisotropiclinearelasticitySmallStrain_hh

// End of file
