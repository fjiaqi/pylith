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

#include <portinfo>

#include "TestIsotropicLinearElasticitySmallStrain.hh" // Implementation of class methods

#include "pylith/materials/Elasticity.hh" // USES Elasticity
#include "pylith/materials/IsotropicLinearElasticitySmallStrain.hh" // USES IsotropicLinearElasticitySmallStrain

// ---------------------------------------------------------------------------------------------------------------------
// Setup testing data.
void
pylith::mmstests::TestIsotropicLinearElasticitySmallStrain::setUp(void) {
    TestElasticity::setUp();

    _rheology = new pylith::materials::IsotropicLinearElasticitySmallStrain();CPPUNIT_ASSERT(_rheology);

    CPPUNIT_ASSERT(_material);
    _material->setBulkRheology(_rheology);
} // setUp


// ---------------------------------------------------------------------------------------------------------------------
// Deallocate testing data.
void
pylith::mmstests::TestIsotropicLinearElasticitySmallStrain::tearDown(void) {
    delete _rheology;_rheology = NULL;

    TestElasticity::tearDown();
} // tearDown


// End of file
