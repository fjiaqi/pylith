# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#
# @file pylith/materials/IsotropicLinearElasticity.py
#
# @brief Python material for isotropic, linearly elastic, plane
# strain material.
#
# Factory: elasticity_rheology

from .RheologyElasticity import RheologyElasticity
from .materials import IsotropicLinearElasticity as ModuleLinearElasticity


class IsotropicLinearElasticity(RheologyElasticity, ModuleLinearElasticity):
    """Python material for isotropic, linearly elastic plane strain.

    FACTORY: elasticity_rheology
    """

    import pythia.pyre.inventory

    useReferenceState = pythia.pyre.inventory.bool("use_reference_state", default=False)
    useReferenceState.meta['tip'] = "Use reference stress/strain state."

    useSmallStrain = pythia.pyre.inventory.bool("use_small_strain", default=False)
    useSmallStrain.meta['tip'] = "Use small strain formulation."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="isotropiclinearelasticity"):
        """Constructor.
        """
        RheologyElasticity.__init__(self, name)
        return

    def _defaults(self):
        from .AuxSubfieldsIsotropicLinearElasticity import AuxSubfieldsIsotropicLinearElasticity
        self.auxiliarySubfields = AuxSubfieldsIsotropicLinearElasticity("auxiliary_subfields")

    def preinitialize(self, problem):
        RheologyElasticity.preinitialize(self, problem)

        ModuleLinearElasticity.useReferenceState(self, self.useReferenceState)
        ModuleLinearElasticity.useSmallStrain(self, self.useSmallStrain)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleLinearElasticity.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def elasticity_rheology():
    """Factory associated with IsotropicLinearElasticity.
    """
    return IsotropicLinearElasticity()


# End of file
