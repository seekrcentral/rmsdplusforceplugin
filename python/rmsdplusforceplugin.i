%module rmsdplusforceplugin
%include "factory.i"

%{
#include "RMSDPlusForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/Vec3.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%{
#include <numpy/arrayobject.h>
%}

%include "std_vector.i"
%include "header.i"
%include "typemaps.i"

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */
using namespace OpenMM;
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%pythoncode %{
import openmm as mm
import openmm.unit as unit
import sys
import numpy
%}

/*
 * Add units to function outputs.
*/

%pythonappend RMSDPlusForcePlugin::RMSDPlusForce::getReferencePositions() const %{
    val = unit.Quantity(val, unit.nanometer)
%}

%pythonappend RMSDPlusForcePlugin::RMSDPlusForce::RMSDPlusForce(const std::vector<Vec3>& referencePositions, const std::vector<int>& alignParticles=std::vector<int>(), const std::vector<int>& rmsdParticles=std::vector<int>())%{
    for arg in args:
        if 'numpy' in sys.modules and isinstance(arg, numpy.ndarray):
            arg = arg.tolist()
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}
namespace RMSDPlusForcePlugin {

class RMSDPlusForce : public OpenMM::Force {
public:
    RMSDPlusForce(const std::vector<Vec3>& referencePositions, const std::vector<int>& alignParticles=std::vector<int>(), const std::vector<int>& rmsdParticles=std::vector<int>());
        
    const std::vector<int>& getAlignParticles() const;
    
    const std::vector<int>& getRMSDParticles() const;
    
    const std::vector<Vec3>& getReferencePositions() const;
        
    void setAlignParticles(const std::vector<int>& particles);
    
    void setRMSDParticles(const std::vector<int>& particles);
    
    void setReferencePositions(const std::vector<Vec3>& positions);
    
    void updateParametersInContext(OpenMM::Context& context);
    
    bool usesPeriodicBoundaryConditions() const;

    /*
     * Add methods for casting a Force to an RMSDPlusForce.
    */
    %extend {
        static RMSDPlusForcePlugin::RMSDPlusForce& cast(OpenMM::Force& force) {
            return dynamic_cast<RMSDPlusForcePlugin::RMSDPlusForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<RMSDPlusForcePlugin::RMSDPlusForce*>(&force) != NULL);
        }
    }
};

}
