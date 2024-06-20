%module rmsdplusforceplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%import(module="simtk.openmm") "swig/typemaps.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "factory.i"
%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "RMSDPlusForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/Vec3.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}

/*
 * Add units to function outputs.
*/

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

%pythonprepend RMSDPlusForcePlugin::RMSDPlusForce(const std::vector<double>& referencePosition, const std::vector<int>& alignParticles, const std::vector<int>& rmsdParticles) %{
    print("CONVERTING")
    if unit.is_quantity(referencePosition):
        referencePosition = referencePosition.value_in_unit(unit.nanometer)
%}

%pythonprepend RMSDPlusForcePlugin::RMSDPlusForce::setReferencePositions(const std::vector< Vec3 > &positions) %{
    if unit.is_quantity(positions):
        positions = positions.value_in_unit(unit.nanometer)
%}

namespace RMSDPlusForcePlugin {

class RMSDPlusForce : public OpenMM::Force {
public:
    //RMSDPlusForce(const std::vector<OpenMM::Vec3>& referencePositions, const std::vector<int>& alignParticles, const std::vector<int>& rmsdParticles);
    
    RMSDPlusForce(const std::vector<double>& referencePosition, const std::vector<int>& alignParticles, const std::vector<int>& rmsdParticles);
    
    int getForceGroup() const;
    
    std::string getName() const;
    
    const std::vector<int> getAlignParticles() const;
    
    const std::vector<int> getRMSDParticles() const;
    
    const std::vector<OpenMM::Vec3> getReferencePositions() const;
        
    void setForceGroup(int group);
    
    void setName(std::string name);
    
    void setAlignParticles(const std::vector<int>& particles);
    
    void setRMSDParticles(const std::vector<int>& particles);
    
    void setReferencePositions(const std::vector<OpenMM::Vec3>& positions);
    
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
