%module rmsdplusforceplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
%include "Vec3.h"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
  %template(vectorv) vector<Vec3>;
};

%{
#include "RMSDPlusForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
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

/* TODO: remove
%pythonappend RMSDPlusForcePlugin::RMSDPlusForce::getBondParameters(int index, int& particle1, int& particle2,
                                                             double& length, double& k) const %{
    val[2] = unit.Quantity(val[2], unit.nanometer)
    val[3] = unit.Quantity(val[3], unit.kilojoule_per_mole/unit.nanometer**4)
%}
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


namespace RMSDPlusForcePlugin {

class RMSDPlusForce : public OpenMM::Force {
public:
    RMSDPlusForce(std::vector<Vec3> referencePositions, std::vector<int> alignParticles, std::vector<int> rmsdParticles);

    int getForceGroup() const;
    
    std::string getName() const;
    
    std::vector<int> getAlignParticles() const;
    
    std::vector<int> getRMSDParticles() const;
    
    std::vector<Vec3> getReferencePositions() const;
    
    void setForceGroup(int group);
    
    void setName(std::string name);
    
    void setAlignParticles(std::vector<int> particles);
    
    void setRMSDParticles(std::vector<int> particles);
    
    void setReferencePositions(std::vector<Vec3> positions);
    
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
