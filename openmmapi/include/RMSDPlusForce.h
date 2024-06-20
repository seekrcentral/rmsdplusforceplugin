#ifndef OPENMM_RMSDPLUSFORCE_H_
#define OPENMM_RMSDPLUSFORCE_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMMRMSDPlusForce                            *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2024                                                *
 * Authors: Anson Noland and Lane Votapka                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Context.h"
#include "openmm/Force.h"
#include "openmm/Vec3.h"
#include <vector>
#include "internal/windowsExportRMSDPlusForce.h"

namespace RMSDPlusForcePlugin {

/**
 * This class implements a RMSD force that has different selections
 * between the alignment and the RMSD calculation.
 */

class OPENMM_EXPORT_RMSDPLUSFORCE RMSDPlusForce : public OpenMM::Force {
public:
    /**
     * Create an RMSDPlusForce.
     */
    //RMSDPlusForce(const std::vector<OpenMM::Vec3>& referencePositions,
    //              const std::vector<int>& alignParticles,
    //              const std::vector<int>& rmsdParticles);

    RMSDPlusForce(const std::vector<double>& referencePosition,
    		      const std::vector<int>& alignParticles,
                  const std::vector<int>& rmsdParticles);

    /**
     * Get the force group this force belongs to
     */
    
    const std::vector<int> getAlignParticles() const {
        return alignParticles;
    }
    
    const std::vector<int> getRMSDParticles() const {
        return rmsdParticles;
    }
    
    const std::vector<OpenMM::Vec3> getReferencePositions() const {
        return referencePositions;
    }
    
    int getNumAlignParticles() const {
    	return alignParticles.size();
    }

    int getNumRMSDParticles() const {
		return rmsdParticles.size();
	}

    int getNumReferencePositions() const {
		return referencePositions.size();
	}

    void getRMSDPlusAlignParameters(int index, int& particle) const;

    void getRMSDPlusRMSDParameters(int index, int& particle) const;

    void getRMSDPlusReferencePosition(int index, OpenMM::Vec3& position) const;

    void setAlignParticles(const std::vector<int>& particles);
    
    void setRMSDParticles(const std::vector<int>& particles);
    
    void setReferencePositions(const std::vector<OpenMM::Vec3>& positions);
    
    void updateParametersInContext(OpenMM::Context& context);
    
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
    
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    std::vector<OpenMM::Vec3> referencePositions;
    std::vector<int> alignParticles;
    std::vector<int> rmsdParticles;
};

} // namespace RMSDPlusForcePlugin

#endif /*OPENMM_RMSDPLUSFORCE_H_*/
