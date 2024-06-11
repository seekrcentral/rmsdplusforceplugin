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
#include "Vec3.h"
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
    RMSDPlusForce(std::vector<Vec3> referencePositions, 
                  std::vector<int> alignParticles,
                  std::vector<int> rmsdParticles);
    /**
     * Get the force group this force belongs to
     */
    
    std::vector<int> getAlignParticles() const {
        return alignParticles;
    }
    
    std::vector<int> getRMSDParticles() const {
        return rmsdParticles;
    }
    
    std::vector<Vec3> getReferencePositions() const {
        return referencePositions
    }
    
    void setAlignParticles(std::vector<int> particles);
    
    void setRMSDParticles(std::vector<int> particles);
    
    void setReferencePositions(std::vector<Vec3> positions);
    
    void updateParametersInContext(OpenMM::Context& context);
    
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
    
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    std::vector<Vec3> referencePositions;
    std::vector<int> alignParticles;
    std::vector<int> rmsdParticles;
};

} // namespace RMSDPlusForcePlugin

#endif /*OPENMM_RMSDPLUSFORCE_H_*/
