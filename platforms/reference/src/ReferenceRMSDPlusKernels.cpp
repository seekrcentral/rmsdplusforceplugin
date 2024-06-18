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

#include "ReferenceRMSDPlusKernels.h"
#include "RMSDPlusForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"

using namespace RMSDPlusForcePlugin;
using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

void ReferenceCalcRMSDPlusForceKernel::initialize(const System& system, const RMSDPlusForce& force) {
    // Initialize bond parameters.
	int numAlignParticles = force.getNumAlignParticles();
	int numRMSDParticles = force.getNumRMSDParticles();
	int numReferencePositions = force.getNumReferencePositions();

    alignParticles.resize(numAlignParticles);
    rmsdParticles.resize(numRMSDParticles);
    referencePositions.resize(numReferencePositions);
    for (int i = 0; i < numAlignParticles; i++)
    	force.getRMSDPlusAlignParameters(i, alignParticles[i]);

	for (int i = 0; i < numRMSDParticles; i++)
		force.getRMSDPlusRMSDParameters(i, rmsdParticles[i]);

	for (int i = 0; i < numRMSDParticles; i++)
		force.getRMSDPlusReferencePosition(i, referencePositions[i]);

}

double ReferenceCalcRMSDPlusForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
    double energy = 0;

    // Compute the interactions.

    // Fill this out later

    return energy;
}

void ReferenceCalcRMSDPlusForceKernel::copyParametersToContext(ContextImpl& context, const RMSDPlusForce& force) {
    if (force.getNumAlignParticles() != alignParticles.size())
        throw OpenMMException("updateParametersInContext: The number of align particles has changed");
    if (force.getNumRMSDParticles() != rmsdParticles.size())
        throw OpenMMException("updateParametersInContext: The number of RMSD particles has changed");
    if (force.getNumReferencePositions() != referencePositions.size())
        throw OpenMMException("updateParametersInContext: The number of reference positions has changed");
    for (int i = 0; i < force.getNumAlignParticles(); i++) {
        int p1;
        force.getRMSDPlusAlignParameters(i, p1);
        if (p1 != alignParticles[i])
            throw OpenMMException("updateParametersInContext: An align particle index has changed");
    }
    for (int i = 0; i < force.getNumRMSDParticles(); i++) {
		int p1;
		force.getRMSDPlusRMSDParameters(i, p1);
		if (p1 != rmsdParticles[i])
			throw OpenMMException("updateParametersInContext: An RMSD particle index has changed");
	}

    for (int i = 0; i < force.getNumReferencePositions(); i++) {
		Vec3 pos;
		force.getRMSDPlusReferencePosition(i, pos);
		if (pos != referencePositions[i])
			throw OpenMMException("updateParametersInContext: A reference position has changed");
	}


}
