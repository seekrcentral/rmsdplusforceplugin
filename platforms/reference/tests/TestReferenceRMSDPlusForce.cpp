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

/**
 * This tests the Reference implementation of RMSDPlusForce.
 */

#include "RMSDPlusForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "sfmt/SFMT.h"

using namespace RMSDPlusForcePlugin;
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerRMSDPlusReferenceKernelFactories();

double estimateRMSDPlusCV(vector<OpenMM::Vec3>& positions, vector<OpenMM::Vec3>& referencePos, vector<int>& alignParticles, vector<int>& rmsdParticles) {
    // Estimate the RMSDCV.  For simplicity we omit the orientation alignment, but they should
    // already be almost perfectly aligned.

    Vec3 center1, center2;
    for (int i : alignParticles) {
        center1 += referencePos[i];
        center2 += positions[i];
    }
    center1 /= alignParticles.size();
    center2 /= alignParticles.size();
    double estimate = 0.0;
    for (int i : rmsdParticles) {
        Vec3 delta = (referencePos[i]-center1) - (positions[i]-center2);
        estimate += delta.dot(delta);
    }
    return sqrt(estimate/rmsdParticles.size());
}

void testRMSDCV() {
    Platform& platform = Platform::getPlatformByName("Reference");
    const int numParticles = 20;
    System system;
    vector<Vec3> referencePos(numParticles);
    vector<Vec3> positions(numParticles);
    vector<int> alignParticles;
    vector<int> rmsdParticles;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; ++i) {
        system.addParticle(1.0);
        referencePos[i] = Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*10;
        positions[i] = referencePos[i] + Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt))*0.2;
        if (i%5 != 0)
            alignParticles.push_back(i);
        if (i%5 == 0)
            rmsdParticles.push_back(i);
    }
    cout << "before\n";
    RMSDPlusForce* force = new RMSDPlusForce(referencePos, alignParticles, rmsdParticles);
    force->setAlignParticles(alignParticles);
    cout << "after\n";
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    double estimate = estimateRMSDPlusCV(positions, referencePos, alignParticles, rmsdParticles);

    // Have the force compute the RMSDPlusCV.  It should be very slightly less than
    // what we calculated above (since that omitted the rotation).

    State state1 = context.getState(State::Energy);
    double RMSDPlusCV = state1.getPotentialEnergy();
    ASSERT(RMSDPlusCV <= estimate);
    cout << "estimate:" << estimate << "\n";
    cout << "RMSDPlusCV:" << RMSDPlusCV << "\n";
    ASSERT(RMSDPlusCV > 0.8*estimate);

    // Translate and rotate all the particles.  This should have no effect on the RMSDCV.

    vector<Vec3> transformedPos(numParticles);
    double cs = cos(1.1), sn = sin(1.1);
    for (int i = 0; i < numParticles; i++) {
        Vec3 p = positions[i];
        transformedPos[i] = Vec3( cs*p[0] + sn*p[1] + 0.1,
                                 -sn*p[0] + cs*p[1] - 11.3,
                                  p[2] + 1.5);
    }
    context.setPositions(transformedPos);
    state1 = context.getState(State::Energy | State::Forces);
    ASSERT_EQUAL_TOL(RMSDPlusCV, state1.getPotentialEnergy(), 1e-4);

    // Take a small step in the direction of the energy gradient and see whether the potential energy changes by the expected amount.

    const vector<Vec3>& forces = state1.getForces();
    double norm = 0.0;
    for (int i = 0; i < (int) forces.size(); ++i)
        norm += forces[i].dot(forces[i]);
    norm = std::sqrt(norm);
    const double stepSize = 0.1;
    double step = 0.5*stepSize/norm;
    vector<Vec3> positions2(numParticles), positions3(numParticles);
    for (int i = 0; i < (int) positions.size(); ++i) {
        Vec3 p = transformedPos[i];
        Vec3 f = forces[i];
        positions2[i] = Vec3(p[0]-f[0]*step, p[1]-f[1]*step, p[2]-f[2]*step);
        positions3[i] = Vec3(p[0]+f[0]*step, p[1]+f[1]*step, p[2]+f[2]*step);
    }
    context.setPositions(positions2);
    State state2 = context.getState(State::Energy);
    context.setPositions(positions3);
    State state3 = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(norm, (state2.getPotentialEnergy()-state3.getPotentialEnergy())/stepSize, 1e-3);

    // Check that updateParametersInContext() works correctly.

    context.setPositions(transformedPos);
    force->setReferencePositions(transformedPos);
    force->updateParametersInContext(context);
    ASSERT_EQUAL_TOL(0.0, context.getState(State::Energy).getPotentialEnergy(), 1e-2);
    context.setPositions(referencePos);
    ASSERT_EQUAL_TOL(RMSDPlusCV, context.getState(State::Energy).getPotentialEnergy(), 1e-4);

}

int main() {
    try {
        registerRMSDPlusReferenceKernelFactories();
        testRMSDCV();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
