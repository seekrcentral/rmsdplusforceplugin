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
 * This tests the CUDA implementation of RMSDPlusForce.
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

extern "C" OPENMM_EXPORT void registerRMSDPlusCudaKernelFactories();

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

void testRMSDPlusCV() {
    cout << "mark0\n";
    Platform& platform = Platform::getPlatformByName("CUDA");
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
        //if (i%5 != 0)
        //    alignParticles.push_back(i);
        //if (i%5 == 0)
        //    rmsdParticles.push_back(i);
    }
    positions[0][0] = 1.1935; positions[0][1] =  1.1972; positions[0][2] = 1.2004;
    positions[1][0] = 1.3461; positions[1][1] =  2.1199; positions[1][2] = 2.0782;
    positions[2][0] = 1.7296; positions[2][1] =  0.2857; positions[2][2] = 2.4467;
    positions[3][0] = 0.5572; positions[3][1] =  1.2005; positions[3][2] = 1.1523;
    positions[4][0] = 1.5434; positions[4][1] =  1.2144; positions[4][2] = 1.6613;
    positions[5][0] = 1.9985; positions[5][1] =  1.5946; positions[5][2] = 0.2312;
    positions[6][0] = 1.9237; positions[6][1] =  1.2023; positions[6][2] = 2.1215;
    positions[7][0] = 1.1097; positions[7][1] =  1.1192; positions[7][2] = 0.4674;
    positions[8][0] = 1.4092; positions[8][1] =  2.0307; positions[8][2] = 0.4728;
    positions[9][0] = -0.2170; positions[9][1] = 2.0325; positions[9][2] = 2.2613;

    referencePos[0][0] = 1.1835; referencePos[0][1] =  1.1872; referencePos[0][2] = 1.2204;
    referencePos[1][0] = 1.3361; referencePos[1][1] =  2.1299; referencePos[1][2] = 2.0882;
    referencePos[2][0] = 1.7096; referencePos[2][1] =  0.2657; referencePos[2][2] = 2.4567;
    referencePos[3][0] = 0.5672; referencePos[3][1] =  1.2205; referencePos[3][2] = 1.1623;
    referencePos[4][0] = 1.5734; referencePos[4][1] =  1.2344; referencePos[4][2] = 1.6513;
    referencePos[5][0] = 1.9885; referencePos[5][1] =  1.5746; referencePos[5][2] = 0.2212;
    referencePos[6][0] = 1.9337; referencePos[6][1] =  1.2323; referencePos[6][2] = 2.1115;
    referencePos[7][0] = 1.1197; referencePos[7][1] =  1.1092; referencePos[7][2] = 0.4574;
    referencePos[8][0] = 1.4292; referencePos[8][1] =  2.0207; referencePos[8][2] = 0.4628;
    referencePos[9][0] = -0.2270; referencePos[9][1] = 2.0225; referencePos[9][2] = 2.2513;

    alignParticles.push_back(0);
    alignParticles.push_back(1);
    alignParticles.push_back(3);
    alignParticles.push_back(7);
    rmsdParticles.push_back(2);
    rmsdParticles.push_back(4);
    rmsdParticles.push_back(5);

    RMSDPlusForce* force = new RMSDPlusForce(referencePos, alignParticles, rmsdParticles);
    force->setAlignParticles(alignParticles);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);

    // TODO; remove
    double estimate = estimateRMSDPlusCV(positions, referencePos, alignParticles, rmsdParticles);
    cout << "mark50\n";
    // Have the force compute the RMSDPlusCV.  It should be very slightly less than
    // what we calculated above (since that omitted the rotation).
    double expected_rmsd = 0.050743;
    double TOLERANCE = 0.0001;
    State state1 = context.getState(State::Energy);
    double RMSDPlusCV = state1.getPotentialEnergy();
    ASSERT(RMSDPlusCV < expected_rmsd + TOLERANCE);
    ASSERT(RMSDPlusCV > expected_rmsd - TOLERANCE);

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
    cout << "mark100\n";
}

int main() {
    try {
        registerRMSDPlusCudaKernelFactories();
        testRMSDPlusCV();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
