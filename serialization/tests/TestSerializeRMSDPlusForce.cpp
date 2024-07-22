/* -------------------------------------------------------------------------- *
 *                                OpenMMRMSDPlusForce                         *
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

#include "RMSDPlusForce.h"
#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace RMSDPlusForcePlugin;
using namespace OpenMM;
using namespace std;

extern "C" void registerRMSDPlusForceSerializationProxies();

void testSerialization() {
	// Create a Force.

	vector<Vec3> refPos;
	for (int i = 0; i < 10; i++)
		refPos.push_back(Vec3(i/5.0, i*1.2, i*i/3.5));
	vector<int> align_particles;
	for (int i = 0; i < 5; i++)
		align_particles.push_back(i*i);
	vector<int> rmsd_particles;
	for (int i = 5; i < 10; i++)
		rmsd_particles.push_back(i*i);
	RMSDPlusForce force(refPos, align_particles, rmsd_particles);
	force.setForceGroup(3);

	// Serialize and then deserialize it.

	stringstream buffer;
	XmlSerializer::serialize<RMSDPlusForce>(&force, "Force", buffer);
	RMSDPlusForce* copy = XmlSerializer::deserialize<RMSDPlusForce>(buffer);

	// Compare the two forces to see if they are identical.

	RMSDPlusForce& force2 = *copy;
	ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
	ASSERT_EQUAL(force.getReferencePositions().size(), force2.getReferencePositions().size());
	for (int i = 0; i < force.getReferencePositions().size(); i++)
		ASSERT_EQUAL_VEC(force.getReferencePositions()[i], force2.getReferencePositions()[i], 0.0);
	ASSERT_EQUAL(force.getAlignParticles().size(), force2.getAlignParticles().size());
	for (int i = 0; i < force.getAlignParticles().size(); i++)
		ASSERT_EQUAL(force.getAlignParticles()[i], force2.getAlignParticles()[i]);
	ASSERT_EQUAL(force.getRMSDParticles().size(), force2.getRMSDParticles().size());
	for (int i = 0; i < force.getRMSDParticles().size(); i++)
		ASSERT_EQUAL(force.getRMSDParticles()[i], force2.getRMSDParticles()[i]);
}

int main() {
    try {
        registerRMSDPlusForceSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
