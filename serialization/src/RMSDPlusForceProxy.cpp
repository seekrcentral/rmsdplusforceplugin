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

#include "RMSDPlusForceProxy.h"
#include "RMSDPlusForce.h"
#include "openmm/Vec3.h"
#include "openmm/serialization/SerializationNode.h"
#include <sstream>
#include <vector>

using namespace RMSDPlusForcePlugin;
using namespace OpenMM;
using namespace std;

RMSDPlusForceProxy::RMSDPlusForceProxy() : SerializationProxy("RMSDPlusForce") {
}

void RMSDPlusForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const RMSDPlusForce& force = *reinterpret_cast<const RMSDPlusForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    SerializationNode& alignParticlesNode = node.createChildNode("alignParticles");
    int particleIndex;
    for (int i = 0; i < force.getNumAlignParticles(); i++) {
    	force.getRMSDPlusAlignParameters(i, particleIndex);
    	alignParticlesNode.createChildNode("alignParticle").setIntProperty("particleIndex", particleIndex);
    }
    SerializationNode& rmsdParticlesNode = node.createChildNode("rmsdParticles");
	for (int i = 0; i < force.getNumRMSDParticles(); i++) {
		force.getRMSDPlusRMSDParameters(i, particleIndex);
		rmsdParticlesNode.createChildNode("rmsdParticle").setIntProperty("particleIndex", particleIndex);
	}

	SerializationNode& referencePositionsNode = node.createChildNode("referencePositions");
	Vec3 referencePosition;
	for (int i = 0; i < force.getNumReferencePositions(); i++) {
		SerializationNode& referencePositionNode = referencePositionsNode.createChildNode("referencePosition");
		force.getRMSDPlusReferencePosition(i, referencePosition);
		referencePositionNode.setDoubleProperty("x", referencePosition[0]);
		referencePositionNode.setDoubleProperty("y", referencePosition[1]);
		referencePositionNode.setDoubleProperty("z", referencePosition[2]);
	}


}

void* RMSDPlusForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");

    RMSDPlusForce* force = NULL;
    try {
		vector<Vec3> referencePositions;
		vector<int> alignParticles;
		vector<int> rmsdParticles;

		const SerializationNode& alignParticlesNode = node.getChildNode("alignParticles");
		for (int i = 0; i < (int) alignParticlesNode.getChildren().size(); i++) {
			const SerializationNode& alignParticleNode = alignParticlesNode.getChildren()[i];
			alignParticles.push_back(alignParticleNode.getIntProperty("particleIndex"));
		}
		const SerializationNode& rmsdParticlesNode = node.getChildNode("rmsdParticles");
		for (int i = 0; i < (int) rmsdParticlesNode.getChildren().size(); i++) {
			const SerializationNode& rmsdParticleNode = rmsdParticlesNode.getChildren()[i];
			rmsdParticles.push_back(rmsdParticleNode.getIntProperty("particleIndex"));
		}

		const SerializationNode& referencePositionsNode = node.getChildNode("referencePositions");
		for (int i = 0; i < (int) referencePositionsNode.getChildren().size(); i++) {
			const SerializationNode& referencePositionNode = referencePositionsNode.getChildren()[i];
			Vec3 position;
			position[0] = referencePositionNode.getDoubleProperty("x");
			position[1] = referencePositionNode.getDoubleProperty("y");
			position[2] = referencePositionNode.getDoubleProperty("z");
			referencePositions.push_back(position);
		}
		force = new RMSDPlusForce(referencePositions, alignParticles, rmsdParticles);
		force->setForceGroup(node.getIntProperty("forceGroup", 0));
		return force;
    }
    catch (...) {
    	if (force != NULL)
			delete force;
		throw;
    }


}
