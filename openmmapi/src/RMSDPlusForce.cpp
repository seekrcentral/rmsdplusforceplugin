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

#include "RMSDPlusForce.h"
#include "internal/RMSDPlusForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Vec3.h"
#include <iostream>

using namespace RMSDPlusForcePlugin;
using namespace OpenMM;
using namespace std;

/*
RMSDPlusForce::RMSDPlusForce(const vector<Vec3>& referencePositions,
                             const vector<int>& alignParticles,
                             const vector<int>& rmsdParticles) {
	setReferencePositions(referencePositions);
	setAlignParticles(alignParticles);
	setRMSDParticles(rmsdParticles);
*/

RMSDPlusForce::RMSDPlusForce(const vector<double>& referencePosition,
		                     const vector<int>& alignParticles,
							 const vector<int>& rmsdParticles) {
	//setReferencePositions(referencePositions);
	setAlignParticles(alignParticles);
	setRMSDParticles(rmsdParticles);

	// DEBUG
	/*
	cout << "DEBUG: initialize. referencePositions: \n";
	for (int i = 0; i < referencePositions.size(); i++) {
		cout << "x: " << referencePositions[i][0] << ", ";
		cout << "y: " << referencePositions[i][1] << ", ";
		cout << "z: " << referencePositions[i][2] << "\n";
	}
	*/

	cout << "DEBUG: initialize. referencePosition: \n";
	cout << "x: " << referencePosition[0] << ", ";
	cout << "y: " << referencePosition[1] << ", ";
	cout << "z: " << referencePosition[2] << "\n";

	cout << "DEBUG: initialize. alignParticles: [";
	for (int i = 0; i < alignParticles.size(); i++) {
		cout << alignParticles[i] << ", ";
	}
	cout << "]\n";
	cout << "DEBUG: initialize. rmsdParticles: [";
	for (int i = 0; i < rmsdParticles.size(); i++) {
		cout << rmsdParticles[i] << ", ";
	}
	cout << "]\n";
}

void RMSDPlusForce::setAlignParticles(const vector<int>& particles) {
    alignParticles = particles;
}

void RMSDPlusForce::setRMSDParticles(const vector<int>& particles) {
    rmsdParticles = particles;
}

void RMSDPlusForce::setReferencePositions(const vector<Vec3>& positions) {
    referencePositions = positions;
}

void RMSDPlusForce::getRMSDPlusAlignParameters(int index, int& particle) const {
	particle = alignParticles[index];
}

void RMSDPlusForce::getRMSDPlusRMSDParameters(int index, int& particle) const{
	particle = rmsdParticles[index];
}

void RMSDPlusForce::getRMSDPlusReferencePosition(int index, Vec3& position) const{
	position = referencePositions[index];
}

ForceImpl* RMSDPlusForce::createImpl() const {
    return new RMSDPlusForceImpl(*this);
}

void RMSDPlusForce::updateParametersInContext(Context& context) {
    dynamic_cast<RMSDPlusForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}
