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

#include "openmm/Vec3.h"
#include "jama_eig.h"

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
	cout << "mark1\n";
	int numAlignParticles = force.getNumAlignParticles();
	int numRMSDParticles = force.getNumRMSDParticles();
	int numReferencePositions = force.getNumReferencePositions();

	if (numAlignParticles == 0)
		throw OpenMMException("alignParticles may not be empty.");

	if (numRMSDParticles == 0)
		throw OpenMMException("rmsdParticles may not be empty.");

	if (numReferencePositions != system.getNumParticles())
		throw OpenMMException("referencePositions must have the same number of particles as the system.");

    alignParticles.resize(numAlignParticles);
    rmsdParticles.resize(numRMSDParticles);
    referencePositions.resize(numReferencePositions);
    for (int i = 0; i < numAlignParticles; i++)
    	force.getRMSDPlusAlignParameters(i, alignParticles[i]);

	for (int i = 0; i < numRMSDParticles; i++)
		force.getRMSDPlusRMSDParameters(i, rmsdParticles[i]);

	for (int i = 0; i < numRMSDParticles; i++)
		force.getRMSDPlusReferencePosition(i, referencePositions[i]);

	Vec3 alignCenter;
	for (int i : alignParticles)
		alignCenter += referencePositions[i];
	alignCenter /= alignParticles.size();
	for (Vec3& p : referencePositions)
		p -= alignCenter;
}

double ReferenceCalcRMSDPlusForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
	cout << "mark2\n";
	vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& forces = extractForces(context);

    int numAlignParticles = alignParticles.size();
	Vec3 alignCenter;
	alignCenter[0] = 0;
	alignCenter[1] = 0;
	alignCenter[2] = 0;
	for (int i : alignParticles)
		alignCenter += pos[i];
	alignCenter /= numAlignParticles;
	vector<Vec3> alignPositions(numAlignParticles);
	for (int i = 0; i < numAlignParticles; i++)
		alignPositions[i] = pos[alignParticles[i]]-alignCenter;

	// Compute the correlation matrix.

	double R[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < numAlignParticles; k++) {
				int index = alignParticles[k];
				R[i][j] += alignPositions[k][i]*referencePositions[index][j];
			}

	// Compute the F matrix.

	Array2D<double> F(4, 4);
	F[0][0] =  R[0][0] + R[1][1] + R[2][2];
	F[1][0] =  R[1][2] - R[2][1];
	F[2][0] =  R[2][0] - R[0][2];
	F[3][0] =  R[0][1] - R[1][0];

	F[0][1] =  R[1][2] - R[2][1];
	F[1][1] =  R[0][0] - R[1][1] - R[2][2];
	F[2][1] =  R[0][1] + R[1][0];
	F[3][1] =  R[0][2] + R[2][0];

	F[0][2] =  R[2][0] - R[0][2];
	F[1][2] =  R[0][1] + R[1][0];
	F[2][2] = -R[0][0] + R[1][1] - R[2][2];
	F[3][2] =  R[1][2] + R[2][1];

	F[0][3] =  R[0][1] - R[1][0];
	F[1][3] =  R[0][2] + R[2][0];
	F[2][3] =  R[1][2] + R[2][1];
	F[3][3] = -R[0][0] - R[1][1] + R[2][2];

	// Find the maximum eigenvalue and eigenvector.

	JAMA::Eigenvalue<double> eigen(F);
	Array1D<double> values;
	eigen.getRealEigenvalues(values);
	Array2D<double> vectors;
	eigen.getV(vectors);

	// Compute the RMSD.

	int numRMSDParticles = rmsdParticles.size();
	vector<Vec3> rmsdPositions(numRMSDParticles);
	for (int i = 0; i < numRMSDParticles; i++)
		rmsdPositions[i] = pos[rmsdParticles[i]]-alignCenter;

	double sum = 0.0;
	for (int i = 0; i < numRMSDParticles; i++) {
		int index = rmsdParticles[i];
		sum += rmsdPositions[i].dot(rmsdPositions[i]) + referencePositions[index].dot(referencePositions[index]);
	    cout << "sum:" << sum << "\n";
	}
	cout << "numRMSDParticles:" << numRMSDParticles << "\n";
	cout << "values[3]:" << values[3] << "\n";
	double msd = (sum-2*values[3])/numRMSDParticles;
	if (msd < 1e-20) {
		// The particles are perfectly aligned, so all the forces should be zero.
		// Numerical error can lead to NaNs, so just return 0 now.
		cout << "msd too small.\n";
		return 0.0;
	}
	double rmsd = sqrt(msd);

	// Compute the rotation matrix.

	double q[] = {vectors[0][3], vectors[1][3], vectors[2][3], vectors[3][3]};
	double q00 = q[0]*q[0], q01 = q[0]*q[1], q02 = q[0]*q[2], q03 = q[0]*q[3];
	double q11 = q[1]*q[1], q12 = q[1]*q[2], q13 = q[1]*q[3];
	double q22 = q[2]*q[2], q23 = q[2]*q[3];
	double q33 = q[3]*q[3];
	double U[3][3] = {{q00+q11-q22-q33, 2*(q12-q03), 2*(q13+q02)},
					  {2*(q12+q03), q00-q11+q22-q33, 2*(q23-q01)},
					  {2*(q13-q02), 2*(q23+q01), q00-q11-q22+q33}};

	// Rotate the reference positions and compute the forces.

	for (int i = 0; i < numRMSDParticles; i++) {
		const Vec3& p = referencePositions[rmsdParticles[i]];
		Vec3 rotatedRef(U[0][0]*p[0] + U[1][0]*p[1] + U[2][0]*p[2],
						U[0][1]*p[0] + U[1][1]*p[1] + U[2][1]*p[2],
						U[0][2]*p[0] + U[1][2]*p[1] + U[2][2]*p[2]);
		forces[rmsdParticles[i]] -= (rmsdPositions[i]-rotatedRef) / (rmsd*numRMSDParticles);
	}
	cout << "rmsd:" << rmsd << "\n";
	return rmsd;
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
