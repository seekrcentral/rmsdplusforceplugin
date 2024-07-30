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

#include "CommonRMSDPlusKernels.h"
#include "CommonRMSDPlusKernelSources.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeForceInfo.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/internal/ContextImpl.h"
#include "jama_eig.h"
#include <set>

using namespace RMSDPlusForcePlugin;
using namespace OpenMM;
using namespace std;


class CommonCalcRMSDPlusForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const RMSDPlusForce& force) : force(force) {
        updateBothParticles();
    }
    void updateBothParticles() {
        alignParticles.clear();
        for (int i=0;  i < force.getNumAlignParticles(); i++) {
        	int alignParticle;
        	force.getRMSDPlusAlignParameters(i, alignParticle);
            alignParticles.insert(alignParticle);
        }
        for (int i=0; i < force.getNumRMSDParticles(); i++) {
			int rmsdParticle;
			force.getRMSDPlusRMSDParameters(i, rmsdParticle);
			rmsdParticles.insert(rmsdParticle);
		}
    }
    bool areAlignParticlesIdentical(int particle1, int particle2) {
        bool include1 = (alignParticles.find(particle1) != alignParticles.end());
        bool include2 = (alignParticles.find(particle2) != alignParticles.end());
        return (include1 == include2);
    }
private:
    const RMSDPlusForce& force;
    set<int> alignParticles;
    set<int> rmsdParticles;
};

void CommonCalcRMSDPlusForceKernel::initialize(const System& system, const RMSDPlusForce& force) {
    // Create data structures.

	ContextSelector selector(cc);
	bool useDouble = cc.getUseDoublePrecision();
	int elementSize = (useDouble ? sizeof(double) : sizeof(float));
	int numAlignParticles = force.getNumAlignParticles();
	int numRMSDParticles = force.getNumRMSDParticles();
	int numReferencePositions = force.getNumReferencePositions();
	if (numAlignParticles == 0)
		throw OpenMMException("alignParticles may not be empty.");

	if (numRMSDParticles == 0)
		throw OpenMMException("rmsdParticles may not be empty.");

	if (numReferencePositions == 0)
		throw OpenMMException("referencePositions may not be empty.");

	referencePos.initialize(cc, system.getNumParticles(), 4*elementSize, "referencePos");
	alignParticles.initialize<int>(cc, numAlignParticles, "alignParticles");
	rmsdParticles.initialize<int>(cc, numRMSDParticles, "rmsdParticles");
	buffer.initialize(cc, 13, elementSize, "buffer");
	recordParameters(force);
    info = new ForceInfo(force);
	cc.addForce(info);

	// Create the kernels.
	blockSize = min(256, cc.getMaxThreadBlockSize());
	map<string, string> defines;
	defines["THREAD_BLOCK_SIZE"] = cc.intToString(blockSize);
	ComputeProgram program = cc.compileProgram(CommonRMSDPlusKernelSources::RMSDPlusForce, defines);
	kernel1 = program->createKernel("computeRMSDPart1");
	kernel2 = program->createKernel("computeRMSDPart2");
	kernel3 = program->createKernel("computeRMSDForces");
	
	kernel1->addArg();
	kernel1->addArg(cc.getPosq());
	kernel1->addArg(referencePos);
	kernel1->addArg(alignParticles);
	kernel1->addArg(buffer);
	
	kernel2->addArg();
    kernel2->addArg(cc.getPosq());
    kernel2->addArg(referencePos);
    kernel2->addArg(rmsdParticles);
    kernel2->addArg(buffer);
	
	kernel3->addArg();
	kernel3->addArg(cc.getPaddedNumAtoms());
	kernel3->addArg(cc.getPosq());
	kernel3->addArg(referencePos);
	kernel3->addArg(rmsdParticles);
	kernel3->addArg(buffer);
	kernel3->addArg(cc.getLongForceBuffer());
}

void CommonCalcRMSDPlusForceKernel::recordParameters(const RMSDPlusForce& force) {
    // Record the parameters and center the reference positions.

	vector<int> alignParticleVec;
	int numAlignParticles = force.getNumAlignParticles();
	if (numAlignParticles == 0)
		throw OpenMMException("alignParticles may not be empty.");
	alignParticleVec.resize(numAlignParticles);
	for (int i = 0; i < numAlignParticles; i++)
	    force.getRMSDPlusAlignParameters(i, alignParticleVec[i]);

	vector<int> rmsdParticleVec;
	int numRMSDParticles = force.getNumRMSDParticles();
	if (numRMSDParticles == 0)
		throw OpenMMException("rmsdParticles may not be empty.");
	rmsdParticleVec.resize(numRMSDParticles);
	for (int i = 0; i < numRMSDParticles; i++)
		force.getRMSDPlusRMSDParameters(i, rmsdParticleVec[i]);

	vector<Vec3> centeredPositions;
	int numReferencePositions = force.getNumReferencePositions();
	if (numReferencePositions != system.getNumParticles())
		throw OpenMMException("referencePositions must have the same number of particles as the system.");
	centeredPositions.resize(numReferencePositions);
	for (int i = 0; i < numReferencePositions; i++)
			force.getRMSDPlusReferencePosition(i, centeredPositions[i]);

    Vec3 center;
    for (int i : alignParticleVec)
        center += centeredPositions[i];
    center /= alignParticleVec.size();
    for (Vec3& p : centeredPositions)
        p -= center;

    // Upload them to the device.

    alignParticles.upload(alignParticleVec);
    rmsdParticles.upload(rmsdParticleVec);
    vector<mm_double4> pos;
    for (Vec3 p : centeredPositions)
        pos.push_back(mm_double4(p[0], p[1], p[2], 0));
    referencePos.upload(pos, true);

    // Record the sum of the norms of the reference positions.

    sumNormRef = 0.0;
    for (int i : alignParticleVec) {
        Vec3 p = centeredPositions[i];
        sumNormRef += p.dot(p);
    }
}

double CommonCalcRMSDPlusForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    if (cc.getUseDoublePrecision())
        return executeImpl<double>(context);
    return executeImpl<float>(context);
}

template <class REAL>
double CommonCalcRMSDPlusForceKernel::executeImpl(OpenMM::ContextImpl& context) {
    // Execute the first kernel.

    int numAlignParticles = alignParticles.getSize();
    int numRMSDParticles = rmsdParticles.getSize();
    kernel1->setArg(0, numAlignParticles);
    kernel1->execute(blockSize, blockSize);
    
    // Download the results, build the F matrix, and find the maximum eigenvalue
    // and eigenvector.

    vector<REAL> b;
    buffer.download(b);

    // JAMA::Eigenvalue may run into an infinite loop if we have any NaN
    for (int i = 0; i < 9; i++) {
        if (b[i] != b[i])
            throw OpenMMException("NaN encountered during RMSDPlus force calculation");
    }
    
    Array2D<double> F(4, 4);
    F[0][0] =  b[0*3+0] + b[1*3+1] + b[2*3+2];
    F[1][0] =  b[1*3+2] - b[2*3+1];
    F[2][0] =  b[2*3+0] - b[0*3+2];
    F[3][0] =  b[0*3+1] - b[1*3+0];
    F[0][1] =  b[1*3+2] - b[2*3+1];
    F[1][1] =  b[0*3+0] - b[1*3+1] - b[2*3+2];
    F[2][1] =  b[0*3+1] + b[1*3+0];
    F[3][1] =  b[0*3+2] + b[2*3+0];
    F[0][2] =  b[2*3+0] - b[0*3+2];
    F[1][2] =  b[0*3+1] + b[1*3+0];
    F[2][2] = -b[0*3+0] + b[1*3+1] - b[2*3+2];
    F[3][2] =  b[1*3+2] + b[2*3+1];
    F[0][3] =  b[0*3+1] - b[1*3+0];
    F[1][3] =  b[0*3+2] + b[2*3+0];
    F[2][3] =  b[1*3+2] + b[2*3+1];
    F[3][3] = -b[0*3+0] - b[1*3+1] + b[2*3+2];
    JAMA::Eigenvalue<double> eigen(F);
    Array1D<double> values;
    eigen.getRealEigenvalues(values);
    Array2D<double> vectors;
    eigen.getV(vectors);

    // Compute the RMSD.
    
    double align_msd = (sumNormRef+b[9]-2*values[3])/numAlignParticles;
    if (align_msd < 1e-20) {
        // The particles are perfectly aligned, so all the forces should be zero.
        // Numerical error can lead to NaNs, so just return 0 now.
        return 0.0;
    }
    
    //double RMSDCV = sqrt(msd);
    //b[9] = RMSDCV;

    // Compute the rotation matrix.

    double q[] = {vectors[0][3], vectors[1][3], vectors[2][3], vectors[3][3]};
    double q00 = q[0]*q[0], q01 = q[0]*q[1], q02 = q[0]*q[2], q03 = q[0]*q[3];
    double q11 = q[1]*q[1], q12 = q[1]*q[2], q13 = q[1]*q[3];
    double q22 = q[2]*q[2], q23 = q[2]*q[3];
    double q33 = q[3]*q[3];
    b[0] = q00+q11-q22-q33;
    b[1] = 2*(q12-q03);
    b[2] = 2*(q13+q02);
    b[3] = 2*(q12+q03);
    b[4] = q00-q11+q22-q33;
    b[5] = 2*(q23-q01);
    b[6] = 2*(q13-q02);
    b[7] = 2*(q23+q01);
    b[8] = q00-q11-q22+q33;

    // Upload it to the device and invoke the kernel to apply forces.
    
    buffer.upload(b);

    // Create a kernel that will compute the RMSD of the rmsdParticles
    kernel2->setArg(0, numRMSDParticles);
    kernel2->execute(numRMSDParticles);
    buffer.download(b);
    double RMSD = b[9];

    kernel3->setArg(0, numRMSDParticles);
    kernel3->execute(numRMSDParticles);
    return RMSD;
}

//TO DO: clean up and refer back to record Parameters
void CommonCalcRMSDPlusForceKernel::copyParametersToContext(ContextImpl& context, const RMSDPlusForce& force) {
    ContextSelector selector(cc);
    vector<int> alignParticleVec;
    int numAlignParticles = force.getNumAlignParticles();
    if (numAlignParticles == 0)
        throw OpenMMException("alignParticles may not be empty.");
    alignParticleVec.resize(numAlignParticles);
    for (int i = 0; i < numAlignParticles; i++)
        force.getRMSDPlusAlignParameters(i, alignParticleVec[i]);
    
    vector<int> rmsdParticleVec;
    int numRMSDParticles = force.getNumRMSDParticles();
    if (numRMSDParticles == 0)
        throw OpenMMException("rmsdParticles may not be empty.");
    rmsdParticleVec.resize(numRMSDParticles);
    for (int i = 0; i < numRMSDParticles; i++)
        force.getRMSDPlusRMSDParameters(i, rmsdParticleVec[i]);

    vector<Vec3> centeredPositions;
    int numReferencePositions = force.getNumReferencePositions();
    if (numReferencePositions != system.getNumParticles())
        throw OpenMMException("referencePositions must have the same number of particles as the system.");
    centeredPositions.resize(numReferencePositions);
    for (int i = 0; i < numReferencePositions; i++)
            force.getRMSDPlusReferencePosition(i, centeredPositions[i]);

    Vec3 center;
    for (int i : alignParticleVec)
        center += centeredPositions[i];
    center /= alignParticleVec.size();
    for (Vec3& p : centeredPositions)
        p -= center;
    
    alignParticles.upload(alignParticleVec);
    rmsdParticles.upload(rmsdParticleVec);
    vector<mm_double4> pos;
    for (Vec3 p : centeredPositions)
        pos.push_back(mm_double4(p[0], p[1], p[2], 0));
    referencePos.upload(pos, true);

    // Record the sum of the norms of the reference positions.

    sumNormRef = 0.0;
    for (int i : alignParticleVec) {
        Vec3 p = centeredPositions[i];
        sumNormRef += p.dot(p);
    }

    info->updateBothParticles();
    cc.invalidateMolecules(info);

}

