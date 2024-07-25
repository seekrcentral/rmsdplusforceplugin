// This file contains kernels to compute the RMSD and its gradient using the algorithm described
// in Coutsias et al, "Using quaternions to calculate RMSD" (doi: 10.1002/jcc.20110).

/**
 * Sum a value over all threads.
 */
DEVICE real reduceValue(real value, LOCAL_ARG volatile real* temp) {
    const int thread = LOCAL_ID;
    SYNC_THREADS;
    temp[thread] = value;
    SYNC_THREADS;
    for (int step = 1; step < 32; step *= 2) {
        if (thread+step < LOCAL_SIZE && thread%(2*step) == 0)
            temp[thread] = temp[thread] + temp[thread+step];
        SYNC_WARPS;
    }
    for (int step = 32; step < LOCAL_SIZE; step *= 2) {
        if (thread+step < LOCAL_SIZE && thread%(2*step) == 0)
            temp[thread] = temp[thread] + temp[thread+step];
        SYNC_THREADS;
    }
    return temp[0];
}

/**
 * Perform the first step of computing the RMSD.  This is executed as a single work group.
 */
KERNEL void computeRMSDPart1(int numAlignParticles, 
        GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT referencePos,
        GLOBAL const int* RESTRICT alignParticles, GLOBAL real* buffer) {
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];

    // Compute the center of the align particle positions.

    real3 center = make_real3(0);
    for (int i = LOCAL_ID; i < numAlignParticles; i += LOCAL_SIZE)
        center += trimTo3(posq[alignParticles[i]]);
    center.x = reduceValue(center.x, temp)/numAlignParticles;
    center.y = reduceValue(center.y, temp)/numAlignParticles;
    center.z = reduceValue(center.z, temp)/numAlignParticles;

    // Compute the correlation matrix.

    real R[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    real sum = 0;
    for (int i = LOCAL_ID; i < numAlignParticles; i += LOCAL_SIZE) {
        int index = alignParticles[i];
        real3 pos = trimTo3(posq[index]) - center;
        real3 refPos = trimTo3(referencePos[index]);
        R[0][0] += pos.x*refPos.x;
        R[0][1] += pos.x*refPos.y;
        R[0][2] += pos.x*refPos.z;
        R[1][0] += pos.y*refPos.x;
        R[1][1] += pos.y*refPos.y;
        R[1][2] += pos.y*refPos.z;
        R[2][0] += pos.z*refPos.x;
        R[2][1] += pos.z*refPos.y;
        R[2][2] += pos.z*refPos.z;
        sum += dot(pos, pos);
    }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            R[i][j] = reduceValue(R[i][j], temp);
    sum = reduceValue(sum, temp);
    
    // Copy everything into the output buffer to send back to the host.

    if (LOCAL_ID == 0) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                buffer[3*i+j] = R[i][j];
        buffer[9] = sum;
        buffer[10] = center.x;
        buffer[11] = center.y;
        buffer[12] = center.z;
    }
}

/**
 * Compute the RMSD of the rmsdParticles
 */
KERNEL void computeRMSDPart2(int numRMSDParticles, 
        GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT referencePos,
        GLOBAL const int* RESTRICT rmsdParticles, GLOBAL real* buffer) {
    LOCAL volatile real temp[THREAD_BLOCK_SIZE];
    
    real3 center = make_real3(buffer[10], buffer[11], buffer[12]);
    real sum = 0;
    for (int i = LOCAL_ID; i < numRMSDParticles; i += LOCAL_SIZE) {
        int index = rmsdParticles[i];
        real3 pos = trimTo3(posq[index]) - center;
        real3 refPos = trimTo3(referencePos[index]);
        sum += (pos.x - refPos.x)*(pos.x - refPos.x);
        sum += (pos.y - refPos.y)*(pos.y - refPos.y);
        sum += (pos.z - refPos.z)*(pos.z - refPos.z);
    }
    sum = reduceValue(sum, temp);
    
    if (LOCAL_ID == 0) {
        double msd = sum / numRMSDParticles;
        double rmsd = SQRT(msd);
        buffer[9] = rmsd;
    }
}

/**
 * Apply forces based on the RMSD.
 */
KERNEL void computeRMSDForces(int numRMSDParticles, int paddedNumAtoms, 
        GLOBAL const real4* RESTRICT posq, GLOBAL const real4* RESTRICT referencePos,
        GLOBAL const int* RESTRICT rmsdParticles, GLOBAL const real* buffer, GLOBAL mm_long* RESTRICT forceBuffers) {
    real3 center = make_real3(buffer[10], buffer[11], buffer[12]);
    real scale = 1 / (real) (buffer[9]*numRMSDParticles);
    for (int i = GLOBAL_ID; i < numRMSDParticles; i += GLOBAL_SIZE) {
        int index = rmsdParticles[i];
        real3 pos = trimTo3(posq[index]) - center;
        real3 refPos = trimTo3(referencePos[index]);
        real3 rotatedRef = make_real3(buffer[0]*refPos.x + buffer[3]*refPos.y + buffer[6]*refPos.z,
                                      buffer[1]*refPos.x + buffer[4]*refPos.y + buffer[7]*refPos.z,
                                      buffer[2]*refPos.x + buffer[5]*refPos.y + buffer[8]*refPos.z);
        real3 force = (rotatedRef-pos)*scale;
        forceBuffers[index] += realToFixedPoint(force.x);
        forceBuffers[index+paddedNumAtoms] += realToFixedPoint(force.y);
        forceBuffers[index+2*paddedNumAtoms] += realToFixedPoint(force.z);
    }
}
