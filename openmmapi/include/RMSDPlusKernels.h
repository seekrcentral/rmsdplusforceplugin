#ifndef RMSDPLUSFORCE_KERNELS_H_
#define RMSDPLUSFORCE_KERNELS_H_

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
#include "openmm/KernelImpl.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <string>

namespace RMSDPlusForcePlugin {

/**
 * This kernel is invoked by RMSDPlusForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcRMSDPlusForceKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "CalcRMSDPlusForce";
    }
    CalcRMSDPlusForceKernel(std::string name, const OpenMM::Platform& platform) : OpenMM::KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the RMSDPlusForce this kernel will be used for
     */
    virtual void initialize(const OpenMM::System& system, const RMSDPlusForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the RMSDPlusForce to copy the parameters from
     */
    virtual void copyParametersToContext(OpenMM::ContextImpl& context, const RMSDPlusForce& force) = 0;
};

} // namespace RMSDPlusForcePlugin

#endif /*RMSDPLUSFORCE_KERNELS_H_*/
