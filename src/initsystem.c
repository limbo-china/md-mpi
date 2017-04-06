#include "initsystem.h"
#include "potential.h"
#include "lattice.h"
#include "parameter.h"

//初始化模拟体系
System* initSystem(Parameter* para){

	System* sys = malloc(sizeof(System));
	sys->energy = NULL;
	sys->potential = NULL;
	sys->lattice = NULL;
   	sys->space = NULL;
   	sys->cells = NULL;
   	sys->atoms = NULL;
   	sys->atomExchange = NULL;

   	initPotInfo(sys->potential);
   	initLatticeInfo(sys->lattice);


   real3 globalExtent;
   globalExtent[0] = cmd.nx * latticeConstant;
   globalExtent[1] = cmd.ny * latticeConstant;
   globalExtent[2] = cmd.nz * latticeConstant;

   sim->domain = initDecomposition(
      cmd.xproc, cmd.yproc, cmd.zproc, globalExtent);

   sim->boxes = initLinkCells(sim->domain, sim->pot->cutoff);
   sim->atoms = initAtoms(sim->boxes);

   // create lattice with desired temperature and displacement.
   createFccLattice(cmd.nx, cmd.ny, cmd.nz, latticeConstant, sim);
   setTemperature(sim, cmd.temperature);
   randomDisplacements(sim, cmd.initialDelta);

   sim->atomExchange = initAtomHaloExchange(sim->domain, sim->boxes);

   // Forces must be computed before we call the time stepper.
   startTimer(redistributeTimer);
   redistributeAtoms(sim);
   stopTimer(redistributeTimer);

   startTimer(computeForceTimer);
   computeForce(sim);
   stopTimer(computeForceTimer);

   kineticEnergy(sim);

   return sim;
}