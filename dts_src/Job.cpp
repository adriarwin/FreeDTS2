/*
Author: Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
Copyright (c) Weria Pezeshkian

Description:
    This class checks the name of the executable, although it's not a critical task for the current version.
*/

#ifdef _OPENMP
    #include <omp.h>
#endif

#ifdef MPI_VERSION
    #include <mpi.h> // Include the MPI header file
#endif

#include <vector>
#include <string>
#include "SimDef.h"
#include "Job.h"
#include "State.h"
#include "RNG.h"
#include "AbstractParallelTemperingSharedMemory.h"
#include "AbstractParallelTemperingDistributedMemory.h"
#include "ParallelTemperingSharedMemory.h"
#include  "ParallelTemperingDistributedMemory.h"


/*
Description:
    This class handles task distribution and allows for execution.

Parameters:
    argument (std::vector<std::string>&): Vector containing input arguments.

*/
Job::Job(const std::vector<std::string> &argument) {
    // Extract executable name from the argument list
    std::string ex_name = Nfunction::SubstringFromRight(argument[0], '/');
    
    // Check if the executable name matches the expected name
    if (ex_name != EXE_NAME) { // EXE_NAME is defined in the SimDef.h file
        std::cout << "--> unrecognized executable name ---> " << ex_name << " :( " << " it should be " << EXE_NAME << std::endl;
        exit(0);
    }
    #if !defined(_OPENMP) && !defined(MPI_VERSION)
        std::cout << "Running simulation on single CPU" << std::endl;
    
        State T_state(argument);
        T_state.Initialize();
        T_state.GetSimulation()->do_Simulation();
    #elif defined(_OPENMP)
        // Code for when either _OPENMP
    // Perform a normal simulation on a single CPU if OpenMP is not enabled
//---> constract an State object
    std::cout<<"OpenMP has been detected. Initializing parallel tempering routine"<<std::endl;
    State T_state(argument);
//---> get parallel tempering data in the input.dts file
    ParallelReplicaData    PRD = T_state.GetParallelReplicaData();
    //---> To be added. Depending on the type of PRD, we will run with either parallel_tempering with shared memory
    // or parallel_tempering with distributed memory.
    //Basically, take into account the type...
    
    
    
//---> here is one openmp is on but still want to perform one single simulation
    if (!PRD.State) {
        std::cout<<"OpenMP has been detected. We run on a single CPU."<<std::endl;
        T_state.Initialize();
        T_state.GetSimulation()->do_Simulation();
    }
else { // run parallel tempering simulations
    std::cout<<"OpenMP has been detected. Initializing parallel tempering routine."<<std::endl;
    AbstractParallelTemperingSharedMemory *pParallelReplicaRun;
    
        
    if (PRD.Type == ParallelTemperingSharedMemory::GetDefaultReadName()){
            
        pParallelReplicaRun = new ParallelTemperingSharedMemory(argument);
        if(pParallelReplicaRun->Initialize(PRD)){
            pParallelReplicaRun->Run();
        }
        else{
            std::cout<<"---> error: faild .... "<<"\n";
        }
    }
    else{
        std::cout<<"---> error: unknow type for "<<AbstractParallelTemperingSharedMemory::GetBaseDefaultReadName()<<"\n";
        exit(0);
    }
}

#elif defined(MPI_VERSION)
std::cout<<"MPI has been detected. Initializing parallel tempering routine."<<std::endl;
#endif
}
Job::~Job() {
    
}




