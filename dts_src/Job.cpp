/*
Author: Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
Copyright (c) Weria Pezeshkian

Description:
    This class checks the name of the executable, although it's not a critical task for the current version.
*/

#ifdef _OPENMP
    #include <omp.h>
#endif

#ifdef MPI_DETECTED
    #include <mpi.h> 
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
//Just a little test
//Second test
Job::Job(const std::vector<std::string> &argument) {
    // Extract executable name from the argument list
    
    
    // Check if the executable name matches the expected name

    #if !defined(_OPENMP) && !defined(MPI_DETECTED)

    std::string ex_name = Nfunction::SubstringFromRight(argument[0], '/');
    if (ex_name != EXE_NAME) { // EXE_NAME is defined in the SimDef.h file
        std::cout << "--> unrecognized executable name ---> " << ex_name << " :( " << " it should be " << EXE_NAME << std::endl;
        exit(0);
    }

    std::cout << "MPI flag not detected. " << std::endl;
    State T_state(argument);
    T_state.Initialize();
    T_state.GetSimulation()->do_Simulation();
    
    #elif defined(_OPENMP)

    std::string ex_name = Nfunction::SubstringFromRight(argument[0], '/');
    if (ex_name != EXE_NAME) { // EXE_NAME is defined in the SimDef.h file
        std::cout << "--> unrecognized executable name ---> " << ex_name << " :( " << " it should be " << EXE_NAME << std::endl;
        exit(0);
    }

    std::cout<<"OpenMP is detected in Job.cpp"<<std::endl;
    State T_state(argument);
//---> get parallel tempering data in the input.dts file
    ParallelReplicaData    PRD = T_state.GetParallelReplicaData();
    //---> To be added. Depending on the type of PRD, we will run with either parallel_tempering with shared memory
    // or parallel_tempering with distributed memory.
    //Basically, take into account the type...
    
    
    
//---> here is one openmp is on but still want to perform one single simulation
    if (!PRD.State) {

    std::cout<<"OpenMP has been detected but we run on a single CPU as stated in input file"<<std::endl;
    T_state.Initialize();
    T_state.GetSimulation()->do_Simulation();

    }
    else{ // run parallel tempering simulations
    std::cout<<"OpenMP has been detected and we initialize parallel tempering routine with shared memory"<<std::endl;
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

    #elif defined(MPI_DETECTED)
    std::cout << "MPI flag detected. " << std::endl;
    State T_state(argument);
    T_state.Initialize();
    //T_state.GetSimulation()->do_Simulation();
    
/*
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //Only run this with one CPU
    

    State T_state(argument);
//---> get parallel tempering data in the input.dts file
    ParallelReplicaData    PRD = T_state.GetParallelReplicaData();

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout<<"Rank: "<<world_rank<<", PRD.Type: "<<PRD.Type<<std::endl;

    std::string ex_name = Nfunction::SubstringFromRight(argument[0], '/');
    if (ex_name != EXE_NAME) { // EXE_NAME is defined in the SimDef.h file
        std::cout << "--> unrecognized executable name ---> " << ex_name << " :( " << " it should be " << EXE_NAME << std::endl;
        exit(0);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank == 0) {
    std::cout<<"MPI has been detected in Job.cpp for all ranks"<<std::endl;}

    if (!PRD.State) {

    if (world_rank == 0) {
    std::cout<<"MPI has been detected but we run on a single CPU as stated in input file"<<std::endl;
    T_state.Initialize();
    T_state.GetSimulation()->do_Simulation();}

    }
    else{ // run parallel tempering simulations
    std::cout<<"We got here all of us!"<<std::endl;
    if (world_rank == 0) {
    std::cout<<"MPI has been detected and we initialize parallel tempering routine with distributed memory"<<std::endl;}
    AbstractParallelTemperingDistributedMemory *pParallelReplicaRun;
    
    std::cout<<"ParallelTemperingDistributedMemory::GetDefaultReadName()"<<ParallelTemperingDistributedMemory::GetDefaultReadName()<<", Rank = "<<world_rank<<std::endl;
    if (PRD.Type == ParallelTemperingDistributedMemory::GetDefaultReadName()){
        if (world_rank == 0) {
        std::cout<<"About to make the call for distributed Memory!"<<std::endl;}

        pParallelReplicaRun = new ParallelTemperingDistributedMemory(argument);
        if(pParallelReplicaRun->Initialize(PRD)){
            if (world_rank == 0) {
            std::cout<<"Parallel tempering routine initialized!"<<std::endl;
            std::cout<<"Running parallel tempering routine with distributed memory"<<std::endl;}
            pParallelReplicaRun->Run();
        //}
        //else{
        //    std::cout<<"---> error: faild .... "<<"\n";
        //}
    }
    else{
        std::cout<<"---> error: unknow type for "<<AbstractParallelTemperingDistributedMemory::GetBaseDefaultReadName()<<"\n";
        exit(0);
    }

    }
    }

*/
#endif
}
Job::~Job() {
    
}




