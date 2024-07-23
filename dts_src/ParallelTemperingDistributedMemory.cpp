/*
Author: Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
Copyright (c) Weria Pezeshkian

Description:
    This class checks the name of the executable, although it's not a critical task for the current version.
*/

#ifdef MPI_DETECTED
#include <mpi.h>


#include <vector>
#include <string>
#include "SimDef.h"
#include "ParallelTemperingDistributedMemory.h"
#include "RNG.h"
#include "State.h"

/*
Description:
    This class handles task distribution and allows for execution.

Parameters:
    argument (std::vector<std::string>&): Vector containing input arguments.

*/
ParallelTemperingDistributedMemory::ParallelTemperingDistributedMemory(std::vector<std::string> Argument) : m_Argument(Argument) {


}
ParallelTemperingDistributedMemory::~ParallelTemperingDistributedMemory() {
    
}
bool ParallelTemperingDistributedMemory::Initialize(ParallelReplicaData PRD){


    std::vector<std::string> data = Nfunction::Split(PRD.Data);
    
    if(data.size() < 4){
        return false;
    }
    m_Rate = Nfunction::String_to_Int(data[0]);
    m_Bins = Nfunction::String_to_Int(data[1]);
    m_minBeta = Nfunction::String_to_Double(data[2]);
    m_maxBeta = Nfunction::String_to_Double(data[3]);
    
    return true;
}
bool ParallelTemperingDistributedMemory::Run() {

    //Comparison of number of threads asked for (m_Bins) and number of threads available
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout<<"Rank: "<<rank<<std::endl;

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout<<"Size: "<<size<<std::endl;

    
    
    State ReplicaState(m_Argument);
        //--> set the temprature
    double beta = m_minBeta + double(rank) * (m_maxBeta - m_minBeta)/double(m_Bins-1);
    ReplicaState.GetSimulation()->SetBeta(beta, 0);
    //--> set the run tag id, we need to update this ID each time that the processor changes its temprature. The id should be temprature dependent
        //std::string gfile = ReplicaState.GetRunTag() + Nfunction::Int_to_String(beta); // general output file name
        //ReplicaState.UpdateRunTag(gfile);
    //ReplicaState.GetNonbinaryTrajectory()->SetFolderName(Nfunction::Int_to_String(rank));
    ReplicaState.Initialize();
        

        //Learn how to access tsi folder name and VTU folder!



   // T_state.GetVisualization()
    // T_state.GetSimulation()->UpdateInitialStep(int ini_step)
   // T_state.GetSimulation()->UpdateFinalStep(int final_step)
    
   // T_state.GetVisualization() = new WritevtuFiles(&T_state, period, foldername);
   // 
    
    ReplicaState.GetSimulation()->do_Simulation();
    //return true;}

    //else{
    //    std::cout<<"I am not doing anything for rank: "<<rank<<std::endl;
    //    return true;
    //}
    std::cout<<"I am done for rank: "<<rank<<std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
}

std::string ParallelTemperingDistributedMemory::CurrentState(){
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName() + " "+ Nfunction::D2S(m_Rate)+" "+Nfunction::D2S(m_Bins);
        state = state +" "+Nfunction::D2S(m_minBeta) +" "+Nfunction::D2S(m_maxBeta);
        return state;
}


#endif

