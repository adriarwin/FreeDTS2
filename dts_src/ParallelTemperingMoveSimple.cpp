

#include <stdio.h>
#include "ParallelTemperingMoveSimple.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "Voxelization.h"
#ifdef MPI_DETECTED
#include <mpi.h>
#endif


/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "ChangeBoxSize" function.
=================================================================================================================
*/
ParallelTemperingMoveSimple::ParallelTemperingMoveSimple(State *pState,int period, int nprocessors, double minbeta, double maxbeta) :
        m_pState(pState),
        m_Period(period),
        m_Nprocessors (nprocessors),
        m_MinBeta (minbeta),
        m_MaxBeta (maxbeta)
        {
        }
            //-- convert the type string into the direction vector

ParallelTemperingMoveSimple::~ParallelTemperingMoveSimple() {
    
}
void ParallelTemperingMoveSimple::Initialize() {
    
    std::cout<<"---> the algorithm for Parallel Tempering involves applying this: "<< GetBaseDefaultReadName()<<" \n";
    #ifdef MPI_DETECTED
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout<<"Rank: "<<rank<<std::endl;
    m_Rank=rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout<<"Size: "<<size<<std::endl;
    m_Size=size;
    m_tempID=rank;
    double beta = m_MinBeta + double(rank) * (m_MaxBeta - m_MinBeta)/double(m_Nprocessors-1);
    m_pState->GetSimulation()->SetBeta(beta, 0);
    //--> set the run tag id, we need to update this ID each time that the processor changes its temprature. The id should be temprature dependent
        //std::string gfile = ReplicaState.GetRunTag() + Nfunction::Int_to_String(beta); // general output file name
        //ReplicaState.UpdateRunTag(gfile);

    //In here, I need to think how to implement the restart in the right way (will be hard but it is extremly important to do it right
    m_pState->GetNonbinaryTrajectory()->SetFolderName( m_pState->GetNonbinaryTrajectory()->GetFolderName() +"_" + Nfunction::Int_to_String(rank));
    m_pState->GetTimeSeriesDataOutput()->SetCustomFileName( m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(rank)+TimeSeriDataExt);

    #endif

}
bool ParallelTemperingMoveSimple::EvolveOneStep(int step){
    /**
     * @brief a call function to change the simulation box size at a given step.
     *
     * This function changes the size of the simulation box based on the current step.
     * It first checks if the current step matches the defined period. If the voxel size
     * is below a threshold, it updates the voxel size and re-voxelizes. It then computes
     * the size change for the box using, and attempts to change the
     * box size based on these calculations.
     *
     * @param step The current step in the simulation.
     * @return true if the box size was changed, false otherwise.
     */
//---> if does not match the preiod, return false
    if(step%m_Period != 0)
        return false;
    //--- first check if the voxel size is fine
    //Now we do everything needed to change a step!

    
    return true;
}

std::string ParallelTemperingMoveSimple::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    //What is this used for?
    //state = state +" "+ Nfunction::D2S(m_Period)+" "+Nfunction::D2S(m_SigmaP)+" "+  m_Type;
    return state;
}




