

#include <stdio.h>
#include "ParallelTemperingMoveSimple.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "Voxelization.h"
#ifdef MPI_DETECTED
#include <mpi.h>



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
    m_TempID=rank;
    m_Counter=0;
    m_SizeIsEven=(m_Size%2==0);
    m_ReceiveBroadcast.resize(m_Size);
    m_RequestBroadcast.resize(m_Size);
    //m_RankWithUpTempID=rank;
    //m_RankWithDownTempID=rank;
    double beta;
    for (int i=0;i<m_Size;i++)
        {
            beta = m_MinBeta + double(i) * (m_MaxBeta - m_MinBeta)/double(m_Nprocessors-1);
            m_BetaVec.push_back(beta);
            m_RankAtTempID.push_back(i);
        }

    m_pState->GetSimulation()->SetBeta(m_BetaVec[m_TempID], 0);
    //--> set the run tag id, we need to update this ID each time that the processor changes its temprature. The id should be temprature dependent
        //std::string gfile = ReplicaState.GetRunTag() + Nfunction::Int_to_String(beta); // general output file name
        //ReplicaState.UpdateRunTag(gfile);

    //In here, I need to think how to implement the restart in the right way (will be hard but it is extremly important to do it right
    m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(m_TempID));
    m_pState->GetTimeSeriesDataOutput()->SetCustomFileName(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(m_TempID)+TimeSeriDataExt);
    
    #endif

}
bool ParallelTemperingMoveSimple::EvolveOneStep(int step){
    /**
     * @brief a call function to exchange temperature between replicas.
     *
     * This function attempts a temperature swap between different replicas.
     * The algorithm is a decentralized one, in the sense that each couple of
     * neighbouring replicas that attempt the exchange do not need of the other ones
     * to continue with the simulation.
     * 
     * O --> Odd rank
     * E --> Even rank
     * N--> Number of processors
     * 
     * If O or E are written, it means that a conditional is imposed for TempID.
     * If your TempID is odd or even, then...
     * @param step The current step in the simulation.
     * @return true if the box size was changed, false otherwise.
     */
//---> if does not match the preiod, return false
    if(step%m_Period != 0)
        return false;
    //--- first check if the voxel size is fine
    //Now we do everything needed to change a step!
    
    bool CountIsEven=(m_Counter%2==0);

    //TEST::
    CountIsEven=true;

    if (m_Counter>0){
        for(int i=0; i<m_Size; i++){
            if(i%2==0){
                MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
            }
        }

        

        for(int i=0; i<m_Size; i++){
            if(i%2==0){
                if (m_ReceiveBroadcast[i]== 1){
                    std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                }
            }
        }
        
        for(int i=0; i<m_Size; i++){
            if(m_RankAtTempID[i]==m_Rank){
                m_TempID=i;
                break;
            }
        }

        std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", RankAtTempID:";
        for (const int& element : m_RankAtTempID) {
            std::cout << element << ' ';
        }
        std::cout << std::endl; 

     


        
    }

    bool TempIDIsEven=(m_TempID%2==0);

    if (CountIsEven){

    if (!TempIDIsEven){

        double sending_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        int DestTempID=m_TempID-1;
        int DestRank=m_RankAtTempID[DestTempID];
        //std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Sending Energy="<<sending_energy<<", Target="<<DestTempID<<std::endl;
        int ReceiveExchangeAcceptedInt;

        

        MPI_Send(&sending_energy, 1, MPI_DOUBLE, DestRank, 0, MPI_COMM_WORLD);
        MPI_Recv(&ReceiveExchangeAcceptedInt, 1, MPI_INT, DestRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        bool ReceiveExchangeAccepted=(ReceiveExchangeAcceptedInt==1);
        //std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Receive Exchange Accepted="<<ReceiveExchangeAccepted<<std::endl;

        //If it is true, first of all prepare the grounds for the exchange.
        //std::vector<int> OldRankAtTempID=m_RankAtTempID;
        if (ReceiveExchangeAccepted){
            int NewTempID=DestTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID)); 
            //m_RankAtTempID[m_TempID]=m_RankAtTempID[NewTempID];
            //m_RankAtTempID[NewTempID]=m_Rank;
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt=0;
            MPI_Send(&EmptyBlockingInt, 1, MPI_INT, DestRank, 3, MPI_COMM_WORLD);
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            //m_TempID=NewTempID;
            }
            
            //Close the file I was writing to regarding energy
            //Open the file I was writing to... (will be a bit hard?)

        for(int i=0;i<m_Size;i++){
            if(i%2==0){
                //std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Broadcasting receive from "<< m_RankAtTempID[i]<<std::endl;
                MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
        }



        /*
        for(int i=0; i<m_Size; i++){
            if(i%2==0){
                MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
            }
        }

        std::cout<<"completed"<<std::endl;
        
        for(int i=0; i<m_Size; i++){
            if(i%2==0){
                if (m_ReceiveBroadcast[i]== 1){
                    std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                }
            }
        }*/
            
        
        //Set a SEND BLOCKING call to m_TempID-1 with sending_energy.
    }

    else{
        //Set a conditional. If m_size%2==0, then it is okay. If m_size%2!=0, then, 
        //do not calculate the energy for the m_TempID==m_size-1 case.
        if (!m_SizeIsEven && m_TempID==m_Size-1){
            //This rank does not need to do anything
            //Before exiting, always add +1 to the counter
            m_Counter++;
            /*std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Exiting."<<std::endl;
            for(int i=0;i<m_Size;i++){
            if(i%2==0){
                MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
            }

            for(int i=0; i<m_Size; i++){
            if(i%2==0){
                MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
            }
            }

            for(   int i=0; i<m_Size; i++){
            if(i%2==0){
                if (m_ReceiveBroadcast[i]== 1){
                    std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                }
            }
            }*/


            return false;
        }

        double own_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        double receive_energy;
        int SourceTempID=m_TempID+1;
        int SourceRank=m_RankAtTempID[SourceTempID];
        MPI_Recv(&receive_energy, 1, MPI_DOUBLE, SourceRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Receiving Energy="<<receive_energy<<", From="<<SourceRank<<std::endl;

        double dE=receive_energy-own_energy;
        double dBeta=m_BetaVec[SourceTempID] - m_BetaVec[m_TempID];
        double RandomNumber=m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        bool ExchangeAccepted=(RandomNumber<exp(dBeta*dE));

        //TEST
        ExchangeAccepted=true;
        // Convert boolean to integer
        int ExchangeAcceptedInt = ExchangeAccepted ? 1 : 0;

        //If it is True, set the ground for the exchange?


        MPI_Send(&ExchangeAcceptedInt, 1, MPI_INT, SourceRank, 1, MPI_COMM_WORLD);
        
        

        //int OldTempID=m_TempID;

        //std::vector<int> OldRankAtTempID=m_RankAtTempID;
        
        if (ExchangeAccepted){
            int NewTempID=SourceTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID));
            //m_RankAtTempID[m_TempID]=m_RankAtTempID[NewTempID];
            //m_RankAtTempID[NewTempID]=m_Rank;
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            //Send
            int EmptyBlockingInt;
            MPI_Recv(&EmptyBlockingInt, 1, MPI_INT, SourceRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            //m_TempID=NewTempID;
            //std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Beta="<<m_pState->GetSimulation()->GetBeta()<<std::endl;
        }
        
        //Now, broadcast to all ranks the result of this...

        //Set a RECEIVING BLOCKING call to m_TempID+1 with receive_energy as buffer.
        for(int i=0;i<m_Size;i++){
            if(i%2==0 && i!=m_TempID){
                //std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Broadcasting receive from "<< m_RankAtTempID[i]<<std::endl;
                MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
            else if(i==m_TempID){
                //std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Broadcasting sending."<<std::endl;
                MPI_Ibcast(&ExchangeAcceptedInt, 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
        }

        m_ReceiveBroadcast[m_TempID]=ExchangeAcceptedInt;

        /*for(int i=0; i<m_Size; i++){
            if(i%2==0){
                MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
            }
        }

        std::cout<<"completed"<<std::endl;
        
        m_ReceiveBroadcast[m_TempID]=ExchangeAcceptedInt;

        for(   int i=0; i<m_Size; i++){
            if(i%2==0){
                if (m_ReceiveBroadcast[i]== 1){
                    std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                }
            }
        }*/


    }

    }
    else{

    if (!TempIDIsEven){

        if (m_SizeIsEven && m_TempID==m_Size-1){
            //This rank does not need to do anything
            //Before exiting, always add +1 to the counter
            m_Counter++;
            std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Exiting."<<std::endl;
            return false;
        }

        double own_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        double receive_energy;
        int SourceTempID=m_TempID+1;
        int SourceRank=m_RankAtTempID[SourceTempID];
        MPI_Recv(&receive_energy, 1, MPI_DOUBLE, SourceRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Receiving Energy="<<receive_energy<<", From="<<SourceRank<<std::endl;

        //Set a RECEIVE BLOCKING call to m_TempID+1 with receive_energy as buffer.
    }

    else{
        if (m_TempID==0){
            m_Counter++;
            std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Exiting."<<std::endl;
            return false;
        }
        //Set a conditional. If m_size%2==0, then it is okay. If m_size%2!=0, then, 
        //do not calculate the energy for the m_TempID==m_size-1 case.
        double sending_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        int DestTempID=m_TempID-1;
        int DestRank=m_RankAtTempID[DestTempID];
        std::cout<<"Rank="<<m_Rank<<", TempID="<<m_TempID<<". Sending Energy="<<sending_energy<<", Target="<<DestTempID<<std::endl;
        MPI_Send(&sending_energy, 1, MPI_DOUBLE, DestRank, 0, MPI_COMM_WORLD);
        //Set a RECEIVING BLOCKING call to m_TempID+1 with receive_energy as buffer.
    }


    }
    

    //O: 
    //E: Attempts exchange. (dE=own_energy - receive_energy, dB=BetaVec[TempID-1]-BetaVec[TempID]. Get a random number **2 and attempt exchange.)

    //O: Blocking receive for the result of the attempt (define AcceptedAttemp and blocking receive it).
    //E: Result of the attempt is a bool. (define AcceptedAttempt). Result must be a blocking sent. 

    //O:
    //E: Non-blocking broadcast of the attempt to all other ranks so that they can reconfigure their m_RankAtTempID vectors. This part might be
    //a bit challenging.

    //O: Non-blocking broadcast receiving. 
    //E: Non-blocking broadcast receiving.

    //O: Continue their execution. Swap a variable by which next time the function is called, all this is repeated but for:
    //E: Continue their execution. Swap a variable by which next time the function is called, all this is repetead but for:

    //O: Blocking recive energy from TempID+1
    //E: Blocking send energy to TempID-1.

    //Repeat but changing one from another. One needs to define one variable for when we have one or other option.

    m_Counter++;
    return true;
}

//double own_energy=m_pState->GetEnergyCalculator()->GetEnergy(); **1
//double other_energy;
//std::cout<<"Energy Correctly Obtained: "<<own_energy<<std::endl;
//m_pState->GetRandomNumberGenerator()->UniformRNG(1.0) **2

bool ParallelTemperingMoveSimple::ChangeToNewTemperatureID(int NewTempID){
    //All that one must do when it is time to change Temperature ID
    //Set New Temp ID
    
    m_TempID=NewTempID;
    //Set New Beta
    m_pState->GetSimulation()->SetBeta(m_BetaVec[m_TempID], 0);
    m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetFolderName() +"_" + Nfunction::Int_to_String(m_TempID));

    //Here, I close the file I was writing to
    //I move to the file I want to write to

    //I have to look into this. Do I have to close the file and open it again?
    //I haveto close the file and open it again...

    //In here I have to look at... what happens if one rank closes one file and opens another.
    //
    //m_pState->GetTimeSeriesDataOutput()->SetCustomFileName(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(m_TempID)+TimeSeriDataExt);
    



    return true;
}

std::string ParallelTemperingMoveSimple::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    //What is this used for?
    //state = state +" "+ Nfunction::D2S(m_Period)+" "+Nfunction::D2S(m_SigmaP)+" "+  m_Type;
    return state;
}




#endif