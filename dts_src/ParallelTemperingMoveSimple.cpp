#include <stdio.h>
#include "ParallelTemperingMoveSimple.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "Voxelization.h"
#ifdef MPI_DETECTED
#include <mpi.h>
#include <sstream>
#include <iostream>


/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "ChangeBoxSize" function.
=================================================================================================================
*/
ParallelTemperingMoveSimple::ParallelTemperingMoveSimple(State *pState,int period, int nprocessors, double minbeta, double maxbeta,int initialsteps) :
        m_pState(pState),
        m_Period(period),
        m_Nprocessors (nprocessors),
        m_MinBeta (minbeta),
        m_MaxBeta (maxbeta),
        m_InitialSteps(initialsteps)
        {
        }
            //-- convert the type string into the direction vector

ParallelTemperingMoveSimple::~ParallelTemperingMoveSimple() {
    
}
void ParallelTemperingMoveSimple::Initialize() {
    
    
    #ifdef MPI_DETECTED
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //std::cout<<"Rank: "<<rank<<std::endl;
    m_Rank=rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //std::cout<<"Size: "<<size<<std::endl;
    m_Size=size;
    m_TempID=rank;
    m_Counter=0;
    m_SizeIsEven=(m_Size%2==0);
    m_ReceiveBroadcast.resize(m_Size);
    m_RequestBroadcast.resize(m_Size);

    if(m_Rank==0){
    std::cout<<"---> the algorithm for Parallel Tempering involves applying this: "<< GetBaseDefaultReadName()<<" \n";}
    //m_RankWithUpTempID=rank;
    //m_RankWithDownTempID=rank;
    if (m_Restart==false){
    for (int i=0;i<m_Size;i++)
        {
            m_RankAtTempID.push_back(i);
        }
    }
    else{
        m_RankAtTempID=ReadLastLineOutput(GetOutputFileName(),m_Size);
        MPI_Barrier(MPI_COMM_WORLD);
        for(int i=0; i<m_Size; i++){
                if(m_RankAtTempID[i]==m_Rank){
                    m_TempID=i;
                    if (m_TempID==0){m_TargetState=true;}
                    break;
                }
        }
    }

    std::vector<double> localBetaVec;
    if (rank == 0) {
        localBetaVec = ReadTemperatures();
    }

    // Broadcast the size of the vector to all ranks
    int vecSize = (rank == 0) ? localBetaVec.size() : 0;
    MPI_Bcast(&vecSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize the local vector based on the size received
    std::vector<double> globalBetaVec(vecSize);
    if (rank == 0) {
        // Ensure localBetaVec.data() is the correct pointer for broadcasting
        MPI_Bcast(localBetaVec.data(), vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Bcast(globalBetaVec.data(), vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Now each rank has a copy of the vector
    m_BetaVec = globalBetaVec;
    if (rank==0){m_BetaVec=localBetaVec;}

    if (m_Size !=m_BetaVec.size()){
        std::cout<<"The selected temperatures do not match the number of processors."<<std::endl;
        exit(0);
    }

    //Getting the number of PT stepts that are going to be run
    int t_ini=m_pState->GetSimulation()->GetInitialStep();
    int t_fin=m_pState->GetSimulation()->GetFinalStep();
    
    m_NumberOfPTSteps= (t_fin - t_ini + 1)/m_Period;

    if (m_Rank==0){
        std::cout<<"Number of parallel tempering steps expected: "<< m_NumberOfPTSteps -1<< std::endl;
    }
    


    m_pState->GetSimulation()->SetBeta(m_BetaVec[m_TempID], 0);
    //--> set the run tag id, we need to update this ID each time that the processor changes its temprature. The id should be temprature dependent
        //std::string gfile = ReplicaState.GetRunTag() + Nfunction::Int_to_String(beta); // general output file name
        //ReplicaState.UpdateRunTag(gfile);

    //In here, I need to think how to implement the restart in the right way (will be hard but it is extremly important to do it right)
    m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(m_TempID));
    m_pState->GetTimeSeriesDataOutput()->SetCustomFileName(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(m_TempID)+TimeSeriDataExt);
    m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(m_TempID));
    
    if (m_Restart==false){
        if (m_Rank == 0) {
            m_TimeSeriesFile.open(GetOutputFileName(), std::ios_base::out);
            m_TimeSeriesFile << "Rank-Temperature ID Mapping: "<<std::endl;
            // Call the new function to write m_RankAtTempID to m_TimeSeriesFile
        }
    }
    else{
        if (m_Rank == 0) {
            m_TimeSeriesFile.open(GetOutputFileName(), std::ios_base::app);
            // Call the new function to write m_RankAtTempID to m_TimeSeriesFile
        }

    }


    #endif
}

void ParallelTemperingMoveSimple::WriteRankAtTempIDToFile() {
    if (!m_TimeSeriesFile.is_open()) {
        std::cerr << "Error: Time series file is not open!" << std::endl;
        return;
    }

    for (const auto& rankID : m_RankAtTempID) {
        m_TimeSeriesFile << rankID << " ";
    }
    m_TimeSeriesFile << std::endl;
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
     * @param step The current step in the simulation.
     * @return true if the box size was changed, false otherwise.
     */
//---> if does not match the preiod, return false
    if(step%m_Period != 0 || step<m_InitialSteps)
        return false;

    MPI_Barrier(MPI_COMM_WORLD);

    bool CountIsEven=(m_Counter%2==0);

    //TEST::



    //RECEIVE BROADCAST FROM LAST EXCHANGE, SWAP RANKS, UPDATE RANKATTEMPID
    if (m_Counter>0){
        if (!CountIsEven){

            #if DEBUG_MODE_PT==Enabled
            std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Setting a wait"<<std::endl;
            #endif 
            for(int i=0; i<m_Size; i++){
                if(i%2==0){
                    MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
                }
            }
            #if DEBUG_MODE_PT==Enabled
            std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Wait finished"<<std::endl;
            #endif 

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
                    if (m_TempID==0){m_TargetState=true;}
                    break;
                }
            }

            //std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", RankAtTempID:";
            //for (const int& element : m_RankAtTempID) {
            //    std::cout << element << ' ';
            //}
            //std::cout << std::endl; 

        }

        else if (CountIsEven){

            //std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<"Checkpoint 1"<<std::endl;     
            #if DEBUG_MODE_PT==Enabled
            std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Setting a wait"<<std::endl;
            #endif 
            for(int i=0; i<m_Size; i++){
                if (m_SizeIsEven){
                    if(i%1==0 && i!=m_Size-1){
                        MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
                    }
                }

                else if (!m_SizeIsEven){
                    if(i%1==0){
                        MPI_Wait(&m_RequestBroadcast[i], MPI_STATUS_IGNORE);
                    }
                }
            }

            //std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<"Checkpoint 2"<<std::endl;
            #if DEBUG_MODE_PT==Enabled
            std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Wait finished"<<std::endl;
            #endif 

            for(int i=0; i<m_Size; i++){

                if (m_SizeIsEven){
                    if(i%1==0 && i!=m_Size-1){
                        if (m_ReceiveBroadcast[i]== 1){
                            std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                        }
                    }
                }
                else if (!m_SizeIsEven){
                    if(i%1==0){
                        if (m_ReceiveBroadcast[i]== 1){
                            std::swap(m_RankAtTempID[i], m_RankAtTempID[i+1]);
                        }
                    }

                }
            }
            
            for(int i=0; i<m_Size; i++){
                if(m_RankAtTempID[i]==m_Rank){
                    m_TempID=i;
                    if (m_TempID==0){m_TargetState=true;}
                    break;
                }
            }

            //std::cout<<"Counter: "<<m_Counter<<", Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", RankAtTempID:";
            //for (const int& element : m_RankAtTempID) {
            //    std::cout << element << ' ';
            //}
            //std::cout << std::endl; 

        }
        
    }

    //Breaking the last communication to do not let any pending comunication
    if (step/m_Period==m_NumberOfPTSteps){
        return true;
    }
    //Just something
    std::fill(m_ReceiveBroadcast.begin(), m_ReceiveBroadcast.end(), 0);

    std::fill(m_RequestBroadcast.begin(), m_RequestBroadcast.end(), MPI_REQUEST_NULL);

    bool TempIDIsEven=(m_TempID%2==0);

    if (m_Rank==0){
        WriteRankAtTempIDToFile();
    }

    //ATTEMPT EXCHANGE
    //CountIsEven=true;
    if (CountIsEven){

    //IF COUNT IS EVEN, ODD TEMPIDS SEND ENERGY AND RECEIVE RESULT OF ATTEMPTED EXCHANGE
    if (!TempIDIsEven){

        //CALCULATE AND SEND ENERGY
        double sending_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        int DestTempID=m_TempID-1;
        int DestRank=m_RankAtTempID[DestTempID];
        int ReceiveExchangeAcceptedInt;

        MPI_Request request_send;
        MPI_Request request_receive;
        MPI_Request request_send_1;
        MPI_Request request_receive_1;
        
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" ,Rank: " << m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " reached stage prior to send energy" << std::endl;
#endif
        MPI_Isend(&sending_energy, 1, MPI_DOUBLE, DestRank, 0, MPI_COMM_WORLD, &request_send);
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" ,Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " sent energy to" <<DestTempID <<std::endl;
#endif
        //MPI_Wait(&request_send, MPI_STATUS_IGNORE);

        //RECEIVE RESULT OF EXCHANGE ATTEMPTS
        MPI_Irecv(&ReceiveExchangeAcceptedInt, 1, MPI_INT, DestRank, 1, MPI_COMM_WORLD, &request_receive);
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step << " ,Rank: " << m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " received energy from" <<DestRank <<std::endl;
#endif
        MPI_Wait(&request_receive, MPI_STATUS_IGNORE);
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step << " ,Rank: " << m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " Completed wait after receive " <<DestRank <<std::endl;
#endif

        bool ReceiveExchangeAccepted=(ReceiveExchangeAcceptedInt==1);

        //IF ACCEPTED, EFFECTIVELY SWAP RANKS
        if (ReceiveExchangeAccepted){
            int NewTempID=DestTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID)); 
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt_1=0;
            int EmptyBlockingInt_2;
            MPI_Isend(&EmptyBlockingInt_1, 1, MPI_INT, DestRank, 3, MPI_COMM_WORLD,&request_send_1);
#if DEBUG_MODE_PT==Enabled
            std::cout << "step:"<< step <<" ,Rank: " << m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " sent blocking INT to " <<DestTempID <<std::endl;
#endif
            //MPI_Wait(&request_send_1,MPI_STATUS_IGNORE);
            MPI_Irecv(&EmptyBlockingInt_2, 1, MPI_INT, DestRank, 4, MPI_COMM_WORLD,&request_receive_1);
#if DEBUG_MODE_PT==Enabled
            std::cout << "step:"<< step <<" ,Rank: " << m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " Waiting " <<DestRank <<std::endl;
#endif
            MPI_Wait(&request_receive_1,MPI_STATUS_IGNORE);
#if DEBUG_MODE_PT==Enabled
            std::cout << "step:"<< step <<" ,Rank: " << m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " Waiting completed from " <<DestRank <<std::endl;
#endif
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(NewTempID));
            }
            
        //BROADCAST TO UPDATE ALL TEMPERATURES ABOUT RESULT OF THE EXCHANGE ATTEMPT
        for(int i=0;i<m_Size;i++){
            if(i%2==0){
                MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
        }

        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" ,Rank: " << m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << ". The received non-sync broadcast is set." <<std::endl;
#endif


    }
    //IF COUNT IS EVEN, EVEN TEMPIDS RECEIVE ENERGY AND ATTEMPT EXCHANGE
    else if (TempIDIsEven){
        //DEALING WITH SPECIAL CASE. IF SIZE IS ODD, THE LAST TEMP ID IS OUT OF THE ATTEMPT.
        if (!m_SizeIsEven && m_TempID==m_Size-1){
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: " << m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << ". Setting special case broadcast for last temp ID." <<std::endl;
#endif
            for(int i=0;i<m_Size;i++){
                if(i%2==0){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << ". Set special case broadcast for last temp ID." <<std::endl;
#endif

            m_Counter++;

            return false;
        }


        //CALCULATING OWN ENERGY AND RECEIVING ENERGY FROM HIGHER CLOSER ODD RANK (TEMPID +1)
        double own_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        double receive_energy;
        int SourceTempID=m_TempID+1;
        int SourceRank=m_RankAtTempID[SourceTempID];
        MPI_Request request_receive;
        MPI_Request request_send_exchange_result;
        MPI_Request request_receive_last_1;
        MPI_Request request_receive_last_2;
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << ". 0 Receive from " <<SourceTempID<<std::endl;
#endif
        MPI_Irecv(&receive_energy, 1, MPI_DOUBLE, SourceRank, 0, MPI_COMM_WORLD, &request_receive);
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << ". 1 Receive from " <<SourceTempID<<std::endl;
#endif 
        MPI_Wait(&request_receive, MPI_STATUS_IGNORE);
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << ". 1 Wait from this completed, Receive from " <<SourceTempID<<std::endl;
#endif 
        //ATTEMPT EXCHANGE
        double dE=receive_energy-own_energy;
        double dBeta=m_BetaVec[SourceTempID] - m_BetaVec[m_TempID];
        double RandomNumber=m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        //std::cout<<RandomNumber<<std::endl;
        bool ExchangeAccepted=(RandomNumber<exp(dBeta*dE));

        //ExchangeAccepted=true;
        int ExchangeAcceptedInt = ExchangeAccepted ? 1 : 0;
        //SENDING RESULT OF ATTEMPT
        MPI_Isend(&ExchangeAcceptedInt, 1, MPI_INT, SourceRank, 1, MPI_COMM_WORLD, &request_send_exchange_result);
#if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << ". Send to " <<SourceTempID<<std::endl;
#endif         
        //MPI_Wait(&request_send_exchange_result, MPI_STATUS_IGNORE);
        //EFFECTIVE SWAP OF TEMPERATURES
        if (ExchangeAccepted){
            int NewTempID=SourceTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID));
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt_1;
            int EmptyBlockingInt_2=0;
            #if DEBUG_MODE_PT==Enabled
             std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " Receiving from "<<SourceRank<<std::endl;
            #endif     
            MPI_Irecv(&EmptyBlockingInt_1, 1, MPI_INT, SourceRank, 3, MPI_COMM_WORLD, &request_receive_last_1);
            #if DEBUG_MODE_PT==Enabled
            std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " Waiting from "<<SourceRank<<std::endl;
            #endif 
            MPI_Wait(&request_receive_last_1,MPI_STATUS_IGNORE);
            #if DEBUG_MODE_PT==Enabled
            std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << " Waiting completed from "<<SourceRank<<std::endl;
            #endif 
            MPI_Isend(&EmptyBlockingInt_2, 1, MPI_INT, SourceRank, 4, MPI_COMM_WORLD, &request_receive_last_2);
            #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Blocking reception from "<<SourceTempID<<std::endl;
#endif        
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(NewTempID));
        }
        

        //BROADCAST TO UPDATE ALL TEMPERATURES ABOUT RESULT OF THE EXCHANGE ATTEMPT

        
        for(int i=0;i<m_Size;i++){
            if(i%2==0 && i!=m_TempID){
                MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
            else if(i==m_TempID){
                MPI_Ibcast(&ExchangeAcceptedInt, 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
            }
        }
        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Setting and receiving broadcast "<<std::endl;
#endif       

        m_ReceiveBroadcast[m_TempID]=ExchangeAcceptedInt;

    }

    }
    //IF COUNT IS ODD, ODD TEMPIDS RECEIVE ENERGY FROM TEMPID +1 AND ATTEMPT EXCHANGE.
    else if(!CountIsEven){

    if (!TempIDIsEven){
        //SPECIAL CASE: IF SIZE IS EVEN AND COUNT IS ODD, THE LAST TEMP ID IS OUT OF THE ATTEMPT.
        if (m_SizeIsEven && m_TempID==m_Size-1){

            #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Setting and receiving broadcast special case "<<std::endl;
#endif       
            for(int i=0;i<m_Size;i++){
                if(i%1==0 && i!=m_TempID){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
            m_Counter++;
            //std::cout<<"Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", Counter: "<<m_Counter<<". End Of Evolve One Step."<<std::endl;
            return false;
        }

        //CALCULATING OWN ENERGY AND RECEIVING ENERGY FROM HIGHER CLOSER ODD RANK (TEMPID +1)
        double own_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        double receive_energy;
        int SourceTempID=m_TempID+1;
        int SourceRank=m_RankAtTempID[SourceTempID];
        MPI_Request requestA,requestB,requestC,requestD;
        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Receiving from "<<SourceTempID<<std::endl;
#endif      
        MPI_Irecv(&receive_energy, 1, MPI_DOUBLE, SourceRank, 0, MPI_COMM_WORLD, &requestA);
        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Received from  "<<SourceTempID<<std::endl;
#endif   
        MPI_Wait(&requestA,MPI_STATUS_IGNORE);
        //ATTEMPT EXCHANGE
        double dE=receive_energy-own_energy;
        double dBeta=m_BetaVec[SourceTempID] - m_BetaVec[m_TempID];
        double RandomNumber=m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        //std::cout<<RandomNumber<<std::endl;
        bool ExchangeAccepted=(RandomNumber<exp(dBeta*dE));


        //ExchangeAccepted=true;
        int ExchangeAcceptedInt = ExchangeAccepted ? 1 : 0;

        //SENDING RESULT OF ATTEMPT
        MPI_Isend(&ExchangeAcceptedInt, 1, MPI_INT, SourceRank, 1, MPI_COMM_WORLD,&requestB);
        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Sent exchange attempt result to "<<SourceTempID<<std::endl;
#endif   
        
        //EFFECTIVE SWAP OF TEMPERATURES
        if (ExchangeAccepted){
            int NewTempID=SourceTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID));
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt_1;
            int EmptyBlockingInt_2=0;

            MPI_Irecv(&EmptyBlockingInt_1, 1, MPI_INT, SourceRank, 3, MPI_COMM_WORLD, &requestC);
            MPI_Wait(&requestC, MPI_STATUS_IGNORE);
            MPI_Isend(&EmptyBlockingInt_2, 1, MPI_INT, SourceRank, 4, MPI_COMM_WORLD, &requestD);
            #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Received blocking from"<<SourceTempID<<std::endl;
#endif   
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(NewTempID));
        }
        

        //BROADCAST TO UPDATE ALL TEMPERATURES ABOUT RESULT OF THE EXCHANGE ATTEMPT
        
        for(int i=0;i<m_Size;i++){
            if (!m_SizeIsEven){
                if(i%1==0 && i!=m_TempID ){

                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
                else if(i==m_TempID){
                    MPI_Ibcast(&ExchangeAcceptedInt, 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
            else if(m_SizeIsEven){
                if(i%1==0 && i!=m_TempID && i!=m_Size-1){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
                else if(i==m_TempID){
                    MPI_Ibcast(&ExchangeAcceptedInt, 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
        }

        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Broadcast ready"<<SourceTempID<<std::endl;
#endif   

        m_ReceiveBroadcast[m_TempID]=ExchangeAcceptedInt;

        
    }
    //ODD COUNT, EVEN TEMPID SEND ENERGY TO TEMPID -1.
    else if (TempIDIsEven){
        if (m_TempID==0){
            #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Special case when size is even..."<<std::endl;
#endif   
            for(int i=0;i<m_Size;i++){
                if (m_SizeIsEven){
                    if(i%1==0 && i!=m_Size-1){
                        MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                    }
                }
                else if (!m_SizeIsEven){
                    if(i%1==0){
                        MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                    }
                }
            }
            //std::cout<<"Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", Counter: "<<m_Counter<<". End Of Evolve One Step"<<std::endl;
            m_Counter++;
            return false;
        }


        //CALCULATE AND SEND ENERGY
        double sending_energy=m_pState->GetEnergyCalculator()->GetEnergy();
        int DestTempID=m_TempID-1;
        int DestRank=m_RankAtTempID[DestTempID];
        int ReceiveExchangeAcceptedInt;
        MPI_Request requestA,requestB,requestC,requestD;
        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Send to "<<DestTempID<<std::endl;
#endif 
        MPI_Isend(&sending_energy, 1, MPI_DOUBLE, DestRank, 0, MPI_COMM_WORLD,&requestA);
        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Sent to "<<DestTempID<<std::endl;
#endif 
        //RECEIVE RESULT OF EXCHANGE ATTEMPTS
        MPI_Irecv(&ReceiveExchangeAcceptedInt, 1, MPI_INT, DestRank, 1, MPI_COMM_WORLD, &requestB);

        MPI_Wait(&requestB,MPI_STATUS_IGNORE);
        bool ReceiveExchangeAccepted=(ReceiveExchangeAcceptedInt==1);

        #if DEBUG_MODE_PT==Enabled
        std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Received from "<<DestTempID<<std::endl;
#endif 

        //IF ACCEPTED, EFFECTIVELY SWAP RANKS
        if (ReceiveExchangeAccepted){
            int NewTempID=DestTempID;
            m_pState->GetSimulation()->SetBeta(m_BetaVec[NewTempID], 0);
            m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(NewTempID)); 
            m_pState->GetTimeSeriesDataOutput()->CloseFile();
            int EmptyBlockingInt_1=0;
            int EmptyBlockingInt_2;
       
            MPI_Isend(&EmptyBlockingInt_1, 1, MPI_INT, DestRank, 3, MPI_COMM_WORLD,&requestC);
            MPI_Irecv(&EmptyBlockingInt_2, 1, MPI_INT, DestRank, 4, MPI_COMM_WORLD,&requestD);
            MPI_Wait(&requestD,MPI_STATUS_IGNORE);
            #if DEBUG_MODE_PT==Enabled
            std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Blocking sent to "<<DestTempID<<std::endl;
#endif 
            m_pState->GetTimeSeriesDataOutput()->OpenFileWithoutHeader(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(NewTempID)+TimeSeriDataExt);
            m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(NewTempID));
            }
            
        //BROADCAST TO UPDATE ALL TEMPERATURES ABOUT RESULT OF THE EXCHANGE ATTEMPT


        for(int i=0;i<m_Size;i++){
            if (!m_SizeIsEven){
                if(i%1==0){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
            else if(m_SizeIsEven){
                if(i%1==0 && i!=m_Size-1){
                    MPI_Ibcast(&m_ReceiveBroadcast[i], 1, MPI_INT, m_RankAtTempID[i], MPI_COMM_WORLD, &m_RequestBroadcast[i]);
                }
            }
        }
        #if DEBUG_MODE_PT==Enabled
            std::cout <<"step:"<< step <<" , Rank: "<< m_Rank << ", TempID: " << m_TempID 
        << ", Counter: " << m_Counter << "Broadcasting"<<std::endl;
#endif 
        
    }


    }
    //std::cout<<"Rank: "<<m_Rank<<", TempID: "<<m_TempID<<", Counter: "<<m_Counter<<". End Of Evolve One Step"<<std::endl;
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

std::vector<double> ParallelTemperingMoveSimple::ReadTemperatures(){

    Nfunction f;
    std::string filename=ParallelTemperingMoveSimple::GetInitialTemperaturesFileName();

    if(f.FileExist(filename)==false)
    {
        std::cout<<"-----> Error: the temperature file name with the name "<<filename<<" does not exist"<<std::endl;
        exit(0);
    }

    std::ifstream TemperatureFile(filename);
    
    // Check if the file was opened successfully
    if (!TemperatureFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(0);
    }

    std::vector<double> betavec;

    // File processing code goes here

    std::string line;
    while (std::getline(TemperatureFile, line)) {
        // Skip empty or whitespace-only lines
        if (is_line_empty(line)) {
            continue;
        }

        std::stringstream ss(line);  // Use a stringstream to parse the line
        double value;
        // Read the double value from the line
        if (ss >> value) {
            // Ensure there is no extra data in the line
            if (ss >> std::ws && !ss.eof()) {
                std::cerr << "Error: More than one value or extra data in line: \"" << line << "\"" << std::endl;
                exit(0);  // Exit with error code
            }

            
            betavec.push_back(value);  // Add each double to the data vector
        } else {
            std::cerr << "Error: Could not read a double value from line: \"" << line << "\"" << std::endl;
            exit(0);  // Exit with error code
        }
    }

    // Close the file
    TemperatureFile.close();

    for (size_t i = 1; i < betavec.size(); ++i) {
        if (betavec[i] < betavec[i - 1]) {
            std::cerr << "Error: Vector is not sorted in ascending order." << std::endl;
            exit(0);  // Exit with error code
        }
    }

    for (size_t i = 0; i < betavec.size(); ++i) {
        betavec[i]=1/betavec[i];
    }

    return betavec;

}

bool ParallelTemperingMoveSimple::is_line_empty(const std::string& line) {
    // Check if the line is empty or contains only whitespace
    return line.empty() || std::all_of(line.begin(), line.end(), [](unsigned char c) { return std::isspace(c); });
}

bool ParallelTemperingMoveSimple::GetTargetState(){
    return m_TargetState;
}

void ParallelTemperingMoveSimple::SetRestart(){
    m_Restart=true;
}

// Function to check if a string is a valid integer
bool ParallelTemperingMoveSimple::isInteger(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

// Function to split a string by spaces and convert to integers
std::vector<int> ParallelTemperingMoveSimple::splitStringToIntVector(const std::string& line) {
    std::vector<int> result;
    std::stringstream ss(line);
    std::string temp;

    while (ss >> temp) {
        if (isInteger(temp)) {
            result.push_back(std::stoi(temp));
        } else {
            // If any part is not an integer, return an empty vector
            return {};
        }
    }

    return result;
}

// Function to read the last line of a file and parse it into a vector
std::vector<int> ParallelTemperingMoveSimple::ReadLastLineOutput(const std::string& output, int N) {
    std::ifstream file(output);
    std::string lastLine;
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << output << std::endl;
        // Return a default vector if the file cannot be opened
        std::vector<int> defaultVec;
        for (int i = 0; i < N; ++i) {
            defaultVec.push_back(i);
        }
        return defaultVec;
    }
    std::cout<<"Here"<<std::endl;
    // Read the file to find the last non-empty line
    while (std::getline(file, line)) {
        if (!line.empty()) {
            lastLine = line;
        }
    }
    std::cout<<"Here 2"<<std::endl;

    file.close();

    // Try to convert the last line to a vector of integers
    std::vector<int> result = splitStringToIntVector(lastLine);
    if (!result.empty() && result.size() == N) {
        return result;
    }

    // If the last line is not a valid integer vector, return a default vector
    std::vector<int> defaultVec;
    for (int i = 0; i < N; ++i) {
        defaultVec.push_back(i);
    }
    return defaultVec;
}


std::string ParallelTemperingMoveSimple::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    //What is this used for?
    //state = state +" "+ Nfunction::D2S(m_Period)+" "+Nfunction::D2S(m_SigmaP)+" "+  m_Type;
    return state;
}


#endif