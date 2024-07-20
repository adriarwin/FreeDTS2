

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include "TimeSeriesDataOutput.h"
#include "State.h"
#include "Nfunction.h"
#include "SimDef.h"

TimeSeriesDataOutput::TimeSeriesDataOutput(){
    m_Periodic = 0;
}
TimeSeriesDataOutput::TimeSeriesDataOutput(State *pState){
    m_pState = pState;
    m_Periodic = 100;
}
TimeSeriesDataOutput::~TimeSeriesDataOutput(){
    if (m_TimeSeriesFile.is_open()) {
        m_TimeSeriesFile.close();
     }
}
void TimeSeriesDataOutput::UpdatePeriod(int period){
    
    m_Periodic = period;
    return;
}
bool TimeSeriesDataOutput::WrireTimeSeriesDataOutput(int step){
    if( m_Periodic == 0 || step%m_Periodic!=0)
        return false;
    m_TimeSeriesFile<<std::fixed;
    m_TimeSeriesFile<<std::setprecision(Precision);
//--> write step and energy
    m_TimeSeriesFile<<step<<"   "<<m_pState->GetSystemEnergy()<<"   ";
//--->write nematic force energy
    if((m_pState->m_pConstant_NematicForce)->m_F0!=0){
        m_TimeSeriesFile<<(m_pState->m_pConstant_NematicForce)->m_ActiveEnergy<<"  ";
    }
//--->write box size
    if(m_pState->GetDynamicBox()->GetTau()!=0){
        m_TimeSeriesFile<<(*(m_pState->m_pMesh->GetBox()))(0)<<"  ";
        m_TimeSeriesFile<<(*(m_pState->m_pMesh->GetBox()))(1)<<"  ";
        m_TimeSeriesFile<<(*(m_pState->m_pMesh->GetBox()))(2)<<"  ";
    }
//---> global curvature
    if (m_pState->GetGlobalCurvature()->GetState() == true){
        m_TimeSeriesFile<<m_pState->GetGlobalCurvature()->GetEnergy()<<"  ";
    }
//----> harmonic force
    if(m_pState->Get2GroupHarmonic()->GetState()==true){
        m_TimeSeriesFile<<m_pState->Get2GroupHarmonic()->GetEnergy()<<"  ";
        m_TimeSeriesFile<<m_pState->Get2GroupHarmonic()->GetForce()<<"  ";
        m_TimeSeriesFile<<m_pState->Get2GroupHarmonic()->GetDistance()<<"  ";
    }
//--> volume coupling
    if(m_pState->GetVolumeCoupling()->GetState()==true){
        m_TimeSeriesFile<<m_pState->GetVolumeCoupling()->GetTotalVolume()<<"   ";
        m_TimeSeriesFile<<m_pState->GetVolumeCoupling()->GetTotalArea()<<"  ";
    }
//--> constant area
    if((m_pState->GetApply_Constant_Area())->GetState()==true){
        m_TimeSeriesFile<<"   "<<(m_pState->GetApply_Constant_Area())->GetTotalArea()<<"  ";
    }
//--> active exchange
    if(m_pState->GetActiveTwoStateInclusion()->GetState() == true)
    {
        m_TimeSeriesFile<<*(m_pState->GetActiveTwoStateInclusion()->GetActiveEnergy())<<"   ";
        m_TimeSeriesFile<<m_pState->GetActiveTwoStateInclusion()->GetDeltaN()<<"   ";
    }
    m_TimeSeriesFile<<std::endl;
    
    return true;
}
bool TimeSeriesDataOutput::OpenTimeSeriesDataOutput(){
    
    std::string filename = m_pState->m_GeneralOutputFilename;
    filename = filename + "-en.xvg";
    
    //-- if it is a restart simulation, then do not clear the file and check for the last line in the file
    if((m_pState->m_RESTART).restartState==true){
        if(!CheckTimeSeriesFile(m_pState->m_Initial_Step)){
            std::cerr << "---> error: energy file does not match with the restart: " << filename << std::endl;
        }
        m_TimeSeriesFile.open(filename.c_str(),std::fstream::app);
    }
    else{
        m_TimeSeriesFile.open(filename.c_str());
    }
    
    if (!m_TimeSeriesFile.is_open()) {
        std::cerr << "---> error: Unable to open file: " << filename << std::endl;
        return false;
    }
//-- write the header if not a restart
    if((m_pState->m_RESTART).restartState==false){
        
            m_TimeSeriesFile<<" ## mcstep  energy ";
            if((m_pState->m_pConstant_NematicForce)->m_F0!=0)
                    m_TimeSeriesFile<<"  active_nematic_energy ";
            if(m_pState->GetDynamicBox()->GetTau()!=0)
                    m_TimeSeriesFile<<" Lx  Ly  Lz ";
            if (m_pState->GetGlobalCurvature()->GetState() == true)
                    m_TimeSeriesFile<<" global_curvature_energy  ";
            if(m_pState->Get2GroupHarmonic()->GetState()==true)
                    m_TimeSeriesFile<<" harmonic_energy harmonic_force distance ";
            if(m_pState->GetVolumeCoupling()->GetState()==true)
                    m_TimeSeriesFile<<" volume area ";
            if((m_pState->GetApply_Constant_Area())->GetState()==true)
                    m_TimeSeriesFile<<" area_constantArea ";
            if(m_pState->GetActiveTwoStateInclusion()->GetState() == true)
                    m_TimeSeriesFile<<" ActiveEnergy DeltaN ";
        
        m_TimeSeriesFile<<std::endl;
    }
        
    return true;
}
bool TimeSeriesDataOutput::FlushTimeSeriesFile(){ // the energy file should be flushed first

    m_TimeSeriesFile.flush(); //

    return true;
}
bool TimeSeriesDataOutput::CheckTimeSeriesFile(int endstep){ // this checks if the energy file matches the restart end step

   // This task is already is being done in the restart class
    
    /*std::string filename = m_pState->m_GeneralOutputFilename;
    filename = filename + "-en.xvg";
    energyfile.open(filename.c_str());

    while(true)
    {
        ren>>step;
        if(ren.eof())
        {
            std::cout<<"----> warning: the energy file does not contains enough output, some is missing "<<std::endl;
            break;
        }
        if(step>=m_pState->m_Initial_Step-1)
            break;
        getline(ren,enstr);
        energyfile<<step<<enstr<<"\n";

    }
    energyfile.close();
    ren.close();
    remove("en.txt");
    */
    
    return true;
}
