

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This class distrubutes the tasks based on inputs provided and makes all the initials variables
 It also receives a blue print of the mesh and the generates the reall mesh, then pass it on to Simualtion object to start ..
 */
#include <fstream>
#include "State.h"
#include <time.h>
#include "MESH.h"
#include "CreateMashBluePrint.h"

State::State()
{
}
State::State(std::vector <std::string> argument)
{
    Nfunction f;
    f.CleanFiles();
    Restart re(this);
    m_Restart = re;
#if TEST_MODE == Enabled
    std::cout<<"----> We have reached the State Class -- "<<std::endl;
#endif
    //============ Initialization of all inputs and data structures for input

    m_pConstant_NematicForce = &m_Constant_NematicForce;
    m_Argument = argument;
    m_Targeted_State = true;
    m_Healthy =true;
    m_Integrator = "MC";
    m_Total_no_Threads = 1;
    m_Initial_Step = 0;
    m_Final_Step = 0;
    m_Seed =36723;
    m_R_Vertex=0.05;   // Move Vertex  within a box with this size
    m_R_Box=0.04;   // box change within this range
    m_TotEnergy = 0;
    m_Mem_Spontaneous_Curvature = 0.0;
    m_TopologyFile = "TrjTSI";
    m_InputFileName = "input.dts";
    m_IndexFileName = "Index.inx";
    m_IndexFile = false;
    m_GeneralOutputFilename = "output";
    m_FreezGroupName = "";
    m_Beta = 1;
    m_MinFaceAngle = -0.5;
    m_MinVerticesDistanceSquare = 1.0;
    m_MaxLinkLengthSquare = 3.0;
    m_Display_periodic  = 1000;
    m_Centering         = 0;
    m_OutPutEnergy_periodic = 100;
    m_CNTCELL(0) = 2;m_CNTCELL(1) = 2;m_CNTCELL(2) = 2;
    m_Parallel_Tempering.State = false;
    m_TRJTSI.tsiPeriod = 0;
    m_TRJTSI.tsiPrecision = 1;
    m_TRJTSI.tsiFolder_name = "tsi";
    m_TRJTSI.tsiState = false;
    m_TRJBTS.btsPeriod = 0;
    m_TRJBTS.btsPrecision = 1;
    m_TRJBTS.btsFile_name = "dts.bts";
    m_TRJBTS.btsState = false;
    m_RESTART.restartState = false;
    m_RESTART.restartPeriod = 1000;
    m_MCMove.VertexMove = true;
    m_MCMove.EdgeVertexMove = true;
    m_MCMove.LinkFlip = true;
    m_MCMove.InclusionMove = true;
    m_MCMove.VertexMoveRate = 1;
    m_MCMove.LinkFlipRate = 1;
    m_MCMove.InclusionMoveRate = 1;
    m_MCMove.EdgeVertexMoveRate = 1;
    m_FrameTension.State = false;
    m_FrameTension.Type = "Position_Rescale";
    m_FrameTension.Tau = 0.0;
    m_FrameTension.updatePeriod = 5;
    m_VolumeConstraint.State = false;
    m_VolumeConstraint.EQSteps = 1000;
    m_VolumeConstraint.DeltaP = 0;
    m_VolumeConstraint.K = 0;
    m_VolumeConstraint.targetV = 1 ;
    m_STRUC_OSMOTIC.State = false;
    m_STRUC_ConstantArea.State = false;
    m_STRUC_ConstantVertexArea.State = false;
    m_STRUC_ActiveTwoStateInclusion.state = false;
    m_STRUC_ActiveTwoStateInclusion.nametype1 = " ";
    m_STRUC_ActiveTwoStateInclusion.nametype2 = " ";
    m_STRUC_ActiveTwoStateInclusion.ep1 = 0;
    m_STRUC_ActiveTwoStateInclusion.ep2 = 0;
    m_STRUC_ActiveTwoStateInclusion.persentage = 0;
    m_STRUC_ActiveTwoStateInclusion.gama = 0;
    m_pMesh = &m_Mesh;

    //====== Updating the input data from provided inputs
    ExploreArguments();     // Find input file name
    
    // this should be changed at some point. Now it needed to be exactly in this place because need the name of the input file but it also needs updates the bending rigidity and etc. so the copy becomes a problem if it is not here.
    Inclusion_Interaction_Map inc_ForceField(m_InputFileName);
    m_inc_ForceField = inc_ForceField;
    m_pinc_ForceField = &m_inc_ForceField;
    
    
    ReadInputFile(m_InputFileName);  // Reading from input file
    ExploreArguments();     // Update last state
    




//
}
State::~State()
{
    
}
void State::ExploreArguments()
{
    Nfunction f;
    for (long i=1;i<m_Argument.size();i=i+2)
    {
        std::string Arg1 = m_Argument.at(i);
        if(Arg1=="-h")
        {
            HelpMessage();
            exit(0);
            break;
        }
        else if(Arg1=="-in")
        {
            m_InputFileName = m_Argument.at(i+1);
        }
        else if(Arg1=="-top")
        {
            m_TopologyFile = m_Argument.at(i+1);
        }
        else if(Arg1=="-defout")
        {
            m_GeneralOutputFilename = m_Argument.at(i+1);
        }
        else if(Arg1=="-ndx")
        {
            m_IndexFileName = m_Argument.at(i+1);
            m_IndexFile = true;
        }
        else if(Arg1=="-b")
        {
            m_Initial_Step = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-nt")
        {
            m_Total_no_Threads = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-e")
        {
            m_Final_Step = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-seed")
        {
            m_Seed = f.String_to_Int(m_Argument.at(i+1));
        }
        else if(Arg1=="-angle")
        {
            m_MinFaceAngle = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-maxDist")
        {
            m_MaxLinkLengthSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-minDist")
        {
            m_MinVerticesDistanceSquare = f.String_to_Double(m_Argument.at(i+1));
        }
        else if(Arg1=="-restart")
        {
            m_RESTART.restartState = true;
            m_RESTART.restartFilename = m_Argument.at(i+1);
        }
        else
        {
            std::cout << "Error: simulation will stop!!:"<<Arg1;
            std::cout<<"\n"<<"For more information and tips run "<< m_Argument.at(0) <<" -h"<<"\n";
            m_Healthy =false;
            exit(0);
            break;
        }
    }
    // check if the files declared in the input file exist
    if (f.FileExist (m_TopologyFile)!=true && m_RESTART.restartState !=true)
    {
        std::cout<<" Error: the topology file  does not exist "<<std::endl;
        m_Healthy =false;
        exit(0);
    }
}
void State::ReadInputFile(std::string file)
{
    Nfunction f;
    // enforcing correct file extension:  check simdef file for value of InExt
    std::string ext = file.substr(file.find_last_of(".") + 1);
    if(ext!=InExt)
    file = file + "." + InExt;
    if (f.FileExist(file)!=true)
    {
        std::cout<<"----> Error: the input file with the name "<<file<< " does not exist "<<std::endl;
        m_Healthy =false;
        exit(0);
    }

    std::ifstream input;
    input.open(file.c_str());
    std::string firstword,rest,str;

    while (true)
    {
        input>>firstword;
        if(input.eof())
            break;
        
        
        if(firstword == "Initial_Step")
        {
            input>>str>>m_Initial_Step;
            getline(input,rest);
        }
        else if(firstword == "Integrator")
        {
            input>>str>>m_Integrator;
            getline(input,rest);
        }
        else if(firstword == "MC_Moves")
        {
            input>>str>>(m_MCMove.VertexMove)>>(m_MCMove.LinkFlip)>>(m_MCMove.EdgeVertexMove)>>(m_MCMove.InclusionMove);
            getline(input,rest);
        }
        else if(firstword == "MC_MovesRate")
        {
            input>>str>>(m_MCMove.VertexMoveRate)>>(m_MCMove.LinkFlipRate)>>(m_MCMove.EdgeVertexMoveRate)>>(m_MCMove.InclusionMoveRate);
            getline(input,rest);
            
            if(m_MCMove.VertexMoveRate<0)
            {
                std::cout<<"---> Error: move rate is negative "<<std::endl;
                exit(0);
            }
            if(m_MCMove.LinkFlipRate<0)
            {
                std::cout<<"---> Error: move rate is negative "<<std::endl;
                exit(0);
            }
            if(m_MCMove.InclusionMoveRate<0)
            {
                std::cout<<"---> Error: move rate is negative "<<std::endl;
                exit(0);
            }
            if(m_MCMove.EdgeVertexMoveRate<0)
            {
                std::cout<<"---> Error: move rate is negative "<<std::endl;
                exit(0);
            }
        }
        else if(firstword == "FreezingAGroup")
        {
            input>>str>>m_FreezGroupName;
            getline(input,rest);
        }
        else if(firstword == "Constant_NematicForce")
        {
            double f0,fd,fn; //(d,p,n)
            input>>str>>f0>>fd>>fn;
            m_Constant_NematicForce.m_F0 = f0;
            m_Constant_NematicForce.m_Fd = fd;
            m_Constant_NematicForce.m_Fn = fn;
            getline(input,rest);
        }
        else if(firstword=="ActiveTwoStateInclusion")
        {
            // ActiveTwoStateInclusion = on p1 p2 ep1 ep2 0.1 gama
            double ep1,ep2,persentage,gama;
            std::string state;
            input>>str>>state;
            if(state == "On" || state=="ON" || state =="on" || state =="yes" || state =="Yes")
            m_STRUC_ActiveTwoStateInclusion.state = true;
            input>>m_STRUC_ActiveTwoStateInclusion.nametype1;
            input>> m_STRUC_ActiveTwoStateInclusion.nametype2;
            input>>ep1>>ep2>>persentage>>gama;
            m_STRUC_ActiveTwoStateInclusion.ep1 = ep1;
            m_STRUC_ActiveTwoStateInclusion.ep2 = ep2;
            m_STRUC_ActiveTwoStateInclusion.persentage = persentage;
            m_STRUC_ActiveTwoStateInclusion.gama = gama;
            getline(input,rest);
        }
        else if(firstword == "Frame_Tension")
        {
            std::string state;
            (m_FrameTension.State)= false;
            input>>str>>state>>(m_FrameTension.Type)>>(m_FrameTension.Tau)>>(m_FrameTension.updatePeriod);
            if(state == "on" || state == "ON" || state == "On" )
                (m_FrameTension.State) = true;
            getline(input,rest);
            
            
        }
        else if(firstword == "Volume_Constraint")
        {
            std::string state;
            m_VolumeConstraint.State = false;
        input>>str>>state>>(m_VolumeConstraint.EQSteps)>>(m_VolumeConstraint.DeltaP)>>(m_VolumeConstraint.K)>>(m_VolumeConstraint.targetV);
            if(state == "on" || state == "ON" || state == "On" )
               m_VolumeConstraint.State = true;
            getline(input,rest);
            
            
            CmdVolumeCouplingSecondOrder C(m_VolumeConstraint.State, m_VolumeConstraint.EQSteps, m_VolumeConstraint.DeltaP,  m_VolumeConstraint.K, m_VolumeConstraint.targetV);
            m_VolumeCouplingSecondOrder = C;
        }
        else if(firstword == "Osmotic_Pressure")
        {
                    std::string state;
          input>>str>>state>>(m_STRUC_OSMOTIC.Type)>>(m_STRUC_OSMOTIC.EQSteps)>>(m_STRUC_OSMOTIC.Gamma)>>(m_STRUC_OSMOTIC.P0);
          if(state=="on" || state=="ON" || state=="On" || state=="yes" || state=="Yes" )
          m_STRUC_OSMOTIC.State = true;
          Apply_Osmotic_Pressure C(m_STRUC_OSMOTIC.State,m_STRUC_OSMOTIC.Type,m_STRUC_OSMOTIC.EQSteps,m_STRUC_OSMOTIC.Gamma,m_STRUC_OSMOTIC.P0);
          m_Apply_Osmotic_Pressure = C;
          
        }
        else if(firstword == "Apply_Constant_Area")
        {
            std::string state;
            m_STRUC_ConstantArea.State = false;
            input>>str>>state>>(m_STRUC_ConstantArea.EQSteps)>>(m_STRUC_ConstantArea.Gamma)>>(m_STRUC_ConstantArea.K0);
            if(state=="on" || state=="ON" || state=="On" || state=="yes" || state=="Yes" )
                m_STRUC_ConstantArea.State = true;
            Apply_Constant_Area C(m_STRUC_ConstantArea.State,m_STRUC_ConstantArea.EQSteps,m_STRUC_ConstantArea.Gamma,m_STRUC_ConstantArea.K0);
            m_Apply_Constant_Area = C;
            
        }
        else if(firstword == "Final_Step")
        {
            input>>str>>m_Final_Step;
            getline(input,rest);
        }
        else if(firstword == "Box_Centering_F")
        {
            input>>str>>m_Centering;
            getline(input,rest);
        }
        else if(firstword == "MinfaceAngle")
        {
            input>>str>>m_MinFaceAngle;
            getline(input,rest);
        }
        else if(firstword == "OutPutEnergy_periodic")
        {
            input>>str>>m_OutPutEnergy_periodic;
            getline(input,rest);
        }
        else if(firstword == "Restart_periodic")
        {
            input>>str>>(m_RESTART.restartPeriod);
            getline(input,rest);
        }
        else if(firstword == "TopologyFile")
        {
            input>>str>>m_TopologyFile;
            getline(input,rest);
        }
        else if(firstword == "CouplingtoFixedGlobalCurvature")
        {
            std::string state;
            double k,c0;
            input>>str>>state>>k>>c0;
            if(state == "yes" || state == "on")
            {
                CouplingtoFixedGlobalCurvature  CoupleGCurvature(true,k,c0);
                m_CoupleGCurvature = CoupleGCurvature;
            }
            else
            {
                CouplingtoFixedGlobalCurvature  CoupleGCurvature(false,k,c0);
                m_CoupleGCurvature = CoupleGCurvature;
            }
            getline(input,rest);
        }
        else if(firstword == "HarmonicPotentialBetweenTwoGroups")
        {
            //HarmonicPotentialBetweenTwoGroups = on 10 0.1 2000 Group1 Group2 0 1 1
            std::string state,g1,g2;
            double k,dx;
            int nx,ny,nz,rate;
            input>>str>>state>>k>>dx>>rate>>g1>>g2>>nx>>ny>>nz;
            if(state == "yes" || state == "on")
            {
                SpringPotentialBetweenTwoGroups  tem(true,k, dx,rate,g1,g2,nx,ny,nz);
                m_SpringPotentialBetweenTwoGroups = tem;
            }
            else
            {
                SpringPotentialBetweenTwoGroups  tem(false,k, dx,rate,g1,g2,nx,ny,nz);
                m_SpringPotentialBetweenTwoGroups = tem;
            }
        }
        else if(firstword == "Seed")
        {
            input>>str>>m_Seed;
            getline(input,rest);
        }
        else if(firstword == "Kappa")
        {
            double a,b;
            input>>str>>a>>b;
            m_inc_ForceField.m_BendingRigidity = a/2;
            m_inc_ForceField.m_GaussianRigidity = b;
            getline(input,rest);
        }
        else if(firstword == "Spont_C")
        {
            double a;
            input>>str>>m_Mem_Spontaneous_Curvature;
            m_inc_ForceField.m_Spontaneous_Curvature = a;
            getline(input,rest);
        }
        else if(firstword == "Edge_Parameters")
        {
            double a,b,c;
            input>>str>>a>>b>>c;
            
            //=== send it to force field class
            m_inc_ForceField.m_Lambda = a;
            m_inc_ForceField.m_KgEdge = b;
            m_inc_ForceField.m_KnEdge = c;

            getline(input,rest);
        }
        else if(firstword == "VertexArea")
        {
            double a,b,c,d;
            input>>str>>a>>b>>c>>d;
            
            double A01 = (1+2*b)*sqrt(3)/2.0;   // selecting b between 0-1
            double A02 = sqrt(1+2*d);   // selecting d between 0-1

            if(b<0 || b>1 || d<0 || d>1)
            {
            std::cout<<"---> error in constant vertex area; gamma is bad; it should be a double number between 0-1 \n";
            exit(0);
            }
            m_inc_ForceField.m_Kva = a/2;
            m_inc_ForceField.m_KvaEdge = c/2;

            m_inc_ForceField.m_av0 = A01;
            m_inc_ForceField.m_av0Edge = A02;

        }
        else if(firstword == "Mem_model_para") ///
        {
            input>>str;
            getline(input,rest);
            std::vector<std::string> memdata = f.split(rest);
            std::vector <double> modelparam;
            for (std::vector<std::string>::iterator it = memdata.begin() ; it != memdata.end(); ++it)
            {
                if((*it)!=";")
                {
                    double value = f.String_to_Double(*it);
                    modelparam.push_back(value);
                }
            }

            m_inc_ForceField.m_Membrane_model_parameters = modelparam;
        }
        else if(firstword == "Display_periodic")
        {
            input>>str>>m_Display_periodic;
            getline(input,rest);
        }
        else if(firstword == "Cell_Size")
        {
            input>>str>>m_CNTCELL(0)>>m_CNTCELL(1)>>m_CNTCELL(2);
            getline(input,rest);
        }
        else if(firstword == "GeneralOutputFilename")
        {
            input>>str>>m_GeneralOutputFilename;
            getline(input,rest);
        }
        else if(firstword == "Min_Max_LinkLenghtsSquare")
        {
            input>>str>>m_MinVerticesDistanceSquare>>m_MaxLinkLengthSquare;
            getline(input,rest);
        }
        else if(firstword == "CoupleToRigidWalls")
        {
            input>>str>>str;
            bool state = false;
            if(str == "on" || str == "yes" || str == "Yes")
            state = true;
            getline(input,rest);
            
            CoupleToWallPotential wall(state,rest);
            m_RigidWallCoupling = wall;
        }
        else if(firstword == "OutPutTRJ_TSI")
        {
            input>>str>>(m_TRJTSI.tsiPeriod)>>(m_TRJTSI.tsiPrecision)>>(m_TRJTSI.tsiFolder_name);
            getline(input,rest);
            if((m_TRJTSI.tsiPeriod)>0)
                m_TRJTSI.tsiState = true;
            else
                m_TRJTSI.tsiState = false;
        }
        else if(firstword == "Parallel_Tempering")
        {
            // Parallel_Tempering  = on PT_steps  PT_minbeta    PT_maxbeta
            std::string state;
            input>>str>>state>>(m_Parallel_Tempering.PT_steps)>>(m_Parallel_Tempering.PT_minbeta)>>(m_Parallel_Tempering.PT_maxbeta);
            getline(input,rest);
            if(state=="on"|| state=="yes"|| state=="On"|| state=="ON"|| state=="Yes"|| state=="YES")
                m_Parallel_Tempering.State = true;
            else
                m_Parallel_Tempering.State = false;
        }
        else if(firstword == "OpenEdgeEvolutionWithConstantVertex")
        {
            // OpenEdgeEvolutionWithConstantVertex  = on PT_steps
            std::string state;
            int rate;
            input>>str>>state>>rate;
            getline(input,rest);
            bool edgestate = false;
            if(state=="on"|| state=="yes"|| state=="On"|| state=="ON"|| state=="Yes"|| state=="YES")
                edgestate = true;
            
            
        }
        else if(firstword == "OutPutTRJ_BTS")
        {
            input>>str>>(m_TRJBTS.btsPeriod)>>(m_TRJBTS.btsPrecision)>>(m_TRJBTS.btsFile_name);
            getline(input,rest);
            if((m_TRJBTS.btsPeriod)>0)
                m_TRJBTS.btsState = true;
            else
                m_TRJBTS.btsState = false;
        }
        else if(firstword == "INCLUSION")
        {
            break;
        }
        else
        {
            if(firstword.at(0)!=';')
            {
                std::cout<<"Error: bad keyword in the input file *** "<<firstword<<" ***\n";
                exit(0);
            }
            getline(input,rest);
        }
    }
    input.close();
    
}
void State::HelpMessage()
{
    std::cout<<"----------------------copyright: Weria Pezeshkian------------------"<<"\n";
    std::cout<<"---------------------- weria.pezeshkian@gmail.com ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"---------------------version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"------------simple example for exacuting MC simulation -------------------"<<"\n";
    std::cout<< "   DTS -f Input.dts  -top topol.q -restart res.res"<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  option    type        default            description "<<"\n";
    std::cout<<"-------------------------------------------------------------------------------"<<"\n";
    std::cout<<"  -in        string       Input.dts         input file name "<<"\n";
    std::cout<<"  -top       string       topology.top      Topology file name "<<"\n";
    std::cout<<"  -b         int          1                 initial time step "<<"\n";
    std::cout<<"  -e         int          100               final time step  "<<"\n";
    std::cout<<"  -seed      int          36723             Random number seed  "<<"\n";
    std::cout<<"  -defout    string       output            an string for run out put files   "<<"\n";
    std::cout<<"  -ndx       string       Index.inx         Index file name  "<<"\n";
    std::cout<<"  -nt        int          1                 Total number of Threads  "<<"\n";
    std::cout<<"  -restart   string       NO                restart file name  "<<"\n";
    std::cout<<"  -angle     double       -0.5              minimum of cos of the angle between two faces  "<<"\n";
    std::cout<<"  -minDist   double       1                 Square of minimum distance between two vertices "<<"\n";
    std::cout<<"  -maxDist   double       3                 Square of maximum link length "<<"\n";
    std::cout<<"=========================================================================="<<"\n";
    std::cout<<"------------------ version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    
}
void State::WriteStateLog()
{
    std::ofstream statelog;
    if (m_RESTART.restartState == true)
    statelog.open("statelog.log",std::fstream::app);
    else
    statelog.open("statelog.log");

    statelog<<"; this file was generated by below command "<<std::endl;
    statelog<<";  ";
    for (std::vector<std::string>::iterator it = m_Argument.begin() ; it != m_Argument.end(); ++it)
    {
        statelog<<(*it)<<"   ";
    }
    statelog<<std::endl;
    statelog<<";------------------------------------------  "<<std::endl;
    statelog<<"Integrator = "<<m_Integrator<<std::endl;
    statelog<<"MC_Moves = "<<(m_MCMove.VertexMove)<<"  "<<(m_MCMove.LinkFlip)<<"  "<<(m_MCMove.InclusionMove)<<std::endl;
    statelog<<"MC_MovesRate = "<<(m_MCMove.VertexMoveRate)<<"  "<<(m_MCMove.LinkFlipRate)<<"  "<<(m_MCMove.InclusionMoveRate)<<std::endl;
    statelog<<"Initial_Step = "<<m_Initial_Step<<std::endl;
    statelog<<"Final_Step = "<<m_Final_Step<<std::endl;
    statelog<<"MinfaceAngle = "<<m_MinFaceAngle<<std::endl;
    statelog<<"OutPutEnergy_periodic = "<<m_OutPutEnergy_periodic<<std::endl;
    statelog<<"Restart_periodic = "<<m_RESTART.restartPeriod<<std::endl;
    statelog<<"TopologyFile = "<<m_TopologyFile<<std::endl;
    statelog<<"Seed = "<<m_Seed<<std::endl;
    statelog<<"Kappa = "<<m_inc_ForceField.m_BendingRigidity<<"  "<<m_inc_ForceField.m_GaussianRigidity<<std::endl;
    statelog<<"Display_periodic = "<<m_Display_periodic<<std::endl;
    statelog<<"CNTCELL = "<<m_CNTCELL(0)<<" "<<m_CNTCELL(1)<<" "<<m_CNTCELL(2)<<" "<<std::endl;
    statelog<<"GeneralOutputFilename = "<<m_GeneralOutputFilename<<std::endl;
    statelog<<"Min_Max_LinkLenghtsSquare = "<<m_MinVerticesDistanceSquare<<"  "<<m_MaxLinkLengthSquare<<std::endl;
    statelog<<"OutPutTRJ_TSI = "<<(m_TRJTSI.tsiPeriod)<<"  "<<(m_TRJTSI.tsiPrecision)<<"  "<<(m_TRJTSI.tsiFolder_name)<<"  "<<std::endl;
    statelog<<"OutPutTRJ_BTS = "<<(m_TRJBTS.btsPeriod)<<"  "<<(m_TRJBTS.btsPrecision)<<"  "<<(m_TRJBTS.btsFile_name)<<"  "<<std::endl;
    std::string state = "off";
    if (m_FrameTension.State == true)
    state = "on";
    statelog<<"Frame_Tension = "<<state<<"  "<<m_FrameTension.Type<<" "<<m_FrameTension.Tau<<"  "<<m_FrameTension.updatePeriod<<std::endl;
    state = "off";
    if (m_VolumeConstraint.State == true)
    state = "on";
    statelog<<"Volume_Constraint = "<<state<<"  "<<(m_VolumeConstraint.EQSteps)<<" "<<(m_VolumeConstraint.DeltaP)<<"  "<<(m_VolumeConstraint.K)<<"  "<<(m_VolumeConstraint.targetV)<<std::endl;
    
    //statelog<<" OpenEdgeEvolutionWithConstantVertex = .... "<<std::endl;
     //if(m_SpringPotentialBetweenTwoGroups->GetState()==true)
    //if (m_CoupleGCurvature->GetState() == true)
   // if (m_RigidWallCoupling->GetState() == true)
   
    statelog.close();
}

