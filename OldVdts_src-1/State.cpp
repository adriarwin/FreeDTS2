

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
    m_pTimeSeriesDataOutput             = new TimeSeriesDataOutput(this);
    m_pDynamicBox                       = new NoBoxChange;
    m_pDynamicTopology                  = new ConstantTopology;
    m_pOpenEdgeEvolution                = new NoEvolution;
    m_pVolumeCoupling                   = new NoCoupling;
    m_pCoupleGlobalCurvature            = new NoGlobalCurvature;
    m_pTotalAreaCoupling                = new NoTotalAreaCoupling;
    m_pBoundary                         = new PBCBoundary;
    m_pConstant_NematicForce = &m_Constant_NematicForce;
    m_Argument = argument;
    m_Targeted_State = true;
    m_Healthy =true;
    m_Integrator = "MC_B";
    m_Total_no_Threads = 1;
    m_Initial_Step = 0;
    m_Final_Step = 0;
    m_Seed =36723;
    m_R_Vertex=0.05;   // Move Vertex  within a box with this size
    m_R_Box=0.04;   // box change within this range
    m_TotEnergy = 0;
    m_Mem_Spontaneous_Curvature = 0.0;
    m_TopologyFile = "topology.top";
    m_InputFileName = "Input.dts";
    m_IndexFileName = "Index.inx";
    m_IndexFile = false;
    m_GeneralOutputFilename = "dts";
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
    m_MCMove.VertexMove = 1;
    m_MCMove.EdgeVertexMove = 1;
    m_MCMove.LinkFlip = 1;
    m_MCMove.InclusionMove_Angle = 1;
    m_MCMove.InclusionMove_Kawasaki = 1;
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
    WriteStateLog();        // writing a log file
    

    //===============================================================================================
    //============= Call different objects for their tasks based on the updated inputs ==============
    //===============================================================================================
    
    //------> Creating the blueprint of the mesh and inclusions and etc .....
    CreateMashBluePrint BluePrint;
    // if restart is false
    bool readingfromrestart = false;
    MeshBluePrint meshblueprint;
    if(m_RESTART.restartState == true)  // reading restart file and update the [State] object to that state
    {
        std::cout<<m_RESTART.restartFilename<<" restart file name is "<<std::endl;
        meshblueprint = m_Restart.ReadRestart(m_RESTART.restartFilename,&readingfromrestart);
        if(readingfromrestart==true)
            std::cout<<"-----> Note: Restart file was successfully read  "<<std::endl;
        else
            m_RESTART.restartState = false;
    }
    if(readingfromrestart==false) // if reading restart fails, we just start from begining
    {
        meshblueprint = BluePrint.MashBluePrintFromInput_Top(m_InputFileName,m_TopologyFile);
        std::cout<<"---> note: Mesh was taken from the topology  "<<std::endl;
    }
        m_Mesh.GenerateMesh(meshblueprint);
    
    // Generating couplings
    LinkFlipMC LPTEM (this);
    m_LinkFlipMC = LPTEM;
    VertexMCMove TVM (this);
    m_VertexMoveMC = TVM;
    EdgeVertexMCMove EVM (this);
    m_EdgeVertexMoveMC = EVM;
    InclusionMCMove TIM(this);
    m_IncMove = TIM;
    Energy Ten(&m_inc_ForceField);  // perhaps later could have a state pointer as well
    m_EnergyCalculator = Ten;
    // this should be called after energy and curvature is called 

    if(m_STRUC_ActiveTwoStateInclusion.state==true)
    {
         double ep1 = m_STRUC_ActiveTwoStateInclusion.ep1;
         double ep2 = m_STRUC_ActiveTwoStateInclusion.ep2;
         double per = m_STRUC_ActiveTwoStateInclusion.persentage;
         double gama = m_STRUC_ActiveTwoStateInclusion.gama;
         std::string t1 = m_STRUC_ActiveTwoStateInclusion.nametype1;
         std::string t2 = m_STRUC_ActiveTwoStateInclusion.nametype2;
         ActiveTwoStateInclusion TA2I(true,ep1,ep2, per,gama, t1, t2);
         m_ActiveTwoStateInclusion = TA2I;

    }
    



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
        else if(Arg1=="-ndx"){   // get the index file and also check if the file exist and it has inx extension
            m_IndexFileName = m_Argument[i+1];
            m_IndexFile = true;
            if(f.SubstringFromRight(m_IndexFileName,'.')!="inx"){
                m_IndexFileName = m_IndexFileName+".inx";
            }
            if(f.FileExist(m_IndexFileName)==false){
                        std::cout<<"---> error: the index file name with the name "<<m_IndexFileName<<" does not exist"<<std::endl;
                        exit(1);
            }

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
        else if(firstword == "Integrator")
        {
            input>>str>>m_Integrator;
            getline(input,rest);
        }
        else if(firstword == "MC_Moves")
        {
            getline(input, str);
            std::vector<std::string> data = f.split(str);

            // Check if the correct number of moves is provided
            const int expectedMoves = 5;
            if (data.size() != expectedMoves + 1) { // +1 to account for the first element is just = sign
                std::cout << "---> error: the defined MC moves does not cover all the moves. It should be " << expectedMoves << " numbers" << std::endl;
                exit(0);
            }

            // Assign normalized move rates
            m_MCMove.VertexMove = f.String_to_Double(data[1]);
            m_MCMove.LinkFlip = f.String_to_Double(data[2]);
            m_MCMove.EdgeVertexMove = f.String_to_Double(data[3]);
            m_MCMove.InclusionMove_Kawasaki = f.String_to_Double(data[4]);
            m_MCMove.InclusionMove_Angle = f.String_to_Double(data[5]);

            // Check if move rates are non-negative
            if (m_MCMove.VertexMove < 0 || m_MCMove.LinkFlip < 0 || m_MCMove.InclusionMove_Angle < 0 || m_MCMove.InclusionMove_Kawasaki < 0 || m_MCMove.EdgeVertexMove < 0 ) {
                std::cout << "---> error: move rates cannot be negative" << std::endl;
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
        else if(firstword == "Dynamic_Box")
        {
            
            std::string type;
            int period = 0;
            double force = 0;

            input >> str >> type >> period >> force;

            if (type == "SingleSideConstantForce") {
                m_pDynamicBox = new DynamicBoxSide(period, force, this);
            }
            else if (type == "PositionRescale_FrameTension") {
                m_pDynamicBox = new PositionRescaleFrameTensionCoupling(period, force,this);
            }

            // Consume remaining input line
            getline(input, rest);
        }
        else if(firstword == "Dynamic_Topology")
        {
            
            std::string type;
            int period = 0;

            input >> str >> type >> period;

            if (type == "Scission_By3Edges") {
                
                m_pDynamicTopology = new Three_Edge_Scission(period, this);
            }
            // Consume remaining input line
            getline(input, rest);
        }
        else if(firstword == "OpenEdgeEvolution")
        {
            // OpenEdgeEvolution =  EvolutionWithConstantVertex PT_steps
            std::string type;
            int period = 0;
            input >> str >> type >> period;
            getline(input,rest);
            if (type == "EvolutionWithConstantVertex") {
                
                m_pOpenEdgeEvolution = new OpenEdgeEvolutionWithConstantVertex(period, this);
            }
            else if (type == "EvolutionWithConstantVertex_B") {
                
                std::cout<<"---> warning: this algorithm for edge treatment has errors. It is up to you to use it \n";
                m_pOpenEdgeEvolution = new OpenEdgeEvolutionWithConstantVertex_B(period, this);
            }
            else {
                
                std::cout<<"---> error: unknown Open Edge Evolution type: "<<type<<"\n";
            }
        }
        else if(firstword == "GlobalCurvature")
        {
            std::string type;
            double k,gc0;
            input>>str>>type>>k>>gc0;
            if(type == "CoupledToHarmonicPotential") {
                m_pCoupleGlobalCurvature = new CouplingtoFixedGlobalCurvature(k,gc0);
            }
            else {
                std::cout<<"---> error: this algorithm for global curvature is unknown \n";
                std::cout<<"----------- no global curvature coupling with be applied \n";
            }
            getline(input,rest);
            
        }
        else if(firstword == "Volume_Constraint")
        {
            
            std::string type;
            input >> str >> type;
            if(type == "SecondOrderCoupling"){
                int eqsteps = 0;
                double dp,k,vt;
                input>>eqsteps>>dp>>k>>vt;
                m_pVolumeCoupling = new CmdVolumeCouplingSecondOrder(eqsteps, dp,  k, vt);
            }
            else if(type == "Osmotic_Pressure"){
                int eqsteps = 0;
                double gamma,p0;
                input>>eqsteps>>gamma>>p0;
                m_pVolumeCoupling = new Apply_Osmotic_Pressure(eqsteps,gamma,p0);
            }
            getline(input,rest);
        }
        else if(firstword == "TotalAreaCoupling"){
            
            std::string type;
            int eq_steps;
            double gamma , k0;
            input>>str>>type;

            if(type=="CouplingTotalAreaToHarmonicPotential"){
                input>>eq_steps>>gamma>>k0;
                m_pTotalAreaCoupling = new CouplingTotalAreaToHarmonicPotential(eq_steps,gamma,k0);
            }
            getline(input,rest);
        }
        else if(firstword == "TotalAreaCoupling"){
            std::string type;
            input>>str>>type;
            if(type == "Spherical"){
                //m_pBoundary
            }
        }
        else if(firstword == "MinfaceAngle")
        {
            input>>str>>m_MinFaceAngle;
            getline(input,rest);
        }
        else if(firstword == "OutPutEnergy_periodic")
        {
            int period;
            input>>str>>period;
            m_pTimeSeriesDataOutput->UpdatePeriod(period);
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
        else if(firstword == "ConstantField")
        {
            double k,x,y,z;
            input>>str>>k>>x>>y>>z;
            getline(input,rest);
            double norm = sqrt(x*x+y*y+z*z);
            Vec3D tem(x,y,z);
            tem = tem*(1.0/norm);
            m_inc_ForceField.m_FieldDirection = tem;
            m_inc_ForceField.m_FieldStrength = k;
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
    statelog<<"MC_Moves = "<<(m_MCMove.VertexMove)<<"  "<<(m_MCMove.LinkFlip)<<"  "<<m_MCMove.EdgeVertexMove<<"  "<<(m_MCMove.InclusionMove_Angle)<<"  "<<(m_MCMove.InclusionMove_Kawasaki)<<std::endl;
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

    
    //statelog<<" OpenEdgeEvolutionWithConstantVertex = .... "<<std::endl;
     //if(m_SpringPotentialBetweenTwoGroups->GetState()==true)
    //if (m_CoupleGCurvature->GetState() == true)
   // if (m_RigidWallCoupling->GetState() == true)
   
    statelog.close();
}

