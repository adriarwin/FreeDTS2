

#ifdef _OPENMP
# include <omp.h>
#endif
#include <thread>
#include "MC_Simulation_B.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "RNG.h"
#include "GenerateCNTCells.h"
#include "Curvature.h"
#include "Energy.h"
#include "LinkFlipMC.h"
#include "WritevtuFiles.h"
#include "Restart.h"
#include "BTSFile.h"
#include "Traj_XXX.h"
#include "CoupleToWallPotential.h"
#include "ActiveTwoStateInclusion.h"


/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 MC simulation class, runs mc simulation if it is defined in the input file.
 */
MC_Simulation_B::MC_Simulation_B(State *pState): m_pMESH(pState->m_pMesh), m_pActiveT((pState->m_pMesh)->m_pActiveT), m_pActiveL((pState->m_pMesh)->m_pActiveL), m_pInclusions((pState->m_pMesh)->m_pInclusion),    m_pHalfLinks1((pState->m_pMesh)->m_pHL),m_pHalfLinks2((pState->m_pMesh)->m_pMHL), m_pActiveV((pState->m_pMesh)->m_pActiveV),m_pSurfV((pState->m_pMesh)->m_pSurfV),m_pEdgeV((pState->m_pMesh)->m_pEdgeV),m_pEdgeL((pState->m_pMesh)->m_pEdgeL)

{
    
    
    
    double simtime = clock();
    m_pState =  pState,
#if TEST_MODE == Enabled
    std::cout<<"----> Note: Simulation has started: MC type B  -- : "<<std::endl;
#endif


///=======
//== Read MC_Simulation_B variables
//===========
std::cout<<std::fixed;
std::cout<<std::setprecision( Precision );
Nfunction f;
double R=pState->m_R_Vertex;   /// Move Vertex  size
double RB=pState->m_R_Box;   /// box move size
int ini=pState->m_Initial_Step;  	// initial step for the Simulation; usually zero
int final=pState->m_Final_Step;  	// final step for the Simulation
int box_centering_f = pState->m_Centering; // box centering
    
m_minAngle = pState->m_MinFaceAngle;
m_Lmin2    = pState->m_MinVerticesDistanceSquare;
m_Lmax2    = pState->m_MaxLinkLengthSquare;
int displaywrite=pState->m_Display_periodic;
std::string gfilename=pState->m_GeneralOutputFilename;
bool Targeted_State = pState->m_Targeted_State;
#if DEBUG_MODE == Enabled
    std::cout<<" simulation: reading varaible \n";
#endif
    
    
#ifdef _OPENMP
    if(Targeted_State==true)
    {
        std::cout<<"targeted thread---> id "<<omp_get_thread_num()<<" beta "<<(pState->m_Beta)<<"\n";
    }
#endif
    
//== Reading from the mesh
    m_pBox     = (pState->m_pMesh)->m_pBox;
    
    
   // std::cout<<" active trinagles sim class "<<m_pActiveT<<"\n";
    
    
   

    
    //========= VTU files at the begining. In case there will be some error later, so we can see how the file looklike.
        WritevtuFiles VTU(pState);
    
        std::string file="conf-1.vtu";
        VTU.Writevtu((m_pMESH->m_pActiveV),(m_pMESH->m_pActiveT),(m_pMESH->m_pHL),file);
    
    if(pState->m_IndexFile == true)
        if(!ReadIndexFile(pState->m_IndexFileName)){
            exit(1);

        }
//=========================================================================================
//======== We now create CNT object so that we can generate the cells anytime we want
//============================================================================================
GenerateCNTCells CNT((m_pMESH->m_pActiveV),pState->m_CNTCELL,m_pBox);
CNT.Generate();
#if DEBUG_MODE == Enabled
    std::cout<<" simulation: creating cnt \n";
#endif
RNG Random1(pState->m_Seed);

//========== Update the Mesh geometry variables ===============
//======== Prepare the trinagles: calculate area and normal vector
for (std::vector<triangle *>::iterator it = (m_pMESH->m_pActiveT).begin() ; it != (m_pMESH->m_pActiveT).end(); ++it)
(*it)->UpdateNormal_Area(m_pBox);
    
    
#if DEBUG_MODE == Enabled
    std::cout<<" simulation: tri normal \n";
#endif
//===== Prepare links:  normal vector and shape operator
for (std::vector<links *>::iterator it = (m_pMESH->m_pHL).begin() ; it != (m_pMESH->m_pHL).end(); ++it)
{
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
}
    
#if DEBUG_MODE == Enabled
    std::cout<<" simulation: link normal \n";
#endif
//======= Prepare vertex:  area and normal vector and curvature of surface vertices not the edge one
for (std::vector<vertex *>::iterator it = (m_pMESH->m_pSurfV).begin() ; it != (m_pMESH->m_pSurfV).end(); ++it)
    (pState->CurvatureCalculator())->SurfVertexCurvature(*it);
    
//====== edge links should be updated
for (std::vector<links *>::iterator it = (m_pMESH->m_pEdgeL).begin() ; it != (m_pMESH->m_pEdgeL).end(); ++it)
        (*it)->UpdateEdgeVector(m_pBox);

for (std::vector<vertex *>::iterator it = (m_pMESH->m_pEdgeV).begin() ; it != (m_pMESH->m_pEdgeV).end(); ++it)
        (pState->CurvatureCalculator())->EdgeVertexCurvature(*it);
    
#if DEBUG_MODE == Enabled
    std::cout<<" simulation: v curvature \n";
#endif
//=========================================================================================
// =================   Calculate Energy at the start of the simulation
//========================================================================================
Energy*   pEnergyCalculator = pState->GetEnergyCalculator();
double *tot_Energy = &(pState->m_TotEnergy);
*tot_Energy=pEnergyCalculator->TotalEnergy(m_pMESH->m_pSurfV,m_pMESH->m_pHL);
*tot_Energy=*tot_Energy+ pEnergyCalculator->TotalEnergy(m_pMESH->m_pEdgeV,m_pMESH->m_pEdgeL);

#if DEBUG_MODE == Enabled
    std::cout<<" simulation: system energy \n";
#endif
#if TEST_MODE == Enabled
#pragma omp critical
std::cout<<"----> Note: Total energy of the starting configuration is: "<<*tot_Energy<<std::endl;
#endif


#if TEST_MODE == Enabled
#pragma omp critical
std::cout<<"----> We printed our first configuration into a vtu file "<<std::endl;
#endif
    
    //========= write a vtu file after the update of the curvature and etc...
    // ===== this becames important for when the inclusion has local orinatation and it will differ from conf-1
    if(Targeted_State==true)
    {
    std::string file0="conf0.vtu";
    VTU.Writevtu((m_pMESH->m_pActiveV),(m_pMESH->m_pActiveT),(m_pMESH->m_pHL),file0);
    }
//=====================================================================================
//=====================================================================================
//================================ MC_Simulation_B Starts ================================
//=====================================================================================
//=====================================================================================
double mc_linkflip = (pState->m_MCMove).LinkFlip;
double mc_vertexmove = (pState->m_MCMove).VertexMove;
double mc_inclusionmove_kawa = (pState->m_MCMove).InclusionMove_Kawasaki;
double mc_inclusionmove_angle = (pState->m_MCMove).InclusionMove_Angle;

#if DEBUG_MODE == Enabled
    std::cout<<" simulation: getting move objects \n";
#endif
if(mc_linkflip == true || mc_vertexmove == true)
{
    if (CheckMesh(pState->m_pMesh)==false)
    {
        #pragma omp critical
        std::cout<<"----> Error: The mesh is bad for vertex and link flip move "<<std::endl;
        std::cout<<"either turn off this moves or use a better mesh "<<std::endl;
        exit(0);
    }
    else
    {
        #pragma omp critical
        std::cout<<"----> note: The mesh looks good, link flip and vertex move can be performed  "<<std::endl;
    }
}
int TotNoVertex=m_pSurfV.size();
    
//------> dynamic box algorithm
    DynamicBox *mc_box = pState->GetDynamicBox();
    mc_box->initialize();
    int BoxChangePeriodTau = mc_box->GetTau();
//----> dynamic topology algorithm -----------
    DynamicTopology *mc_topology = pState->GetDynamicTopology();
    mc_topology->initialize();
//------> edge treatment -----
    OpenEdgeEvolution * mc_edge_evo = pState->GetOpenEdgeEvolution();
    mc_edge_evo->Initialize();
//----> voulme coupling
    VolumeCoupling *volumecoupling = pState->GetVolumeCoupling();
    volumecoupling->Initialize(m_pActiveT);
//---> global curvature coupling
    GlobalCurvature * gc_globalcurvature = pState->GetGlobalCurvature();
    gc_globalcurvature->Initialize(m_pActiveV);
    
    CoupleToWallPotential * mc_rigidwall = pState->GetRigidWallCoupling();
    if(mc_rigidwall->GetState()==true)
    mc_rigidwall->Initialize(m_pActiveV);

    SpringPotentialBetweenTwoGroups * harmonic2groups = pState->Get2GroupHarmonic();
    if(harmonic2groups->GetState()==true)
    harmonic2groups->MakeGroups(m_pActiveV, pState->m_IndexFileName);
    if((pState->m_pConstant_NematicForce)->m_F0!=0)
    (pState->m_pConstant_NematicForce)->Initialize();
    if((pState->GetApply_Constant_Area())->GetState()==true)
    (pState->GetApply_Constant_Area())->Initialize(m_pActiveT);
    ActiveTwoStateInclusion *ActiveTwoState = pState->GetActiveTwoStateInclusion();
    if(ActiveTwoState->GetState()==true)
    ActiveTwoState->Initialize(m_pInclusions, ((pState->m_pMesh)->m_pInclusionType), pState->m_pinc_ForceField, &Random1);
    LinkFlipMC *mc_LFlip = pState->GetMCMoveLinkFlip();
    VertexMCMove *mc_VMove = pState->GetMCAVertexMove();

    EdgeVertexMCMove *mc_EdgeVMove = pState->GetMCEdgeVertexMove();
    InclusionMCMove *mc_IncMove = pState->GetInclusionMCMove();
    Restart *pRestart = pState->GetRestart();
    BTSFile btsFile ((pState->m_TRJBTS).btsFile_name, (pState->m_RESTART).restartState, "w");
    Traj_XXX TSIFile (m_pBox,(pState->m_TRJTSI).tsiFolder_name);

//=========== Some rate variables
    int totallmove = 0;
    int totalvmove = 0;
    int totalboxmove = 0;
    int InRate = 0;
    int VRate = 0;
    int LRate = 0;
    int boxrate = 0;
    int totalinmove = 0;
//============================= Energy Files =================================
//---> time seri file open: energy file
    m_pState->m_pTimeSeriesDataOutput->OpenTimeSeriesDataOutput();


//======================
#if TEST_MODE == Enabled
    #pragma omp critical
    std::cout<<"----> note: next is the first step in the sim loop "<<std::endl;
#endif
    #pragma omp critical
    std::cout<<"==============================================================  \n";
    std::cout<<"---> Running simulation from  "<<ini<<" to "<<final<<std::endl;
for (int mcstep=ini;mcstep<final+1;mcstep++)
{
    double VerexORbox=0;
    if( BoxChangePeriodTau!=0){
        VerexORbox=Random1.UniformRNG(1.0)*double(BoxChangePeriodTau);
        VerexORbox=1.0/VerexORbox;
    }
    if(VerexORbox<1){

            //== do link flip if it is allowed
            int no_link = m_pMESH->m_pHL.size();
            int no_link_iter = static_cast<int>(mc_linkflip * m_pMESH->m_pHL.size());
            if(no_link_iter!=0)
            for(int t=0;t<no_link;++t){
                int m=Random1.IntRNG(no_link);
                links *Tlinks = (m_pMESH->m_pHL)[m];  // chose a link randomly
                if(Tlinks->GetMirrorFlag()){
                    double thermal=Random1.UniformRNG(1.0);
                    mc_LFlip->MC_FlipALink(mcstep,Tlinks,thermal);
                    LRate+=mc_LFlip->GetMoveValidity();
                    ++totallmove;
                }
            }
        
        //== do vertex move
            int no_surfV = (m_pMESH->m_pSurfV).size();
            int no_surfV_iter = int(mc_vertexmove*(m_pMESH->m_pSurfV).size());

            for(int t=0;t<no_surfV_iter;++t)
            {
                int n=Random1.IntRNG(no_surfV);
                vertex *lpvertex = (m_pMESH->m_pSurfV)[n];   //
                if(lpvertex->GetGroupName()!=pState->m_FreezGroupName)
                {
                double dx=1-2*Random1.UniformRNG(1.0);            // Inside a cube with the side length of R
                double dy=1-2*Random1.UniformRNG(1.0);
                double dz=1-2*Random1.UniformRNG(1.0);
                double thermal=Random1.UniformRNG(1.0);
                bool cwp = true;
                if(mc_rigidwall->GetState()==true)
                cwp=mc_rigidwall->CheckVertexMoveWithinWalls(mcstep,R*dx,R*dy,R*dz,lpvertex);

                if(cwp==true ){
                    mc_VMove->MC_MoveAVertex(mcstep,lpvertex,R*dx,R*dy,R*dz,thermal);
                    VRate+=mc_VMove->GetMoveValidity();
                    totalvmove++;
                }
            }
        }// end of if(mc_vertexmove == true)
        
        //== do edge v move
        int no_edgeV = int(((pState->m_MCMove).EdgeVertexMove)*((m_pMESH->m_pEdgeV).size()));
        for(int t=0;t<no_edgeV;t++){
            int n=Random1.IntRNG((m_pMESH->m_pEdgeV).size());
            vertex *lpvertex = (m_pMESH->m_pEdgeV)[n];   //
            if(lpvertex->GetGroupName()!=pState->m_FreezGroupName){
                double dx=1-2*Random1.UniformRNG(1.0);            // Inside a cube with the side length of R
                double dy=1-2*Random1.UniformRNG(1.0);
                double dz=1-2*Random1.UniformRNG(1.0);
                double thermal=Random1.UniformRNG(1.0);
                bool cwp = true;
                if(mc_rigidwall->GetState()==true) // check if the move is valid within the defined wall type
                cwp=mc_rigidwall->CheckVertexMoveWithinWalls(mcstep,R*dx,R*dy,R*dz,lpvertex);
                if(cwp==true )
                {
                    mc_EdgeVMove->MC_MoveAVertex(mcstep,lpvertex,R*dx,R*dy,R*dz,thermal);
                   // VRate =mc_EdgeVMove->GetMoveValidity();
                    //totalvmove++;
                }
           }
        
        }// end if(edgevertexmove == true)
        
        int no_kawa = int(mc_inclusionmove_kawa*m_pInclusions.size());
        for(int t=0;t<no_kawa;t++){
            int n=Random1.IntRNG(m_pInclusions.size());
            inclusion *linclusion = m_pInclusions.at(n);   //
            mc_IncMove->MC_Move_AnInclusion(linclusion,&Random1, 1);
            InRate+=mc_IncMove->GetMoveValidity();
            totalinmove++;
        }
        int no_angle = int(mc_inclusionmove_angle*m_pInclusions.size());
        for(int t=0;t<no_angle;t++){
            int n=Random1.IntRNG(m_pInclusions.size());
            inclusion *linclusion = m_pInclusions.at(n);   //
            mc_IncMove->MC_Move_AnInclusion(linclusion,&Random1, 2);
            InRate+=mc_IncMove->GetMoveValidity();
            totalinmove++;
        }

    }
    else
    { // box move
        
        double dr=1-2*Random1.UniformRNG(1.0);
        double thermal=Random1.UniformRNG(1.0);
        bool move = mc_box->MCMoveBoxChange(RB*dr, tot_Energy, thermal, mcstep, (&CNT));
        if(move==true)
            boxrate=boxrate+1;
        if(mc_box->GetCNTCondition()==false)
            CNT.Generate();
           totalboxmove++;
        
        
        //=== attempt to change topology

    }
//==============  Edge evolotion
    if(mc_edge_evo->GetRate()!=0 && mcstep%(mc_edge_evo->GetRate())==0)
    {
        mc_edge_evo->MC_Move(&Random1,m_Lmin2,m_Lmax2,m_minAngle);
    }
    
    //== change topology
    {
        double thermal=Random1.UniformRNG(1.0);
        bool top_move = mc_topology->MCMove(mcstep, tot_Energy, &Random1, (&CNT));
    }

//==== Active Inclsuion exchange two state; as these movie are indepenednt of the energy based move, the have no conditions
//===================================
    if(ActiveTwoState->GetState() == true)
    ActiveTwoState->ActiveExchange(tot_Energy);

//==================== centering the mesh in the box =============================================
    if (box_centering_f != 0 && mcstep % box_centering_f == 0) {
        (pState->m_pMesh)->CenterMesh(); // Center the mesh in the box
        CNT.Generate(); // Update CNT cells due to vertex location changes
    }
//==================== We are now writing some dynamics files
if(Targeted_State==true)
if(mcstep%displaywrite==0 && displaywrite!=0)
{
        file="conf"+f.Int_to_String(int(mcstep/displaywrite))+".vtu";
        VTU.Writevtu((m_pMESH->m_pActiveV),(m_pMESH->m_pActiveT),(m_pMESH->m_pHL),file);
}
if((pState->m_TRJBTS).btsPeriod!=0 && mcstep%((pState->m_TRJBTS).btsPeriod)==0 )
{
    btsFile.WrireBTSFile(mcstep,  (pState->m_pMesh));
}
if((pState->m_TRJTSI).tsiPeriod!=0 && mcstep%((pState->m_TRJTSI).tsiPeriod)==0 )
{
    std::string file=gfilename+f.Int_to_String(int(mcstep/(pState->m_TRJTSI).tsiPeriod))+"."+TSIExt;
   TSIFile.WriteTSI(mcstep,file, m_pActiveV, m_pActiveT, m_pInclusions);
}
//--> write energy and other info into the file output-en.xvg
    m_pState->m_pTimeSeriesDataOutput->WrireTimeSeriesDataOutput(mcstep);

if( pState->m_pTimeSeriesDataOutput->GetPeriodic()!=0  && mcstep%(pState->m_pTimeSeriesDataOutput->GetPeriodic())==0)
{
    if(ActiveTwoState->GetState() == true){
        // I do not know what this suppose to do
        *(ActiveTwoState->GetActiveEnergy()) = 0;
    }
    //======= write acceptance rate
    if(Targeted_State==true)
        {
            std::cout<<"Step = "<<mcstep<<"/"<<final<<std::flush;
            std::cout << std::fixed << std::setprecision(2);
            std::cout<<"; Accpetance rates="<<std::flush;
        if(mc_vertexmove >0)
            std::cout<<" Postion update: "<<double(VRate)/double(totalvmove)*100.0<<" % "<<std::flush;
        if(mc_linkflip >0)
            std::cout<<"; Alexander move: "<<LRate/double(totallmove)*100.0<<" % "<<std::flush;
        if(m_pInclusions.size()!=0)
            std::cout<<"Incluions move: "<<InRate/double(totalinmove)*100.0<<" %. "<<std::flush;
        if(BoxChangePeriodTau!=0)
            std::cout<<"; Box size update:  "<<boxrate/double(totalboxmove)*100.0<<" % "<<std::flush;
            std::cout << '\r';
            std::cout << "\033[K";

        }
}
if(Targeted_State==true)
if((pState->m_RESTART).restartPeriod!=0 && mcstep%((pState->m_RESTART).restartPeriod)==0){
    
        m_pState->m_pTimeSeriesDataOutput->FlushTimeSeriesFile();
        pRestart->WrireRestart(mcstep,gfilename,(pState->m_pMesh),R,RB);
}
//========== Optimiszing R and RB
if(mcstep%500==0) // Optimize R and RB
{
    double vrate=double(VRate)/double(totalvmove);
    double Boxrate=double(boxrate)/double(totalboxmove);
        if(mcstep<double(final)*0.01 && mcstep<10000){
            if(mc_vertexmove == true){
                if(vrate<0.4){
                    if(R>0.01)
                        R=R/1.1;
                }
                else if(vrate>0.6){
                    if(R<0.2)
                        R=1.1*R;
                }
            }
        if(BoxChangePeriodTau !=0){
            if(Boxrate<0.4){
                if(RB>0.01)
                RB=RB/1.1;
            }
            else if(Boxrate>0.6){
                if(RB<0.2)
                RB=1.1*RB;
            }
         }
        }
        
        VRate=0;
        totallmove=0;
        LRate=0;
        totalvmove=0;
        InRate=0;
        totalinmove=0;
        boxrate=0;
        totalboxmove=0;
    
}// Optimize R and RB
}//End of Sim Loop for (int i=ini;i<final+1;i++)
//================ checking for energy leak, both should be equal

    
    DetailedSystemEnergy();
    double fen = SystemEnergy();
    std::cout<<std::setprecision(8)<<fen<<"  "<<*tot_Energy<<std::endl;
    simtime=clock()-simtime;
    #pragma omp critical
    std::cout<<" total time "<<((float)simtime)/CLOCKS_PER_SEC<<" second \n";

    

}//End of object constructor
//======================================================================================================
//======================================================================================================
//=========== The MC simulation constructor is over, the rest are accessory functions ==================
//======================================================================================================
//======================================================================================================
MC_Simulation_B::~MC_Simulation_B()
{
}
//================A function to calculate energy of the mesh from scratch  ==================================
double  MC_Simulation_B::SystemEnergy()
{
    double en = 0;
    for (std::vector<vertex *>::iterator it1 = (m_pMESH->m_pEdgeV).begin() ; it1 != (m_pMESH->m_pEdgeV).end(); ++it1)
    {
        double x=(*it1)->GetVXPos();
        double y=(*it1)->GetVYPos();
        double z=(*it1)->GetVZPos();
        (*it1)->UpdateVXPos(x);
        (*it1)->UpdateVYPos(y);
        (*it1)->UpdateVZPos(z);
    }
    
    for (std::vector<triangle *>::iterator it = (m_pMESH->m_pActiveT).begin() ; it != (m_pMESH->m_pActiveT).end(); ++it){
        double oldarea = (*it)->GetArea();
        (*it)->UpdateNormal_Area(m_pBox);
        double newarea = (*it)->GetArea();
        if(fabs(newarea-oldarea)>0.0002)
            std::cout<<"  area is different \n";
    }

    //===== Prepare links:  normal vector and shape operator
         for (std::vector<links *>::iterator it = (m_pMESH->m_pHL).begin() ; it != (m_pMESH->m_pHL).end(); ++it)
    {
        Vec3D oldN = (*it)->GetNormal();
              (*it)->UpdateNormal();
              (*it)->UpdateShapeOperator(m_pBox);
        Vec3D newN = (*it)->GetNormal()-oldN;
         if(fabs(newN(0))>0.001 || fabs(newN(1))>0.001 || fabs(newN(2))>0.001)
             std::cout<<"  normal is different "<<(*it)->GetID()<<"  "<<(*it)->GetMirrorLink()->GetID()<<"\n";

    }

    //======= Prepare vertex:  area and normal vector and curvature of surface vertices not the edge one
      for (std::vector<vertex *>::iterator it = (m_pMESH->m_pSurfV).begin() ; it != (m_pMESH->m_pSurfV).end(); ++it)
        (m_pState->CurvatureCalculator())->SurfVertexCurvature(*it);
        
    //====== edge links should be updated
      for (std::vector<links *>::iterator it = (m_pMESH->m_pEdgeL).begin() ; it != (m_pMESH->m_pEdgeL).end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = (m_pMESH->m_pEdgeV).begin() ; it != (m_pMESH->m_pEdgeV).end(); ++it)
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(*it);

    Energy*   pEnergyCalculator = m_pState->GetEnergyCalculator();
    en=pEnergyCalculator->TotalEnergy(m_pMESH->m_pSurfV,m_pMESH->m_pHL);
    en+=pEnergyCalculator->TotalEnergy(m_pMESH->m_pEdgeV,m_pMESH->m_pEdgeL);
    
    return en;
}
bool MC_Simulation_B::CheckMesh(MESH *pMesh)
{
    Vec3D *pBox     = pMesh->m_pBox;
    std::vector<vertex*> pV= pMesh->m_pActiveV;
    std::vector<links*> pL=pMesh->m_pHL;
    // Check if there are any pair of vertices that are to close

    for (int i=0;i<pV.size();i++)
    for (int j=i+1;j<pV.size();j++)
    {
        double l2 = CheckLengthBetweenTwoVertex(pV[i],pV[j], pBox);
        if(l2<m_Lmin2)
        return false;
    }
    // Check if there are any pair of connected vertices that are to far
    for (std::vector<vertex *>::iterator it = pV.begin() ; it != pV.end(); ++it)
    {
        std::vector <vertex *> NV = (*it)->GetVNeighbourVertex();
        for (std::vector<vertex *>::iterator it2 = NV.begin() ; it2 != NV.end(); ++it2)
        {
            double l2 = CheckLengthBetweenTwoVertex(*it, *it2, pBox);
            if(l2>m_Lmax2)
            return false;
        }
    }

    // Check the angle of the faces and see if the are bent to much
    for (std::vector<links *>::iterator it = pL.begin() ; it != pL.end(); ++it)
    {
        double face  = CheckFaceAngle(*it, pBox);
        if(face<m_minAngle)
            return false;
    }
    return true;
}
double MC_Simulation_B::CheckLengthBetweenTwoVertex(vertex* v1, vertex* v2, Vec3D *pBox)
{
    double x2=v2->GetVXPos();
    double y2=v2->GetVYPos();
    double z2=v2->GetVZPos();
    double x1=v1->GetVXPos();
    double y1=v1->GetVYPos();
    double z1=v1->GetVZPos();
    double dr[3];
    dr[0]=x2-x1;
    dr[1]=y2-y1;
    dr[2]=z2-z1;
    for (int i=0;i<3;i++)
    if(fabs(dr[i])>(*pBox)(i)/2.0)
    {
        if(dr[i]<0)
            dr[i]=(*pBox)(i)+dr[i];
        else if(dr[i]>0)
            dr[i]=dr[i]-(*pBox)(i);
    }
    double l2 = 0;
    for (int i=0;i<3;i++)
    l2+=dr[i]*dr[i];
    return l2;
}
double   MC_Simulation_B::CheckFaceAngle(links * l, Vec3D *pBox)
{
    vertex* v1 = l->GetV1();
    vertex* v2 = l->GetV2();
    vertex* v3 = l->GetV3();
    vertex* v4 = (l->GetMirrorLink())->GetV3();
    Vec3D   N1=CalculateNormal(v1,v2,v3, pBox);
    Vec3D   N2=CalculateNormal(v1,v4,v2, pBox);
    
    double faceangle =  N1.dot(N1,N2) ;
    return faceangle;
}
Vec3D   MC_Simulation_B::CalculateNormal(vertex* v1 ,vertex* v2 ,vertex* v3, Vec3D *pBox)
{
    Vec3D X1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
    Vec3D X2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
    Vec3D X3(v3->GetVXPos(),v3->GetVYPos(),v3->GetVZPos());

   Vec3D dX1=X2-X1;
    for (int i=0;i<3;i++)
    if(fabs(dX1(i))>(*pBox)(i)/2.0)
    {
        if(dX1(i)<0)
            dX1(i)=(*pBox)(i)+dX1(i);
        else if(dX1(i)>0)
            dX1(i)=dX1(i)-(*pBox)(i);
    }
    Vec3D dX2=X3-X1;
    for (int i=0;i<3;i++)
        if(fabs(dX2(i))>(*pBox)(i)/2.0)
        {
            if(dX2(i)<0)
                dX2(i) = (*pBox)(i)+dX2(i);
            else if(dX2(i)>0)
                dX2(i)=dX2(i)-(*pBox)(i);
        }
    Vec3D N=dX1*dX2;
    double area=N.norm();
    N=N*(1.0/area);
    
    return N;
}
//--- this function reads the index file and assigned the group name to the vertex
bool MC_Simulation_B::ReadIndexFile(std::string filename) {
    std::ifstream indexfile(filename);
    if (!indexfile) {
        std::cerr << "---> error: Unable to open index file " << filename << std::endl;
        return false;
    }

    std::string name;
    int NAtom;
    int gid = 1;
    while (indexfile >> name >> NAtom) {
        int vid;
        for (int i = 0; i < NAtom; ++i) {
            if (!(indexfile >> vid)) {
                std::cerr << "Error: Failed to read ID from index file" << std::endl;
                return false;
            }
            if (vid < m_pActiveV.size()) {
                m_pActiveV[vid]->UpdateGroupName(name);
                m_pActiveV[vid]->UpdateGroup(gid);
            } else {
                std::cerr << "Error: ID " << vid << " exceeds active vector size" << std::endl;
                return false;
            }
        }
        gid++;
    }
    indexfile.close();
    
    return true;
}
//--- end of reading index file
void MC_Simulation_B::DetailedSystemEnergy(){
    std::vector<double> energyv,energyl;
    for (int i=0;i<m_pActiveL.size();i++){
        energyl.push_back((m_pActiveL[i])->GetIntEnergy());
    }
    for (int i=0;i<m_pActiveV.size();i++){
        energyv.push_back((m_pActiveV[i])->GetEnergy());
    }
double fen = SystemEnergy();
for (int i=0;i<(m_pMESH->m_pActiveL).size();i++){
    if(fabs((m_pActiveL[i])->GetIntEnergy()-energyl[i])>0.000001){
        std::cout<<"int energy: "<<std::setprecision(5)<<(m_pActiveL[i])->GetIntEnergy()<<"  "<<energyl[i]<<" vid "<<((m_pActiveL[i])->GetV1())->GetVID()<<"  ";
    std::cout<<((m_pActiveL[i])->GetV2())->GetVID()<<"  "<<((m_pActiveL[i])->GetV3())->GetVID()<<" \n ";
    }
}
for (int i=0;i<(m_pMESH->m_pActiveV).size();i++){
        if(fabs((m_pActiveV[i])->GetEnergy()-energyv[i])>0.000001){
            std::cout<<"v energy: "<<std::setprecision(5)<<(m_pActiveV[i])->GetEnergy()<<"  "<<energyv[i]<<" vid "<<((m_pActiveV[i]))->GetVID()<<"  ";
        }
    }

    return;
}












