#if !defined(AFX_State_H_INCLUDED_)
#define AFX_State_H_INCLUDED_

#include "MESH.h"
#include "SimDef.h"
#include "CreateMashBluePrint.h"
#include "Inclusion_Interaction_Map.h"
#include "CouplingtoFixedGlobalCurvature.h"
#include "SpringPotentialBetweenTwoGroups.h"
#include "CmdVolumeCouplingSecondOrder.h"
#include "CoupleToWallPotential.h"
#include "LinkFlipMC.h"
#include "VertexMCMove.h"
#include "EdgeVertexMCMove.h"
#include "InclusionMCMove.h"
#include "Restart.h"
#include "Nfunction.h"
#include "ActiveTwoStateInclusion.h"
#include "Apply_Osmotic_Pressure.h"
#include "Apply_Constant_Area.h"
#include "Curvature.h"
#include "Energy.h"
#include "Constant_NematicForce.h"
#include "DynamicBox.h"
#include "PositionRescaleFrameTensionCoupling.h"
#include "DynamicBoxSide.h"
#include "DynamicTopology.h"
#include "Three_Edge_Scission.h"
#include "OpenEdgeEvolution.h"
#include "OpenEdgeEvolutionWithConstantVertex.h"
#include "OpenEdgeEvolutionWithConstantVertex_B.h"
#include "VolumeCoupling.h"

/*#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"*/


/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This class reads the inputs and distrubutes the tasks based on inputs provided and makes all the initials variables
 It is called by Job class: 
 */
struct STRUC_Membrane_Parameters {
    double lambda;
    double kappa_geo;
    double kappa_normal;
    double kappa; // this is read by another parameters for now; later may change
    double kappa_g; // this is read by another parameters for now; later may change
};
struct STRUC_TRJTSI {    // data structure for tsi trajectory file
    int tsiPeriod;
    int tsiPrecision;
    std::string tsiFolder_name;
    bool tsiState;
};
struct STRUC_RESTART {  // restart, how often write a file and check if this is a restart sim
    int  restartPeriod;
    bool restartState;
    std::string restartFilename;

};
struct STRUC_TRJBTS { //data structure for bts trajectory file (binary file format)
    int btsPeriod;
    int btsPrecision;
    std::string btsFile_name;
    bool btsState;
};
struct STRUC_ConstantArea {  // data structure for inputs of constant area algorithm
    bool State;
    int EQSteps;
    double Gamma;
    double K0;
    
};
struct STRUC_MCMOVES {  // data structure for turning on and off certain moves and specify the rate
    double VertexMove;
    double LinkFlip;
    double InclusionMove_Angle;
    double InclusionMove_Kawasaki;
    double EdgeVertexMove;


};
struct Parallel_Tempering {  // data structure for turning on and off certain moves
    bool State;
    int PT_steps;
    double PT_minbeta;
    double PT_maxbeta;
};
struct STRUC_ActiveTwoStateInclusion {
    bool state;
    std::string nametype1;
    std::string nametype2;
    double ep1 ;
    double ep2 ;
    double persentage;
    double gama;
};

class State
{
public:
    
	State(std::vector <std::string> argument);
    State();

	 ~State();
inline Curvature *CurvatureCalculator()                    {return &m_CurvatureCalculations;}
inline CouplingtoFixedGlobalCurvature *GetGlobalCurvature()                      {return &m_CoupleGCurvature;}
inline SpringPotentialBetweenTwoGroups *Get2GroupHarmonic()                      {return &m_SpringPotentialBetweenTwoGroups;}
inline PositionRescaleFrameTensionCoupling *GetRescaleTension()                  {return &m_RescaleTenCoupl;}
inline Apply_Constant_Area *GetApply_Constant_Area()                           {return &m_Apply_Constant_Area;}
inline LinkFlipMC *GetMCMoveLinkFlip()                         {return &m_LinkFlipMC;}
inline VertexMCMove *GetMCAVertexMove()                         {return &m_VertexMoveMC;}
inline EdgeVertexMCMove *GetMCEdgeVertexMove()                         {return &m_EdgeVertexMoveMC;}
inline InclusionMCMove *GetInclusionMCMove()                     {return &m_IncMove;}
inline Restart *GetRestart()                                    {return &m_Restart;}
inline CoupleToWallPotential *GetRigidWallCoupling()                                    {return &m_RigidWallCoupling;}
inline ActiveTwoStateInclusion *GetActiveTwoStateInclusion()                                    {return &m_ActiveTwoStateInclusion;}
inline Energy *GetEnergyCalculator()                                    {return &m_EnergyCalculator;}
inline DynamicBox *GetDynamicBox()                                    {return m_pDynamicBox;}
inline DynamicTopology *GetDynamicTopology()                    {return m_pDynamicTopology;}
inline OpenEdgeEvolution     *GetOpenEdgeEvolution()            {return m_pOpenEdgeEvolution;}
inline VolumeCoupling *GetVolumeCoupling()                         {return m_pVolumeCoupling;}

public:
    
    bool m_Healthy;   // To check if the input data are read correctly
    bool m_Targeted_State; // Only relavant for Parallel Tempering by 2022; Which state carries the target temparature 
    std::vector <std::string> m_Argument;
    int m_Initial_Step; // initial step
    int m_Final_Step; // final step
    int m_Centering; // box centering frequency; important for pressure coupling
    int m_Seed;       // seed for random number generator
    int m_Total_no_Threads;       // Total no of Threads
    std::string   m_TopologyFile; // name of the topology file, *.top, *.dat *.tsi *.bts
    std::string   m_InputFileName; // name of the topology file, *.top, *.dat *.tsi *.bts
    double m_Beta ;        // system temperature only applicable if Parallel Tempering Methods is on
    double m_Mem_Spontaneous_Curvature ; //Spontaneous Curvature of the membrane, same effect if the membrane is covered fully by proteins but faster
    double m_MinFaceAngle;          //  minimum angle between the face (smaller will results in error), this is the value of the cos
    double m_R_Vertex;   // Move Vertex  within a box with this size
    double m_R_Box;   // box change within this range
    double m_MinVerticesDistanceSquare; //  minimum distance allowed between two vertices  (smaller will results in error)
    double m_MaxLinkLengthSquare;       //  maximum distance allowed between two nighbouring vertices  (larger will results in error)
    std::string m_GeneralOutputFilename; //  a general file flag for specific run
    std::string m_IndexFileName;            // Name of the index file for group specification
    bool m_IndexFile;                       // to check if the index file is provided
    std::string m_FreezGroupName ;            // Name of a group to be frozen
    STRUC_TRJTSI m_TRJTSI;                  //  an object for tsi file format
    STRUC_TRJBTS m_TRJBTS;                  //  an object for binary  trajectory
    int  m_Display_periodic ;               //  periodic for paraview output
    int m_OutPutEnergy_periodic;            //  periodic for energy file
    Vec3D m_CNTCELL;                    // for domain decomposition
    std::string m_Integrator;               //  Type of integrator (for now only mc exist)
    STRUC_RESTART m_RESTART;                // To check if this is a restart simulation of fresh start
    STRUC_ConstantArea m_STRUC_ConstantArea;
    STRUC_MCMOVES m_MCMove;                 // data structure for turning on and off certain moves
    STRUC_ActiveTwoStateInclusion m_STRUC_ActiveTwoStateInclusion;  // input data to start active two state membrane
    Parallel_Tempering m_Parallel_Tempering; // an object that includes info about Parallel Tempering method that we are applying
    MESH          *m_pMesh;
    Inclusion_Interaction_Map *m_pinc_ForceField;

    double m_TotEnergy;
    CouplingtoFixedGlobalCurvature  m_CoupleGCurvature;
    SpringPotentialBetweenTwoGroups m_SpringPotentialBetweenTwoGroups;
    PositionRescaleFrameTensionCoupling m_RescaleTenCoupl;
    Apply_Constant_Area m_Apply_Constant_Area;
    ActiveTwoStateInclusion m_ActiveTwoStateInclusion;
    LinkFlipMC m_LinkFlipMC;
    VertexMCMove m_VertexMoveMC;
    EdgeVertexMCMove m_EdgeVertexMoveMC;
    InclusionMCMove m_IncMove;
    Restart m_Restart;
    CoupleToWallPotential m_RigidWallCoupling;
    Curvature m_CurvatureCalculations;
    Energy m_EnergyCalculator;
    Constant_NematicForce *m_pConstant_NematicForce;
    DynamicBox            *m_pDynamicBox;
    DynamicTopology       *m_pDynamicTopology; 
    OpenEdgeEvolution     *m_pOpenEdgeEvolution;
    VolumeCoupling        *m_pVolumeCoupling;
private:
    Inclusion_Interaction_Map m_inc_ForceField;
    Constant_NematicForce m_Constant_NematicForce;
private:
    MESH          m_Mesh;
    void ExploreArguments();         // updates variables based on the command line arguments
    void HelpMessage();              // writes a help message
    void ReadInputFile(std::string);    // updates variables based on data in the inputfile
    void WriteStateLog();
   
    


};

#endif
