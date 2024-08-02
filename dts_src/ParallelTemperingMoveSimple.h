#if !defined(AFX_ParallelTemperingMoveSimple_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_)
#define AFX_ParallelTemperingMoveSimple_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "HarmonicPotentialBetweenTwoGroups.h"
#include "triangle.h"
#include "MESH.h"
#include "AbstractParallelTemperingMove.h"
#ifdef MPI_DETECTED
#include <mpi.h>
#endif


class State;
class ParallelTemperingMoveSimple : public AbstractParallelTemperingMove {
public:
    
    ParallelTemperingMoveSimple(State *pState, int period, int nprocessors, double minbeta, double maxbeta);
	~ParallelTemperingMoveSimple();

    void Initialize();
    bool EvolveOneStep(int step);

    inline  std::string GetDerivedDefaultReadName()  {return "ParallelTemperingMoveSimple";}
    inline  static std::string GetBaseDefaultReadName()  {return "ParallelTemperingMoveSimple";}
    std::string CurrentState();

private:
    bool ChangeToNewTemperatureID(int NewTempID);
//Here I can insert things no one will ever know...
    
    
private:
    int m_Nprocessors;
    double m_MinBeta;
    double m_MaxBeta;
    int m_Period;
    int m_Rank;
    int m_Size;
    int m_TempID;
    int m_Counter;
    bool m_SizeIsEven;
    std::vector<int>  m_RankAtTempID;   // a map for temp_id to thread_id, temp = temperatures.at(temp_id)
    std::vector<double> m_BetaVec;        // the 1/temprature of each temp_id
    State *m_pState;
    
    std::vector<int> m_ReceiveBroadcast;
    #ifdef MPI_DETECTED
    std::vector<MPI_Request> m_RequestBroadcast;
    #endif
    


};


#endif