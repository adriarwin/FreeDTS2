#if !defined(AFX_AnnealedImportanceSamplingMove_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_)
#define AFX_AnnealedImportanceSamplingMove_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "HarmonicPotentialBetweenTwoGroups.h"
#include "triangle.h"
#include "MESH.h"
#include "AbstractAnnealedImportanceSamplingMove.h"
#include <fstream>
#ifdef MPI_DETECTED
#include <mpi.h>
#endif


class State;
class AnnealedImportanceSamplingMove : public AbstractAnnealedImportanceSamplingMove {
public:
    
    AnnealedImportanceSamplingMove(State *pState, std::string period_file, std::string temperature_file);
	~AnnealedImportanceSamplingMove();

    void Initialize();
    bool EvolveOneStep(int step);

    inline  std::string GetDerivedDefaultReadName()  {return "AnnealedImportanceSampling";}
    inline  static std::string GetBaseDefaultReadName()  {return "AnnealedImportanceSampling";}

    
    std::string CurrentState();

    bool AnnealedImportanceSamplingMoveOn();
    void SetRestart();


private:  
    bool ChangeToNewTemperatureID(int NewTempID);
    bool isInteger(const std::string& s);
    std::vector<int> splitStringToIntVector(const std::string& line);
    void WriteWeight(double weight);
//Here I can insert things no one will ever know...
    
    
private:

    std::string m_pPeriodFile; //File with all periods
    std::string m_pTemperatureFile; //File with all temperatures



    std::vector<double> m_pBetaVector;
    std::vector<int> m_pPeriodVector;

    int m_pPeriod;
    int m_pCounter;
    State *m_pState;
    
    inline  std::string GetInitialTemperaturesFileName() {return "initial_temperatures.txt";}
    inline  std::string GetOutputFileName() {return "weights.txt";}
    std::vector<double> ReadTemperatures();
    std::vector<int> ReadPeriods();
    bool is_line_empty(const std::string& line);

    std::ofstream m_TimeSeriesFile;
    void WriteParentalDescendingMappingToFile();
    bool m_pRestart=false;    

    

};


#endif