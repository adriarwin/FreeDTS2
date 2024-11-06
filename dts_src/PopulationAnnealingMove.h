#if !defined(AFX_PopulationAnnealingMove_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_)
#define AFX_PopulationAnnealingMove_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "HarmonicPotentialBetweenTwoGroups.h"
#include "triangle.h"
#include "MESH.h"
#include "AbstractPopulationAnnealingMove.h"
#include <fstream>
#ifdef MPI_DETECTED
#include <mpi.h>
#endif


class State;
class PopulationAnnealingMove : public AbstractPopulationAnnealingMove {
public:
    
    PopulationAnnealingMove(State *pState, std::string period_file, std::string temperature_file, std::string input_files, int n_processors);
	~PopulationAnnealingMove();

    void Initialize();
    bool EvolveOneStep(int step);

    inline  std::string GetDerivedDefaultReadName()  {return "PopulationAnnealingMoveSimple";}
    inline  static std::string GetBaseDefaultReadName()  {return "PopulationAnnealingMoveSimple";}

    
    std::string CurrentState();

    bool PopulationAnnealingMoveOn();
    void SetRestart();


private:
    std::vector<int> ReorderNewRanks(std::vector<int> new_ranks);   
    bool ChangeToNewTemperatureID(int NewTempID);
    bool isInteger(const std::string& s);
    std::vector<int> splitStringToIntVector(const std::string& line);
    std::vector<int> ReadLastLineOutput(const std::string& output, int N);
    std::string GetTopologyFile(){return m_pTopologyFileString;};
//Here I can insert things no one will ever know...
    
    
private:

    std::vector<std::string> m_pRestartFilesVector;
    std::string m_pPeriodFile; //File with all periods
    std::string m_pTemperatureFile; //File with all temperatures
    std::string m_pTopologyFile; //File with path to all input files
    int m_pInputSize; //Number of processors that you input
    int m_pSize; //Number of processors detected

    std::vector<int> m_pNewRank;

    std::vector<double> m_pBetaVector;
    std::vector<int> m_pPeriodVector;
    std::vector<std::string> m_pTopologyFilesVector;
    std::string m_pTopologyFileString;
    int m_pRank; //What is the rank of current processor.
    int m_pPeriod;
    int m_pCounter;
    int m_pNextConfiguration;
    State *m_pState;
    
    inline  std::string GetInitialTemperaturesFileName() {return "initial_temperatures.txt";}
    inline  std::string GetOutputFileName() {return "output_descending_mapping.txt";}
    std::vector<double> ReadTemperatures();
    std::vector<int> ReadPeriods();
    bool is_line_empty(const std::string& line);
    void WriteRankAtTempIDToFile();

    std::vector<double> m_pEnergyVector;
    std::ofstream m_TimeSeriesFile;
    void WriteParentalDescendingMappingToFile();
    bool m_pRestart=false;    

    

};


#endif