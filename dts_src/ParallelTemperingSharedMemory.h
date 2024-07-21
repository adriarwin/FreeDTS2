#ifndef ParallelTemperingSharedMemory_FREEDTS_H_INCLUDED
#define ParallelTemperingSharedMemory_FREEDTS_H_INCLUDED

#include "AbstractParallelTemperingSharedMemory.h"
#include "SimDef.h"
#include <vector>
#include <string>

class ParallelReplicaData;
class ParallelTemperingSharedMemory : public AbstractParallelTemperingSharedMemory {
public:
    // Constructor
    ParallelTemperingSharedMemory(std::vector<std::string> Argument);
    
    // Destructor
    ~ParallelTemperingSharedMemory();
    inline std::string GetDerivedDefaultReadName()  {return "ParallelTemperingSharedMemory";}
    inline static std::string GetDefaultReadName()  {return "ParallelTemperingSharedMemory";}
    std::string CurrentState();
    
    bool Initialize(ParallelReplicaData PRD);
    bool Run();
    
private:
    // Private member variables and functions can be added here
    int m_Rate;
    int m_Bins;
    double m_minBeta;
    double m_maxBeta;
    std::vector<std::string> m_Argument;
    
};

#endif // ParallelTemperingSharedMemory_H_INCLUDED
