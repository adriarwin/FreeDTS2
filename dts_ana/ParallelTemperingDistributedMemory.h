#ifndef ParallelTemperingDistributedMemory_FREEDTS_H_INCLUDED
#define ParallelTemperingDistributedMemory_FREEDTS_H_INCLUDED

#include "AbstractParallelTemperingDistributedMemory.h"
#include "SimDef.h"
#include <vector>
#include <string>

class ParallelReplicaData;
class ParallelTemperingDistributedMemory : public AbstractParallelTemperingDistributedMemory {
public:
    // Constructor
    ParallelTemperingDistributedMemory(std::vector<std::string> Argument);
    
    // Destructor
    ~ParallelTemperingDistributedMemory();
    inline std::string GetDerivedDefaultReadName()  {return "ParallelTemperingDistributedMemory";}
    inline static std::string GetDefaultReadName()  {return "ParallelTemperingDistributedMemory";}
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

#endif // ParallelTemperingDistributedMemory_H_INCLUDED
