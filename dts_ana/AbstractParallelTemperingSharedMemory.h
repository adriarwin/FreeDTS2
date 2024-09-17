#ifndef AbstractParallelTemperingSharedMemory_FREEDTS_H_INCLUDED
#define AbstractParallelTemperingSharedMemory_FREEDTS_H_INCLUDED

#include "SimDef.h"
#include <vector>
#include <string>

class ParallelReplicaData;
class AbstractParallelTemperingSharedMemory {
public:
    // Constructor
    AbstractParallelTemperingSharedMemory(){
        
    }
    // Destructor
    virtual ~AbstractParallelTemperingSharedMemory(){
        
    }
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "ParallelTemperingSharedMemory";}
    
    virtual bool Initialize(ParallelReplicaData PRD) = 0;
    virtual bool Run() = 0;
    
private:

};


#endif // AbstractParallelReplicaRun_H_INCLUDED
