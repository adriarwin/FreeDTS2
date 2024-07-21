#ifndef AbstractParallelTemperingDistributedMemory_FREEDTS_H_INCLUDED
#define AbstractParallelTemperingDistributedMemory_FREEDTS_H_INCLUDED

#include "SimDef.h"
#include <vector>
#include <string>

class ParallelReplicaData;
class AbstractParallelTemperingDistributedMemory {
public:
    // Constructor
    AbstractParallelTemperingDistributedMemory(){
        
    }
    // Destructor
    virtual ~AbstractParallelTemperingDistributedMemory(){
        
    }
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "ParallelTemperingDistributedMemory";}
    
    virtual bool Initialize(ParallelReplicaData PRD) = 0;
    virtual bool Run() = 0;
    
private:

};


#endif // AbstractParallelTemperingDistributedMemory_H_INCLUDED
