

#include <stdio.h>
#include "DynamicTopology.h"
#include "Nfunction.h"
#include "vertex.h"
#include "Curvature.h"
#include "State.h"

/*
===============================================================================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
=================================================================================================================
*/
DynamicTopology::DynamicTopology() {
    // Constructor implementation

}

DynamicTopology::~DynamicTopology() {
    // Destructor implementation
}


//---- a class for no box change
ConstantTopology::ConstantTopology()
{

}
ConstantTopology::~ConstantTopology(){
    
}
void ConstantTopology::initialize(){
    std::cout<<"---> note: there is no algorithm for topology change, the surface topology will remain constant \n";
}
bool ConstantTopology::MCMove(double * tot_Energy, RNG *rng, GenerateCNTCells *pGenCNT) {

    return false;
}






