

#include <stdio.h>
#include "DynamicBox.h"
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
DynamicBox::DynamicBox() {
    // Constructor implementation

}

DynamicBox::~DynamicBox() {
    // Destructor implementation
}


//---- a class for no box change
NoBoxChange::NoBoxChange()
{
    m_F0 = 0;
    m_Tau = 0;
}
NoBoxChange::~NoBoxChange(){
    
}
void NoBoxChange::initialize(){
    std::cout<<"---> note: there is no algorithm for box change, the box size will remain constant \n";
    m_Tau = 0;
    m_F0 = 0;
}
bool NoBoxChange::GetCNTCondition(){
    return true;
}
int NoBoxChange::GetTau() {
    return m_Tau;
}
bool NoBoxChange::MCMoveBoxChange(double dx, double * tot_Energy, double temp, int step, GenerateCNTCells *pGenCNT) {

    return false;
}






