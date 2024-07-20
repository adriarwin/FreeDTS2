

#include <stdio.h>
#include "VolumeCoupling.h"
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
VolumeCoupling::VolumeCoupling() {
    // Constructor implementation

}

VolumeCoupling::~VolumeCoupling() {
    // Destructor implementation
}


//---- a class for no box change
NoCoupling::NoCoupling()
{
    m_TotalVolume = 0;
    m_TotalArea = 0;
    m_NoEQStep = 0;
    m_KV = 0;
    m_TargetV = 0;
    m_DeltaP = 0;
}
NoCoupling::~NoCoupling(){
    
}
void NoCoupling::Initialize(std::vector<triangle *> pTriangle){
    std::cout<<"---> note: there is no algorithm for volume coupling \n";
    m_TotalVolume = 0;
    m_TotalArea = 0;
    m_NoEQStep = 0;
    m_KV = 0;
    m_TargetV = 0;
    m_DeltaP = 0;
}






