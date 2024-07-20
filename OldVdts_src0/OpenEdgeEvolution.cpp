

#include <stdio.h>
#include "OpenEdgeEvolution.h"
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
OpenEdgeEvolution::OpenEdgeEvolution() {
    // Constructor implementation
}

OpenEdgeEvolution::~OpenEdgeEvolution() {
    // Destructor implementation
}


//---- a class for no box change
NoEvolution::NoEvolution()
{
}
NoEvolution::~NoEvolution(){
    
}
void NoEvolution::Initialize(){
    std::cout<<"---> note: there is no defined algorithm for edge treatment, edge size remains constant \n";
}
void NoEvolution::MC_Move(RNG* rng, double lmin, double lmax, double maxangle){

    return ;
}






