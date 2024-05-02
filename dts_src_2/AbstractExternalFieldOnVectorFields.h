#if !defined(AFX_AbstractExternalFieldOnVectorFields_H)
#define AFX_AbstractExternalFieldOnVectorFields_H
#include <iostream>
#include "vertex.h"

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a Abstract class for curvature calculations.
========================================================
*/
class State;
class  AbstractExternalFieldOnVectorFields {
public:
    AbstractExternalFieldOnVectorFields(){
        
    }
    virtual ~ AbstractExternalFieldOnVectorFields(){
        
    }
    virtual double GetCouplingEnergy(vertex *pvertex) = 0;
    
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "ExternalFieldOnVectorFields";}

    
private:
    
};
//---- a class for no box change
class NoExternalField : public AbstractExternalFieldOnVectorFields {
public:
    NoExternalField(){
        
    }
    ~NoExternalField(){
        
    }
    inline std::string GetDerivedDefaultReadName()  {return "NoExternalField";};
    double GetCouplingEnergy(vertex *pvertex){
        return 0;
    }
};
#endif
