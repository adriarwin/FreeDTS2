#if !defined(AFX_Boundary_H)
#define AFX_Boundary_H
#include <iostream>

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for rigid wall boundry
========================================================
*/
class  Boundary {
public:
    Boundary(){
        
    }
    virtual ~ Boundary(){
        
    }
    virtual void Initialize(std::vector <vertex *> &Apv);
    virtual bool MoveHappensWithinTheBoundary(int step, double x, double y, double z, vertex* v) = 0;

    

};
//---- a class for no box change
class PBCBoundary : public  Boundary {
public:
    PBCBoundary(){
        
    }
    ~PBCBoundary(){
        
    }

    
    void Initialize(std::vector <vertex *> &Apv){
     
        return;
    }
    bool MoveHappensWithinTheBoundary(int step, double x, double y, double z, vertex* v){
        return true;
    }
    
};

#endif
