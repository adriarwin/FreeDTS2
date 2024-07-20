#if !defined(AFX_TotalAreaCoupling_H)
#define AFX_TotalAreaCoupling_H
#include <iostream>
#include "GenerateCNTCells.h"

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for global Curvature and how energy should change
========================================================
*/
class  TotalAreaCoupling {
public:
    TotalAreaCoupling(){
        
    }
    virtual ~ TotalAreaCoupling(){
        
    }
    virtual  double GetTotalArea() = 0;
    virtual  bool GetState() = 0;
    virtual  void Initialize(std::vector<triangle *> &pTriangle) = 0;
    virtual  void UpdateArea(double oldarea, double newarea) = 0;
    virtual  double CalculateEnergyChange(int step, double oldarea,  double newarea) = 0;
    

};
//---- a class for no box change
class NoTotalAreaCoupling : public  TotalAreaCoupling {
public:
    NoTotalAreaCoupling(){
        
    }
    ~NoTotalAreaCoupling(){
        
    }
    
    inline bool GetState()                           {return false;} // if the coupling is active
    inline double GetTotalArea()                        {return 0;} // Only part of energy asscoiated with this class
    
    void Initialize(std::vector<triangle *> &pTriangle){
     
        return;
    }
    void UpdateArea(double oldarea, double newarea){
        return;
    }
    double CalculateEnergyChange(int step, double oldarea,  double newarea){
        return 0;
    }
    
};

#endif
