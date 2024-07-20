#if !defined(AFX_GlobalCurvature_H)
#define AFX_GlobalCurvature_H
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
class  GlobalCurvature {
public:
    GlobalCurvature(){
        
    }
    virtual ~ GlobalCurvature(){
        
    }
    virtual  double GetEnergy() = 0;
    virtual  bool GetState() = 0;
    virtual  void Initialize(std::vector<vertex *> &Ver) = 0;
    virtual  void UpdateEnergyChange(double delta_area, double delta_curvature) = 0;
    virtual  double CalculateEnergyChange(double delta_area, double delta_curvature) = 0;
    

};
//---- a class for no box change
class NoGlobalCurvature : public  GlobalCurvature {
public:
    NoGlobalCurvature(){
        
    }
    ~NoGlobalCurvature(){
        
    }
    
    inline bool GetState()                           {return false;} // if the coupling is active
    inline double GetEnergy()                        {return 0;} // Only part of energy asscoiated with this class
    
    void Initialize(std::vector<vertex *> &Ver){
     
        return;
    }
    void UpdateEnergyChange(double delta_area, double delta_curvature){
        return;
    }
    double CalculateEnergyChange(double delta_area, double delta_curvature){
        return 0;
    }
    
};

#endif
