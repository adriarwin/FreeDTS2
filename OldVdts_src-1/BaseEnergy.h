#if !defined(AFX_BaseEnergy_H)
#define AFX_BaseEnergy_H
#include <iostream>

// Define a base class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for system energy. Only includes the individual vertex and inclsuion energy
========================================================
*/
class  BaseEnergy {
public:
    BaseEnergy(){
        
    }
    virtual ~ BaseEnergy(){
        
    }
    virtual  double CalaculateEnergyofSingleVertex(vertex * pvertex) = 0;
    
    
//---> real functions
    inline double GetEnergy()               const                 {return m_TotalEnergy;}
    double UpdateTotalEnergy(double en){
       
        m_TotalEnergy = en;
        return;
    }
    
private:
    double m_TotalEnergy;
    
};

#endif
