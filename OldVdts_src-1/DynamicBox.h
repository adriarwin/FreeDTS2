#if !defined(AFX_DynamicBox_H)
#define AFX_DynamicBox_H
#include <iostream>
#include "GenerateCNTCells.h"

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class DynamicBox {
public:
    DynamicBox();
    virtual ~DynamicBox();

    virtual bool MCMoveBoxChange(double dx, double* totalenergy, double temp, int step, GenerateCNTCells* pGenCNT) = 0;
    virtual void initialize() = 0;
    virtual  bool GetCNTCondition() = 0;
    virtual int GetTau() = 0;
    

};
//---- a class for no box change
class NoBoxChange : public DynamicBox {
public:
    NoBoxChange();
    ~NoBoxChange();

    
    bool GetCNTCondition();
    int GetTau();
    void initialize();
    bool MCMoveBoxChange(double dx, double * TotalEnergy, double temp, int step, GenerateCNTCells *pGenCNT );


private:

    int m_Tau;
    double m_F0;
    bool m_UpdateCNT;
};

#endif
