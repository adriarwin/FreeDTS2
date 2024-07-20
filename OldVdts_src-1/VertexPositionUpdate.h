#if !defined(AFX_VertexPositionUpdate_H)
#define AFX_VertexPositionUpdate_H
#include <iostream>

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for moving a vertex, this could include, MD, MD, SD,
========================================================
*/
class VertexPositionUpdate {
public:
    VertexPositionUpdate();
    virtual ~VertexPositionUpdate();

    virtual bool MCMoveBoxChange(double dx, double* totalenergy, double temp, int step, GenerateCNTCells* pGenCNT) = 0;
    virtual void initialize() = 0;
    virtual  bool GetCNTCondition() = 0;
    virtual int GetTau() = 0;
    

};


#endif
