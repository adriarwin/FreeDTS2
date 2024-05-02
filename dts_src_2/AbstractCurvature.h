#if !defined(AFX_AbstractCurvature_H)
#define AFX_AbstractCurvature_H
#include <iostream>

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a Abstract class for curvature calculations.
========================================================
*/
class  AbstractCurvature {
public:
    AbstractCurvature(){
        
    }
    virtual ~ AbstractCurvature(){
        
    }
    virtual  bool SurfVertexCurvature(vertex *pvertex) = 0;
    virtual  bool EdgeVertexCurvature(vertex *pvertex) = 0;
    virtual  bool VertexCurvature(vertex *pvertex) = 0;
    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    
    inline static std::string GetBaseDefaultReadName() {return "CurvatureMethod";}
    
private:
    
};

#endif
