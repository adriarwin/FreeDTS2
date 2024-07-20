#if !defined(AFX_VolumeCoupling_H)
#define AFX_VolumeCoupling_H
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
class VolumeCoupling {
public:
    VolumeCoupling();
    virtual ~VolumeCoupling();

    virtual inline bool GetState()= 0;
    virtual inline double GetTotalVolume() = 0;
    virtual inline double GetTotalArea() = 0;
    virtual void Initialize(std::vector<triangle *> pTriangle) = 0;
    virtual double SingleTriangleVolume(triangle * ptriangle) = 0;
    virtual double GetEnergyChange(int s, double oa, double ov, double na, double nv) = 0;
    virtual void UpdateArea_Volume(double oa, double ov, double na, double nv) = 0;
    

};
//---- a class for no box change
class NoCoupling : public VolumeCoupling {
public:
    NoCoupling();
    ~NoCoupling();
    inline bool GetState()                   {return false;}
    inline double GetTotalVolume()                  {return 0;}
    inline double GetTotalArea()                  {return 0;}
    
    void    Initialize(std::vector<triangle *> pTriangle);
    double  VolumeofTrianglesAroundVertex(vertex * pVeretx)  {return 0;}
    double  GetEnergyChange(int s, double oa, double ov, double na, double nv) {return 0;}
    double  Energy(double volume, double area, double a) {return 0;}
    void    UpdateArea_Volume(double oa, double ov, double na, double nv) {return ;}
    double  SingleTriangleVolume(triangle * ptriangle)  {return 0;}


private:

    double m_TotalVolume;
    double m_TotalArea;
    int m_NoEQStep;
    double m_KV;
    double m_TargetV;
    double m_DeltaP;
};

#endif
