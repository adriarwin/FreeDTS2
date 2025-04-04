

#include <stdio.h>
#include "CoupleToWallPotential.h"
#include "Nfunction.h"
CoupleToWallPotential::CoupleToWallPotential()
{
    m_State = false;    
}
CoupleToWallPotential::CoupleToWallPotential(bool state, std::string info)
{
    m_ReachTargetWall = false;
    m_State = state;
    Nfunction f;
    m_EQTime = 0;
    std::vector<std::string> data = f.split(info);
    if(data.size()<2)
    {
        std::cout<<"----> Error: provided data for wall potential is not enough "<<std::endl;
        exit(0);
    }
    else
    {
        m_PotentialType = data[0];
        m_EQTime = f.String_to_Int(data[1]);
        for (int i=2;i<data.size();i++)
         m_Data.push_back(data.at(i));
    }
    m_H = 0;
    m_Lx = 0;
    m_Ly = 0;
    m_Lz = 0;
    m_A = 0;
    m_B = 0;
    m_C = 0;
    m_h = 0;
    
}
CoupleToWallPotential::~CoupleToWallPotential()
{
    
}
void CoupleToWallPotential::Initialize(std::vector <vertex *> Apv)
{

    Nfunction f;
    m_pAllV = Apv;
    // Obtain the center of geometry
    m_COG(0) = 0; m_COG(1) = 0; m_COG(2) = 0;
    for (std::vector<vertex*>::iterator it = Apv.begin() ; it != Apv.end(); ++it)
    {
        m_COG(0)= m_COG(0) + (*it)->GetVXPos()/double(Apv.size());
        m_COG(1)= m_COG(1) + (*it)->GetVYPos()/double(Apv.size());
        m_COG(2)= m_COG(2) + (*it)->GetVZPos()/double(Apv.size());
    }
    
    if(m_PotentialType=="TwoFlatParallelWall")
    {
        if(m_Data.size()== 0)
        {
            std::cout<<"----> Error: provided data for wall potential is not enough "<<std::endl;
            exit(0);
        }
        else
        {
        m_H = f.String_to_Double(m_Data[0]);
            if(m_H<=0)
            std::cout<<"----> Warning: The wall will never reaches the targeted value "<<std::endl;
        }
        // Finding m_h, we first see how much m_h can be, if it smaller then m_H, then we set it to m_H
        m_h = 0;
        for (std::vector<vertex*>::iterator it = Apv.begin() ; it != Apv.end(); ++it)
        {
            if(m_h<fabs((*it)->GetVZPos()-m_COG(2)))
            m_h=fabs((*it)->GetVZPos()-m_COG(2));
        }
        if(m_h<m_H)
        m_h = m_H;
    }
    else if(m_PotentialType=="Cuboid")
    {
        if(m_Data.size()< 3)
        {
            std::cout<<"----> Error: provided data for wall potential is not enough "<<std::endl;
            exit(0);
        }
        else
        m_Lx = f.String_to_Double(m_Data[0]);
        m_Ly = f.String_to_Double(m_Data[1]);
        m_Lz = f.String_to_Double(m_Data[2]);
        
        // Finding m_lx, m_ly, m_lz, we first see how much m_lx ... can be, if it smaller then m_Lx, then we set it to m_Lx
        m_lx = 0; m_ly = 0; m_lz = 0;
        for (std::vector<vertex*>::iterator it = Apv.begin() ; it != Apv.end(); ++it)
        {
            double dx = fabs((*it)->GetVXPos()-m_COG(0));
            double dy = fabs((*it)->GetVYPos()-m_COG(1));
            double dz = fabs((*it)->GetVZPos()-m_COG(2));

            if(m_lx<dx)
            m_lx=dx;
            if(m_ly<dy)
            m_ly=dy;
            if(m_lz<dz)
            m_lz=dz;
        }
        if(m_lx<m_Lx)
        m_lx = m_Lx;
        if(m_ly<m_Ly)
        m_ly = m_Ly;
        if(m_lz<m_Lz)
        m_lz = m_Lz;

    }
    else if(m_PotentialType=="Ellipsoid")
    {
        m_A = f.String_to_Double(m_Data[0]);
        m_B = f.String_to_Double(m_Data[1]);
        m_C = f.String_to_Double(m_Data[2]);
        m_r = 0;
        for (std::vector<vertex*>::iterator it = Apv.begin() ; it != Apv.end(); ++it)
        {
            double dx = (*it)->GetVXPos()-m_COG(0);
            double dy = (*it)->GetVYPos()-m_COG(1);
            double dz = (*it)->GetVZPos()-m_COG(2);
            
            double r2 = dx*dx+m_A*m_A*dy*dy/(m_B*m_B)+m_A*m_A*dz*dz/(m_C*m_C);
                if(r2>m_r*m_r)
                 m_r = sqrt(r2);
        }
    }
    else if(m_PotentialType=="EllipsoidalShell")
    {
        m_A = f.String_to_Double(m_Data[0]);
        m_B = f.String_to_Double(m_Data[1]);
        m_C = f.String_to_Double(m_Data[2]);
        m_H = f.String_to_Double(m_Data[3]);
        m_B = m_B/m_A; m_C = m_C/m_A;
        
        double rmax = 0; double rmin = 10000000000;
        for (std::vector<vertex*>::iterator it = m_pAllV.begin() ; it != m_pAllV.end(); ++it)
        {
            double dx = (*it)->GetVXPos()-m_COG(0);
            double dy = (*it)->GetVYPos()-m_COG(1);
            double dz = (*it)->GetVZPos()-m_COG(2);
            
            double r2 = dx*dx+dy*dy/(m_B*m_B)+dz*dz/(m_C*m_C);
            if(r2>rmax*rmax)
            rmax = sqrt(r2);
            if(r2<rmin*rmin)
            rmin = sqrt(r2);
        }
        m_h = rmax - rmin;
        m_r = (rmax + rmin)/2;
    }
    else
    {
        std::cout<<"----> Error:  wall potential type "<<m_PotentialType<<" has not been defined yet. "<<std::endl;
        exit(0);
    }
    
}
bool CoupleToWallPotential::CheckVertexMoveWithinWalls(int step, double dx, double dy, double dz, vertex* v)
{
    bool accept=true;
    
    // check if the step is larger then m_EQTime
    if(m_ReachTargetWall==false)
    {
        if(step<m_EQTime)
        {
            if(AllVerticesAreInsideTheBound()==true)
            MoveTheWallsTowardTheTarget(step);
        }
        else
        {
            m_ReachTargetWall = true; // tem solution
            std::cout<<"---> Warning, the simulation step is larger equilibration step but we have not reach the targeted wall "<<std::endl;
            std::cout<<" better to increase the equilibration steps "<<std::endl;
        }
    }
    Vec3D X(v->GetVXPos()+dx,v->GetVYPos()+dy,v->GetVZPos()+dz);
    if(APointIsInsideTheBound(X)==false)
        accept = false;
    
    return accept;
}
void CoupleToWallPotential::MoveTheWallsTowardTheTarget(int step)
{
    double DR = 0.02;
    if(m_PotentialType=="TwoFlatParallelWall")
    {
        double dx = fabs(m_h-m_H);
        if(dx<DR)
        {
            m_h =m_H;
            m_ReachTargetWall = true;
        }
        else
        {
            m_h = m_h-DR;
        }
    }
    else if(m_PotentialType=="Cuboid")
    {
        double dx = fabs(m_lx-m_Lx);
        double dy = fabs(m_ly-m_Ly);
        double dz = fabs(m_lz-m_Lz);

        if(dx<DR)
        {
            m_lx =m_Lx;
        }
        else
        {
            m_lx = m_lx-DR;
        }
        if(dy<DR)
        {
            m_ly =m_Ly;
        }
        else
        {
            m_ly = m_ly-DR;
        }
        if(dz<DR)
        {
            m_lz =m_Lz;
        }
        else
        {
            m_lz = m_lz-DR;
        }
        if(m_lz == m_Lz && m_lx == m_Lx && m_ly == m_Ly)
        {
            m_ReachTargetWall = true;
        }

    }
    else if(m_PotentialType=="Ellipsoid")
    {
        double dx = fabs(m_r-m_A);
        if(dx<DR)
        {
            m_r =m_A;
            m_ReachTargetWall = true;
        }
        else
        {
            m_r = m_r-DR;
        }
    }
    else if(m_PotentialType=="EllipsoidalShell")
    {
        double dx = fabs(m_h-m_H);
        if(dx<DR)
        {
            m_h =m_H;
            m_ReachTargetWall = true;
        }
        else
        {
            m_h = m_h-DR;
        }
        
        double rmax = 0; double rmin = 10000000000;
        for (std::vector<vertex*>::iterator it = m_pAllV.begin() ; it != m_pAllV.end(); ++it)
        {
            double dx = (*it)->GetVXPos()-m_COG(0);
            double dy = (*it)->GetVYPos()-m_COG(1);
            double dz = (*it)->GetVZPos()-m_COG(2);
            
            double r2 = dx*dx+dy*dy/(m_B*m_B)+dz*dz/(m_C*m_C);
            if(r2>rmax*rmax)
                rmax = sqrt(r2);
            if(r2<rmin*rmin)
                rmin = sqrt(r2);
        }
        m_r = (rmax + rmin)/2;
    }
        
}
bool CoupleToWallPotential::AllVerticesAreInsideTheBound()
{
    bool Itis = true;
    // When we make the wall more confined, it is better to check if all the vertices moved inside the bound before making the wall more confinded again
    for (std::vector<vertex*>::iterator it = m_pAllV.begin() ; it != m_pAllV.end(); ++it)
    {
        Vec3D X((*it)->GetVXPos(),(*it)->GetVYPos(),(*it)->GetVZPos());
        if(APointIsInsideTheBound(X)==false)
        {
            Itis = false;
            break;
        }
    }
    
    return Itis;
}
bool CoupleToWallPotential::APointIsInsideTheBound(Vec3D X)
{
    // checking if a point is inside the wall
    bool Itis = true;
    if(m_PotentialType=="TwoFlatParallelWall")
    {
        double dx = fabs(X(2)-m_COG(2));
        if(dx>m_h)
        {
        Itis = false;
        }
    }
    else if(m_PotentialType=="Cuboid")
    {
        double dx = fabs(X(0)-m_COG(0));
        double dy = fabs(X(1)-m_COG(1));
        double dz = fabs(X(2)-m_COG(2));

        if(dx>m_lx)
        Itis = false;
        if(dy>m_ly)
        Itis = false;
        if(dz>m_lz)
        Itis = false;
    }
    else if(m_PotentialType=="Ellipsoid")
    {
        double dx = fabs(X(0)-m_COG(0));
        double dy = fabs(X(1)-m_COG(1));
        double dz = fabs(X(2)-m_COG(2));
        
        double r2 = dx*dx+m_A*m_A*dy*dy/(m_B*m_B)+m_A*m_A*dz*dz/(m_C*m_C);
        if(r2>m_r*m_r)
        Itis = false;
    }
    else if(m_PotentialType=="EllipsoidalShell")
    {
        double dx = X(0)-m_COG(0);
        double dy = X(1)-m_COG(1);
        double dz = X(2)-m_COG(2);
        double r2 = dx*dx+dy*dy/(m_B*m_B)+dz*dz/(m_C*m_C);
        
        if(r2>(m_r+m_h/2)*(m_r+m_h/2) || r2<(m_r-m_h/2)*(m_r-m_h/2) )
            Itis = false;
    }
    
    return Itis;
}




