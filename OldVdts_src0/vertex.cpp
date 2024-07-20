 #if !defined(AFX_vertex_CPP_7F4A21C7_C13C_1223_BF2E_124095086234__INCLUDED_)
#define AFX_vertex_CPP_7F4A21C7_C13C_1223_BF2E_124095086234__INCLUDED_

#include <stdio.h>
#include "vertex.h"
#include "links.h"
#include "triangle.h"
#include "Nfunction.h"
vertex::vertex()
{
}
vertex::vertex(int id, double x, double y, double z)
{

m_X = x;
m_Y = y;
m_Z = z;
m_ID = id;
    m_SimTimeStep=-1;
    m_Group = 0;
    m_OwnInclusion = false;
    m_GroupName = "system";
    
    
    m_Geodesic_Curvature = 0;
    m_Normal_Curvature = 0;
    m_VertexType = 0;
    m_VLength = 0;                       // length of the vertex
    m_Lambda = 0;                   // line tension
    
    
}
vertex::vertex(int id)
{

m_ID=id;
m_X=0;
m_Y=0;
m_Z=0;
    m_SimTimeStep=-1;
    m_Group = 0;
    m_OwnInclusion = false;
    m_GroupName = "system";
    
    
    m_Geodesic_Curvature = 0;
    m_Normal_Curvature = 0;
    m_VertexType = 0;
    m_VLength = 0;                       // length of the vertex
    m_Lambda = 0;                   // line tension

}
vertex::~vertex()
{
    
}
void vertex::UpdateGroupName(std::string x)
{
    m_GroupName=x;
}
void vertex::UpdateBox(Vec3D *x)
{
m_pBox=x;
}
void vertex::UpdateVID(int i)
{
    m_ID=i;
}
void vertex::UpdateVXPos(double x) {
    double boxX = (*m_pBox)(0);

    // Adjust X position if it's outside the box bounds
    if (x >= boxX) {
        m_X = x - boxX;
    } else if (x < 0) {
        m_X = x + boxX;
    } else {
        m_X = x;
    }
}
void vertex::UpdateVYPos(double y) {
    double boxY = (*m_pBox)(1);

    // Adjust Y position if it's outside the box bounds
    if (y >= boxY) {
        m_Y = y - boxY;
    } else if (y < 0) {
        m_Y = y + boxY;
    } else {
        m_Y = y;
    }
}
void vertex::UpdateVZPos(double z) {
    double boxZ = (*m_pBox)(2);

    // Adjust Z position if it's outside the box bounds
    if (z >= boxZ) {
        m_Z = z - boxZ;
    } else if (z < 0) {
        m_Z = z + boxZ;
    } else {
        m_Z = z;
    }
}
void vertex::AddtoLinkList(links* z)
{
m_VLinkList.push_back(z);
}
void vertex::RemoveFromLinkList(links* z)
{

m_VLinkList.erase(std::remove(m_VLinkList.begin(), m_VLinkList.end(), z), m_VLinkList.end());
}
void vertex::AddtoTraingleList(triangle* z)
{
m_VTraingleList.push_back(z);
}
void vertex::RemoveFromTraingleList(triangle* z)
{
m_VTraingleList.erase(std::remove(m_VTraingleList.begin(), m_VTraingleList.end(), z), m_VTraingleList.end());
}
void vertex::AddtoNeighbourVertex(vertex* z)
{
m_VNeighbourVertex.push_back(z);
}
void vertex::RemoveFromNeighbourVertex(vertex* z)
{
m_VNeighbourVertex.erase(std::remove(m_VNeighbourVertex.begin(), m_VNeighbourVertex.end(), z), m_VNeighbourVertex.end());
}
void vertex::UpdateInclusion(inclusion * inc)
{
m_pInclusion=inc;
}
void vertex::UpdateOwnInclusion(bool own)
{
m_OwnInclusion=own;
}
void vertex::UpdateVCNTCell(CNTCell * z)
{
m_CNTCell=z;
}
void vertex::UpdateGroup(int z)
{
    m_Group=z;
}
void vertex::UpdateEnergy(double z)
{
m_Energy=z;
}
//////////// Functions related to Curvature energy and normal
void vertex::UpdateCurvature(double x,double y)
{
    m_Curvature.clear();
    m_Curvature.push_back(x);
    m_Curvature.push_back(y);
}
void vertex::UpdateNormal_Area(Vec3D v,double a)
{
    m_Normal=v;
    m_Area=a;
}
void vertex::UpdateL2GTransferMatrix(Tensor2 v)
{
    m_T_Local_2_Global = v;
}
void vertex::UpdateG2LTransferMatrix(Tensor2 v)
{
    m_T_Global_2_Local = v;
}
void vertex::UpdateSimTimeStep(int v)
{
    m_SimTimeStep=v;
}
bool vertex::CheckCNT() {
    Vec3D min = m_CNTCell->GetCNTSidemin();
    Vec3D max = m_CNTCell->GetCNTSidemax();

    // Check X coordinate
    if (m_X < min(0) || m_X >= max(0)) {
        std::cout << m_X << "  x " << min(0) << "  " << max(0) << "\n";
        return false;
    }

    // Check Y coordinate
    if (m_Y < min(1) || m_Y >= max(1)) {
        std::cout << m_Y << "  y " << min(1) << "  " << max(1) << "\n";
        return false;
    }

    // Check Z coordinate
    if (m_Z < min(2) || m_Z >= max(2)) {
        std::cout << m_Z << "  z " << min(2) << "  " << max(2) << "\n";
        return false;
    }

    // All coordinates are within bounds
    return true;
}
#endif



