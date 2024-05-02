

#include "links.h"

links::links(int id, vertex *v1, vertex *v2, triangle *t1)
{
    m_IntEnergy = 0;
    m_T1=t1;
    m_V1=v1;
    m_V2=v2;
    m_SimTimeStep=-1;
    m_ID=id;
    m_LinkSide = 0;
    m_mirorflag=false;
    m_Show=true;
    m_EdgeSize = 0;
    m_LinkType = 0;
}
links::links(int id)
{
    m_IntEnergy = 0;
    m_ID=id;
    m_LinkSide = 0;
    m_mirorflag=false;
    m_Show=true;
    m_SimTimeStep=-1;
    m_EdgeSize = 0;
    m_LinkType = 0;
}

links::~links()
{
    
}
void links::UpdateTriangle(triangle *v){
    m_T1=v;
    return;
}
void links::UpdateV(vertex *v1,vertex *v2,vertex *v3){
    m_V3=v3;
    m_V2=v2;
    m_V1=v1;
    return;
}
void links::UpdateSimTimeStep(int v){
    m_SimTimeStep=v;
    return;
}
void links::UpdateV3(vertex *v3){
    m_V3=v3;
    return;
}
void links::UpdateMirrorLink(links* v){
    m_mirorlink=v;
    return;
}
bool links::SetCopy(){           // Copies the key ellements into the old type
    m_OldT1 = m_T1;     //
    m_OldV1 = m_V1;
    m_OldV2 = m_V2;
    m_OldV3 = m_V3;
    m_Oldmirorlink = m_mirorlink;
    m_Oldneighborlink1 = m_neighborlink1;
    m_Oldneighborlink2  = m_neighborlink2;
    m_Oldmirorflag = m_mirorflag;
    m_OldLinkSide = m_LinkSide;
    m_OldNormal = m_Normal;
    m_OldBe = m_Be;
    m_OldHe = m_He;
    m_OldIntEnergy = m_IntEnergy;
    m_OldEdgeVector = m_EdgeVector;
    m_OldEdgeSize = m_EdgeSize;
    m_OldLinkType = m_LinkType;

    return true;
}
bool links::Reverse2PreviousCopy(){           // reverse to the last copy and part of the mirror variables
    
    m_T1 = m_OldT1;
    m_V1 = m_OldV1;
    m_V2 = m_OldV2;
    m_V3 = m_OldV3;
    m_mirorlink = m_Oldmirorlink;
    m_neighborlink1 = m_Oldneighborlink1;
    m_neighborlink2 = m_Oldneighborlink2;
    m_mirorflag = m_Oldmirorflag;
    m_LinkSide = m_OldLinkSide;
    m_Normal = m_OldNormal;
    m_Be = m_OldBe;
    m_He = m_OldHe;
    m_IntEnergy = m_OldIntEnergy;
    m_EdgeVector = m_OldEdgeVector;
    m_EdgeSize = m_OldEdgeSize;
    m_LinkType = m_OldLinkType;
    
    //---- update some for mirror
    m_mirorlink->PutShapeOperator(m_Be,m_He);
    m_mirorlink->UpdateIntEnergy(m_IntEnergy);
    m_mirorlink->PutNormal(m_Normal);

    return true;
}
void links::UpdateNeighborLink1(links* v){
    m_neighborlink1=v;
}
void links::UpdateNeighborLink2(links* v){
    m_neighborlink2=v;
}
void links::UpdateVisualize(bool v){
    m_Show=v;
}
void links::UpdateMirrorFlag(bool v){
    m_mirorflag=v;
    return;
}
void links::UpdateNormal()
{
   if(this->GetMirrorFlag()==true){
       Vec3D v2=(m_mirorlink->GetTriangle())->GetNormalVector();
       Vec3D v1=m_T1->GetNormalVector();
       m_Normal=v1+v2;
       double norm=m_Normal.norm();
       if(norm==0){
           std::cout<<"error 2022----> one of the normals has zero size; normal link cannot be defined  \n";
           exit(0);
       }
       m_Normal=m_Normal*(1.0/norm);
       m_mirorlink->PutNormal(m_Normal);
    }
    else{
        // this is an edge link
        std::cout<<" developer error: link type and id "<<m_LinkType<<"  "<<m_ID<<" \n";
        std::cout<<"error ----> normal vector for edge links has not been defined   \n";
        exit(0);
    }
    return;
}
void links::PutNormal(Vec3D n)
{
    m_Normal=n;
    
}
void links::PutShapeOperator(Vec3D Be,double He)
{
    m_Be=Be;
    m_He=He;
    
}
void links::UpdateIntEnergy(double en)
{
    m_IntEnergy = en;
}
void links::UpdateShapeOperator(Vec3D *pBox)
{
    UpdateEdgeVector(pBox);
   if(this->GetMirrorFlag()==true)
   {
       Vec3D Re = m_EdgeVector;
       Re=Re*(1.0/m_EdgeSize);
       Vec3D Be=m_Normal*Re;

    //====== Finding the size of the Be vector to make it normaized; this should not be needed
    // just have it so in case.
       double size=Be.norm();
       if(size!=0)
       {
           size=1.0/size;
       }
       else
       {
           std::cout<<" error 7634---> this should not happen \n";
           exit(0);
       }
       Be=Be*size;
       Vec3D Nf1=(m_mirorlink->GetTriangle())->GetNormalVector();
       Vec3D Nf2=m_T1->GetNormalVector();

       
//=== this is different from the orginal paper; it is faster
//==========
       double sign=Re.dot(Nf1*Nf2,Re);
       double tangle=Re.dot(Nf1,Nf2);
       double He=0;

       if(tangle<1)
       {
           if(sign>0)
           {
               He=-m_EdgeSize*sqrt(0.5*(1.0-tangle));
           }
           else if(sign<0)
           {
            He=m_EdgeSize*sqrt(0.5*(1.0-tangle));  //He=2*cos(m_Dihedral/2.0)*renorm;
           }
           else
           {
            He=0;
           }
       }
       else if(tangle>=1 && tangle<1.01)
       {
           // in case some numerical probelm happens; 1.01 is too large however,
           He=0;
   
       }
       else if(tangle>1.01)
       {
           std::cout<<"error--->: somthing wrong with this link \n";
           exit(0);
       }
	
       m_Be=Be;
       m_He=He;
       m_mirorlink->PutShapeOperator(m_Be,m_He);


  }
  else
  {

  }


}
void links::UpdateEdgeVector(Vec3D *pBox)
{
        double x1=m_V1->GetVXPos();
        double y1=m_V1->GetVYPos();
        double z1=m_V1->GetVZPos();
        double x2=m_V2->GetVXPos();
        double y2=m_V2->GetVYPos();
        double z2=m_V2->GetVZPos();
        
        double dx1=x2-x1;
        if(fabs(dx1)>(*pBox)(0)/2.0)
        {
            if(dx1<0)
            dx1=(*pBox)(0)+dx1;
            else if(dx1>0)
            dx1=dx1-(*pBox)(0);
        }
        double dy1=y2-y1;
        if(fabs(dy1)>(*pBox)(1)/2.0)
        {
            if(dy1<0)
            dy1=(*pBox)(1)+dy1;
            else if(dy1>0)
            dy1=dy1-(*pBox)(1);
        }
        double dz1=z2-z1;
        if(fabs(dz1)>(*pBox)(2)/2.0)
        {
            if(dz1<0)
            dz1=(*pBox)(2)+dz1;
            else if(dz1>0)
            dz1=dz1-(*pBox)(2);
        }
      
        Vec3D Re(dx1,dy1,dz1);
        m_EdgeVector = Re;
        m_EdgeSize=(m_EdgeVector.norm());
    
        if(this->GetMirrorFlag()==true)
        {
            m_mirorlink->PutEdgeVector(m_EdgeVector*(-1),m_EdgeSize);
        }

}
void links::PutEdgeVector(Vec3D v, double l)
{
    m_EdgeSize   = l;
    m_EdgeVector = v;
}
void links::Flip()
{
    
   if(this->GetMirrorFlag()==true){
    
    triangle *T2 = m_mirorlink->GetTriangle();
    vertex  *V4 = m_mirorlink->GetV3();
    vertex  *v1 = m_V1;
    vertex  *v2 = m_V2;
    vertex  *v3 = m_V3;
    vertex  *v4 = V4;
    links *l1=this->GetNeighborLink1();
    links *l2=this->GetNeighborLink2();
    links *l3=m_mirorlink->GetNeighborLink1();
    links *l4=m_mirorlink->GetNeighborLink2();

       m_V1->RemoveFromNeighbourVertex(m_V2);
       m_V2->RemoveFromNeighbourVertex(m_V1);
       V4->AddtoNeighbourVertex(m_V3);
       m_V3->AddtoNeighbourVertex(V4);

       m_V1->RemoveFromLinkList(this);
       m_V2->RemoveFromLinkList(m_mirorlink);
       V4->AddtoLinkList(this);
       m_V3->AddtoLinkList(m_mirorlink);
       m_V1->RemoveFromTraingleList(T2);
       m_V2->RemoveFromTraingleList(m_T1);
       m_V3->AddtoTraingleList(T2);
       V4->AddtoTraingleList(m_T1);

       this->UpdateNeighborLink1(l2);
       this->UpdateNeighborLink2(l3);
       m_mirorlink->UpdateNeighborLink1(l4);
       m_mirorlink->UpdateNeighborLink2(l1);
       l1->UpdateV3(V4);
       l2->UpdateV3(V4);
       l3->UpdateV3(m_V3);
       l4->UpdateV3(m_V3);
       l1->UpdateNeighborLink1(m_mirorlink);
       l1->UpdateNeighborLink2(l4);
       l2->UpdateNeighborLink1(l3);
       l2->UpdateNeighborLink2(this);
    
    l3->UpdateNeighborLink1(this);
    l3->UpdateNeighborLink2(l2);
    
    l4->UpdateNeighborLink1(l1);
    l4->UpdateNeighborLink2(m_mirorlink);
    

    l1->UpdateTriangle(T2);
    l3->UpdateTriangle(m_T1);
    l4->UpdateTriangle(T2);
    l2->UpdateTriangle(m_T1);
 
  
    m_V1=v4;
    m_V2=v3;
    m_V3=v1;
    V4=v2;
    m_mirorlink->UpdateV(m_V2,m_V1,V4);



    
    int id2=T2->GetTriID();
    int id1=m_T1->GetTriID();
    triangle tm1(id1,m_V1,m_V2,m_V3);
    triangle tm2(id2,m_V2,m_V1,V4);
      *m_T1 = tm1;
   
      *(m_mirorlink->GetTriangle()) = tm2;
    
   
   }
    else
    {
        std::cout<<"error---> a link without a mirror, possibly an edge link, is asked to be flipped, such an action is not possible \n";
        exit(0);
    }
    
}
//calaculates the dot(n1,n2) of this edge trinagle with the face of the mirror. this should not be smaller
// the minangle value
bool links::CheckFaceAngleWithMirrorFace(double &minangle) {
    // Check if mirror link exists
    if (!m_mirorflag) {
        return true;
    }
    // Calculate normal vectors of the faces
    Vec3D n1 = m_T1->GetNormalVector();
    Vec3D n2 = m_mirorlink->GetTriangle()->GetNormalVector();
    
    if(n1.dot(n1,n2)<minangle){
        return false;
    }
    
    // Calculate and return the dot product of the normal vectors
    return true;
}
//calaculates the dot(n1,n2) of this edge trinagle with the face of the next edge. this should not be smaller
// the minangle value
bool links::CheckFaceAngleWithNextEdgeFace(double &minangle){
    
    if (!m_neighborlink1->GetMirrorFlag()) {
        std::cout << "---> error: This operation is not defined for such an edge." << std::endl;
        return 0.0; // Return 0.0 as default value
    }
    Vec3D n1 = m_T1->GetNormalVector();
    Vec3D n2 = m_neighborlink1->GetMirrorLink()->GetTriangle()->GetNormalVector();
    
    if(n1.dot(n1,n2)<minangle){
        return false;
    }
    
    // Calculate and return the dot product of the normal vectors
    return true;
}



