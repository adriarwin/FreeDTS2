

#include <stdio.h>
#include "EvolveVerticesByMetropolisAlgorithm.h"
#include "State.h"

EvolveVerticesByMetropolisAlgorithm::EvolveVerticesByMetropolisAlgorithm(State *pState)
    : m_pState(pState),
      m_pSurfV(pState->GetMesh()->GetSurfV()),
      m_pEdgeV(pState->GetMesh()->GetEdgeV()),
      m_Beta(pState->GetSimulation()->GetBeta()),
      m_DBeta(pState->GetSimulation()->GetBeta()),
      m_MinLength2(pState->GetSimulation()->GetMinL2()),
      m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
      m_MinAngle(pState->GetSimulation()->GetMinAngle()){
          
          m_NumberOfMovePerStep_Surf = 1.0;
          m_NumberOfMovePerStep_Edge = 1.0;
          m_DR = 0.05;
          
      }
EvolveVerticesByMetropolisAlgorithm::EvolveVerticesByMetropolisAlgorithm(State *pState, double rate_surf, double rate_edge, double dr)
    : m_pState(pState),
      m_pSurfV(pState->GetMesh()->GetSurfV()),
      m_pEdgeV(pState->GetMesh()->GetEdgeV()),
      m_Beta(pState->GetSimulation()->GetBeta()),
      m_DBeta(pState->GetSimulation()->GetBeta()),
      m_MinLength2(pState->GetSimulation()->GetMinL2()),
      m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
      m_MinAngle(pState->GetSimulation()->GetMinAngle()) {
          
          m_NumberOfMovePerStep_Surf = rate_surf;
          m_NumberOfMovePerStep_Edge = rate_edge;
          m_DR = dr;
      }
EvolveVerticesByMetropolisAlgorithm::~EvolveVerticesByMetropolisAlgorithm(){

}
void EvolveVerticesByMetropolisAlgorithm::Initialize(){
    m_pBox = m_pState->GetMesh()->GetBox();
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;

    return;
}
bool EvolveVerticesByMetropolisAlgorithm::EvolveOneStep(int step){
 
    int no_surf_v = m_pSurfV.size();
    int no_edge_v = m_pEdgeV.size();
    int no_steps_edge = no_edge_v*m_NumberOfMovePerStep_Edge;
    int no_steps_surf = no_surf_v*m_NumberOfMovePerStep_Surf;

  for (int i = 0; i< no_steps_surf;i++) {
    
      int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_surf_v);
      vertex *pvertex = m_pSurfV[r_vid];
      if( m_FreezGroupName == pvertex->GetGroupName()){
          continue;
      }
      double dx=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
      double dy=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
      double dz=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
      
      if(!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(dx,dy,dz, pvertex)){
          continue;
      }
      
      double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
      if(EvolveOneVertex(step, pvertex, m_DR*dx, m_DR*dy, m_DR*dz,thermal)){
          m_AcceptedMoves++;
      }
      m_NumberOfAttemptedMoves++;
    }
    
    for (int i = 0; i< no_steps_edge;i++) {
      
        int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_edge_v);
        vertex *pvertex = m_pEdgeV[r_vid];
        if( m_FreezGroupName == pvertex->GetGroupName()){
            continue;
        }
        double dx=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        double dy=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        double dz=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        if(!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(dx,dy,dz, pvertex)){
            continue;
        }
        
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

        if(EvolveOneVertex(step, pvertex, m_DR*dx, m_DR*dy, m_DR*dz,thermal)){
            m_AcceptedMoves++;
        }
        m_NumberOfAttemptedMoves++;
    }
    
    return true;
}
bool EvolveVerticesByMetropolisAlgorithm::EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp){
    
    double old_energy = 0;
    double new_energy = 0;

//---> first checking if all the distances will be fine if we move the vertex
    if(!VertexMoveIsFine(pvertex,dx,dy,dz,m_MinLength2,m_MaxLength2))  // this function could get a booling varaible to say, it crossed the voxel
        return 0;
//---> lets get global variables before moving
    //--- obtain vertices energy terms
    old_energy = pvertex->GetEnergy();
    const std::vector<vertex *>& vNeighbourV = pvertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
        
        old_energy+=(*it)->GetEnergy();
    }
    // find the links in which there interaction energy changes
    std::vector<links*> Affected_links = GetEdgesWithInteractionChange(pvertex);
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        old_energy += 2*(*it)->GetIntEnergy();
    }
    // --- obtaining global variables that can change by the move. Note, this is not the total volume, only the one that can change.
     double old_Tvolume = 0;
     double old_Tarea = 0;
     double old_Tcurvature = 0;
     double new_Tvolume = 0;
     double new_Tarea = 0;
     double new_Tcurvature = 0;

//--->
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateAVertexRingContributionToGlobalVariables(pvertex, old_Tvolume, old_Tarea, old_Tcurvature);
    }

//----> Move the vertex;
        pvertex->PositionPlus(dx,dy,dz);
    //--- update triangles normal
    std::vector<triangle *> N_triangles = pvertex->GetVTraingleList();
    for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
        (*it)->Copy();
        (*it)->UpdateNormal_Area(m_pBox);
    }
    //  check new faces angles, if bad, reverse the trinagles
    if(!CheckFacesAfterAVertexMove(pvertex)){
        
        for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
            (*it)->Reverse2PreviousCopy();
        }
        return false;
    }
//---->
    //--> clauclate edge shape operator;
    const std::vector<links *>& v_NLinks = pvertex->GetVLinkList();
    for (std::vector<links *>::const_iterator it = v_NLinks.begin() ; it != v_NLinks.end(); ++it){
        
        (*it)->SetCopy();
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
                
        (*it)->GetNeighborLink1()->SetCopy();
        (*it)->GetNeighborLink1()->UpdateNormal();
        (*it)->GetNeighborLink1()->UpdateShapeOperator(m_pBox);
    }

    // --> calculate vertex shape operator
    pvertex->SetCopy();
    (m_pState->GetCurvatureCalculator())->UpdateVertexCurvature(pvertex);
    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
        (*it)->SetCopy();
        (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(*it);
    }

    //---> calculate new energies
    new_energy = (m_pState->GetEnergyCalculator())->SingleVertexEnergy(pvertex);

    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
        new_energy += (m_pState->GetEnergyCalculator())->SingleVertexEnergy(*it);
    }
    
    //---> get energy for ApplyConstraintBetweenGroups
    Vec3D Dx(dx,dy,dz);
    double dE_Cgroup = m_pState->GetApplyConstraintBetweenGroups()->CalculateEnergyChange(pvertex, Dx);
    
    //---> for now, only active nematic force: ForceonVerticesfromInclusions
    double dE_force_from_inc  = m_pState->GetForceonVerticesfromInclusions()->Energy_of_Force(pvertex, Dx);
    
//---> new global variables
    if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
        m_pState->GetVAHGlobalMeshProperties()->CalculateAVertexRingContributionToGlobalVariables(pvertex, new_Tvolume, new_Tarea, new_Tcurvature);
    }
    //---> energy change of global variables
    double dE_volume =  m_pState->GetVolumeCoupling()->GetEnergyChange(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
    double dE_t_area = m_pState->GetTotalAreaCoupling()->CalculateEnergyChange(old_Tarea, new_Tarea);
    double dE_g_curv = m_pState->GetGlobalCurvature()->CalculateEnergyChange(new_Tarea-old_Tarea, new_Tcurvature-old_Tcurvature);
    
    //--> only elatsic energy
    double diff_energy = new_energy - old_energy;
    //--> sum of all the energies
    double tot_diff_energy = diff_energy + dE_Cgroup + dE_force_from_inc + dE_volume + dE_t_area + dE_g_curv ;

    //---> accept or reject the move
    if(tot_diff_energy <= 0 || exp(-m_Beta * tot_diff_energy + m_DBeta) > temp ) {

        // move is accepted
        (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);
        //---> if vertex is out of the voxel, update its voxel
        if(!pvertex->CheckVoxel()){
            pvertex->UpdateVoxelAfterAVertexMove();
        }
        //---> ApplyConstraintBetweenGroups
        m_pState->GetApplyConstraintBetweenGroups()->AcceptMove();
        
        //---> global variables
        if(m_pState->GetVAHGlobalMeshProperties()->GetCalculateVAH()){
             m_pState->GetVolumeCoupling()->UpdateArea_Volume(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
            m_pState->GetTotalAreaCoupling()->UpdateTotalArea(old_Tarea, new_Tarea);
            m_pState->GetGlobalCurvature()->UpdateEnergyChange(new_Tarea-old_Tarea, new_Tcurvature-old_Tcurvature);
        }
    }
    else {
//---> reverse the changes that has been made to the system
        //---> reverse the triangles
        for (std::vector<triangle *>::iterator it = N_triangles.begin() ; it != N_triangles.end(); ++it){
            
            (*it)->Reverse2PreviousCopy();
        }
        //---> reverse the links
        for (std::vector<links *>::const_iterator it = v_NLinks.begin() ; it != v_NLinks.end(); ++it){
            
            (*it)->Reverse2PreviousCopy();
            (*it)->GetNeighborLink1()->Reverse2PreviousCopy();
        }
        //--> the shape operator of these links has not been affected, therefore we only update the interaction energy
        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            (*it)->Reverse_InteractionEnergy();
        }
        //---> reverse the vertices
        pvertex->Reverse2PreviousCopy();
        for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
            (*it)->Reverse2PreviousCopy();
        }
     }

    return true;
}
//---> this does not check the angle of the faces. Because this should be done after the move:
//waste of calculation if we do ith before the move. Unless, we store the values. That also not good because move could get rejected.
bool EvolveVerticesByMetropolisAlgorithm::VertexMoveIsFine(vertex* pvertex, double dx,double dy, double dz,  double mindist2, double maxdist2){
//--->  vertex new position if get accepted
    
        double new_x = pvertex->GetXPos() + dx;
        double new_y = pvertex->GetYPos() + dy;
        double new_z = pvertex->GetZPos() + dz;
    
//--->  let check the distances with the nighbours
    std::vector <vertex *> npvertex = pvertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::iterator it = npvertex.begin() ; it != npvertex.end(); ++it){
        double dist2 = (*it)->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it);
            if(dist2<mindist2 || dist2>maxdist2)
            return false;
    }
//---> now check it within the voxel cells
//---> lets get the object voxel, note: the voxel could be different from the associated one as it is moving
        //-- obtain the object (vertex) new cell id, with respect to the current cell
        int i = int(new_x/pvertex->GetVoxel()->GetXSideVoxel((*m_pBox)(0)))-pvertex->GetVoxel()->GetXIndex();
        int j = int(new_y/pvertex->GetVoxel()->GetYSideVoxel((*m_pBox)(1)))-pvertex->GetVoxel()->GetYIndex();
        int k = int(new_z/pvertex->GetVoxel()->GetZSideVoxel((*m_pBox)(2)))-pvertex->GetVoxel()->GetZIndex();
        //-- check if it has moved too far
        if(i > 1 || i < -1 || j > 1 || j < -1 || k > 1 || k < -1) {
            std::cout << " ---> warning: the object might moved more than one voxel " << std::endl;
            return false;
        }

        Voxel<vertex>* new_pvox = pvertex->GetVoxel()->GetANeighbourCell(i, j, k);
    
        for(int n=-1;n<2;n++)
        for(int m=-1;m<2;m++)
        for(int s=-1;s<2;s++){
            std::vector <vertex *> CV = new_pvox->GetANeighbourCell(n, m, s)->GetContentObjects();
            for (std::vector<vertex *>::iterator it = CV.begin() ; it != CV.end(); ++it){
                if(*it!=pvertex){
                    if((*it)->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it)<mindist2)
                        return false;
                }
                
            }
        } ///   for(int s=-1;s<2;s++){
    return true;
}
// finding the distance of the current vertex from the pv2; also considering the pbc conditions
bool EvolveVerticesByMetropolisAlgorithm::CheckFacesAfterAVertexMove(vertex* p_vertex) {
    std::vector<links*> linkList = p_vertex->GetVLinkList();
    for (std::vector<links*>::iterator it = linkList.begin(); it != linkList.end(); ++it) {
        links* link = *it;
        if (!link->CheckFaceAngleWithMirrorFace(m_MinAngle) || !link->CheckFaceAngleWithNextEdgeFace(m_MinAngle)) {
            return false;
        }
    }
    return true;
}

//========
// this function can be deleted any time; it is for test cases only
double  EvolveVerticesByMetropolisAlgorithm::SystemEnergy()
{
    /*
    MESH* m_pMESH = m_pState->m_pMesh;
    std::vector<vertex *> ActiveV = m_pMESH->m_pActiveV;
    std::vector<triangle *> pActiveT = m_pMESH->m_pActiveT;
    std::vector<links *> mLink = m_pMESH->m_pHL;
    std::vector<links *>  pEdgeL = m_pMESH->m_pEdgeL;
    std::vector<vertex *> EdgeV  = m_pMESH->m_pEdgeV;
    double en = 0;
    

    for (std::vector<triangle *>::iterator it = pActiveT.begin() ; it != pActiveT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    
    
    for (std::vector<links *>::iterator it = mLink.begin() ; it != mLink.end(); ++it)
    {
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    
    
    for (std::vector<vertex *>::iterator it = ActiveV.begin() ; it != ActiveV.end(); ++it)
        (m_pState->CurvatureCalculator())->SurfVertexCurvature(*it);

    //====== edge links should be updated
    for (std::vector<links *>::iterator it = pEdgeL.begin() ; it != pEdgeL.end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = EdgeV.begin() ; it != EdgeV.end(); ++it)
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(*it);
    
    en=m_pEnergyCalculator->TotalEnergy(ActiveV,mLink);
    //en=en+ m_pEnergyCalculator->TotalEnergy(EdgeV,pEdgeL);
   
    
    return en;
     */
    return 0;
}
std::string EvolveVerticesByMetropolisAlgorithm::CurrentState(){
    
    std::string state = AbstractVertexPositionIntegrator::GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();
    state = state +" "+ Nfunction::D2S(m_NumberOfMovePerStep_Surf) +" "+Nfunction::D2S(m_NumberOfMovePerStep_Edge);
    state = state +" "+ Nfunction::D2S(m_DR);
    return state;
}
std::vector<links*> EvolveVerticesByMetropolisAlgorithm::GetEdgesWithInteractionChange(vertex* p_vertex){
    
    std::vector<links*> edge_with_interaction_change;
    std::vector<links *> all_temlinks;
    
    if( p_vertex->GetVertexType() == 1 )
        all_temlinks.push_back(p_vertex->GetPrecedingEdgeLink());
    
    
    std::vector<vertex *> neighbor_vertices = p_vertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::iterator it = neighbor_vertices.begin() ; it != neighbor_vertices.end(); ++it)
    {
            if((*it)->VertexOwnInclusion()==true) {
                std::vector<links *> ltem=(*it)->GetVLinkList();
                all_temlinks.insert(all_temlinks.end(), ltem.begin(), ltem.end());
                if((*it)->m_VertexType==1)
                    all_temlinks.push_back((*it)->m_pPrecedingEdgeLink);
            }
    }
    //-- now we remove the repeated links
    for (std::vector<links *>::iterator it = all_temlinks.begin() ; it != all_temlinks.end(); ++it)
    {
        bool Should_be_added = true;
        for (std::vector<links *>::iterator it2 = edge_with_interaction_change.begin() ; it2 != edge_with_interaction_change.end(); ++it2)
        {
            if((*it2)==(*it)){
                Should_be_added = false;
            }
            else if((*it2)->GetMirrorFlag() && (*it2)->GetMirrorLink() == (*it)){
                Should_be_added = false;
            }
        }
        if(Should_be_added == true)
            edge_with_interaction_change.push_back((*it));

    }
    
    return edge_with_interaction_change;
}
/*
std::vector<links*> GetEdgesWithInteractionChange(vertex* p_vertex){
    
    std::vector<links*> edge_with_interaction_change;
    std::unordered_set<links*> all_temlinks_set; // Change to unordered_set

    if( p_vertex->GetVertexType() == 1 )
        all_temlinks_set.insert(p_vertex->GetPrecedingEdgeLink());

    std::vector<vertex*> neighbor_vertices = p_vertex->GetVNeighbourVertex();
    for (std::vector<vertex*>::iterator it = neighbor_vertices.begin(); it != neighbor_vertices.end(); ++it) {
        if((*it)->VertexOwnInclusion()) {
            std::vector<links*> ltem = (*it)->GetVLinkList();
            all_temlinks_set.insert(ltem.begin(), ltem.end());
            if((*it)->m_VertexType == 1)
                all_temlinks_set.insert((*it)->m_pPrecedingEdgeLink);
        }
    }

    // Copy unique links from unordered_set to vector
    std::vector<links*> all_temlinks(all_temlinks_set.begin(), all_temlinks_set.end());

    
    for (std::vector<links *>::iterator it = all_temlinks.begin() ; it != all_temlinks.end(); ++it)
    {
        bool addit=true;
        for (std::vector<links *>::iterator it2 = edge_with_interaction_change.begin() ; it2 != edge_with_interaction_change.end(); ++it2)
        {
                if((*it2)->GetMirrorFlag()==true && ((*it2)->GetMirrorLink())==(*it)){
                    addit=false;
                }
        }
        if(addit==true)
            edge_with_interaction_change.push_back((*it));
    }
    
    return edge_with_interaction_change;
}
*/

