#if !defined(AFX_MESH_H_INCLUDED_)
#define AFX_MESH_H_INCLUDED_
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "CreateMashBluePrint.h"

class MESH
{
public:

    MESH();
    ~MESH();

    inline double GetMinLength()                    const      {return m_MinLength;}
    inline double GetMaxLength()                    const      {return m_MaxLength;}
    inline double GetMinAngle()                     const      {return m_MinAngle;}

    
private:
    std::vector<vertex>         m_Vertex;
    std::vector<triangle>       m_Triangle;
    std::vector<links>          m_Links;
    std::vector<inclusion>      m_Inclusion;
    Vec3D                       m_Box;
    double m_MinLength;
    double m_MaxLength;
    double m_MinAngle;
public:
    std::vector <InclusionType> m_InclusionType;
    std::vector <InclusionType*> m_pInclusionType;
    

    
    std::vector<vertex*>        m_pActiveV; // all the active vertices edge + surf
    std::vector<vertex*>        m_pSurfV; // all the active vertices  surf
    std::vector<vertex*>        m_pEdgeV;  // edge
    
    
    std::vector<links*>         m_pActiveL;   // all the links
    std::vector<links*>         m_pHL;
    std::vector<links*>         m_pMHL;
    std::vector<links*>         m_pEdgeL;
    
    
    std::vector<triangle*>      m_pActiveT;

    std::vector<inclusion*>     m_pInclusion;
    Vec3D                       *m_pBox;
    
    void GenerateMesh(MeshBluePrint meshblueprint);
    MeshBluePrint Convert_Mesh_2_BluePrint(MESH *mesh);
    


    
    
};



#endif
