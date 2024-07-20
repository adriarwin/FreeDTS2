#if !defined(AFX_vertex_H_9P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_)
#define AFX_vertex_H_9P4B21B8_C13C_5648_BF23_124095086234__INCLUDED_
#include "Voxel.h"
#include "SimDef.h"
#include "CNTCell.h"
#include "Vec3D.h"
#include "Tensor2.h"
#include "inclusion.h"
/*******************
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Vertex object: can be defined by and integer and three double numbers for the positions of the vertices
 Also there are a few lists pointing to the nighbouring links and trinagles and also the inclusion if the vertex own one.
 *******************/
class links;
class triangle;
class vertex
{
public:
	vertex(int id, double x, double y, double z);
	vertex(int id);
    vertex();

	 ~vertex();

// A set of function to get vertex variables (since all of them are private)
	    inline int GetVID()                                 {return m_ID;}
        inline double GetVXPos()                            {return m_X;}
        inline double GetVYPos()                            {return m_Y;}
        inline double GetVZPos()                            {return m_Z;}
        inline double GetXPos()                             {return m_X;}
        inline double GetYPos()                             {return m_Y;}
        inline double GetZPos()                             {return m_Z;}
        inline double GetArea()                             {return m_Area;}
        inline Tensor2  GetL2GTransferMatrix()              {return m_T_Local_2_Global;}
        inline Tensor2  GetG2LTransferMatrix()              {return m_T_Global_2_Local;}
        inline Vec3D GetNormalVector()                      {return m_Normal;}
        inline std::vector <double> GetCurvature()          {return m_Curvature;}// surface curvature
        inline double GetEnergy()                           {return m_Energy;}
        inline CNTCell * GetVCNTCell()                      {return m_CNTCell;}
        inline inclusion* GetInclusion()                    {return m_pInclusion;}
        inline bool VertexOwnInclusion()                    {return m_OwnInclusion;}
        inline Vec3D *GetBox()                              {return m_pBox;}
        inline int GetSimTimeStep()                         {return m_SimTimeStep;}
        inline int GetGroup()                               {return m_Group;}
        inline std::string GetGroupName()                   {return m_GroupName;}
        inline std::vector <links *> GetVLinkList()             {return m_VLinkList;}
        inline std::vector <triangle *> GetVTraingleList()         {return m_VTraingleList;}
        inline std::vector <vertex *> GetVNeighbourVertex()     {return m_VNeighbourVertex;}
    
public:
    
    bool SetCopy();            // Copies the key ellements into the old type
    bool Reverse2PreviousCopy();  // reverse the edge to the value set at the time of MakeCopy()
  // A set of functions to update vertex variables
  void UpdateVXPos(double x);    // a function for update the x position of a vertex
  void UpdateVYPos(double y);   // a function for update the y position of a vertex
  void UpdateVZPos(double z);   // a function for updates the z position of a vertex
  void UpdateGroupName(std::string z); // A vertex can only have one group name, is different from group id
  void UpdateVCNTCell(CNTCell * z);
  void UpdateBox(Vec3D *z);
  void UpdateCurvature(double,double); // vis
  void UpdateEnergy(double); 
  void UpdateNormal_Area(Vec3D,double); // vis
  void UpdateL2GTransferMatrix(Tensor2 v);
  void UpdateG2LTransferMatrix(Tensor2 v);
  void UpdateInclusion(inclusion* );
  void UpdateOwnInclusion(bool );
  void UpdateVID(int); // this should never be called, only for temporarily vis functions
  void AddtoLinkList(links* z);
  bool AddtoLinkListCarefully(links* z);
  void AddtoTraingleList(triangle * z);
  bool AddtoTriangleListCarefully(triangle* t);
  void AddtoNeighbourVertex(vertex* z);
  bool AddtoNeighbourVertexCarefully(vertex* v);
  void RemoveFromLinkList(links* z);
  bool RemoveFromLinkListCarefully(links* l);
  void RemoveFromTraingleList(triangle * z);
  void RemoveFromNeighbourVertex(vertex* z);
  void UpdateGroup(int z);
  void UpdateVoxel(Voxel<vertex> * pVoxel);

  void UpdateSimTimeStep(int v);   // we should remove this function at some point

public:
    bool CheckCNT();
    bool VertexMoveIsFine(double dx,double dy, double dz,  double mindist, double maxdist);
    bool CheckVoxel();
    bool UpdateVoxelAfterAVertexMove(); // update the voxel after the move has happened.
    //-- checks the face angles of trinagle around the vertex and with the next triangles
    bool CheckFacesAfterAVertexMove(double &minangle);  // this checks if the faces are fine, if not, the move need to be rejected.

private:
    double SquareDistanceFromAVertex(vertex* pv2);
    double SquareDistanceOfAVertexFromAPoint(double X, double Y, double Z, vertex* pv2);
private:

    int m_ID;         // ID of the vertex, a unique number
    double m_X;       // X coordinate
    double m_Y;       // Y coordinate
    double m_Z;       // Z coordinate
    std::vector <triangle *> m_VTraingleList;  // A list holding all the nighbouring triangles
    std::vector <links *> m_VLinkList;      // A list holding all the nighbouring edges (linkss)
    std::vector <vertex *> m_VNeighbourVertex; // A list holding all the nighbouring vertexes
    inclusion *m_pInclusion;                    // pointer to an inclusion that the vertex hold (could be empty)
    bool m_OwnInclusion;                        // to check if the vertex own any inclusion
    double m_Area;                              // area of the vertex
    CNTCell * m_CNTCell;                        // a unitcell that the vertex belong to at any point of the simulation, it will be chnage during a simulation
    int m_SimTimeStep;                          // some extera access (should be removed )
    int m_Group;            // Id of a group that the vertex belong too
    private:
    Vec3D m_Normal;
    Vec3D *m_pBox;
    std::vector<double> m_Curvature;
    double m_Energy;
    Tensor2  m_T_Local_2_Global;         //  Local to global transformation matrix
    Tensor2  m_T_Global_2_Local;        //  global to local transformation matrix
    std::string m_GroupName;
     Voxel<vertex> * m_pVoxel;
   // Voxel * m_pVoxel;

    
    
public:
    // lets have them public for now
    //================================================
    // development of Aug 2023; for membranes with hole
    //a*kg^2+b*k^2    note: any spontaneous one will come from inclusions
    
    // we might need to define n and t and b as well, for now lets ignore it
    //================================================
    double m_Geodesic_Curvature;          // Edge Vertex Curvature
    double m_Normal_Curvature;          // Edge Vertex Curvature
    double m_VLength;                       // length of the vertex
    double m_Lambda;                   // line tension
    int m_VertexType;                   // 0 surface vertex; 1 edge vertex;
    links * m_pEdgeLink;
    links * m_pPrecedingEdgeLink;// preceding link at the edge
    
    
// members for copying.
private:
    double m_OldX;       // X coordinate
    double m_OldY;       // Y coordinate
    double m_OldZ;       // Z coordinate
    std::vector <triangle *> m_OldVTraingleList;  // A list holding all the nighbouring triangles
    std::vector <links *> m_OldVLinkList;      // A list holding all the nighbouring edges (linkss)
    std::vector <vertex *> m_OldVNeighbourVertex; // A list holding all the nighbouring vertexes
    inclusion *m_OldpInclusion;                    // pointer to an inclusion that the vertex hold (could be empty)
    bool m_OldOwnInclusion;                        // to check if the vertex own any inclusion
    double m_OldArea;                              // area of the vertex
    CNTCell * m_OldCNTCell;                        // a unitcell that the vertex belong to at any point of the simulation, it will be chnage during a simulation
    int m_OldGroup;            // Id of a group that the vertex belong too
    Vec3D m_OldNormal;
    std::vector<double> m_OldCurvature;
    double m_OldEnergy;
    Tensor2  m_OldT_Local_2_Global;         //  Local to global transformation matrix
    Tensor2  m_OldT_Global_2_Local;        //  global to local transformation matrix
    std::string m_OldGroupName;
    Voxel<vertex> * m_OldpVoxel;
    double m_OldGeodesic_Curvature;          // Edge Vertex Curvature
    double m_OldNormal_Curvature;          // Edge Vertex Curvature
    double m_OldVLength;                       // length of the vertex
    double m_OldLambda;                   // line tension
    int m_OldVertexType;                   // 0 surface vertex; 1 edge vertex;
    links * m_OldpEdgeLink;
    links * m_OldpPrecedingEdgeLink;// preceding link at the edge
    
};


#endif
