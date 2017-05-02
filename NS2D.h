#include <AFEPack/EasyMesh.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/Operator.h>
#include <AFEPack/BilinearOperator.h>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>


#include <AFEPack/Functional.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/solver_bicgstab.h>
#include <lac/solver_gmres.h>
#include <lac/solver_minres.h>
#include <lac/sparse_ilu.h>
#include <lac/sparse_mic.h>
#include <base/exceptions.h>
#include <lac/vector.h>

#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <set>
#include<time.h>
#include <fstream>


#define DIM  2

#define PI (atan(1.0) * 4.0)
using namespace std;

/**
 * NS class.
 * 
 */
class NS2D
{

private:
	/// Projection method, only one finite element space. 
	EasyMesh mesh;
	/// geometry template;
	TemplateGeometry<DIM> template_geometry; 
	TemplateGeometry<DIM> twin_template_geometry; 
	/// coordinate transform;
	CoordTransform<DIM,DIM> coord_transform;
	CoordTransform<DIM,DIM> twin_coord_transform;
	/// template of the degree of freedom;
	TemplateDOF<DIM> template_dof;
	TemplateDOF<DIM> twin_template_dof;
	/// the administrator of the basis function;
	BasisFunctionAdmin<double,DIM,DIM> basis_function;
	BasisFunctionAdmin<double,DIM,DIM> twin_basis_function;
	/// element template;
	std::vector<TemplateElement<double,DIM,DIM> > template_element;
	/// finite element space;
	FEMSpace<double,DIM> fem_space;
	/// mesh file;
	std::string mesh_file;
	/// current time of the simulation;
	double t;
	/// current time step of the simulation;
	double dt;
	/// a parameter from outside;
	double a;
	/// current numerical solution;


	std::vector < FEMFunction<double, DIM> > u_h;
	HGeometryTree<DIM> h_tree;/**< 网格树. */
	
	std::vector < IrregularMesh<DIM> * >irregular_mesh_v;  /**< V 空间不规则网格. (链表) */
	std::vector < SparseMatrix<double> > stiff_matrix;//稀疏矩阵模版
	std::vector < SparseMatrix<double> * > matrix;//稀疏矩阵

	std::vector < FEMSpace<double, DIM> > fem_space_v; /**< 速度空间. */
	
	std::vector <SparsityPattern> sp_stiff_matrix;//细网格稀疏矩阵模板
	std::vector < Vector<double> > rhs;//方程的右端项；
	
	std::vector< std::vector <int> > locate;//粗网格上序号为i的点在细网格的编号；第i个locate代表i层到i+1层的关系；
	std::vector< std::vector <int> > flag;//网格编号，说明该点属于哪层网格;
	std::vector< std::vector< std::set<int> > > rconnect;//第i个集合存储与i点相邻的粗网格点的编号，i为细网格上的细网格点编号，i为细网格上的粗网格点时值为1
	std::vector< std::vector< std::set<int> > > connect;//第i个集合存储与i点相邻的细网格点的编号，i为粗网格上的粗网格点编号;
	
	int mesh_l;//多重网格层数,网格编号由粗到细为0-mesh_l-1
	double tol;			// 收敛精度. 
	int V1,V2;
	
	//指示器，用于表明element需要加密还是粗化
	std::vector< Indicator<DIM> > indicator;
	
public:
	/** 
	 * Constructor.
	 * 
	 * @param file mesh file name.
	 */
	NS2D(const std::string& file);
	
	/** 
	 * A standard destructor.
	 * 
	 */
	virtual ~ NS2D();
public:
	/** 
	 * Main procedure.
	 * 
	 */
	void run();
/** 
	 * 读入配置文件.
	 * 
	 * @param _config_file 配置文件名.
	 */
	void config(std::string _config_file);

	/** 
	 * The simulation develop for one time step.
	 * 
	 */
	void stepForward();

	/** 
	 * All the prepare things.
	 * 
	 */
	void initialize();
	
	void getindex(int n);
	
	void Construct_Matrix(int n_dof_p,int flag);
	
	void Matrix_build(int n);
	
	void Matrix_project(int n , SparseMatrix<double>& A);
	
	void GaussSidel(const SparseMatrix<double> & M, Vector<double>& x, const Vector<double>& r,const int& s) const ; 
	
	Vector<double>  matrix_multi(const SparseMatrix<double>& P,Vector<double>& rhs); 	
	
	void solve();
	
	void VMG(int n,int v1,int v2);
	
	/** 
	 * Arrange the initial values of the simulation.
	 * 
	 */
	void initialValue();

	/** 
	 * Set the boundary values here, but may useless.
	 * 
	 */
	void boundaryValue();

	/** 
	 * Output the current numerical solution.
	 * 
	 */
	virtual void outputSolution();
	
};

//
// end of file
//////////////////////////////////////////////////////////////////////////////
