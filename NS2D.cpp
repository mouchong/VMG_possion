
/**
 * @file   NS2D.cpp
 * @author MC
 * @date  
 *
 * @brief VMG for 2D possion equations 
 * solver.
 *
 *
 */

#include "NS2D.h"
#include <sstream>
#include <AFEPack/Functional.h>
#define PI (atan(1.0) * 4.0)

/**
 ＊ 维数在NS2D.h中声明改变,
 ＊ 以下维数不同时,用if语句选择在地方有三处,
 ＊ 	1:l63.读入初始网格文件;
 ＊ 	2:l97.读入几何信息文件;
 ＊ 3:l416.获取粗细网格间索引时,由粗网格总自由度个数确定细网格总自由度个数时;
 ＊初步时间测算三维VMG结果(4次加密):cost=G_S迭代部分(主要耗时是在最细网格上的部分)+非G_S迭代部分;
 ＊若以非G_S迭代部分耗时为100
 ＊G_S迭代部分,V2每+1,耗时增加10-15
 ＊加密次数增加,G_S迭代部分占总时间的比例也在增加;
 ＊
*/
double g_t;
double g_a;

double u(const double * p)//
{
        return sin(PI * p[0]) * sin(2 * PI * p[1]) * sin(PI * p[2]);
};

double f(const double * p)
{
        return 6 * PI * PI * u(p);
};

double boundary_value(const double *);

double boundary_value(const double *p)
{
	return 0;
};

double zero_bnd(const double *);

double zero_bnd(const double *p)
{
	return 0;
};

double real_solution(const double *);

double real_solution(const double *p)
{
	double result = (1.0 - p[0] * p[0]) / 2.0;
	double s = 0.0;
	/// formula (1.5), in page 11.
	for (int k = 1; k < 100; k = k + 2)
	{
		s += sin(k * PI * (1.0 + p[0]) / 2.0) 
			/ (k * k * k * sinh(k * PI)) 
			* (sinh(k * PI * (1.0 + p[1]) / 2.0) 
			   + sinh(k * PI * (1.0 - p[1]) / 2.0)); 
	}

	return (result - 16.0 / (PI * PI * PI) * s);
};

///////////////////////////////////////////////////////////////////////////////

NS2D::NS2D(const std::string& file) :
	mesh_file(file), a(0.005), t(0.25)
{
	g_a = a;
	g_t = t;
};


/**
 * Standard decontructor.
 *
 */
NS2D::~NS2D()
{};

/**
 * All the prepare things here.
 *	
 */
 //主函数
 //考虑全局加密，则一次性将所有层的网格和有限元空间准备好；
void NS2D::initialize()
{
	/// 读取配置文件,多重网格层数和收敛精度,VMG(V1,V2)的V1，V2
	config("config");
	std::cout << "The mesh tree is " << mesh_file << std::endl;
	/// 读入网格.
	//mesh_l本次网格层数，由粗到细编号0-mesh_l-1
 	//n该函数细网格所在层数
 	if(DIM==3)
	h_tree.readMesh(mesh_file);
	else if(DIM==2)
	h_tree.readEasyMesh(mesh_file);
	
	//给出参数向量大小.
	irregular_mesh_v.resize(mesh_l);
	sp_stiff_matrix.resize(mesh_l);
	fem_space_v.resize(mesh_l);
	u_h.resize(mesh_l);
	stiff_matrix.resize(mesh_l);
	matrix.resize(mesh_l);
	
	locate.resize(mesh_l-1);
	flag.resize(mesh_l-1);
	rconnect.resize(mesh_l-1);
	connect.resize(mesh_l-1);
	
	for(int i=0;i<mesh_l;i++)
		irregular_mesh_v[i] = new IrregularMesh<DIM>;
	
	//读取网格并建立两层网格之间的索引.
	irregular_mesh_v[0]->reinit(h_tree);
	for(int i=0;i<mesh_l-1;i++)
		getindex(i);
	
	std::cout << "index settled,mesh complete." << std::endl;
	
	/// 基于维数，读入几何信息;
	//要用到局部加密网格的话，二维时需要三角形和双生三角形的配置信息来生成有限元空间；
	if(DIM==3)
	{
		template_geometry.readData("tetrahedron.tmp_geo");
		coord_transform.readData("tetrahedron.crd_trs");
		template_dof.reinit(template_geometry);
		template_dof.readData("tetrahedron.1.tmp_dof");
		basis_function.reinit(template_dof);
		basis_function.readData("tetrahedron.1.bas_fun");
	}
	else if(DIM==2)
	{
		template_geometry.readData("triangle.tmp_geo");
		coord_transform.readData("triangle.crd_trs");
		template_dof.reinit(template_geometry);
		template_dof.readData("triangle.1.tmp_dof");
		basis_function.reinit(template_dof);
		basis_function.readData("triangle.1.bas_fun");
	}
	
	template_element.resize(1);
	template_element[0].reinit(template_geometry, template_dof,coord_transform, basis_function);
	
	
	/// 根据不同层的网格, 建立各自的空间.
	for(int i=0;i<mesh_l;i++)
	fem_space_v[i].reinit(irregular_mesh_v[i]->regularMesh(), template_element);
	
	/// 建立各自的有限元空间.
	std::vector<int> n_ele(mesh_l);
	for(int i=0;i<mesh_l;i++)
	{
		n_ele[i] = irregular_mesh_v[i]->regularMesh().n_geometry(DIM);
		fem_space_v[i].element().resize(n_ele[i]);
		for (int j = 0; j < n_ele[i]; ++j)
			fem_space_v[i].element(j).reinit(fem_space_v[i], j, 0);
		fem_space_v[i].buildElement();
		fem_space_v[i].buildDof();
		fem_space_v[i].buildDofBoundaryMark();
		u_h[i].reinit(fem_space_v[i]);
	}
	
	std::cout << "space complete." << std::endl;
	
	//建立每层网格的稀疏矩阵模板A
	std::vector <int > n_dof_v;
	n_dof_v.resize(mesh_l);
	
	for(int i=0;i<mesh_l;++i)
	{
		n_dof_v[i]=fem_space_v[i].n_dof();
		Construct_Matrix(n_dof_v[i],i);
		stiff_matrix[i].reinit(sp_stiff_matrix[i]);
		matrix[i]=&stiff_matrix[i];
	}

	rhs.resize(mesh_l);
	for(int i=0;i<mesh_l;++i)
		rhs[i].reinit(n_dof_v[i]);
	
	//VMG只要拼装出最细网格的矩阵;
	Matrix_build(mesh_l-1);
	//通过几何关系映射得到所有粗网格的矩阵;
	//通过几何关系映射得到粗网格的矩阵,层数n->n-1;
	for(int n=mesh_l-1;n>0;--n)
	{
		const SparsityPattern& spMf = matrix[n]->get_sparsity_pattern();
		const std::size_t * f_rowstart = spMf.get_rowstart_indices();
		const unsigned int * f_colnum = spMf.get_column_numbers();
		
		const SparsityPattern& spMc = matrix[n-1]->get_sparsity_pattern();
		const std::size_t * c_rowstart = spMc.get_rowstart_indices();
		const unsigned int * c_colnum = spMc.get_column_numbers();
		
		//粗网格的a(n-1)(ii)=a(n)(ii)+0.25*求和(a(n)(jk)|j,k为i在细网格上的零点)
		//按行遍历粗网格
		for (int i = 0;i < matrix[n-1]->n();i ++)
		{
			//找到第i行的非零元，开始遍历
			for (int j = c_rowstart[i];j < c_rowstart[i+1];j ++)
			{
				//a是j点的纵坐标
				const int& a = c_colnum[j];
				double temp=0;
				//基于几何信息计算a(i,j),遍历所有细网格上的非零元,k为细网格上与i点相邻的点(包括本身),l为细网格上与a点相邻的点(包括本身)
				for(set<int>::iterator k=connect[n-1][i].begin();k!=connect[n-1][i].end();k++)
				{
					for(set<int>::iterator l=connect[n-1][a].begin();l!=connect[n-1][a].end();l++)
					{
					//基于几何信息，设定映射权重
						double p1=0.5,p2=0.5;
						if((i==locate[n-1][*k])&&(i!=0))
							p1=1.0;
						if((i==0)&&(*k==0))
							p1=1.0;
							
						if((a==locate[n-1][*l])&&(a!=0))
							p2=1.0;
						if((a==0)&&(*l==0))
							p2=1.0;
						temp+=(matrix[n]->el(*k,*l)*p1*p2);
					}
				}
				matrix[n-1]->add(i, a , temp);
			}
		}
	}
	std::cout << "initialize complete, press a key to continue." << std::endl;
	getchar();
	
};

void NS2D::solve()
{
	//以上准备工作完成，下面开始求解线性方程组，Au=f；
	double ehs=1;
	int iter_t=0;
	double ehs_t;
	while(ehs>tol)
	{
		VMG(mesh_l-1,V1,V2);
		//GaussSidel(*matrix[mesh_l-1],u_h[mesh_l-1],rhs[mesh_l-1],1);
		//下面基于改变了的u_h[mesh_l-1]计算残差
		Vector<double> temp(matrix[mesh_l-1]->n());
		temp=matrix_multi(*matrix[mesh_l-1],u_h[mesh_l-1]);
		for(int i=0;i<matrix[mesh_l-1]->n();i++)
			temp(i)=rhs[mesh_l-1](i)-temp(i);
		ehs_t=ehs;
		double error1=0;
		for(int i=0;i<matrix[mesh_l-1]->n();i++)
			error1+=temp(i)*temp(i);
		ehs=sqrt(error1);
		//std::cout<<"ehs= "<<ehs<<std::endl;
		//std::cout<<"con rate= "<<ehs/ehs_t<<std::endl;
		std::cout<<ehs<<std::endl;
		iter_t++;
	}
	std::cout<<"converge with residual= "<<ehs<<" iteration times= "<<iter_t<<std::endl;

//	double l2error = Functional::L2Error(u_h[mesh_l-1], FunctionFunction<double>(&real_solution), 3);
//	std::cout << "L2 error  = " << l2error << std::endl;
};

//V-cycle，一次VMG,n为最细的网格编号;
void NS2D::VMG(int n,int v1,int v2)
{
	double cost;
	int level=n;
	std::vector< Vector<double> > tem;//存储用于restrict的残差
	std::vector< Vector<double> > zero;//存储用于Interpolation的逼近
	tem.resize(mesh_l);
	zero.resize(mesh_l);
	for(int i=0;i<=level;i++)
	{
		tem[i].reinit(matrix[i]->n());
		zero[i].reinit(matrix[i]->n());
	}
	//1.在最细的网格Gi上smooth v1次
	GaussSidel(*matrix[level],u_h[level],rhs[level],v1);
	for (int i = 0;i < matrix[level]->n(); i ++)
		zero[level](i)=u_h[level](i);
	//restrict
	while(level!=0)
	{
		//将初值修正到右端项
		tem[level]=matrix_multi(*matrix[level],zero[level]);
		for (int i = 0;i < matrix[level]->n(); i ++)
			tem[level](i)=rhs[level](i)-tem[level](i);
		//2.将残差传递到粗一级的网格上Gi->Gi-1,ri-1=p*ri
		for(int i=0;i<tem[level-1].size();i++)
		{
			tem[level-1](i)=0.0;
			for(set<int>::iterator k=connect[level-1][i].begin();k!=connect[level-1][i].end();k++)
			{
				double p1=0.5;
				if((i==locate[level-1][*k])&&(i!=0))
					p1=1.0;
				if((i==0)&&(*k==0))
					p1=1.0;
				
				tem[level-1](i)+=p1*tem[level](*k);
			}
		}
		//存储此时的残差
		for (int i = 0;i < matrix[level-1]->n(); i ++)
				rhs[level-1](i)=tem[level-1](i);
		//3.取ui-1=0,在Gi-1上smooth v1次
		GaussSidel(*matrix[level-1],zero[level-1],rhs[level-1],v1);
		//存储逼近
		for (int i = 0;i < matrix[level-1]->n(); i ++)
				u_h[level-1](i)=zero[level-1](i);
		level--;
		//4.在最粗的网格G0上解出确解,因为矩阵很小,GS迭代到确值
		if(level==0)
		{
			GaussSidel(*matrix[0],zero[0],rhs[0],100);
			for (int i = 0;i < rhs[0].size(); i ++)
				u_h[0](i)=zero[0](i);
			break;
		}
	}
	//interpolation
	while(level!=n)
	{
		//5.将粗网格上的解扩展到细一级的网格
		for(int i=0;i<zero[level+1].size();i++)
		{
			zero[level+1](i)=0.0;
			//1.处理旧点，等价过去即可
			if(flag[level][i]==1)
				zero[level+1](i)=u_h[level](locate[level][i]);
			//2.处理新点，相邻粗点线性延拓
			if(flag[level][i]==0)
			{
				for(set<int>::iterator it=rconnect[level][i].begin();it!=rconnect[level][i].end();it++)
					zero[level+1](i)+=0.5*u_h[level](locate[level][*it]);
			}		
		}
		//6.更新细网格上的解
		for (int i = 0;i < matrix[level+1]->n(); i ++)
			zero[level+1](i)=u_h[level+1](i)+zero[level+1](i);
		//6.在细网格上smoothv2次
		GaussSidel(*matrix[level+1],zero[level+1],rhs[level+1],v2);
		for (int i = 0;i < matrix[level+1]->n(); i ++)
			u_h[level+1](i)=zero[level+1](i);
		level++;
	}
};

//基于AFEPack加密过程中的信息，建立粗网格与细网格间索引关系，如果是多重网格，相邻两层网格都要做一次；
void  NS2D::getindex(int n)
{
	//产生粗网格.
	if(n==0)
	{
	irregular_mesh_v[n]->semiregularize();
	irregular_mesh_v[n]->regularize(false);
	}

	RegularMesh<DIM> &mesh_p = irregular_mesh_v[n]->regularMesh();
	//产生对应细网格
	irregular_mesh_v[n+1]=new IrregularMesh<DIM> (*irregular_mesh_v[n]);
	//建立网格间点编号索引,2维和3维时粗网格和细网格的自由度比值不同,故n_nex_vtx=当前维度单位元个数＊其上的自由度个数不同；
	int n_of_point;
	if(DIM==3)
		n_of_point=32;
	else if(DIM==2)
		n_of_point=12;
	int n_nex_vtx =	mesh_p.n_geometry(DIM)*n_of_point;
	
	std::vector<int> idx_ir(n_nex_vtx);
	std::vector<int> idx_it(n_nex_vtx);
	
	//全局加密,生成细网格先不编号来确认粗网格位置
	irregular_mesh_v[n+1]->globalRefine(1);
	
	//基于AFEPack的几何数构建两个索引，编号前细网格（实际上内部的是粗网格信息）,编号后细网格到遍历element时顺序
	//以下三次基于单元顺序遍历几何树
	IrregularMesh<DIM, DIM>::ActiveIterator v_iterator = irregular_mesh_v[n+1]->beginActiveElement();
	IrregularMesh<DIM, DIM>::ActiveIterator v_end = irregular_mesh_v[n+1]->endActiveElement();
	for (int k = 0; v_iterator != v_end; ++v_iterator)
	{
		int n_vtx = v_iterator->h_element->n_vertex;
		for (int i = 0; i < n_vtx; ++i)
		{
			int idx = v_iterator->h_element->vertex[i]->index;
			idx_ir[k] = idx;
			k++;
		}
	}
	//对细网格编号
	irregular_mesh_v[n+1]->semiregularize();
	irregular_mesh_v[n+1]->regularize(false);
	
	int k;
	for (v_iterator = irregular_mesh_v[n+1]->beginActiveElement(), k = 0;
	v_iterator != v_end; ++v_iterator)
	{
		int n_vtx = v_iterator->h_element->n_vertex;
		for (int i = 0; i < n_vtx; ++i)
		{
			int idx = v_iterator->h_element->vertex[i]->index;
			idx_it[idx]=k;
			k++;
		}
	}
	
	//至此，建立了细网格点到粗网格点的索引关系，即细网格上点i->该点在粗网格上的编号idx_ir[idx_it[i]]，若为新产生的点，则为0(此处要注意存在编号为0的点在两层网格上编号都为0）
	//最后得到locate函数,n代表索引是网格编号n-1和网格n之间的对应关系;
	//locate[n][i]为细网格上序号为i的点在粗网格上的编号;
	//标记点所属的网格,n代表flag是网格编号n-1和网格n之间的区别;
	//flag[n][i]即细网格编号下粗细网格点的值，粗为1，细为0;
	int pointn=irregular_mesh_v[n+1]->regularMesh().n_geometry(0);
	locate[n].resize(pointn);
	flag[n].resize(pointn);
	for(int i=0;i<pointn;++i)
	{
		flag[n][i]=1;
		locate[n][i]=idx_ir[idx_it[i]];
		if((locate[n][i]==0)&&(i!=0))
			flag[n][i]=0;
	}
	
	rconnect[n].resize(pointn);
	connect[n].resize(irregular_mesh_v[n]->regularMesh().n_geometry(0));
	//再遍历一次细网格的几何树，建立rconnect
	//通过flag函数,当单位元里有粗网格点i时，将与其相邻的细网格点的编号存到connect[i]中,对细网格点i,将与其相邻的粗网格点的编号存到rconnect[i]中;
	//CMG下，只需要知道每个细网格点相邻的粗网格点编号，初值线性延拓过去即可；
	v_iterator = irregular_mesh_v[n+1]->beginActiveElement();
	for (int k = 0; v_iterator != v_end; ++v_iterator)
	{
		int n_vtx = v_iterator->h_element->n_vertex;
		for (int i = 0; i < n_vtx; ++i)
		{
			int idxi = v_iterator->h_element->vertex[i]->index;
			for (int j = 0; j < n_vtx; ++j)
			{
				int idxj = v_iterator->h_element->vertex[j]->index;
				if((flag[n][idxi]==0)&&(flag[n][idxj]==1))
					rconnect[n][idxi].insert(idxj);
					
				if(flag[n][idxi]==1)
					connect[n][locate[n][idxi]].insert(idxj);
			}
		}
	}
};

//构建矩阵A的系数矩阵模板
void NS2D::Construct_Matrix(int n_dof_p,int n)
{
	/// 准备统计系数矩阵的每一行有多少个非零元.
	std::vector<unsigned int> n_non_zero_per_row_p(n_dof_p);

	FEMSpace<double,DIM>::ElementIterator the_element_p= fem_space_v[n].beginElement();
	FEMSpace<double,DIM>::ElementIterator end_element_p= fem_space_v[n].endElement();

	//拼装粗网格模板
	//第一次遍历，统计非零元个数
	for (; the_element_p != end_element_p; ++the_element_p)
	{
		const std::vector<int>& element_dof_p = the_element_p->dof();
		int n_element_dof_p = the_element_p->n_dof();
		for (int j = 0; j < n_element_dof_p; ++j)
			for (int k = 0; k < n_element_dof_p; ++k)
				n_non_zero_per_row_p[element_dof_p[j]] += 1;
	}
	sp_stiff_matrix[n].reinit(n_dof_p, n_dof_p, n_non_zero_per_row_p, true);
	
	/// 第二次遍历, 指定每个非零元的坐标.
	for (the_element_p = fem_space_v[n].beginElement();
			the_element_p != end_element_p; ++the_element_p)
	{
		const std::vector<int>& element_dof_p = the_element_p->dof();
		int n_element_dof_p = the_element_p->n_dof();
		for (int j = 0; j < n_element_dof_p; ++j)
			for (int k = 0; k < n_element_dof_p; ++k)
				sp_stiff_matrix[n].add(element_dof_p[j], element_dof_p[k]);
	}
	sp_stiff_matrix[n].compress();
	
};

//矩阵A的拼装
void NS2D::Matrix_build(int n)
{
	//拼装当前网格上的矩阵
	FEMSpace<double,DIM>::ElementIterator the_element_v = fem_space_v[n].beginElement();
	FEMSpace<double,DIM>::ElementIterator end_element_v = fem_space_v[n].endElement();

	for (the_element_v = fem_space_v[n].beginElement();the_element_v != end_element_v; ++the_element_v)
	{
		double volume = the_element_v->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(0);
		std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<std::vector<double> > > basis_gradient = the_element_v->basis_function_gradient(q_point);
		std::vector<std::vector<double> >  basis_value = the_element_v->basis_function_value(q_point);
		const std::vector<int>& element_dof = the_element_v->dof();
		int n_element_dof = the_element_v->n_dof();
		for (int l = 0; l < n_quadrature_point; ++l)
		{
			double Jxw = quad_info.weight(l) * jacobian[l] * volume;
			for (int i = 0; i < n_element_dof; ++i)
			{
				for (int j = 0; j < n_element_dof; ++j)
				{
					double cont = Jxw * innerProduct(basis_gradient[i][l], basis_gradient[j][l]);
					matrix[n]->add(element_dof[i], element_dof[j], cont);
				}
				rhs[n](element_dof[i]) += Jxw * 1.0 * basis_value[i][l];
			}
		}
	}
	
	//边界条件处理
	//ex1.1 2D
//	BoundaryFunction<double,2> boundary(BoundaryConditionInfo::DIRICHLET, 1, &zero_bnd);
//	BoundaryConditionAdmin<double,2> boundary_admin(fem_space_v[n]);
//        boundary_admin.add(boundary);
//        boundary_admin.apply(*matrix[n], u_h[n], rhs[n]);
        
    //维下边界条件设定
	const std::size_t * rowstart = sp_stiff_matrix[n].get_rowstart_indices();
	const unsigned int * colnum = sp_stiff_matrix[n].get_column_numbers();
	FEMFunction<double, DIM> v_h_a;
	v_h_a.reinit(fem_space_v[n]);
	/// 遍历全部自由度。
	for (unsigned int i = 0; i < matrix[n]->n(); ++i)
	{
		unsigned int bm = fem_space_v[n].dofInfo(i).boundary_mark;
		/// 若该自由度标识为 1，这个是在生成网格时手工指定的。
		//std::cout << "build." <<i<< std::endl;
			if (bm == 1)
			{
				/// 计算并代入真解，因为这里是 Dirichlet 条件。
				v_h_a(i) = u(fem_space_v[n].dofInfo(i).interp_point);
				/// 对应右端项设为真解乘以对角元。这样可以避免过度修改矩阵条件。 
				rhs[n](i) = stiff_matrix[n].diag_element(i) * v_h_a(i); 
				/// 遍历非对角非零元.
				for (unsigned int j = rowstart[i] + 1; j < rowstart[i + 1]; ++j) 
				{ 
					/// 全部按成 0 .
					matrix[n]->global_entry(j) -= stiff_matrix[n].global_entry(j);
					/// 记下列号.
					unsigned int k = colnum[j];
					/// 在第 k 行寻找第 n_total_dof_v + i 列. (不可能是对角元) .
					const unsigned int *p = std::find(&colnum[rowstart[k] + 1], 
					&colnum[rowstart[k + 1]], i);
					/// 如果是非零元.
					if (p != &colnum[rowstart[k + 1]]) 
					{
						unsigned int l = p - &colnum[rowstart[0]];
						/// 把对应项移到右端.
						rhs[n](k) -= stiff_matrix[n].global_entry(l) 
											* rhs[n](i) 
											/ stiff_matrix[n].diag_element(i);
						matrix[n]->global_entry(l) -= stiff_matrix[n].global_entry(l);
					}
				}
			}
	}
};

//G_S 迭代
void NS2D::GaussSidel(const SparseMatrix<double> & M,Vector<double>& x,const Vector<double>& r,const int& s)const
{
	const SparsityPattern& spM = M.get_sparsity_pattern();
	const std::size_t * rowstart = spM.get_rowstart_indices();
	const unsigned int * colnums = spM.get_column_numbers();
	for (int i = 0;i < s;i ++)//iteration times
	{
		for (int j = 0;j < M.m();j ++)
		{
			double r0 = r(j);
			for (int k = rowstart[j] + 1;k < rowstart[j+1];k ++)
				r0 -= x(colnums[k]) * M.global_entry(k);
			x(j) = r0 / M.global_entry(rowstart[j]);
		}
	}
}

//矩阵＊向量
Vector<double> NS2D::matrix_multi(const SparseMatrix<double>& P,Vector<double>& rhs)
{
	const SparsityPattern& spP = P.get_sparsity_pattern();
	const std::size_t * P_rowstart = spP.get_rowstart_indices();
	const unsigned int * P_colnum = spP.get_column_numbers();
	Vector<double> A(P.m());

	for (int i = 0;i < P.m();i ++)
	{
		for (int j = P_rowstart[i];j < P_rowstart[i+1];j ++)
		{
			const int& a = P_colnum[j];
			A(i)+=P.el(i,a)*rhs(a);
		}
	}
	return A;
}

void NS2D::outputSolution()
{
	std::stringstream result;
	result.setf(std::ios::fixed);
	result.precision(4);
	result << "v_h_" << int((t - 0.25) / dt) << ".dx";
	u_h[mesh_l-1].writeOpenDXData(result.str());
	result.str(std::string());
};

void NS2D::run()
{
	initialize();
	solve();
	outputSolution();
};

void NS2D::initialValue()
{
};

void NS2D::boundaryValue()
{
};
void NS2D::stepForward()
{
};

//
// end of file
///////////////////////////////////////////////////////////////////////////////
