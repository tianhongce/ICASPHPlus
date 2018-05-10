#pragma once
#include "Vector3f.h"
#include "IISPHParticle.h"
#include "StaticRigidObject.h"
#include "ParticleObject.h"
#include "Boundary.h"
#include "FunctionKit.h"
#include "SPHKernel.h"
#include "ShareData.h"
#include <vector>



class IISPHFluidObject
{
public:
	//时间步骤的初始信息
	double IISPH_TimeStep = 0.001f;           //只有在流体模拟的时候，使用时间步骤来进行
										      //FluidObject将粒子的固定属性提取出来，减少每个粒子的对象大小
						
	//粒子的固定属性
	double IISPH_RestDensity = 1000.f;        //流体的静止密度，总会用到的参数
	int    IISPH_ParticleNum;                 //粒子数量需要计算出来

	//水体渲染的初始信息
	double IISPH_FluidWitdth = 1.0f;
	double IISPH_FluidLength = 1.0f;
	double IISPH_FluidHeight = 1.0f;          //水体的大小
	double IISPH_OriginX = 0.25f;
	double IISPH_OriginY = 0.25f;
	double IISPH_OriginZ = 0.075f;              //水体的起始渲染位置									

	//边界相关信息
	Boundary boundary;                                             //设置边界的属性：用于计算边界粒子的位置初始化
	RigidBodyParticleObject boundaryObj;                           //粒子组装成的静态刚体边界对象

	//临近表
	vector<IISPHParticle*> particleList;                             //当前水体的粒子列表
	//vector<vector<IISPHParticle>> neighborList;                     //水体的每个粒子的邻居表：用于流体内部进行属性计算
	vector<vector<StaticRigidParticle*>> boundaryNeighborList;       //边界粒子的邻居：用于边界PSI计算
	//vector<vector<StaticRigidParticle>> fluidBoundaryNeighborList;  //流体的边界粒子邻居：用于流体密度补正

	//水体初始化的信息
	Vector3i FluidWLH;          //这个是水体初始化的时候，长宽高方向上有多少个粒子,只在初始化时候有用
	double particleInitialPad = 0.075f;  //水体粒子的初始化半径，由于之后粒子的大小自适应，因此这个值只在初始化的时候有用

	//流体表中含有流体粒子的个数
	int numParticleObj;


public:
	//对当前水体进行初始化
	void Initialise();                //根据默认属性进行初始化

	void ComputeFluidWLH();
	void InitialiseParticlePosition();//初始化水体的粒子位置
	void InitialiseParticleDensity(); //初始化水体的密度

	void UpdateParticleHashValue();   //更新每个粒子的哈希值
	void ComputeParticleNum(int x, int y, int z);
	void ClearNeighborList();
	//清除流体粒子周围的边界粒子邻居
	void ClearFluidBoundaryNeighborList();

	//////////////////////////////////////////////////////////////////////////
	//***************粒子操作***********************************************//
	//////////////////////////////////////////////////////////////////////////
	//对所有边界粒子进行初始化
	void InitializeBoundary();
	//通过刚体计算从而完成所有的边界粒子的更新
	void AddRigidBodyObject(StaticRigidObject& rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles);
	void ComputeBoundaryPsi(RigidBodyParticleObject& bound);
	//计算流体中粒子对象的个数
	void ComputeNumParticleObj();
};