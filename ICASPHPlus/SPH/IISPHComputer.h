#pragma once
#include "HashGridList.h"
#include "GridCell.h"
#include "Vector3f.h"
#include "Vector3i.h"
#include "ShareData.h"
#include "FunctionKit.h"
#include "SPHKernel.h"
#include "WCSPHFluidObject.h"
#include "IISPHParticle.h"
#include "IISPHFluidObject.h"
#include "Boundary.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
using namespace std;


class IISPHComputer
{
public:
	HashGridList hashGridList;       //哈希网格表
	IISPHFluidObject fluidModel;     //流体表面模型
	Boundary  boundary;              //流体的边界模型
	double viscosityConstant = 0.03; //计算黏度所需要的常量
	double gravity = -9.81;          //流体的重力
	//关于迭代
	int m_iterations = 0;
	double m_maxError = 0.01;
	int m_maxIterations = 100;

	//////////////////////////////////////////////////////////////////////////
	// 无限自适应扩展
	//////////////////////////////////////////////////////////////////////////
	double baseDistance;       //最远的距离
	double fineDistance;//最近距离
	double baseMass = 0.2f;           //基础的质量
	double fineMass = 0.025f;           //最小的质量
	double adaptivityFactor = fineMass/baseMass;   //自适应因子



public:
	//初始化
	void Initialization(); //在SPH计算之前调用一次，对粒子的属性、位置进行初始化

	//IISPH的主要步骤函数：以process函数开头的阶段计算函数
	//计算部分的阶段函数函数
	void Computation();  //SPH计算，该函数封装下列方法
	//1.计算临近粒子列表
	void ProcessNeighborList();
	//2.处理密度计算
	void ProcessDensity();
	//3.处理非压力的力
	void ProcessViscosity();
	//4.更新时间步骤长
	void ProcessUpdateTimeStep();
	//5.预测平流
	void ProcessPredictAdvection();
	//6.压力解算
	void ProcessPressureSolve();
	//7.积分计算
	void ProcessIntegration();

	//IISPH的计算细节函数，下面的函数是为了配合上面的阶段函数进行计算的
	//重置当前所有粒子的加速度
	void ClearAcceleration();
    //计算粒子周围邻居粒子
	vector<IISPHParticle*> ComputeNeighborParticle(IISPHParticle* origin);
	//计算粒子周围的边界粒子：用于补正粒子密度
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(IISPHParticle* origin);
	//计算边界粒子周围的边界粒子：用于进行边界粒子PSI计算
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(StaticRigidParticle* origin);
	//根据粒子周围的邻居，计算粒子的密度，并且所有粒子的压强已经计算成功
	double ComputeDensity(vector<IISPHParticle>& neighborList, vector<StaticRigidParticle>& boundaryNeighborList, IISPHParticle& origin);
    //根据粒子周围的邻居，计算粒子在粘性力下的加速度
	Vector3f ComputeViscostyAcceleration(vector<IISPHParticle>& neighborList, IISPHParticle& origin);
	//在计算好压力的情况下，计算粒子的压力加速度
	void ComputePressureAcceleration();
	//计算流体的邻居表
	void UpdateFluidNeighborList();
	//计算流体的边界邻居表，并且需要不断更新
	void UpdateFluidBoundaryNeighborList();
	//计算边界的边界邻居表，只需要初始化调用一次
	void ComputeBoundaryNeighborList();
	//根据粒子的位置进行网格分类
	void MapParticleToGrid(); 
	//将边界粒子映射到网格当中
	void MapBoundaryParticleToGrid(); 
	//基于压力计算得到每个粒子的压力加速度
	void ComputePressureAccels();
	//论文中的方法
	void ComputePressureAccels2();
	
	//////////////////////////////////////////////////////////////////////////
	// 无限自适应扩展
	//////////////////////////////////////////////////////////////////////////
	//计算自适应因子
	void ComputeAdaptivityFactor();
	//计算所有粒子到表面的距离
	void ComputeDistance();
	//根据粒子的mark来对粒子进行分类
	void ClassifyParticle();
	//计算所有粒子的最优粒子质量
	void ComputeOptMass();
	//依据粒子的印记 计算所有粒子的O
	void ComputeO();
	//隐性融合密度
	void BlendDensity();
	//隐性融合速度
	void BlendVelocity();
	//更新所有粒子的融合因子
	void UpdateBlendFactor();

	//add
	//计算每个粒子的距离域DistanceField
	void ComputePtoSandMopt();
	void ProcessMergeandSplit();
};