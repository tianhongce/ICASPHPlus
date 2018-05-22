#include "IISPHComputer.h"

void IISPHComputer::Initialization()
{
	//对哈希网格以及流体模型都使用默认参数进行初始化
	hashGridList.boundary = boundary;
	hashGridList.Intialization();
	fluidModel.boundary = boundary;
	fluidModel.Initialise();
	//将所有的粒子都装入对应的网格当中
	MapParticleToGrid();
	//将所有的边界粒子装入对应的网格当中
	MapBoundaryParticleToGrid();
	//根据网格当中的邻居边界粒子，计算得到边界粒子的边界粒子邻居列表
	ComputeBoundaryNeighborList();
	//计算边界粒子的PSI值
	fluidModel.ComputeBoundaryPsi(fluidModel.boundaryObj);
	//需要再一次将边界粒子装入网格当中，来更新BoundaryPSI值
	MapBoundaryParticleToGrid();
	//在所有的流体粒子以及边界粒子完成网格装填之后
	//更新当前流体粒子的流体粒子邻居表
	UpdateFluidNeighborList();
	//更新当前流体粒自的边界粒子邻居表
	UpdateFluidBoundaryNeighborList();
	//对SPHComputation当中的参数装入全局变量当中
	SharedData::SetGravity(gravity);
	SharedData::SetViscosityConstant(viscosityConstant);
}

void IISPHComputer::Computation()
{
	cout << "time_frame" << endl;

	//将所有的粒子的加速度初始化为重力加速度
	//粒子的β值每一步的递减放在函数ClearAcceleration()中
	ClearAcceleration();
	//将所有粒子的邻居表进行更新
	ProcessNeighborList();
	//计算所有粒子的密度
	ProcessDensity();
	//计算所有粒子的黏度
	ProcessViscosity();
	//使用CFL条件约束，更新TimeStep（待定）
	ProcessUpdateTimeStep();
	//预测平流
	ProcessPredictAdvection();
	//压力计算
	ProcessPressureSolve();	
	//计算所有的位置
	ProcessIntegration();
	//完成一次计算，根据粒子的新位置，计算得到粒子的新哈希值
	fluidModel.UpdateParticleHashValue();
	//将所有粒子的位置都更新到网格当中
	MapParticleToGrid();
	//完成一次IISPH计算
	//自适应分裂
	ComputePtoSandMopt();

	//ProcessMergeandSplit();

	fluidModel.ComputeNumParticleObj();
}

void IISPHComputer::ProcessNeighborList()
{
	cout << "ProcessNeighborList()" << endl;
	//更新流体的流体粒子列表
	UpdateFluidNeighborList();
	//更新流体的边界粒子列表0-+
	UpdateFluidBoundaryNeighborList();
}

void IISPHComputer::ProcessDensity()
{
	cout << "ProcessDensity()" << endl;
	//得到当前粒子个数
	const int numParticles = fluidModel.particleList.size();
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			double &density = fluidModel.particleList[i]->density;
			double massi = fluidModel.particleList[i]->particleMass;
			double supportRad = fluidModel.particleList[i]->particleSupportRad;
			// Compute current density for particle i
			density = massi * Poly6KernelOne::W_zero(supportRad);
			Vector3f &xi = fluidModel.particleList[i]->position;

			vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
			vector<StaticRigidParticle*> boundaryNeighbors = fluidModel.particleList[i]->boundaryNeighbors;
			//////////////////////////////////////////////////////////////////////////
			// 流体计算
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
			{
				if (fluidNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = fluidNeighbors[j]->position;
				double massj = fluidNeighbors[j]->particleMass;
				density += massj * Poly6KernelOne::W(supportRad,xi - xj);
			}
			double test1 = density;


			//////////////////////////////////////////////////////////////////////////
			// 边界计算
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < boundaryNeighbors.size(); j++)
			{
				Vector3f &xj = boundaryNeighbors[j]->position;
				density += boundaryNeighbors[j]->boundaryPsi* Poly6KernelOne::W(supportRad, xi - xj);
			}

			//得到当前粒子的混合因子
			const double beita = fluidModel.particleList[i]->β;
			//得到当前粒子的lastDensity
			double& density0 = fluidModel.particleList[i]->lastDensity;
			//ρblended i = (1 − βi)*ρi + βi* ρO 刷新当前的密度
			//cout << i << " : density0: " << density0 << endl;
			//cout << "beita: " << beita <<" : "<<i<<" : density: "<<density<< endl;
			density = (1 - beita)*density + beita*density0;
			//cout << i << " : density: " << density << endl;
			density0 = density;
			//cout << fluidModel.particleList[i]->β << endl;

			//if (i == 75)
			//{
			//	cout << "Density Two: " << density << endl;
			//}
			//cout<<"Density Two: "<<density<<endl;
		}
	}

	////得到当前粒子个数
	//const unsigned int numParticles = fluidModel.particleList.size();
 //   #pragma omp parallel default(shared)
	//{
 //       #pragma omp for schedule(static)  
	//	for (int i = 0; i < (int)numParticles; i++)
	//	{
	//		//计算得到流体的邻居
	//		vector<IISPHParticle> fluidNeighbor = fluidModel.neighborList[i];
	//		//计算得到流体的边界邻居
	//		vector<StaticRigidParticle> boundaryNeighbor = fluidModel.fluidBoundaryNeighborList[i];
	//		//计算当前粒子的密度，并且对粒子的密度进行赋值
	//		fluidModel.particleList[i].density = ComputeDensity(fluidNeighbor, boundaryNeighbor, fluidModel.particleList[i]);
	//		//cout << fluidModel.particleList[i].density << endl;
	//	}
	//}
}

void IISPHComputer::ProcessViscosity()
{
	cout << "ProcessViscosity()" << endl;
	const unsigned int numParticles = fluidModel.particleList.size();
	const double h = fluidModel.IISPH_TimeStep;
	const double invH = (1.0 / h);

	// Compute viscosity forces (XSPH)
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			Vector3f &xi = fluidModel.particleList[i]->position;
			Vector3f &vi = fluidModel.particleList[i]->velocity;
			//直接对加速度进行赋值
			Vector3f &ai = fluidModel.particleList[i]->acceleration;
			const double density_i = fluidModel.particleList[i]->density;
			const double supportRad = fluidModel.particleList[i]->particleSupportRad;
			vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
			{
				if (fluidNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = fluidNeighbors[j]->position;
				Vector3f &vj = fluidNeighbors[j]->velocity;

				// Viscosity
				const double density_j = fluidNeighbors[j]->density;
				double massj = fluidNeighbors[j]->particleMass;
				ai -= invH * viscosityConstant * ( massj / density_j) * (vi - vj) * CubicKernelOne::W(supportRad, xi - xj);
			}
		}
	}
	//// 使用XSPH的方法对粒子黏度进行计算
 //   #pragma omp parallel default(shared)
	//{
 //       #pragma omp for schedule(static)  
	//	for (int i = 0; i < (int)fluidModel.IISPH_ParticleNum; i++)
	//	{
	//		vector<IISPHParticle> fluidNeibors = fluidModel.neighborList[i];
	//		fluidModel.particleList[i].acceleration	 = ComputeViscostyAcceleration(fluidNeibors, fluidModel.particleList[i]);
	//	}
	//}
}

//处理时间步长
void IISPHComputer::ProcessUpdateTimeStep()
{
	cout << "ProcessUpdateTimeStep()......待定	" << endl;

}

void IISPHComputer::ProcessIntegration()
{
	cout << "ProcessIntegration()" << endl;
	int numParticles = fluidModel.particleList.size();
	//计算当前的所有粒子的压力加速度
	//ComputePressureAcceleration();
	ComputeO();
	ComputePressureAccels2();

	//为每一个粒子计算其在时间步骤当中的位移以及速度变化
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			fluidModel.particleList[i]->Intigration(boundary, fluidModel.IISPH_TimeStep);
		}
	}
}

void IISPHComputer::ClearAcceleration()
{
	cout << "ClearAcceleration()" << endl;
	//将所有的粒子的加速度都初始化为G 
	const unsigned int count = fluidModel.particleList.size();
	const Vector3f G(0.0, 0.0, gravity);
	for (unsigned int i = 0; i < count; i++)
	{
		if (fluidModel.particleList[i] == NULL)
		{
			continue;
		}
		fluidModel.particleList[i]->acceleration = G;
		//cout << "fluidModel.particleList[i]->β : " << i <<" " << fluidModel.particleList[i]->β << endl;
		//此处添加粒子分裂、合并中的 β
		if (fluidModel.particleList[i]->β > 0)
		{
			fluidModel.particleList[i]->β = fluidModel.particleList[i]->β - 0.1;
		}
		
	}
}

vector<IISPHParticle*> IISPHComputer::ComputeNeighborParticle(IISPHParticle* origin)
{
	double h = SharedData::GetCellLength(); //得到粒子的支撑半径，也就是网格单元的长度
	Vector3f originPosition = origin->position;
	//得到角落里的粒子表格
	Vector3f tempPosition; //用于对27个网格进行遍历的工具
	//int gridCellOffserNum = origin->gridCellNum;
	int gridCellOffserNum = 1;
	gridCellOffserNum = origin->gridCellNum;
	//cout << "origin->gridCellNum: " << origin->gridCellNum << endl;
	//cout << "gridCellOffserNum: : " << gridCellOffserNum << endl;
	if (gridCellOffserNum == 8101)
	{
		cout << "aaa" << endl;
	}
	tempPosition.SetX(originPosition.GetX() - h*gridCellOffserNum);
	tempPosition.SetY(originPosition.GetY() - h*gridCellOffserNum);
	tempPosition.SetZ(originPosition.GetZ() - h*gridCellOffserNum);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	int offsetXNum = gridCellOffserNum * 2 + 1;
	int offsetYNum = gridCellOffserNum * 2 + 1;
	int offsetZNum = gridCellOffserNum * 2 + 1;

	vector<IISPHParticle*> neighborList; //粒子的邻居粒子列表
	int cellNum = hashGridList.GetCellNum(); //得到当前网格的数量
	for (int i = 0; i < offsetXNum; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < offsetYNum; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < offsetZNum; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				vector<IISPHParticle*> tempParticleList = targetGridCell.GetIISPHParticleList();
				//将目标粒子列表当中的粒子一个一个的塞进邻居列表
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//如果在邻域搜索当中出现自己，那么就忽略掉
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//只压入在h范围之内（在粒子的影响范围之内）的邻居粒子
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < origin->particleSupportRad)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//其他的舍弃
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

vector<StaticRigidParticle*> IISPHComputer::ComputeNeighborBoundaryParticle(IISPHParticle* origin)
{
	double h = SharedData::GetCellLength(); //得到粒子的支撑半径，也就是网格单元的长度
	Vector3f originPosition = origin->position;
	//得到角落里的粒子表格
	Vector3f tempPosition; //用于对27个网格进行遍历的工具
	int gridCellNum = 1;
	gridCellNum = origin->gridCellNum;
	tempPosition.SetX(originPosition.GetX() - h*gridCellNum);
	tempPosition.SetY(originPosition.GetY() - h*gridCellNum);
	tempPosition.SetZ(originPosition.GetZ() - h*gridCellNum);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	int offsetXNum = gridCellNum * 2 + 1;
	int offsetYNum = gridCellNum * 2 + 1;
	int offsetZNum = gridCellNum * 2 + 1;
	vector<StaticRigidParticle*> neighborList; //粒子的邻居粒子列表
	int cellNum = hashGridList.GetCellNum(); //得到当前网格的数量

	for (int i = 0; i < offsetXNum; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < offsetYNum; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < offsetZNum; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				//运行笔记：由于找不到对应的网格，返回的网格是错误的，所以，这个地方的代码是错误的
				//解决原因：其中一个函数调用参数用的是网格数量而不是哈希表格的数量
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				//根据目标粒子的位置，直接查找到对应的哈希列表当中的网格，并且获取其中的粒子列表
				//运行笔记：找到对应网格之后，无法找到网格当中的粒子
				//解决方式：粒子映射网格拷贝错误，导致无法找到网格当中的粒子
				vector<StaticRigidParticle*> tempParticleList = targetGridCell.GetBoundaryParticleList();
				//将目标粒子列表当中的粒子一个一个的塞进邻居列表
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//如果在邻域搜索当中出现自己，那么就忽略掉
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//只压入在h范围之内（在粒子的影响范围之内）的邻居粒子
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < origin->particleSupportRad)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//其他的舍弃
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

vector<StaticRigidParticle*> IISPHComputer::ComputeNeighborBoundaryParticle(StaticRigidParticle* origin)
{
	double h = SharedData::GetCellLength(); //得到粒子的支撑半径，也就是网格单元的长度
	Vector3f originPosition = origin->position;
	//得到角落里的粒子表格
	Vector3f tempPosition; //用于对27个网格进行遍历的工具
	tempPosition.SetX(originPosition.GetX() - h);
	tempPosition.SetY(originPosition.GetY() - h);
	tempPosition.SetZ(originPosition.GetZ() - h);
	double originX = tempPosition.GetX();
	double originY = tempPosition.GetY();
	double originZ = tempPosition.GetZ();
	vector<StaticRigidParticle*> neighborList; //粒子的邻居粒子列表
	int cellNum = hashGridList.GetCellNum(); //得到当前网格的数量

	for (int i = 0; i < 3; i++, tempPosition.SetX(tempPosition.GetX() + h))
	{
		for (int j = 0; j < 3; j++, tempPosition.SetY(tempPosition.GetY() + h))
		{
			for (int k = 0; k < 3; k++, tempPosition.SetZ(tempPosition.GetZ() + h))
			{
				//运行笔记：由于找不到对应的网格，返回的网格是错误的，所以，这个地方的代码是错误的
				//解决原因：其中一个函数调用参数用的是网格数量而不是哈希表格的数量
				GridCell targetGridCell = hashGridList.GetGridCell(tempPosition);
				//根据目标粒子的位置，直接查找到对应的哈希列表当中的网格，并且获取其中的粒子列表
				//运行笔记：找到对应网格之后，无法找到网格当中的粒子
				//解决方式：粒子映射网格拷贝错误，导致无法找到网格当中的粒子
				vector<StaticRigidParticle*> tempParticleList = targetGridCell.GetBoundaryParticleList();
				//将目标粒子列表当中的粒子一个一个的塞进邻居列表
				for (int q = 0; q < tempParticleList.size(); q++)
				{
					//如果在邻域搜索当中出现自己，那么就忽略掉
					if (origin->position == tempParticleList[q]->position)
					{
						continue;
					}
					//只压入在h范围之内（在粒子的影响范围之内）的邻居粒子
					if (FunctionKit::DistanceComputation(originPosition,
						tempParticleList[q]->position) < h)
					{
						neighborList.push_back(tempParticleList[q]);
					}
					//其他的舍弃
				}
			}
			tempPosition.SetZ(originZ);
		}
		tempPosition.SetY(originY);
	}
	return neighborList;
}

double IISPHComputer::ComputeDensity(vector<IISPHParticle>& neighborList, vector<StaticRigidParticle>& boundaryNeighborList, IISPHParticle & origin)
{
	int i;
	double mass = origin.particleMass;
	double result = CubicKernelOne::W_zero(origin.particleSupportRad)*mass;
	double h = hashGridList.GetCellLength(); //得到当前单元的长度
	Vector3f xi = origin.position;
	//得到粒子的支撑半径
	double supportRad = origin.particleSupportRad;

	//计算流体粒子对流体粒子之间的影响
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (i = 0; i < neighborList.size(); i++)
		{
			Vector3f xj = neighborList[i].position;
			double massj = neighborList[i].particleMass;
			//在核函数的范围之内的话
			Vector3f temp123 = xi - xj;
			double temp = CubicKernelOne::W(supportRad, xi - xj);
			result += CubicKernelOne::W(supportRad, xi - xj)*massj;
			//result += SPHKernal::Poly6Kernel(h, distance);
		}
	}
	double test1 = result;
	//计算边界粒子对流体粒子的影响，计算密度补正
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int j = 0; j < boundaryNeighborList.size(); j++)
		{
			Vector3f xj = boundaryNeighborList[j].position;
			double kernelResult = CubicKernelOne::W(supportRad, xi - xj);
			double boundPSI = boundaryNeighborList[j].boundaryPsi;
			//调用得到对应的边界邻居粒子
			result += boundaryNeighborList[j].boundaryPsi * CubicKernelOne::W(supportRad, xi - xj);
		}
	}
	
	//得到当前粒子的混合因子
	const double beita = origin.β;
	//得到当前粒子的lastDensity
	double& density0 = origin.lastDensity;
	//ρblended i = (1 − βi)*ρi + βi* ρO 刷新当前的密度
	result = (1 - beita)*result + beita*density0;
	density0 = result;
	//cout << "Density One：" << test1 << "Density Two: " << result << endl;
	return result;
}

Vector3f IISPHComputer::ComputeViscostyAcceleration(vector<IISPHParticle>& neighborList,  IISPHParticle & origin)
{
	double viscosityCons = viscosityConstant;
	//此处调用的是XSPH当中的Viscosity的计算方法
	Vector3f vi = origin.velocity;
	Vector3f xi = origin.position;
	double mi = origin.particleMass;
	double hi = origin.particleSupportRad;
	double di = origin.density;
	Vector3f result = origin.acceleration;
	double invH = 1.0 / SharedData::GetTimeStep();
	double supportRad = origin.particleSupportRad;
	//////////////////////////////////////////////////////////////////////////
	///流体自身的黏度计算
	//////////////////////////////////////////////////////////////////////////
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (int j = 0; j < neighborList.size(); j++)
		{
			Vector3f vj = neighborList[j].velocity;
			Vector3f xj = neighborList[j].position;
			double dj = neighborList[j].density;
			double mj = neighborList[j].particleMass;
			result -= invH*viscosityCons*(mj / dj)*(vi - vj)*CubicKernelOne::W(supportRad,xi - xj);
		}
	}
	return result;
}

void IISPHComputer::ComputePressureAcceleration()
{
	int numParticles = fluidModel.particleList.size();
	//计算压力
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			//对临近列表进行设置
			vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
			vector<StaticRigidParticle*> boundaryNeighbors = fluidModel.particleList[i]->boundaryNeighbors;

			//初始化位置
			Vector3f xi = fluidModel.particleList[i]->position;
			double density_i = fluidModel.particleList[i]->density;
			Vector3f ai = fluidModel.particleList[i]->pressureAcceleration;
			ai.Zero();
			const double dpi = fluidModel.particleList[i]->pressure / (density_i*density_i);
			//当前粒子的支撑半径
			double supportRad = fluidModel.particleList[i]->particleSupportRad;
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (int j = 0; j < fluidNeighbors.size(); j++)
			{
				Vector3f xj = fluidNeighbors[j]->position;
				// 对压力进行计算
				const double density_j = fluidNeighbors[j]->density;
				const double dpj = fluidNeighbors[j]->pressure / (density_j*density_j);
				ai -= fluidNeighbors[j]->particleMass * (dpi + dpj) * CubicKernelOne::gradW(supportRad, xi - xj);
			}
			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			for (int j = 0; j < boundaryNeighbors.size(); j++)
			{
				Vector3f xj = boundaryNeighbors[j]->position;
				Vector3f a = boundaryNeighbors[j]->boundaryPsi * (dpi)* CubicKernelOne::gradW(supportRad, xi - xj);
				ai -= a;
			}
			fluidModel.particleList[i]->pressureAcceleration = ai;
			//if (i == 75)
			//{
			//	cout << fluidModel.particleList[i]->pressureAcceleration << endl;
			//}
		}
	}
}

void IISPHComputer::UpdateFluidNeighborList()
{//在刷新新的粒子列表的时候，清除当前的邻居列表
	fluidModel.ClearNeighborList();
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		//遍历当前所有粒子的列表
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			//计算得到所有的临近的粒子
			vector<IISPHParticle*> temp = ComputeNeighborParticle(fluidModel.particleList[i]);
			//fluidModel.particleList[i]->fluidNeighbors =temp;
			for (int j = 0; j < temp.size(); j++)
			{
				if (temp[j] == NULL)
				{
					continue;
				}
				//将临近的所有粒子装入对应的粒子邻居表当中
				fluidModel.particleList[i]->fluidNeighbors.push_back(temp[j]);
			}
		}
	}
}

void IISPHComputer::UpdateFluidBoundaryNeighborList()
{
	//清理流体周围的边界粒子
	fluidModel.ClearFluidBoundaryNeighborList();
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			//得到对应的流体粒子中的临近静态粒子列表
			vector<StaticRigidParticle*> temp = ComputeNeighborBoundaryParticle(fluidModel.particleList[i]);
			for (int j = 0; j < temp.size(); j++)
			{
				if (temp[j] == NULL)
				{
					continue;
				}
				//将对应的静态流体粒子装入到流体周围的静态粒子列表当中
				fluidModel.particleList[i]->boundaryNeighbors.push_back(temp[j]);
			}
		}
	}
}

//初始化的时候进行调用
void IISPHComputer::ComputeBoundaryNeighborList()
{
	for (int i = 0; i < fluidModel.boundaryObj.staticRigidParticleList.size(); i++)
	{
		StaticRigidParticle* test2 = fluidModel.boundaryObj.staticRigidParticleList[i];
		//得到对应的流体粒子中的临近静态粒子列表
		vector<StaticRigidParticle*> temp = ComputeNeighborBoundaryParticle(fluidModel.boundaryObj.staticRigidParticleList[i]);
		for (int j = 0; j < temp.size(); j++)
		{
			//将对应的静态流体粒子装入到静态粒子周围的静态粒子列表当中
			fluidModel.boundaryNeighborList[i].push_back(temp[j]);
		}
	}
}

void IISPHComputer::MapParticleToGrid()
{
	cout << "MapParticleToGrid()" << endl;
	//清除网格当中的粒子之后再向其中添加新的粒子
	//将流体粒子装入对应位置
	static int a = 0;
	hashGridList.ClearParticle();
	for (int i = 0; i < fluidModel.particleList.size(); i++)
	{
		if (fluidModel.particleList[i] == NULL)
		{
			continue;
		}
		//输出每个粒子的信息
		//{
		//	cout << "i ：" << i << endl;
		//	cout << "Times :" << a << endl;
		//	cout << "Position: " << fluidModel.particleList[i]->position << endl;
		//	cout << "Pressure: " << fluidModel.particleList[i]->pressure << endl;
		//	cout << "Pressure Acceleration: " << fluidModel.particleList[i]->pressureAcceleration << endl;
		//	cout << "Velocity: " << fluidModel.particleList[i]->velocity << endl;
		//}

		//如果超过当前的值，那么输出当前问题粒子的诸多问题
		if (fluidModel.particleList[i]->position.x < -0.6f || fluidModel.particleList[i]->position.y < -0.6f || fluidModel.particleList[i]->position.z < -0.6f)
		//if(i==0)
		{
			a++;
			cout << "i ：" << i << endl;
			cout << "Times :" << a << endl;
			cout << "Position: " << fluidModel.particleList[i]->position << endl;
			cout << "Pressure: " << fluidModel.particleList[i]->pressure << endl;
			cout << "Pressure Acceleration: " << fluidModel.particleList[i]->pressureAcceleration << endl;
			cout << "Velocity: " << fluidModel.particleList[i]->velocity << endl;
		}

		hashGridList.PushParticle(fluidModel.particleList[i]);
	
		
	}
}

void IISPHComputer::MapBoundaryParticleToGrid()
{
	for (int i = 0; i < fluidModel.boundaryObj.GetStaticParticleNum(); i++)
	{
		StaticRigidParticle* tempStaticParticle = fluidModel.boundaryObj.staticRigidParticleList[i];
		hashGridList.PushParticle(tempStaticParticle);
	}
	//把所有的边界粒子存入边界网格当中
}

void IISPHComputer::ProcessPredictAdvection()
{
	cout << "ProcessPredictAdvection()" << endl;
	const unsigned int numParticles = fluidModel.particleList.size();
	const double h = fluidModel.IISPH_TimeStep;
	//预测对应粒子的属性 v_adv
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			//这里增加速度的值是粘性力施加的速度值
			Vector3f &vel = fluidModel.particleList[i]->velocity;
			const Vector3f &accel = fluidModel.particleList[i]->acceleration;
			Vector3f &dii = fluidModel.particleList[i]->dii;

			vel += h * accel;

			// Compute d_ii
			dii.Zero();
			double density = fluidModel.particleList[i]->density;
			const double density2 = density*density;
			Vector3f &xi = fluidModel.particleList[i]->position;
			double supportRad = fluidModel.particleList[i]->particleSupportRad;

			vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
			vector<StaticRigidParticle*> boundaryNeighbors = fluidModel.particleList[i]->boundaryNeighbors;
			//////////////////////////////////////////////////////////////////////////
			// 流体部分的计算
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
			{
				if (fluidNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = fluidNeighbors[j]->position;
				dii -= fluidNeighbors[j]->particleMass / density2 * CubicKernelOne::gradW(supportRad, xi - xj);
			}

			//////////////////////////////////////////////////////////////////////////
			// 边界粒子的计算
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < boundaryNeighbors.size(); j++)
			{
				if (boundaryNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = boundaryNeighbors[j]->position;
				dii -= boundaryNeighbors[j]->boundaryPsi / density2 * CubicKernelOne::gradW(supportRad, xi - xj);
			}
			//cout << "DII:   " << dii(0,0)<<"  "<<dii(1,0)<<"  "<<dii(2,0)<< endl;
		}
	}

	// 计算属性 rho_adv
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			const double &density = fluidModel.particleList[i]->density;
			double &densityAdv = fluidModel.particleList[i]->density_adv;
			densityAdv = density;
			Vector3f &xi = fluidModel.particleList[i]->position;
			Vector3f &vi = fluidModel.particleList[i]->velocity;
			double supportRad = fluidModel.particleList[i]->particleSupportRad;

			vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
			vector<StaticRigidParticle*> boundaryNeighbors = fluidModel.particleList[i]->boundaryNeighbors;
			//////////////////////////////////////////////////////////////////////////
			// 流体邻居部分计算
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
			{
				if (fluidNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = fluidNeighbors[j]->position;
				Vector3f &vj = fluidNeighbors[j]->velocity;
				densityAdv += h* fluidNeighbors[j]->particleMass * (vi - vj).Dot(CubicKernelOne::gradW(supportRad,xi - xj));
			}

			//////////////////////////////////////////////////////////////////////////
			// 边界邻居部分计算
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < boundaryNeighbors.size(); j++)
			{
				if (boundaryNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = boundaryNeighbors[j]->position;
				Vector3f &vj = boundaryNeighbors[j]->velocity;
				densityAdv += h* boundaryNeighbors[j]->boundaryPsi * (vi - vj).Dot(CubicKernelOne::gradW(supportRad,xi - xj));
			}
			//cout << densityAdv << endl;

			const double &pressure = fluidModel.particleList[i]->pressure;
			double &lastPressure = fluidModel.particleList[i]->lastPressure;
			lastPressure = 0.5*pressure;

			// 计算粒子 a_ii
			double &aii = fluidModel.particleList[i]->aii;
			aii = 0.0;
			Vector3f &dii = fluidModel.particleList[i]->dii;

			const double dpi = fluidModel.particleList[i]->particleMass / (density*density);

			//////////////////////////////////////////////////////////////////////////
			// 流体部分计算
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
			{
				if (fluidNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = fluidNeighbors[j]->position;
				// Compute d_ji
				const Vector3f kernel = CubicKernelOne::gradW(supportRad, xi - xj);
				Vector3f dji = dpi * kernel;
				aii += fluidNeighbors[j]->particleMass * (dii - dji).Dot(kernel);
			}

			//////////////////////////////////////////////////////////////////////////
			// 边界粒子计算
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < boundaryNeighbors.size(); j++)
			{
				if (boundaryNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = boundaryNeighbors[j]->position;
				const Vector3f kernel = CubicKernelOne::gradW(supportRad, xi - xj);
				Vector3f dji = dpi * kernel;
				aii += boundaryNeighbors[j]->boundaryPsi * (dii - dji).Dot(kernel);
			}
			//cout << aii << endl;
		}
	}
}

//处理压力解算器
void IISPHComputer::ProcessPressureSolve()
{
	cout << "ProcessPressureSolve()" << endl;
	const unsigned int numParticles = fluidModel.particleList.size();

	const double density0 = fluidModel.IISPH_RestDensity;
	const double h = fluidModel.IISPH_TimeStep;
	const double h2 = h*h;

	m_iterations = 0;
	const double omega = 0.5;

	const double eta = m_maxError * 0.01 * density0;  

	double avg_density = 0.0;
	//((avg_density - density0) > eta) ||
	while (( (m_iterations < 2)) && (m_iterations < m_maxIterations))
	{
		avg_density = 0.0;

		// 计算 dij_pj
       #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				if (fluidModel.particleList[i] == NULL)
				{
					continue;
				}
				Vector3f &dij_pj = fluidModel.particleList[i]->dij_pj;
				dij_pj.Zero();
				Vector3f &xi = fluidModel.particleList[i]->position;
				double supportRad = fluidModel.particleList[i]->particleSupportRad;

				vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
				vector<StaticRigidParticle*> boundaryNeighbors = fluidModel.particleList[i]->boundaryNeighbors;

				//////////////////////////////////////////////////////////////////////////
				// 流体部分计算
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
				{
					if (fluidNeighbors[j] == NULL)
					{
						continue;
					}
					Vector3f &xj = fluidNeighbors[j]->position;
					const double &densityj = fluidNeighbors[j]->density;
					dij_pj -= fluidNeighbors[j]->particleMass / (densityj*densityj) * fluidNeighbors[j]->lastPressure * CubicKernelOne::gradW(supportRad, xi - xj);
				}
			}
		}

		// 计算新的压力
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				if (fluidModel.particleList[i] == NULL)
				{
					continue;
				}
				const double &aii = fluidModel.particleList[i]->aii;
				const double density = fluidModel.particleList[i]->density;
				Vector3f &xi = fluidModel.particleList[i]->position;
				double dpi =  fluidModel.particleList[i]->particleMass / (density*density);
				double sum = 0.0;
				double supportRad = fluidModel.particleList[i]->particleSupportRad;

				vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
				vector<StaticRigidParticle*> boundaryNeighbors = fluidModel.particleList[i]->boundaryNeighbors;
				//////////////////////////////////////////////////////////////////////////
				// 流体粒子模拟
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
				{
					if (fluidNeighbors[j] == NULL)
					{
						continue;
					}
					Vector3f &xj = fluidNeighbors[j]->position;
					Vector3f &d_jk_pk = fluidNeighbors[j]->dij_pj;
					// Compute \sum_{k \neq i} djk*pk
					// Compute d_ji
					const Vector3f kernel = CubicKernelOne::gradW(supportRad, xi - xj);
					const Vector3f dji = dpi * kernel;
					Vector3f d_ji_pi = dji * fluidNeighbors[j]->lastPressure;

					// \sum ( mj * (\sum dij*pj - djj*pj - \sum_{k \neq i} djk*pk) * model->gradW)
					sum += fluidNeighbors[j]->particleMass * (fluidNeighbors[j]->dij_pj - fluidNeighbors[j]->dii * fluidNeighbors[j]->lastPressure - (d_jk_pk - d_ji_pi)).Dot(kernel);
				}

				//////////////////////////////////////////////////////////////////////////
				// 边界粒子
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < boundaryNeighbors.size(); j++)
				{
					if (boundaryNeighbors[j] == NULL)
					{
						continue;
					}
					Vector3f &xj = boundaryNeighbors[j]->position;
					sum += boundaryNeighbors[j]->boundaryPsi * fluidModel.particleList[i]->dij_pj.Dot(CubicKernelOne::gradW(supportRad, xi - xj));
				}

				const double b = density0 - fluidModel.particleList[i]->density_adv;
				double &pi = fluidModel.particleList[i]->pressure;
				const double &lastPi = fluidModel.particleList[i]->lastPressure;
				const double denom = aii*h2;
				if (fabs(denom) > 1.0e-9)
					pi = max((1.0 - omega)*lastPi + omega / denom * (b - h2*sum), 0.0);
				else
					pi = 0.0;

				if (pi != 0.0)
				{
					const double newDensity = (aii*pi + sum)*h2 - b + density0;
                    #pragma omp atomic
					avg_density += newDensity;
				}
				else
				{
                    #pragma omp atomic
					avg_density += density0;
				}
			}
		}

		for (int i = 0; i < (int)numParticles; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			const double &pi = fluidModel.particleList[i]->pressure;
			double &lastPi = fluidModel.particleList[i]->lastPressure;
			lastPi = pi;
		}

		avg_density /= numParticles;

		m_iterations++;
	}
	/*for (int i = 0; i < numParticles; i++)
	{
		cout << m_simulationData.getDij_pj(i)(0,0)<<" "<< m_simulationData.getDij_pj(i)(1, 0)<<" "<< m_simulationData.getDij_pj(i)(2, 0) << endl;
	}*/
}

void IISPHComputer::ComputePressureAccels()
{
	const unsigned int numParticles = fluidModel.particleList.size();
	//计算压力
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			Vector3f &xi = fluidModel.particleList[i]->position;
			const double &density_i = fluidModel.particleList[i]->density;

			Vector3f &ai = fluidModel.particleList[i]->pressureAcceleration;
			ai.Zero();
			double dpi = fluidModel.particleList[i]->pressure / (density_i*density_i);
			double supportRad = fluidModel.particleList[i]->particleSupportRad;

			vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
			vector<StaticRigidParticle*> boundaryNeighbors = fluidModel.particleList[i]->boundaryNeighbors;

			//////////////////////////////////////////////////////////////////////////
			// 计算边界粒子
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
			{
				if (fluidNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = fluidNeighbors[j]->position;
				// Pressure 
				const double &density_j = fluidNeighbors[j]->density;
				double dpj = fluidNeighbors[j]->pressure / (density_j*density_j);
				ai -= fluidNeighbors[j]->particleMass * (dpi + dpj) * CubicKernelOne::gradW(supportRad,xi - xj);
			}

			//////////////////////////////////////////////////////////////////////////
			// 计算边界粒子
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < boundaryNeighbors.size(); j++)
			{
				if (boundaryNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = boundaryNeighbors[j]->position;
				const Vector3f a = boundaryNeighbors[j]->boundaryPsi * (dpi)* CubicKernelOne::gradW(supportRad,xi - xj);
				ai -= a;
				//model->getForce(pid, neighborIndex) += model->getMass(i) * a;
			}
		}
	}
}

void IISPHComputer::ComputePressureAccels2()
{
	const unsigned int numParticles = fluidModel.particleList.size();
	//计算压力
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			Vector3f &xi = fluidModel.particleList[i]->position;
			const double &density_i = fluidModel.particleList[i]->density;

			Vector3f &ai = fluidModel.particleList[i]->pressureAcceleration;
			ai.Zero();
			double Ωi = fluidModel.particleList[i]->O;
			//cout << Ωi << endl;
			double dpi = fluidModel.particleList[i]->pressure / (density_i*density_i*Ωi);
			double supportRad = fluidModel.particleList[i]->particleSupportRad;

			vector<IISPHParticle*> fluidNeighbors = fluidModel.particleList[i]->fluidNeighbors;
			vector<StaticRigidParticle*> boundaryNeighbors = fluidModel.particleList[i]->boundaryNeighbors;

			//////////////////////////////////////////////////////////////////////////
			// 计算流体粒子
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < fluidNeighbors.size(); j++)
			{
				if (fluidNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = fluidNeighbors[j]->position;
				// Pressure 
				const double &density_j = fluidNeighbors[j]->density;
				//!!!!!!!!!!!!!!!!!!!!!!!!!!
				//double Ωj = fluidModel.particleList[j]->O;
				double Ωj = fluidNeighbors[j]->O;
				double dpj = fluidNeighbors[j]->pressure / (density_j*density_j*Ωj);
				ai -= fluidNeighbors[j]->particleMass * (dpi + dpj) * CubicKernelOne::gradW(supportRad, xi - xj);
			}

			//////////////////////////////////////////////////////////////////////////
			// 计算边界粒子
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < boundaryNeighbors.size(); j++)
			{
				if (boundaryNeighbors[j] == NULL)
				{
					continue;
				}
				Vector3f &xj = boundaryNeighbors[j]->position;
				const Vector3f a = boundaryNeighbors[j]->boundaryPsi * (dpi)* CubicKernelOne::gradW(supportRad, xi - xj);
				ai -= a;
				//model->getForce(pid, neighborIndex) += model->getMass(i) * a;
			}
		}
	}
}

//计算自适应粒子
void IISPHComputer::ComputeAdaptivityFactor()
{
	//计算得到一个自适应因子
	adaptivityFactor = fineMass / baseMass;
}

void IISPHComputer::ComputeDistance()
{
}

//为例子进行分类
void IISPHComputer::ClassifyParticle()
{
	//依据当前的质量、质量比值，进而完成对应的标记
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		for (int i=0; i< fluidModel.particleList.size(); i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			//得到当前粒子的质量
			double massi = fluidModel.particleList[i]->particleMass;
			//得到当前粒子的最佳质量
			double mopt = fluidModel.particleList[i]->mopt;
			//然后得到当前的粒子质量的比值
			double mrel = massi / mopt;
			//得到当前粒子的标记
			PARTICLE_TYPE& particleT = fluidModel.particleList[i]->pt;
			//开始进行标签化分类
			if (mrel < 0.5f)
			{
				particleT = PARTICLE_TYPE::S;
			}
			else if (mrel >= 0.5 && mrel <= 0.9)
			{
				particleT = PARTICLE_TYPE::s;
			}
			else if (mrel > 0.9 && mrel < 1.1)
			{
				particleT = PARTICLE_TYPE::o;
			}
			else if (mrel >= 1.1&&mrel <= 2.0)
			{
				particleT = PARTICLE_TYPE::l;
			}
			else if(mrel>2.0)
			{
				particleT = PARTICLE_TYPE::L;
			}
		}
	}
}

//计算粒子最优质量
void IISPHComputer::ComputeOptMass()
{
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)  
		//遍历当前所有粒子的列表，计算得到每一个粒子的最佳质量
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			double& mopt = fluidModel.particleList[i]->mopt;
			double distance = fluidModel.particleList[i]->distance;
			//计算得到当前粒子的最优质量
			mopt = baseMass*((min(distance, baseDistance) / baseDistance)*(1 - adaptivityFactor) + adaptivityFactor);
		}
	}
}

//计算欧米伽
void IISPHComputer::ComputeO()
{
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		for (int i = 0; i < fluidModel.particleList.size(); i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			//得到Ω的引用
			double& O = fluidModel.particleList[i]->O;
			//得到当前粒子的位置
			Vector3f xi = fluidModel.particleList[i]->position;
			//得到当前粒子的密度
			const double di = fluidModel.particleList[i]->density;
			//得到当前粒子的粒子支撑半径
			const double hi = fluidModel.particleList[i]->particleSupportRad;
			//得到当前粒子的粒子质量
			const double massi = fluidModel.particleList[i]->particleMass;
			//得到当前粒子的类型
			PARTICLE_TYPE particleT = fluidModel.particleList[i]->pt;
			if (particleT == PARTICLE_TYPE::l)
			{
				//这里求得是WZero的偏导数
				O = 1 + hi / (3 * di)*massi*Poly6KernelOne::WZeroPartialDerivativeH(hi);
			}
			else
			{
				//得到当前粒子的流体邻居
				vector<IISPHParticle*> fluidNeibors = fluidModel.particleList[i]->fluidNeighbors;
		        //在此基础上，运算W针对邻居的偏导数
				for (int j = 0; j < fluidNeibors.size(); j++)
				{
					if (fluidNeibors[j] == NULL)
					{
						continue;
					}
					//得到邻居的粒子质量
					Vector3f xj = fluidNeibors[j]->position;
					const double mj = fluidNeibors[j]->particleMass;
					//W函数对h的偏导数
					O += mj * Poly6KernelOne::WPartialDerivativeH(hi, xi - xj);
				}
				O = O*(hi / (3 * di)) + 1;
			}
		}
	}
}

//混合当前密度
void IISPHComputer::BlendDensity()
{
	//注：其中贝塔为0.5 则为第一次分裂，贝塔初始化为0.2，则为第一次合并
	int particleNum = fluidModel.particleList.size();
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (int i = 0; i < particleNum; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			//得到当前粒子是否是更新后的第一帧
			const bool newUpdate = fluidModel.particleList[i]->NEW_UPDATE;
			if (newUpdate)
			{
				//得到当前粒子速度
				double& density = fluidModel.particleList[i]->density;
				//得到当前粒子的混合因子
				const double beita = fluidModel.particleList[i]->β;
				//得到当前粒子的lastDensity
				const double density0 = fluidModel.particleList[i]->lastDensity;
				//ρblended i = (1 − βi)*ρi + βi* ρO 刷新当前的密度
				density = (1 - beita)*density + beita*density0;
			}
		}
	}
}

//混合当前速度
void IISPHComputer::BlendVelocity()
{
	int particleNum = fluidModel.particleList.size();
     #pragma omp parallel default(shared)
	{
		for (int i = 0; i < particleNum; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			//得到当前是否是第一次搞的
			bool& newUpdate = fluidModel.particleList[i]->NEW_UPDATE;
			if (newUpdate)
			{
				//在第一次更新之后，将第一个标志位设置为false
				newUpdate = false;
				//得到当前粒子速度
				Vector3f& velocity = fluidModel.particleList[i]->velocity;
				//得到当前粒子的混合因子
				const double beita = fluidModel.particleList[i]->β;
				//得到当前粒子的lastDensity
				const Vector3f velocity0 = fluidModel.particleList[i]->lastVelocity;
				//ρblended i = (1 − βi)*ρi + βi* ρO 刷新当前的密度
				velocity = (1 - beita)*velocity + beita*velocity0;
			}
		}
	}
}

//每一个粒子都在遍历当中刷新当前的混合因子
void IISPHComputer::UpdateBlendFactor()
{
	int particleNum = fluidModel.particleList.size();
	double d = 0.0001f;
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static)
		for (int i = 0; i < particleNum; i++)
		{
			if (fluidModel.particleList[i] == NULL)
			{
				continue;
			}
			//如果当前粒子>0 那么就减0.1
			if (fluidModel.particleList[i]->β > 0 + d)
			{
				fluidModel.particleList[i]->β -= 0.1;
			}
			//否则不对当前粒子做操作
		}
	}
}


//add


void IISPHComputer::ComputePtoSandMopt()
{
	cout << "ComputePtoSandMopt()" << endl;
	baseDistance = -10.0;
	fineDistance = 10.0;
	//得到当前粒子个数
	const int numParticles = fluidModel.particleList.size();
	//遍历所有粒子
	for (int i = 0; i < (int)numParticles; i++)
	{
		if (fluidModel.particleList[i] == NULL)
		{
			continue;
		}
		double _x = 0.0;
		double _y = 0.0;
		double _z = 0.0;
		double _r = 0.0;
		double wj = 0.0;//权值
						//double distValue = 0.0;
						//double disToSurface = 0.0;

		IISPHParticle *pi;
		IISPHParticle *pj;
		IISPHParticle *pk;

		pi = fluidModel.particleList[i];
		//计算distance field
		for (int j = 0; j < fluidModel.particleList[i]->fluidNeighbors.size(); j++)
		{
			if (fluidModel.particleList[i]->fluidNeighbors[j] == NULL)
			{
				continue;
			}
			pj = fluidModel.particleList[i]->fluidNeighbors[j];
			double rj = fluidModel.particleList[i]->fluidNeighbors[j]->particleRad;
			double Rj = fluidModel.particleList[i]->fluidNeighbors[j]->particleSupportRad;

			double uws = abs(sqrt(pow(pj->position.x - pi->position.x, 2) + pow(pj->position.y - pi->position.y, 2) + pow(pj->position.z - pi->position.z, 2))) / Rj;
			double uw = max(0.0, pow((1 - uws*uws), 3));

			double dw = 0.0;
			for (int k = 0; k < fluidModel.particleList[i]->fluidNeighbors.size(); k++)
			{
				if (fluidModel.particleList[i]->fluidNeighbors[k] == NULL)
				{
					continue;
				}
				pk = fluidModel.particleList[i]->fluidNeighbors[k];

				double dws = abs(sqrt(pow(pk->position.x - pi->position.x, 2) + pow(pk->position.y - pi->position.y, 2) + pow(pk->position.z - pi->position.z, 2))) / Rj;
				dw += max(0.0, pow((1 - dws*dws), 3));

			}

			wj = ((uw == 0 && dw == 0) ? 1 : (uw / dw));

			_x += wj*pj->position.x;
			_y += wj*pj->position.y;
			_z += wj*pj->position.z;

			_r += wj * pj->particleRad;


		}
		//粒子i的距离域值
		pi->distValue = abs(sqrt(pow((pi->position.x - _x), 2) + pow((pi->position.y - _y), 2) + pow((pi->position.z - _z), 2))) - _r;
		//cout << "pi->distValue: " << pi->distValue << endl;
		//求粒子i的到表面距离
		double minDIST = -1.5*pi->particleRad;
		double maxDIST = +1.5*pi->particleRad;
		if (pi->fluidNeighbors.size() < 5)
		{
			pi->disToSurface = -0.85*pi->particleRad;
		}
		else if (pi->distValue < minDIST)
		{
			pi->disToSurface = minDIST;
		}
		else if (pi->distValue > maxDIST)
		{
			pi->disToSurface = maxDIST;
		}
		else
		{
			pi->disToSurface = pi->distValue;
		}
		//cout << pi->disToSurface << endl;

		//当距离域值大于-0.85*pi->particleRad时，将粒子标记为表面粒子
		if (pi->disToSurface >= -0.85*pi->particleRad)
		{
			pi->mark = 0;
		}
		//将其他粒子标记为内部粒子
		else
		{
			pi->mark = 1;
		}

		if (pi->disToSurface < -0.85*pi->particleRad)
		{
			double tmpDIST = -18.0*pi->particleRad;
			for (int m = 0; m < pi->fluidNeighbors.size(); m++)
			{
				if (pi->fluidNeighbors[m] == NULL)
				{
					continue;
				}
				IISPHParticle *pm = pi->fluidNeighbors[m];

				if (0 == pm->mark || 10 == pm->mark)
				{
					double t = pi->disToSurface - sqrt(pow(pi->position.x - pm->position.x, 2) + pow(pi->position.y - pm->position.y, 2) + pow(pi->position.z - pm->position.z, 2));
					if (tmpDIST < t && t < -0.85*pi->particleRad)
					{
						tmpDIST = t;
					}
					//cout << "aaaa: " << t << endl;
				}


			}
			//cout << "tmpDIST: " << tmpDIST << endl;
			//求出来的disToSurface值为负，表面粒子在流体内部
			pi->disToSurface = tmpDIST;
			//将粒子标记为接近表面
			pi->mark = 10;
		}
		pi->disToSurface = abs(pi->disToSurface);
		//cout << "pi->disToSurface: " << pi->disToSurface << endl;

		if (fineDistance > pi->disToSurface)
		{
			fineDistance = pi->disToSurface;
		}

		if (baseDistance < pi->disToSurface)
		{
			baseDistance = pi->disToSurface;
		}

	}
	//cout << "maxxxxx: " << baseDistance << endl;
	//cout << "minnnnn: " << fineDistance << endl;
	for (int i = 0; i < fluidModel.particleList.size(); i++)
	{
		if (fluidModel.particleList[i] == NULL)
		{
			continue;
		}
		IISPHParticle *pi = fluidModel.particleList[i];
		if (pi->disToSurface <= 0.85*pi->particleRad)
		{
			pi->mopt = fineMass;
			//cout << particles[i].mopt << endl;
		}
		else
		{

			pi->mopt = baseMass*((pi->disToSurface / baseDistance)*(1.0 - fineMass / baseMass) + adaptivityFactor);

		}
		//cout <<i <<" "<< fluidModel.particleList[i]->mopt << endl;

	}
	//cout << fluidModel.particleList[512]->mopt << endl;
}
void IISPHComputer::ProcessMergeandSplit()
{
	cout << "ProcessMergeandSplit()" << endl;
	for (int i = 0; i < fluidModel.particleList.size(); i++)
	{
		if (fluidModel.particleList[i] == NULL)
		{
			continue;
		}

		//ProcessNeighborList();//搜索邻居方法需要改，加判断，只有指针指向不为空的 才添加到粒子i的邻居内部；或者重写一个方法ProcessNeighborList(int i),
							  //在每次遍历粒子列表的时候都要添加指针不为空的判断；
		IISPHParticle *pi = fluidModel.particleList[i];

		double m_rel = pi->particleMass / pi->mopt;
		//cout << "mrel: " << m_rel <<"--"<< (m_rel == 0.5) << endl;
		//对S粒子进行处理
		if (m_rel <= 0.5 && pi->β<=0)
		{
			cout << "S粒子sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss" << endl;
			//S粒子
			int num_S_neighbor = 0;
			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
			{
				if (pi->fluidNeighbors[j] == NULL)
				{
					continue;
				}

				IISPHParticle *pj = pi->fluidNeighbors[j];
				if ((pj->particleMass / pj->mopt >= 0.5&&pj->particleMass / pj->mopt <= 0.9) || (pj->particleMass / pj->mopt < 0.5))
				{
					num_S_neighbor = num_S_neighbor + 1;
				}

			}
			double mn = pi->particleMass / num_S_neighbor;

			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
			{
				if (pi->fluidNeighbors[j] == NULL)
				{
					continue;
				}

				IISPHParticle *pj = pi->fluidNeighbors[j];
				if ((pj->particleMass / pj->mopt >= 0.5&&pj->particleMass / pj->mopt <= 0.9) || (pj->particleMass / pj->mopt < 0.5))
				{
					pj->position = (pi->position*mn + pj->position*pj->particleMass) / (mn + pj->particleMass);
					pj->velocity = (pi->velocity*mn + pj->velocity*pj->particleMass) / (mn + pj->particleMass);
					pj->particleMass = pj->particleMass + mn;

					pj->β = 0.2;
				}
				int qq = 11;

			}

			//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
			{
				//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，pk指向设为空
				for (int j = 0; j < pi->fluidNeighbors.size(); j++)
				{
					if (pi->fluidNeighbors[j] == NULL)
					{
						continue;
					}

					IISPHParticle *pj = pi->fluidNeighbors[j];
					for (int k = 0; k < pj->fluidNeighbors.size(); k++)
					{
						if (pj->fluidNeighbors[k] == NULL)
						{
							continue;
						}

						IISPHParticle *pk = pj->fluidNeighbors[k];
						if (pk == pi)
						{
							pj->fluidNeighbors[k] = NULL;
						}
					}

				}
				delete fluidModel.particleList[i];
				fluidModel.particleList[i] = NULL;
				//int qqq = 43;

			}
			//pi->β = 0.2;
			//cout << "pi->β ："<< i <<" " << pi->β << endl;
		}

		//-------------------------------------------------------------------------------
		//l 粒子
		if (m_rel >= 1.1 && m_rel < 2 && pi->β <= 0)
		{
			cout << "l粒子lllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll" << endl;
			//S粒子
			double mex = pi->particleMass - pi->mopt;

			int num_S_neighbor = 0;
			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
			{
				if (pi->fluidNeighbors[j] == NULL)
				{
					continue;
				}

				IISPHParticle *pj = pi->fluidNeighbors[j];
				if ((pj->particleMass / pj->mopt >= 0.5&&pj->particleMass / pj->mopt <= 0.9) || (pj->particleMass / pj->mopt < 0.5))
				{
					num_S_neighbor = num_S_neighbor + 1;
				}

			}
			double mn = mex / num_S_neighbor;

			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
			{
				if (pi->fluidNeighbors[j] == NULL)
				{
					continue;
				}

				IISPHParticle *pj = pi->fluidNeighbors[j];
				if ((pj->particleMass / pj->mopt >= 0.5&&pj->particleMass / pj->mopt <= 0.9) || (pj->particleMass / pj->mopt < 0.5))
				{
					pj->position = (pi->position*mn + pj->position*pj->particleMass) / (mn + pj->particleMass);
					pj->velocity = (pi->velocity*mn + pj->velocity*pj->particleMass) / (mn + pj->particleMass);
					pj->particleMass = pj->particleMass + mn;

					pj->β = 0.2;
				}
				//int qq = 11;

			}

			//pi->β = 0.2;
			//cout << "pi->β ："<< i <<" " << pi->β << endl;
		}

		//-------------------------------------------------------------------------------
		//对L粒子
		//此方法中，生成的新粒子的position有待优化……

		//if (isnormal(m_rel) && m_rel > 2 && pi->β <= 0)
		//{
		//	if (pi->density == 0)
		//	{
		//		cout << "pi->density == 0" << endl;
		//	}

		//	cout << "pi->denstiy" << pi->density << endl;
		//	cout << "L粒子LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL" << endl;
		//	int num_np = ceil(pi->particleMass/pi->mopt);//分裂成的新粒子的个数。结果向上取整

		//	if (num_np == 2)
		//	{
		//		cout << "2222222222222222222222222222222222222222222222222" << endl;

		//		double r_Spliting = 0.3*pi->particleSupportRad;//新生成的粒子所在球面的半径，0.3为取得自定义的值，此处用来测试
		//		double mass_np = pi->particleMass / num_np;
		//		

		//		IISPHParticle *np1 = new IISPHParticle(*pi);
		//		IISPHParticle *np2 = new IISPHParticle(*pi);

		//		np1->lastPressure = pi->pressure;
		//		np2->lastPressure = pi->pressure;

		//		Vector3f p1, p2;
		//		p1.SetValue(pi->position.x, pi->position.y, pi->position.z+ r_Spliting);

		//		int hashValueP1 = FunctionKit::PositionMapHash(p1,
		//			SharedData::GetCellLength(),
		//			SharedData::GetHashCellNum(),
		//			boundary.offset); //得到当前position所对应的哈希值
		//		cout << "pi->denstiy" << pi->density << endl;
		//		cout << "massssssssssssssssssssssssss_npppppppppppppppppppppppp" << mass_np << endl;
		//		cout << "massssssssssssssssssssssssss_npppppppppppppppppppppppphashValueP1" << hashValueP1 << endl;
		//		cout << "massssssssssssssssssssssssss_nppppppppppppppppppppppppdensity" << pi->density << endl;
		//		cout << "massssssssssssssssssssssssss_nppppppppppppppppppppppppp1" << p1 << endl;

		//		np1->InitializationPVS(p1, hashValueP1, mass_np, pi->density);//需要添加构造函数（至少）：mass，position，velocity, density
		//		
		//		//vector<IISPHParticle*> tempnp1nei = ComputeNeighborParticle(np1);
		//		////fluidModel.particleList[i]->fluidNeighbors =temp;
		//		//for (int j = 0; j < tempnp1nei.size(); j++)
		//		//{
		//		//	if (tempnp1nei[j] == NULL)
		//		//	{
		//		//		continue;
		//		//	}
		//		//	//将临近的所有粒子装入对应的粒子邻居表当中
		//		//	np1->fluidNeighbors.push_back(tempnp1nei[j]);
		//		//}


		//		int hashValueP2 = FunctionKit::PositionMapHash(p2,
		//			SharedData::GetCellLength(),
		//			SharedData::GetHashCellNum(),
		//			boundary.offset); //得到当前position所对应的哈希值
		//		p2.SetValue(pi->position.x, pi->position.y, pi->position.z - r_Spliting);
		//		np2->InitializationPVS(p2, hashValueP2, mass_np, pi->density);

		//		np1->β = 0.5;
		//		np2->β = 0.5;
		//		fluidModel.particleList.push_back(np1);
		//		fluidModel.particleList.push_back(np2);

		//		//删除粒子部分代码开始标志--------------------------------------------------------------------------------------------
		//		//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
		//		{
		//			//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，新增其指向生成粒子的指针
		//			//aa
		//			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
		//			{
		//				if (pi->fluidNeighbors[j] == NULL)
		//				{
		//					continue;
		//				}

		//				IISPHParticle *pj = pi->fluidNeighbors[j];
		//				for (int k = 0; k < pj->fluidNeighbors.size(); k++)
		//				{
		//					if (pj->fluidNeighbors[k] == NULL)
		//					{
		//						continue;
		//					}

		//					IISPHParticle *pk = pj->fluidNeighbors[k];
		//					if (pk == pi)
		//					{
		//						pj->fluidNeighbors[k] = NULL;
		//						pj->fluidNeighbors.push_back(np1);
		//						pj->fluidNeighbors.push_back(np2);
		//						
		//					}
		//				}

		//			}
		//			delete fluidModel.particleList[i];
		//			fluidModel.particleList[i] = NULL;
		//			//int qqq = 43;

		//		}

		//		//删除粒子部分代码结束标志----------------------------------------------------------------------------------------------------
		//	}

		//	if (num_np == 3)
		//	{
		//		cout << "333333333333333333333333333333333333333333333333333" << endl;

		//		double r_Spliting = 0.3*pi->particleSupportRad;//新生成的粒子所在球面的半径，0.3为取得自定义的值，此处用来测试
		//		double mass_np = pi->particleMass / num_np;
		//		double c = sqrt(3.0)*r_Spliting*0.5;

		//		IISPHParticle *np1 = new IISPHParticle(*pi);
		//		IISPHParticle *np2 = new IISPHParticle(*pi);
		//		IISPHParticle *np3 = new IISPHParticle(*pi);

		//		np1->lastPressure = pi->pressure;
		//		np2->lastPressure = pi->pressure;
		//		np3->lastPressure = pi->pressure;

		//		Vector3f p1, p2, p3;
		//		p1.SetValue(pi->position.x, r_Spliting+pi->position.y, pi->position.z);
		//		p2.SetValue(-c+pi->position.x, -0.5*r_Spliting+pi->position.y, pi->position.z);
		//		p3.SetValue(c+pi->position.x, -0.5*r_Spliting + pi->position.y, pi->position.z);

		//		int hashValueP1 = FunctionKit::PositionMapHash(p1,SharedData::GetCellLength(),SharedData::GetHashCellNum(),boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP2 = FunctionKit::PositionMapHash(p2,SharedData::GetCellLength(),SharedData::GetHashCellNum(),boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP3 = FunctionKit::PositionMapHash(p3, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		
		//		np1->InitializationPVS(p1, hashValueP1, mass_np, pi->density);//需要添加构造函数（至少）：mass，position，velocity, density
		//		np2->InitializationPVS(p2, hashValueP2, mass_np, pi->density);
		//		np3->InitializationPVS(p3, hashValueP3, mass_np, pi->density);
		//		
		//		np1->β = 0.5;
		//		np2->β = 0.5;
		//		np3->β = 0.5;

		//		fluidModel.particleList.push_back(np1);
		//		fluidModel.particleList.push_back(np2);
		//		fluidModel.particleList.push_back(np3);

		//		//删除粒子部分代码开始标志--------------------------------------------------------------------------------------------
		//		//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
		//		{
		//			//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，新增其指向生成粒子的指针
		//			//aa
		//			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
		//			{
		//				if (pi->fluidNeighbors[j] == NULL)
		//				{
		//					continue;
		//				}

		//				IISPHParticle *pj = pi->fluidNeighbors[j];
		//				for (int k = 0; k < pj->fluidNeighbors.size(); k++)
		//				{
		//					if (pj->fluidNeighbors[k] == NULL)
		//					{
		//						continue;
		//					}

		//					IISPHParticle *pk = pj->fluidNeighbors[k];
		//					if (pk == pi)
		//					{
		//						pj->fluidNeighbors[k] = NULL;
		//						pj->fluidNeighbors.push_back(np1);
		//						pj->fluidNeighbors.push_back(np2);
		//						pj->fluidNeighbors.push_back(np3);


		//					}
		//				}

		//			}
		//			delete fluidModel.particleList[i];
		//			fluidModel.particleList[i] = NULL;
		//			//int qqq = 43;

		//		}

		//		//删除粒子部分代码结束标志----------------------------------------------------------------------------------------------------
		//	}

		//	if (num_np == 4)
		//	{
		//		cout << "4444444444444444444444444444444444444444444444444" << endl;

		//		double r_Spliting = 0.3*pi->particleSupportRad;//新生成的粒子所在球面的半径，0.3为取得自定义的值，此处用来测试
		//		double mass_np = pi->particleMass / num_np;
		//		double c = 1.0 / sqrt(3.0)*r_Spliting;

		//		IISPHParticle *np1 = new IISPHParticle(*pi);
		//		IISPHParticle *np2 = new IISPHParticle(*pi);
		//		IISPHParticle *np3 = new IISPHParticle(*pi);
		//		IISPHParticle *np4 = new IISPHParticle(*pi);

		//		np1->lastPressure = pi->pressure;
		//		np2->lastPressure = pi->pressure;
		//		np3->lastPressure = pi->pressure;
		//		np4->lastPressure = pi->pressure;

		//		Vector3f p1, p2, p3, p4;
		//		p1.SetValue(c + pi->position.x, c + pi->position.y, c + pi->position.z);
		//		p2.SetValue(-c + pi->position.x, c + pi->position.y, -c + pi->position.z);
		//		p3.SetValue(-c + pi->position.x, -c + pi->position.y, c + pi->position.z);
		//		p4.SetValue(c + pi->position.x, -c + pi->position.y, -c + pi->position.z);

		//		int hashValueP1 = FunctionKit::PositionMapHash(p1, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP2 = FunctionKit::PositionMapHash(p2, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP3 = FunctionKit::PositionMapHash(p3, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP4 = FunctionKit::PositionMapHash(p4, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		
		//		np1->InitializationPVS(p1, hashValueP1, mass_np, pi->density);//需要添加构造函数（至少）：mass，position，velocity, density
		//		np2->InitializationPVS(p2, hashValueP2, mass_np, pi->density);
		//		np3->InitializationPVS(p3, hashValueP3, mass_np, pi->density);
		//		np4->InitializationPVS(p4, hashValueP4, mass_np, pi->density);

		//		np1->β = 0.5;
		//		np2->β = 0.5;
		//		np3->β = 0.5;
		//		np4->β = 0.5;

		//		fluidModel.particleList.push_back(np1);
		//		fluidModel.particleList.push_back(np2);
		//		fluidModel.particleList.push_back(np3);
		//		fluidModel.particleList.push_back(np4);


		//		//删除粒子部分代码开始标志--------------------------------------------------------------------------------------------
		//		//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
		//		{
		//			//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，新增其指向生成粒子的指针
		//			//aa
		//			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
		//			{
		//				if (pi->fluidNeighbors[j] == NULL)
		//				{
		//					continue;
		//				}

		//				IISPHParticle *pj = pi->fluidNeighbors[j];
		//				for (int k = 0; k < pj->fluidNeighbors.size(); k++)
		//				{
		//					if (pj->fluidNeighbors[k] == NULL)
		//					{
		//						continue;
		//					}

		//					IISPHParticle *pk = pj->fluidNeighbors[k];
		//					if (pk == pi)
		//					{
		//						pj->fluidNeighbors[k] = NULL;
		//						pj->fluidNeighbors.push_back(np1);
		//						pj->fluidNeighbors.push_back(np2);
		//						pj->fluidNeighbors.push_back(np3);
		//						pj->fluidNeighbors.push_back(np4);

		//					}
		//				}

		//			}
		//			delete fluidModel.particleList[i];
		//			fluidModel.particleList[i] = NULL;
		//		}

		//		//删除粒子部分代码结束标志----------------------------------------------------------------------------------------------------
		//	}

		//	if (num_np == 5)
		//	{
		//		cout << "5555555555555555555555555555555555555555555555" << endl;

		//		double r_Spliting = 0.3*pi->particleSupportRad;//新生成的粒子所在球面的半径，0.3为取得自定义的值，此处用来测试
		//		double mass_np = pi->particleMass / num_np;
		//		double c = 1.0 / sqrt(3.0)*r_Spliting;

		//		IISPHParticle *np1 = new IISPHParticle(*pi);
		//		IISPHParticle *np2 = new IISPHParticle(*pi);
		//		IISPHParticle *np3 = new IISPHParticle(*pi);
		//		IISPHParticle *np4 = new IISPHParticle(*pi);
		//		IISPHParticle *np5 = new IISPHParticle(*pi);

		//		np1->lastPressure = pi->pressure;
		//		np2->lastPressure = pi->pressure;
		//		np3->lastPressure = pi->pressure;
		//		np4->lastPressure = pi->pressure;
		//		np5->lastPressure = pi->pressure;

		//		Vector3f p1, p2, p3, p4, p5;
		//		p1.SetValue(c + pi->position.x, c + pi->position.y, c + pi->position.z);
		//		p2.SetValue(-c + pi->position.x, c + pi->position.y, -c + pi->position.z);
		//		p3.SetValue(-c + pi->position.x, -c + pi->position.y, c + pi->position.z);
		//		p4.SetValue(c + pi->position.x, -c + pi->position.y, -c + pi->position.z);
		//		p5.SetValue(pi->position.x, pi->position.y, pi->position.z);

		//		int hashValueP1 = FunctionKit::PositionMapHash(p1, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP2 = FunctionKit::PositionMapHash(p2, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP3 = FunctionKit::PositionMapHash(p3, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP4 = FunctionKit::PositionMapHash(p4, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP5 = FunctionKit::PositionMapHash(p5, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值

		//		np1->InitializationPVS(p1, hashValueP1, mass_np, pi->density);//需要添加构造函数（至少）：mass，position，velocity, density
		//		np2->InitializationPVS(p2, hashValueP2, mass_np, pi->density);
		//		np3->InitializationPVS(p3, hashValueP3, mass_np, pi->density);
		//		np4->InitializationPVS(p4, hashValueP4, mass_np, pi->density);
		//		np5->InitializationPVS(p5, hashValueP5, mass_np, pi->density);

		//		np1->β = 0.5;
		//		np2->β = 0.5;
		//		np3->β = 0.5;
		//		np4->β = 0.5;
		//		np5->β = 0.5;

		//		fluidModel.particleList.push_back(np1);
		//		fluidModel.particleList.push_back(np2);
		//		fluidModel.particleList.push_back(np3);
		//		fluidModel.particleList.push_back(np4);
		//		fluidModel.particleList.push_back(np5);


		//		//删除粒子部分代码开始标志--------------------------------------------------------------------------------------------
		//		//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
		//		{
		//			//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，新增其指向生成粒子的指针
		//			//aa
		//			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
		//			{
		//				if (pi->fluidNeighbors[j] == NULL)
		//				{
		//					continue;
		//				}

		//				IISPHParticle *pj = pi->fluidNeighbors[j];
		//				for (int k = 0; k < pj->fluidNeighbors.size(); k++)
		//				{
		//					if (pj->fluidNeighbors[k] == NULL)
		//					{
		//						continue;
		//					}

		//					IISPHParticle *pk = pj->fluidNeighbors[k];
		//					if (pk == pi)
		//					{
		//						pj->fluidNeighbors[k] = NULL;
		//						pj->fluidNeighbors.push_back(np1);
		//						pj->fluidNeighbors.push_back(np2);
		//						pj->fluidNeighbors.push_back(np3);
		//						pj->fluidNeighbors.push_back(np4);
		//						pj->fluidNeighbors.push_back(np5);

		//					}
		//				}

		//			}
		//			delete fluidModel.particleList[i];
		//			fluidModel.particleList[i] = NULL;
		//		}

		//		//删除粒子部分代码结束标志----------------------------------------------------------------------------------------------------
		//	}

		//	if (num_np == 6)
		//	{
		//		cout << "66666666666666666666666666666666666666666666666" << endl;

		//		double r_Spliting = 0.3*pi->particleSupportRad;//新生成的粒子所在球面的半径，0.3为取得自定义的值，此处用来测试
		//		double mass_np = pi->particleMass / num_np;
		//		double c = r_Spliting;

		//		IISPHParticle *np1 = new IISPHParticle(*pi);
		//		IISPHParticle *np2 = new IISPHParticle(*pi);
		//		IISPHParticle *np3 = new IISPHParticle(*pi);
		//		IISPHParticle *np4 = new IISPHParticle(*pi);
		//		IISPHParticle *np5 = new IISPHParticle(*pi);
		//		IISPHParticle *np6 = new IISPHParticle(*pi);

		//		np1->lastPressure = pi->pressure;
		//		np2->lastPressure = pi->pressure;
		//		np3->lastPressure = pi->pressure;
		//		np4->lastPressure = pi->pressure;
		//		np5->lastPressure = pi->pressure;
		//		np6->lastPressure = pi->pressure;

		//		Vector3f p1, p2, p3, p4, p5, p6;
		//		p1.SetValue(0 + pi->position.x, 0 + pi->position.y, c + pi->position.z);
		//		p2.SetValue(0 + pi->position.x, 0 + pi->position.y, -c + pi->position.z);
		//		p3.SetValue(c + pi->position.x, 0 + pi->position.y, 0 + pi->position.z);
		//		p4.SetValue(-c + pi->position.x,0 + pi->position.y, 0 + pi->position.z);
		//		p5.SetValue(0 + pi->position.x, -c + pi->position.y, 0 + pi->position.z);
		//		p6.SetValue(0 + pi->position.x, c + pi->position.y, 0 + pi->position.z);

		//		int hashValueP1 = FunctionKit::PositionMapHash(p1, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP2 = FunctionKit::PositionMapHash(p2, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP3 = FunctionKit::PositionMapHash(p3, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP4 = FunctionKit::PositionMapHash(p4, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP5 = FunctionKit::PositionMapHash(p5, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP6 = FunctionKit::PositionMapHash(p6, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值

		//		np1->InitializationPVS(p1, hashValueP1, mass_np, pi->density);//需要添加构造函数（至少）：mass，position，velocity, density
		//		np2->InitializationPVS(p2, hashValueP2, mass_np, pi->density);
		//		np3->InitializationPVS(p3, hashValueP3, mass_np, pi->density);
		//		np4->InitializationPVS(p4, hashValueP4, mass_np, pi->density);
		//		np5->InitializationPVS(p5, hashValueP5, mass_np, pi->density);
		//		np6->InitializationPVS(p6, hashValueP6, mass_np, pi->density);

		//		np1->β = 0.5;
		//		np2->β = 0.5;
		//		np3->β = 0.5;
		//		np4->β = 0.5;
		//		np5->β = 0.5;
		//		np6->β = 0.5;

		//		fluidModel.particleList.push_back(np1);
		//		fluidModel.particleList.push_back(np2);
		//		fluidModel.particleList.push_back(np3);
		//		fluidModel.particleList.push_back(np4);
		//		fluidModel.particleList.push_back(np5);
		//		fluidModel.particleList.push_back(np6);


		//		//删除粒子部分代码开始标志--------------------------------------------------------------------------------------------
		//		//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
		//		{
		//			//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，新增其指向生成粒子的指针
		//			//aa
		//			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
		//			{
		//				if (pi->fluidNeighbors[j] == NULL)
		//				{
		//					continue;
		//				}

		//				IISPHParticle *pj = pi->fluidNeighbors[j];
		//				for (int k = 0; k < pj->fluidNeighbors.size(); k++)
		//				{
		//					if (pj->fluidNeighbors[k] == NULL)
		//					{
		//						continue;
		//					}

		//					IISPHParticle *pk = pj->fluidNeighbors[k];
		//					if (pk == pi)
		//					{
		//						pj->fluidNeighbors[k] = NULL;
		//						pj->fluidNeighbors.push_back(np1);
		//						pj->fluidNeighbors.push_back(np2);
		//						pj->fluidNeighbors.push_back(np3);
		//						pj->fluidNeighbors.push_back(np4);
		//						pj->fluidNeighbors.push_back(np5);
		//						pj->fluidNeighbors.push_back(np6);

		//					}
		//				}

		//			}
		//			delete fluidModel.particleList[i];
		//			fluidModel.particleList[i] = NULL;
		//		}

		//		//删除粒子部分代码结束标志----------------------------------------------------------------------------------------------------
		//	}

		//	if (num_np == 7)
		//	{
		//		cout << "77777777777777777777777777777777777777777777777" << endl;

		//		double r_Spliting = 0.3*pi->particleSupportRad;//新生成的粒子所在球面的半径，0.3为取得自定义的值，此处用来测试
		//		double mass_np = pi->particleMass / num_np;
		//		double c = r_Spliting;

		//		IISPHParticle *np1 = new IISPHParticle(*pi);
		//		IISPHParticle *np2 = new IISPHParticle(*pi);
		//		IISPHParticle *np3 = new IISPHParticle(*pi);
		//		IISPHParticle *np4 = new IISPHParticle(*pi);
		//		IISPHParticle *np5 = new IISPHParticle(*pi);
		//		IISPHParticle *np6 = new IISPHParticle(*pi);
		//		IISPHParticle *np7 = new IISPHParticle(*pi);

		//		np1->lastPressure = pi->pressure;
		//		np2->lastPressure = pi->pressure;
		//		np3->lastPressure = pi->pressure;
		//		np4->lastPressure = pi->pressure;
		//		np5->lastPressure = pi->pressure;
		//		np6->lastPressure = pi->pressure;
		//		np7->lastPressure = pi->pressure;

		//		Vector3f p1, p2, p3, p4, p5, p6, p7;
		//		p1.SetValue(0 + pi->position.x, 0 + pi->position.y, c + pi->position.z);
		//		p2.SetValue(0 + pi->position.x, 0 + pi->position.y, -c + pi->position.z);
		//		p3.SetValue(c + pi->position.x, 0 + pi->position.y, 0 + pi->position.z);
		//		p4.SetValue(-c + pi->position.x, 0 + pi->position.y, 0 + pi->position.z);
		//		p5.SetValue(0 + pi->position.x, -c + pi->position.y, 0 + pi->position.z);
		//		p6.SetValue(0 + pi->position.x, c + pi->position.y, 0 + pi->position.z);
		//		p7.SetValue(0 + pi->position.x, 0 + pi->position.y, 0 + pi->position.z);

		//		int hashValueP1 = FunctionKit::PositionMapHash(p1, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP2 = FunctionKit::PositionMapHash(p2, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP3 = FunctionKit::PositionMapHash(p3, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP4 = FunctionKit::PositionMapHash(p4, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP5 = FunctionKit::PositionMapHash(p5, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP6 = FunctionKit::PositionMapHash(p6, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP7 = FunctionKit::PositionMapHash(p7, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值

		//		np1->InitializationPVS(p1, hashValueP1, mass_np, pi->density);//需要添加构造函数（至少）：mass，position，velocity, density
		//		np2->InitializationPVS(p2, hashValueP2, mass_np, pi->density);
		//		np3->InitializationPVS(p3, hashValueP3, mass_np, pi->density);
		//		np4->InitializationPVS(p4, hashValueP4, mass_np, pi->density);
		//		np5->InitializationPVS(p5, hashValueP5, mass_np, pi->density);
		//		np6->InitializationPVS(p6, hashValueP6, mass_np, pi->density);
		//		np7->InitializationPVS(p7, hashValueP7, mass_np, pi->density);

		//		np1->β = 0.5;
		//		np2->β = 0.5;
		//		np3->β = 0.5;
		//		np4->β = 0.5;
		//		np5->β = 0.5;
		//		np6->β = 0.5;
		//		np7->β = 0.5;

		//		fluidModel.particleList.push_back(np1);
		//		fluidModel.particleList.push_back(np2);
		//		fluidModel.particleList.push_back(np3);
		//		fluidModel.particleList.push_back(np4);
		//		fluidModel.particleList.push_back(np5);
		//		fluidModel.particleList.push_back(np6);
		//		fluidModel.particleList.push_back(np7);


		//		//删除粒子部分代码开始标志--------------------------------------------------------------------------------------------
		//		//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
		//		{
		//			//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，新增其指向生成粒子的指针
		//			//aa
		//			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
		//			{
		//				if (pi->fluidNeighbors[j] == NULL)
		//				{
		//					continue;
		//				}

		//				IISPHParticle *pj = pi->fluidNeighbors[j];
		//				for (int k = 0; k < pj->fluidNeighbors.size(); k++)
		//				{
		//					if (pj->fluidNeighbors[k] == NULL)
		//					{
		//						continue;
		//					}

		//					IISPHParticle *pk = pj->fluidNeighbors[k];
		//					if (pk == pi)
		//					{
		//						pj->fluidNeighbors[k] = NULL;
		//						pj->fluidNeighbors.push_back(np1);
		//						pj->fluidNeighbors.push_back(np2);
		//						pj->fluidNeighbors.push_back(np3);
		//						pj->fluidNeighbors.push_back(np4);
		//						pj->fluidNeighbors.push_back(np5);
		//						pj->fluidNeighbors.push_back(np6);
		//						pj->fluidNeighbors.push_back(np7);

		//					}
		//				}

		//			}
		//			delete fluidModel.particleList[i];
		//			fluidModel.particleList[i] = NULL;
		//		}

		//		//删除粒子部分代码结束标志----------------------------------------------------------------------------------------------------
		//	}
		//	if (num_np == 8)
		//	{
		//		cout << "888888888888888888888888888888888888888888888888888888888" << endl;

		//		double r_Spliting = 0.3*pi->particleSupportRad;//新生成的粒子所在球面的半径，0.3为取得自定义的值，此处用来测试
		//		double mass_np = pi->particleMass / num_np;
		//		double c = 1.0 / sqrt(3.0)*r_Spliting;

		//		IISPHParticle *np1 = new IISPHParticle(*pi);
		//		IISPHParticle *np2 = new IISPHParticle(*pi);
		//		IISPHParticle *np3 = new IISPHParticle(*pi);
		//		IISPHParticle *np4 = new IISPHParticle(*pi);
		//		IISPHParticle *np5 = new IISPHParticle(*pi);
		//		IISPHParticle *np6 = new IISPHParticle(*pi);
		//		IISPHParticle *np7 = new IISPHParticle(*pi);
		//		IISPHParticle *np8 = new IISPHParticle(*pi);

		//		np1->lastPressure = pi->pressure;
		//		np2->lastPressure = pi->pressure;
		//		np3->lastPressure = pi->pressure;
		//		np4->lastPressure = pi->pressure;
		//		np5->lastPressure = pi->pressure;
		//		np6->lastPressure = pi->pressure;
		//		np7->lastPressure = pi->pressure;
		//		np8->lastPressure = pi->pressure;

		//		Vector3f p1, p2, p3, p4, p5, p6, p7, p8;
		//		p1.SetValue(c + pi->position.x, c + pi->position.y, c + pi->position.z);
		//		p2.SetValue(-c + pi->position.x, c + pi->position.y, c + pi->position.z);
		//		p3.SetValue(c + pi->position.x, -c + pi->position.y, c + pi->position.z);
		//		p4.SetValue(c + pi->position.x, c + pi->position.y, -c + pi->position.z);
		//		p5.SetValue(-c + pi->position.x, -c + pi->position.y, c + pi->position.z);
		//		p6.SetValue(-c + pi->position.x, c + pi->position.y, -c + pi->position.z);
		//		p7.SetValue(c + pi->position.x, -c + pi->position.y, -c + pi->position.z);
		//		p8.SetValue(-c + pi->position.x, -c + pi->position.y, -c + pi->position.z);

		//		int hashValueP1 = FunctionKit::PositionMapHash(p1, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP2 = FunctionKit::PositionMapHash(p2, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP3 = FunctionKit::PositionMapHash(p3, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP4 = FunctionKit::PositionMapHash(p4, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP5 = FunctionKit::PositionMapHash(p5, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP6 = FunctionKit::PositionMapHash(p6, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP7 = FunctionKit::PositionMapHash(p7, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP8 = FunctionKit::PositionMapHash(p8, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值

		//		np1->InitializationPVS(p1, hashValueP1, mass_np, pi->density);//需要添加构造函数（至少）：mass，position，velocity, density
		//		np2->InitializationPVS(p2, hashValueP2, mass_np, pi->density);
		//		np3->InitializationPVS(p3, hashValueP3, mass_np, pi->density);
		//		np4->InitializationPVS(p4, hashValueP4, mass_np, pi->density);
		//		np5->InitializationPVS(p5, hashValueP5, mass_np, pi->density);
		//		np6->InitializationPVS(p6, hashValueP6, mass_np, pi->density);
		//		np7->InitializationPVS(p7, hashValueP7, mass_np, pi->density);
		//		np8->InitializationPVS(p8, hashValueP8, mass_np, pi->density);

		//		np1->β = 0.5;
		//		np2->β = 0.5;
		//		np3->β = 0.5;
		//		np4->β = 0.5;
		//		np5->β = 0.5;
		//		np6->β = 0.5;
		//		np7->β = 0.5;
		//		np8->β = 0.5;

		//		fluidModel.particleList.push_back(np1);
		//		fluidModel.particleList.push_back(np2);
		//		fluidModel.particleList.push_back(np3);
		//		fluidModel.particleList.push_back(np4);
		//		fluidModel.particleList.push_back(np5);
		//		fluidModel.particleList.push_back(np6);
		//		fluidModel.particleList.push_back(np7);
		//		fluidModel.particleList.push_back(np8);


		//		//删除粒子部分代码开始标志--------------------------------------------------------------------------------------------
		//		//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
		//		{
		//			//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，新增其指向生成粒子的指针
		//			//aa
		//			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
		//			{
		//				if (pi->fluidNeighbors[j] == NULL)
		//				{
		//					continue;
		//				}

		//				IISPHParticle *pj = pi->fluidNeighbors[j];
		//				for (int k = 0; k < pj->fluidNeighbors.size(); k++)
		//				{
		//					if (pj->fluidNeighbors[k] == NULL)
		//					{
		//						continue;
		//					}

		//					IISPHParticle *pk = pj->fluidNeighbors[k];
		//					if (pk == pi)
		//					{
		//						pj->fluidNeighbors[k] = NULL;
		//						pj->fluidNeighbors.push_back(np1);
		//						pj->fluidNeighbors.push_back(np2);
		//						pj->fluidNeighbors.push_back(np3);
		//						pj->fluidNeighbors.push_back(np4);
		//						pj->fluidNeighbors.push_back(np5);
		//						pj->fluidNeighbors.push_back(np6);
		//						pj->fluidNeighbors.push_back(np7);
		//						pj->fluidNeighbors.push_back(np8);

		//					}
		//				}

		//			}
		//			delete fluidModel.particleList[i];
		//			fluidModel.particleList[i] = NULL;
		//		}

		//		//删除粒子部分代码结束标志----------------------------------------------------------------------------------------------------
		//	}

		//	if (num_np == 9)
		//	{
		//		cout << "99999999999999999999999999999999999999999999" << endl;

		//		double r_Spliting = 0.3*pi->particleSupportRad;//新生成的粒子所在球面的半径，0.3为取得自定义的值，此处用来测试
		//		double mass_np = pi->particleMass / num_np;
		//		double c = 1.0 / sqrt(3.0)*r_Spliting;

		//		IISPHParticle *np1 = new IISPHParticle(*pi);
		//		IISPHParticle *np2 = new IISPHParticle(*pi);
		//		IISPHParticle *np3 = new IISPHParticle(*pi);
		//		IISPHParticle *np4 = new IISPHParticle(*pi);
		//		IISPHParticle *np5 = new IISPHParticle(*pi);
		//		IISPHParticle *np6 = new IISPHParticle(*pi);
		//		IISPHParticle *np7 = new IISPHParticle(*pi);
		//		IISPHParticle *np8 = new IISPHParticle(*pi);
		//		IISPHParticle *np9 = new IISPHParticle(*pi);

		//		np1->lastPressure = pi->pressure;
		//		np2->lastPressure = pi->pressure;
		//		np3->lastPressure = pi->pressure;
		//		np4->lastPressure = pi->pressure;
		//		np5->lastPressure = pi->pressure;
		//		np6->lastPressure = pi->pressure;
		//		np7->lastPressure = pi->pressure;
		//		np8->lastPressure = pi->pressure;
		//		np9->lastPressure = pi->pressure;

		//		Vector3f p1, p2, p3, p4, p5, p6, p7, p8, p9;
		//		p1.SetValue(c + pi->position.x, c + pi->position.y, c + pi->position.z);
		//		p2.SetValue(-c + pi->position.x, c + pi->position.y, c + pi->position.z);
		//		p3.SetValue(c + pi->position.x, -c + pi->position.y, c + pi->position.z);
		//		p4.SetValue(c + pi->position.x, c + pi->position.y, -c + pi->position.z);
		//		p5.SetValue(-c + pi->position.x, -c + pi->position.y, c + pi->position.z);
		//		p6.SetValue(-c + pi->position.x, c + pi->position.y, -c + pi->position.z);
		//		p7.SetValue(c + pi->position.x, -c + pi->position.y, -c + pi->position.z);
		//		p8.SetValue(-c + pi->position.x, -c + pi->position.y, -c + pi->position.z);
		//		p9.SetValue(0 + pi->position.x, 0 + pi->position.y, 0 + pi->position.z);

		//		int hashValueP1 = FunctionKit::PositionMapHash(p1, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP2 = FunctionKit::PositionMapHash(p2, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP3 = FunctionKit::PositionMapHash(p3, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP4 = FunctionKit::PositionMapHash(p4, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP5 = FunctionKit::PositionMapHash(p5, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP6 = FunctionKit::PositionMapHash(p6, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP7 = FunctionKit::PositionMapHash(p7, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP8 = FunctionKit::PositionMapHash(p8, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值
		//		int hashValueP9 = FunctionKit::PositionMapHash(p9, SharedData::GetCellLength(), SharedData::GetHashCellNum(), boundary.offset); //得到当前position所对应的哈希值

		//		np1->InitializationPVS(p1, hashValueP1, mass_np, pi->density);//需要添加构造函数（至少）：mass，position，velocity, density
		//		np2->InitializationPVS(p2, hashValueP2, mass_np, pi->density);
		//		np3->InitializationPVS(p3, hashValueP3, mass_np, pi->density);
		//		np4->InitializationPVS(p4, hashValueP4, mass_np, pi->density);
		//		np5->InitializationPVS(p5, hashValueP5, mass_np, pi->density);
		//		np6->InitializationPVS(p6, hashValueP6, mass_np, pi->density);
		//		np7->InitializationPVS(p7, hashValueP7, mass_np, pi->density);
		//		np8->InitializationPVS(p8, hashValueP8, mass_np, pi->density);
		//		np9->InitializationPVS(p9, hashValueP9, mass_np, pi->density);

		//		np1->β = 0.5;
		//		np2->β = 0.5;
		//		np3->β = 0.5;
		//		np4->β = 0.5;
		//		np5->β = 0.5;
		//		np6->β = 0.5;
		//		np7->β = 0.5;
		//		np8->β = 0.5;
		//		np9->β = 0.5;

		//		fluidModel.particleList.push_back(np1);
		//		fluidModel.particleList.push_back(np2);
		//		fluidModel.particleList.push_back(np3);
		//		fluidModel.particleList.push_back(np4);
		//		fluidModel.particleList.push_back(np5);
		//		fluidModel.particleList.push_back(np6);
		//		fluidModel.particleList.push_back(np7);
		//		fluidModel.particleList.push_back(np8);
		//		fluidModel.particleList.push_back(np9);


		//		//删除粒子部分代码开始标志--------------------------------------------------------------------------------------------
		//		//删除i粒子，需要先删除所有指向i粒子的指针，1、粒子表的指针；2、其邻居指向的指针。然后删除i粒子，pi指向设为空
		//		{
		//			//遍历i粒子的j个邻居，对于邻居j（邻居是相互的）：j粒子的邻居中有i粒子的指针，遍历j粒子的邻居k，删除其指向i的指针，新增其指向生成粒子的指针
		//			//aa
		//			for (int j = 0; j < pi->fluidNeighbors.size(); j++)
		//			{
		//				if (pi->fluidNeighbors[j] == NULL)
		//				{
		//					continue;
		//				}

		//				IISPHParticle *pj = pi->fluidNeighbors[j];
		//				for (int k = 0; k < pj->fluidNeighbors.size(); k++)
		//				{
		//					if (pj->fluidNeighbors[k] == NULL)
		//					{
		//						continue;
		//					}

		//					IISPHParticle *pk = pj->fluidNeighbors[k];
		//					if (pk == pi)
		//					{
		//						pj->fluidNeighbors[k] = NULL;
		//						pj->fluidNeighbors.push_back(np1);
		//						pj->fluidNeighbors.push_back(np2);
		//						pj->fluidNeighbors.push_back(np3);
		//						pj->fluidNeighbors.push_back(np4);
		//						pj->fluidNeighbors.push_back(np5);
		//						pj->fluidNeighbors.push_back(np6);
		//						pj->fluidNeighbors.push_back(np7);
		//						pj->fluidNeighbors.push_back(np8);
		//						pj->fluidNeighbors.push_back(np9);

		//					}
		//				}

		//			}
		//			delete fluidModel.particleList[i];
		//			fluidModel.particleList[i] = NULL;
		//		}

		//		//删除粒子部分代码结束标志----------------------------------------------------------------------------------------------------
		//	}

		//}



	}

	//
	//	//对L粒子的操作：


	//	//*****************************************************************************************************************
	//	//论文中粒子分裂方法：
	//	Particle &p = particles[0];
	//	int n = 9;
	//
	//	//粒子分裂(n=2)
	//	if (n == 2)
	//	{
	//		double r_Spliting = 0.3*4.0*p.r;
	//
	//		double mass = p.mass / double(n);
	//		Particle np1(p.pos.x, p.pos.y, r_Spliting + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np2(p.pos.x, p.pos.y, -r_Spliting + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//
	//		particles.push_back(np1);
	//		particles.push_back(np2);
	//
	//		particles.erase(particles.begin() + 0);
	//	}
	//
	//	//粒子分裂(n=3)
	//	if (n == 3)
	//	{
	//		double r_Spliting = 0.3*4.0*p.r;
	//
	//		double c = sqrt(3.0)*r_Spliting*0.5;
	//		double mass = p.mass / double(n);
	//		Particle np1(p.pos.x, r_Spliting + p.pos.y, p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np2(-c + p.pos.x, -0.5*r_Spliting + p.pos.y, p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np3(c + p.pos.x, -0.5*r_Spliting + p.pos.y, p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//
	//		particles.push_back(np1);
	//		particles.push_back(np2);
	//		particles.push_back(np3);
	//
	//		particles.erase(particles.begin() + 0);
	//	}
	//
	//
	//	//粒子分裂(n=4)
	//	if (n == 4)
	//	{
	//		double r_Spliting = 0.3*4.0*p.r;
	//
	//		double c = 1.0 / sqrt(3.0)*r_Spliting;
	//		double mass = p.mass / double(n);
	//		Particle np1(c + p.pos.x, c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np2(-c + p.pos.x, c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np3(-c + p.pos.x, -c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np4(c + p.pos.x, -c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//
	//		particles.push_back(np1);
	//		particles.push_back(np2);
	//		particles.push_back(np3);
	//		particles.push_back(np4);
	//
	//		particles.erase(particles.begin() + 0);
	//	}
	//
	//
	//	//粒子分裂(n=5)
	//	if (n == 5)
	//	{
	//		double r_Spliting = 0.3*4.0*p.r;
	//
	//		double c = 1.0 / sqrt(3.0)*r_Spliting;
	//		double mass = p.mass / double(n);
	//		Particle np0(p.pos.x, p.pos.y, p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np1(c + p.pos.x, c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np2(-c + p.pos.x, c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np3(-c + p.pos.x, -c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np4(c + p.pos.x, -c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		particles.push_back(np0);
	//		particles.push_back(np1);
	//		particles.push_back(np2);
	//		particles.push_back(np3);
	//		particles.push_back(np4);
	//
	//		particles.erase(particles.begin() + 0);
	//
	//	}
	//
	//	//粒子分裂(n=6)
	//	if (n == 6)
	//	{
	//		double r_Spliting = 0.3*4.0*p.r;
	//
	//		double c = r_Spliting;
	//		double mass = p.mass / double(n);
	//		Particle np1(0 + p.pos.x, 0 + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np2(0 + p.pos.x, 0 + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np3(c + p.pos.x, 0 + p.pos.y, 0 + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np4(-c + p.pos.x, 0 + p.pos.y, 0 + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np5(0 + p.pos.x, -c + p.pos.y, 0 + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np6(0 + p.pos.x, c + p.pos.y, 0 + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//
	//		particles.push_back(np1);
	//		particles.push_back(np2);
	//		particles.push_back(np3);
	//		particles.push_back(np4);
	//		particles.push_back(np5);
	//		particles.push_back(np6);
	//
	//		particles.erase(particles.begin() + 0);
	//	}
	//
	//
	//	//粒子分裂(n=7)
	//	if (n == 7)
	//	{
	//		double r_Spliting = 0.3*4.0*p.r;
	//
	//		double c = r_Spliting;
	//		double mass = p.mass / double(n);
	//		Particle np0(p.pos.x, p.pos.y, p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np1(0 + p.pos.x, 0 + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np2(0 + p.pos.x, 0 + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np3(c + p.pos.x, 0 + p.pos.y, 0 + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np4(-c + p.pos.x, 0 + p.pos.y, 0 + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np5(0 + p.pos.x, -c + p.pos.y, 0 + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np6(0 + p.pos.x, c + p.pos.y, 0 + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		particles.push_back(np0);
	//		particles.push_back(np1);
	//		particles.push_back(np2);
	//		particles.push_back(np3);
	//		particles.push_back(np4);
	//		particles.push_back(np5);
	//		particles.push_back(np6);
	//
	//		particles.erase(particles.begin() + 0);
	//	}
	//
	//
	//	//粒子分裂(n=8)
	//
	//	if (n == 8)
	//	{
	//		double r_Spliting = 0.3*4.0*p.r;
	//
	//		double c = 1.0 / sqrt(3.0)*r_Spliting;
	//		double mass = p.mass / double(n);
	//		Particle np1(c + p.pos.x, c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np2(-c + p.pos.x, c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np3(c + p.pos.x, -c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np4(c + p.pos.x, c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np5(-c + p.pos.x, -c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np6(-c + p.pos.x, c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np7(c + p.pos.x, -c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np8(-c + p.pos.x, -c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//
	//		particles.push_back(np1);
	//		particles.push_back(np2);
	//		particles.push_back(np3);
	//		particles.push_back(np4);
	//		particles.push_back(np5);
	//		particles.push_back(np6);
	//		particles.push_back(np7);
	//		particles.push_back(np8);
	//
	//		particles.erase(particles.begin() + 0);
	//	}
	//
	//
	//	//粒子分裂(n=9)
	//	if (n == 9)
	//	{
	//		double r_Spliting = 0.3*4.0*p.r;
	//
	//		double c = 1.0 / sqrt(3.0)*r_Spliting;
	//		double mass = p.mass / double(n);
	//		Particle np0(p.pos.x, p.pos.y, p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np1(c + p.pos.x, c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np2(-c + p.pos.x, c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np3(c + p.pos.x, -c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np4(c + p.pos.x, c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np5(-c + p.pos.x, -c + p.pos.y, c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np6(-c + p.pos.x, c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np7(c + p.pos.x, -c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		Particle np8(-c + p.pos.x, -c + p.pos.y, -c + p.pos.z, mass, FLUID_PARTICLE, Point3D(0, 0, 0));
	//		particles.push_back(np0);
	//		particles.push_back(np1);
	//		particles.push_back(np2);
	//		particles.push_back(np3);
	//		particles.push_back(np4);
	//		particles.push_back(np5);
	//		particles.push_back(np6);
	//		particles.push_back(np7);
	//		particles.push_back(np8);
	//
	//		particles.erase(particles.begin() + 0);
	//	}
	//
	//	//***********************************************************************************************************
	//
	//
	//	int num0 = particles.size();
	//	cout << "改变之后粒子的数量：" << num0 << endl;
}