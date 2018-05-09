#include "IISPHFluidObject.h"

void IISPHFluidObject::ClearNeighborList()
{
    #pragma omp parallel default(shared)
	{
       #pragma omp for schedule(static) 
		//ע����������ǰ�ھ��б�����б�
		for (int i = 0; i < IISPH_ParticleNum; i++)
		{
			if (particleList[i] == NULL)
			{
				continue;
			}
			particleList[i]->fluidNeighbors.clear();
		}
	}
}

void IISPHFluidObject::ClearFluidBoundaryNeighborList()
{
    #pragma omp parallel default(shared)
	{
       #pragma omp for schedule(static) 
		for (int i = 0; i < IISPH_ParticleNum; i++)
		{
			if (particleList[i] == NULL)
			{
				continue;
			}
			particleList[i]->boundaryNeighbors.clear();
		}
	}
}

void IISPHFluidObject::InitializeBoundary()
{
	//���ݱ߽����ӣ������Ӧ�ı߽����ӵ�λ������
	vector<Vector3f> boundaryParticlePosition = boundary.initBoundaryData();
	//���ݼ���õ��ı߽����ӵĸ����Ա߽��ھ������б��������
	boundaryNeighborList.resize(boundaryParticlePosition.size());
	//��ʼ����̬���������
	StaticRigidObject rb;
	//ʹ�þ�̬������߽����ӵ�λ�öԾ�̬���ӱ߽���г�ʼ��
	AddRigidBodyObject(rb, boundaryParticlePosition.size(), boundaryParticlePosition);

	//������������SPHComputation������
	//���ݾ�̬���ӵ�λ�ã�������ӳ�䵽������
	//����ÿһ����̬�߽����ӵ�BoundaryPSI
}

void IISPHFluidObject::AddRigidBodyObject(StaticRigidObject & rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles)
{
	//�ӹ�ϣ�����б��еõ������ؼ��ľ�ֵ̬,���õĶ��Ǿ�ֵ̬
	int cellNum = SharedData::GetHashCellNum();
	double cellLength = SharedData::GetCellLength();

	std::cout << "��Ӹ������" << std::endl;
	boundaryObj.staticRigidParticleList.resize(numBoundaryParticles);    //��ʼ����ǰ��̬���������б�Ĵ�С
	//������г�ʼ��
	//ÿһ��ָ�붼ָ��һ����̬��������
	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		boundaryObj.staticRigidParticleList[i] = new StaticRigidParticle();
	}

	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		boundaryObj.staticRigidParticleList[i]->x0 = boundaryParticles[i];
		boundaryObj.staticRigidParticleList[i]->position = boundaryParticles[i];
		boundaryObj.staticRigidParticleList[i]->velocity.Zero();
		boundaryObj.staticRigidParticleList[i]->m_f.Zero();
		//���ݵ�ǰλ�ü���õ���Ӧ���ӵĹ�ϣֵ
		boundaryObj.staticRigidParticleList[i]->hashValue =
			FunctionKit::PositionMapHash(boundaryParticles[i], cellLength,
				cellNum, boundary.offset);
	}
	//�Ա߽����ĸ������Խ��и�ֵ
	boundaryObj.rd = rbo;
}

void IISPHFluidObject::ComputeBoundaryPsi(RigidBodyParticleObject & bound)
{
	double radius = boundary.pRad * 4;
	const double density0 = IISPH_RestDensity; //�õ������ʼ�ܶ�
	//�õ���̬���Ӹ���
	const unsigned int numBoundaryParticles = boundaryObj.GetStaticParticleNum();
	//������ǰ�ı߽�����
	for (int i = 0; i < (int)numBoundaryParticles; i++)
	{
		double delta = CubicKernelOne::W_zero(radius);
		//�õ���Ӧ������bound���еľ�̬��������
		StaticRigidParticle* a = boundaryObj.staticRigidParticleList[i];
		for (unsigned int j = 0; j < boundaryNeighborList[i].size(); j++)
		{
			//�õ���Ӧ�ٽ����ӵ�ǰ���ٽ���
			StaticRigidParticle* b = boundaryNeighborList[i][j];
			//����õ���ǰ������������õ��Ĳ�ֵ
			Vector3f positionSub = a->position - b->position;
			delta += CubicKernelOne::W(radius, positionSub);
		}
		const double volume = 1.0 / delta;
		//�Ե�ǰ��̬���ӵ�boundaryPsi���и�ֵ
		boundaryObj.staticRigidParticleList[i]->boundaryPsi = density0 * volume;
		//std::cout << "  BoundaryPSI : " << boundaryObj.staticRigidParticleList[i]->boundaryPsi << std::endl;
	}
}

void IISPHFluidObject::Initialise()
{
	ComputeFluidWLH();
	ComputeParticleNum(FluidWLH.GetX(), FluidWLH.GetY(), FluidWLH.GetZ()); 
	InitialiseParticlePosition(); 
	InitialiseParticleDensity();
	InitializeBoundary();

	//���ݳ�ʼ�������Ӱ뾶ȷ����ǰ��֧�Ű뾶
	Poly6Kernel::setRadius(4 * particleInitialPad);
	SpikyKernel::setRadius(4 * particleInitialPad);
	CubicKernel::setRadius(4 * particleInitialPad);
}

void IISPHFluidObject::ComputeFluidWLH()
{
	int numX = floor(double(IISPH_FluidWitdth / (2 * particleInitialPad)));
	int numY = floor(double(IISPH_FluidLength / (2 * particleInitialPad)));
	int numZ = floor(double(IISPH_FluidHeight / (2 * particleInitialPad)));
 	FluidWLH.SetValue(numX, numY, numZ);
}

void IISPHFluidObject::InitialiseParticlePosition()
{
	//��ʼ�����ӵ�λ��
	//�����б��е����
	int particleIndex = 0;
	//���ӵ�ֱ��
	double particleLength = particleInitialPad * 2;
	for (int i = 0; i < FluidWLH.GetX(); i++)
	{
		for (int j = 0; j < FluidWLH.GetY(); j++)
		{
			for (int k = 0; k < FluidWLH.GetZ(); k++)
			{
				Vector3f position(0, 0, 0); //��ʼ����ǰ���������ӵ�λ��
				position.SetX(IISPH_OriginX + particleInitialPad + particleLength*i); //���õ�ǰ���ӵ�x����
				position.SetY(IISPH_OriginY + particleInitialPad + particleLength*j); //���õ�ǰ���ӵ�y����
				position.SetZ(IISPH_OriginZ + particleInitialPad + particleLength*k); //���õ�ǰ���ӵ�z����
				int hashValue = FunctionKit::PositionMapHash(position,
					SharedData::GetCellLength(),
					SharedData::GetHashCellNum(),
					boundary.offset); //�õ���ǰposition����Ӧ�Ĺ�ϣֵ
				IISPHParticle* newParticle=new IISPHParticle(); //��ʼ��һ���µ�����
				newParticle->Initialization(position, hashValue, particleInitialPad); //�������ӽ��г�ʼ��
                //����ʼ���õ�������װ�������б���
				particleList[particleIndex] = newParticle;
				particleIndex++;
			}
		}
	}
}

void IISPHFluidObject::InitialiseParticleDensity()
{
	//��ˮ���������ӵ��ܶ�����Ϊˮ��ĳ�ʼ�ܶ�
	for (int i = 0; i < particleList.size(); i++)
	{
		if (particleList[i] == NULL)
		{
			continue;
		}
		particleList[i]->density = IISPH_RestDensity;
		//����������
		particleList[i]->lastDensity = IISPH_RestDensity;
	}
}

void IISPHFluidObject::UpdateParticleHashValue()
{
	//ˢ�µ�ǰ���ӵĹ�ϣֵ
    #pragma omp parallel default(shared)
	{
        #pragma omp for schedule(static) 
		for (int i = 0; i < particleList.size(); i++)
		{
			if (particleList[i] == NULL)
			{
				continue;
			}
			double cellLength = SharedData::GetCellLength(); //�õ�����Ĵ�С
			int hashCellNum = SharedData::GetHashCellNum();  //�õ����й�ϣ��������
			int hashValue = FunctionKit::PositionMapHash(particleList[i]->position, cellLength, hashCellNum, boundary.offset);
			particleList[i]->hashIndex = hashValue;
		}
	}
}

void IISPHFluidObject::ComputeParticleNum(int x, int y, int z)
{
	IISPH_ParticleNum = x*y*z;               //����x,y,z�᲻ͬ�ķ�����������������
	particleList.resize(IISPH_ParticleNum);  //����������������ˢ�����������ھ��б�����
	//neighborList.resize(IISPH_ParticleNum);
	//fluidBoundaryNeighborList.resize(IISPH_ParticleNum);
}
