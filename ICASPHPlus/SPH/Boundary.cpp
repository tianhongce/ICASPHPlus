#include "Boundary.h"

Boundary::Boundary()
{
	GenerateNewBoundary(boundaryOriginX, boundaryOriginY, boundaryOriginZ, boundaryWidth, boundaryLength, boundaryHeight, offset);
}

Boundary::Boundary(Boundary & a)
{
	operator=(a);
}

Boundary::Boundary(double x, double y, double z, double w, double l, double h, double os)
{
	GenerateNewBoundary(x, y, z, w, l, h, os);

}

Boundary::~Boundary()
{
}

void Boundary::GenerateNewBoundary(double x, double y, double z, double w, double l, double h, double os)
{
	boundaryOriginX = x;
	boundaryOriginY = y;
	boundaryOriginZ = z;
	boundaryWidth = w;
	boundaryLength = l;
	boundaryHeight = h;
	offset = os;
	gridBoundaryOriginX = boundaryOriginX - offset;
	gridBoundaryOriginY = boundaryOriginY - offset;
	gridBoundaryOriginZ = boundaryOriginZ - offset;
	gridBoundaryWidth = boundaryWidth + 2 * offset;
	gridBoundaryLength = boundaryLength + 2 * offset;
	gridBoundaryHeight = boundaryHeight + 2 * offset;

	//������ֵ
	Vector3f origin(x, y, z);
	Vector3f wlh(w, l, h);
	//����ʼ��һ��boundaryʱ��������SharedData��������
	SharedData::SetBoundaryOrigin(origin);
	SharedData::SetBoundaryWLH(wlh);
	SharedData::SetOffset(os);
}

void Boundary::GenerateNewBoundary(Vector3f position, Vector3f fluidWLH, double offs)
{
	boundaryOriginX = position.GetX();
	boundaryOriginY = position.GetX();
	boundaryOriginZ = position.GetZ();
	boundaryWidth = fluidWLH.GetX();
	boundaryLength = fluidWLH.GetY();
	boundaryHeight = fluidWLH.GetZ();
	offset = offs;
	gridBoundaryOriginX = boundaryOriginX - offset;
	gridBoundaryOriginY = boundaryOriginY - offset;
	gridBoundaryOriginZ = boundaryOriginZ - offset;
	gridBoundaryWidth = boundaryWidth + 2 * offset;
	gridBoundaryLength = boundaryLength + 2 * offset;
	gridBoundaryHeight = boundaryHeight + 2 * offset;
}

Boundary & Boundary::operator=(Boundary & a)
{
	boundaryOriginX = a.boundaryOriginX;
	boundaryOriginY = a.boundaryOriginY;
	boundaryOriginZ = a.boundaryOriginZ;
	boundaryWidth = a.boundaryWidth;
	boundaryLength = a.boundaryLength;
	boundaryHeight = a.boundaryHeight;
	offset = a.offset;
	gridBoundaryOriginX = boundaryOriginX - offset;
	gridBoundaryOriginY = boundaryOriginY - offset;
	gridBoundaryOriginZ = boundaryOriginZ - offset;
	gridBoundaryWidth = boundaryWidth + 2 * offset;
	gridBoundaryLength = boundaryLength + 2 * offset;
	gridBoundaryHeight = boundaryHeight + 2 * offset;
	pRad = a.pRad;
	return *this;
}

//���ݱ߽�ľ�����Ϣ����������ǽ�������ÿ�����ӵ�λ��
vector<Vector3f> Boundary::initBoundaryData()
{
	std::cout << "��ʼ���߽�����" << std::endl;
	std::vector<Vector3f> boundaryParticles;
	//StaticRigidObject *rb = new StaticRigidObject();
	const double particleRad = pRad;   //Ӳ�Թ涨Ϊ��С�뾶
	const double diam = particleRad * 2; //�õ���ǰ���ӵ�ֱ��

	Vector3f boundaryMinP(boundaryOriginX-particleRad, boundaryOriginY- particleRad, boundaryOriginZ- particleRad);
	Vector3f boundaryMaxP(boundaryOriginX + boundaryWidth + particleRad, boundaryOriginY + boundaryLength + particleRad, boundaryOriginZ + boundaryHeight + particleRad);

	//ʹ��ֱ������Ϊÿһ����ǰ�ı仯��
	double xshift = diam;
	double yshift = diam;
	double zshift = diam;

	const int stepsX = (int)round((boundaryWidth + 2 * diam) / xshift);
	const int stepsY = (int)round((boundaryLength + 2 * diam) / yshift);
	const int stepsZ = (int)round((boundaryHeight + 2 * diam) / zshift);

	Vector3f start = boundaryMinP;

	//����ǽ��ʹ�����ӽ���ǽ�Ĵ���
	for (int i = 0; i < stepsX; i++)
	{
		for (int j = 0; j < stepsY; j++)
		{
			Vector3f currPos = Vector3f(i*xshift, j*yshift, 0) + boundaryMinP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = Vector3f(i*xshift, 0, k*zshift) + boundaryMinP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int j = 0; j < stepsY; j++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = Vector3f(0, j*yshift, k*zshift) + boundaryMinP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int j = 0; j < stepsY; j++)
		{
			Vector3f currPos = -Vector3f(i*xshift, j*yshift, 0) + boundaryMaxP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = -Vector3f(i*xshift, 0, k*zshift) + boundaryMaxP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int j = 0; j < stepsY; j++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = -Vector3f(0, j*yshift, k*zshift) + boundaryMaxP;
			boundaryParticles.push_back(currPos);
		}
	}
	//std::cout << "boundaryParticles.size():" << boundaryParticles.size() << std::endl;
	//model.addRigidBodyObject(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
	std::cout << "��ʼ���߽��������" << std::endl;
	return boundaryParticles;
}
