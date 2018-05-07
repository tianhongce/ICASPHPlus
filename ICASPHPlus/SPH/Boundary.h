#pragma once

#include "Vector3f.h"
#include "Vector3i.h"
#include "StaticRigidObject.h"
#include "ShareData.h"

//////////////////////////////////////////////////////////////////////////
//Boundary��
//���ã�����һЩ�߽����ԡ����ݱ߽����Բ����õ���ǰ�߽������λ������
//////////////////////////////////////////////////////////////////////////
class Boundary
{
public:
	//�߽������
	double boundaryOriginX = 0;
	double boundaryOriginY = 0;
	double boundaryOriginZ = 0; //���ڱ߽���Եı߽���ʼ��
	double boundaryWidth = 2;
	double boundaryLength = 4;
	double boundaryHeight = 2;  //���ڱ߽���Եı߽糤���
	//������ʷǷ��ڴ�Ľ������
	double  offset =  0.6f;               //������������Ĳ��ֵ�ƫ�ƣ�һ��Ϊ��������ĳ���
	double gridBoundaryOriginX = -0.5f;
	double gridBoundaryOriginY = -0.5f;
	double gridBoundaryOriginZ = -0.5f; //������������ı߽���ʼ��
	double gridBoundaryWidth = 3.0f;
	double gridBoundaryLength = 5.0f;
	double gridBoundaryHeight = 3.0f;   //������������ĳ����

	double pRad = 0.075f;
public:
	//������������
	Boundary();
	Boundary(Boundary& a);
	Boundary(double x, double y, double z, double w, double l, double h, double os);
	~Boundary();
	//�����ض���ֵ�����µ���������
	void GenerateNewBoundary(double x, double y, double z, double w, double l, double h, double os); //�����µ�����
	void GenerateNewBoundary(Vector3f position, Vector3f fluidWLH, double offs);
	Boundary& operator=(Boundary& a);
	//���ݲ����õ��ı߽����ݲ�������õ���ǰ�ı߽���������
	vector<Vector3f> initBoundaryData();
};