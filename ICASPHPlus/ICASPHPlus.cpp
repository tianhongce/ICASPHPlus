#include "ICASPHPlus.h"

ICASPHPlus::ICASPHPlus(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	init();
}
ICASPHPlus::~ICASPHPlus()
{

}

void ICASPHPlus::init()
{
	glwidget = new GLWidget(&ui);
	ui.BaseLayout->addWidget(glwidget);

}
void ICASPHPlus::onStart()
{
	ui.actionStart->setChecked(true);
	ui.actionStop->setChecked(false);
	glwidget->StartRun();

}

void ICASPHPlus::onStop()
{
	ui.actionStop->setChecked(true);
	ui.actionStart->setChecked(false);
	glwidget->StopRun();
}

void ICASPHPlus::isSavePic()
{
	if (ui.actionSavePic->isChecked() == true)
	{
		glwidget->setOutputFramePicture(true);
	}
	else
	{
		glwidget->setOutputFramePicture(false);
	}
	
}
void ICASPHPlus::changeIISPH()
{
	//��ֹͣ
	glwidget->StopRun();
	//����ΪIISPH����Ⱦ��ʽ
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//��ʼ��
	glwidget->sphManager.Initialize();
	//ˢ�½���
	glwidget->update();
	ui.SPHTypeName->setText("IISPH");
}
void ICASPHPlus::changeWCSPH()
{
	//��ֹͣ
	glwidget->StopRun();
	//����ΪIISPH����Ⱦ��ʽ
	glwidget->sphManager.SetSPHType(SPH_TYPE::WCSPH);
	//��ʼ��
	glwidget->sphManager.Initialize();
	//ˢ�½���
	glwidget->update();
	ui.SPHTypeName->setText("WCSPH");
}
void ICASPHPlus::isShowBoundaryPartilces()
{
	if (ui.actionshowBoundaryPar->isChecked() == true)
	{
		glwidget->setShowBoundaryParticles(true);
	}
	else
	{
		glwidget->setShowBoundaryParticles(false);
	}
}