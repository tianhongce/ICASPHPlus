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
	//先停止
	glwidget->StopRun();
	//设置为IISPH的渲染方式
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//初始化
	glwidget->sphManager.Initialize();
	//刷新界面
	glwidget->update();
	ui.SPHTypeName->setText("IISPH");
}
void ICASPHPlus::changeWCSPH()
{
	//先停止
	glwidget->StopRun();
	//设置为IISPH的渲染方式
	glwidget->sphManager.SetSPHType(SPH_TYPE::WCSPH);
	//初始化
	glwidget->sphManager.Initialize();
	//刷新界面
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