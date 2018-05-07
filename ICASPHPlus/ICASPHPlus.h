#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_ICASPHPlus.h"
#include "glwidget.h"

class ICASPHPlus : public QMainWindow
{
	Q_OBJECT

public:
	ICASPHPlus(QWidget *parent = 0);
	~ICASPHPlus();

private:
	Ui::ICASPHPlusClass ui;

	void init();
	GLWidget *glwidget;
private slots :
	void onStart();
	void onStop();
	void isSavePic();
	void changeIISPH();
	void changeWCSPH();
	void isShowBoundaryPartilces();
};
