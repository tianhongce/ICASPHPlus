/********************************************************************************
** Form generated from reading UI file 'ICASPHPlus.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ICASPHPLUS_H
#define UI_ICASPHPLUS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ICASPHPlusClass
{
public:
    QAction *actionStart;
    QAction *actionStop;
    QAction *actionSavePic;
    QAction *actionshowBoundaryPar;
    QAction *actionIISPH;
    QAction *actionWCSPH;
    QWidget *centralWidget;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *BaseLayout;
    QLabel *SPHTypeName;
    QFrame *line;
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout_down;
    QLabel *label_width;
    QDoubleSpinBox *SpinBox_width;
    QLabel *label_length;
    QDoubleSpinBox *SpinBox_length;
    QLabel *label_height;
    QDoubleSpinBox *SpinBox_height;
    QWidget *horizontalLayoutWidget_2;
    QHBoxLayout *horizontalLayout_top;
    QLabel *label_radius;
    QDoubleSpinBox *RadiusSpinBox;
    QMenuBar *menuBar;
    QMenu *menuSPH;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *ICASPHPlusClass)
    {
        if (ICASPHPlusClass->objectName().isEmpty())
            ICASPHPlusClass->setObjectName(QStringLiteral("ICASPHPlusClass"));
        ICASPHPlusClass->resize(1060, 795);
        actionStart = new QAction(ICASPHPlusClass);
        actionStart->setObjectName(QStringLiteral("actionStart"));
        actionStart->setCheckable(true);
        QIcon icon;
        icon.addFile(QStringLiteral("Resources/start.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionStart->setIcon(icon);
        actionStop = new QAction(ICASPHPlusClass);
        actionStop->setObjectName(QStringLiteral("actionStop"));
        actionStop->setCheckable(true);
        QIcon icon1;
        icon1.addFile(QStringLiteral("Resources/stop.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionStop->setIcon(icon1);
        actionSavePic = new QAction(ICASPHPlusClass);
        actionSavePic->setObjectName(QStringLiteral("actionSavePic"));
        actionSavePic->setCheckable(true);
        QIcon icon2;
        icon2.addFile(QStringLiteral("Resources/savepic.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSavePic->setIcon(icon2);
        actionshowBoundaryPar = new QAction(ICASPHPlusClass);
        actionshowBoundaryPar->setObjectName(QStringLiteral("actionshowBoundaryPar"));
        actionshowBoundaryPar->setCheckable(true);
        QIcon icon3;
        icon3.addFile(QStringLiteral("Resources/showBoundaryP.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionshowBoundaryPar->setIcon(icon3);
        actionIISPH = new QAction(ICASPHPlusClass);
        actionIISPH->setObjectName(QStringLiteral("actionIISPH"));
        actionWCSPH = new QAction(ICASPHPlusClass);
        actionWCSPH->setObjectName(QStringLiteral("actionWCSPH"));
        centralWidget = new QWidget(ICASPHPlusClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        verticalLayoutWidget = new QWidget(centralWidget);
        verticalLayoutWidget->setObjectName(QStringLiteral("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(0, 80, 1061, 661));
        BaseLayout = new QVBoxLayout(verticalLayoutWidget);
        BaseLayout->setSpacing(6);
        BaseLayout->setContentsMargins(11, 11, 11, 11);
        BaseLayout->setObjectName(QStringLiteral("BaseLayout"));
        BaseLayout->setContentsMargins(0, 0, 0, 0);
        SPHTypeName = new QLabel(centralWidget);
        SPHTypeName->setObjectName(QStringLiteral("SPHTypeName"));
        SPHTypeName->setGeometry(QRect(20, 10, 131, 51));
        QFont font;
        font.setFamily(QStringLiteral("Segoe Print"));
        font.setPointSize(22);
        SPHTypeName->setFont(font);
        line = new QFrame(centralWidget);
        line->setObjectName(QStringLiteral("line"));
        line->setGeometry(QRect(220, 0, 20, 81));
        line->setFrameShape(QFrame::VLine);
        line->setFrameShadow(QFrame::Sunken);
        horizontalLayoutWidget = new QWidget(centralWidget);
        horizontalLayoutWidget->setObjectName(QStringLiteral("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(230, 40, 821, 41));
        horizontalLayout_down = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout_down->setSpacing(6);
        horizontalLayout_down->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_down->setObjectName(QStringLiteral("horizontalLayout_down"));
        horizontalLayout_down->setContentsMargins(0, 0, 0, 0);
        label_width = new QLabel(horizontalLayoutWidget);
        label_width->setObjectName(QStringLiteral("label_width"));
        QFont font1;
        font1.setFamily(QStringLiteral("Segoe Print"));
        font1.setPointSize(16);
        label_width->setFont(font1);

        horizontalLayout_down->addWidget(label_width);

        SpinBox_width = new QDoubleSpinBox(horizontalLayoutWidget);
        SpinBox_width->setObjectName(QStringLiteral("SpinBox_width"));
        SpinBox_width->setDecimals(3);
        SpinBox_width->setMinimum(0.1);
        SpinBox_width->setMaximum(1.5);
        SpinBox_width->setSingleStep(0.025);
        SpinBox_width->setValue(0.5);

        horizontalLayout_down->addWidget(SpinBox_width);

        label_length = new QLabel(horizontalLayoutWidget);
        label_length->setObjectName(QStringLiteral("label_length"));
        label_length->setFont(font1);

        horizontalLayout_down->addWidget(label_length);

        SpinBox_length = new QDoubleSpinBox(horizontalLayoutWidget);
        SpinBox_length->setObjectName(QStringLiteral("SpinBox_length"));
        SpinBox_length->setDecimals(3);
        SpinBox_length->setMinimum(0.1);
        SpinBox_length->setMaximum(1.5);
        SpinBox_length->setSingleStep(0.025);
        SpinBox_length->setValue(0.5);

        horizontalLayout_down->addWidget(SpinBox_length);

        label_height = new QLabel(horizontalLayoutWidget);
        label_height->setObjectName(QStringLiteral("label_height"));
        label_height->setFont(font1);

        horizontalLayout_down->addWidget(label_height);

        SpinBox_height = new QDoubleSpinBox(horizontalLayoutWidget);
        SpinBox_height->setObjectName(QStringLiteral("SpinBox_height"));
        SpinBox_height->setDecimals(3);
        SpinBox_height->setMinimum(0.1);
        SpinBox_height->setMaximum(1.5);
        SpinBox_height->setSingleStep(0.025);
        SpinBox_height->setValue(0.5);

        horizontalLayout_down->addWidget(SpinBox_height);

        horizontalLayoutWidget_2 = new QWidget(centralWidget);
        horizontalLayoutWidget_2->setObjectName(QStringLiteral("horizontalLayoutWidget_2"));
        horizontalLayoutWidget_2->setGeometry(QRect(230, 0, 821, 41));
        horizontalLayout_top = new QHBoxLayout(horizontalLayoutWidget_2);
        horizontalLayout_top->setSpacing(6);
        horizontalLayout_top->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_top->setObjectName(QStringLiteral("horizontalLayout_top"));
        horizontalLayout_top->setContentsMargins(0, 0, 0, 0);
        label_radius = new QLabel(horizontalLayoutWidget_2);
        label_radius->setObjectName(QStringLiteral("label_radius"));
        label_radius->setFont(font1);

        horizontalLayout_top->addWidget(label_radius);

        RadiusSpinBox = new QDoubleSpinBox(horizontalLayoutWidget_2);
        RadiusSpinBox->setObjectName(QStringLiteral("RadiusSpinBox"));
        RadiusSpinBox->setDecimals(3);
        RadiusSpinBox->setMinimum(0.025);
        RadiusSpinBox->setMaximum(0.25);
        RadiusSpinBox->setSingleStep(0.001);
        RadiusSpinBox->setValue(0.025);

        horizontalLayout_top->addWidget(RadiusSpinBox);

        ICASPHPlusClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(ICASPHPlusClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1060, 23));
        menuSPH = new QMenu(menuBar);
        menuSPH->setObjectName(QStringLiteral("menuSPH"));
        ICASPHPlusClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(ICASPHPlusClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        ICASPHPlusClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(ICASPHPlusClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        ICASPHPlusClass->setStatusBar(statusBar);

        menuBar->addAction(menuSPH->menuAction());
        menuSPH->addAction(actionIISPH);
        menuSPH->addAction(actionWCSPH);
        mainToolBar->addAction(actionStart);
        mainToolBar->addAction(actionStop);
        mainToolBar->addSeparator();
        mainToolBar->addAction(actionSavePic);
        mainToolBar->addAction(actionshowBoundaryPar);

        retranslateUi(ICASPHPlusClass);
        QObject::connect(actionStart, SIGNAL(triggered()), ICASPHPlusClass, SLOT(onStart()));
        QObject::connect(actionStop, SIGNAL(triggered()), ICASPHPlusClass, SLOT(onStop()));
        QObject::connect(actionSavePic, SIGNAL(triggered()), ICASPHPlusClass, SLOT(isSavePic()));
        QObject::connect(actionshowBoundaryPar, SIGNAL(triggered()), ICASPHPlusClass, SLOT(isShowBoundaryPartilces()));

        QMetaObject::connectSlotsByName(ICASPHPlusClass);
    } // setupUi

    void retranslateUi(QMainWindow *ICASPHPlusClass)
    {
        ICASPHPlusClass->setWindowTitle(QApplication::translate("ICASPHPlusClass", "ICASPHPlus", Q_NULLPTR));
        actionStart->setText(QApplication::translate("ICASPHPlusClass", "Start", Q_NULLPTR));
        actionStop->setText(QApplication::translate("ICASPHPlusClass", "Stop", Q_NULLPTR));
        actionSavePic->setText(QApplication::translate("ICASPHPlusClass", "SavePic", Q_NULLPTR));
        actionshowBoundaryPar->setText(QApplication::translate("ICASPHPlusClass", "showBoundaryPar", Q_NULLPTR));
        actionIISPH->setText(QApplication::translate("ICASPHPlusClass", "IISPH", Q_NULLPTR));
        actionWCSPH->setText(QApplication::translate("ICASPHPlusClass", "WCSPH", Q_NULLPTR));
        SPHTypeName->setText(QApplication::translate("ICASPHPlusClass", "IISPH", Q_NULLPTR));
        label_width->setText(QApplication::translate("ICASPHPlusClass", "WIDTH:", Q_NULLPTR));
        label_length->setText(QApplication::translate("ICASPHPlusClass", "LENGTH:", Q_NULLPTR));
        label_height->setText(QApplication::translate("ICASPHPlusClass", "HEIGHT:", Q_NULLPTR));
        label_radius->setText(QApplication::translate("ICASPHPlusClass", "RADIUS:", Q_NULLPTR));
        menuSPH->setTitle(QApplication::translate("ICASPHPlusClass", "SPH Type", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class ICASPHPlusClass: public Ui_ICASPHPlusClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ICASPHPLUS_H
