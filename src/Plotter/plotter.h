#pragma once

#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtCore/QFile>
#include <QtCore/QTextStream>

class Plotter {
public:
    Plotter(class System* system);
    int plotEnergy(int argc, char *argv[]);
    int plotOneBodyDensity(int argc, char *argv[]);
    ~Plotter();
protected:
    class System* m_system = nullptr;
};
