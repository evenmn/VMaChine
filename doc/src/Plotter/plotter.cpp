#include "plotter.h"
#include "../system.h"
#include "../sampler.h"

#include <string>
#include <iostream>

Plotter::Plotter(System* system) {
    m_system = system;
}

QT_CHARTS_USE_NAMESPACE

int Plotter::plotEnergy(int argc, char *argv[]) {
    QApplication a(argc, argv);

    //![1]
    QLineSeries *series = new QLineSeries();
    //![1]

    //![2]
    std::string filename = m_system->getSampler()->generateFileName("energy", ".dat");
    QFile energies(QString::fromStdString(filename));
    if (!energies.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return 1;
    }

    QTextStream stream(&energies);
    double maxElement = -1e5;
    double minElement = +1e5;
    int i = 1;
    while (!stream.atEnd()) {
        QString line = stream.readLine();
        QStringList values = line.split(" ", QString::SkipEmptyParts);
        if(values[0].toDouble() > maxElement) { maxElement = values[0].toDouble(); }
        if(values[0].toDouble() < minElement) { minElement = values[0].toDouble(); }
        series->append(i, values[0].toDouble());
        i++;
    }
    energies.close();
    //![2]

    double diff = maxElement - minElement;

    //![3]
    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->legend()->hide();
    chart->setTitle("Average energy");
    //![3]

    //![4]
    QValueAxis *axisX = new QValueAxis;
    axisX->setLabelFormat("%i");
    axisX->setTitleText("Iteration");
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis;
    axisY->setTitleText("Average energy");
    axisY->setRange(minElement - diff/10, maxElement + diff/10);
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);
    //![4]

    //![5]
    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    //![5]

    //![6]
    QMainWindow window;
    window.setCentralWidget(chartView);
    //window.setWhatsThis("HHe");
    window.resize(820, 600);
    window.show();
    //![6]

    return a.exec();
}

int Plotter::plotOneBodyDensity(int argc, char *argv[]) {
    QApplication a(argc, argv);

    //![1]
    QLineSeries *series = new QLineSeries();
    //![1]

    //![2]
    std::string filename = m_system->getSampler()->generateFileName("onebody", ".dat");
    QFile densities(QString::fromStdString(filename));
    if (!densities.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return 1;
    }

    QTextStream stream(&densities);
    double maxElement = -1e5;
    double minElement = +1e5;
    int i = 1;
    while (!stream.atEnd()) {
        QString line = stream.readLine();
        QStringList values = line.split(" ", QString::SkipEmptyParts);
        if(values[0].toDouble() > maxElement) { maxElement = values[0].toDouble(); }
        if(values[0].toDouble() < minElement) { minElement = values[0].toDouble(); }
        series->append(i, values[0].toDouble());
        i++;
    }
    densities.close();
    //![2]

    double diff = maxElement - minElement;

    //![3]
    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->legend()->hide();
    chart->setTitle("One-body Density");
    //![3]

    //![4]
    QValueAxis *axisX = new QValueAxis;
    axisX->setLabelFormat("%i");
    axisX->setTitleText("Radial distance");
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis;
    axisY->setTitleText("Density");
    axisY->setRange(minElement - diff/10, maxElement + diff/10);
    chart->addAxis(axisY, Qt::AlignLeft);
    series->attachAxis(axisY);
    //![4]

    //![5]
    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);
    //![5]

    //![6]
    QMainWindow window;
    window.setCentralWidget(chartView);
    //window.setWhatsThis("HHe");
    window.resize(820, 600);
    window.show();
    //![6]

    return a.exec();
}

Plotter::~Plotter() {};
