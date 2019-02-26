#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtCore/QFile>
#include <QtCore/QTextStream>

QT_CHARTS_USE_NAMESPACE

int plotter(int argc, char *argv[])
{
    QApplication a(argc, argv);

    //![1]
    QLineSeries *series = new QLineSeries();
    //![1]

    //![2]
    QFile sunSpots("../data/energy.dat");
    if (!sunSpots.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return 1;
    }

    QTextStream stream(&sunSpots);
    int i = 1;
    while (!stream.atEnd()) {
        QString line = stream.readLine();
        QStringList values = line.split(" ", QString::SkipEmptyParts);
        series->append(i, values[0].toDouble());
        i++;
    }
    sunSpots.close();
    //![2]

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
    axisY->setLabelFormat("%i");
    axisY->setTitleText("Average energy");
    axisY->setRange(2, 17);
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
    window.setWhatsThis("HHe");
    window.resize(820, 600);
    window.show();
    //![6]

    return a.exec();
}
