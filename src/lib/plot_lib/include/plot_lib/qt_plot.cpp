#include "qt_plot.h"
#include <QDebug>
#include <QCoreApplication>

#include <thread>


// ===========================================================

namespace as64_
{

namespace pl_
{

QtPlot *QtPlot::QtPlot_;
bool QtPlot::initialized = false;
int QtPlot::fig_count = 0;
std::map<Color, QColor> QtPlot::color;
std::map<pl_::PROPERTY, std::string> QtPlot::property_name;
const int QtPlot::N_COLORS = 13;

void QtPlot::init(QWidget *parent)
{
  if (!initialized)
  {
    // if (!QCoreApplication::instance()) throw std::runtime_error("[QtPlot::init]: There is no QApplication...");

    if (!QCoreApplication::instance())
    {
      Semaphore sem;

      std::thread([&sem]()
      {
        int argc = 0;
        char **argv = 0;
        QApplication app(argc, argv);
        QThread::currentThread()->setPriority(QThread::LowestPriority);
        app.setQuitOnLastWindowClosed(false);
//        QMainWindow *win = new QMainWindow();
//        win->setVisible(false);
        QtPlot::QtPlot_ = new QtPlot();
        QObject::connect(QtPlot_, SIGNAL(terminateSignal()), &app, SLOT(quit()));
        sem.notify();
        app.exec();

        std::cerr << "[QtPlot::init::thread]: Finished exec!\n";

        //delete (win);
        delete (QtPlot_);
      }).detach();

       sem.wait();
    }
    else
    {
      if (QCoreApplication::instance()->thread()->currentThreadId() != QThread::currentThread()->currentThreadId())
      {
        throw std::runtime_error("[QtPlot::init]: Must be called from the QApplication thread!");
      }
      QtPlot::QtPlot_ = new QtPlot(parent);
      // QtPlot::QtPlot_->moveToThread( QApplication::instance()->thread() );
    }

    initialized = true;
  }
  else
  {
    std::cerr << "\033[1m" << "\033[33m" << "[WARNING]: QtPlot has already been initialized!\033[0m\n";
  }
}

void QtPlot::terminate()
{
  QtPlot::QtPlot_->terminateSignal();
}

QColor QtPlot::getQColor(Color c)
{
  return (QtPlot::color.find(static_cast<Color>(c)))->second;
}

std::string QtPlot::getPropertyName(pl_::PROPERTY p)
{
  return QtPlot::property_name.find(p)->second;
}

Figure *QtPlot::figure(const std::string &title_, const std::vector<int> position)
{
  if (!QtPlot::initialized) throw std::runtime_error("[QtPlot::figure]: QtPlot has not been initialized...");

  Figure *fig;
  QtPlot::QtPlot_->figureSignal(&fig);

  if (!title_.empty()) fig->setTitle(title_);

  if (position.size() == 2) fig->resize(position[0], position[1]);
  else if (position.size() == 4) fig->setPosition(position[0], position[1], position[2], position[3]);
  else throw std::runtime_error("[QtPlot::figure]: \"position\" must have size 2 or 4...");

  return fig;
}

QtPlot::QtPlot(QWidget *parent)
{
  QtPlot::color[BLUE] = QColor(0, 0, 255);
  QtPlot::color[GREEN] = QColor(0, 255, 0);
  QtPlot::color[BROWN] = QColor(153, 51, 0);
  QtPlot::color[MAGENTA] = QColor(255, 0, 255);
  QtPlot::color[CYAN] = QColor(0, 255, 255);
  QtPlot::color[RED] = QColor(255, 0, 0);
  QtPlot::color[YELLOW] = QColor(230, 230, 0);
  QtPlot::color[LIGHT_BROWN] = QColor(217, 84, 26);
  QtPlot::color[PURPLE] = QColor(125, 46, 143);
  QtPlot::color[MUSTARD] = QColor(237, 176, 33);
  QtPlot::color[PINK] = QColor(255, 153, 199);
  QtPlot::color[BLACK] = QColor(0, 0, 0);
  QtPlot::color[GREY] = QColor(200, 200, 200);

  QtPlot::property_name[Color_] = "Color";
  QtPlot::property_name[LineWidth_] = "LineWidth";
  QtPlot::property_name[LineStyle_] = "LineStyle";
  QtPlot::property_name[MarkerSize_] = "MarkerSize";
  QtPlot::property_name[MarkerStyle_] = "MarkerStyle";
  QtPlot::property_name[FontSize_] = "FontSize";
  QtPlot::property_name[FontWeight_] = "FontWeight";
  QtPlot::property_name[FontFamily_] = "FontFamily";

  QObject::connect(this, SIGNAL(figureSignal(Figure **)), this, SLOT(figureSlot(Figure **)), Qt::BlockingQueuedConnection);
}

QtPlot::~QtPlot()
{
   // std::cerr << "[QtPlot::~QtPlot]: deleting self...\n";
}

void QtPlot::figureSlot(Figure **fig)
{
  *fig = new Figure(this);
  QtPlot::fig_count++;
}


// ===========================================
// =============   Figure   ==================
// ===========================================


Figure::Figure(QWidget *parent) : QMainWindow(parent)
{
  QObject::connect(this, &Figure::setAxesSignal, this, &Figure::setAxesSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Figure::clearAxesSignal, this, &Figure::clearAxesSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Figure::setTitleSignal, this, &Figure::setTitleSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Figure::resizeSignal, this, &Figure::resizeSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Figure::setPositionSignal, this, &Figure::setPositionSlot, Qt::BlockingQueuedConnection);

  // this->moveToThread( QApplication::instance()->thread() );

  central_widget = new QWidget(this);
  this->setCentralWidget(central_widget);

  QPalette pal = palette();
  pal.setColor(QPalette::Background, Qt::white);
  central_widget->setAutoFillBackground(true);
  central_widget->setPalette(pal);

  if (this->objectName().isEmpty()) this->setObjectName(QStringLiteral("Figure"));
  QMainWindow::resize(500, 400);
  grid_layout = new QGridLayout(central_widget);
  grid_layout->setObjectName(QStringLiteral("gridLayout"));

  this->setWindowTitle(QString("Figure ") + QString().setNum(QtPlot::fig_count+1));
  this->setAttribute(Qt::WA_QuitOnClose, true);

  setAxesSlot(1,1);

  // this->setWindowTitle(QApplication::translate("PlotDialog", "Dialog", 0));
  show();
}

Figure::~Figure()
{
  // std::cerr << "[Figure::~Figure]: deleting self...\n";
  clearAxesSlot();
  // delete grid_layout;
  QtPlot::fig_count--;
}

void Figure::closeEvent(QCloseEvent *event)
{
  delete this;
}

Axes *Figure::getAxes(int k)
{
  return axes[k];
}

Axes *Figure::getAxes(int row, int col)
{
  return axes[col + row*n2];
}

void Figure::setAxes(int n1, int n2)
{
  setAxesSignal(n1,n2);
}

void Figure::clearAxes(int k)
{
  clearAxesSignal(k);
}

void Figure::setTitle(const std::string &title_)
{
  setTitleSignal(QString(title_.c_str()));
}

void Figure::resize(int w, int h)
{
  resizeSignal(w,h);
}

void Figure::setPosition(int i1, int i2, int w, int h)
{
  setPositionSignal(i1, i2, w, h);
}


void Figure::setAxesSlot(int n1, int n2)
{
  this->n1 = n1;
  this->n2 = n2;

  clearAxesSlot();
  delete grid_layout;
  grid_layout = new QGridLayout(central_widget);
  // this->setLayout(grid_layout);

  axes.resize(n1*n2);
  for (int i=0; i<n1; i++)
  {
    for (int j=0; j<n2; j++)
    {
      int k = getAxesIndex(i,j);
      axes[k] = new Axes;
      grid_layout->addWidget(axes[k], i, j, 1, 1);
    }
  }
}

void Figure::clearAxesSlot(int k)
{
  if (k < 0)
  {
    for (int i=0;i<axes.size();i++) delete axes[i];
  }
  else
  {
    axes.erase(axes.begin()+k);
  }
}

void Figure::setTitleSlot(const QString &title_)
{
  this->setWindowTitle(title_);
}

void Figure::resizeSlot(int w, int h)
{
  QMainWindow::resize(w,h);
}

void Figure::setPositionSlot(int i1, int i2, int w, int h)
{
  setGeometry(i1, i2, w, h);
}

// =========================================
// =============   Axes   ==================
// =========================================

Axes::Axes(Figure *parent): QCustomPlot(parent)
{
  // qDebug() << "[Axes::Figure]: " << QThread::currentThread() << "\n";

  this->parent = parent;

  legend_ = new Legend(this->QCustomPlot::legend, this);
  title_ = 0;
  x_ax_label = new AxisLabel(this->xAxis);
  y_ax_label = new AxisLabel(this->yAxis);

  qRegisterMetaType<QVector<QString>>("QVector<QString>");

  QObject::connect(this, &Axes::holdSignal, this, &Axes::holdSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Axes::gridSignal, this, &Axes::gridSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Axes::setTitleSignal, this, &Axes::setTitleSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Axes::setXLabelSignal, this, &Axes::setXLabelSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Axes::setYLabelSignal, this, &Axes::setYLabelSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Axes::setLegendSignal, this, &Axes::setLegendSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, SIGNAL(plotSignal(const void *, const void *)), this, SLOT(plotSlot(const void *, const void *)), Qt::BlockingQueuedConnection);
  QObject::connect(this, &Axes::drawnowSignal, this, &Axes::drawnowSlot, Qt::BlockingQueuedConnection);

  this->resize(667, 452);

  this->setAttribute(Qt::WA_DeleteOnClose);

  // cplot = new QCustomPlot;
  // set locale to english, so we get english decimal separator:
  this->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom));
  // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
  this->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

  color_ind = 0;

  holdSlot(false);
  gridSlot(false);

  this->show();
}

Axes::~Axes()
{
  // std::cerr << "[Axes::~Axes]: deleting self...\n";
    // delete this;
}

void Axes::hold(bool set)
{
  holdSignal(set);
}

void Axes::grid(bool set)
{
  gridSignal(set);
}

Graph *Axes::plot(const arma::vec &data)
{
  arma::vec x_data(data.size());
  for (int i=0; i<data.size(); i++) x_data(i) = i;
  return plot(x_data, data);
}

Graph *Axes::plot(const arma::vec &x_data, const arma::vec &y_data)
{
  int n_data = x_data.size();

  if (y_data.size() != n_data) throw std::runtime_error("[Axes::plot]: Incompatible sizes between x_data and y_data");

  QVector<double> qx_data(n_data);
  QVector<double> qy_data(n_data);
  for (int i=0; i<n_data; i++)
  {
    qx_data[i] = x_data(i);
    qy_data[i] = y_data(i);
  }

  plotSignal(reinterpret_cast<const void *>(&qx_data), reinterpret_cast<const void *>(&qy_data));

  return last_graph;
}

Graph *Axes::plot(const arma::rowvec &data)
{
  arma::vec data2 = data.t();
  return plot(data2);
}

Graph *Axes::plot(const arma::rowvec &x_data, const arma::rowvec &y_data)
{
  arma::vec x_data2 = x_data.t();
  arma::vec y_data2 = y_data.t();
  return plot(x_data2, y_data2);
}

TextLabel *Axes::title(const std::string &title_text)
{
  setTitleSignal(QString(title_text.c_str()));
  return this->title_;
}

AxisLabel *Axes::xlabel(const std::string &label)
{
  setXLabelSignal(QString(label.c_str()));
  return x_ax_label;
}

AxisLabel *Axes::ylabel(const std::string &label)
{
  setYLabelSignal(QString(label.c_str()));
  return y_ax_label;
}

Legend *Axes::legend(const std::vector<std::string> &legend_labels)
{
  QVector<QString> qlegend_labels(legend_labels.size());
  for (int i=0; i<legend_labels.size(); i++) qlegend_labels[i] = legend_labels[i].c_str();
  setLegendSignal(qlegend_labels);

  return this->legend_;
}

void Axes::drawnow()
{
  drawnowSignal();
}

void Axes::plotSlot(const void *x_data, const void *y_data)
{
  // qDebug() << "[Axes::plotSlot]: " << QThread::currentThread() << "\n";
  const QVector<double> &qx_data = *reinterpret_cast<const QVector<double> *>(x_data);
  const QVector<double> &qy_data = *reinterpret_cast<const QVector<double> *>(y_data);

  if (!hold_on)
  {
    this->clearGraphs();
    color_ind = 0;
  }

  QCPGraph *graph = new QCPGraph(this->xAxis, this->yAxis);
  // graph->setLineStyle(QCPGraph::lsLine); // connect points with straight lines
  graph->setData(qx_data, qy_data);

  graphs.push_back( new Graph(graph, this) );
  last_graph = graphs.back();

  last_graph->setColorSlot(QtPlot::getQColor(static_cast<Color>(color_ind)));
  color_ind = (color_ind+1)%QtPlot::N_COLORS;

  last_graph->setLineStyleSlot(SolidLine);
  last_graph->setLineWidthSlot(2.0);
  last_graph->setMarkerStyleSlot(ssDot);
  last_graph->setMarkerSizeSlot(2.0);
  //graph->setPen(QPen(QBrush(c), 2.0, Qt::SolidLine)); // line color blue for first graph
  // graph->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDot, 4)); // set marker type and size. The width is determined by the linewidth
  this->rescaleAxes();
}

void Axes::holdSlot(bool set)
{
  hold_on = set;
}

void Axes::gridSlot(bool set)
{
  this->xAxis->grid()->setVisible(set);
  this->yAxis->grid()->setVisible(set);
}

void Axes::setTitleSlot(const QString &title)
{
  int fontsize = 16;
  QString font_family = "Arial";
  this->plotLayout()->insertRow(0);
  QCPTextElement *qcp_title = new QCPTextElement(this, title, QFont(font_family, fontsize, QFont::Normal));
  this->plotLayout()->addElement(0, 0, qcp_title);

  this->title_ = new TextLabel(qcp_title);
}

void Axes::setXLabelSlot(const QString &label)
{
  x_ax_label->setTextSlot(label);
  x_ax_label->setColorSlot(QColor(0,0,0));
  // x_axis->setLabelFont(QFont(font_family, fontsize, QFont::Normal));
}

void Axes::setYLabelSlot(const QString &label)
{
  y_ax_label->setTextSlot(label);
  y_ax_label->setColorSlot(QColor(0,0,0));
}

void Axes::setLegendSlot(const QVector<QString> &legend_labels)
{
  // qDebug() << "[Axes::setLegendSlot]: " << QThread::currentThread() << "\n";

  legend_->setLabelsSlot(legend_labels);
  legend_->setVisibleSlot(true);
  legend_->setBgColorSlot(QColor(255,255,255,230));
  legend_->setAlignmentSlot(Qt::AlignBottom|Qt::AlignRight);

  // int fontsize = 14;
  // QString font_family = "Arial";
  // qcp_legend->setFont(QFont(font_family, fontsize, QFont::Normal));
  // qcp_legend->setBrush(QBrush(QColor(255,255,255,230)));
  // this->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);
}

void Axes::drawnowSlot()
{
  this->replot();
}


// ==========================================
// =============   Graph   ==================
// ==========================================


Graph::Graph(QCPGraph *qcp_graph, Axes *parent)
{
  QObject::connect(this, &Graph::setColorSignal, this, &Graph::setColorSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Graph::setLineStyleSignal, this, &Graph::setLineStyleSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Graph::setLineWidthSignal, this, &Graph::setLineWidthSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Graph::setMarkerStyleSignal, this, &Graph::setMarkerStyleSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Graph::setMarkerSizeSignal, this, &Graph::setMarkerSizeSlot, Qt::BlockingQueuedConnection);

  this->qcp_graph = qcp_graph;
}

Graph::~Graph()
{

}

void Graph::setProperty() {}

void Graph::setPropertyHelper(pl_::PROPERTY p, pl_::Color p_value)
{
  if (p == pl_::Color_) setColor(p_value);
  else std::cerr << "[Graph::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Graph::setPropertyHelper(pl_::PROPERTY p, const QColor &p_value)
{
  if (p == pl_::Color_) setColor(p_value);
  else std::cerr << "[Graph::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Graph::setPropertyHelper(pl_::PROPERTY p, double p_value)
{
  if (p == pl_::LineWidth_) setLineWidth(p_value);
  else if (p == pl_::MarkerSize_) setMarkerSize(p_value);
  else std::cerr << "[Graph::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Graph::setPropertyHelper(pl_::PROPERTY p, pl_::LineStyle p_value)
{
  if (p == pl_::LineStyle_) setLineStyle(p_value);
  else std::cerr << "[Graph::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Graph::setPropertyHelper(pl_::PROPERTY p, pl_::MarkerStyle p_value)
{
  if (p == pl_::MarkerStyle_) setMarkerStyle(p_value);
  else std::cerr << "[Graph::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Graph::setColor(Color color)
{
  setColor(QtPlot::getQColor(static_cast<Color>(color)));
}

void Graph::setColor(const QColor &color)
{
  setColorSignal(color);
}

void Graph::setLineStyle(LineStyle style)
{
  setLineStyleSignal(style);
}

void Graph::setLineWidth(double width)
{
  setLineWidthSignal(width);
}

void Graph::setMarkerStyle(MarkerStyle type)
{
  setMarkerStyleSignal(type);
}

void Graph::setMarkerSize(double size)
{
  setMarkerSizeSignal(size);
}

void Graph::setColorSlot(const QColor &color)
{
  QPen pen = qcp_graph->pen();
  pen.setColor(color);
  qcp_graph->setPen(pen);
}

void Graph::setLineStyleSlot(int style)
{
  if (style == NoLine)
  {
    qcp_graph->setLineStyle(QCPGraph::lsNone);
  }
  else
  {
    qcp_graph->setLineStyle(QCPGraph::lsLine);
    QPen pen = qcp_graph->pen();
    pen.setStyle(static_cast<Qt::PenStyle>(style));
    qcp_graph->setPen(pen);
  }
}

void Graph::setLineWidthSlot(double width)
{
  QPen pen = qcp_graph->pen();
  pen.setWidth(width);
  qcp_graph->setPen(pen);
}

void Graph::setMarkerStyleSlot(int type)
{
  QCPScatterStyle marker = qcp_graph->scatterStyle();
  marker.setShape( static_cast<QCPScatterStyle::ScatterShape>(type));
  qcp_graph->setScatterStyle(marker);
}

void Graph::setMarkerSizeSlot(double size)
{
  QCPScatterStyle marker = qcp_graph->scatterStyle();
  marker.setSize(size);
  qcp_graph->setScatterStyle(marker);
}


// ===========================================
// =============   Legend   ==================
// ===========================================

Legend::Legend(QCPLegend *qcp_legend, Axes *axes)
{
  this->qcp_legend = qcp_legend;
  this->axes = axes;

  QObject::connect(this, &Legend::setAlignmentSignal, this, &Legend::setAlignmentSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Legend::setBgColorSignal, this, &Legend::setBgColorSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Legend::setVisibleSignal, this, &Legend::setVisibleSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Legend::setLabelsSignal, this, &Legend::setLabelsSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Legend::setColorSignal, this, &Legend::setColorSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Legend::setFontSizeSignal, this, &Legend::setFontSizeSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Legend::setFontFamilySignal, this, &Legend::setFontFamilySlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &Legend::setFontWeightSignal, this, &Legend::setFontWeightSlot, Qt::BlockingQueuedConnection);
}

void Legend::setLabels(const std::vector<std::string> &labels)
{
  QVector<QString> qlabels(labels.size());
  for (int i=0; i<labels.size(); i++) qlabels[i] = labels[i].c_str();
  setLabelsSignal(qlabels);
}

void Legend::setColor(Color c)
{
  setColor(QtPlot::getQColor(static_cast<Color>(c)));
}

void Legend::setColor(const QColor &c)
{
  setColorSignal(c);
}

void Legend::setFontSize(int size)
{
  setFontSizeSignal(size);
}

void Legend::setFontFamily(const std::string &family)
{
  setFontFamilySignal(QString(family.c_str()));
}

void Legend::setFontWeight(FontWeight fweight)
{
  setFontWeightSignal(fweight);
}

void Legend::setVisible(bool set)
{
  setVisibleSignal(set);
}

void Legend::setBgColor(const QColor &color)
{
  setBgColorSignal(color);
}

void Legend::setAlignment(Qt::Alignment alignment)
{
  setAlignmentSignal(alignment);
}


void Legend::setLabelsSlot(const QVector<QString> &labels)
{
  int n_graphs = axes->graphCount();
  int n_labels = labels.size();

  if (n_labels > n_graphs)
  {
    std::cerr << "[Legend::setLabelsSlot]: Extra legend labels will be ignored...\n";
    n_labels = n_graphs;
  }

  for (int i=0; i<n_labels; i++) axes->graph(i)->setName(labels[i]);

  // qcp_legend->setVisible(true);
  // qcp_legend->setBrush(QBrush(QColor(255,255,255,230)));
  // int fontsize = 14;
  // QString font_family = "Arial";
  // qcp_legend->setFont(QFont(font_family, fontsize, QFont::Normal));
  // axes->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);
}

void Legend::setColorSlot(const QColor &c)
{
  qcp_legend->setTextColor(c);
}

void Legend::setFontSizeSlot(int size)
{
  QFont font = qcp_legend->font();
  font.setPointSize(size);
  qcp_legend->setFont(font);
}

void Legend::setFontFamilySlot(const QString &family)
{
  QFont font = qcp_legend->font();
  font.setFamily(family);
  qcp_legend->setFont(font);
}

void Legend::setFontWeightSlot(FontWeight fweight)
{
  QFont font = qcp_legend->font();
  font.setWeight(fweight);
  qcp_legend->setFont(font);
}

void Legend::setVisibleSlot(bool set)
{
  this->qcp_legend->setVisible(set);
}

void Legend::setBgColorSlot(const QColor &color)
{
  this->qcp_legend->setBrush(color);
}

void Legend::setAlignmentSlot(Qt::Alignment alignment)
{
  axes->axisRect()->insetLayout()->setInsetAlignment(0, alignment);
}

void Legend::setProperty()
{}

void Legend::setPropertyHelper(pl_::PROPERTY p, pl_::Color p_value)
{
  if (p == pl_::Color_) setColor(p_value);
  else std::cerr << "[Legend::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Legend::setPropertyHelper(pl_::PROPERTY p, const QColor &p_value)
{
  if (p == pl_::Color_) setColor(p_value);
  else std::cerr << "[Legend::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Legend::setPropertyHelper(pl_::PROPERTY p, int p_value)
{
  if (p == pl_::FontSize_) setFontSize(p_value);
  else std::cerr << "[Legend::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Legend::setPropertyHelper(pl_::PROPERTY p, const std::string &p_value)
{
  if (p == pl_::FontFamily_) setFontFamily(p_value);
  else std::cerr << "[Legend::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void Legend::setPropertyHelper(pl_::PROPERTY p, FontWeight p_value)
{
  if (p == pl_::FontWeight_) setFontWeight(p_value);
  else std::cerr << "[Legend::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}


// ===============================================
// =============   Axis Label   ==================
// ===============================================

AxisLabel::AxisLabel(QCPAxis *qcp_axis)
{
  this->qcp_axis = qcp_axis;

  QObject::connect(this, &AxisLabel::setTextSignal, this, &AxisLabel::setTextSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &AxisLabel::setColorSignal, this, &AxisLabel::setColorSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &AxisLabel::setFontSizeSignal, this, &AxisLabel::setFontSizeSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &AxisLabel::setFontFamilySignal, this, &AxisLabel::setFontFamilySlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &AxisLabel::setFontWeightSignal, this, &AxisLabel::setFontWeightSlot, Qt::BlockingQueuedConnection);

}

void AxisLabel::setText(const std::string &s)
{
  setTextSignal(QString(s.c_str()));
}

void AxisLabel::setColor(Color c)
{
  setColor(QtPlot::getQColor(static_cast<Color>(c)));
}

void AxisLabel::setColor(const QColor &c)
{
  setColorSignal(c);
}

void AxisLabel::setFontSize(int size)
{
  setFontSizeSignal(size);
}

void AxisLabel::setFontFamily(const std::string &family)
{
  setFontFamilySignal(QString(family.c_str()));
}

void AxisLabel::setFontWeight(FontWeight fweight)
{
  setFontWeightSignal(fweight);
}

void AxisLabel::setProperty() {}

void AxisLabel::setPropertyHelper(pl_::PROPERTY p, pl_::Color p_value)
{
  if (p == pl_::Color_) setColor(p_value);
  else std::cerr << "[AxisLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void AxisLabel::setPropertyHelper(pl_::PROPERTY p, const QColor &p_value)
{
  if (p == pl_::Color_) setColor(p_value);
  else std::cerr << "[AxisLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void AxisLabel::setPropertyHelper(pl_::PROPERTY p, int p_value)
{
  if (p == pl_::FontSize_) setFontSize(p_value);
  else std::cerr << "[AxisLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void AxisLabel::setPropertyHelper(pl_::PROPERTY p, const std::string &p_value)
{
  if (p == pl_::FontFamily_) setFontFamily(p_value);
  else std::cerr << "[AxisLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void AxisLabel::setPropertyHelper(pl_::PROPERTY p, FontWeight p_value)
{
  if (p == pl_::FontWeight_) setFontWeight(p_value);
  else std::cerr << "[AxisLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void AxisLabel::setTextSlot(const QString &s)
{
  qcp_axis->setLabel(s);
}

void AxisLabel::setColorSlot(const QColor &c)
{
  qcp_axis->setLabelColor(c);
}

void AxisLabel::setFontSizeSlot(int size)
{
  QFont font = qcp_axis->labelFont();
  font.setPointSize(size);
  qcp_axis->setLabelFont(font);
}

void AxisLabel::setFontFamilySlot(const QString &family)
{
  QFont font = qcp_axis->labelFont();
  font.setFamily(family);
  qcp_axis->setLabelFont(font);
}

void AxisLabel::setFontWeightSlot(FontWeight fweight)
{
  QFont font = qcp_axis->labelFont();
  font.setWeight(fweight);
  qcp_axis->setLabelFont(font);
}



// ===============================================
// =============   Text Label   ==================
// ===============================================

TextLabel::TextLabel(QCPTextElement *qcp_text)
{
  this->qcp_text = qcp_text;

  QObject::connect(this, &TextLabel::setTextSignal, this, &TextLabel::setTextSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &TextLabel::setColorSignal, this, &TextLabel::setColorSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &TextLabel::setFontSizeSignal, this, &TextLabel::setFontSizeSlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &TextLabel::setFontFamilySignal, this, &TextLabel::setFontFamilySlot, Qt::BlockingQueuedConnection);
  QObject::connect(this, &TextLabel::setFontWeightSignal, this, &TextLabel::setFontWeightSlot, Qt::BlockingQueuedConnection);

}

void TextLabel::setText(const std::string &s)
{
  setTextSignal(QString(s.c_str()));
}

void TextLabel::setColor(Color c)
{
  setColor(QtPlot::getQColor(static_cast<Color>(c)));
}

void TextLabel::setColor(const QColor &c)
{
  setColorSignal(c);
}

void TextLabel::setFontSize(int size)
{
  setFontSizeSignal(size);
}

void TextLabel::setFontFamily(const std::string &family)
{
  setFontFamilySignal(QString(family.c_str()));
}

void TextLabel::setFontWeight(FontWeight fweight)
{
  setFontWeightSignal(fweight);
}

void TextLabel::setProperty() {}

void TextLabel::setPropertyHelper(pl_::PROPERTY p, pl_::Color p_value)
{
  if (p == pl_::Color_) setColor(p_value);
  else std::cerr << "[TextLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void TextLabel::setPropertyHelper(pl_::PROPERTY p, const QColor &p_value)
{
  if (p == pl_::Color_) setColor(p_value);
  else std::cerr << "[TextLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void TextLabel::setPropertyHelper(pl_::PROPERTY p, int p_value)
{
  if (p == pl_::FontSize_) setFontSize(p_value);
  else std::cerr << "[TextLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void TextLabel::setPropertyHelper(pl_::PROPERTY p, const std::string &p_value)
{
  if (p == pl_::FontFamily_) setFontFamily(p_value);
  else std::cerr << "[TextLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void TextLabel::setPropertyHelper(pl_::PROPERTY p, FontWeight p_value)
{
  if (p == pl_::FontWeight_) setFontWeight(p_value);
  else std::cerr << "[TextLabel::setPropertyHelper]: ** Invalid property \"" << QtPlot::getPropertyName(p) << "\" **\n";
}

void TextLabel::setTextSlot(const QString &s)
{
  qcp_text->setText(s);
}

void TextLabel::setColorSlot(const QColor &c)
{
  qcp_text->setTextColor(c);
}

void TextLabel::setFontSizeSlot(int size)
{
  QFont font = qcp_text->font();
  font.setPointSize(size);
  qcp_text->setFont(font);
}

void TextLabel::setFontFamilySlot(const QString &family)
{
  QFont font = qcp_text->font();
  font.setFamily(family);
  qcp_text->setFont(font);
}

void TextLabel::setFontWeightSlot(FontWeight fweight)
{
  QFont font = qcp_text->font();
  font.setWeight(fweight);
  qcp_text->setFont(font);
}



} // namespace pl_

} // namespace as64_
