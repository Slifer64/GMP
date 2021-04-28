#ifndef AS64_QTPLOT_H
#define AS64_QTPLOT_H

#include <QDialog>
#include <QLabel>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QCloseEvent>
#include <QVector>
#include <QApplication>
#include <QMainWindow>

#include <vector>
#include <cstring>
#include <exception>
#include <memory>
#include <armadillo>

#include "qcustomplot.h"

#include <plot_lib/utils.h>

// Q_DECLARE_METATYPE(QVector<QString>);

namespace as64_
{

namespace pl_
{

class Figure; // forward declaration
class Axes; // forward declaration

enum PROPERTY
{
  Color_,
  LineWidth_,
  LineStyle_,
  MarkerSize_,
  MarkerStyle_,
  FontSize_,
  FontWeight_,
  FontFamily_,
};

enum Color {
  BLUE = 0,        //   0   0   255
  GREEN = 1,       //   0  255    0
  MAGENTA = 2,     // 255    0  255
  BROWN = 3,       // 153   51    0
  CYAN = 4,        //   0  255  255
  RED = 5,         // 255    0    0
  YELLOW = 6,      // 230  230    0
  LIGHT_BROWN = 7, // 217   84   26
  PURPLE = 8,      // 125   46  143
  MUSTARD = 9,     // 237  176   33
  PINK = 10,        // 255  153  199
  BLACK = 11,       //   0    0    0
  GREY = 12,        // 200  200  200
};

enum LineStyle {
  NoLine = 0,
  SolidLine = 1,
  DashLine = 2,
  DotLine = 3,
  DashDotLine = 4,
  DashDotDotLine = 5
};

enum MarkerStyle {
  ssNone,
  ssDot,
  ssCross,
  ssPlus,
  ssCircle,
  ssDisc,
  ssSquare,
  ssDiamond,
  ssStar,
  ssTriangle,
  ssTriangleInverted,
  ssCrossSquare,
  ssPlusSquare,
  ssCrossCircle,
  ssPlusCircle,
  ssPeace
};

enum FontWeight {
  Light = 25,
  Normal = 50,
  DemiBold = 63,
  Bold = 75,
  Black = 87
};


// ===============================================
// =============   Text Label   ==================
// ===============================================

class TextLabel : public QWidget
{

Q_OBJECT

public:
  TextLabel(QCPTextElement *qcp_text);

  template<typename T, typename... Arguments>
  void setProperty(PROPERTY p, T p_value, Arguments... parameters)
  {
    setPropertyHelper(p, p_value);
    setProperty(parameters...);
  }

  void setText(const std::string &s);
  void setColor(Color c);
  void setColor(const QColor &c);
  void setFontSize(int size);
  void setFontFamily(const std::string &family);
  void setFontWeight(FontWeight fweight);

signals:
  void setTextSignal(const QString &s);
  void setColorSignal(const QColor &c);
  void setFontSizeSignal(int size);
  void setFontFamilySignal(const QString &family);
  void setFontWeightSignal(FontWeight fweight);

private slots:
  void setTextSlot(const QString &s);
  void setColorSlot(const QColor &c);
  void setFontSizeSlot(int size);
  void setFontFamilySlot(const QString &family);
  void setFontWeightSlot(FontWeight fweight);

private:
  void setProperty();
  void setPropertyHelper(PROPERTY p, Color p_value);
  void setPropertyHelper(pl_::PROPERTY p, const QColor &p_value);
  void setPropertyHelper(PROPERTY p, int p_value);
  void setPropertyHelper(PROPERTY p, const std::string &p_value);
  void setPropertyHelper(PROPERTY p, FontWeight p_value);

  QCPTextElement *qcp_text;
};


// ===============================================
// =============   Axis Label   ==================
// ===============================================

class AxisLabel : public QWidget
{

Q_OBJECT

friend Axes;

public:
  AxisLabel(QCPAxis *qcp_axis);

  template<typename T, typename... Arguments>
  void setProperty(PROPERTY p, T p_value, Arguments... parameters)
  {
    setPropertyHelper(p, p_value);
    setProperty(parameters...);
  }

  void setText(const std::string &s);
  void setColor(Color c);
  void setColor(const QColor &c);
  void setFontSize(int size);
  void setFontFamily(const std::string &family);
  void setFontWeight(FontWeight fweight);

signals:
  void setTextSignal(const QString &s);
  void setColorSignal(const QColor &c);
  void setFontSizeSignal(int size);
  void setFontFamilySignal(const QString &family);
  void setFontWeightSignal(FontWeight fweight);

private slots:
  void setTextSlot(const QString &s);
  void setColorSlot(const QColor &c);
  void setFontSizeSlot(int size);
  void setFontFamilySlot(const QString &family);
  void setFontWeightSlot(FontWeight fweight);

private:
  void setProperty();
  void setPropertyHelper(PROPERTY p, Color p_value);
  void setPropertyHelper(pl_::PROPERTY p, const QColor &p_value);
  void setPropertyHelper(PROPERTY p, int p_value);
  void setPropertyHelper(PROPERTY p, const std::string &p_value);
  void setPropertyHelper(PROPERTY p, FontWeight p_value);

  QCPAxis *qcp_axis;
};


// ==========================================
// =============   Legend  ==================
// ==========================================

class Legend : public QWidget
{

Q_OBJECT

  friend Axes;

public:
  Legend(QCPLegend *qcp_legend, Axes *axes);

  template<typename T, typename... Arguments>
  void setProperty(PROPERTY p, T p_value, Arguments... parameters)
  {
    setPropertyHelper(p, p_value);
    setProperty(parameters...);
  }

  void setLabels(const std::vector<std::string> &labels);

  void setColor(Color c);
  void setColor(const QColor &c);
  void setFontSize(int size);
  void setFontFamily(const std::string &family);
  void setFontWeight(FontWeight fweight);

  void setVisible(bool set);
  void setBgColor(const QColor &color);
  void setAlignment(Qt::Alignment alignment);

//  void setColor(Color c, int i);
//  void setColor(const QColor &c, int i);
//  void setFontSize(int size, int i);
//  void setFontFamily(const std::string &family, int i);
//  void setFontWeight(FontWeight fweight, int i);

signals:
  void setLabelsSignal(const QVector<QString> &labels);
  void setColorSignal(const QColor &c);
  void setFontSizeSignal(int size);
  void setFontFamilySignal(const QString &family);
  void setFontWeightSignal(FontWeight fweight);

  void setVisibleSignal(bool set);
  void setBgColorSignal(const QColor &color);
  void setAlignmentSignal(Qt::Alignment alignment);

private slots:
  void setLabelsSlot(const QVector<QString> &labels);
  void setColorSlot(const QColor &c);
  void setFontSizeSlot(int size);
  void setFontFamilySlot(const QString &family);
  void setFontWeightSlot(FontWeight fweight);

  void setVisibleSlot(bool set);
  void setBgColorSlot(const QColor &color);
  void setAlignmentSlot(Qt::Alignment alignment);

private:
  void setProperty();
  void setPropertyHelper(PROPERTY p, Color p_value);
  void setPropertyHelper(pl_::PROPERTY p, const QColor &p_value);
  void setPropertyHelper(PROPERTY p, int p_value);
  void setPropertyHelper(PROPERTY p, const std::string &p_value);
  void setPropertyHelper(PROPERTY p, FontWeight p_value);

  QCPLegend *qcp_legend;
  Axes *axes;
};

// ==========================================
// =============   Graph   ==================
// ==========================================

class Graph : public QWidget {

Q_OBJECT

friend Axes;

public:
  Graph(QCPGraph *qcp_graph, Axes *parent = 0);

  ~Graph();

  template<typename T, typename... Arguments>
  void setProperty(pl_::PROPERTY p, T p_value, Arguments... parameters)
  {
    setPropertyHelper(p, p_value);
    setProperty(parameters...);
  }

  void setColor(Color color);

  void setColor(const QColor &color);

  void setLineStyle(LineStyle style);

  void setLineWidth(double width);

  void setMarkerStyle(MarkerStyle type);

  void setMarkerSize(double size);

signals:

  void setColorSignal(const QColor &color);

  void setLineStyleSignal(int style);

  void setLineWidthSignal(double width);

  void setMarkerStyleSignal(int type);

  void setMarkerSizeSignal(double size);

private slots:

  void setColorSlot(const QColor &color);

  void setLineStyleSlot(int style);

  void setLineWidthSlot(double width);

  void setMarkerStyleSlot(int type);

  void setMarkerSizeSlot(double size);

private:
  void setProperty();
  void setPropertyHelper(pl_::PROPERTY p, pl_::Color p_value);
  void setPropertyHelper(pl_::PROPERTY p, const QColor &p_value);
  void setPropertyHelper(pl_::PROPERTY p, double p_value);
  void setPropertyHelper(pl_::PROPERTY p, pl_::LineStyle p_value);
  void setPropertyHelper(pl_::PROPERTY p, pl_::MarkerStyle p_value);

  Axes *parent;
  QCPGraph *qcp_graph;
};


// =========================================
// =============   Axes   ==================
// =========================================

class Axes : public QCustomPlot {
Q_OBJECT

friend Legend;

public:

  Axes(Figure *parent=0);

  ~Axes();

  void hold(bool set);

  void grid(bool set);

  Graph *plot(const arma::vec &data);
  Graph *plot(const arma::rowvec &data);

  Graph *plot(const arma::vec &x_data, const arma::vec &y_data);
  Graph *plot(const arma::rowvec &x_data, const arma::rowvec &y_data);

  template<typename... Arguments>
  Graph *plot(const arma::rowvec  &x, const arma::rowvec  &y, Arguments... properties)
  {
    Graph *graph = plot(x, y);
    graph->setProperty(properties...);
    return graph;
  }

  template<typename... Arguments>
  Graph *plot(const arma::vec  &x, const arma::vec  &y, Arguments... properties)
  {
    Graph *graph = plot(x, y);
    graph->setProperty(properties...);
    return graph;
  }

  template<typename... Arguments>
  TextLabel *title(const std::string &title_text, Arguments... properties)
  {
    TextLabel *title_label = this->title(title_text);
    title_label->setProperty(properties...);
    return title_label;
  }

  TextLabel *title(const std::string &title_text);

  AxisLabel *xlabel(const std::string &label);

  template<typename... Arguments>
  AxisLabel *xlabel(const std::string &label, Arguments... properties)
  {
    AxisLabel *axis_label = this->xlabel(label);
    axis_label->setProperty(properties...);
    return axis_label;
  }

  AxisLabel *ylabel(const std::string &label);

  template<typename... Arguments>
  AxisLabel *ylabel(const std::string &label, Arguments... properties)
  {
    AxisLabel *axis_label = this->ylabel(label);
    axis_label->setProperty(properties...);
    return axis_label;
  }

  template<typename... Arguments>
  Legend *legend(const std::vector<std::string> &legend_labels, Arguments... properties)
  {
    Legend *legend_elem = this->legend(legend_labels);
    legend_elem->setProperty(properties...);
    return legend_elem;
  }

  Legend *legend(const std::vector<std::string> &legend_labels);

  void drawnow();

signals:

  void holdSignal(bool set);

  void gridSignal(bool set);

  void plotSignal(const void *x_data, const void *y_data);

  void setTitleSignal(const QString &title);

  void setXLabelSignal(const QString &label);

  void setYLabelSignal(const QString &label);

  void setLegendSignal(const QVector<QString> &legend_labels);

  void drawnowSignal();

private slots:

  void holdSlot(bool set);

  void gridSlot(bool set);

  void plotSlot(const void *x_data, const void *y_data);

  void setTitleSlot(const QString &title);

  void setXLabelSlot(const QString &label);

  void setYLabelSlot(const QString &label);

  void setLegendSlot(const QVector<QString> &legend_labels);

  void drawnowSlot();


private:
  bool hold_on;
  Figure *parent;
  std::vector<Graph *> graphs;

  QCPAxisRect *axes;
  QCPLayoutGrid *axes_grid;

  int color_ind;

  Graph *last_graph;
  TextLabel *title_;
  Legend *legend_;
  AxisLabel *x_ax_label;
  AxisLabel *y_ax_label;
};

// ===========================================
// =============   Figure   ==================
// ===========================================

class Figure : public QMainWindow {
Q_OBJECT

public:

  Figure(QWidget *parent = 0);

  ~Figure();

  Axes *getAxes(int k);
  Axes *getAxes(int row, int col);
  void setAxes(int n1, int n2);
  void clearAxes(int k = -1);
  void setTitle(const std::string &title_);
  void resize(int w, int h);
  void setPosition(int i1, int i2, int w, int h);

signals:

  void setAxesSignal(int n1, int n2);
  void clearAxesSignal(int k = -1);
  void setTitleSignal(const QString &title_);
  void resizeSignal(int w, int h);
  void setPositionSignal(int i1, int i2, int w, int h);

private slots:

  void setAxesSlot(int n1, int n2);
  void clearAxesSlot(int k = -1);
  void setTitleSlot(const QString &title_);
  void resizeSlot(int w, int h);
  void setPositionSlot(int i1, int i2, int w, int h);

  void closeEvent(QCloseEvent *event) override;

private:
  QGridLayout *grid_layout;
  std::vector<Axes *> axes;
  QWidget *central_widget;
  int n1;
  int n2;

  int getAxesIndex(int i, int j) { return j + i * n2; }
};


// ===========================================
// ===========================================

class QtPlot : public QWidget {
Q_OBJECT

public:
  QtPlot(QWidget *parent = 0);

  ~QtPlot();

  static void init(QWidget *parent = 0);
  static void terminate();

  static Figure *figure(const std::string &title_="", const std::vector<int> position={500, 400});

  static int fig_count;

  static QColor getQColor(Color c);
  static std::string getPropertyName(pl_::PROPERTY p);

  static const int N_COLORS;

signals:

  void figureSignal(Figure **);
  void terminateSignal();

private slots:

  void figureSlot(Figure **);

private:
  static bool initialized;
  static QtPlot *QtPlot_;
  static std::map<Color, QColor> color;
  static std::map<pl_::PROPERTY, std::string> property_name;
};

} // namespace pl_

} // namespace as64_

#endif // AS64_QTPLOT_H
