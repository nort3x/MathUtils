#include "MathUtils.h"
#include "gnuplot.h"
#include "DataType.h"

typedef Algebric::MultiDimPoint<double> Vec3;
const double G = 1;

class Body {
private:
    Vec3 position;
    Vec3 velocity;
    double mass;
    std::function<Vec3(Vec3, Vec3, double)> Field_v_r_t_get = [](const Vec3 &v, const Vec3 &r, double t) -> Vec3 {
        return {0, 0, 0};
    };
    std::function<Vec3(Vec3, Vec3, double)> Field_v_r_t_send = [this](const Vec3 &v, const Vec3 &r, double t) -> Vec3 {
        return -(G * mass / std::pow((r + (-1 * position)).Norm(), 3)) * (r + (-1 * position));

    };
public:
    Body(double _mass, const Vec3 &pos, const Vec3 &vel) {
        mass = _mass;
        position = pos;
        velocity = vel;
    }

    Vec3 &acc_Position() {
        return position;
    }

    Vec3 &acc_getVelocity() {
        return velocity;
    }

    double &acc_Mass() {
        return mass;
    }

    void AcceptForce(const std::function<Vec3(Vec3, Vec3, double)> &Force) {
        Field_v_r_t_get = [&](const Vec3 &v, const Vec3 &r, double t) -> Vec3 {
            Vec3 v1 = Field_v_r_t_get(v, r, t);
            Vec3 v2 = mass * Force(v, r, t);
            return v1 + v2;
        };
    }

    std::function<Vec3(Vec3, Vec3, double)> &getForce_byref() {
        return Field_v_r_t_send;
    }

    Vec3 UnderWhatForce_question_mark() {
        return Field_v_r_t_get(velocity, position, 0);
    }
};

class World {
private:
    std::vector<Body> Bodies;
    double tp;
public:
    World() {}

    World(const std::vector<Body> &b) {
        Bodies = b;
    }

    void AddBody(Body &b) {
        for (auto bs: Bodies) {
            bs.AcceptForce(b.getForce_byref());
        }
        Bodies.push_back(b);
    }

};


int main() {
    Body b1(1, {1, 0, 0}, {0, 0, 0});
    Body b2(1, {0, 1, 0}, {0, 0, 0});
    World w;
    w.AddBody(b1);
    w.AddBody(b2);

    std::cout << b1.UnderWhatForce_question_mark();

    std::function<double(double)> potential = [](double x) {
        return x * x;
    };

    plt::TwoDim::Gnuplot2D<double, double> g("/root/Desktop/", "set yr [-3:3]");
    for (int i = 0; i <= 10; ++i) {
        double E = 0.1 * i;
        std::function<double(double, double, double)> shro = [&potential, E](double si_p, double si, double x) {
            double h = 1;
            double m = 1;
            return (2 * m / (h * h)) * (potential(x) - E) * si;
        };
        g.addFunctionPlot(Calculus::SingleVar::ODE<double, double>::EulerMethodSecondOrder(shro, 0, 0, 0.5, 0.1),
                          -2 * M_PI, 2 * M_PI, 100, "w l notitle");
    }
    double E = 2*(3+0.5);
    std::function<double(double, double, double)> shro = [&potential, E](double si_p, double si, double x) {
        double h = 1;
        double m = 1;
        return (2 * m / (h * h)) * (potential(x) - E) * si;};
    g.addFunctionPlot(Calculus::SingleVar::ODE<double, double>::EulerMethodSecondOrder(shro, 0, 0, 0.5, 0.001),
                      -2 * M_PI, 2 * M_PI, 100, "w l title 'E=0.7025' lw 3");
    g.plot();

    return 0;
}