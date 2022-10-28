#ifndef SPLINE_INTERPOLATION_CUBICSLPINE_H
#define SPLINE_INTERPOLATION_CUBICSLPINE_H

#include <vector>
#include <array>

std::vector<double> solveTriagonalSlae(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c, const std::vector<double> &d) {
    std::vector<double> p(d.size() - 1);
    std::vector<double> q(d.size() - 1);
    std::vector<double> x(d.size());

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];
    for (int i = 1; i < p.size(); ++i) {
        p[i] = -c[i] / (a[i] * p[i - 1] + b[i]);
        q[i] = (d[i] - a[i] * q[i - 1]) / (a[i] * p[i - 1] + b[i]);
    }
    x.back() = (d.back() - a.back() * q.back()) / (a.back() * p.back() + b.back());
    for (int i = x.size() - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }
    return x;
}

class CubicSpline {
private:
    std::vector<std::array<double, 4>> a;
    std::vector<double> x0;

    void initialize(const std::vector<double> &x, const std::vector<double> &f) {
        for (int i = 0; i < a.size(); ++i) {
            a[i][0] = f[i + 1];
        }
        a.back()[2] = 0;

        std::vector<double> alpha(a.size() - 1);
        std::vector<double> beta(a.size() - 1, 2);
        std::vector<double> gamma(a.size() - 1);
        std::vector<double> delta(a.size() - 1);

        for (int i = 0; i < delta.size(); ++i) {
            if (i != 0)
                alpha[i] = (x[i + 1] - x[i]) / (x[i + 2] - x[i]);
            if (i != delta.size() - 1) {
                gamma[i] = (x[i + 2] - x[i + 1]) / (x[i + 2] - x[i]);
            }
            delta[i] = ((f[i + 2] - f[i + 1]) / (x[i + 2] - x[i + 1]) - (f[i + 1] - f[i]) / (x[i + 1] - x[i])) /
                       (x[i + 2] - x[i]);
        }
        auto c = solveTriagonalSlae(alpha, beta, gamma, delta);

        for (int i = 0; i < c.size(); ++i) {
            a[i][2] = c[i];
        }

        a[0][1] = a[0][2] * (x[1] - x[0]) / 3 + (f[1] - f[0]) / (x[1] - x[0]);
        a[0][3] = a[0][2] / (x[1] - x[0]);

        for (int i = 1; i < a.size(); ++i) {
            a[i][1] = (a[i][2] / 3 + a[i - 1][2] / 6) * (x[i + 1] - x[i]) + (f[i + 1] - f[i]) / (x[i + 1] - x[i]);
            a[i][3] = (a[i][2] - a[i - 1][2]) / (x[i + 1] - x[i]);
        }
    }

public:
    CubicSpline(const std::vector<double> &x, const std::vector<double> &f) : x0(x), a(std::vector<std::array<double, 4>>(x.size() - 1)) {
        initialize(x, f);
    }

    double interpolate(double x) {
        for (int i = 0; i < a.size(); ++i) {
            if (x0[i] <= x && x <= x0[i + 1]) {
                double d = x - x0[i + 1];
                return (a[i][0] + a[i][1] * d + a[i][2] * d * d / 2 + a[i][3] * d * d * d / 6);
            }
        }
        return 0.0;
    }
};


#endif //SPLINE_INTERPOLATION_CUBICSLPINE_H
