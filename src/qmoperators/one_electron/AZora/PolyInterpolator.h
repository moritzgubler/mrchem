#include <vector>
#include <Eigen/Dense>

double polynomialInterpolate5(Eigen::Vector5d &x_in, Eigen::Vector5d &y_in, double x){
    double xm2 = x_in(0);
    double xm1 = x_in(1);
    double x00 = x_in(2);
    double xp1 = x_in(3);
    double xp2 = x_in(4);

    double ym2 = y_in(0);
    double ym1 = y_in(1);
    double y00 = y_in(2);
    double yp1 = y_in(3);
    double yp2 = y_in(4);
    return ym2 + (x - xm2)*((ym1 - ym2)/(xm1 - xm2) +
        (x - xm1)*(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2) +
        (x - x00)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
        ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/
        (-xm2 + xp1)) + (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/
        (-xm1 + xp1)) + ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/
        (-xm1 + xp2)))/(-xm2 + xp2))));
}

double polynomialInterpolate5_deriv(Eigen::Vector5d &x_in, Eigen::Vector5d &y_in, double x) {
    double xm2 = x_in(0);
    double xm1 = x_in(1);
    double x00 = x_in(2);
    double xp1 = x_in(3);
    double xp2 = x_in(4);

    double ym2 = y_in(0);
    double ym1 = y_in(1);
    double y00 = y_in(2);
    double yp1 = y_in(3);
    double yp2 = y_in(4);

    return (ym1 - ym2)/(xm1 - xm2) + (x - xm1)*(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2) +
             (x - x00)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                 ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
              ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1)) +
                   (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1)) +
                      ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/(-xm1 + xp2)))/(-xm2 + xp2))) +
        (x - xm2)*(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2) +
           (x - x00)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                 ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
              ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1)) +
                   (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1)) +
                      ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/(-xm1 + xp2)))/(-xm2 + xp2)) +
           (x - xm1)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                 ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
              ((x - x00)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1)) +
                   (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1)) +
                      ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/(-xm1 + xp2)))/(-xm2 + xp2) +
              ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
                        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1)) +
                   (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1)) +
                      ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/(-xm1 + xp2)))/(-xm2 + xp2)));
}

class PolyInterpolator {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    double :: xmin;
    double :: xmax;
    public:
    PolyInterpolator(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in){
        x = x_in;
        y = y_in;
        xmin = x_in(0);
        xmax = x_in(x_in.size() - 1);
    }
    
    void evalf(const double &x, double &y, double &yp) const {
        if (x < xmin || x > xmax) {
            if (x < xmin) {
                y = y(0);
                yp = 0;
            } else {
                y = y(y.size() - 1);
                yp = 0;
            }
        }

        // Find the interval in which x lies using binary search
        int i = 0;
        int j = x.size() - 1;
        while (j - i > 1) {
            int k = (i + j)/2;
            if (x < x(k)) {
                j = k;
            } else {
                i = k;
            }
        }

        Eigen::Vector5d x_in;
        Eigen::Vector5d y_in;

        if (i == 0) i = 2;
        if (i == 1) i = 2;
        if (i == x.size() - 1) i = x.size() - 3;
        if (i == x.size() - 2) i = x.size() - 3;
        for (int k = 0; k < 5; k++) {
            x_in(k) = x(i - 2 + k);
            y_in(k) = y(i - 2 + k);
        }

        y = polynomialInterpolate5(x_in, y_in, x);
        yp = polynomialInterpolate5_deriv(x_in, y_in, x);

    }
}