#include <iostream>
#include <random>

using namespace std;

random_device rd;
mt19937 gen(rd());

double getRand(double s = 0, double e = 1) {
    uniform_real_distribution<> dist(s, e);
    return dist(gen);
}

double getNormalRand(double m, double sd) {
    normal_distribution<> dist(m, sd);
    return dist(gen);
}

int getIntRand(int s, int e) {
    uniform_int_distribution<> dist(s, e-1);
    return dist(gen);
}

class Point {
public:
    Point(double xc = 0, double yc = 0) {
        x = xc;
        y = yc;
        color = -1;
    }

    static double distance(Point a, Point b) {
        return sqrt((a.getX()-b.getX())*(a.getX()-b.getX())+(a.getY()-b.getY())*(a.getY()-b.getY()));
    }

    int getColor() {
        return color;
    }

    void setColor(int c) {
        color = c;
    }

    double getX() {
        return x;
    }

    double getY() {
        return y;
    }

    Point operator+(Point p) {
        return Point(x+p.getX(), y+p.getY());
    }

    Point operator-(Point p) {
        return Point(x-p.getX(), y-p.getY());
    }

    Point operator*(double n) {
        return Point(x*n, y*n);
    }

    Point operator/(double n) {
        return Point(x/n, y/n);
    }

    void operator+=(Point p) {
        x += p.getX();
        y += p.getY();
    }

    void operator-=(Point p) {
        x -= p.getX();
        y -= p.getY();
    }

    void operator*=(double n) {
        x *= n;
        y *= n;
    }

    void operator/=(double n) {
        x /= n;
        y /= n;
    }
private:
    double x;
    double y;
    int color;
};

void cluster(Point *p, int n, Point *m, int nm) {
    for (int i = 0; i < n; ++i) {
        int min = -1;
        double mind = -1;
        for (int j = 0; j < nm; ++j) {
            double d = Point::distance(p[i], m[j]);
            if (d < mind || min == -1) {
                mind = d;
                min = j;
            }
        }
        p[i].setColor(min);
    }
    for (int i = 0; i < nm; ++i) {
        Point c(0, 0);
        int num = 0;
        for (int j = 0; j < n; ++j) {
            if (p[j].getColor() == i) {
                c += p[j];
                num++;
            }
        }
        m[i] = c/num;
    }
}

void getCluster(Point *p, int n, int c) {
    Point *m = new Point[c];
    for (int i = 0; i < c; ++i) {
        m[i] = p[getIntRand(i*n/c, (i+1)*n/c)];
    }
    Point *pm = new Point[c];
    for (int i = 0; i < c; ++i) {
        pm[i] = Point(m[i].getX(), m[i].getY());
    }
    while (true) {
        cluster(p, n, m, c);
        bool next = false;
        for (int i = 0; i < c; ++i) {
            if (Point::distance(m[i], pm[i]) > 0.0000000001) {
                next = true;
                break;
            }
        }
        if (!next) break;
        for (int i = 0; i < c; ++i) {
            pm[i] = Point(m[i].getX(), m[i].getY());
        }
    }
}

int main(int argc, char const * argv[] ) {
    int n = 100, c = 5;
    if (argc >= 2) {
        n = stoi(argv[1]);
    }
    if (argc == 3) {
        c = stoi(argv[2]);
    }
    Point *p = new Point[n];
    Point *m = new Point[c]; // To generate clustered random points
    for (int i = 0; i < c; ++i) {
        m[i] = Point(getRand(), getRand());
    }
    for (int i = c; i < n; ++i) {
        Point pp = m[getIntRand(0, c)];
        double sd = getRand(0.4, 0.8)/c;
        p[i] = Point(getNormalRand(pp.getX(), sd), getNormalRand(pp.getY(), sd));
    }
    delete m;
    getCluster(p, n, c);
    for (int j = 0; j < n; ++j) {
        cout << p[j].getX() << " " << p[j].getY() << " " << p[j].getColor() << endl;
    }
    delete p;
    return 0;
}
