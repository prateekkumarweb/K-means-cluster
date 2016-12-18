#include <iostream>
#include <fstream>
#include <random>
#include <ctime>

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

int getIntRand(int n, double *d) {
    double r = getRand();
    for (int i=0; i<n; i++) {
        if (i!=0) {
            d[i] += d[i-1];
        }
        if (r < d[i]) {
            return i;
        }
    }
    return n-1;
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

void genRandPoints(Point *p, int n, Point *m, int c) {
    double *d = new double[n];
    double sum = 0;
    int *cc = new int[c];
    cc[0] = getIntRand(0, n);
    for (int i=1; i<c; i++) {
        sum = 0;
        for (int j=0; j<n; j++) {
            d[j] = -1;
            for (int k=0; k<=i; k++) {
                if (cc[k] != j) {
                    double dst = Point::distance(p[cc[k]], p[j]);
                    dst = dst*dst;
                    if (d[j] < 0 || dst < d[j]) {
                        d[j] = dst;
                    }
                }
                else {
                    d[j] = 0;
                }
            }
            sum += d[j];
        }
        for (int j=0; j<n; j++) {
            d[j] = d[j] / sum;
        }
        cc[i] = getIntRand(n, d);
    }
    for (int i=0; i<c; i++) {
        m[i] = p[cc[i]];
    }
}

void printDBIndex(Point *p, int n, Point *m, int c) {
    double *s = new double[c];
    for (int i=0; i<c; i++) {
        s[i] = 0;
        int k = 0;
        for (int j=0; j<n; j++) {
            if (i == p[j].getColor()) {
                s[i] += Point::distance(m[i], p[j]);
                k++;
            }
        }
        s[i] /= k;
    }
    double d;
    for (int i=0; i<c; i++) {
        double max = -1;
        for (int j=0; j<c; j++) {
            if (j != i) {
                double v = (s[i]+s[j])/Point::distance(m[i], m[j]);
                if (max < v) {
                    max = v;
                }
            }
        }
        d += max;
    }
    d = d/c;
    cerr << "Davies-Bouldin Index : " << d << endl;
}

void printDIndex(Point *p, int n, Point *m, int c) {
    double min = -1;
    for (int i=0; i<c; i++) {
        for (int j=i+1; j<c; j++) {
            double d = Point::distance(m[i], m[j]);
            if (min == -1 || min > d) {
                min = d;
            }
        }
    }
    double max = -1;
    for (int k=0; k<c; k++) {
        for (int i=0; i<n; i++) {
            if (k == p[i].getColor()) {
                for (int j=i+1; j<n; j++) {
                    if (k == p[j].getColor()) {
                        double d = Point::distance(p[i], p[j]);
                        if (max < d) {
                            max= d;
                        }
                    }
                }
            }
        }
    }
    cerr << "Dunn Index : " << min/max << endl;
}

void getCluster(Point *p, int n, int c) {
    for (int i=0; i<n; i++) {
        p[i].setColor(-1);
    }
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
    cerr << "Evaluation indices :: kmeans algorithm" << endl;
    printDBIndex(p, n, m, c);
    printDIndex(p, n, m, c);
    delete m;
}

void getClusterpp(Point *p, int n, int c) {
    for (int i=0; i<n; i++) {
        p[i].setColor(-1);
    }
    Point *m = new Point[c];
    genRandPoints(p, n, m, c);
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
    cerr << "Evaluation indices :: kmeans++ algorithm" << endl;
    printDBIndex(p, n, m, c);
    printDIndex(p, n, m, c);
    delete m;
}

int main(int argc, char const * argv[] ) {
    int n, c;
    if (argc < 5) {
        cerr << "Invalid number of algorithms" << endl;
        return 1;
    }
    n = stoi(argv[1]);
    c = stoi(argv[2]);
    ofstream f1(argv[3], ios::trunc);
    ofstream f2(argv[4], ios::trunc);
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
    clock_t s, e;
    s = clock(); 
    getCluster(p, n, c);
    e = clock();
    cout << "Time taken : " << (double)(e-s)/CLOCKS_PER_SEC << endl;
    for (int j = 0; j < n; ++j) {
        f1 << p[j].getX() << " " << p[j].getY() << " " << p[j].getColor() << endl;
    }
    s = clock();
    getClusterpp(p, n, c);
    e = clock();
    cout << "Time taken : " << (double)(e-s)/CLOCKS_PER_SEC << endl;
    for (int j = 0; j < n; ++j) {
        f2 << p[j].getX() << " " << p[j].getY() << " " << p[j].getColor() << endl;
    }
    delete p;
    return 0;
}
