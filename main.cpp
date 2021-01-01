#include <iostream>
#include<fstream>
#include <string>
#include<cmath>
#include<new>

using namespace std;

class Errors {
public:
    virtual void perr(ostream &out) {
    }
};

class CurErr : public Errors {
    double c, t, h;
public:
    CurErr(double c, double t, double h) {
        this->c = c;
        this->h = h;
        this->t = t;
    }

    void perr(ostream &out) override {
        cout << "для c,t и h должно выполняться c*t/h<1";
    }
};

class ErrNew : public Errors {
    bad_alloc x;
public:
    ErrNew(bad_alloc x) {
        this->x = x;
    }

    void perr(ostream &out) override {
        cout << x.what();
    }
};

class Scheme {
protected:
    int T, L;
    double t, h, c;
public:
    Scheme(double t, int T, double h, int L, double c = 1) {
        try {
            this->t = t;
            this->T = T;
            this->h = h;
            this->L = L;
            this->c = c;
            if (c * t / h > 1) {
                CurErr err(c, t, h);
                throw err;
            }
        }
        catch (CurErr err) {
            err.perr(cout);
        }
    }
};

class CIR : public Scheme {
private:
    double *U;
public:
    CIR(double const *U, double t, int T, double h, int L, double c = 1) : Scheme(t, T, h, L, c) {
        try {
            this->U = new double[L];
            for (int i = 0; i < L; i++)
                this->U[i] = U[i];
        }
        catch (bad_alloc x) {
            ErrNew err(x);
            err.perr(cout);
        }
    }

    ~CIR() {
        delete[] U;
    }

    void makefiles(string name = "", string path = "") {
        U[0] = 0;
        double a = 0;
        double *U_old = new double[L];
        for (int i = 0; i < L; i++)
            U_old[i] = U[i];
        ofstream fout;
        fout.open(path + name + ("0.csv"), ios::out);
        fout.clear();
        for (int i = 0; i < L; i++) {
            fout << i * h << ' ' << U[i] << '\n';
        }
        fout.close();
        for (int j = 1; j < T; j++) {
            fout.open(path + name + to_string(j) + ".csv", ios::out);
            fout.clear();
            fout << 0 << ' ' << 0 << '\n';
            for (int i = 1; i < L; i++) {
                U[i] = U_old[i] - c * t / h * (U_old[i] - U_old[i - 1]);
                fout << i * h << ' ' << round(U[i] * 1000) / 1000 << '\n';
            }
            for (int i = 0; i < L; i++)
                U_old[i] = U[i];
            fout.close();
        }
    }
};

class LW : public Scheme {
private:
    double *U;
public:
    LW(double const *U, double t, int T, double h, int L, double c = 1) : Scheme(t, T, h, L, c) {
        try {
            this->U = new double[L];
            for (int i = 0; i < L; i++)
                this->U[i] = U[i];
        }
        catch (bad_alloc x) {
            ErrNew err(x);
            err.perr(cout);
        }
    }

    ~LW() {
        delete[] U;
    }

    void makefiles(string name,string path = "") {
        U[0] = 0;
        U[L - 1] = 0;
        double *U_old = new double[L];
        for (int i = 0; i < L; i++)
            U_old[i] = U[i];
        ofstream fout;
        fout.open(path+name + ("0.csv"), ios::out);
        fout.clear();
        for (int i = 0; i < L; i++) {
            fout << i * h << ' ' << U[i] << '\n';
        }
        fout.close();
        for (int j = 1; j < T; j++) {
            fout.open(path+name + to_string(j) + ".csv", ios::out);
            fout.clear();
            fout << 0 << ' ' << 0 << '\n';
            for (int i = 1; i < L - 1; i++) {
                U[i] = U_old[i] - c * t / h * (U_old[i] - U_old[i - 1]) +
                       c * c * t * t / 2 / h / h * (U_old[i + 1] - 2 * U_old[i] + U_old[i - 1]);
                fout << i * h << ' ' << round(U[i] * 1000) / 1000 << '\n';
            }
            for (int i = 0; i < L; i++)
                U_old[i] = U[i];
            fout << 0;
            fout.close();
        }
    }
};

class TVD : public Scheme {
private:
    double *U;

    double MClimiter(double a0, double a1, double a2) {
        double ans = min(min(abs((a2 - a0) / 2),  2*abs(a2 - a1)),  2*abs(a1 - a0));
        if (a2 - a1 > 0)
            return ans;
        else return -ans;
    }
    double monmodlimiter(double a0, double a1, double a2){
        double ans = min(abs(a2 - a1), abs(a1 - a0));
        if (a2 - a1 > 0)
            return ans;
        else return -ans;
    }
    double Superbeelimiter(double a0, double a1, double a2){
        double ans = max(min(2*abs(a2-a1),abs(a1-a0)),min(abs(a2-a1),2*abs(a1-a0)));
        return ans;
    }
    double vanLeerlimiter(double a0, double a1, double a2){
        if(a2!=a0){
            return 2*(a2-a1)*(a1-a0)/(a2-a0);
        }
        else return MClimiter(a0,a1,a2);
    }

public:
    TVD(double const *U, double t, int T, double h, int L, double c = 1) : Scheme(t, T, h, L, c) {
        try {
            this->U = new double[L];
            for (int i = 0; i < L; i++)
                this->U[i] = U[i];
        }
        catch (bad_alloc x) {
            ErrNew err(x);
            err.perr(cout);
        }
    }

    ~TVD() {
        delete[] U;
    }

    void makefiles(string name = "", string path = "") {
        U[0] = 0;
        U[1] = 0;
        U[L - 1] = 0;
        U[L - 2] = 0;
        double *U_old = new double[L];
        for (int i = 0; i < L; i++)
            U_old[i] = U[i];
        ofstream fout;
        fout.open(path + name + ("0.csv"), ios::out);
        fout.clear();
        for (int i = 0; i < L; i++) {
            fout << i * h << ' ' << U[i] << '\n';
        }
        fout.close();
        for (int j = 1; j < T; j++) {
            fout.open(path + name + to_string(j) + ".csv", ios::out);
            fout.clear();
            fout << 0 << ' ' << 0 << '\n';
            for (int i = 2; i < L - 2; i++) {
                double interstep1, interstep2, lim2, curr;
                curr = c * t / h;
                if (curr > 0) {
                    interstep2 = U_old[i] + (1 - curr) / 2 * MClimiter(U_old[i - 1], U_old[i], U_old[i + 1]);
                    interstep1 = U_old[i - 1] + (1 - curr) / 2 * MClimiter(U_old[i - 2], U_old[i - 1], U_old[i]);
                } else {
                    interstep2 = U_old[i] - (1 - curr) / 2 * MClimiter(U_old[i], U_old[i + 1], U_old[i + 2]);
                    interstep1 = U_old[i - 1] - (1 - curr) / 2 * MClimiter(U_old[i - 1], U_old[i], U_old[i + 1]);
                }
                U[i] = U_old[i] - c * t / h * (interstep2 - interstep1);
                fout << i * h << ' ' << round(U[i] * 1000) / 1000 << '\n';
            }
            for (int i = 0; i < L; i++)
                U_old[i] = U[i];
            fout << 0;
            fout.close();
        }
    }
};

int main() {
    int L = 200;
    double *U;
    double *U1;
    try {
        U = new double[L];
        for (int i = 0; i < 50; i++)
            U[i] = 0;
        for (int i = 50; i < 100; i++)
            U[i] = 1;
        for (int i = 100; i < 200; i++)
            U[i] = 0;
    }
    catch (bad_alloc x) {
        ErrNew err(x);
        err.perr(cout);
        return 1;
    }
    TVD c(U, 0.5, 400, 1, L);
    c.makefiles("TVD");
    delete[] U;
    return 0;
}
