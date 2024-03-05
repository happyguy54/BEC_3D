#ifndef C2_H
#define C2_H

// ROOT
#include<TF3.h>


class C2 {
  public:
    C2();
    virtual void setup(TF3*);

    virtual double operator() (double*, double*) = 0;

    double index;
    size_t nParams;
    double rejFrom;
    double rejTo;
    std::vector<double> rejFromOsl;
    std::vector<double> rejToOsl;
    double rej2From;
    double rej2To;
    bool useEps;
    bool fixC0;
};


class C2_1 : public C2 {
  public:
    C2_1();
    void setup(TF3*);

    double operator() (double*, double*);
};


class C2_2 : public C2 {
  public:
    C2_2();
    void setup(TF3*);

    double operator() (double*, double*);
};


class C2_3 : public C2 {
  public:
    C2_3();
    void setup(TF3*);

    double operator() (double*, double*);
};

class C2_4: public C2 {
  public:
    C2_4();
    void setup(TF3*);

    double operator() (double*, double*);
};

class C2_5: public C2 {
  public:
    C2_5();
    void setup(TF3*);

    double operator() (double*, double*);
};

class C2_6: public C2 {
  public:
    C2_6();
    void setup(TF3*);

    double operator() (double*, double*);
};

#endif /* C2 */
