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
    std::vector<double> alpha;
    bool useEps;
    bool fixC0;
    bool with_errors;
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

class C2_7: public C2 {
  public:
    C2_7();
    void setup(TF3*);

    double operator() (double*, double*);
};

class C2_8: public C2 {
  public:
    C2_8();
    void setup(TF3*);

    double operator() (double*, double*);
};

class C2_9: public C2 {
  public:
    C2_9();
    void setup(TF3*);

    double operator() (double*, double*);
};

class C2_10: public C2 {
  public:
    C2_10();
    void setup(TF3*);

    double operator() (double*, double*);
};

#endif /* C2 */
