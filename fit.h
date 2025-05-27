#ifndef FIT_H
#define FIT_H

// ROOT
#include<TH3.h>
#include<TF3.h>


TH3D* hist;
TH3D* histLarge;
TF3* func;
TF3* funcInit;


void MakeHistogram();
void Fit();
void Output();
void Plot();

#endif /* FIT_H */
