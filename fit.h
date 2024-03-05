#ifndef FIT_H
#define FIT_H

// ROOT
#include<TH3.h>
#include<TF3.h>


TH3D* hist;
TF3* func;


void MakeHistogram();
void Fit();
void Output();
void Plot();

#endif /* FIT_H */
