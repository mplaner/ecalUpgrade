#ifndef utils_h
#define utils_h

#include </usr/include/math.h>

#include "init.h"

float dR(double eta1, double phi1, double eta2, double phi2);
float dRetaphi(double eta1, double phi1, double eta2, double phi2);
float SignalShapeEB(float x);
float SignalShapeEBfit(double *x, double *par);
float SignalShapeEE(float x);
float SignalShapeEEfit(double *x, double *par);
float invMass(float e1, float eta1, float phi1, float e2, float eta2, float phi2);
//vector<float> defineSmallestIntervalToContain(float interval, vector<float> v);
float defineSmallestIntervalToContain(float interval, vector<float> v);
void sortKeepAligned(vector<float> &v1, vector<float> &v2);

float defineAverage(vector<float> v);


#endif 
