#include <algorithm>
#include <fstream>
#include <iostream>

#include "matrix.h"
#include "parameters.h"

using namespace std;

Parameters::Parameters(const string inputfile){
  ifstream in;
  string var, eqs;
  real val;

  in.open(inputfile);
  if (in.is_open()){
    while (in >> var >> eqs >> val){
      if ( var.compare("xlength") == 0 ){
        xlength = val;}
      else if ( var.compare("ylength") == 0 ){
        ylength = val;}
      else if ( var.compare("imax") == 0 ){
        imax = val;}
      else if ( var.compare("jmax") == 0 ){
        jmax = val;}
      else if ( var.compare("dt") == 0 ){
        dt = val;}
      else if ( var.compare("t_end") == 0 ){
        t_end = val;}
      else if ( var.compare("dt_out") == 0 ){
        dt_out = val;}
      else if ( var.compare("tau") == 0 ){
        tau = val;}
      else if ( var.compare("itermax") == 0 ){
        itermax = val;}
      else if ( var.compare("eps") == 0 ){
        eps = val;}
      else if ( var.compare("omega") == 0 ){
        omega = val;}
      else if ( var.compare("Re") == 0 ){
        Re = val;}
      else if ( var.compare("gx") == 0 ){
        gx = val;}
      else if ( var.compare("gy") == 0 ){
        gy = val;}
      else if ( var.compare("UI") == 0 ){
        UI = val;}
      else if ( var.compare("VI") == 0 ){
        VI = val;}
      else if ( var.compare("PI") == 0 ){
        PI = val;}
      else if ( var.compare("rho") == 0 ){
        rho = val;}
      else if ( var.compare("dynvis") == 0 ){
        vis = val;}
      else if ( var.compare("inflow") == 0 ){
        inflow = val;}
      else if ( var.compare("wW") == 0 ){
        wW = val;}
      else if ( var.compare("wE") == 0 ){
        wE = val;}
      else if ( var.compare("wN") == 0 ){
        wN = val;}
      else if ( var.compare("wS") == 0 ){
        wS = val;}
      else{
        cout << "The parameter " + var + " is invalid." << endl;
      }
    };
  in.close();
  }else{
    cout << "Did not manage to open the input file!" << endl;
  }
  dx = xlength/(real)(imax);
  dy = ylength/(real)(jmax);
}

void Parameters::setScale(const real uMax, const real vMax, const real pMax){
  L = std::max(xlength,ylength);
  uInf = std::max(uMax,vMax);
  pInf = pMax;
  rhoInf = rho;
  Re = rhoInf*uInf*L/vis;
  std::cout << "The Reynold number is: " << Re << std::endl;
}

void Parameters::toDimensionless(matrix<real>* U, matrix<real>* V, matrix<real>* P){
  xlength = xlength/L;
  ylength = ylength/L;
  t_end = uInf*t_end/L;
  dt = uInf*dt/L;
  inflow = inflow/uInf;
  for (int i = 0; i<imax+2; ++i){
    for (int j = 1; j<jmax+2; j++){
      U->set(i,j,U->get(i,j)/uInf);
      V->set(i,j,V->get(i,j)/uInf);
      P->set(i,j,(P->get(i,j)-pInf)/(rhoInf*uInf*uInf));
    }
  }
}

void Parameters::toDimensional(matrix<real>* U, matrix<real>* V, matrix<real>* P){
  xlength = xlength*L;
  ylength = ylength*L;
  t_end = t_end*L/uInf;
  dt = dt*L/uInf;
  inflow = inflow*uInf;
  for (int i = 0; i<imax+2; ++i){
    for (int j = 1; j<jmax+2; j++){
      U->set(i,j,U->get(i,j)*uInf);
      V->set(i,j,V->get(i,j)*uInf);
      P->set(i,j,P->get(i,j)*rhoInf*uInf*uInf+pInf);
    }
  }
}
