#include <fstream>
#include <iostream>

#include "parameters.h"

using namespace std;

Parameters::Parameters(string inputfile){
  ifstream in;
  string line;
  string var, eqs;
  Real val;

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
      else if ( var.compare("wW") == 0 ){
        wW = val;}
      else if ( var.compare("wE") == 0 ){
        wE = val;}
      else if ( var.compare("wN") == 0 ){
        wN = val;}
      else if ( var.compare("wS") == 0 ){
        wS = val;}
      else{
        cout << "The parameter " << var << " is invalid." << endl;
      }
    };
  in.close();
  }else{
    cout << "Did not manage to open the input file!" << endl;
  }
  dx = xlength/imax;
  dy = ylength/jmax;
}
