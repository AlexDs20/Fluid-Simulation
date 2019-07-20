#include <fstream>
#include <iostream>
#include "Parameters.h"

using namespace std;

Parameters::Parameters(string inputfile){
  ifstream in;
  string line;
  string var, eqs;
  Real val;

  in.open(inputfile);
  if (in.is_open()){
    while (in >> var >> eqs >> val){
      if ( var.compare("xlength") == 0 )
        xlength = val;
      if ( var.compare("ylength") == 0 )
        ylength = val;
      if ( var.compare("imax") == 0 )
        imax = val;
      if ( var.compare("jmax") == 0 )
        jmax = val;
      if ( var.compare("dt") == 0 )
        dt = val;
      if ( var.compare("t_end") == 0 )
        t_end = val;
      if ( var.compare("tau") == 0 )
        tau = val;
      if ( var.compare("itermax") == 0 )
        itermax = val;
      if ( var.compare("eps") == 0 )
        eps = val;
      if ( var.compare("omega") == 0 )
        omega = val;
      if ( var.compare("Re") == 0 )
        Re = val;
      if ( var.compare("gx") == 0 )
        gx = val;
      if ( var.compare("gy") == 0 )
        gy = val;
      if ( var.compare("UI") == 0 )
        UI = val;
      if ( var.compare("VI") == 0 )
        VI = val;
      if ( var.compare("PI") == 0 )
        PI = val;
      if ( var.compare("wW") == 0 )
        wW = val;
      if ( var.compare("wE") == 0 )
        wE = val;
      if ( var.compare("wN") == 0 )
        wN = val;
      if ( var.compare("wS") == 0 )
        wS = val;
    }
  in.close();
  }
  dx = xlength/imax;
  dy = ylength/jmax;
}
