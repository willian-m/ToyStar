#ifndef EOSPOLYTROPIC_H
#define EOSPOLYTROPIC_H

#include<math.h>
#include "eos_base.h"

class EOSPolytropic : public EOSBase{

private:
    EOSPolytropic();
    double k;
    double n;

public:
    EOSPolytropic(double lk, double ln);
    double  get_pressure(double density);
};

inline double EOSPolytropic::get_pressure(double rho){ return k*std::pow(rho,1.+1./n); };

#endif