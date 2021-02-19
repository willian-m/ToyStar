//Abstract class. Other EOS equations must inherit it and implement actual code.

#ifndef EOSBASE_H
#define EOSBASE_H

class EOSBase{

    
public:
    virtual double get_pressure(double density) = 0;
};

#endif