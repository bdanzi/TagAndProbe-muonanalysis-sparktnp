#include "RooErfExp.h"

ClassImp(RooErfExp)

RooErfExp::RooErfExp(const char *name, const char *title, 
		     RooAbsReal& _x,
		     RooAbsReal& _alpha,
		     RooAbsReal& _gamma,
		     RooAbsReal& _beta,
		     RooAbsReal& _n):
RooAbsPdf(name,title), 
  x("x","x",this,_x),
  alpha("alpha","alpha",this,_alpha),
  gamma("gamma","gamma",this,_gamma),
  beta("beta","beta",this,_beta),
  n("n","n",this,_n)
{ }

RooErfExp::RooErfExp(const RooErfExp& other, const char* name):
RooAbsPdf(other,name),
  x("x",this,other.x),
  alpha("alpha",this,other.alpha),
  gamma("gamma",this,other.gamma),
  beta("beta",this,other.beta),
  n("n",this,other.n)
{ }

Double_t RooErfExp::evaluate() const
{
   return RooMath::erfc((x-alpha)*beta) * exp(-gamma*(pow(x, n) - pow(alpha, n)));
}

