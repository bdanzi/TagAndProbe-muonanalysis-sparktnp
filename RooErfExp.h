#ifndef ROOERFEXP
#define ROOERFEXP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "TMath.h"
#include "RooMath.h"

class RooErfExp : public RooAbsPdf
{
 public:
   RooErfExp() {};
   RooErfExp(const char *name, const char *title,
	     RooAbsReal& _x,
	     RooAbsReal& _alpha,
	     RooAbsReal& _gamma,
	     RooAbsReal& _beta,
	     RooAbsReal& _n);
   
   RooErfExp(const RooErfExp& other, const char* name);
   inline virtual TObject* clone(const char* newname) const { return new RooErfExp(*this,newname); }
   inline ~RooErfExp() {}
   Double_t evaluate() const;
   
   ClassDef(RooErfExp, 2);
   
 protected:
   
   RooRealProxy x;
   RooRealProxy alpha;
   RooRealProxy gamma;
   RooRealProxy beta;
   RooRealProxy n;
};

#endif
