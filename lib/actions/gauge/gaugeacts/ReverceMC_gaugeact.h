// -*- C++ -*-
/*! \file
 *  \brief  gauge action for using continuous beta method of ReverceMC
 */

#ifndef __ReverceMC_gaugeact_h__
#define __ReverceMC_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace ReverceMCGaugeActEnv 
  { 
    bool registerAll();

    //! Parameter structure
    /*! @ingroup gaugeacts */
    struct Params 
    {
      //! Base Constructor
      Params() {}
	//I is index matrix for parsing beta lattice.
	//Values assigned based on symmetric matrix.
      void InitI(){
	I[0][1] = 0;
	I[0][2] = 1;
	I[0][3] = 2;
	I[1][2] = 3;
	I[1][3] = 4;
	I[2][3] = 5;
	I[1][0] = I[0][1];
	I[2][0] = I[0][2];
	I[3][0] = I[0][3];
	I[2][1] = I[1][2];
	I[3][1] = I[1][3];
	I[3][2] = I[2][3];
      }
      
      //Declare lattice of betas.
      //multi1d<LatticeReal> beta(Nd*(Nd-1)/2);
      //Nd not recognized at compile time.
      multi1d<LatticeReal> beta;
      bool Normalization;
      //Flag for determining whether normalization term is turned on/off.
      int I[Nd][Nd];
      
      //! Read params from some root path
      Params(XMLReader& xml_in, const std::string& path);

      //These are probably no longer needed.
      //Real  beta_F;   // Coupling for fundamental plaquette
      //Real  beta_A;   // Coupling for adjoint plaquette
      //Real  beta_S;   // Coupling for sextet plaquette
      //Only betas relating to the ReverceMC will be used.
      Real  beta_initial;   // An initial beta that may be restored.
      Real  alpha;          // Lower limit of integration for probability distribution.
      //alpha is a paramater that may be tuned for testing.

    };
  

    //! ReverceMC  gauge action
    /*! \ingroup gaugeacts
     *
     * The standard gauge action.
     */
    class GaugeAct : public LinearGaugeAction
    {
    public:
      // Typedefs to save typing
      typedef multi1d<LatticeColorMatrix>  P;
      typedef multi1d<LatticeColorMatrix>  Q;

      //! General CreateGaugeState<P,Q>
      //! Read coeff from a param struct
      GaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, const Params& p) :
	cgs(cgs_), param(p) {}

      //! Return the set on which the gauge action is defined
      /*! Defined on the even-off (red/black) set */
      const Set& getSet() const {return rb;}

      //! Compute staple
      /*! Default version. Derived class should override this if needed. */
      void staple(LatticeColorMatrix& result,
		  const Handle< GaugeState<P,Q> >& state,
		  int mu, int cb) const;

      //This function will be called in the monomials and set the betas.
      void setbeta(const Handle< GaugeState<P,Q> >& state, const multi1d<LatticeReal> & input);
      //void setbeta(multi1d<LatticeColorMatrix>& result,
      		 //const Handle< GaugeState<P,Q> >& state) const;
      //Also used in the monomials, sets beta lattice to one value.
      void restorebeta();
      //! Compute dS/dU
      void deriv(multi1d<LatticeColorMatrix>& result,
		 const Handle< GaugeState<P,Q> >& state) const;

      //! Compute the actions
      Double S(const Handle< GaugeState<P,Q> >& state) const;

      //! Produce a gauge create state object
      const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

      //! Destructor is automatic
      ~GaugeAct() {}

    protected:
      //! Hide assignment
      void operator=(const GaugeAct& a) {}
      
      //! Compute the site-level action
      void siteAction(multi2d<LatticeColorMatrix>& site_act, 
		      const Handle< GaugeState<P,Q> >& state) const;

    private:
      Handle< CreateGaugeState<P,Q> >  cgs;  /*!< Create Gauge State */
      Params               param;            /*!< The parameters */
    };

  }

}


#endif
