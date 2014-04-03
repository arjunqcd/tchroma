/*! \file
 *  \Plaquette gauge action as extension of Reverce Monte Carlo method
 */

#include "chromabase.h"
#include "actions/gauge/gaugeacts/ReverceMC_gaugeact.h"
#include "actions/gauge/gaugeacts/gaugeact_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"
#include "io/aniso_io.h"
#include <math.h>
//Include math so exponents and logs may be computed.

namespace Chroma
{
 
  namespace ReverceMCGaugeActEnv
  { 
    namespace
    {
      GaugeAction< multi1d<LatticeColorMatrix>, 
		   multi1d<LatticeColorMatrix> >* createGaugeAct(XMLReader& xml, 
								 const std::string& path) 
      {
	return new GaugeAct(CreateGaugeStateEnv::reader(xml, path), 
			    Params(xml, path));
      }
      
      const std::string name = "REVERCEMC_GAUGEACT";

      //! Local registration flag
      static bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeActFactory::Instance().registerObject(name, createGaugeAct);
	registered = true;
      }
      return success;
    }


    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      XMLReader paramtop(xml_in, path);
      InitI();
      /*param.I[0][1] = 0;
      param.I[0][2] = 1;
      param.I[0][3] = 2;
      param.I[1][2] = 3;
      param.I[1][3] = 4;
      param.I[2][3] = 5;
      param.I[1][0] = I[0][1];
      param.I[2][0] = I[0][2];
      param.I[3][0] = I[0][3];
      param.I[2][1] = I[1][2];
      param.I[3][1] = I[1][3];
      param.I[3][2] = I[2][3];*/
      Normalization = false;
      beta = multi1d<LatticeReal> (Nd*(Nd-1)/2);
      //beta = multi1d<LatticeReal> (6);
      try 
      {	
	//read(paramtop, "beta_F", beta_F);
	//read(paramtop, "beta_A", beta_A);
	//read(paramtop, "beta_S", beta_S);
    	  //New inputs for to read from XML.
    	  read(paramtop, "beta_initial", beta_initial);
    	  read(paramtop, "alpha", alpha);
	  //Initialize beta matrix as ordinary action. 
	  beta = beta_initial;
      }
      catch( const std::string& e ) { 
	QDPIO::cerr << "Error reading XML: " <<  e << endl;
	QDP_abort(1);
      }
    }


    //! Compute the action
    Double GaugeAct::S(const Handle< GaugeState<P,Q> >& state) const
    {
      // Action at the site level
      multi2d<LatticeColorMatrix> plq;
      this->siteAction(plq, state);

      // Total action
      Double act_F  = zero;
      Double LogTerm = zero ;
      Double three = Nc;
      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  //LatticeComplex plaq=trace(plq[mu][nu]);
	  //Used as template for new action.
	    LatticeReal plaq=real(trace(plq[mu][nu])-three);
	  // Sum over plaquettes
	  //New action below, including switch for normalization term.
		act_F += sum(param.beta[param.I[mu][nu]]*plaq);
	  if(param.Normalization == true)
	    LogTerm += sum(log(1/plaq*(exp((param.beta_initial - param.alpha)*plaq/Real(Nc)) - 1)));
	}
      }

      // Normalize
      Real act = -(1 / Real(Nc)) * act_F + LogTerm ;
      //The normalization is edited to adjust for Reverce MC.

      return act;
    }
 

    //! Compute the plaquette
    void GaugeAct::siteAction(multi2d<LatticeColorMatrix>& plq, const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      // Initialize
      plq.resize(Nd,Nd);
      plq = zero;

      // Handle< const GaugeState<P,Q> > u_bc(createState(u));
      // Apply boundaries
      const multi1d<LatticeColorMatrix>& u = state->getLinks();

      // Compute the average plaquettes
      for(int mu=1; mu < Nd; ++mu)
      {
	for(int nu=0; nu < mu; ++nu)
	{
	  /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
	  /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu) */
	  /* wplaq_tmp = tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)) */
		plq[mu][nu] += (u[mu]*shift(u[nu],FORWARD,mu)*adj(shift(u[mu],FORWARD,nu))*adj(u[nu])) ;
		//Plaquette should still be calculated normally.
	  
	  // Keep a copy
	  plq[nu][mu] = plq[mu][nu];
	}
      }


      END_CODE();
    }
 

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void GaugeAct::staple(LatticeColorMatrix& result,
			  const Handle< GaugeState<P,Q> >& state,
			  int mu, int cb) const
    {
      QDPIO::cerr << __func__ << ": staple not possible\n";
      QDP_abort(1);
    }

    //Set the betas for Reverce MC.
    //void setbeta(const multi1d<LatticeReal> & input)
    void GaugeAct::setbeta(const Handle< GaugeState<P,Q> >& state, const multi1d<LatticeReal> & input)
    {
    	param.beta = input;
    	param.Normalization = true;
    }
    //void restorebeta()
    void GaugeAct::restorebeta()
    {
    	param.beta = param.beta_initial;
    	param.Normalization = false;
    }
    //! Compute dS/dU
    void GaugeAct::deriv(multi1d<LatticeColorMatrix>& ds_u,
			 const Handle< GaugeState<P,Q> >& state) const
    {
      START_CODE();

      LatticeColorMatrix tmp_1;
      LatticeColorMatrix tmp_2;

      const multi1d<LatticeColorMatrix>& u = state->getLinks();
      multi1d<LatticeColorMatrix> deriv_fun(Nd); 
      //multi1d<LatticeColorMatrix> deriv_adj(Nd);
      //multi1d<LatticeColorMatrix> deriv_sex(Nd);
      ds_u.resize(Nd);
      for(int mu=0; mu < Nd; mu++) 
      {
	deriv_fun[mu] = zero ;
	//deriv_adj[mu] = zero ;
	//deriv_sex[mu] = zero ;
	for(int nu = 0; nu < Nd; nu++) 
	{ 
	  if (mu == nu) continue;
	  //Here will be the cases to call the beta lattice.
	  LatticeColorMatrix tmp_1 = shift(u[nu], FORWARD, mu);
	  LatticeColorMatrix tmp_2 = shift(u[mu], FORWARD, nu);


	  //These lines have been modified to include multiplying by beta for ReverceMC.
	  LatticeColorMatrix up_plq   = (-1/Real(2*Nc))*param.beta[param.I[mu][nu]]*u[mu]*tmp_1*adj(tmp_2)*adj(u[nu]);
	  LatticeColorMatrix adjtmp_1 = adj(tmp_1);
	  LatticeColorMatrix down_plq = u[mu]*shift((-1/Real(2*Nc))*param.beta[param.I[mu][nu]]*adjtmp_1/*adj(tmp_1)*/*adj(u[mu])*u[nu],BACKWARD,nu);

	  //The original derivative is replaced with dS/dU of the ReverceMC action.
	  deriv_fun[mu] += up_plq + down_plq ;
	  if(param.Normalization == true)
	  {
		  LatticeReal Action = real(trace(up_plq));
		  deriv_fun[mu] += up_plq*(1/Action - (param.beta_initial - param.alpha)/(exp((param.alpha - param.beta_initial)*Action/((-1/Real(2*Nc))*param.beta[param.I[mu][nu]])) - 1)/param.beta[param.I[mu][nu]]);
		  Action = real(trace(down_plq));
		  deriv_fun[mu] += down_plq*(1/Action - (param.beta_initial - param.alpha)/(exp((param.alpha - param.beta_initial)*Action/((-1/Real(2*Nc))*param.beta[param.I[mu][nu]])) - 1)/param.beta[param.I[mu][nu]]);
	  }

	  //Adjoint and sextet not used.
	  //deriv_adj[mu] += up_plq*conj(trace(up_plq)) +
	  //                 down_plq*conj(trace(down_plq)) ;

	  //deriv_sex[mu] += up_plq*up_plq + down_plq*down_plq +
	    //               up_plq*trace(up_plq) +
	      //             down_plq*trace(down_plq) ;
	
	}// nu
	// Fold in the normalization from the action
	//Only fundamental used below.
	ds_u[mu]  = deriv_fun[mu];
	//ds_u[mu]  = (-1/Real(2*Nc)        ) * deriv_fun[mu];
	//ds_u[mu] += (-param.beta_A/Real(Nc*Nc)       ) * deriv_adj[mu];
	//ds_u[mu] += (-param.beta_S/Real(2*(Nc*Nc+Nc))) * deriv_sex[mu];
      }// mu

      // Zero the force on any fixed boundaries
      getGaugeBC().zero(ds_u);

      END_CODE();
    }


  }//ReverceMCGaugeActEnv

} // Chroma
