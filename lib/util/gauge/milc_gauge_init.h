// -*- C++ -*-
// $Id: milc_gauge_init.h,v 3.1 2007-02-04 22:06:42 edwards Exp $
/*! \file
 *  \brief Read a MILC config
 */

#ifndef __milc_gauge_init_h__
#define __milc_gauge_init_h__

#include "util/gauge/gauge_init.h"

namespace Chroma
{

  //! Name and registration
  namespace MILCGaugeInitEnv
  {
    extern const std::string name;
    bool registerAll();
  

    //! Params for initializing config
    /*! @ingroup gauge */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      string cfg_file;		/*!< File name */
    };



    //! Gauge initialization
    /*! @ingroup gauge
     *
     * MILC reader
     */
    class GaugeIniter : public GaugeInit
    {
    public:
      //! Full constructor
      GaugeIniter(const Params& p) : params(p) {}

      //! Initialize the gauge field
      void operator()(XMLReader& gauge_file_xml,
		      XMLReader& gauge_xml,
		      multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      GaugeIniter() {}

    private:
      Params  params;
    };

  }  // end namespace


  //! Reader
  /*! @ingroup gauge */
  void read(XMLReader& xml, const string& path, MILCGaugeInitEnv::Params& param);

  //! Writer
  /*! @ingroup gauge */
  void write(XMLWriter& xml, const string& path, const MILCGaugeInitEnv::Params& param);

}  // end namespace Chroma


#endif
