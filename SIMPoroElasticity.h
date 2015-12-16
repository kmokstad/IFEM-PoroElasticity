// $Id$
//==============================================================================
//!
//! \file SIMPoroElasticity.h
//!
//! \date April 16 2015
//!
//! \author Yarded Bekele
//!
//! \brief Simulation driver for poroelasticity problems
//!
//==============================================================================

#ifndef SIMPOROELASTICITY_H_
#define SIMPOROELASTICITY_H_


#include "PoroElasticity.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "ASMbase.h"
#include "AlgEqSystem.h"
#include "ASMstruct.h"
#include "Functions.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Property.h"
#include "DataExporter.h"
#include "tinyxml.h"
#include "Vec3Oper.h"
#include "TimeStep.h"
#include "InitialConditionHandler.h"
#include "SIMSolver.h"
#include "Profiler.h"
#include "PoroMaterial.h"
#include <memory>


/*!
  \brief Driver class for poro-elastic simulators.
*/

template<class Dim> class SIMPoroElasticity : public Dim
{
public:
  //! \brief Dummy declaration, no setup properties needed
  typedef bool SetupProps;

  //! \brief The constructor initializes the references to the integrand.
  SIMPoroElasticity() : Dim({Dim::dimension,1}), poroel(Dim::dimension)
  {
    Dim::myProblem = &poroel;
  }

  //! \brief Destructor.
  virtual ~SIMPoroElasticity()
  {
    Dim::myProblem = nullptr;
    Dim::myInts.clear();
  }

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"poroelasticity"))
      return this->Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
    {
      if (!strcasecmp(child->Value(),"gravity"))
        this->parseGravity(child);
      else if (!strcasecmp(child->Value(),"scaling"))
        this->parseScaling(child);
      else if (!strcasecmp(child->Value(),"isotropic")) {
        int code = this->parseMaterialSet(child,mVec.size());
        std::cout <<"\tMaterial code "<< code <<":" << std::endl;
        mVec.push_back(std::unique_ptr<PoroMaterial>(new PoroMaterial));
        mVec.back()->parse(child);
        mVec.back()->printLog();
      }
      else
        this->Dim::parse(child);
    }

    if (!mVec.empty())
      poroel.setMaterial(mVec.front().get());

    return true;
  }

  //! \brief Parse gravity from an XML element.
  void parseGravity(const TiXmlElement* elem)
  {
    const char* value = utl::getValue(elem,"gravity");
    Vec3 grav;
    if (value)
    {
      std::stringstream str;
      str << value;
      str >> grav[0] >> grav[1] >> grav[2];
      poroel.setGravity(grav);
    }
    IFEM::cout << "Gravity vector: "
               << "\n\tg = " << grav << std::endl;
  }

  //! \brief Parse scaling parameter from an XML element.
  void parseScaling(const TiXmlElement* elem)
  {
    double sc = 1e14;
    utl::getAttribute(elem, "value", sc);

    poroel.setScaling(sc);
    IFEM::cout << "Scaling:\n\tsc = " << sc << std::endl;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);
    typename Dim::TracFuncMap::const_iterator tit = Dim::myTracs.find(propInd);

    if (vit != Dim::myVectors.end())
      poroel.setTraction(vit->second);
    else if (tit != Dim::myTracs.end())
      poroel.setTraction(tit->second);
    else
      return false;

    return true;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const
  {
    return "PoroElasticity";
  }

  //! \brief Obtain const reference to solution vector.
  //! \param[in] i Solution vector to get reference to.
  //! \return Const reference to requested vector.
  const Vector& getSolution(int i)
  {
    return solution[i];
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (Dim::opt.format < 0)
      return true;

    nBlock = 0;
    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] tp Time step identifier
  //! \param[in] nBlock Running VTF block counter
  bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (tp.step%Dim::opt.saveInc > 0 || Dim::opt.format < 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;

    // Write solution fields
    bool result = this->writeGlvS(solution.front(), iDump, nBlock,
                                  tp.time.t, true, "vector", 89);

    return result && this->writeGlvStep(iDump, tp.time.t);
  }

  //! \brief Initializes the simulator time stepping loop.
  bool init()
  {
    solution.resize(this->getNoSolutions());
    for (size_t i=0;i<solution.size();++i)
      solution[i].resize(this->getNoDOFs(),true);

    return true;
  }

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep& tp)
  {
    // Update vectors between time steps
    const int nNusols = solution.size();
    for (int n = nNusols-1; n > 0; n--)
      solution[n] = solution[n-1];

    return true;
  }

  //! \brief Computes the solution for the current time step.
  bool solveStep(TimeStep& tp)
  {
    if (Dim::msgLevel >= 0)
      IFEM::cout <<"\n  step = "<< tp.step <<"  time = "<< tp.time.t << std::endl;

    if (!this->assembleSystem(tp.time, solution))
      return false;

    if (!this->solveSystem(solution.front(), Dim::msgLevel-1,"displacement+pressure"))
      return false;

    return true;
  }

  //! \brief Register fields for data export
  void registerFields(DataExporter& exporter)
  {
    exporter.registerField("u,p","primary",DataExporter::SIM,DataExporter::PRIMARY);
    exporter.setFieldValue("u,p",this,&getSolution(0));
  }

  //! \brief Sets initial conditions.
  void setInitialConditions()
  {
    SIM::setInitialConditions(*this);
  }

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd)
  {
    if (propInd >= mVec.size())
      propInd = mVec.size()-1;

    poroel.setMaterial(mVec[propInd].get());
    return true;
  }

  //! \brief Prints a summary of the calculated solution to std::cout.
  //! \param[in] solution The solution vector
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual void printSolutionSummary(const Vector& solution, int, const char*,
                                    std::streamsize outPrec)
  {
    const size_t nsd = this->getNoSpaceDim();
    size_t iMax[nsd+1];
    double dMax[nsd+1];
    double dNorm = this->solutionNorms(solution,dMax,iMax,this->getNoFields(1));
    if (this->getNoFields(2) > 0) {
      double dpNorm = this->solutionNorms(solution,dMax+nsd,
                                          iMax+nsd,this->getNoFields(2), 'P');
      dNorm = sqrt(pow(dNorm,2.0)+pow(dpNorm,2.0));
    }

    std::stringstream str;
    if (this->adm.getProcId() == 0)
    {
      if (outPrec > 0) str.precision(outPrec);

      str <<"  Primary solution summary: L2-norm            : "<< utl::trunc(dNorm);

      char D = 'X';
      for (size_t d = 0; d < nsd; d++, D++)
        if (utl::trunc(dMax[d]) != 0.0)
          str <<"\n                            Max "<< char('X'+d)
              <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];
        str <<"\n                            Max pressure       : "
            << dMax[nsd] <<" node "<< iMax[nsd] << std::endl;
    }

    utl::printSyncronized(std::cout,str,this->adm.getProcId());
  }

private:
  PoroElasticity poroel;            //!< Poroelasticity integrand
  Vectors solution;                 //!< Solution vectors
  std::vector<std::unique_ptr<PoroMaterial>> mVec;  //!< Material data
};


//! \brief Partial specialization for configurator
template<class Dim>
struct SolverConfigurator< SIMPoroElasticity<Dim> > {
  int setup(SIMPoroElasticity<Dim>& pe,
            const typename SIMPoroElasticity<Dim>::SetupProps& props,
            char* infile)
  {
    utl::profiler->start("Model input");

    // Read input file
    if(!pe.read(infile))
      return 1;

    // Configure finite element library
    if(!pe.preprocess())
      return 2;

    // Setup integration
    pe.setQuadratureRule(pe.opt.nGauss[0],true);
    pe.initSystem(pe.opt.solver,1,1);
    pe.setAssociatedRHS(0,0);
    pe.setMode(SIM::DYNAMIC);

    // Time-step loop
    if (!pe.init())
      return 3;

    pe.setInitialConditions();

    return 0;
  }
};

#endif /* SIMPOROELASTICITY_H_ */
