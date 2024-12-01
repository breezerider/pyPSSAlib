#include <memory>     // STL C-style pointer wrapper
#include <atomic>     // STL atomic definitions
#include <thread>

// pSSAlib
#include <pssalib/PSSA.h>
#include <pssalib/util/Timing.h>

#define PYBIND11_DETAILED_ERROR_MESSAGES

// Python bindings
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "utils.h"
#include "models.h"


namespace py = pybind11;


class pSSAlibWrapper
{
public:
    enum tagTestCase
    {
        tcCLC,
        tcCA,
        tcHomo,
        tcSBD,
        tcTCS,
        tcED,
        tcNone
    };

protected:

  /**
   * Initialize the data model for a test case
   *
   * @param model reference to a @link pssalib::datamodel::detail::Model object
   * @return @true if initialization succeeds, @false otherwise
   */
  bool initializeModel(pssalib::datamodel::detail::Model &model,
                       enum tagTestCase testCase,
                       double * arParams, size_t szParams) const
  {
    switch(testCase)
    {
    case tcCLC:
      generateCyclicLinearChain(model, arParams, szParams);
      break;
    case tcCA:
      generateColloidalAggregation(model, arParams, szParams);
      break;
    case tcHomo:
      generateHomoreaction(model, arParams, szParams);
      break;
    case tcSBD:
      generateSingleBirthDeath(model, arParams, szParams);
      break;
    case tcTCS:
      generateTwoComponentSystem(model, arParams, szParams);
      break;
    case tcED:
      generateEnzymaticDegradation(model, arParams, szParams);
      break;
    default:
      PY_ERRMSG("Invalid test case identifier")
      return false;
    }

    return true;
  }

  /**
   * Simulation thread main function
   *
   * @param ptrSimInfo pointer to a @link pssalib::datamodel::SimulationInfo
   * object
   * @return void
   */
  void runSimulation(pssalib::datamodel::SimulationInfo *ptrSimInfo)
  {
    pssalib::PSSA simEngine;

    simEngine.setMethod(method);

    simEngine.run(ptrSimInfo);

    bProcessing = false;
  }

  /**
   * Dispatcher function that starts the simulation thread.
   *
   * @param testCase one of @link tagTestCase
   * @param params one-dimensioanl array of model parameters
   * @param samples number of samples to be generated
   * @param timeEnd final simulation time
   * @param timeStart initial output time
   * @param timeStep time interval between subsequent outputs
   * @return a @link py::array_t<unsigned> containing the results
   */
  py::array_t<unsigned int>
  run(enum tagTestCase testCase,
      py::array_t<double, py::array::c_style | py::array::forcecast> params,
      size_t samples, double timeEnd, double timeStart, double timeStep) {
    if (bProcessing)
      throw std::runtime_error("Another simualtion is running already");
    bProcessing = true;

    try {
      if (samples < 1)
        throw std::runtime_error(
            "Number of samples must be a positive integer");

      py::buffer_info buf_params = params.request();

      if (buf_params.ndim != 1)
        throw std::runtime_error("Parameters must be a float vector");

      size_t numTimepoints;
      if (timeStart < 0.0)
        numTimepoints = 1;
      else {
        numTimepoints =
            pssalib::timing::getNumTimePoints(timeStart, timeEnd, timeStep);

        if (numTimepoints <= 0)
          throw std::runtime_error(
              "Time interval must include at least one time point");
      }

      // input parameters vector
      double *ptr_params = static_cast<double *>(buf_params.ptr);
      size_t num_params = buf_params.size;

      // intialize simulation info
      pssalib::datamodel::SimulationInfo simInfo;

      if (!initializeModel(simInfo.getModel(), testCase, ptr_params,
                           num_params))
        throw std::runtime_error("Failed to intialize test case model");

      // setup the simulation params
      simInfo.unSamplesTotal = samples;
      simInfo.dTimeEnd = timeEnd;
      simInfo.dTimeStep = timeStep;

      if (timeStart < 0.0) {
        simInfo.dTimeStart = timeEnd;
        simInfo.unOutputFlags =
            pssalib::datamodel::SimulationInfo::ofRawFinalPops;
      } else {
        simInfo.dTimeStart = timeStart;
        simInfo.unOutputFlags =
            pssalib::datamodel::SimulationInfo::ofRawTrajectory;
      }

      // setup output
      simInfo.unOutputFlags |=
          pssalib::datamodel::SimulationInfo::ofNone
          // | pssalib::datamodel::SimulationInfo::ofTrajectory
          // | pssalib::datamodel::SimulationInfo::ofStatus
          | pssalib::datamodel::SimulationInfo::ofLog
          // | pssalib::datamodel::SimulationInfo::ofTrace
          // | pssalib::datamodel::SimulationInfo::ofInfo
          // | pssalib::datamodel::SimulationInfo::ofWarning
          | pssalib::datamodel::SimulationInfo::ofError
          // | pssalib::datamodel::SimulationInfo::eofModuleGrouping
          // | pssalib::datamodel::SimulationInfo::eofModuleSampling
          // | pssalib::datamodel::SimulationInfo::eofModuleUpdate
          ;

      simInfo.resetOutput();

      simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog,
                                 std::cerr.rdbuf());

      // output array (# species x # time points)
      size_t numSpecies = simInfo.getModel().getSpeciesCount();
      auto result = py::array_t<unsigned int>(
          {samples, numTimepoints, numSpecies}, // shape
          {sizeof(unsigned int) * numTimepoints * numSpecies,
           sizeof(unsigned int) * numSpecies, sizeof(unsigned int)} // strides
      );
      py::buffer_info buf_result = result.request();
      simInfo.ptrarRawPopulations = static_cast<unsigned int *>(buf_result.ptr);

      // run simulation
      std::thread simThread(&pSSAlibWrapper::runSimulation, this, &simInfo);

      while (bProcessing) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        if (PyErr_CheckSignals() == -1) {
          PY_ERRMSG("simulation cancelled by keyboard interrupt")

          simInfo.bInterruptRequested = true;
        }
      }

      simThread.join();

      // done
      bProcessing = false;

      return result;

    } catch (const std::exception &e) {
      bProcessing = false;
      throw;
    }
  }

public:
  //! Simulation method
  pssalib::PSSA::EMethod method;

  //! Atomic flag to interrupt execution
  std::atomic<bool> bProcessing;

  //! Default constructor
  pSSAlibWrapper() : bProcessing(false), method(pssalib::PSSA::M_PDM) {
    // do nothing
  }

  //! Constructor
  pSSAlibWrapper(pssalib::PSSA::EMethod m) : bProcessing(false), method(m) {
    // do nothing
  }

  //! Copy Constructor
  pSSAlibWrapper(const pSSAlibWrapper &other) : bProcessing(false), method(other.method) {
    // do nothing
  }

  //! Convert to a string
  std::string str() const {
    return (boost::format("method='%1%'") %
            pssalib::PSSA::getMethodName(method))
        .str();
  }

  //! Stringified represantation of the object
  std::string repr() const {
    return (boost::format("pSSAlib(method='%1%')") %
            pssalib::PSSA::getMethodName(method))
        .str();
  }

  //! Sample individual trajectories
  py::array_t<unsigned int> sample_trajectories(
      enum tagTestCase testCase,
      py::array_t<double, py::array::c_style | py::array::forcecast> params,
      double timeEnd, size_t samples, double timeStart, double timeStep) {
    return run(testCase, params, samples, timeEnd, timeStart, timeStep);
  }

  //! Sample species population at a final time
  py::array_t<unsigned int> sample_finalpops(
      enum tagTestCase testCase,
      py::array_t<double, py::array::c_style | py::array::forcecast> params,
      double timeEnd, size_t samples) {
    return run(testCase, params, samples, timeEnd, -1.0, 0.0);
  }

  //! Compute mass-action ODEs for a given model
  py::array_t<double> compute_odes(enum tagTestCase testCase,
      py::array_t<double, py::array::c_style | py::array::forcecast> params,
      py::array_t<double, py::array::c_style | py::array::forcecast> population)
  {
      // check params
      py::buffer_info buf_params = params.request();

      if (buf_params.ndim != 1)
          throw std::runtime_error("Parameters must be a float vector");

      // check population
      py::buffer_info buf_population = population.request();

      if (buf_population.ndim != 1)
          throw std::runtime_error("Populations must be a float vector");

      // input parameters vector
      double *ptr_params = static_cast<double *>(buf_params.ptr);
      size_t num_params = buf_params.size;

      // initalize the model
      pssalib::datamodel::detail::Model model;
      if (!initializeModel(model, testCase, ptr_params,
                            num_params))
          throw std::runtime_error("Failed to intialize test case model");

      // normalize the model
      model.normalize();

      // input population vector
      double *ptr_population = static_cast<double *>(buf_population.ptr);
      size_t num_population = buf_population.size;

      if(num_population != model.getSpeciesCount())
          throw std::length_error("population array size must have same size as number of species in the model");

      // intialize output vector
      size_t numSpecies = model.getSpeciesCount();
      auto result = py::array_t<double>(
          {numSpecies}, // shape
          {sizeof(double)} // strides
      );
      py::buffer_info buf_result = result.request();
      double *ptr_result = static_cast<double *>(buf_result.ptr);

      // compute the ODEs
      computeODEs(model, ptr_population, ptr_result);

      return result;
  }

  //! Compute Jacobian of mass-action ODEs for a given model
  py::array_t<double> compute_jacobian(enum tagTestCase testCase,
      py::array_t<double, py::array::c_style | py::array::forcecast> params,
      py::array_t<double, py::array::c_style | py::array::forcecast> population)
  {
      // check params
      py::buffer_info buf_params = params.request();

      if (buf_params.ndim != 1)
          throw std::runtime_error("Parameters must be a float vector");

      // check population
      py::buffer_info buf_population = population.request();

      if (buf_population.ndim != 1)
          throw std::runtime_error("Populations must be a float vector");

      // input parameters vector
      double *ptr_params = static_cast<double *>(buf_params.ptr);
      size_t num_params = buf_params.size;

      // initalize the model
      pssalib::datamodel::detail::Model model;
      if (!initializeModel(model, testCase, ptr_params,
                            num_params))
          throw std::runtime_error("Failed to intialize test case model");

      // normalize the model
      model.normalize();

      // input population vector
      double *ptr_population = static_cast<double *>(buf_population.ptr);
      size_t num_population = buf_population.size;

      if(num_population != model.getSpeciesCount())
          throw std::length_error("population array size must have same size as number of species in the model");

      // intialize output vector
      size_t numSpecies = model.getSpeciesCount();
      auto result = py::array_t<double>(
          {numSpecies, numSpecies}, // shape
          {numSpecies * sizeof(double), sizeof(double)} // strides
      );
      py::buffer_info buf_result = result.request();
      double *ptr_result = static_cast<double *>(buf_result.ptr);

      // compute the ODEs
      computeJacobian(model, ptr_population, ptr_result);

      return result;
  }

  //! Compute Q-term for Lyapunov Equation based on mass-action ODEs for a given model
  py::array_t<double> compute_lyapunovQ(enum tagTestCase testCase,
      py::array_t<double, py::array::c_style | py::array::forcecast> params,
      py::array_t<double, py::array::c_style | py::array::forcecast> population)
  {
      // check params
      py::buffer_info buf_params = params.request();

      if (buf_params.ndim != 1)
          throw std::runtime_error("Parameters must be a float vector");

      // check population
      py::buffer_info buf_population = population.request();

      if (buf_population.ndim != 1)
          throw std::runtime_error("Populations must be a float vector");

      // input parameters vector
      double *ptr_params = static_cast<double *>(buf_params.ptr);
      size_t num_params = buf_params.size;

      // initalize the model
      pssalib::datamodel::detail::Model model;
      if (!initializeModel(model, testCase, ptr_params,
                            num_params))
          throw std::runtime_error("Failed to intialize test case model");

      // normalize the model
      model.normalize();

      // input population vector
      double *ptr_population = static_cast<double *>(buf_population.ptr);
      size_t num_population = buf_population.size;

      if(num_population != model.getSpeciesCount())
          throw std::length_error("population array size must have same size as number of species in the model");

      // intialize output vector
      size_t numSpecies = model.getSpeciesCount();
      auto result = py::array_t<double>(
          {numSpecies, numSpecies}, // shape
          {numSpecies * sizeof(double), sizeof(double)} // strides
      );
      py::buffer_info buf_result = result.request();
      double *ptr_result = static_cast<double *>(buf_result.ptr);

      // compute the ODEs
      computeQ(model, ptr_population, ptr_result);

      return result;
  }

  //! Print reaction network for a given model
  std::string print_test_case(enum tagTestCase testCase,
      py::array_t<double, py::array::c_style | py::array::forcecast> params, bool odes = false)
  {
      // check params
      py::buffer_info buf_params = params.request();

      if (buf_params.ndim != 1)
          throw std::runtime_error("Parameters must be a float vector");

      // input parameters vector
      double *ptr_params = static_cast<double *>(buf_params.ptr);
      size_t num_params = buf_params.size;

      // initalize the model
      pssalib::datamodel::detail::Model model;
      if (!initializeModel(model, testCase, ptr_params,
                            num_params))
          throw std::runtime_error("Failed to intialize test case model");

      // normalize the model
      model.normalize();

      // print reaction network
      std::stringstream ssTemp;
      printReactionNetwork(model, ssTemp);
      if(odes)
        printODEs(model, ssTemp);

      return ssTemp.str();
  }
};

//! Compute analytical PDF for Homoreaction model
py::array_t<double> homoreactionPDF(
    py::array_t<unsigned int, py::array::c_style | py::array::forcecast> n,
    double k1, double k2, double omega) {
  py::buffer_info buf_n = n.request();

  if (buf_n.ndim != 1)
    throw std::runtime_error("Species population must a vector");

  auto result = py::array_t<double>(buf_n.size);

  py::buffer_info buf_result = result.request();

  unsigned int *ptr_n = static_cast<unsigned int *>(buf_n.ptr);
  double *ptr_result = static_cast<double *>(buf_result.ptr);

  computeHomoreactionPDF(ptr_n, buf_n.shape[0], k1, k2, omega, ptr_result);

  return result;
}

PYBIND11_MODULE(pypssalib, m) {
  m.doc() = R"pbdoc(
      Python bindings for PSSAlib
      ---------------------------

      .. currentmodule:: pypssalib

      .. autosummary::
          :toctree: _generate

          Testcase
          pSSAlib
          homoreaction_pdf
  )pbdoc";

  py::enum_<pSSAlibWrapper::tagTestCase>(m, "Testcase")
      .value("none", pSSAlibWrapper::tcNone)
      .value("clc", pSSAlibWrapper::tcCLC)
      .value("CLC", pSSAlibWrapper::tcCLC)
      .value("CyclicLinearChain", pSSAlibWrapper::tcCLC)
      .value("ca", pSSAlibWrapper::tcCA)
      .value("CA", pSSAlibWrapper::tcCA)
      .value("ColloidalAggregation", pSSAlibWrapper::tcCA)
      .value("homoreaction", pSSAlibWrapper::tcHomo)
      .value("Homoreaction", pSSAlibWrapper::tcHomo)
      .value("sbd", pSSAlibWrapper::tcSBD)
      .value("SBD", pSSAlibWrapper::tcSBD)
      .value("SingleBirthDeath", pSSAlibWrapper::tcSBD)
      .value("tcs", pSSAlibWrapper::tcTCS)
      .value("TCS", pSSAlibWrapper::tcTCS)
      .value("TwoComponentSystem", pSSAlibWrapper::tcTCS)
      .value("ed", pSSAlibWrapper::tcED)
      .value("ED", pSSAlibWrapper::tcED)
      .value("EnzymaticDegration", pSSAlibWrapper::tcED);

  py::enum_<pssalib::PSSA::EMethod>(m, "SSA")
      .value("none", pssalib::PSSA::M_Invalid)
      .value("dm", pssalib::PSSA::M_DM)
      .value("DM", pssalib::PSSA::M_DM)
      .value("DirectMethod", pssalib::PSSA::M_DM)
      .value("pdm", pssalib::PSSA::M_PDM)
      .value("PDM", pssalib::PSSA::M_PDM)
      .value("PartialPropensityDirectMethod", pssalib::PSSA::M_PDM)
      .value("spdm", pssalib::PSSA::M_SPDM)
      .value("SPDM", pssalib::PSSA::M_SPDM)
      .value("SortingPartialPropensityDirectMethod", pssalib::PSSA::M_SPDM)
      .value("pssacr", pssalib::PSSA::M_PSSACR)
      .value("PSSACR", pssalib::PSSA::M_PSSACR)
      .value("CompositionRejectionSampling", pssalib::PSSA::M_PSSACR);

  py::class_<pSSAlibWrapper>(m, "pSSAlib", py::buffer_protocol(), R"pbdoc(
      Interface to pSSAlib simulation engine
      --------------------------------------

      .. currentclass:: pSSAlib

      .. autosummary::
          :toctree: _generate

          method
          sample_testcase_trajectory
          sample_testcase_population
  )pbdoc")
      .def(py::init<>())
      .def(py::init<const pssalib::PSSA::EMethod>())
      .def("__str__", &pSSAlibWrapper::str)
      .def("__repr__", &pSSAlibWrapper::repr)
      .def_readwrite("method", &pSSAlibWrapper::method, R"pbdoc(
        SSA used by the simulation backend.

    )pbdoc")
      .def("sample_testcase_trajectory", &pSSAlibWrapper::sample_trajectories,
           R"pbdoc(
        Generate simulated trajectories of a given test case
        model using an SSA set by `method` attribute.

        Parameters
        ----------

        :param test_case: test case identifier
        :type test_case: Testcase
        :param params: model parameters for the test case
        :type params: ndarray[float]
        :param time_end: final simulation time
        :type time_end: float
        :param time_start: simulation time when output begins
        :type time_start: float
        :param samples: number of trajectories to generate
        :type samples: int
        :param time_step: output time step
        :type time_step: float
        :return: trajectories in an array [samples x time_points x species]
        :rtype: ndarray[float]
    )pbdoc",
           py::arg("test_case"), py::arg("params"), py::arg("time_end"),
           py::arg("samples") = 1, py::arg("time_start") = 0,
           py::arg("time_step") = 1e-1)
      .def("sample_testcase_population", &pSSAlibWrapper::sample_finalpops,
           R"pbdoc(
        Generate populations of a given test case model at a certain time point
        using an SSA set by `method` attribute.

        Parameters
        ----------

        :param test_case: test case identifier
        :type test_case: Testcase
        :param params: model parameters for the test case
        :type params: ndarray[float]
        :param time_end: final simulation time
        :type time_end: float
        :param samples: number of trajectories to generate
        :type samples: int
        :return: species populations in an array [samples x 1 x species]
        :rtype: ndarray[float]
    )pbdoc",
           py::arg("test_case"), py::arg("params"), py::arg("time_end"),
           py::arg("samples") = 1)
      .def("odes", &pSSAlibWrapper::compute_odes, R"pbdoc(
        Compute mass-action ODEs for a given model

        Parameters
        ----------

        :param test_case: test case identifier
        :type test_case: Testcase
        :param params: model parameters for the test case
        :type params: ndarray[float]
        :param population: species population
        :type population: ndarray[float]
        :return: mass-action ODEs value for respective species population
        :rtype: ndarray[float]
    )pbdoc",
        py::arg("test_case"), py::arg("params"), py::arg("population"))
      .def("jacobian", &pSSAlibWrapper::compute_jacobian, R"pbdoc(
        Compute Jacobian of mass-action ODEs for a given model

        Parameters
        ----------

        :param test_case: test case identifier
        :type test_case: Testcase
        :param params: model parameters for the test case
        :type params: ndarray[float]
        :param population: species population
        :type population: ndarray[float]
        :return: Jacobian of mass-action ODEs value for respective species population
        :rtype: ndarray[float]
    )pbdoc",
        py::arg("test_case"), py::arg("params"), py::arg("population"))
      .def("lyapunovQ", &pSSAlibWrapper::compute_lyapunovQ, R"pbdoc(
        Compute Q-term for Lyapunov Equation based on mass-action ODEs for a given model

        Parameters
        ----------

        :param test_case: test case identifier
        :type test_case: Testcase
        :param params: model parameters for the test case
        :type params: ndarray[float]
        :param population: species population
        :type population: ndarray[float]
        :return: Q-term for Lyapunov Equation value for respective species population
        :rtype: ndarray[float]
    )pbdoc",
        py::arg("test_case"), py::arg("params"), py::arg("population"))
      .def("print_test_case", &pSSAlibWrapper::print_test_case, R"pbdoc(
        Get reaction network for a given model as a string

        Parameters
        ----------

        :param test_case: test case identifier
        :type test_case: Testcase
        :param params: model parameters for the test case
        :type params: ndarray[float]
        :param odes: if True, prints a system of ordinary differential equations
                     for temporal evolution of species concentrations in the model
        :type params: boolean
        :return: string containing reaction network for a given model
        :rtype: str
    )pbdoc",
        py::arg("test_case"), py::arg("params"), py::arg("odes") = false)
      .def(py::pickle(
        [](const pSSAlibWrapper &instance) { // __getstate__
            /* Return a tuple that fully encodes the state of the object */
            return py::make_tuple(instance.method);
        },
        [](py::tuple t) { // __setstate__
            if (t.size() != 1)
                throw std::runtime_error("Invalid state!");

            /* Create a new C++ instance */
            return pSSAlibWrapper(t[0].cast<pssalib::PSSA::EMethod>());
        }
    ));

  m.def("homoreaction_pdf", homoreactionPDF, R"pbdoc(
        Compute analytic PDF for Homoreaction model

        Parameters
        ----------

        :param n: species populations.
        :type n: ndarray[int]
        :param k1: first reaction rate
        :type k1: float
        :param k2: second reaction rate
        :type k2: float
        :param omega: subvolume size
        :type omega: float
        :return: analytical PDF values for respective species populations
        :rtype: ndarray[float]
    )pbdoc",
        py::arg("n"), py::arg("k1"), py::arg("k2"), py::arg("omega"));

#ifdef PYPSSALIB_VERSION
    m.attr("__version__") = MACRO_STRINGIFY(PYPSSALIB_VERSION);
#else
    m.attr("__version__") = "dev";
#endif
}
