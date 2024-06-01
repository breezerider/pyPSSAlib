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

  void runSimulation(pssalib::datamodel::SimulationInfo *ptrSimInfo)
  {
    pssalib::PSSA simEngine;

    simEngine.setMethod(pssalib::PSSA::M_PDM);

    simEngine.run(ptrSimInfo);

    bProcessing = false;
  }

public:


    //! Atomic flag to interrupt execution
    std::atomic<bool> bProcessing;

    //! Default constructor
    pSSAlibWrapper():
        bProcessing(false)
    {
        // do nothing
    }

    py::array_t<unsigned int> run(enum tagTestCase testCase,
                                  py::array_t<double, py::array::c_style | py::array::forcecast> params,
                                  double timeStart, double timeEnd, double timeStep)
    {
        if(bProcessing)
            throw std::runtime_error("Another simualtion is running already");
        bProcessing = true;

        size_t numTimepoints = pssalib::timing::getNumTimePoints(timeStart, timeEnd, timeStep);

        if(numTimepoints <= 0)
            throw std::runtime_error("Time interval must include art least one time point");

        py::buffer_info buf_params = params.request();

        if (buf_params.ndim != 1)
            throw std::runtime_error("Number of dimensions must be one");

        // input parameters vector
        double *ptr_params = static_cast<double *>(buf_params.ptr);
        size_t  num_params = buf_params.size;

        // intialize simulation info
        pssalib::datamodel::SimulationInfo simInfo;

        if(!initializeModel(simInfo.getModel(), testCase, ptr_params, num_params))
            throw std::runtime_error("Failed to intialize test case model");

        // setup the simulation params
        simInfo.unSamplesTotal = 1;
        simInfo.dTimeStart = timeStart;
        simInfo.dTimeEnd = timeEnd;
        simInfo.dTimeStep = timeStep;

        // setup output
        simInfo.unOutputFlags = pssalib::datamodel::SimulationInfo::ofNone
          | pssalib::datamodel::SimulationInfo::ofRawTrajectory
          | pssalib::datamodel::SimulationInfo::ofTrajectory
    //       | pssalib::datamodel::SimulationInfo::ofStatus
          | pssalib::datamodel::SimulationInfo::ofLog
    //       | pssalib::datamodel::SimulationInfo::ofTrace
    //       | pssalib::datamodel::SimulationInfo::ofInfo
    //       | pssalib::datamodel::SimulationInfo::ofWarning
          | pssalib::datamodel::SimulationInfo::ofError
    //       | pssalib::datamodel::SimulationInfo::eofModuleGrouping
    //       | pssalib::datamodel::SimulationInfo::eofModuleSampling
    //       | pssalib::datamodel::SimulationInfo::eofModuleUpdate
          ;
        simInfo.resetOutput();

        simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());

        // output array (# species x # time points)
        auto result = py::array_t<unsigned int>(simInfo.getModel().getSpeciesCount()*numTimepoints);
        py::buffer_info buf_result = result.request();
        simInfo.ptrarRawTrajectory = static_cast<unsigned int *>(buf_result.ptr);

        // run simulation
        std::thread simThread(&pSSAlibWrapper::runSimulation, this, &simInfo);

        while(bProcessing)
        {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            if(PyErr_CheckSignals() == -1)
            {
                PY_ERRMSG("simulation cancelled by keyboard interrupt")

                simInfo.bInterruptRequested = true;
            }
        }

        simThread.join();

        // done
        bProcessing = false;

        return result;
    }
};

PYBIND11_MODULE(pypssalib, m) {
    m.doc() = R"pbdoc(
        Python bindings for PSSAlib
        ---------------------------

        .. currentmodule:: pypssalib

        .. autosummary::
           :toctree: _generate

           run_testcase
           num_timepoints
           num_species
           test_case
           time_start
           time_step
           time_end
    )pbdoc";

    py::enum_<pSSAlibWrapper::tagTestCase>(m, "Testcase")
        .value("none", pSSAlibWrapper::tcNone)
        .value("sbd", pSSAlibWrapper::tcSBD)
        .value("SBD", pSSAlibWrapper::tcSBD)
        .value("SingleBirthDeath", pSSAlibWrapper::tcSBD);

    py::class_<pSSAlibWrapper>(m, "pSSAlib", py::buffer_protocol())
        .def(py::init<>())
        //.def("__str__", &pSSAlibWrapper::str)
        //.def("__repr__", &pSSAlibWrapper::repr)
        .def("run_testcase", &pSSAlibWrapper::run, "Generate simulated trajectory of a given test case model using an SSA",
             py::arg("test_case") = pSSAlibWrapper::tcNone, py::arg("params") = py::array_t<double>(),
             py::arg("time_start") = 0, py::arg("time_end") = 1.0, py::arg("time_step") = 1e-1);
        // .def_readwrite("test_case", &pSSAlibWrapper::testCase)
        // .def_readonly("num_timepoints", &pSSAlibWrapper::szTimePointsCount)
        // .def_readonly("num_species", &pSSAlibWrapper::szSpeciesCount)
        // .def_readwrite("time_start", &pSSAlibWrapper::dTimeStart)
        // .def_readwrite("time_step", &pSSAlibWrapper::dTimeStep)
        // .def_readwrite("time_end", &pSSAlibWrapper::dTimeEnd)
        // .def_buffer([](pSSAlibWrapper &w) -> py::buffer_info {
        //         return py::buffer_info(
        //             w.arRawTrajectories.get(),                     /* Pointer to buffer */
        //             sizeof(unsigned int),                          /* Size of one scalar */
        //             py::format_descriptor<unsigned int>::format(), /* Python struct-style format descriptor */
        //             2,                                             /* Number of dimensions */
        //             { w.szTimePointsCount, w.szSpeciesCount },     /* Buffer dimensions */
        //             { sizeof(unsigned int) * w.szSpeciesCount,     /* Strides (in bytes) for each index */
        //             sizeof(unsigned int) }
        //         );
        //     });



#ifdef PYPSSALIB_VERSION
    m.attr("__version__") = MACRO_STRINGIFY(PYPSSALIB_VERSION);
#else
    m.attr("__version__") = "dev";
#endif
}
