// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_tools.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
#include <boost/multi_array.hpp>
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
// Problem parameters.
  static double mu_stem = 1.0e-2;
  static double mu_leaf = 1.0e-2;
  static const double PI = acos(-1);
  static double beta_stem = 1.0e1;
  static double beta_leaf = 1.0e1;
  static double J_maximum=1;
  static double J_minimum=1;
  static double t_set=0;
  static double kappa_s = 5.0e11;
  static double pulse_period = 1.89;                //period of a pulse in seconds
  static double delta_time = 0.54;           //contraction time
  static double expand_time = 0.71;           //expansion time
  static double rest_time = 0.64;           //rest time

  static double phase_difference = PI/4.0;           //phase difference between polyps in radians
  static double disp_x1 = 0.0;                      // displacement of coral 1 in x-direction
  static double disp_y1 = 0.0;                      // displacement of coral 1 in y-direction
  static double disp_z1 = 0.0;                      // displacement of coral 1 in z-direction
  static double rc1 = 0.0009;                        //radius of the stalk in meters
  static double xc1 = 0.00;                       //x-coordinate of first polyp center
  static double zc1 = 0.0;                         //y-coordinate of first polyp center
  static double phi1 = PI/16.0;            // angles for test motion
  static double phi2 = PI/8.0;
  static double phi3 = PI/12.0;
  static double init_time1 = pulse_period/10.0;      // Time to move tentacles to initial position
  static double betao = 4.0;
  static double betam = 0.2;
  static double to = 0.4;
  static double gamma = 2.0;
  static double tentacle_length = 0.0041;
  static double prefactor = 1.9;

// The following function outputs data

void
output_data(
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
    Mesh& mesh,
    EquationSystems* equation_systems,
    const int iteration_num,
    const double loop_time,
    const string& data_dump_dirname);



// The following two functions move the corals //////////////////////////////////////
void
target1_force_function(
    VectorValue<double>& F,
	const TensorValue<double>& /*FF*/,
	const libMesh::Point& X,
	const libMesh::Point& s,
	Elem* const /*elem*/,
	const std::vector<const std::vector<double>*>& /*var_data*/,
	const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
	double time,
    void* /*ctx*/)
{
  libMesh::Point s_dump;
  double kappa;

  // s_dump is the target point position 
  // s is the reference configuration. 
  // The current position is X

  // parameters for the coral///////////////////////////////////////////

  // parameters used for calculations
  double rc;             // distance from the coral center to lagrangian point
  double radb;             // distance from tentacle base to Lagrangian point
  double xb;             // x- position of lagrangian point relative to coral center
  double zb;             // y- position of lagrangian point relative to coral center
  double thetab;         // angle of lagrangian point relative to positive x-axis
  double pulse_no;      // pulse number of the coral 
  double curr_time;     // time since beginning of current pulse cycle
  double y0;             // height of beginning polynomial fit
  double y1;             // height of ending polynomial fit
  double r0;             // radial position of beginning polynomial fit
  double r1;             // radial position of ending polynomial fit
  double yfit;             // height of current fit
  double rfit;              // radial position of current fit

  double beta = 0.5;
  double zl=0.0;
  double phi = 0.0;
  double ro=0.0;
  double thetaj = 0.0;
  double ds1 = 1.0/128.0;
  
  int s1=0;
  int N_tent=0;

  kappa=kappa_s;

  // Determine radius and angle of Lagrangian point relative to coral center
  
     s_dump = s;                // initialize target point position
     xb = s(0)+disp_x1-xc1;     // determine x-position relative to center
     zb = s(2)+disp_z1-zc1;     // determine y-position relative to center
     rc = sqrt(xb*xb+zb*zb);    // determine distance from center
     thetab=atan2(zb,xb);       // determine angle from x-axis
     radb = rc-rc1;            // determine distance from tentacle base 

        if((time<init_time1) && (rc>=rc1))     //check if in inital moving stage
        {
            // move into initial position
	    phi = 0;
	    beta = betao+(betam-betao)*pow((0.0/delta_time),gamma);

	    thetaj = -prefactor*(1.0-exp(-(s1)*ds1/beta));
            N_tent = floor((radb/tentacle_length)/ds1);
            zl=0.0;
            ro=0.0;
              
            for (s1 = 1; s1 <N_tent; ++s1){
              thetaj = -prefactor*(1-exp(-(s1)*ds1/beta));
              zl += ds1*sin(thetaj);
              ro += ds1*cos(thetaj);
            }
	    
            y0 = s(1)+disp_y1;                 // height at beginning
            y1 = s(1)+disp_y1 - tentacle_length*zl;   // height at end
            r0 = radb;                        // distance from tentacle base at beginning
            r1 = tentacle_length*ro;              // distance from tentacle base at end
            yfit = y0*(1-(time/init_time1))+y1*(time/init_time1);   // interpolated current height
            rfit = r0*(1-(time/init_time1))+r1*(time/init_time1);   // interpolated current distance from base in xy plane
            s_dump(0)=disp_x1+(rc1+rfit)*cos(thetab);          // target x-position    
            s_dump(2)=disp_z1+(rc1+rfit)*sin(thetab);          // target y - xposition
            s_dump(1)=s(1)+disp_y1+yfit;                            // target z-position
            F = kappa*(s_dump-X);                                   // target force
        }
        if((time>=init_time1) && (rc>=rc1)) // we are in the pulsing cycle
        {
            pulse_no = floor((time-init_time1)/pulse_period);         //figure out what pulse number we are on
            curr_time = time-(init_time1+pulse_no*pulse_period);    //figure out time since beginning of pulse
              if (curr_time<delta_time){  //contract tentacle
                phi = 0;
                beta = betao+(betam-betao)*pow((curr_time/delta_time),gamma);
              }
              else if ((curr_time>=delta_time) && (curr_time<= delta_time+expand_time)){ //expand tentacle
                phi = 0; 
                beta = betam+(betao-betam)*pow(((curr_time-delta_time)/(expand_time)),gamma);
              }
			  else {
				phi = 0;
				beta = betao;
			  }
              thetaj = -1.55*(1.0-exp(-(s1)*ds1/beta));
              N_tent = floor((radb/tentacle_length)/ds1);
              zl=0.0;
              ro=0.0;
              
              for (s1 = 1; s1 <N_tent; ++s1){
                thetaj = -1.55*(1-exp(-(s1)*ds1/beta));
                zl += ds1*sin(thetaj);
                ro += ds1*cos(thetaj);
              }

                yfit = -tentacle_length*zl;   // interpolated current height
                rfit = tentacle_length*ro;   // interpolated current distance from base in xy plane
                s_dump(0)=disp_x1+(rc1+rfit)*cos(thetab);          // target x-position    
                s_dump(2)=disp_z1+(rc1+rfit)*sin(thetab);          // target y - xposition
                s_dump(1)=s(1)+disp_y1+yfit;                            // target z-position
                F = kappa*(s_dump-X);                                   // target force
          }

        if(rc<rc1)      //check if on stem
        {
            s_dump(0)=s(0)+disp_x1;
            s_dump(1)=s(1)+disp_y1;
            s_dump(2)=s(2)+disp_z1;
            F = kappa*(s_dump-X);
        } 

    return;
} 

// Stress tensor functions.
void
PK1_dev_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
	const libMesh::Point& /*X*/,
	const libMesh::Point& /*s*/,
	Elem* const /*elem*/,
	const std::vector<const std::vector<double>*>& /*var_data*/,
	const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
	double /*time*/,
    void* /*ctx*/)
{
double mu;
 mu= mu_stem;
  PP = (mu)*FF;
    return;
}// PK1_dev_stress_function

void
PK1_dil_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
	const libMesh::Point& /*X*/,
	const libMesh::Point& /*s*/,
	Elem* const /*elem*/,
	const std::vector<const std::vector<double>*>& /*var_data*/,
	const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
	double /*time*/,
    void* /*ctx*/)
{
double mu;
double beta;

mu= mu_stem;
beta=beta_stem;
 
PP = (-1.0*(mu)+beta*log(FF.det()))*tensor_inverse_transpose(FF,NDIM);

    return;
}// PK1_dil_stress_function
}

using namespace ModelData;

// Function prototypes
static ofstream forcex_stream, forcey_stream, forcez_stream, velx_stream, vely_stream, velz_stream, workx_stream, worky_stream, workz_stream, U_L1_norm_stream, U_L2_norm_stream, U_max_norm_stream;
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(
    int argc,
    char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    {// cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string mesh1_exodus_filename = app_initializer->getExodusIIFilename("mesh1");

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        const string restart_read_dirname = app_initializer->getRestartReadDirectory();
        const int restart_restore_num = app_initializer->getRestartRestoreNumber();
        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Create a simple FE mesh with Dirichlet boundary conditions.
        //
        // Note that boundary condition data must be registered with each FE
        // system before calling IBFEMethod::initializeFEData().

        Mesh mesh1(init.comm(), NDIM);

        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC")*dx;
	const string mesh_filename = input_db->getStringWithDefault("MESH_FILENAME","Leaf.e");
     string elem_type = input_db->getString("ELEM_TYPE");

    mesh1.read("coral_tri_0p1.e");
    mesh1.prepare_for_use();
    std::cout<< "WHAT DOES MY MESH LOOK LIKE"  << std::endl;
    std::cout<< mesh1<< std::endl;
    vector<Mesh*> meshes(1);

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        	mu_leaf = input_db->getDouble("MU_LEAF");
	mu_stem = input_db->getDouble("MU_STEM");
	kappa_s = input_db->getDouble("KAPPA_S");
	beta_stem = input_db->getDouble("BETA_STEM");
	beta_leaf = input_db->getDouble("BETA_LEAF");

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator", app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n" <<
                       "Valid options are: COLLOCATED, STAGGERED");
        }
		Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh1,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true,
                           restart_read_dirname,
                           restart_restore_num);
        Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
            "IBHierarchyIntegrator", app_initializer->getComponentDatabase("IBHierarchyIntegrator"), ib_method_ops, navier_stokes_integrator);
Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = main_db->getStringWithDefault("adv_diff_solver_type", "PREDICTOR_CORRECTOR");
        if (adv_diff_solver_type == "PREDICTOR_CORRECTOR")
        {
            Pointer<AdvectorExplicitPredictorPatchOps> predictor = new AdvectorExplicitPredictorPatchOps("AdvectorExplicitPredictorPatchOps", app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
            adv_diff_integrator = new AdvDiffPredictorCorrectorHierarchyIntegrator("AdvDiffPredictorCorrectorHierarchyIntegrator", app_initializer->getComponentDatabase("AdvDiffPredictorCorrectorHierarchyIntegrator"), predictor);
        }
        else if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator("AdvDiffSemiImplicitHierarchyIntegrator", app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n" <<
                       "Valid options are: PREDICTOR_CORRECTOR, SEMI_IMPLICIT");
        }
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);


        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
            "PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

    const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;

    Pointer<CellVariable<NDIM,double> > T_var = new CellVariable<NDIM,double>("T");
    adv_diff_integrator->registerTransportedQuantity(T_var);
    adv_diff_integrator->setDiffusionCoefficient(T_var, input_db->getDouble("KAPPA"));
    adv_diff_integrator->setInitialConditions(T_var, new muParserCartGridFunction("T_init", app_initializer->getComponentDatabase("TemperatureInitialConditions"), grid_geometry));
    RobinBcCoefStrategy<NDIM>* T_bc_coef = NULL;
    if (!periodic_domain)
      {
        T_bc_coef = new muParserRobinBcCoefs("T_bc_coef", app_initializer->getComponentDatabase("TemperatureBcCoefs"), grid_geometry);
        adv_diff_integrator->setPhysicalBcCoef(T_var, T_bc_coef);
      }
    adv_diff_integrator->setAdvectionVelocity(T_var, navier_stokes_integrator->getAdvectionVelocityVariable());

    // Configure the IBFE solver.
	ib_method_ops->initializeFEEquationSystems();
    	
	FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
    IBFEMethod::LagBodyForceFcnData target1_force_data(target1_force_function);
    IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function);
    IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);

    PK1_dev_stress_data.quad_order = Utility::string_to_enum<libMeshEnums::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER","SEVENTH"));
    PK1_dil_stress_data.quad_order = Utility::string_to_enum<libMeshEnums::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER","SEVENTH"));

    ib_method_ops->registerLagBodyForceFunction(target1_force_data, 0);
    ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data, 0);
    ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data, 0);

	EquationSystems* mesh1_equation_systems = ib_method_ops->getFEDataManager(0)->getEquationSystems();

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction("u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction("p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction("f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        
	    std::unique_ptr<ExodusII_IO> mesh1_exodus_io(uses_exodus ? new ExodusII_IO(mesh1) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus) {
                mesh1_exodus_io->write_timestep(mesh1_exodus_filename, *mesh1_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
            }
        }
        if (SAMRAI_MPI::getRank() == 0)
        {
            forcex_stream.open("F_x.curve", ios_base::out | ios_base::trunc);
            forcey_stream.open("F_y.curve", ios_base::out | ios_base::trunc);
            forcez_stream.open("F_z.curve", ios_base::out | ios_base::trunc);
	        velx_stream.open("U_x.curve", ios_base::out | ios_base::trunc);
            vely_stream.open("U_y.curve", ios_base::out | ios_base::trunc);
            velz_stream.open("U_z.curve", ios_base::out | ios_base::trunc);
	        workx_stream.open("W_x.curve", ios_base::out | ios_base::trunc);
            worky_stream.open("W_y.curve", ios_base::out | ios_base::trunc);
            workz_stream.open("W_z.curve", ios_base::out | ios_base::trunc);
	        U_L1_norm_stream.open("U_L1.curve", ios_base::out | ios_base::trunc);
	        U_L2_norm_stream.open("U_L2.curve", ios_base::out | ios_base::trunc);
	        U_max_norm_stream.open("U_max.curve", ios_base::out | ios_base::trunc);
        }
	    static const double PI = acos(-1);
        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
	    double max_J = 1;
	    double min_J = 1;

        while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();
	    
            pout <<                                                    "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
	    
            time_integrator->advanceHierarchy(dt);
            max_J=SAMRAI_MPI::maxReduction(J_maximum);
            min_J=SAMRAI_MPI::minReduction(J_minimum);
            pout << "Max Jacobian " <<  max_J << "\n";
            pout << "Min Jacobian " <<  min_J << "\n";
            
            loop_time += dt;

            pout <<                                                    "\n";
            pout << "At end       of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout <<                                                    "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num%viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus){
                    mesh1_exodus_io->write_timestep(mesh1_exodus_filename, *mesh1_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num%restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num%timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }


 if (dump_postproc_data && (iteration_num%postproc_data_dump_interval == 0 || last_step))
            {
                pout << "\nWriting state data...\n\n";
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            forcex_stream.close();
            forcey_stream.close();
            forcez_stream.close();
            velx_stream.close();
            vely_stream.close();
            velz_stream.close();
            workx_stream.close();
            worky_stream.close();
            workz_stream.close();
            U_L1_norm_stream.close();
            U_L2_norm_stream.close();
            U_max_norm_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
}// main

void postprocess_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                      Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int /*iteration_num*/,
                      const double loop_time,
                      const string& /*data_dump_dirname*/)
{
    const unsigned int dim = mesh.mesh_dimension();
    {
        double F_integral[NDIM];
	    double U_integral[NDIM];
	    double W_integral[NDIM];

        for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d) W_integral[d] = 0.0;
        System& F_system = equation_systems->get_system<System>(IBFEMethod::FORCE_SYSTEM_NAME);
	    System& U_system = equation_systems->get_system<System>(IBFEMethod::VELOCITY_SYSTEM_NAME);
        NumericVector<double>* F_vec = F_system.solution.get();
        NumericVector<double>* U_vec = U_system.solution.get();
        NumericVector<double>* F_ghost_vec = F_system.current_local_solution.get();
        NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
        F_vec->localize(*F_ghost_vec);
        U_vec->localize(*U_ghost_vec);
        DofMap& F_dof_map = F_system.get_dof_map();
        DofMap& U_dof_map = U_system.get_dof_map();
        std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);
        std::vector<std::vector<unsigned int> > U_dof_indices(NDIM);
        unique_ptr<FEBase> fe(FEBase::build(dim, F_dof_map.variable_type(0)));
        unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, FIFTH);
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<std::vector<double> >& phi = fe->get_phi();
        const std::vector<double>& JxW = fe->get_JxW();
        boost::multi_array<double, 2> F_node;
        boost::multi_array<double, 2> U_node;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_dof_map.dof_indices(elem, F_dof_indices[d], d);
                U_dof_map.dof_indices(elem, U_dof_indices[d], d);
            }
            const int n_qp = qrule->n_points();
            const int n_basis = F_dof_indices[0].size();
            get_values_for_interpolation(F_node, *F_ghost_vec, F_dof_indices);
            get_values_for_interpolation(U_node, *U_ghost_vec, U_dof_indices);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                for (int k = 0; k < n_basis; ++k)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F_integral[d] += F_node[k][d] * phi[k][qp] * JxW[qp];
			            U_integral[d] += U_node[k][d] * phi[k][qp] * JxW[qp];
			            W_integral[d] += U_node[k][d]*F_node[k][d] * phi[k][qp] * JxW[qp];
                    }
                }
            }
        }
        SAMRAI_MPI::sumReduction(F_integral, NDIM);
        SAMRAI_MPI::sumReduction(U_integral, NDIM);
        SAMRAI_MPI::sumReduction(W_integral, NDIM);

        if (SAMRAI_MPI::getRank() == 0)
        {
            forcex_stream << loop_time << " " << F_integral[0] << endl;
            forcey_stream << loop_time << " " << F_integral[1] << endl;
            forcez_stream << loop_time << " " << F_integral[2] << endl;
            velx_stream << loop_time << " " << U_integral[0] << endl;
            vely_stream << loop_time << " " << U_integral[1] << endl;
            velz_stream << loop_time << " " << U_integral[2] << endl;
            workx_stream << loop_time << " " << W_integral[0] << endl;
            worky_stream << loop_time << " " << W_integral[1] << endl;
            workz_stream << loop_time << " " << W_integral[2] << endl;
        }
    }

    {
        double U_L1_norm = 0.0, U_L2_norm = 0.0, U_max_norm = 0.0;
        System& U_system = equation_systems->get_system<System>(IBFEMethod::VELOCITY_SYSTEM_NAME);
        NumericVector<double>* U_vec = U_system.solution.get();
        NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
        U_vec->localize(*U_ghost_vec);
        DofMap& U_dof_map = U_system.get_dof_map();
        std::vector<std::vector<unsigned int> > U_dof_indices(NDIM);
        unique_ptr<FEBase> fe(FEBase::build(dim, U_dof_map.variable_type(0)));
        unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, FIFTH);
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<std::vector<double> >& phi = fe->get_phi();
        const std::vector<double>& JxW = fe->get_JxW();
        VectorValue<double> U_qp;
        boost::multi_array<double, 2> U_node;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_dof_map.dof_indices(elem, U_dof_indices[d], d);
            }
            const int n_qp = qrule->n_points();
            get_values_for_interpolation(U_node, *U_ghost_vec, U_dof_indices);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                interpolate(U_qp, qp, U_node, phi);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    U_L1_norm += std::abs(U_qp(d)) * JxW[qp];
                    U_L2_norm += U_qp(d) * U_qp(d) * JxW[qp];
                    U_max_norm = std::max(U_max_norm, std::abs(U_qp(d)));
                }
            }
        }
        SAMRAI_MPI::sumReduction(&U_L1_norm, 1);
        SAMRAI_MPI::sumReduction(&U_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&U_max_norm, 1);
        U_L2_norm = sqrt(U_L2_norm);
        if (SAMRAI_MPI::getRank() == 0)
        {
            U_L1_norm_stream << loop_time << " " << U_L1_norm << endl;
            U_L2_norm_stream << loop_time << " " << U_L2_norm << endl;
            U_max_norm_stream << loop_time << " " << U_max_norm << endl;
        }
    }
    return;
} // postprocess_datai
void
output_data(
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
    Mesh& mesh,
    EquationSystems* equation_systems,
    const int iteration_num,
    const double loop_time,
    const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(), navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(), navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

   // Write Lagrangian data.
   file_name = data_dump_dirname + "/" + "fe_mesh.";
   sprintf(temp_buf, "%05d", iteration_num);
   file_name += temp_buf;
   file_name += ".xda";
   mesh.write(file_name);
   file_name = data_dump_dirname + "/" + "fe_equation_systems.";
   sprintf(temp_buf, "%05d", iteration_num);
   file_name += temp_buf;
   equation_systems->write(file_name, (EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA));
   return;
}// output_data
   
