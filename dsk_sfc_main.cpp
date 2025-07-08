#include "dsk_sfc_include.h"
using namespace std;

int main(int argc, char *argv[])
{
  //To measure execution time
  using namespace std::chrono;
  auto start = high_resolution_clock::now();

  //Command line parser
  CLI::App app{"DSK_SFC, the Dystopian Schumpeter meeting Keynes Stock Flow Consistent model"};

  //Add command line input: json input file
  string inputstring = "default.json";
  app.add_option("inputfile", inputstring, "A path to the json input file")
    ->required()
    ->check(CLI::ExistingFile);

  //Add command line input: run name
  string str_runname;
  app.add_option("-r,--run", str_runname, "A name for the run")
    ->default_val("test");

  //Add command line input: seed
  int exseed{1};
  app.add_option("-s", exseed, "A seed (positive integer)")
    ->default_val(1)
    ->check(CLI::PositiveNumber);

  //Add command line input: full output dummy
  int reducedout{0};
  app.add_option("-f,--reducedoutput", reducedout, "If set to 1, full output will be saved");

  //Add command line input: error printing dummy
  int cerr_enabled{0};
  app.add_option("-c,--cerr", cerr_enabled, "If set to 1, print error messages to the console");

  //Add command line input: verbose dummy
  int verbose{0};
  app.add_option("-v,--verbose", verbose, "If set to 1, print simulation progress");

  //Parse command line input
  CLI11_PARSE(app, argc, argv);
  //If run name contains spaces, remove them
  auto noSpace = std::remove(str_runname.begin(), str_runname.end(), ' ');
  str_runname.erase(noSpace, str_runname.end());

  //If dummy for error printing to console is set to 0, suppress error messages
  if(cerr_enabled==0){std::cerr.setstate(std::ios_base::failbit);}

  //Initialise variables needed to process console inputs
  //char const* seednumber = to_string(exseed).c_str();
  std::string seedstring = std::to_string(exseed);
  const char* seednumber = seedstring.c_str();
  char* exec_dir;
  char* runname = &str_runname[0];
  reducedoutput=reducedout;

  //Path to executable as char; needed to create output folder
  exec_dir = argv[0];

  //Seed needs to be converted to negative integer to work with functions generating random draws
  exseed = -exseed;

  //Read JSON file and convert it to a document from which to read parameter, initial and flag values
  using namespace rapidjson;
  std::ifstream ifs { inputstring };
  //Exit if unable to open input file
  if ( !ifs.is_open() )
  {
    cout << "Could not open input file for reading!"<< endl;
    return EXIT_FAILURE;
  }
  //Parse
  IStreamWrapper isw { ifs };
  Document inputs {};
  inputs.ParseStream( isw );
  //Exit if encounter error during parsing
  if (inputs.HasParseError())
  {
    cout << "Parsing of input file failed!"<< endl;
    return EXIT_FAILURE;
  }
  if(verbose){cout << "Finished parsing input file" << endl;}

  //Use the parsed input file to set model parameters, initial values, flags, etc.
  SETPARAMS(inputs);
  if(verbose){cout << "Exiting function SETPARAMS" << endl;}

  //If they do not exist yet, create folders to hold simulation output and error logs
  FOLDERS(exec_dir, runname);
  if(verbose){cout << "Finished creating output folders" << endl;}

  char pathname[32];
  char output_folder[80];
  strcpy(output_folder, "output/");
  char* run_output_folder = strcat(output_folder,runname);
  char* filepath=strcpy(pathname,output_folder);
  //Generate name for the error file
  strcpy(errorfilename,"output/");
  char* name_error=strcat(errorfilename,runname);
  name_error=strcat(errorfilename,"/errors/Errors_");
  name_error=strcat(errorfilename,runname);
  name_error=strcat(errorfilename,"_");
  name_error=strcat(errorfilename,seednumber);
  strcat(errorfilename,".txt");

  //Create name for output file
  GENFILEYMC(filepath,"/ymc", runname, seednumber);

  if(verbose){cout << "Finished creating output file names" << endl;}

  //Resize variable arrays to the needed dimensions based on the numbers of agents and periods
  RESIZE();
  if(verbose){cout << "Exiting function RESIZE" << endl;}

  //Initialise endogenous model variables
	INITIALIZE(exseed);
  if(verbose){cout << "Exiting function INITIALIZE; Entering simulation loop" << endl;}

  INITOUTPUTCOLLECTION();

  //This only works on Linux systems (not Windows or Mac)
  //Mostly needed for large batch runs on HPC to avoid getting stuck
  #ifdef __linux__
    //Set the run to time out after T*2 seconds
    signal(SIGALRM, catchAlarm);
    alarm(T*2);
  #endif

  // Open file for printing Errors
  Errors.open(errorfilename,ios::app);

  /* // To measure initialisation time
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);
  cout << "Time taken by simulation " << seednumber << " initialisation: " << duration.count() << " seconds" << endl; */

  //Enter loop over simulation periods - AA
  for (t=1; t<=T; t++)
  {
    /* //To measure execution time of timestep
    if (t==2|t==400|t==401){
      stop = high_resolution_clock::now();
      duration = duration_cast<seconds>(stop - start);
      cout << "Time taken until after simulation " << seednumber << " period " << t-1 << ": " << duration.count() << " seconds" << endl;
    } */

    if(t==t_regime_shifts)
    {
      REGIMESHIFTS();
      if(verbose){cout << "Exiting function REGIMESHIFTS in period " << t << endl;}
    }

    if(verbose){cout << "Entering simulation period " << t << endl;}
    //Resets endogenous variables as needed to begin a new simulation period
    SETVARS();
    if(verbose){cout << "Exiting function SETVARS in period " << t << endl;}
    //Interest on bank deposits is paid
    DEPOSITINTEREST();
    if(verbose){cout << "Exiting function DEPOSITINTEREST in period " << t << endl;}
    //only enter here every freqclim periods

    if(t%freqclim==0)
    {
      //Update climate policy (only carbon tax for now)
      CLIMATE_POLICY();
      if(verbose){cout << "Exiting function CLIMATE_POLICY in period " << t << endl;}
    }

    //Deliver machines ordered in t-1, calculate unit cost for C and K-Firms, update prices
    MACH();
    if(verbose){cout << "Exiting function MACH in period " << t << endl;}

    //Banks determine maximum amount they are willing to lend
    TOTCREDIT();
    if(verbose){cout << "Exiting function TOTCREDIT in period " << t << endl;}

    //Banks determine the loan interest rate to be charged to each C-Firm customer
    LOANRATES();
    if(verbose){cout << "Exiting function LOANRATES in period " << t << endl;}

    //K-Firms send out brochures to attract customers
    BROCHURE();
    if(verbose){cout << "Exiting function BROCHURE in period " << t << endl;}

    //C-Firms set expected demand & desired production and determine desired investment
    INVEST();
    if(verbose){cout << "Exiting function INVEST in period " << t << endl;}

    ALLOCATECREDIT();
    if(verbose){cout << "Exiting function ALLOCATECREDIT in period " << t << endl;}

    PRODMACH();
    if(verbose){cout << "Exiting function PRODMACH in period " << t << endl;}

    EN_DEM_HOUSEHOLDS();
    if(verbose){cout << "Exiting function EN_DEM_HOUSEHOLDS in period " << t << endl;}

    EN_DEM_TOT();
    if(verbose){cout << "Exiting function EN_DEM in period " << t << endl;}

    EMISS_IND(); //Calculate total emissions from C-firms and K-firms
    if(verbose){cout << "Exiting function EMISS_IND in period " << t << endl;}

    ENERGY_INV_PROD();
    if(verbose){cout << "Exiting function ENERGY_INV_PROD in period " << t << endl;}

    ENERGY_RandD();
    if(verbose){cout << "Exiting function ENERGY_RandD in period " << t << endl;}

    PAY_LAB_INV();
    if(verbose){cout << "Exiting function PAY_LAB_INV in period " << t << endl;}

    COMPET2();
    if(verbose){cout << "Exiting function COMPET2 in period " << t << endl;}

    PROFIT();
    if(verbose){cout << "Exiting function PROFIT in period " << t << endl;}

    MACRO();
    if(verbose){cout << "Exiting function MACRO in period " << t << endl;}

    ENTRYEXIT();
    if(verbose){cout << "Exiting function ENTRYEXIT in period " << t << endl;}

    BANKING();
    if(verbose){cout << "Exiting function BANKING in period " << t << endl;}

    BAILOUT();
    if(verbose){cout << "Exiting function BAILOUT in period " << t << endl;}

    GOV_BUDGET();
    if(verbose){cout << "Exiting function GOV_BUDGET in period " << t << endl;}

    TAYLOR();
    if(verbose){cout << "Exiting function TAYLOR in period " << t << endl;}

    SETTLEMENT();
    if(verbose){cout << "Exiting function SETTLEMENT in period " << t << endl;}

    TECHANGEND();
    if(verbose){cout << "Exiting function TECHANGEND in period " << t << endl;}

    ENERGYINEQUALITY();
    if(verbose){cout << "Exiting function ENERGYINEQUALITY in period " << t << endl;}

    CARBONINEQUALITY();
    if(verbose){cout << "Exiting function CARBONINEQUALITY in period " << t << endl;}

    if(t>t_start_climbox && t%freqclim==0)
    {
      if (flag_cum_emissions==0)
      {

        CLIMATEBOX();
        if(iterations==niterclim)
        {
          // std::cerr << "\n\n ERROR: Carbon content loop did not converge in period " << t << endl;
          // Errors << "\n Carbon content loop did not converge in period " << t << endl;
        }

        if(verbose){cout << "Exiting function CLIMATEBOX in period " << t << endl;}
      }
      else
      {
        CLIMATEBOX_CUM_EMISS();
        if(verbose){cout << "Exiting function CLIMATEBOX_CUM_EMISS in period " << t << endl;}
      }
    }

    if(flag_exogenousshocks==1)
    {
      SINGLESHOCK();
      if(verbose){cout << "Exiting function SINGLESHOCK in period " << t << endl;}
    }
    else if(a2_nord>0)
    {
      SHOCKSNORD();
      if(verbose){cout << "Exiting function SHOCKSNORD in period " << t << endl;}
    }
    else
    {
      SHOCKS();
      if(verbose){cout << "Exiting function SHOCKS in period " << t << endl;}
    }

    DEPOSITCHECK();
    if(verbose){cout << "Exiting function DEPOSITCHECK in period " << t << endl;}

    NEGATIVITYCHECK();
    if(verbose){cout << "Exiting function NEGATIVITYCHECK in period " << t << endl;}

    CHECKSUMS();
    if(verbose){cout << "Exiting function CHECKSUMS in period " << t << endl;}

    ADJUSTSTOCKS();
    if(verbose){cout << "Exiting function ADJUSTSTOCKS in period " << t << endl;}

    SFC_CHECK();
    if(verbose){cout << "Exiting function SFC_CHECK in period " << t << endl;}

    SAVE();
    if(verbose){cout << "Exiting function SAVE in period " << t << endl;}

    UPDATE();
    if(verbose){cout << "Exiting function UPDATE in period " << t << endl;}

    UPDATECLIMATE();
    if(verbose){cout << "Exiting function UPDATECLIMATE in period " << t << endl;}

    OVERBOOST();

    if(verbose){cout << "Exiting function OVERBOOST; end of period " << t << endl;}
  }
  if(verbose){cout << "End of simulation loop"<< endl;}

  //Write output files
  WRITEFILE(nomefile2, output_string);

  //Close file for printing Errors
  Errors.close();

  //To measure execution time of entire run
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<seconds>(stop - start);
  cout << "Time taken by simulation " << seednumber << ": " << duration.count() << " seconds" << endl;

  return(EXIT_SUCCESS);
}

///////////MODEL FUNCTIONS/////////////////////////////

void SETPARAMS(const rapidjson::Document& inputs){

      //initialise households classes
      classes_mh.push_back("wr");
      classes_mh.push_back("pr");
      classes_mh.push_back("ma");
      n_classes_mh = classes_mh.size();

      //Read the values of parameters, flags, initial values and one-off shocks from the document "inputs" and set the respective variables to those values
      N1=inputs["params"][0]["N1"].GetInt();
      N2=inputs["params"][0]["N2"].GetInt();
      NB=inputs["params"][0]["NB"].GetInt();
      T=inputs["params"][0]["T"].GetInt();
      varphi=inputs["params"][0]["varphi"].GetDouble();
      nu=inputs["params"][0]["nu"].GetDouble();
      xi=inputs["params"][0]["xi"].GetDouble();
      o1=inputs["params"][0]["o1"].GetDouble();
      o2=inputs["params"][0]["o2"].GetDouble();
      uu11=inputs["params"][0]["uu11"].GetDouble();
      uu21=inputs["params"][0]["uu21"].GetDouble();
      uu12=inputs["params"][0]["uu12"].GetDouble();
      uu22=inputs["params"][0]["uu22"].GetDouble();
      uu31=inputs["params"][0]["uu31"].GetDouble();
      uu41=inputs["params"][0]["uu41"].GetDouble();
      uu32=inputs["params"][0]["uu32"].GetDouble();
      uu42=inputs["params"][0]["uu42"].GetDouble();
      uu51=inputs["params"][0]["uu51"].GetDouble();
      uu61=inputs["params"][0]["uu61"].GetDouble();
      uu52=inputs["params"][0]["uu52"].GetDouble();
      uu62=inputs["params"][0]["uu62"].GetDouble();
      uinf=inputs["params"][0]["uinf"].GetDouble();
      usup=inputs["params"][0]["usup"].GetDouble();
      b_a11=inputs["params"][0]["b_a11"].GetDouble();
      b_b11=inputs["params"][0]["b_b11"].GetDouble();
      b_a12=inputs["params"][0]["b_a12"].GetDouble();
      b_b12=inputs["params"][0]["b_b12"].GetDouble();
      b_a2=inputs["params"][0]["b_a2"].GetDouble();
      b_b2=inputs["params"][0]["b_b2"].GetDouble();
      b_a3=inputs["params"][0]["b_a3"].GetDouble();
      b_b3=inputs["params"][0]["b_b3"].GetDouble();
      mi1=inputs["params"][0]["mi1"].GetDouble();
      mi2=inputs["params"][0]["mi2"].GetDouble();
      Gamma=inputs["params"][0]["Gamma"].GetDouble();
      chi=inputs["params"][0]["chi"].GetDouble();
      omega1=inputs["params"][0]["omega1"].GetDouble();
      omega2=inputs["params"][0]["omega2"].GetDouble();
      omega4=inputs["params"][0]["omega4"].GetDouble();
      psi1=inputs["params"][0]["psi1"].GetDouble();
      psi2=inputs["params"][0]["psi2"].GetDouble();
      psi3=inputs["params"][0]["psi3"].GetDouble();
      deltami2=inputs["params"][0]["deltami2"].GetDouble();
      w_min=inputs["params"][0]["w_min"].GetDouble();
      pmin=inputs["params"][0]["pmin"].GetDouble();
      theta=inputs["params"][0]["theta"].GetDouble();
      u=inputs["params"][0]["u"].GetDouble();
      alfa=inputs["params"][0]["alfa"].GetDouble();
      b=inputs["params"][0]["b"].GetDouble();
      dim_mach=inputs["params"][0]["dim_mach"].GetDouble();
      agemax=inputs["params"][0]["agemax"].GetDouble();
      credit_multiplier=inputs["params"][0]["credit_multiplier"].GetDouble();
      beta_basel=inputs["params"][0]["beta_basel"].GetDouble();
      floor_default_probability=inputs["params"][0]["floor_default_probability"].GetDouble();
      upsilon=inputs["params"][0]["upsilon"].GetDouble();
      lambdaB1=inputs["params"][0]["lambdaB1"].GetDouble();
      lambdaB2=inputs["params"][0]["lambdaB2"].GetDouble();
      riskWeightLoans=inputs["params"][0]["riskWeightLoans"].GetDouble();
      riskWeightGovBonds=inputs["params"][0]["riskWeightGovBonds"].GetDouble();
      capitalAdequacyRatioTarget=inputs["params"][0]["capitalAdequacyRatioTarget"].GetDouble();
      bankmarkdown=inputs["params"][0]["bankmarkdown"].GetDouble();
      centralbankmarkdown=inputs["params"][0]["centralbankmarkdown"].GetDouble();
      d1=inputs["params"][0]["d1"].GetDouble();
      d2=inputs["params"][0]["d2"].GetDouble();
      db=inputs["params"][0]["db"].GetDouble();
      repayment_share=inputs["params"][0]["repayment_share"].GetDouble();
      bonds_share=inputs["params"][0]["bonds_share"].GetDouble();
      pareto_a=inputs["params"][0]["pareto_a"].GetDouble();
      pareto_k=inputs["params"][0]["pareto_k"].GetDouble();
      pareto_p=inputs["params"][0]["pareto_p"].GetDouble();
      d_cpi_target=inputs["params"][0]["d_cpi_target"].GetDouble();
      ustar=inputs["params"][0]["ustar"].GetDouble();
      w1sup=inputs["params"][0]["w1sup"].GetDouble();
      w1inf=inputs["params"][0]["w1inf"].GetDouble();
      w2sup=inputs["params"][0]["w2sup"].GetDouble();
      w2inf=inputs["params"][0]["w2inf"].GetDouble();
      k_const=inputs["params"][0]["k_const"].GetDouble();
      k_const2=inputs["params"][0]["k_const2"].GetDouble();
      aliqw_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["aliqw_mh"]);
      aliqwealth_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["aliqwealth_mh"]);      
      aliqdiv=inputs["params"][0]["aliqdiv"].GetDouble();
      taylor1=inputs["params"][0]["taylor1"].GetDouble();
      taylor2=inputs["params"][0]["taylor2"].GetDouble();
      bondsmarkdown=inputs["params"][0]["bondsmarkdown"].GetDouble();
      mdw=inputs["params"][0]["mdw"].GetDouble();
      phi2=inputs["params"][0]["phi2"].GetDouble();
      b1sup=inputs["params"][0]["b1sup"].GetDouble();
      b1inf=inputs["params"][0]["b1inf"].GetDouble();
      b2sup=inputs["params"][0]["b2sup"].GetDouble();
      b2inf=inputs["params"][0]["b2inf"].GetDouble();
      aliq=inputs["params"][0]["aliq"].GetDouble();
      aliqb=inputs["params"][0]["aliqb"].GetDouble();
      wu=inputs["params"][0]["wu"].GetDouble();
      de=inputs["params"][0]["de"].GetDouble();
      a1_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["a1_mh"]);
      a2_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["a2_mh"]);
      a3_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["a3_mh"]);
      a4_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["a4_mh"]);
      a5_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["a5_mh"]);
      f2_entry_min=inputs["params"][0]["f2_entry_min"].GetDouble();
      kappa=inputs["params"][0]["kappa"].GetDouble();
      taylor=inputs["params"][0]["taylor"].GetDouble();
      omicron=inputs["params"][0]["omicron"].GetDouble();
      I_max=inputs["params"][0]["I_max"].GetDouble();
      persistence=inputs["params"][0]["persistence"].GetDouble();
      omega3=inputs["params"][0]["omega3"].GetDouble();
      profit_share_energy_inv=inputs["params"][0]["profit_share_energy_inv"].GetDouble();
      redistribute_co2TaxRev=inputs["params"][0]["redistribute_co2TaxRev"].GetDouble();
      d_f=inputs["params"][0]["d_f"].GetDouble();
      d_f_sh_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["d_f_sh_mh"]);
      g_ls=inputs["params"][0]["g_ls"].GetDouble();
      ld_ratios_mh=read_input_vector_into_umap(classes_mh,inputs["params"][0]["ld_ratios_mh"]);
      aliqee=inputs["params"][0]["aliqee"].GetDouble();
      aliqef=inputs["params"][0]["aliqef"].GetDouble();
      tre=inputs["params"][0]["tre"].GetDouble();
      tref=inputs["params"][0]["tref"].GetDouble();
      passthrough=inputs["params"][0]["passthrough"].GetDouble();
      r_base=inputs["params"][0]["r_base"].GetDouble();
      g_cons_GDP_ratio=inputs["params"][0]["g_cons_GDP_ratio"].GetDouble();
      LS_sh_mh=read_input_vector_into_umap(classes_mh, inputs["params"][0]["LS_sh_mh"]);
      w_ratios_mh=read_input_vector_into_umap(classes_mh, inputs["params"][0]["w_ratios_mh"]);
      Ownership_sh_e_mh = read_input_vector_into_umap(classes_mh, inputs["params"][0]["Ownership_sh_e_mh"]);
      Ownership_sh_ff_mh = read_input_vector_into_umap(classes_mh, inputs["params"][0]["Ownership_sh_ff_mh"]);
      Ownership_sh_b_mh = read_input_vector_into_umap(classes_mh, inputs["params"][0]["Ownership_sh_b_mh"]);
      Ownership_sh_1_mh0 = read_input_vector_into_umap(classes_mh, inputs["params"][0]["Ownership_sh_1_mh0"]);
      Ownership_sh_2_mh0 = read_input_vector_into_umap(classes_mh, inputs["params"][0]["Ownership_sh_2_mh0"]);
      energy_expenditure_sh_mh = read_input_vector_into_umap(classes_mh, inputs["params"][0]["energy_expenditure_sh_mh"]);
      ratio_c_en_h_firms=inputs["params"][0]["ratio_c_en_h_firms"].GetDouble();
      demand_persistency_h=inputs["params"][0]["demand_persistency_h"].GetDouble();

      //Parameters for experiments
      t_regime_shifts = inputs["expparams"][0]["t_regime_shifts"].GetInt();
      income_sh_regime_shift_mh = read_input_vector_into_umap(classes_mh, inputs["expparams"][0]["income_sh_regime_shift_mh"]);
      a1_mh_regime_shift = read_input_vector_into_umap(classes_mh, inputs["expparams"][0]["a1_mh_regime_shift"]);
      a2_mh_regime_shift = read_input_vector_into_umap(classes_mh, inputs["expparams"][0]["a2_mh_regime_shift"]);
      a3_mh_regime_shift = read_input_vector_into_umap(classes_mh, inputs["expparams"][0]["a3_mh_regime_shift"]);
      en_exp_sh_mh_regime_shift = read_input_vector_into_umap(classes_mh, inputs["expparams"][0]["en_exp_sh_mh_regime_shift"]);
      c_price_increase_perc = inputs["expparams"][0]["c_price_increase_perc"].GetDouble();
      ratio_mi_en_shock = inputs["expparams"][0]["ratio_mi_en_shock"].GetDouble();
      aliqdiv_regime_shift = inputs["expparams"][0]["aliqdiv_regime_shift"].GetDouble();
      aliqw_regime_shift_mh = read_input_vector_into_umap(classes_mh, inputs["expparams"][0]["aliqw_regime_shift_mh"]);
      aliqw_progressivity_regime_shift = inputs["expparams"][0]["aliqw_progressivity_regime_shift"].GetDouble();
      aliq_households_increase = inputs["expparams"][0]["aliq_households_increase"].GetDouble();
      aliqwealth_regime_shift_mh = read_input_vector_into_umap(classes_mh, inputs["expparams"][0]["aliqwealth_regime_shift_mh"]);
      g_cons_GDP_ratio_regime_shift = inputs["expparams"][0]["g_cons_GDP_ratio_regime_shift"].GetDouble();
      env_subsidy_per_machine = inputs["expparams"][0]["env_subsidy_per_machine"].GetDouble();
      correlation_prod_and_green_tech = inputs["expparams"][0]["correlation_prod_and_green_tech"].GetDouble();
      rs_uu_lp = inputs["expparams"][0]["rs_uu_lp"].GetDouble();
      rs_uu_ee = inputs["expparams"][0]["rs_uu_ee"].GetDouble();
      rs_uu_ef = inputs["expparams"][0]["rs_uu_ef"].GetDouble();
      K_ge_END_perc = inputs["expparams"][0]["K_ge_END_perc"].GetDouble();
      t_length_energy_transition = inputs["expparams"][0]["t_length_energy_transition"].GetInt();
      renew_impact_on_p_e = inputs["expparams"][0]["renew_impact_on_p_e"].GetDouble();
      bonuses_share = inputs["expparams"][0]["bonuses_share"].GetDouble();


      //Energy & climate parameters
      share_RD_en=inputs["climparams"][0]["share_RD_en"].GetDouble();
      share_de_0=inputs["climparams"][0]["share_de_0"].GetDouble();
      payback_en=inputs["climparams"][0]["payback_en"].GetInt();
      life_plant=inputs["climparams"][0]["life_plant"].GetInt();
      exp_quota=inputs["climparams"][0]["exp_quota"].GetDouble();
      o1_en=inputs["climparams"][0]["o1_en"].GetDouble();
      uu1_en=inputs["climparams"][0]["uu1_en"].GetDouble();
      uu2_en=inputs["climparams"][0]["uu2_en"].GetDouble();
      Transfer_shock_sh_mh=read_input_vector_into_umap(classes_mh, inputs["climparams"][0]["Transfer_shock_sh_mh"]);
      exp_quota_param=inputs["params"][0]["exp_quota_param"].GetDouble();
      t_start_climbox=inputs["climparams"][0]["t_start_climbox"].GetInt();
      T_pre=inputs["climparams"][0]["T_pre"].GetDouble();
      intercept_temp=inputs["climparams"][0]["intercept_temp"].GetDouble();
      slope_temp=inputs["climparams"][0]["slope_temp"].GetDouble();
      tc1=inputs["climparams"][0]["tc1"].GetDouble();
      tc2=inputs["climparams"][0]["tc2"].GetDouble();
      ndep=inputs["climparams"][0]["ndep"].GetInt();
      laydep.ReSize(ndep);
      laydep(1)=inputs["climparams"][0]["laydep1"].GetDouble();
      laydep(2)=inputs["climparams"][0]["laydep2"].GetDouble();
      laydep(3)=inputs["climparams"][0]["laydep3"].GetDouble();
      laydep(4)=inputs["climparams"][0]["laydep4"].GetDouble();
      laydep(5)=inputs["climparams"][0]["laydep5"].GetDouble();
      fertil=inputs["climparams"][0]["fertil"].GetDouble();
      heatstress=inputs["climparams"][0]["heatstress"].GetDouble();
      humtime=inputs["climparams"][0]["humtime"].GetDouble();
      biotime=inputs["climparams"][0]["biotime"].GetDouble();
      humfrac=inputs["climparams"][0]["humfrac"].GetDouble();
      eddydif=inputs["climparams"][0]["eddydif"].GetDouble();
      ConrefT=inputs["climparams"][0]["ConrefT"].GetDouble();
      rev0=inputs["climparams"][0]["rev0"].GetDouble();
      revC=inputs["climparams"][0]["revC"].GetDouble();
      niterclim=inputs["climparams"][0]["niterclim"].GetInt();
      forCO2=inputs["climparams"][0]["forCO2"].GetDouble();
      otherforcefac=inputs["climparams"][0]["otherforcefac"].GetDouble();
      outrad=inputs["climparams"][0]["outrad"].GetDouble();
      secyr=inputs["climparams"][0]["secyr"].GetDouble();
      seasurf=inputs["climparams"][0]["seasurf"].GetDouble();
      heatcap=inputs["climparams"][0]["heatcap"].GetDouble();
      freqclim=inputs["climparams"][0]["freqclim"].GetInt();
      g_emiss_global=inputs["climparams"][0]["g_emiss_global"].GetDouble();
      emiss_share=inputs["climparams"][0]["emiss_share"].GetDouble();

      nshocks=inputs["climshockparams"][0]["nshocks"].GetInt();
      a_0.ReSize(nshocks);
      a_0(1)=inputs["climshockparams"][0]["a1_0"].GetDouble();
      a_0(2)=inputs["climshockparams"][0]["a2_0"].GetDouble();
      a_0(3)=inputs["climshockparams"][0]["a3_0"].GetDouble();
      a_0(4)=inputs["climshockparams"][0]["a4_0"].GetDouble();
      a_0(5)=inputs["climshockparams"][0]["a5_0"].GetDouble();
      a_0(6)=inputs["climshockparams"][0]["a6_0"].GetDouble();
      a_0(7)=inputs["climshockparams"][0]["a7_0"].GetDouble();
      a_0(8)=inputs["climshockparams"][0]["a8_0"].GetDouble();
      a_0(9)=inputs["climshockparams"][0]["a9_0"].GetDouble();
      b_0.ReSize(nshocks);
      b_0(1)=inputs["climshockparams"][0]["b1_0"].GetDouble();
      b_0(2)=inputs["climshockparams"][0]["b2_0"].GetDouble();
      b_0(3)=inputs["climshockparams"][0]["b3_0"].GetDouble();
      b_0(4)=inputs["climshockparams"][0]["b4_0"].GetDouble();
      b_0(5)=inputs["climshockparams"][0]["b5_0"].GetDouble();
      b_0(6)=inputs["climshockparams"][0]["b6_0"].GetDouble();
      b_0(7)=inputs["climshockparams"][0]["b7_0"].GetDouble();
      b_0(8)=inputs["climshockparams"][0]["b8_0"].GetDouble();
      b_0(9)=inputs["climshockparams"][0]["b9_0"].GetDouble();
      shockexponent1.ReSize(nshocks);
      shockexponent1(1)=inputs["climshockparams"][0]["shockexponent1_1"].GetDouble();
      shockexponent1(2)=inputs["climshockparams"][0]["shockexponent2_1"].GetDouble();
      shockexponent1(3)=inputs["climshockparams"][0]["shockexponent3_1"].GetDouble();
      shockexponent1(4)=inputs["climshockparams"][0]["shockexponent4_1"].GetDouble();
      shockexponent1(5)=inputs["climshockparams"][0]["shockexponent5_1"].GetDouble();
      shockexponent1(6)=inputs["climshockparams"][0]["shockexponent6_1"].GetDouble();
      shockexponent1(7)=inputs["climshockparams"][0]["shockexponent7_1"].GetDouble();
      shockexponent1(8)=inputs["climshockparams"][0]["shockexponent8_1"].GetDouble();
      shockexponent1(9)=inputs["climshockparams"][0]["shockexponent9_1"].GetDouble();
      shockexponent2.ReSize(nshocks);
      shockexponent2(1)=inputs["climshockparams"][0]["shockexponent1_2"].GetDouble();
      shockexponent2(2)=inputs["climshockparams"][0]["shockexponent2_2"].GetDouble();
      shockexponent2(3)=inputs["climshockparams"][0]["shockexponent3_2"].GetDouble();
      shockexponent2(4)=inputs["climshockparams"][0]["shockexponent4_2"].GetDouble();
      shockexponent2(5)=inputs["climshockparams"][0]["shockexponent5_2"].GetDouble();
      shockexponent2(6)=inputs["climshockparams"][0]["shockexponent6_2"].GetDouble();
      shockexponent2(7)=inputs["climshockparams"][0]["shockexponent7_2"].GetDouble();
      shockexponent2(8)=inputs["climshockparams"][0]["shockexponent8_2"].GetDouble();
      shockexponent2(9)=inputs["climshockparams"][0]["shockexponent9_2"].GetDouble();
      a2_nord=inputs["climshockparams"][0]["a2_nord"].GetDouble();
      sd_nord=inputs["climshockparams"][0]["sd_nord"].GetDouble();

      //Inits
      A0=inputs["inits"][0]["A0"].GetDouble();
      LS0_wr=inputs["inits"][0]["LS0_wr"].GetDouble();
      W10=inputs["inits"][0]["W10"].GetDouble();
      W20=inputs["inits"][0]["W20"].GetDouble();
      L0=inputs["inits"][0]["L0"].GetDouble();
      w0_wr=inputs["inits"][0]["w0_wr"].GetDouble();
      K0=inputs["inits"][0]["K0"].GetDouble();
      bankmarkup_init=inputs["inits"][0]["bankmarkup_init"].GetDouble();
      FirmDefaultProbability_init=inputs["inits"][0]["FirmDefaultProbability_init"].GetDouble();
      A0_en=inputs["inits"][0]["A0_en"].GetDouble();
      A0_ef=inputs["inits"][0]["A0_ef"].GetInt();
      A_ef_ratio0=inputs["inits"][0]["A_ef_ratio0"].GetDouble();
      K_ge0_perc=inputs["inits"][0]["K_ge0_perc"].GetDouble();
      pf0=inputs["inits"][0]["pf0"].GetDouble();
      mi_en0=inputs["inits"][0]["mi_en0"].GetDouble();
      A_de0=inputs["inits"][0]["A_de0"].GetDouble();
      EM0=inputs["inits"][0]["EM0"].GetDouble();
      CF_ge0=inputs["inits"][0]["CF_ge0"].GetDouble();
      t_CO2_0=inputs["inits"][0]["t_CO2_0"].GetDouble();
      t_CO2_en_0=inputs["inits"][0]["t_CO2_en_0"].GetDouble();
      d_Am_init=inputs["inits"][0]["d_Am_init"].GetDouble();
      D_h0=inputs["inits"][0]["D_h0"].GetDouble();
      D_h0_sh_mh=read_input_vector(D_h0_sh_mh, inputs["inits"][0]["D_h0_sh_mh"], n_classes_mh);
      NW_b0=inputs["inits"][0]["NW_b0"].GetDouble();
      pm=inputs["inits"][0]["pm"].GetDouble();
      D_e0=inputs["inits"][0]["D_e0"].GetDouble();

      Emiss_yearly_0=inputs["climinits"][0]["Emiss_yearly_0"].GetDouble();
      Cum_emissions_0=inputs["climinits"][0]["Cum_emissions_0"].GetDouble();
      T_0_cumemiss=inputs["climinits"][0]["T_0_cumemiss"].GetDouble();
      Con00=inputs["climinits"][0]["Con00"].GetDouble();
      Conref=Con00*laydep(1);
      NPP0=inputs["climinits"][0]["NPP0"].GetDouble();
      Cat0=inputs["climinits"][0]["Cat0"].GetDouble();
      Catinit0=inputs["climinits"][0]["Catinit0"].GetDouble();
      Honinit0.ReSize(ndep);
      Coninit0.ReSize(ndep);
      Honinit0(1)=inputs["climinits"][0]["Honinit01"].GetDouble();
      Honinit0(2)=inputs["climinits"][0]["Honinit02"].GetDouble();
      Honinit0(3)=inputs["climinits"][0]["Honinit03"].GetDouble();
      Honinit0(4)=inputs["climinits"][0]["Honinit04"].GetDouble();
      Honinit0(5)=inputs["climinits"][0]["Honinit05"].GetDouble();
      Coninit0(1)=inputs["climinits"][0]["Coninit01"].GetDouble();
      Coninit0(2)=inputs["climinits"][0]["Coninit02"].GetDouble();
      Coninit0(3)=inputs["climinits"][0]["Coninit03"].GetDouble();
      Coninit0(4)=inputs["climinits"][0]["Coninit04"].GetDouble();
      Coninit0(5)=inputs["climinits"][0]["Coninit05"].GetDouble();
      Tmixedinit0=inputs["climinits"][0]["Tmixedinit0"].GetDouble();
      biominit0=inputs["climinits"][0]["biominit0"].GetDouble();
      huminit0=inputs["climinits"][0]["huminit0"].GetDouble();
      Catinit1=inputs["climinits"][0]["Catinit1"].GetDouble();
      Honinit1.ReSize(ndep);
      Coninit1.ReSize(ndep);
      Honinit1(1)=inputs["climinits"][0]["Honinit11"].GetDouble();
      Honinit1(2)=inputs["climinits"][0]["Honinit12"].GetDouble();
      Honinit1(3)=inputs["climinits"][0]["Honinit13"].GetDouble();
      Honinit1(4)=inputs["climinits"][0]["Honinit14"].GetDouble();
      Honinit1(5)=inputs["climinits"][0]["Honinit15"].GetDouble();
      Coninit1(1)=inputs["climinits"][0]["Coninit11"].GetDouble();
      Coninit1(2)=inputs["climinits"][0]["Coninit12"].GetDouble();
      Coninit1(3)=inputs["climinits"][0]["Coninit13"].GetDouble();
      Coninit1(4)=inputs["climinits"][0]["Coninit14"].GetDouble();
      Coninit1(5)=inputs["climinits"][0]["Coninit15"].GetDouble();
      Tmixedinit1=inputs["climinits"][0]["Tmixedinit1"].GetDouble();
      biominit1=inputs["climinits"][0]["biominit1"].GetDouble();
      huminit1=inputs["climinits"][0]["huminit1"].GetDouble();

      //flags
      flag_exogenousshocks=inputs["flags"][0]["flag_exogenousshocks"].GetInt();
      flag_cum_emissions=inputs["flags"][0]["flag_cum_emissions"].GetInt();
      flag_tax_CO2=inputs["flags"][0]["flag_tax_CO2"].GetInt();
      flag_capshocks=inputs["flags"][0]["flag_capshocks"].GetInt();
      flag_outputshocks=inputs["flags"][0]["flag_outputshocks"].GetInt();
      flag_inventshocks=inputs["flags"][0]["flag_inventshocks"].GetInt();
      flag_encapshocks=inputs["flags"][0]["flag_encapshocks"].GetInt();
      flag_popshocks=inputs["flags"][0]["flag_popshocks"].GetInt();
      flag_energyshocks=inputs["flags"][0]["flag_energyshocks"].GetInt();
      flag_energyshocks_MP=inputs["flags"][0]["flag_energyshocks_MP"].GetInt();
      flag_demandshocks=inputs["flags"][0]["flag_demandshocks"].GetInt();
      flag_RDshocks=inputs["flags"][0]["flag_RDshocks"].GetInt();
      flag_prodshocks1=inputs["flags"][0]["flag_prodshocks1"].GetInt();
      flag_prodshocks2=inputs["flags"][0]["flag_prodshocks2"].GetInt();
      flag_share_END=inputs["flags"][0]["flag_share_END"].GetInt();
      flag_energy_exp=inputs["flags"][0]["flag_energy_exp"].GetInt();
      flagbailout=inputs["flags"][0]["flagbailout"].GetInt();
      flag_entry=inputs["flags"][0]["flag_entry"].GetInt();
      flag_nonCO2_force=inputs["flags"][0]["flag_nonCO2_force"].GetInt();
      flag_inventories=inputs["flags"][0]["flag_inventories"].GetInt();
      flag_uniformshocks=inputs["flags"][0]["flag_uniformshocks"].GetInt();
      flag_rate_setting_markup=inputs["flags"][0]["flag_rate_setting_markup"].GetInt();
      flag_rate_setting_loans=inputs["flags"][0]["flag_rate_setting_loans"].GetInt();
      flag_endogenous_exp_quota=inputs["flags"][0]["flag_endogenous_exp_quota"].GetInt();
      flag_c_price_increase=inputs["flags"][0]["flag_c_price_increase"].GetInt();
      flag_change_income_shares=inputs["flags"][0]["flag_change_income_shares"].GetInt();
      flag_change_apc=inputs["flags"][0]["flag_change_apc"].GetInt();
      flag_change_en_exp_shares=inputs["flags"][0]["flag_change_en_exp_shares"].GetInt();
      flag_change_dividends_tax_rate=inputs["flags"][0]["flag_change_dividends_tax_rate"].GetInt();
      flag_change_wage_tax_rate=inputs["flags"][0]["flag_change_wage_tax_rate"].GetInt();
      flag_highly_progressive_taxation=inputs["flags"][0]["flag_highly_progressive_taxation"].GetInt();
      flag_change_wealth_tax_rate=inputs["flags"][0]["flag_change_wealth_tax_rate"].GetInt();
      flag_change_public_spending=inputs["flags"][0]["flag_change_public_spending"].GetInt();
      flag_env_friendl_in_market_shares=inputs["flags"][0]["flag_env_friendl_in_market_shares"].GetInt();
      flag_environmental_subsidies_C_firms=inputs["flags"][0]["flag_environmental_subsidies_C_firms"].GetInt();
      flag_correlate_prod_and_green_tech=inputs["flags"][0]["flag_correlate_prod_and_green_tech"].GetInt();

      shocks_cfirms.ReSize(N2);
      shocks_kfirms.ReSize(N1);
      for(i=1; i<=N1; i++){
        shocks_kfirms(i)=inputs["shocks_kfirms"][(i-1)].GetDouble();
      }
      for(j=1; j<=N2; j++){
        shocks_cfirms(j)=inputs["shocks_cfirms"][(j-1)].GetDouble();
      }
      shock_scalar=inputs["shock_scalar"][0].GetDouble();
}

void RESIZE(void)
{
  //Resize all vectors, matrices and arrays to the correct dimension
  //Required because the correct dimensions depend on the number of agents, which are parameters
  X_a.ReSize(nshocks);
  X_b.ReSize(nshocks);
  A.ReSize(T,N1);
  C.ReSize(T,N1);
  C_secondhand.ReSize(T,N1);
  A_en.ReSize(T,N1);
  A_ef.ReSize(T,N1);
  A_de.ReSize(T);
  EM_de.ReSize(T);
  C_de.ReSize(T);
  G_de.ReSize(T);
  G_ge.ReSize(T);
  G_ge_n.ReSize(T);
  CF_ge.ReSize(T);
  IC_en_quota.ReSize(T);
  G_de_temp.ReSize(T);
  shocks_encapstock_de.ReSize(T);
  shocks_encapstock_ge.ReSize(T);

  //C-firms
  D2.ReSize(2,N2);
  De.ReSize(N2);
  f2.ReSize(3,N2);
  E2.ReSize(N2);
  c2.ReSize(N2);
  c2p.ReSize(N2);
  Q2.ReSize(N2);
  N.ReSize(2,N2);
  Inventories.ReSize(2,N2);
  Deposits_2.ReSize(2,N2);
  NW_2.ReSize(2,N2);
  NW_2_c.ReSize(N2);
  CapitalStock.ReSize(2,N2);
  deltaCapitalStock.ReSize(2,N2);
  DebtServiceToSales2.ReSize(N2);
  scrap_age.ReSize(N2);
  I.ReSize(N2);
  EI.ReSize(2,N2);
  EI_n.ReSize(N2);
  SI.ReSize(N2);
  SI_n.ReSize(N2);
  S2.ReSize(2,N2);
  Sales2.ReSize(N2);
  Loans_2.ReSize(2,N2);
  DebtService_2.ReSize(2,N2);
  CreditDemand.ReSize(N2);
  p2.ReSize(N2);
  DebtServiceToSales2_temp.ReSize(N2);
  mu2.ReSize(2,N2);
  LoanInterest_2.ReSize(N2);
  DebtRemittances2.ReSize(N2);
  baddebt_2.ReSize(N2);
  EId.ReSize(N2);
  SId.ReSize(N2);
  SIp.ReSize(N2);
  EIp.ReSize(N2);
  Ip.ReSize(N2);
  A2.ReSize(N2);
  A2_mprod.ReSize(N2);
  A2e.ReSize(N2);
  c2e.ReSize(N2);
  Ld2_wr.ReSize(N2);
  l2.ReSize(N2);
  n_mach.ReSize(N2);
  Qd.ReSize(N2);
  K.ReSize(N2);
  K_cur.ReSize(N2);
  Kd.ReSize(N2);
  Ktrig.ReSize(N2);
  Pi2.ReSize(2,N2);
  Q2temp.ReSize(N2);
  Wages_2_i.ReSize(N2);
  Investment_2.ReSize(N2);
  EnergyPayments_2.ReSize(N2);
  InterestDeposits_2.ReSize(N2);
  Taxes_2.ReSize(N2);
  Taxes_CO2_2.ReSize(N2);
  f_temp2.ReSize(N2);
  D_temp2.ReSize(N2);
  dN.ReSize(N2);
  dNm.ReSize(N2);
  fornit.ReSize(N2);
  Cmach.ReSize(N2);
  CmachEI.ReSize(N2);
  CmachSI.ReSize(N2);
  Ne.ReSize(N2);
  mol.ReSize(N2);
  Dividends_2_i.ReSize(N2);
  Bonuses_2_i.ReSize(N2);
  Match.ReSize(N2,N1);
  BankingSupplier_2.ReSize(N2);
  r_deb_h.ReSize(N2);
  FirmDefaultProbability.ReSize(N2);
  k.ReSize(N2);
  D2_en.ReSize(N2);
  A2e_en.ReSize(N2);
  A2e_ef.ReSize(N2);
  A2e2.ReSize(N2);
  A2e_en2.ReSize(N2);
  A2e_ef2.ReSize(N2);
  A2_en.ReSize(N2);
  A2_ef.ReSize(N2);
  A2_ci.ReSize(N2);
  Emiss2.ReSize(N2);
  Injection_2.ReSize(N2);
  exiting_2.ReSize(N2);
  exit_payments2.ReSize(N2);
  exit_equity2.ReSize(N2);
  exit_marketshare2.ReSize(N2);
  n_mach_entry.ReSize(N2);
  shocks_capstock.ReSize(N2);
  shocks_invent.ReSize(N2);
  Loss_Capital.ReSize(N2);
  Loss_Capital_mol.ReSize(N2);
  Loss_Inventories.ReSize(N2);
  k_entry.ReSize(N2);
  EntryShare.ReSize(N2);
  CompEntry.ReSize(N2);
  K_temp.ReSize(N2);
  K_loss.ReSize(N2);
  C_loss.ReSize(N2);
  I_loss.ReSize(N1);
  marker_age.ReSize(N2);
  pass_2.ReSize(N2);
  Transfer_shock_f2.ReSize(N2);
  EnvSubsidies_2_i.ReSize(N2);

  //K-firms
  pass_1.ReSize(N1);
  Transfer_shock_f1.ReSize(N1);
  Deposits_1.ReSize(2,N1);
  NW_1.ReSize(2,N1);
  NW_1_c.ReSize(N1);
  Wages_1_i.ReSize(N1);
  S1.ReSize(N1);
  S1_pre.ReSize(N1);
  S1_post.ReSize(N1);
  Sales1.ReSize(N1);
  p1.ReSize(N1);
  RD.ReSize(2,N1);
  f1.ReSize(2,N1);
  Q1.ReSize(N1);
  D1.ReSize(N1);
  A1.ReSize(N1);
  Anew.ReSize(N1);
  A1inn.ReSize(N1);
  A1pinn.ReSize(N1);
  A1imm.ReSize(N1);
  A1pimm.ReSize(N1);
  Pi1.ReSize(2,N1);
  Ld1_prod_wr.ReSize(N1);
  ee1.ReSize(N1);
  nclient.ReSize(N1);
  Dividends_1_i.ReSize(N1);
  Bonuses_1_i.ReSize(N1);
  EnergyPayments_1.ReSize(N1);
  InterestDeposits_1.ReSize(N1);
  Taxes_1.ReSize(N1);
  Taxes_CO2_1.ReSize(N1);
  BankingSupplier_1.ReSize(N1);
  c1.ReSize(N1);
  A1p.ReSize(N1);
  RDin.ReSize(N1);
  RDim.ReSize(N1);
  Inn.ReSize(N1);
  Imm.ReSize(N1);
  A1f.ReSize(N1);
  A1pf.ReSize(N1);
  Td.ReSize(N1+1);
  D1_en.ReSize(N1);
  A1_en.ReSize(N1);
  A1_ef.ReSize(N1);
  A1_ci.ReSize(N1);
  A1p_en.ReSize(N1);
  A1p_ef.ReSize(N1);
  EE_inn.ReSize(N1);
  EEp_inn.ReSize(N1);
  EE_imm.ReSize(N1);
  EEp_imm.ReSize(N1);
  EF_inn.ReSize(N1);
  EFp_inn.ReSize(N1);
  EF_imm.ReSize(N1);
  EFp_imm.ReSize(N1);
  Emiss1.ReSize(N1);
  Injection_1.ReSize(N1);
  Balances_1.ReSize(N1);
  baddebt_1.ReSize(N1);
  exiting_1.ReSize(N1);
  exiting_1_payments.ReSize(N1);
  c1p.ReSize(N1);
  EnvSubsidiesPerMachine_1_i.ReSize(N1);
  Attractivess_1_i.ReSize(N1);


  shocks_machprod.ReSize(N1);
  shocks_techprod.ReSize(N1);
  shocks_machprod_max.ReSize(N1);
  shocks_techprod_max.ReSize(N1);
  shocks_rd.ReSize(N1);
  shocks_labprod1.ReSize(N1);
  shocks_labprod2.ReSize(N2);
  shocks_eneff1.ReSize(N1);
  shocks_eneff2.ReSize(N2);
  shocks_output1.ReSize(N1);
  shocks_output2.ReSize(N2);
  S1_temp.ReSize(2,N1);
  S2_temp.ReSize(2,N2);

  // Resize household classes variables
  for (const string& cl:classes_mh){
    //With lag
    Income_gross_mh[cl].ReSize(2);
    Taxes_div_mh[cl].ReSize(2);
    Taxes_bon_mh[cl].ReSize(2);
    Bonuses_mh[cl].ReSize(2);

    //For C- and K-firms
    Ld2_i_mh[cl].ReSize(N2);
    Wages_2_i_mh[cl].ReSize(N2);
    Dividends_2_i_mh[cl].ReSize(N2);
    Ownership_sh_2_i_mh[cl].ReSize(N2);
    Ownership_2_i_mh[cl].ReSize(N2);

    Ld1_i_mh[cl].ReSize(N1);
    Ld1_prod_i_mh[cl].ReSize(N1);
    Ld1_rd_i_mh[cl].ReSize(N1);
    Wages_1_i_mh[cl].ReSize(N1);
    Dividends_1_i_mh[cl].ReSize(N1);
    Ownership_sh_1_i_mh[cl].ReSize(N1);
    Ownership_1_i_mh[cl].ReSize(N1);
  }

  // Define number of banks
  if(N2>=200){NB=10;}
  else{NB=1;}

  NbClient_1.ReSize(NB);
  NbClient_2.ReSize(NB);
  bonds_dem.ReSize(NB);
  DebtServiceToSales2_bank.ReSize(N2,NB);
  DS2_rating.ReSize(N2,NB);
  fB.ReSize(2,NB);
  NW_b.ReSize(2,NB);
  NW_b_c.ReSize(NB);
  Loans_b.ReSize(2,NB);
  BankMatch_1.ReSize(N1,NB);
  BankMatch_2.ReSize(N2,NB);
  BankCredit.ReSize(NB);
  Taxes_b.ReSize(NB);
  BaselBankCredit.ReSize(NB);
  InterestDeposits.ReSize(NB);
  Deposits_b.ReSize(2,NB);
  BankProfits.ReSize(2,NB);
  Dividends_b_i.ReSize(NB);
  Bonuses_b_i.ReSize(NB);
  Bailout_b.ReSize(NB);
  LossAbsorbed.ReSize(NB);
  BankEquity_temp.ReSize(NB);
  Bank_active.ReSize(NB);
  GB_b.ReSize(2,NB);
  bonds_purchased.ReSize(NB);
  BondRepayments_b.ReSize(NB);
  NL_1.ReSize(NB);
  NL_2.ReSize(NB);
  BankProfits_temp.ReSize(NB);
  r_deb.ReSize(NB);
  bankmarkup.ReSize(NB);
  riskWeightedAssets.ReSize(NB);
  capitalAdequacyRatio.ReSize(NB);
  buffer.ReSize(NB);
  Bond_share.ReSize(NB);
  Deposits_hb.ReSize(2,NB);
  Deposits_eb.ReSize(2,NB);
  Advances_b.ReSize(2,NB);
  Reserves_b.ReSize(2,NB);
  InterestBonds_b.ReSize(NB);
  LoanInterest.ReSize(NB);
  InterestReserves_b.ReSize(NB);
  InterestAdvances_b.ReSize(NB);
  Outflows.ReSize(NB);
  Inflows.ReSize(NB);
  DepositShare_e.ReSize(NB);
  DepositShare_h.ReSize(NB);
  baddebt_b.ReSize(NB);
  capitalRecovered.ReSize(NB);
  capitalRecovered2.ReSize(NB);
  capitalRecoveredShare.ReSize(NB);
  LossEntry_b.ReSize(NB);
  ReserveBalance.ReSize(NB);
  ShareBonds.ReSize(NB);
  ShareReserves.ReSize(NB);
  ShareAdvances.ReSize(NB);
  Adjustment.ReSize(NB);
  prior.ReSize(NB);

  fluxC.ReSize((ndep-1));
  Con.ReSize(2,ndep);
  Hon.ReSize(2,ndep);
  Ton.ReSize(2,ndep);
  Cax.ReSize(niterclim);
  Caxx.ReSize(niterclim);
  Cay.ReSize(niterclim);
  Cayy.ReSize(niterclim);
  Caa.ReSize(niterclim);
  fluxH.ReSize((ndep-1));
  Emiss_TOT.ReSize((freqclim*2));

  //Initialising all resized variables to zero to avoid memory issues
  X_a=0;
  X_b=0;
  A=0;
  C=0;
  C_secondhand=0;
  A_en=0;
  A_ef=0;
  A_de=0;
  EM_de=0;
  C_de=0;
  G_de=0;
  G_ge=0;
  G_ge_n=0;
  CF_ge=0;
  IC_en_quota=0;
  G_de_temp=0;
  shocks_encapstock_de=0;
  shocks_encapstock_ge=0;
  D2=0;
  De=0;
  f2=0;
  E2=0;
  c2=0;
  c2p=0;
  Q2=0;
  N=0;
  Inventories=0;
  Deposits_2=0;
  NW_2=0;
  NW_2_c=0;
  CapitalStock=0;
  deltaCapitalStock=0;
  DebtServiceToSales2=0;
  scrap_age=0;
  I=0;
  EI=0;
  EI_n=0;
  SI=0;
  SI_n=0;
  S2=0;
  Sales2=0;
  Loans_2=0;
  DebtService_2=0;
  CreditDemand=0;
  p2=0;
  DebtServiceToSales2_temp=0;
  mu2=0;
  LoanInterest_2=0;
  DebtRemittances2=0;
  baddebt_2=0;
  EId=0;
  SId=0;
  SIp=0;
  EIp=0;
  Ip=0;
  A2=0;
  A2_mprod=0;
  A2e=0;
  c2e=0;
  Ld2_wr=0;
  Ld2_control=0;
  l2=0;
  n_mach=0;
  Qd=0;
  K=0;
  K_cur=0;
  Kd=0;
  Ktrig=0;
  Pi2=0;
  Q2temp=0;
  Investment_2=0;
  EnergyPayments_2=0;
  InterestDeposits_2=0;
  Taxes_2=0;
  Taxes_CO2_2=0;
  f_temp2=0;
  D_temp2=0;
  dN=0;
  dNm=0;
  fornit=0;
  Cmach=0;
  CmachEI=0;
  CmachSI=0;
  Ne=0;
  mol=0;
  Dividends_2=0;
  Match=0;
  BankingSupplier_2=0;
  r_deb_h=0;
  FirmDefaultProbability=0;
  k=0;
  D2_en=0;
  A2e_en=0;
  A2e_ef=0;
  A2e2=0;
  A2e_en2=0;
  A2e_ef2=0;
  A2_en=0;
  A2_ef=0;
  Emiss2=0;
  Injection_2=0;
  exiting_2=0;
  exit_payments2=0;
  exit_equity2=0;
  exit_marketshare2=0;
  n_mach_entry=0;
  shocks_capstock=0;
  shocks_invent=0;
  Loss_Capital=0;
  Loss_Capital_mol=0;
  Loss_Inventories=0;
  k_entry=0;
  EntryShare=0;
  CompEntry=0;
  K_temp=0;
  K_loss=0;
  C_loss=0;
  I_loss=0;
  marker_age=0;
  pass_2=0;
  Transfer_shock_f2=0;
  pass_1=0;
  Transfer_shock_f1=0;
  Deposits_1=0;
  NW_1=0;
  NW_1_c=0;
  S1=0;
  S1_pre=0;
  S1_post=0;
  Sales1=0;
  p1=0;
  RD=0;
  f1=0;
  Q1=0;
  D1=0;
  A1=0;
  Anew=0;
  A1inn=0;
  A1pinn=0;
  A1imm=0;
  A1pimm=0;
  Pi1=0;
  ee1=0;
  Dividends_h=0;
  Dividends_1=0;
  Dividends_2=0;
  Dividends_b=0;
  Dividends_e=0;
  Bonuses_1=0;
  Bonuses_2=0;
  Bonuses_b=0;
  Bonuses_e=0;
  EnergyPayments_1=0;
  InterestDeposits_1=0;
  Taxes_1=0;
  Taxes_CO2_1=0;
  BankingSupplier_1=0;
  c1=0;
  A1p=0;
  RDin=0;
  RDim=0;
  Inn=0;
  Imm=0;
  A1f=0;
  A1pf=0;
  Td=0;
  D1_en=0;
  A1_en=0;
  A1_ef=0;
  A1p_en=0;
  A1p_ef=0;
  EE_inn=0;
  EEp_inn=0;
  EE_imm=0;
  EEp_imm=0;
  EF_inn=0;
  EFp_inn=0;
  EF_imm=0;
  EFp_imm=0;
  Emiss1=0;
  Injection_1=0;
  Balances_1=0;
  baddebt_1=0;
  exiting_1=0;
  exiting_1_payments=0;
  c1p=0;
  shocks_machprod=0;
  shocks_techprod=0;
  shocks_machprod_max=0;
  shocks_techprod_max=0;
  shocks_rd=0;
  shocks_labprod1=0;
  shocks_labprod2=0;
  shocks_eneff1=0;
  shocks_eneff2=0;
  shocks_output1=0;
  shocks_output2=0;
  S1_temp=0;
  S2_temp=0;
  NbClient_1=0;
  NbClient_2=0;
  bonds_dem=0;
  DebtServiceToSales2_bank=0;
  DS2_rating=0;
  fB=0;
  NW_b=0;
  NW_b_c=0;
  Loans_b=0;
  BankMatch_1=0;
  BankMatch_2=0;
  BankCredit=0;
  Taxes_b=0;
  BaselBankCredit=0;
  InterestDeposits=0;
  Deposits_b=0;
  BankProfits=0;
  Dividends_b=0;
  Bailout_b=0;
  LossAbsorbed=0;
  BankEquity_temp=0;
  Bank_active=0;
  GB_b=0;
  bonds_purchased=0;
  BondRepayments_b=0;
  NL_1=0;
  NL_2=0;
  BankProfits_temp=0;
  r_deb=0;
  bankmarkup=0;
  riskWeightedAssets=0;
  capitalAdequacyRatio=0;
  buffer=0;
  Bond_share=0;
  Deposits_hb=0;
  Deposits_eb=0;
  Advances_b=0;
  Reserves_b=0;
  InterestBonds_b=0;
  LoanInterest=0;
  InterestReserves_b=0;
  InterestAdvances_b=0;
  Outflows=0;
  Inflows=0;
  DepositShare_e=0;
  DepositShare_h=0;
  baddebt_b=0;
  capitalRecovered=0;
  capitalRecovered2=0;
  capitalRecoveredShare=0;
  LossEntry_b=0;
  ReserveBalance=0;
  ShareBonds=0;
  ShareReserves=0;
  ShareAdvances=0;
  Adjustment=0;
  prior=0;
  fluxC=0;
  Con=0;
  Hon=0;
  Ton=0;
  Cax=0;
  Caxx=0;
  Cay=0;
  Cayy=0;
  Caa=0;
  fluxH=0;
  Emiss_TOT=0;
  EnvSubsidies_2_i=0;
  EnvSubsidies_2=0;

  age.resize(T);
  C_pb.resize(T);
  g_pb.resize(T);
  g_c.resize(T);
  g_c2.resize(T);
  g_c3.resize(T);
  g.resize(T);
  gtemp.resize(T);
  g_price.resize(T);
  g_secondhand.resize(T);
  g_secondhand_p.resize(T);
  age_secondhand.resize(T);
  for (std::size_t i = 0; i < T; i++)
  {
    age[i].resize(N1);
    C_pb[i].resize(N1);
    g_pb[i].resize(N1);
    g_c[i].resize(N1);
    g_c2[i].resize(N1);
    g_c3[i].resize(N1);
    g[i].resize(N1);
    gtemp[i].resize(N1);
    g_price[i].resize(N1);
    g_secondhand[i].resize(N1);
    g_secondhand_p[i].resize(N1);
    age_secondhand[i].resize(N1);
    for (std::size_t j = 0; j < N1; j++)
    {
      age[i][j].resize(N2);
      C_pb[i][j].resize(N2);
      g_pb[i][j].resize(N2);
      g_c[i][j].resize(N2);
      g_c2[i][j].resize(N2);
      g_c3[i][j].resize(N2);
      g[i][j].resize(N2);
      gtemp[i][j].resize(N2);
      g_price[i][j].resize(N2);
      g_secondhand[i][j]=0;
      g_secondhand_p[i][j]=0;
      age_secondhand[i][j]=0;
        for (std::size_t ii = 0; ii < N2; ii++){
          age[i][j][ii] = 0;
          C_pb[i][j][ii] = 0;
          g_pb[i][j][ii] = 0;
          g[i][j][ii] = 0;
          gtemp[i][j][ii] = 0;
          g_c[i][j][ii] = 0;
          g_c2[i][j][ii] = 0;
          g_c3[i][j][ii] = 0;
          g_price[i][j][ii] = 0;
        }
    }
  }

  //Initialising all other global variables to zero to avoid memory issues
  i=0;
  ii=0;
  iii=0;
  j=0;
  jjj=0;
  t=0;
  tt=0;
  rni=0;
  t0=0;
  t00=0;
  n=0;
  iterations=0;
  pareto_rv=0;
  tolerance=0;
  deviation=0;
  parber=0;
  rnd=0;
  N1r=0;
  N2r=0;
  step=0;
  stepbis=0;
  cont=0;
  age0=0;
  Amax=0;
  A1pmax=0;
  A1_en_max=0;
  A1_ef_max=0;
  A1p_en_max=0;
  A1p_ef_max=0;
  D20=0;
  DS2_min_index=0;
  newbroch=0;
  indforn=0;
  flag=0;
  payback=0;
  jmax=0;
  tmax=0;
  imax=0;
  nmachprod=0;
  nmp_temp=0;
  cmin=0;
  imin=0;
  jmin=0;
  tmin=0;
  MaxFunds=0;
  prestmax=0;
  p1prova=0;
  rated_firm_2=0;
  Qpast=0;
  Ipast=0;
  scrapmax=0;
  cmax=0;
  ind_i=0;
  ind_tt=0;
  scrap_n=0;
  sendingBank=0;
  receivingBank=0;
  c_de_min=0;
  cf_min_ge=0;
  Q_de_temp=0;
  idmin=0;
  parber_en_de=0;
  parber_en_ge=0;
  l2m=0;
  p2m=0;
  Cres=0;
  Cresbis=0;
  cpi_temp=0;
  maxbank=0;
  max_equity=0;
  multip_bailout=0;
  min_equity=0;
  ns1=0;
  ns2=0;
  mD1=0;
  mD2=0;
  multip_entry=0;
  injection=0;
  injection2=0;
  n_mach_exit=0;
  n_mach_exit2=0;
  n_mach_needed=0;
  n_mach_resid=0;
  n_mach_resid2=0;
  n_exit2=0;
  cpi_init=0;
  GDP_init=0;
  baddebt_2_temp=0;
  markdownCapital=0;
  post=0;
  prior_cb=0;
  post_cb=0;
  DepositsCheck_1=0;
  DepositsCheck_2=0;
  p2_entry=0;
  f2_exit=0;
  CurrentDemand=0;
  CompEntry_m=0;
  K_gap=0;
  K_top=0;
  loss=0;
  lossj=0;
  rani=0;
  rant=0;
  ranj=0;
  reduction=0;
  K_temp_sum=0;
  mi_en_preshock=0;
  pf_preshock=0;
  mi_en_shock=0;
  c_en_preshock=0;
  pf_shock=0;
  c_infra_t=0;
  ptemp=0;
  Ldtemp=0;
  Deposits_h=0;
  Deposits_e=0;
  GB_cb=0;
  GB=0;
  Deposits_fuel=0;
  Deposits_fuel_cb=0;
  Advances=0;
  Reserves=0;
  CapitalStock_e=0;
  NW_h=0;
  NW_gov=0;
  NW_cb=0;
  NW_e=0;
  NW_f=0;
  NW_h_c=0;
  NW_gov_c=0;
  NW_cb_c=0;
  NW_e_c=0;
  NW_f_c=0;
  NWSum=0;
  RealAssets=0;
  Wages_en=0;
  Wages=0;
  Benefits=0;
  EnergyPayments=0;
  Taxes_h=0;
  Taxes_CO2_e=0;
  Taxes_CO2=0;
  InterestDeposits_h=0;
  InterestDeposits_e=0;
  InterestBonds=0;
  InterestBonds_cb=0;
  BondRepayments_cb=0;
  Taxes_g=0;
  Consumption_h=0;
  Consumption_g=0;
  Consumption=0;
  FirmTransfers=0;
  FirmTransfers_1=0;
  FirmTransfers_2=0;
  InterestReserves=0;
  InterestAdvances=0;
  TransferCB=0;
  FuelCost=0;
  TransferFuel=0;
  Taxes_e_shock=0;
  Taxes_f_ff_shock=0;
  Taxes_e_ff_shock=0;
  Transfer_shock=0;
  Transfer_shock_f=0;
  Balance_h=0;
  Balance_1=0;
  Balance_2=0;
  Balance_e=0;
  Balance_b=0;
  Balance_g=0;
  Balance_cb=0;
  Balance_f=0;
  BalanceSum=0;
  LS=0;
  U=0;
  // Divtot_1=0;
  // Divtot_2=0;
  Cons_h=0;
  Cons = 0;
  Deposits_recovered_1=0;
  Deposits_recovered_2=0;
  A1top=0;
  A1ptop=0;
  A1_en_top=0;
  A1p_en_top=0;
  A1_ef_top=0;
  A1p_ef_top=0;
  r_depo=0;
  bonds_dem_tot=0;
  r_bonds=0;
  Bailout=0;
  G=0;
  Deficit=0;
  PSBR=0;
  NewBonds=0;
  EntryCosts=0;
  BankTransfer=0;
  r_cbreserves=0;
  r_a=0;
  r=0;
  ProfitCB=0;
  Adjustment_cb=0;
  d_cpi_target_a=0;
  inflation_a=0;
  pf=0;
  mi_en=0;
  c_en=0;
  D1_en_TOT=0;
  D2_en_TOT=0;
  D_en_TOT=0;
  K_ge=0;
  K_de=0;
  K_gelag=0;
  K_delag=0;
  Q_ge=0;
  Q_de=0;
  EI_en=0;
  EI_en_de=0;
  EI_en_ge=0;
  IC_en=0;
  LDexp_en=0;
  PC_en=0;
  c_infra=0;
  share_de=0;
  Rev_en=0;
  RD_en_de=0;
  RD_en_ge=0;
  // LDrd_de=0;
  // LDrd_ge=0;
  Inn_en_ge=0;
  Inn_en_de=0;
  A_de_inn=0;
  EM_de_inn=0;
  CF_ge_inn=0;
  ProfitEnergy=0;
  G_de_0=0;
  G_ge_0=0;
  G_ge_n_0=0;
  Tmixed=0;
  Emiss_yearly_calib=0;
  g_rate_em_y=0;
  Emiss_yearly=0;
  Emiss1_TOT=0;
  Emiss2_TOT=0;
  Emiss_en=0;
  NPP=0;
  Cum_emissions=0;
  shock_pop=0;
  shock_cons=0;
  t_CO2=0;
  t_CO2_en=0;
  Emiss_gauge=0;
  Cat=0;
  humrelease=0;
  hum=0;
  biorelease=0;
  biom=0;
  Cat1=0;
  dCat1=0;
  Con1=0;
  Ctot1=0;
  FCO2=0;
  Fin=0;
  Fout=0;
  Emiss_global=0;
  Am=0;
  Am_a=0;
  Am2=0;
  Am1=0;
  ftot=0;
  Em2=0;
  cpi=0;
  kpi=0;
  Am_en=0;
  // LD1rdtot=0;
  // LDentot=0;
  Tdtot=0;
  // LD1tot=0;
  // LD2tot=0;
  // LSe=0;
  LD=0;
  LD2=0;
  Pitot1=0;
  Pitot2=0;
  dNtot=0;
  dNmtot=0;
  ExpansionInvestment_r=0;
  ExpansionInvestment_n=0;
  ReplacementInvestment_r=0;
  ReplacementInvestment_n=0;
  Investment_r=0;
  Investment_n=0;
  Consumption_r=0;
  CreditDemand_all=0;
  CreditSupply_all=0;
  Q2tot=0;
  Q1tot=0;
  Q2dtot=0;
  D2tot=0;
  A_mi=0;
  A1_mi=0;
  A2_en_mi=0;
  A2_ef_mi=0;
  A1_en_mi=0;
  A1_ef_mi=0;
  A_sd=0;
  H1=0;
  H2=0;
  HB=0;
  GDP_r=0;
  GDP_n=0;
  d_U=0;
  d_cpi=0;
  d_Am=0;
  dw=0;
  dw2=0;
  A2scr=0;
  A1scr=0;
  Utilisation=0;
  counter_bankfailure=0;
  Loan_interest_e=0;
  DebtService_e=0;
  DeafaultedDebtRecovered_e=0;
  BadDebt_e=0;
  DebtRemittances_e=0;
  DebtWrittenOff_e=0;
  DefaultedDeposits_e=0;
  Loans_preDefault_e=0;
  EnvSubsidiesPerMachine_1_i = 0;
  Attractivess_1_i = 0;
}

void INITIALIZE(int Exseed)
{
  //Set seed
  seed=Exseed;
  //Pointer to seed
  p_seed=&seed;
  //Tolerance level used to check deviations from stock-flow consistency
  tolerance=1e-06;
  //Numbers of agents as doubles
  N1r=double(N1);
	N2r=double(N2);
  //Generate empty output files
  INTFILE();

  //Initialize energy sector
  Deposits_e=D_e0;

  Loans_e=0;

  //Determines number of firm customers of each bank
  ALLOCATEBANKCUSTOMERS();

  //Set all variables describing firm-bank networks to zero
  BankingSupplier_2=0; //Note: this initialises all values of the matrix as 0
  BankMatch_2=0;
  NbClient_2=0;
  BankingSupplier_1=0;
  BankMatch_1=0;
  NbClient_1=0;
  //Initial deposits and for firms
  Deposits_1=W10;
  Deposits_2=W20;
  Loans_2=L0;
  //Set deposits and loans from point of view of banks to zero
  Deposits_b=0;
  Loans_b=0;
  //First, iterate over all C-Firms
  for(j=1; j<=N2; j++)
  {
    //Match C-firms to banks and give those banks the initial deposits and loans
    //While firm j does not have a supplier of banking services
    while (BankingSupplier_2(j)==0)
    {
      //Randomly draw a bank
      rni=int(ran1(p_seed)*N1*N2)%NB+1;
      //If that bank has not yet reached the number of customers determined above
      if (NbClient_2(rni)< NL_2(rni))
      {
        //j becomes a customer of that bank
        BankMatch_2(j,rni)=1;
        BankingSupplier_2(j)= rni;
        NbClient_2(rni)++;
        Deposits_b(1,rni)+=Deposits_2(1,j);
        Deposits_b(2,rni)+=Deposits_2(2,j);
        Loans_b(1,rni)+=Loans_2(1,j);
        Loans_b(2,rni)+=Loans_2(2,j);
      }
      else
      {
        BankMatch_2(j,rni)=0;
      }
    }
  }

  //Iterate over all K-Firms
  for(i=1; i<=N1; i++)
  {
    //Match K-firms to banks and give those banks the initial deposits
    //While firm i does not have a supplier of banking services
    while (BankingSupplier_1(i)==0)
    {
      //Randomly draw a bank
      rni=int(ran1(p_seed)*N1*N2)%NB+1;
      //If that bank has not yet reached the number of customers determined above
      if (NbClient_1(rni)< NL_1(rni))
      {
        BankMatch_1(i,rni)=1;
        BankingSupplier_1(i)= rni;
        NbClient_1(rni)++;
        Deposits_b(1,rni)+=Deposits_1(1,i);
        Deposits_b(2,rni)+=Deposits_1(2,i);
      }
      else {
      BankMatch_1(i,rni)=0;
      }
    }
  }

  //Energy price
  c_en(2)=mi_en0;
  //Set initial prices
  p2=(1+mi2)*(w_tot_for_1_wr_mh(1)/A0+mi_en0/A0_en);
  p1=(1+mi1)*(w_tot_for_1_wr_mh(1)/(A0*pm)+mi_en0/A0_en);
  //Initialise measure of average labour productivity
  Am(2)=(A0*N2+A0*pm*N1)/(N1+N2);


  //Set initial household deposits
  Deposits_h=D_h0;

  //Initialise household classes variables
    //Deposits
  RowVector D_h0_mh=D_h0_sh_mh*D_h0;
  Deposits_mh=create_umap_1_lag(classes_mh, D_h0_mh);

  LS=0;
  w_tot_for_1_wr_mh=0;
  RowVector cl_w(2);
  RowVector NW_mh_0(2);
  for (const string& cl:classes_mh){
    //Labour supply
    LS_mh[cl]=LS0_wr*LS_sh_mh[cl]/LS_sh_mh["wr"];
    LS+=LS_mh[cl];
    //Wages with a 1-lag
    cl_w=w0_wr*w_ratios_mh[cl];
    w_mh[cl]=cl_w;
    //Total wages to pay for each worker employed (a sort of total unitary cost of labour)
    w_tot_for_1_wr_mh+=w_mh[cl](1)*ld_ratios_mh[cl];
    //Net worth
    NW_mh_0=Deposits_mh[cl];
    NW_mh[cl]=NW_mh_0;
    //Firms' ownership
    Ownership_sh_1_i_mh[cl]=Ownership_sh_1_mh0[cl];
    Ownership_sh_2_i_mh[cl]=Ownership_sh_2_mh0[cl];
    //Deposit shares
    Deposits_sh_mh[cl]=Deposits_mh[cl](1)/Deposits_h(1);
    //Taxes on dividends and bonuses for consumption function
    Taxes_div_mh[cl](2)=0;
    Taxes_bon_mh[cl](2)=0;
	}

  //Finish initialising bank balance sheets, distribute remaining stocks between banks
  for(i=1; i<=NB; i++)
  {
    //Bank's market share is equal to its share of firm customers
    fB(1,i)=(NbClient_2(i)+NbClient_1(i))/(N2+N1);
    fB(2,i)=(NbClient_2(i)+NbClient_1(i))/(N2+N1);
    //Bank's initial stock of government bonds is a % of its loan portfolio
    GB_b(1,i)=varphi*Loans_b(1,i);
    GB_b(2,i)=varphi*Loans_b(2,i);
    //Give the bank a share of household and energy sector deposits
    Deposits_b(1,i)+=fB(1,i)*(Deposits_h(1)+Deposits_e(1));
    Deposits_b(2,i)+=fB(2,i)*(Deposits_h(2)+Deposits_e(2));
    Deposits_hb(1,i)+=fB(1,i)*(Deposits_h(1));
    Deposits_hb(2,i)+=fB(2,i)*(Deposits_h(2));
    Deposits_eb(1,i)+=fB(1,i)*(Deposits_e(1));
    Deposits_eb(2,i)+=fB(2,i)*(Deposits_e(2));
    //Initialise share of bank in aggregate household and energy sector deposits
    DepositShare_e(i)=Deposits_eb(1,i)/(Deposits_e(1));
    DepositShare_h(i)=Deposits_hb(1,i)/(Deposits_h(1));
    //Bank reserves are given as a residual when setting an exogenous initial value for bank net worth (NW_b0)
    //and assuming that initial central bank advances are zero
    Reserves_b(1,i)=Deposits_b(1,i)+NW_b0*fB(1,i)-GB_b(1,i)-Loans_b(1,i);
    Reserves_b(2,i)=Deposits_b(2,i)+NW_b0*fB(2,i)-GB_b(2,i)-Loans_b(2,i);
    NW_b(1,i)=Loans_b(1,i)+Reserves_b(1,i)+GB_b(1,i)-Deposits_b(1,i)-Advances_b(1,i);
    NW_b(2,i)=Loans_b(2,i)+Reserves_b(2,i)+GB_b(2,i)-Deposits_b(2,i)-Advances_b(2,i);
    riskWeightedAssets(i)=riskWeightLoans*Loans_b(1,i)+riskWeightGovBonds*GB_b(1,i);
    capitalAdequacyRatio(i)=NW_b(1,i)/riskWeightedAssets(i);
  }
  //Governmen bonds held by central bank are given by identity from the above
  GB_cb(1)=NW_b0+W10*N1+W20*N2+D_h0+Deposits_e(1)-GB_b.Row(1).Sum()-Loans_b.Row(1).Sum();
  GB_cb(2)=NW_b0+W10*N1+W20*N2+D_h0+Deposits_e(2)-GB_b.Row(2).Sum()-Loans_b.Row(2).Sum();
  //Initial stock of overall government bonds (banks+CB)
  GB(1)=GB_cb(1)+GB_b.Row(1).Sum();
  GB(2)=GB_cb(2)+GB_b.Row(2).Sum();
  //CB advances assumed to be zero initially
  Advances_b=0;
  Advances=0;
  //CB Stock of gov. bonds implies an equal aggregate stock of reserves
  Reserves=GB_cb;
  //Set nominal value of initial capital stock
  CapitalStock=K0/dim_mach*(1+mi1)*(w_tot_for_1_wr_mh(1)/(A0*pm)+c_en(2)/A0_en);
  deltaCapitalStock=0;
  //Initialise net worth of households, K-Firms, government and central bank
  NW_h=Deposits_h;
  NW_1=Deposits_1;
  NW_gov=-GB;
  NW_cb=GB_cb-Reserves;
  //Initial central bank rate is set by converting annual r_base to quarterly rate
  r=pow((1+r_base),0.25)-1;
  //Convert quarterly inflation target to annual one
  d_cpi_target_a=pow((1+d_cpi_target),4)-1;

  //Initialise deposit, CB reserve and gov. bond interest rates
  if(flag_rate_setting_markup==0 || flag_rate_setting_markup==2)
  {
    //additive mark-up case
    r_depo=0;       // For now, interest rate set to zero, change that to "r_depo=r-bankmarkdown;"
                    // when banks compete on the deposit market
    r_cbreserves=0; // For now, interest rate set to zero, change that to "r_cbreserves=r-centralbankmarkdown;"
                    // later
    if (flag_rate_setting_markup==0){
			r_bonds=r;      // For now, interest rate set to zero, change that to "r_bonds=r-bondsmarkdown;"
                    		// when the bond market is fixed
    }
    else if(flag_rate_setting_markup==2){
      r_bonds=bondsmarkdown;
    }
  }
  else if(flag_rate_setting_markup==1)
  {
    //Multiplicative mark-up case
    r_depo=r*(1-bankmarkdown);
    r_cbreserves=r*(1-centralbankmarkdown);
    r_bonds=r*(1-bondsmarkdown);
  }


  //Initial central bank profits
  ProfitCB(2)=r_cbreserves*Reserves(2);
  //Initial wages energy sector
  Wages_en=0;
  for (const string& cl:classes_mh){
    Wages_en_mh[cl]=0;
  }
  //Initial change in mean productivity
  d_Am=d_Am_init;
  //Initial wage gross growth rate
  dw2=1;
  //Initial passthrough rates (for energy price shock)
  pass_1=1;
  pass_2=1;

  //Set initial value for carbon tax
  if(flag_tax_CO2>0 && flag_tax_CO2<5)
  {
    t_CO2_en=t_CO2_en_0;
  }

  //Initial productivities, energy efficiencies, environmental friendliness
	A1=A0;
	A1p=A0*pm;
	A2=A0;
  A2_mprod=A0;
  A=A0;

  A1_en=A0_en;
  A1p_en=A0_en*pm;
  A2_en=A0_en;
  A_en=A0_en;

  A1_ef=A0_ef;
  A1p_ef=A0_ef*A_ef_ratio0;
  A2_ef=A0_ef;
  A_ef=A0_ef;

  //Set initial values for energy sector
  //Brown energy tech.
  A_de=A_de0;
  EM_de=EM0;
  //Fossil fuel price
  pf=pf0;
  //Mark-up
  mi_en=mi_en0;
  //Green energy tech
  CF_ge(1)=CF_ge0;

  //Initialise climate module
  if (flag_nonCO2_force==0)
  {
    //Case of no non-CO2 forcing
    //Initial atmospheric carbon
    Cat(1)=Catinit0;
    //Initial biosphere carbon
    biom(1)=biominit0;
    //Initial humus carbon
    hum(1)=huminit0;
    //Initial temperature anomaly
    Tmixed(1)=Tmixedinit0;
    //Initial carbon and heat content of ocean layers
    for (j=1;  j<=ndep; j++)
    {
      Con(1,j)=Coninit0(j);
      Hon(1,j)=Honinit0(j);
    }
  }
  else
  {
    //Case of non-CO2 forcing
    //Initial atmospheric carbon
    Cat(1)=Catinit1;
    //Initial biosphere carbon
    biom(1)=biominit1;
    //Initial humus carbon
    hum(1)=huminit1;
    //Initial temperature anomaly
    Tmixed(1)=Tmixedinit1;
    //Initial carbon and heat content of ocean layers
    for (j=1;  j<=ndep; j++)
    {
      Con(1,j)=Coninit1(j);
      Hon(1,j)=Honinit1(j);
    }
  }

  if(flag_cum_emissions==1)
  {
    //If using cumulative emissions module, initialise temperature anomaly
    Tmixed(1)=T_0_cumemiss;
  }
  //Initialise cumulative emissions, yearly emissions, etc.
  Cum_emissions=Cum_emissions_0;
  Emiss_yearly_calib=Emiss_yearly_0;
  g_rate_em_y=0;
  Emiss_TOT=0;
  Emiss_yearly=0;

  //Consumer price index
	cpi(2)=(1+mi2)*(w_tot_for_1_wr_mh(1)/A0+c_en(2)/A0_en);
  cpi(3)=(1+mi2)*(w_tot_for_1_wr_mh(1)/A0+c_en(2)/A0_en);
  cpi(4)=(1+mi2)*(w_tot_for_1_wr_mh(1)/A0+c_en(2)/A0_en);
  cpi(5)=(1+mi2)*(w_tot_for_1_wr_mh(1)/A0+c_en(2)/A0_en);

  //Initialise cost vectors
	C=w_tot_for_1_wr_mh(1)/A0+c_en(2)/A0_en;
	c2=w_tot_for_1_wr_mh(1)/A0+c_en(2)/A0_en;
	c1=w_tot_for_1_wr_mh(1)/(A0*pm)+c_en(2)/(A0_en*pm);
  //Competitiveness and ability to serve demand
	Em2=1;
	l2=1;
  E2=1;
	f1=1/N1r;
	f2=1/N2r;
  //Set initial capital stock in terms of productive capacity
	K=K0;
  //t0 needed when iterating over technology arrays
	t0=1;

  for (i=1; i<=N2; i++)
	{
    FirmDefaultProbability(i)=FirmDefaultProbability_init;
  }

  Bank_active=1;
  //Set initial baseline lending rate; either additive or multiplicative mark-up
  //over CB rate
  if(flag_rate_setting_markup==0 || flag_rate_setting_markup==2)
  {
    bankmarkup=bankmarkup_init;
    r_deb=r+bankmarkup;
  }

  if(flag_rate_setting_markup==1)
  {
    //Adjust initial value for multiplicative case
    bankmarkup_init*=100;
    bankmarkup=bankmarkup_init;
    r_deb=r*(1+bankmarkup);
  }

  //initialise debt service
  DebtService_2=L0*(r+bankmarkup_init)+repayment_share*L0;

	EI=0;
  //Initial investment is given by amount necessary to replace average amount of machines reaching max. age
	I=ROUND((((K0)/(agemax+1)))/dim_mach)*dim_mach;
	mu2=mi2;
	Td.element(0)=0;
	step=N2/N1;
	cont=0;
	Match=0;

  //Match C-firms to K-firms & distribute initial stocks of machines
	for (i=1; i<=N1; i++)
	{
		cont+=step;
		for (j=0; j<step; j++)
		{
			Match(cont-j,i)=1;
			fornit(cont-j)=i;
		}
	}

	for (j=1; j<=N2; j++)
	{
		n_mach(j)=K(j)/dim_mach;
		while (n_mach(j) > 0)
		{
			i++;
			if (i > N1)
      {
				i=1;
      }
      //Initial capital stock is assumed to not come from initial supplier
      if (fornit(j) != i)
			{
				//Random age
        age0=int(ran1(p_seed)*(agemax+1))%int((agemax+1))+1;
        g[0][i-1][j-1]++;
        g_price[0][i-1][j-1]=p1(i);
				gtemp[0][i-1][j-1]++;
				g_c[0][i-1][j-1]++;
        g_c2[0][i-1][j-1]++;
        g_c3[0][i-1][j-1]++;
				age[0][i-1][j-1]=age0;
        n_mach(j)--;
			}
		}
	}

  //Set initial demand for consumption goods and energy to a level consistent with unemployment=ustar
  Dtot0 = w_tot_for_1_wr_mh(1)*(LS0_wr*(1-ustar)-I.Sum()/dim_mach/(A0*pm));  //Initial total demand for C goods as wage costs * (tot employment - exmployment in K sector) = wages to employed in C sector
  Den0 = 0;
  for (const string& cl:classes_mh){
    // Split initial expenditure to calculate initial energy expenditure
    Expenditure_tot_mh[cl] = Dtot0 * (w_mh[cl](1) * ld_ratios_mh[cl])/w_tot_for_1_wr_mh(1); // Based on contribution to total wages for each worker
    // Need to calcualte energy expenditure to subtract to total and get initial C-goods demand
    Den0_mh[cl] = energy_expenditure_sh_mh[cl] * Expenditure_tot_mh[cl];
    Den0 += Den0_mh[cl];
  }
  D20 = Dtot0 - Den0;
	D2=D20/N2r;
  De=D2.Row(1);
  N=omicron*D20/N2r;
  //Nominal value of C-firm inventories
  Inventories=(omicron*D20/N2r)*p2(1);

  NW_2=Deposits_2+CapitalStock+Inventories-Loans_2;
  //Assume initial firm sales are homogeneous
	S1=(((I(1)/dim_mach)*N2r)/N1r)*p1;
  S1_pre=(((I(1)/dim_mach)*N2r)/N1r)*p1;
  S1_post=(((I(1)/dim_mach)*N2r)/N1r)*p1;
	S2=D20/N2r*p2(1);
  //Rough measure of initial net revenue
  mol=D20/N2r*p2(1)-D20/N2r/A0*w_tot_for_1_wr_mh(1)-D20/N2r/A0_en*c_en(2);
  //Initial dividends
  Dividends_h(2)=0; //(mol(1)-L0*(r+bankmarkup_init))*N2r*d2+(S1.Sum()-(((I.element(1)/dim_mach)*N2r))*c1(1))*d1+db*(r_bonds*GB_b.Row(1).Sum()+L0*(r+bankmarkup_init)*N2r)+de*c_en(2)*(De.Sum()/A0_en+((I(1)/dim_mach)*N2r)/(A0_en*pm));
  for (const string& cl:classes_mh){
    //Initial ownership
    Ownership_1_i_mh[cl]=SP(Ownership_sh_1_i_mh[cl], NW_1.Row(2));
    Ownership_2_i_mh[cl]=SP(Ownership_sh_2_i_mh[cl], NW_2.Row(2));
    //Initial dividends
    RowVector Dividends_cl_0(2);
    Dividends_cl_0(1) = 0;
    Dividends_cl_0(2) = Ownership_sh_2_mh0[cl]*Dividends_h(2); //Assumed based on C-firms initial shares of ownership, since most part comes from C sector
    Dividends_mh[cl] = Dividends_cl_0;
    //Initial bonuses
    Bonuses_mh[cl]=0;
  }
  Bonuses_h=0;
  Ownership_1_mh=umap_sum_each_key(Ownership_1_i_mh);
  Ownership_2_mh=umap_sum_each_key(Ownership_2_i_mh);

  U(2)=(LS0_wr-((D20)/A0+(((I.element(1)/dim_mach)*N2r))/(A0*pm)))/LS0_wr;

  //Give the energy sector an initial stock of productive capacity
  G_ge(1)=K_ge0_perc*(1.1*(D20/A0_en+((I(1)/dim_mach)*N2r)/(A0_en*pm)));
  G_ge_n(1)=CF_ge(1)*G_ge(1);
  G_de(1)=(1.1*(D20/A0_en+((I(1)/dim_mach)*N2r)/(A0_en*pm)))-G_ge(1);
  C_de(1)=pf/A_de(1)+t_CO2_en*EM_de(1);
  K_gelag=G_ge(1);
  K_delag=G_de(1);
  G_de_0=G_de(1);
  G_ge_0=G_ge(1);
  G_ge_n_0=G_ge_n(1);
  CapitalStock_e=G_ge_n_0;
  NW_e=Deposits_e+CapitalStock_e;

  //Initialise technological change
  A1top=A0;
  A1ptop=A0*pm;
  A1f=A0;
  A1pf=A0*pm;
  A1_en_top=A0_en;
  A1p_en_top=A0_en;
  A1_ef_top=A0_ef;
  A1p_ef_top=A0_ef;

  RD.Row(1)=nu*S1;
  RD.Row(2)=nu*S1;
  t=0;

  //Energy sector initalisation in ENERGY_INV_PROD in energy module
  TECHANGEND();
}

void REGIMESHIFTS(void)
{
  if (flag_change_income_shares==1){
    double w_tot_for_1_wr_mh_check=0;
    double Wages_tot_check = 0;
    for (const string& cl : classes_mh){
      //Update wage and wage ratio for each class
      double new_tot_wage_cl = Income_gross_h(2) * income_sh_regime_shift_mh[cl] - Dividends_mh[cl](2) - Benefits_mh[cl] - InterestDeposits_mh[cl];
      w_mh[cl](1) = new_tot_wage_cl / (LS_mh[cl] * (1-U_mh[cl]));
      w_ratios_mh[cl] = w_mh[cl](1) / w_mh["wr"](1);
      //The one above is the wage based on previous timestep income, before wage update --> update
      w_mh[cl](1) *= (1+dw);
      w_mh[cl](2) = w_mh[cl](1);

      //For checks
      w_tot_for_1_wr_mh_check += (w_mh[cl](1) * ld_ratios_mh[cl]);
      Wages_tot_check += w_mh[cl](1) * (LS_mh[cl] * (1-U_mh[cl]));
    }

    // Check that total wages remain the same
    double Wages_next_timestep = w_tot_for_1_wr_mh(1) * LS_mh["wr"] * (1-U_mh["wr"]);
    deviation = fabs((Wages_tot_check - Wages_next_timestep) / Wages_next_timestep);
    if(deviation>0.001)
    {
      std::cerr<<"Error: total wages changed by " << deviation*100 << "% during income shares change in period " << t << endl;
      Errors << "\n Total wages changed by " << deviation*100 << "% during income shares change in period " << t << endl;
    }
    // Check that wage to be paid for each worker remains the same
    deviation = fabs((w_tot_for_1_wr_mh_check - w_tot_for_1_wr_mh(1)) / w_tot_for_1_wr_mh(1));
    if(deviation>0.001)
    {
      std::cerr<<"Error: total wages to be paid for each worker changed by " << deviation*100 << "% during income shares change in period " << t << endl;
      Errors << "\n Total wages to be paid for each worker changed by " << deviation*100 << "% during income shares change in period " << t << endl;
    }
  } else if(flag_change_income_shares==2){
    w_tot_for_1_wr_mh(1)=0;
    for (const string& cl : classes_mh){
      w_ratios_mh[cl] = income_sh_regime_shift_mh[cl];
      w_mh[cl](1) = w_ratios_mh[cl] * w_mh["wr"](1);
      w_tot_for_1_wr_mh(1) += (w_mh[cl](1) * ld_ratios_mh[cl]);
    }
    w_tot_for_1_wr_mh(2)=w_tot_for_1_wr_mh(1);
  }

  if (flag_change_apc==1){
    a1_mh  = a1_mh_regime_shift;
    a2_mh  = a2_mh_regime_shift;
    a3_mh  = a3_mh_regime_shift;
  }

  if (flag_change_en_exp_shares==1){
    energy_expenditure_sh_mh = en_exp_sh_mh_regime_shift;
  }

  // Exogenous change in C-firms markup
  if (flag_c_price_increase==1){
    mu2.Row(1) = c_price_increase_perc;
    mu2.Row(2) = mu2.Row(1);
  } else if (flag_c_price_increase==2){
    mu2.Row(1)*=(1+c_price_increase_perc);
    mu2.Row(2) = mu2.Row(1);
  }

  //Exogenous change in target public spending to GDP
  if (flag_change_public_spending==1){
    g_cons_GDP_ratio = g_cons_GDP_ratio_regime_shift;
  }else if (flag_change_public_spending==2 && (flag_change_wage_tax_rate!=0 || flag_change_dividends_tax_rate!=0 || flag_change_wealth_tax_rate!=0)){
    //Save tax rates pre-change
    for (const string& cl : classes_mh){
      aliqw_mh_prechange[cl] = aliqw_mh[cl];
      aliqwealth_mh_prechange[cl] = aliqwealth_mh[cl];
    }
    aliqdiv_prechange = aliqdiv;
  }

  
  std::map<string, double> labour_income_tax_sh_mh;

  //Exogenous change in wage tax rate
  if (flag_change_wage_tax_rate==1){ //With set parameter
    for (const string& cl : classes_mh){
      aliqw_mh[cl] = aliqw_regime_shift_mh[cl];
    }
  } else if (flag_change_wage_tax_rate==2){ //Keep same wage tax revenues
    double denominator = 0;
    for (const string& cl : classes_mh){
      denominator += Wages_mh[cl] * aliqw_regime_shift_mh[cl];
    }
    aliqw_mh["wr"] = Taxes_w_h / denominator;
    for (const string& cl : classes_mh){
      aliqw_mh[cl] = aliqw_mh["wr"] * aliqw_regime_shift_mh[cl];
    }
    // Check that total taxes on wages remain the same
    double Taxes_w_h_new = 0;
    for (const string& cl : classes_mh){
      Taxes_w_h_new += Wages_mh[cl] * aliqw_mh[cl];
    }
    deviation = fabs((Taxes_w_h_new - Taxes_w_h) / Taxes_w_h);
    if(deviation>tolerance)
    {
      std::cerr<<"Error: total taxes on wages changed by " << deviation*100 << "% during tax rate change in period " << t << endl;
      Errors << "\n Total taxes on wages changed by " << deviation*100 << "% during tax ratechange in period " << t << endl;
    }
  } else if (flag_change_wage_tax_rate==3){ //Keep same tax revenues on households but change also tax rate on dividends to aliqdiv_regime_shift
    aliqdiv = aliqdiv_regime_shift;
    double Taxes_income_h = Taxes_w_h;
    double denominator = 0;
    for (const string& cl : classes_mh){
      denominator += Wages_mh[cl] * aliqw_regime_shift_mh[cl];
    }
    aliqw_mh["wr"] = (Taxes_income_h - Dividends_h(1) * aliqdiv) / denominator;
    for (const string& cl : classes_mh){
      aliqw_mh[cl] = aliqw_mh["wr"] * aliqw_regime_shift_mh[cl];
    }
    // Check that total taxes on households remain the same
    double Taxes_w_h_new = 0;
    for (const string& cl : classes_mh){
      Taxes_w_h_new += Wages_mh[cl] * aliqw_mh[cl];
    }
    double Taxes_tot_new = Taxes_w_h_new + Dividends_h(1) * aliqdiv;
    deviation = fabs((Taxes_tot_new - Taxes_income_h) / Taxes_tot_new);
    if(deviation>tolerance)
    {
      std::cerr<<"Error: total taxes on households changed by " << deviation*100 << "% during tax rate change in period " << t << endl;
      Errors << "\n Total taxes on households changed by " << deviation*100 << "% during tax rate change in period " << t << endl;
    }
  } else if (flag_change_wage_tax_rate==4){ //Change tax rates based only on one parameter, aliqw_progressivity_regime_shift 
    std::map<string, double> labour_income_share_mh;
    std::map<string, double> correction_factors_mh;
    double tot_correction_factors = 0;
    double average_tax_rate_labour_income = 0;
    for (const string& cl : classes_mh){
      //Calculate labour income tax shares
      labour_income_share_mh[cl] += Wages_mh[cl] / Wages;
      //Calculate correction factors
      correction_factors_mh[cl] = labour_income_share_mh[cl] * pow((labour_income_share_mh[cl] / (LS_mh[cl]/LS)), aliqw_progressivity_regime_shift);
      tot_correction_factors += correction_factors_mh[cl];
      //Calculate average tax rate for following calculation
      average_tax_rate_labour_income += labour_income_share_mh[cl] * aliqw_mh[cl];
    }

    for (const string& cl : classes_mh){
      //Normalise correction factors to obtain target shares of labour income taxes paid by each class
      labour_income_tax_sh_mh[cl] = correction_factors_mh[cl] / tot_correction_factors;
      //Calculate new tax rate for each class
      aliqw_mh[cl] = average_tax_rate_labour_income * labour_income_tax_sh_mh[cl] / labour_income_share_mh[cl];
    }
    
    //Check that total taxes on wages remain the same
    double Taxes_w_h_new = 0;
    for (const string& cl : classes_mh){
      Taxes_w_h_new += Wages_mh[cl] * aliqw_mh[cl];
    }
    deviation = fabs((Taxes_w_h_new - Taxes_w_h) / Taxes_w_h);
    if(deviation>tolerance)
    {
      std::cerr<<"Error: total taxes on wages changed by " << deviation*100 << "% during tax rate change in period " << t << endl;
      Errors << "\n Total taxes on wages changed by " << deviation*100 << "% during tax ratechange in period " << t << endl;
    }
  }

  //Exogenous change in dividends tax rate
  if (flag_change_dividends_tax_rate==1){
    aliqdiv = aliqdiv_regime_shift;
  } else if(flag_change_dividends_tax_rate==2){
    // Get previous shares of personal income taxation between capital and labour income
    double Taxes_income_h = Taxes_w_h + Taxes_div_h;
    double Taxes_w_sh_inc_h = Taxes_w_h / Taxes_income_h;
    double Taxes_div_sh_inc_h = Taxes_div_h / Taxes_income_h;
    // Define new shares based on policy parameter
    Taxes_w_sh_inc_h -= aliqdiv_regime_shift;
    Taxes_div_sh_inc_h += aliqdiv_regime_shift;
    // Define new taxation from capital and labour income 
    double New_Taxes_div_h = Taxes_div_sh_inc_h * Taxes_income_h;
    double New_Taxes_w_h = Taxes_w_sh_inc_h * Taxes_income_h;
    // Define new dividends tax rate
    aliqdiv = New_Taxes_div_h / Dividends_h(1);
    // Define new labour income tax rates
    std::map<string, double> New_Taxes_w_mh;
    for (const string& cl : classes_mh){
      if (flag_change_wage_tax_rate!=4){ //If we are also redistributing labour income taxation, keep new tax shares calculated above
        //Calculate labour income tax shares
        labour_income_tax_sh_mh[cl] = Taxes_w_mh[cl] / Taxes_w_h;
      }
      // Define new labour income taxation to collect from each class
      New_Taxes_w_mh[cl] = labour_income_tax_sh_mh[cl] * New_Taxes_w_h;
      // Define new labour income tax rates
      aliqw_mh[cl] = New_Taxes_w_mh[cl] / Wages_mh[cl];
    }
    // Check that total taxes on households remain the same
    double Taxes_w_h_new = 0;
    for (const string& cl : classes_mh){
      Taxes_w_h_new += Wages_mh[cl] * aliqw_mh[cl];
    }
    double Taxes_tot_new = Taxes_w_h_new + Dividends_h(1) * aliqdiv;
    deviation = fabs((Taxes_tot_new - Taxes_income_h) / Taxes_tot_new);
    if(deviation>tolerance)
    {
      std::cerr<<"Error: total taxes on households changed by " << deviation*100 << "% during tax rate change in period " << t << endl;
      Errors << "\n Total taxes on households changed by " << deviation*100 << "% during tax rate change in period " << t << endl;
    }
    // Print if getting negative values of tax rates
    if(aliqdiv<0)
    {
      std::cerr<<"Error: Negative tax rate on dividends after shifting taxation from capital to labour income";
      Errors << "\n Negative tax rate on dividends after shifting taxation from capital to labour income";      
    }
    if(aliqdiv>1)
    {
      std::cerr<<"Error: Tax rate on dividends higher than 1 after shifting taxation from labour to capital income";
      Errors << "\n Tax rate on dividends higher than 1 after shifting taxation from labour to capital income";      
    }
    for (const string& cl : classes_mh){
      if(aliqw_mh[cl]<0)
      {
        std::cerr<<"Error: Negative tax rate on labour income of class " << cl << " after shifting taxation from labour to capital income";
        Errors << "\n Negative tax rate on labour income of class " << cl << " after shifting taxation from labour to capital income";      
      }
    }
  }

  if (flag_highly_progressive_taxation==1){
    aliqw_mh["ma"] += aliq_households_increase;
    if(aliqw_mh["ma"]<0)
    {
      std::cerr<<"Error: Negative tax rate on labour income of Managers after decrease at " << aliqw_mh["ma"];
      Errors << "\n Negative tax rate on labour income of Managers after decrease at " << aliqw_mh["ma"];      
    }
    if(aliqw_mh["ma"]>1)
    {
      std::cerr<<"Error: Tax rate on labour income of Managers higher than 1 after increase at " << aliqw_mh["ma"];
      Errors << "\n Negative tax rate on labour income of Managers higher than 1 after increase at " << aliqw_mh["ma"];      
    }
  }else if(flag_highly_progressive_taxation==2){
    aliqdiv += aliq_households_increase;
    if(aliqdiv<0)
    {
      std::cerr<<"Error: Negative tax rate on dividends after decrease at " << aliqw_mh["ma"];
      Errors << "\n Negative tax rate on dividends after decrease at " << aliqw_mh["ma"];      
    }
    if(aliqdiv>1)
    {
      std::cerr<<"Error: Tax rate on dividends higher than 1 after increase at " << aliqw_mh["ma"];
      Errors << "\n Negative tax rate on dividends higher than 1 after increase at " << aliqw_mh["ma"];      
    }
  }else if (flag_highly_progressive_taxation==3){
    aliqw_mh["wr"] += aliq_households_increase;
    if(aliqw_mh["wr"]<0)
    {
      std::cerr<<"Error: Negative tax rate on labour income of Workers after decrease at " << aliqw_mh["wr"];
      Errors << "\n Negative tax rate on labour income of Workers after decrease at " << aliqw_mh["wr"];      
    }
    if(aliqw_mh["wr"]>1)
    {
      std::cerr<<"Error: Tax rate on labour income of Workers higher than 1 after increase at " << aliqw_mh["wr"];
      Errors << "\n Negative tax rate on labour income of Workers higher than 1 after increase at " << aliqw_mh["wr"];      
    }
  }

  //Exogenous change in wealth tax rate
  if (flag_change_wealth_tax_rate==1){
    for (const string& cl : classes_mh){
      aliqwealth_mh[cl] = aliqwealth_regime_shift_mh[cl];
    }
  }

  //Carbon tax introduction
  if (flag_tax_CO2==5 || flag_tax_CO2==6){
    t_CO2=t_CO2_0;
    t_CO2_en_0=t_CO2_0;
    t_CO2_en=t_CO2_en_0;
    cpi_t_regime_shifts=cpi(2);
  }
}

void SETVARS(void)
{
  for (j=1; j<=N2; ++j)
  {
    Loans_2(1,j)=Loans_2(2,j);
    Deposits_2(1,j)=Deposits_2(2,j);
    deltaCapitalStock(1,j)=0;
    InterestDeposits_2(j)=0;
    CapitalStock(1,j)=CapitalStock(2,j);
    c2(j)=0;
    c2p(j)=0;
		A2(j)=0;
    A2_mprod(j)=0;
    A2e(j)=0;
    A2e_en(j)=0;
    A2e_ef(j)=0;
    A2e2(j)=0;
    A2e_en2(j)=0;
    A2e_ef2(j)=0;
	  c2e(j)=0;
    EI(1,j)=0;
		SI(j)=0;
		I(j)=0;
    Dividends_2_i(j)=0;
    Bonuses_2_i(j)=0;
    for (const string& cl:classes_mh){
      Dividends_2_i_mh[cl](j)=0;
    }
    Taxes_2(j)=0;
    Injection_2(j)=0;
    DebtRemittances2(j)=0;
    A2_en(j)=0;
    A2_ef(j)=0;
    baddebt_2(j)=0;
    S2(1,j)=0;
    S2_temp(1,j)=0;
    D2(1,j)=0;
    Q2(j)=0;
    k(j)=0;
    EId(j)=0;
    SId(j)=0;
    exiting_2(j)=0;
    exit_payments2(j)=0;
    exit_equity2(j)=0;
    exit_marketshare2(j)=0;
    n_mach_entry(j)=0;
    Loss_Capital(j)=0;
    Loss_Capital_mol(j)=0;
    Loss_Inventories(j)=0;
    k_entry(j)=0;
    marker_age(j)=0;
    K_loss(j)=0;
    C_loss=0;
    I_loss=0;
  }

  for (i=1; i<=N1; i++)
  {
    Deposits_1(1,i)=Deposits_1(2,i);
    InterestDeposits_1(i)=0;
    Q1(i)=0;
	  D1(i)=0;
    S1(i)=0;
    S1_pre(i)=0;
    S1_post(i)=0;
    S1_temp(1,i)=0;
    Dividends_1_i(i)=0;
    Bonuses_1_i(i)=0;
    for (const string& cl:classes_mh){
      Dividends_1_i_mh[cl](i)=0;
    }
    Taxes_1(i)=0;
    Injection_1(i)=0;
    baddebt_1(i)=0;
    exiting_1(i)=0;
    exiting_1_payments(i)=0;
  }

  for (i=1; i<=NB; i++)
  {
    Loans_b(1,i)=Loans_b(2,i);
    Deposits_b(1,i)=Deposits_b(2,i);
    Deposits_hb(1,i)=Deposits_hb(2,i);
    Deposits_eb(1,i)=Deposits_eb(2,i);
    GB_b(1,i)=GB_b(2,i);
    Advances_b(1,i)=Advances_b(2,i);
    Reserves_b(1,i)=Reserves_b(2,i);
    InterestDeposits(i)=0;
    baddebt_b(i)=0;
    Outflows(i)=0;
    Inflows(i)=0;
    NW_b(1,i)=NW_b(2,i);
    Dividends_b_i(i)=0;
    Bonuses_b_i(i)=0;
    Taxes_b(i)=0;
    LoanInterest(i)=0;
    bonds_dem(i)=0;
    Bailout_b(i)=0;
    LossAbsorbed(i)=0;
    capitalRecovered(i)=0;
    capitalRecovered2(i)=0;
    capitalRecoveredShare(i)=0;
    LossEntry_b(i)=0;
  }

  Taxes_g=0;
  Taxes_CO2(1)=0;
  Loan_interest_e=0;
  DebtService_e=0;
  DeafaultedDebtRecovered_e=0;
  BadDebt_e=0;
  DebtRemittances_e=0;
  DebtWrittenOff_e=0;
  DefaultedDeposits_e=0;
  Loans_preDefault_e=0;
  Taxes_e_shock=0;
  Taxes_f_ff_shock=0;
  Taxes_e_ff_shock=0;
  Transfer_shock_f=0;
  Wages=0;
  Deposits_recovered_1=0;
  Deposits_recovered_2=0;
  EntryCosts=0;
  BankTransfer=0;
  EnvSubsidies_2=0;
  EnvSubsidies_2_i=0;
  Deposits_h(1)=Deposits_h(2); //I think useless, since opposite set is done in UPDATE at the end of the period
  Deposits_e(1)=Deposits_e(2);
  GB_cb(1)=GB_cb(2);
  GB(1)=GB(2);
  Advances(1)=Advances(2);
  Reserves(1)=Reserves(2);
  Wages=0;
  Dividends_h(1)=0;
  Dividends_1=0;
  Dividends_2=0;
  Dividends_b=0;
  Bonuses_h=0;
  Bonuses_1=0;
  Bonuses_2=0;
  Bonuses_b=0;
  for (const string& cl:classes_mh){
    Wages_mh[cl]=0;
    Deposits_recovered_1_mh[cl]=0;
    Deposits_recovered_2_mh[cl]=0;
    Deposits_mh[cl](1)=Deposits_mh[cl](2);
    Deposits_sh_mh[cl]=Deposits_mh[cl](1)/Deposits_h(1);
    Dividends_mh[cl](1)=0;
    Dividends_1_mh[cl]=0;
    Dividends_2_mh[cl]=0;
    Dividends_b_mh[cl]=0;
    Dividends_e_mh[cl]=0; //Even if not necessary since assigned (not summed)
    Bonuses_mh[cl]=0;
    Bonuses_1_mh[cl]=0;
    Bonuses_2_mh[cl]=0;
    Bonuses_b_mh[cl]=0;
    Bonuses_e_mh[cl]=0; //Even if not necessary since assigned (not summed)
  }
  //Set contribution to entry of new C-/K- firms based on deposit shares
  Entry_financing_sh_mh["wr"] = 0; //Households never finance entry of new firms
  /* if (Entry_financing_sh_mh["wr"]<tolerance){    //VERSION WITH WORKERS FINANCING IF THEY CAN
    Entry_financing_sh_mh["wr"]=0;
  }
  else{
    cout << "Households will contribute to entry of C-/K-firms since their deposits share is " << Entry_financing_sh_mh["wr"] << endl;
  } */
  Entry_financing_sh_mh["pr"] = Deposits_sh_mh["pr"];
  Entry_financing_sh_mh["ma"] = 1 - Entry_financing_sh_mh["pr"]; //Managers finance all the rest || NOT A PROBLEM UNTIL WORKERS DO NOT HAVE A LOT OF DEPOSITS

  EnergyPayments=0;
  InterestReserves=0;
  InterestAdvances=0;
  FirmTransfers=0;
  FirmTransfers_1=0;
  FirmTransfers_2=0;

  Pitot1=0;
	Pitot2=0;
  /* LD2tot_wr=0;
	LD1tot_wr=0; */
  dNtot=0;
	dNmtot=0;
  DebtServiceToSales2_bank=0;
  Bailout=0;
  NewBonds=0;
  cpi(1)=0;
  kpi=0;
  Em2(1)=0;
  ns1=0;
  ns2=0;
  mD1=0;
  mD2=0;
  n_mach_exit=0;
  n_mach_needed=0;
  Emiss_en=0;
  Am(1)=0;
  Am2=0;
  Am1=0;
  Am_a=0;
  Am_en(1)=0;
  ftot=0;
  Consumption_r=0;
  A_mi=0;
  A1_mi=0;
  A2_en_mi=0;
  A2_ef_mi=0;
  A1_en_mi=0;
  A1_ef_mi=0;
	A_sd=0;
  H1=0;
	H2=0;
	HB=0;
  PC_en=0;
  IC_en=0;
  p2_entry=0;
  f2_exit=0;
  EntryShare=0;
  CompEntry=0;
  C_secondhand=1000000;
  counter_bankfailure=0;
  FuelCost=0;
  Dividends_e=0;
  Wages_en=0;
  for (const string& cl:classes_mh){
    Wages_en_mh[cl]=0;
  }
  Loans_e(1)=0;
}

void DEPOSITINTEREST(void)
{
  //Firms, Households and the Energy sector receive deposit interest
  //C-firms
  for (i=1; i<=N1; i++)
	{
    sendingBank=BankingSupplier_1(i);
    InterestDeposits_1(i)=r_depo*Deposits_1(2,i);
    Deposits_1(1,i)+=InterestDeposits_1(i);
    Deposits_b(1,sendingBank)+=InterestDeposits_1(i);
    InterestDeposits(sendingBank)+=InterestDeposits_1(i);
  }

  //K-firms
  for (j=1; j<=N2; j++)
	{
    sendingBank=BankingSupplier_2(j);
    InterestDeposits_2(j)=r_depo*Deposits_2(2,j);
    InterestDeposits(sendingBank)+=InterestDeposits_2(j);
    Deposits_2(1,j)+=InterestDeposits_2(j);
    Deposits_b(1,sendingBank)+=InterestDeposits_2(j);
  }

  //Households
  InterestDeposits_h=0;
  for (const string& cl:classes_mh){
    InterestDeposits_mh[cl]=r_depo*Deposits_mh[cl](2);
    InterestDeposits_h+=InterestDeposits_mh[cl];

    Deposits_mh[cl](1)+=InterestDeposits_mh[cl];
  }
  Deposits_h(1)+=InterestDeposits_h;
  //Energy sector
  InterestDeposits_e=r_depo*Deposits_e(2);
  Deposits_e(1)+=InterestDeposits_e;

  for(i=1; i<=NB; i++)
  {
    InterestDeposits(i)+=r_depo*Deposits_hb(2,i);
    Deposits_hb(1,i)+=r_depo*Deposits_hb(2,i);
    Deposits_b(1,i)+=r_depo*Deposits_hb(2,i);
    InterestDeposits(i)+=r_depo*Deposits_eb(2,i);
    Deposits_eb(1,i)+=r_depo*Deposits_eb(2,i);
    Deposits_b(1,i)+=r_depo*Deposits_eb(2,i);
  }

  for(i=1; i<=NB; i++)
  {
    if(Deposits_hb.Row(1).Sum()>0)
    {
      DepositShare_h(i)=Deposits_hb(1,i)/Deposits_hb.Row(1).Sum();
    }
    else
    {
      DepositShare_h(i)=(NL_1(i)+NL_2(i))/(N1+N2);
    }

    if(Deposits_eb.Row(1).Sum()>0)
    {
      DepositShare_e(i)=Deposits_eb(1,i)/Deposits_eb.Row(1).Sum();
    }
    else
    {
      DepositShare_e(i)=(NL_1(i)+NL_2(i))/(N1+N2);
    }
  }
}

void MACH(void)
{
  //If simulating an energy price shock scenario, set indicator variables for passthrough
  //to zero for all firms in the period in which the shock hits
  if((flag_energyshocks==1 && t==(t_regime_shifts+1)) || (flag_energyshocks==3 && t==(t_regime_shifts+1)))
  {
    pass_1=0;
    pass_2=0;
  }

  if((flag_energyshocks==1 && t>t_regime_shifts) || (flag_energyshocks==3 && t>t_regime_shifts))
  {
    //Iterate over all K-Firms. If K-Firm i has not passed on the energy price shock yet (pass_1=0)
    //it does so in this period with probability passthrough
    for (i=1; i<=N1; i++)
    {
      if(pass_1(i)==0)
      {
        pass_1(i)=passthrough;
      }
    }
    //Same as above, but iterating over C-Firms
    for (j=1; j<=N2; j++)
	  {
      if(pass_2(j)==0)
      {
        pass_2(j)=passthrough;
      }
    }
  }

  //Iterate over all K-Firms and all machine vintages still in use (newer than t0) to determine production cost and selling price for K-firms
  for (i=1; i<=N1; i++)
  {
    for (tt=t0; tt<=t; tt++)
    {
      if (A(tt,i) > 0 & A_en(tt,i)>0)
      {
        //Calculate the unit cost of production implied by using a machine of vintage tt
        //produced by i to produce consumption goods in the current period
        C(tt,i)=w_tot_for_1_wr_mh(2)/A(tt,i)+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i);
      }
      else
      {
        //Machine labour productivity (A) and energy efficiency (A_en) should be greater than 0
        std::cerr << "\n\n ERROR: A_en(tt,i) or A(tt,i) <= 0 in period " << t << " for K-firm "<< i << endl;
        Errors << "\n A_en(tt,i) or A(tt,i) <= 0 in period " << t << " for K-firm "<< i << endl;
        exit(EXIT_FAILURE);
      }
    }

    //Calculate unit cost of production for K-Firm i (unit labour cost+unit energy cost+emission tax to be paid per unit of output)
    c1(i)=w_tot_for_1_wr_mh(2)/((1-shocks_labprod1(i))*A1p(i))+c_en(2)/((1-shocks_eneff1(i))*A1p_en(i))+t_CO2*A1p_ef(i)/((1-shocks_eneff1(i))*A1p_en(i));
    //In case of energy price shock, cost which will be used in setting the price
    //depends on whether i is already passing on the shock
    if(pass_1(i)==1)
    {
      c1p(i)=c1(i);
    }
    else
    {
      c1p(i)=w_tot_for_1_wr_mh(2)/((1-shocks_labprod1(i))*A1p(i))+(c_en_preshock+pass_1(i)*(c_en(2)-c_en_preshock))/((1-shocks_eneff1(i))*A1p_en(i))+t_CO2*A1p_ef(i)/((1-shocks_eneff1(i))*A1p_en(i));
    }

    //i can only reset its price in t with probability theta
    rnd=ran1(p_seed);
    if(rnd<=theta){
      //Price is set as simple mark-up over cost
      p1(i)=(1+mi1)*c1p(i);
    }

    //Ensure that price does not fall below a lower bound
    if (p1(i) < pmin)
    {
      p1(i)=pmin;
    }
  }


  //C-firms receive capital ordered in the last period
  for (j=1; j<=N2; j++)
	{
    //Productive capacity is increased by amount of expansion investment (EI is expressed in units producible)
    //and decreased by capacity which was scrapped in the previous period
    K(j)+=EI(2,j)-scrap_age(j);
    //Nominal value of the capital stock is updated by changes from last period
    CapitalStock(1,j)+=deltaCapitalStock(2,j);
    //Update the machine frequency arrays
    for (i=1; i<=N1; i++)
		{
		 	for (tt=t0; tt<=t; tt++)
			{
        //gtemp incorporates all changes to stocks of individual machine vintages held by
        //each firm. It is continuously updated throughout the period to record all changes
        //Here, the machine frequency array g is set to be identical to gtemp to incorporate
        //all changes which happened over the previous period.
        g[tt-1][i-1][j-1]=gtemp[tt-1][i-1][j-1];
        //g_c, g_c2 and g_c3 are copies of g which are used when calculating C-Firms' production costs
				g_c[tt-1][i-1][j-1]=g[tt-1][i-1][j-1];
        g_c2[tt-1][i-1][j-1]=g[tt-1][i-1][j-1];
        g_c3[tt-1][i-1][j-1]=g[tt-1][i-1][j-1];
      }
    }
    //Re-initialise the variable recording amount of productive capacity to be scrapped
    scrap_age(j)=0;
  }


  //C-firms determine cost of production, revise mark-up and set their price
  for (j=1; j<=N2; j++)
  {
    //Number of machines owned by j
    n_mach(j)=K(j)/dim_mach;
    for (i=1; i<=N1; i++)
    {
      for (tt=t0; tt<=t; tt++)
      {
        if (n_mach(j)>0)
        {
          //j's unit cost is a weighted average of the unit costs implied by all of the machine vintages of which j owns one or more units
          c2(j)+=(w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/((1-shocks_eneff2(j))*A_en(tt,i))+t_CO2*A_ef(tt,i)/((1-shocks_eneff2(j))*A_en(tt,i)))*g[tt-1][i-1][j-1]/n_mach(j);
          //In energy price shock scenario, the cost j uses to set its price depends on whether or not it has already passed on the price shock
          if(pass_2(j)==1)
          {
            c2p(j)+=(w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/((1-shocks_eneff2(j))*A_en(tt,i))+t_CO2*A_ef(tt,i)/((1-shocks_eneff2(j))*A_en(tt,i)))*g[tt-1][i-1][j-1]/n_mach(j);
          }
          else
          {
            c2p(j)+=(w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en_preshock/((1-shocks_eneff2(j))*A_en(tt,i))+t_CO2*A_ef(tt,i)/((1-shocks_eneff2(j))*A_en(tt,i)))*g[tt-1][i-1][j-1]/n_mach(j);
          }
          //Also record j's labour productivity, energy efficiency and environmental friendliness as weighted averages of those of the machine vintages owned by j
          A2(j)+=(1-shocks_labprod2(j))*A(tt,i)*g[tt-1][i-1][j-1]/n_mach(j);
          A2_mprod(j)+=A(tt,i)*g[tt-1][i-1][j-1]/n_mach(j);
          A2_en(j)+=(1-shocks_eneff2(j))*A_en(tt,i)*g[tt-1][i-1][j-1]/n_mach(j);
          A2_ef(j)+=A_ef(tt,i)*g[tt-1][i-1][j-1]/n_mach(j);
        }
        else
        {
          //Every C-Firm should own at least 1 machine. If j does not, this indicates a bug somewhere
          std::cerr << "\n\n ERROR n_mach = 0 in period " << t << " for C-firm "<< j << endl;
          Errors << "\n n_mach = 0 in period " << t << " for C-firm "<< j << endl;
          exit(EXIT_FAILURE);
        }
      }
    }

    //C-Firm j updates its mark-up based on the previous period's change in its market share
    if (f2(3,j)>0)
    {
      mu2(1,j)=max(0.0,mu2(2,j)*(1+deltami2*((f2(2,j)-f2(3,j))/f2(3,j))));
    }
    else
    {
      mu2(1,j)=mu2(2,j);
    }

    if (mu2(1,j) <= 0)
    {
      //Mark-ups should be positive
      std::cerr << "\n\n ERROR: mark-up out of range in period " << t << " for C-firm "<< j << endl;
      Errors << "\n Mark-up out of range in period " << t << " for C-firm "<< j << endl;
      exit(EXIT_FAILURE);
    }

    //As in the case of K-Firms, C-Firm j can only reset its price in t with probability theta
    rnd=ran1(p_seed);
    if(rnd<=theta){
      p2(j)=(1+mu2(1,j))*c2p(j);
    }

    //Ensure that price does not fall below a lower bound
    if (p2(j) < pmin)
    {
      p2(j)=pmin;
    }



  }
}

void BROCHURE(void)
{
  //K-Firms send brochures to attract customers
  //Re-initialise number of clients of each K-Firm
  nclient=0;

  //If a C-Firm does not have a valid supplier of capital goods, randomly assign one
  //This happens e.g. if j's supplier exited in the previous period
	for (j=1; j<=N2; j++)
	{
    if (fornit(j) < 1 || fornit(j) > N1)
		{
      fornit(j)=int(ran1(p_seed)*N1*N2)%N1+1;
      Match(j,fornit(j))=1;

    }
  }

  //K-firms send brochures to potential customers
	for (i=1; i<=N1; i++)
	{
		//Count number of C-Firms matched to K-Firm i
    for (j=1; j<=N2; j++)
    {
			nclient(i)+=Match(j,i);
    }

    //Number of brochures sent is a function of number of existing clients
		newbroch=int(ROUND(nclient(i)*Gamma));

    //Ensure that every firm sends at least 1 brochure
    if (newbroch==0)
    {
      newbroch++;
    }

    //Brochures are sent to randomly drawn C-Firms
    //Note that this does not ensure that a randomly drawn C-Firm is not already a customer of i!
    while (newbroch > 0)
		{
			rni=int(ran1(p_seed)*N1*N2)%N2+1;
      Match(rni,i)=1;
      newbroch--;
		}
	}

  //Calculate attractiveness of each K-firm
  for (i=1; i<=N1; i++){
    Attractivess_1_i(i) = p1(i) - EnvSubsidiesPerMachine_1_i(i) + (w_tot_for_1_wr_mh(2)/A1(i)+c_en(2)/A1_en(i)+t_CO2*A1_ef(i)/A1_en(i))*b;
  }

  //C-firms choose their preferred supplier of machine tools
  for (j=1; j<=N2; j++)
	{
		//current supplier of j as integer
    indforn=int(fornit(j));
		for (i=1; i<=N1; i++)
		{
			if (A1(i) > 0)
			{
        //If j has received a brochure from i and i's technology is more convenient, j switches
        if (Match(j,i)==1 && Attractivess_1_i(i) < Attractivess_1_i(indforn))
        {
          indforn=i;
        }
			}
			else
      {
        //Technology offered by i should imply a positive labour productivity
        std::cerr << "\n\n ERROR: A1(i) = 0 in period " << t << " for K-firm "<< i << endl;
        Errors << "\n A1(i) = 0 in period " << t << " for K-firm "<< i << endl;
        exit(EXIT_FAILURE);
      }
    }
    //Update supplier of j
		fornit(j)=indforn;

    //Reset the K-Firm C-Firm network by setting the entries of all i from which j has received
    //brochures but which j did not choose as its supplier to zero
    for (i=1; i<=N1; i++)
		{
			if (i != indforn)
      {
				Match(j,i)=0;
      }
		}
	}

  //Update the number of clients of each K-Firm
	nclient=0;
	for (i=1; i<=N1; i++)
	{
		for (j=1; j<=N2; j++)
		{
			if (Match(j,i) == 1)
      {
				nclient(i)++;
      }
		}
	}
}

void INVEST(void)
{
	for (j=1; j<=N2; j++)
	{
    //C-firms determine expected demand, desired production, and demand for investment
    De(j)=alfa*De(j)+(1-alfa)*D2(2,j);
    if (De(j)<=0)
    {
        De(j)=1;
    }

    //Desired inventories
    Ne(j)=De(j)*omicron;
    //Desired output
		Qd(j)=De(j)+Ne(j)-N(2,j);

		if (Qd(j) < 0)
    {
			Qd(j)=0;
    }

    //Desired output implies desired capital stock
    Kd(j)=Qd(j)/u;

    //C-firm determines which machines should be scrapped
    SCRAPPING();

    //Determine the capital stock which the firm will have available in t+1 once aged machines are scrapped
    //Machines scrapped due to age can still be used in t but will be removed at the beginning of t+1
		Ktrig(j)=ceil((K(j)-scrap_age(j))/dim_mach)*dim_mach;

    //If the desired capital stock is larger than what the firm will have available after scrapping, it wants to engage in expansion investment
		if (Kd(j) >= Ktrig(j))
		{
      //If there is a constraint on expansion investment, determine the maximum capital stock which can be achieved
      if(I_max>0)
      {
        K_top=K(j)*(1+I_max);
        K_top=ROUND(K_top/dim_mach)*dim_mach;
        if(K_top<(K(j)+dim_mach))
        {
          K_top+=dim_mach;
        }
      }
      else
      {
        K_top=Kd(j)+1;
      }

      if(Kd(j)>K_top)
      {
        EId(j)=K_top-Ktrig(j);
      }
      else
      {
        EId(j)=ceil((Kd(j)-(K(j)-scrap_age(j)))/dim_mach)*dim_mach;
      }
    }
		else
    {
      EId(j)=0;
    }

    if(SId(j)==0 && EId(j)==0 && marker_age(j)==1)
    {
      EId(j)=dim_mach;
    }

		if (Qd(j) > K(j))
		{
			Qd(j)=K(j);
		}

    //C-firms determine effective production cost
		if (Qd(j) > 0 && Qd(j) < K(j))
		{
			COSTPROD();
		}
		else
		{
			A2e(j)=A2(j);
			c2e(j)=c2(j);
      A2e_en(j)=A2_en(j);
      A2e_ef(j)=A2_ef(j);
		}
	}

	ORD();
}

void SCRAPPING(void)
{
  //Put capital stock of j in temporary storage
  K_temp(j)=K(j)/dim_mach;
  indforn=int(fornit(j));
  //C-firm j determines which machines should be scrapped due to age and/or due to superior tech being available
  for (i=1; i<=N1; i++)
	{
		for (tt=t0; tt<=t; tt++)
		{
			//Initialise array elements to hold unit cost and quantity of machines to be scrapped
      C_pb[tt-1][i-1][j-1]=0;
			g_pb[tt-1][i-1][j-1]=0;

      //If a machine has reached its maximum age, it is scrapped, unless the firm has only 1 machine remaining
      if (g[tt-1][i-1][j-1] > 0 && age[tt-1][i-1][j-1] > (agemax))
			{
        g_pb[tt-1][i-1][j-1]=min(g[tt-1][i-1][j-1],(K_temp(j)-1));
        C_pb[tt-1][i-1][j-1]=C(tt,i);
        scrap_age(j)+=dim_mach*g_pb[tt-1][i-1][j-1];
        K_temp(j)-=g_pb[tt-1][i-1][j-1];
        //If the firm has only one machine left, set the marker_age flag. This will ensure later that
        //its desired investment is at least 1 machine
        if(K_temp(j)==1)
        {
          marker_age(j)=1;
        }
			}
      //REPLACEMENT DECISION
        //Machines which have not reached their maximum age are scrapped if a superior technology is available from the current supplier
        //Only enter here if the relevant machines have not already been scrapped due to age
      if (g[tt-1][i-1][j-1] > 0 && g_pb[tt-1][i-1][j-1]==0)
			{
        if (w_tot_for_1_wr_mh(2) > 0 && A(tt,i) > 0 && A1(indforn) > 0 && A1_en(indforn)>0 && A_en(tt,i)>0)
        {
          //calculate payback variable: (price of machine)/(unit cost difference offered by the machine)
          //Payback can be interpreted as number of units which have to be produced with a new machine such that
          //the savings in unit cost are equal to the purchase price of the machine
          payback=(p1(indforn)-EnvSubsidiesPerMachine_1_i(indforn))/(w_tot_for_1_wr_mh(2)/A(tt,i)+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i)-w_tot_for_1_wr_mh(2)/A1(indforn)-c_en(2)/A1_en(indforn)-t_CO2*A1_ef(indforn)/A1_en(indforn));
        }
        else
        {
          //All of the variables in the denominators of the formula above should never become zero or negative
          std::cerr << "\n\n ERROR: payback division by zero in period " << t << " for C-firm "<< j << endl;
          Errors << "\n Payback division by zero in period " << t << " for C-firm "<< j << endl;
          exit(EXIT_FAILURE);
        }

        //If payback is smaller than an exogenous threshold, firm wants to replace the machine in question with the newer vintage
				if (payback <= b && payback>0)
				{
					g_pb[tt-1][i-1][j-1]=g[tt-1][i-1][j-1];
		 			C_pb[tt-1][i-1][j-1]=C(tt,i);
					SId(j)+=dim_mach*g_pb[tt-1][i-1][j-1];
				}
			}
		}
	}
}

void COSTPROD(void)
{
  //C-firms determine effective production cost based on desired production; most efficient machines used first
  nmachprod=ceil(Qd(j)/dim_mach);
	nmp_temp=nmachprod;

	while (nmp_temp > 0)
	{
    cmin=10000000000000000;
    imin=0;
    jmin=0;
    tmin=0;

    for (i=1; i<=N1; i++)
    {
      for (tt=t0; tt<=t; tt++)
      {
        if (g_c[tt-1][i-1][j-1] > 0 && (w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i)) < cmin)
        {
          cmin=w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i);
          imin=i;
          jmin=j;
          tmin=tt;
        }
      }
    }

    if (nmachprod>0)
    {
      if (g_c[tmin-1][imin-1][jmin-1] >= nmp_temp)
      {
        A2e(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*nmp_temp/nmachprod;
        A2e_en(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*nmp_temp/nmachprod;
        A2e_ef(j)+=A_ef(tmin,imin)*nmp_temp/nmachprod;
        c2e(j)+=(w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tmin,imin))+c_en(2)/((1-shocks_eneff2(j))*A_en(tmin,imin))+t_CO2*A_ef(tmin,imin)/((1-shocks_eneff2(j))*A_en(tmin,imin)))*nmp_temp/nmachprod;
        g_c[tmin-1][imin-1][jmin-1]-= nmp_temp;
        nmp_temp=0;
      }
      else
      {
        A2e(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*g_c[tmin-1][imin-1][jmin-1]/nmachprod;
        A2e_en(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*g_c[tmin-1][imin-1][jmin-1]/nmachprod;
        A2e_ef(j)+=A_ef(tmin,imin)*g_c[tmin-1][imin-1][jmin-1]/nmachprod;
        c2e(j)+=(w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tmin,imin))+c_en(2)/((1-shocks_eneff2(j))*A_en(tmin,imin))+t_CO2*A_ef(tmin,imin)/((1-shocks_eneff2(j))*A_en(tmin,imin)))*g_c[tmin-1][imin-1][jmin-1]/nmachprod;
        nmp_temp-=g_c[tmin-1][imin-1][jmin-1];
        g_c[tmin-1][imin-1][jmin-1]=0;
      }
    }
    else
    {
      std::cerr << "\n\n ERROR: nmachprod = 0!!!" << endl;
      Errors << "\n nmachprod = 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
	}
}

void ORD(void)
{
  if((flag_energyshocks>0 && t>t_regime_shifts))
  //Transfer to C-firms
  {
    for (j=1; j<=N2; ++j)
    {
      Transfer_shock_f2(j)=tref*max(0.0,(c_en(2)-c_en_preshock)*(Qd(j)/A2e(j)));
      receivingBank=BankingSupplier_2(j);
      Deposits_2(1,j)+=Transfer_shock_f2(j);
      Deposits_b(1,receivingBank)+=Transfer_shock_f2(j);
      Inflows(receivingBank)+=Transfer_shock_f2(j);
      Transfer_shock_f+=Transfer_shock_f2(j);
      mol(j)+=Transfer_shock_f2(j);
    }
  }

  //C-firms calculate internal funds & the maximum amount of loans they are willing to take on; based on this, investment is possibly scaled back
  for (j=1; j<=N2; ++j)
  {
		MaxFunds=max(0.0,Deposits_2(1,j)+phi2*mol(j)-Loans_2(1,j)-c2e(j)*Qd(j));

    indforn=int(fornit(j));
		p1prova=p1(indforn);
    if ((EId(j)/dim_mach)*p1prova < MaxFunds)
    {
      EIp(j)=EId(j);
      MaxFunds-=(EId(j)/dim_mach)*p1(indforn);
    }
    else
    {
      p1prova=p1(indforn);
      EIp(j)=floor(MaxFunds/p1(indforn))*dim_mach;
      if(EIp(j)<0)
      {
        EIp(j)=0;
        MaxFunds=0;
      }else
      {
        MaxFunds=0;
      }
    }

    if ((SId(j)/dim_mach)*p1prova < MaxFunds)
    {
      SIp(j)=SId(j);
    }
    else
    {
      SIp(j)=floor(MaxFunds/p1(indforn))*dim_mach;
    }

    Ip(j)=EIp(j)+SIp(j);

    //Determine cost of planned investment; 
		if (Ip(j) > 0)
		{
      CmachEI(j)=p1(indforn)*EIp(j)/dim_mach;
			CmachSI(j)=p1(indforn)*SIp(j)/dim_mach;
			Cmach(j)=p1(indforn)*Ip(j)/dim_mach;
		}
		else
    {
      CmachEI(j)=0;
			CmachSI(j)=0;
			Cmach(j)=0;
    }
  }
}

void ALLOCATECREDIT(void)
{
  //C-firms calculate credit demand based on outstanding loans (need to be rolled over), planned investment and production
  for (j=1; j<=N2; ++j)
  {
		if (Loans_2(1,j)+Cmach(j)+(c2e(j)*Qd(j))<=Deposits_2(1,j))
    {
      CreditDemand(j)=0;
    }
		else
    {
      CreditDemand(j)=Loans_2(1,j)+Cmach(j)+(c2e(j)*Qd(j))-Deposits_2(1,j);
    }
	}

  //Banks allocate credit
  for (i=1; i<=NB; i++)
  {
    for (j=1; j<=NL_2(i); j++)
    {
      DS2_rating.Column(i).Minimum1(rated_firm_2);
      DS2_rating(rated_firm_2,i)=DS2_rating.Column(i).Maximum()+1;
			if (BankMatch_2(rated_firm_2,i)==1)
      {
        //If customer does not need credit, they use deposits to repay their outstanding loans
        if (CreditDemand(rated_firm_2)==0)
        {
			    Q2(rated_firm_2)=Qd(rated_firm_2);
          EI(1,rated_firm_2)=EIp(rated_firm_2);
          SI(rated_firm_2)=SIp(rated_firm_2);
          I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
          Deposits_2(1,rated_firm_2)=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
          Loans_b(1,i)=max(0.0,Loans_b(1,i)-Loans_2(1,rated_firm_2));
          Deposits_b(1,i)=Deposits_b(1,i)-Loans_2(1,rated_firm_2);
          Loans_2(1,rated_firm_2)=0;
        }
        else if(CreditDemand(rated_firm_2)>0)
        {
			    //If remaining credit supply of bank is sufficient, demand is fully satisfied
          if (CreditDemand(rated_firm_2) <= BankCredit(i))
          {
            Deposits_2(1,rated_firm_2)=max(0.0,Deposits_2(1,rated_firm_2)+CreditDemand(rated_firm_2)-Loans_2(1,rated_firm_2));
            Deposits_b(1,i)=Deposits_b(1,i)+CreditDemand(rated_firm_2)-Loans_2(1,rated_firm_2);
            Loans_b(1,i)=Loans_b(1,i)+CreditDemand(rated_firm_2)-Loans_2(1,rated_firm_2);
            Loans_2(1,rated_firm_2)=CreditDemand(rated_firm_2);
				    BankCredit(i)-=CreditDemand(rated_firm_2);
            Q2(rated_firm_2)=Qd(rated_firm_2);
				    EI(1,rated_firm_2)=EIp(rated_firm_2);
    				SI(rated_firm_2)=SIp(rated_firm_2);
    				I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
          }
          else
          {
				    //If remaining credit supply is insufficient, first remove replacement investment
            if(Loans_2(1,rated_firm_2)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2) <= BankCredit(i))
            {
              if (Loans_2(1,rated_firm_2)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2)>=0)
              {
                Loans_b(1,i)=Loans_b(1,i)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
                Loans_2(1,rated_firm_2)=Loans_2(1,rated_firm_2)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
				        Deposits_b(1,i)=Deposits_b(1,i)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
                Deposits_2(1,rated_firm_2)=CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2));
              }
              else
					    {
                Loans_b(1,i)=max(0.0,Loans_b(1,i)-Loans_2(1,rated_firm_2));
                Deposits_b(1,i)=Deposits_b(1,i)-Loans_2(1,rated_firm_2);
                Deposits_2(1,rated_firm_2)=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                Loans_2(1,rated_firm_2)=0;
              }

					    BankCredit(i)-=Loans_2(1,rated_firm_2);
              Q2(rated_firm_2)=Qd(rated_firm_2);
              EI(1,rated_firm_2)=EIp(rated_firm_2);
              SI(rated_firm_2)=0;
              I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
            }
            //If remaining credit supply is still insufficient, remove also expansion investment
            else if (Loans_2(1,rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2)<=BankCredit(i))
            {
              if (Loans_2(1,rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2)>=0)
              {
                Loans_b(1,i)=Loans_b(1,i)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
                Loans_2(1,rated_firm_2)=Loans_2(1,rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
				        Deposits_b(1,i)=Deposits_b(1,i)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
                Deposits_2(1,rated_firm_2)=(c2e(rated_firm_2)*Qd(rated_firm_2));
              }
              else
              {
                Loans_b(1,i)=max(0.0,Loans_b(1,i)-Loans_2(1,rated_firm_2));
                Deposits_b(1,i)=Deposits_b(1,i)-Loans_2(1,rated_firm_2);
                Deposits_2(1,rated_firm_2)=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                Loans_2(1,rated_firm_2)=0;
              }
              BankCredit(i)-=Loans_2(1,rated_firm_2);
              Q2(rated_firm_2)=Qd(rated_firm_2);
              EI(1,rated_firm_2)=0;
              SI(rated_firm_2)=0;
              I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
            }
            else if (Loans_2(1,rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2)>BankCredit(i))
				    {
              //If possible, scale back desired production to a level which can be financed
              if (Loans_2(1,rated_firm_2)-Deposits_2(1,rated_firm_2)<=BankCredit(i))
              {
						    Q2(rated_firm_2)=(BankCredit(i)-Loans_2(1,rated_firm_2)+Deposits_2(1,rated_firm_2))/c2e(rated_firm_2);

                //If production needs to be scaled back too much, firm exits
                if (Q2(rated_firm_2) < 1)
						    {
                  if (Loans_2(1,rated_firm_2)>Deposits_2(1,rated_firm_2))
                  {
                    Loans_b(1,i)-=Deposits_2(1,rated_firm_2);
                    Deposits_b(1,i)-=Deposits_2(1,rated_firm_2);
                    Loans_2(1,rated_firm_2)-=Deposits_2(1,rated_firm_2);
                    Deposits_2(1,rated_firm_2)=0;
								    baddebt_2(rated_firm_2)=Loans_2(1,rated_firm_2);
                    baddebt_b(i)+=Loans_2(1,rated_firm_2);
    								BankCredit(i)-=Loans_2(1,rated_firm_2);
                    Loans_b(1,i)=max(0.0, Loans_b(1,i)-Loans_2(1,rated_firm_2));
                    Loans_2(1,rated_firm_2)=0;
                  }
                  else
                  {
                    Loans_b(1,i)=max(0.0, Loans_b(1,i)-Loans_2(1,rated_firm_2));
                    baddebt_2(rated_firm_2)=Loans_2(1,rated_firm_2)-Deposits_2(1,rated_firm_2);
                    Deposits_recovered_2+=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                    for (const string& cl:classes_mh){
                      Deposits_recovered_2_mh[cl] += (Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2)) * Ownership_sh_2_i_mh[cl](j);
                    }
                    Deposits_b(1,i)-=Deposits_2(1,rated_firm_2);
                    Outflows(i)+=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                    Deposits_2(1,rated_firm_2)=0;
                    Loans_2(1,rated_firm_2)=0;
                  }
                  Q2(rated_firm_2)=0;
                  exiting_2(rated_firm_2)=1;
    							f2(1,rated_firm_2)=0;
    							f2(2,rated_firm_2)=0;
    							f2(3,rated_firm_2)=0;
                }
                else
                {
                  Loans_b(1,i)=Loans_b(1,i)+BankCredit(i)-Loans_2(1,rated_firm_2);
                  Deposits_2(1,rated_firm_2)=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2)+BankCredit(i);
                  Deposits_b(1,i)=Deposits_b(1,i)-Loans_2(1,rated_firm_2)+BankCredit(i);
                  Loans_2(1,rated_firm_2)=BankCredit(i);
                  BankCredit(i)=0;
                }
    						EI(1,rated_firm_2)=0;
    						SI(rated_firm_2)=0;
    						I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
              }
              //Otherwise, firm is forced to exit
              else
              {
  						  Q2(rated_firm_2)=0;
                BankCredit(i)=0;
                if (Loans_2(1,rated_firm_2)>Deposits_2(1,rated_firm_2))
                {
                  Loans_b(1,i)-=Deposits_2(1,rated_firm_2);
                  Deposits_b(1,i)-=Deposits_2(1,rated_firm_2);
                  Loans_2(1,rated_firm_2)-=Deposits_2(1,rated_firm_2);
                  Deposits_2(1,rated_firm_2)=0;
                  baddebt_2(rated_firm_2)=Loans_2(1,rated_firm_2);
                  baddebt_b(i)+=Loans_2(1,rated_firm_2);
                  Loans_b(1,i)=max(0.0, Loans_b(1,i)-Loans_2(1,rated_firm_2));
                  Loans_2(1,rated_firm_2)=0;
                }
                else
                {
                  Loans_b(1,i)=max(0.0, Loans_b(1,i)-Loans_2(1,rated_firm_2));
                  baddebt_2(rated_firm_2)=Loans_2(1,rated_firm_2)-Deposits_2(1,rated_firm_2);
                  Deposits_recovered_2+=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                  for (const string& cl:classes_mh){
                    Deposits_recovered_2_mh[cl] += (Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2)) * Ownership_sh_2_i_mh[cl](j);
                  }
                  Deposits_b(1,i)-=Deposits_2(1,rated_firm_2);
                  Outflows(i)+=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                  Deposits_2(1,rated_firm_2)=0;
                  Loans_2(1,rated_firm_2)=0;
                }
    						EI(1,rated_firm_2)=0;
    						SI(rated_firm_2)=0;
    						I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
                exiting_2(rated_firm_2)=1;
    						f2(1,rated_firm_2)=0;
    						f2(2,rated_firm_2)=0;
    						f2(3,rated_firm_2)=0;
              }
            }
          }
        }
      }
    }
  }


  //C-firms calculate labour demand based on production which can be financed
	for (j=1; j<=N2; ++j)
  {
      if (A2e(j) > 0)
      {
        Ld2_wr(j)=Q2(j)/A2e(j); 
        for (const string& cl:classes_mh){
          double cl_ld=Ld2_wr(j)*ld_ratios_mh[cl];
          Ld2_i_mh[cl](j)=cl_ld;
        }
      }
      else
      {
        std::cerr << "\n\n ERROR: A2e(j) = 0 in period " << t << " for C-firm " << j << endl;
        Errors << "\n ERROR: A2e(j) = 0 in period " << t << " for C-firm " << j << endl;
        exit(EXIT_FAILURE);
      }
  }
  LD2_mh=umap_sum_each_key(Ld2_i_mh); //Total L demand for each household class

}

void PRODMACH(void)
{
  //K-firms receive demand from C-firms and calculate labour demand
  for (i=1; i<=N1; i++)
  {
		for (j=1; j<=N2; j++)
		{
      if (Match(j,i) == 1)
      {
          D1(i)+=I(j)/dim_mach;
      }
	  }

    Q1(i)=D1(i);

    if (A1p(i)>0)
    {
			Ld1_prod_wr(i)=Q1(i)/((1-shocks_labprod1(i))*A1p(i));
      for (const string& cl:classes_mh){
        double cl_ld=Ld1_prod_wr(i)*ld_ratios_mh[cl];
        Ld1_prod_i_mh[cl](i)=cl_ld;
        Ld1_i_mh[cl](i)=Ld1_prod_i_mh[cl](i)+Ld1_rd_i_mh[cl](i); //Sum R&D demand from previous period
      }
    }
    else
    {
      std::cerr << "\n\n ERROR: A1p(i)=0 in period " << t << " for K-firm " << i << endl;
      Errors << "\n A1p(i)=0 in period " << t << " for K-firm " << i << endl;
      exit(EXIT_FAILURE);
    }
  }

  if((flag_energyshocks>0 && t>t_regime_shifts))
  //Transfer to K-firms
  {
    for (i=1; i<=N1; i++)
    {
      Transfer_shock_f1(i)=tref*max(0.0,(c_en(2)-c_en_preshock)*(Q1(i)/A1p_en(i)));
      receivingBank=BankingSupplier_1(i);
      Deposits_1(1,i)+=Transfer_shock_f1(i);
      Deposits_b(1,receivingBank)+=Transfer_shock_f1(i);
      Inflows(receivingBank)+=Transfer_shock_f1(i);
      Transfer_shock_f+=Transfer_shock_f1(i);
    }
  }
  LD1_mh=umap_sum_each_key(Ld1_i_mh); //Total L demand for each household class

  //Determine total labour demand
  LABOR();

  //Calculate Energy demand based on actual production
  for (i=1; i<=N1; i++)
  {
    if (A1p_en(i)>0)
    {
      D1_en(i)=Q1(i)/((1-shocks_eneff1(i))*A1p_en(i));
    }
    else
    {
      std::cerr << "\n\n ERROR: A1p_en(i)=0 in period " << t << " for K-firm " << i << endl;
      Errors << "\n A1p_en(i) = 0 in period " << t << " for K-firm " << i << endl;
      exit(EXIT_FAILURE);
    }

    S1_pre(i)=Q1(i)*p1(i);
  }

  for (j=1; j<=N2; j++)
  {
    if (A2e_en(j) > 0)
    {
      D2_en(j)=Q2(j)/A2e_en(j);
    }
    else
    {
      std::cerr << "\n\n ERROR: A2e_en(j) = 0 in period " << t << " for C-firm " << j <<endl;
      Errors << "\n A2e_en(j) = 0 in period " << t << " for C-firm " << j << endl;
      exit(EXIT_FAILURE);
    }
  }

  if((flag_outputshocks==1))
  {
    for (i=1; i<=N1; i++)
    {
      Qpast=Q1(i);
      if (Qpast > 0)
      {
        Q1(i)=floor(Q1(i)*(1-shocks_output1(i)));
        loss=Qpast-Q1(i);
        I_loss(i)=loss*dim_mach;
        while(loss>0)
        {
          ranj=int(ran1(p_seed)*N1*N2)%N2+1;
          if (Match(ranj,i) == 1 && I(ranj)>0)
          {
            I(ranj)-=dim_mach;
            if (I(ranj)<EI(1,ranj))
            {
              EI(1,ranj)=I(ranj);
            }
            SI(ranj)=I(ranj)-EI(1,ranj);
            loss-=1;
          }
        }
      }
    }
  }

  if((flag_outputshocks==2))
  {
    loss=ceil(Q1.Sum()*shocks_output1(1));
    while(loss>0)
    {
      rani=int(ran1(p_seed)*N1*N2)%N1+1;
      Qpast=Q1(rani);
      if (Qpast > 0)
      {
        if(Q1(rani)>=loss)
        {
          Q1(rani)-=loss;
          loss=0;
        }
        else
        {
          loss-=Q1(rani);
          Q1(rani)=0;
        }
        lossj=Qpast-Q1(rani);
        I_loss(rani)=lossj*dim_mach;
        while(lossj>0)
        {
          ranj=int(ran1(p_seed)*N1*N2)%N2+1;
          if (Match(ranj,rani) == 1 && I(ranj)>0)
          {
            I(ranj)-=dim_mach;
            if (I(ranj)<EI(1,ranj))
            {
              EI(1,ranj)=I(ranj);
            }
            SI(ranj)=I(ranj)-EI(1,ranj);
            lossj-=1;
          }
        }
      }
    }
  }

  if((flag_outputshocks==3))
  {
    loss=ceil(Q1.Sum()*shocks_output1(1));
    while(loss>0)
    {
      rani=0;
      while(rani==0)
      {
        rani=int(ran1(p_seed)*N1*N2)%N1+1;
      }
      Qpast=Q1(rani);
      if (Qpast > 0)
      {
        rnd=ran1(p_seed);
        lossj=ceil(rnd*Qpast);
        if(lossj>=loss)
        {
          Q1(rani)-=loss;
          loss=0;
        }
        else
        {
          loss-=lossj;
          Q1(rani)-=lossj;
        }
        lossj=Qpast-Q1(rani);
        I_loss(rani)=lossj*dim_mach;
        while(lossj>0)
        {
          ranj=int(ran1(p_seed)*N1*N2)%N2+1;
          if (Match(ranj,rani) == 1 && I(ranj)>0)
          {
            I(ranj)-=dim_mach;
            if (I(ranj)<EI(1,ranj))
            {
              EI(1,ranj)=I(ranj);
            }
            SI(ranj)=I(ranj)-EI(1,ranj);
            lossj-=1;
          }
        }
      }
    }
  }

  //Calculate nominal value of actual investment
  for (j=1; j<=N2; j++)
  {
    indforn=int(fornit(j));
    SI_n(j)=p1(indforn)*SI(j)/dim_mach;
    EI_n(j)=p1(indforn)*EI(1,j)/dim_mach;
    S1(indforn)+=p1(indforn)*I(j)/dim_mach;
    S1_post(indforn)+=p1(indforn)*I(j)/dim_mach;
  }

  //Old machines are scrapped; temporary machine frequency arrays are updated based on expansion & replacement investment
	for (j=1; j<=N2; j++)
	{

    CANCMACH();

		for (i=1; i<=N1; i++)
		{
			for (tt=t0; tt<=t; tt++)
			{
				if (gtemp[tt-1][i-1][j-1] == 0)
        {
					age[tt-1][i-1][j-1]=0;
        }
			}
		}

		indforn=int(fornit(j));
		gtemp[t-1][indforn-1][j-1]+=I(j)/dim_mach;
    g_price[t-1][indforn-1][j-1]=p1(indforn);
    deltaCapitalStock(1,j)+=EI_n(j);

		if (I(j) > 0)
		{
			if (gtemp[t-1][indforn-1][j-1] == I(j)/dim_mach){
        age[t-1][indforn-1][j-1]=0;
      }
		}
	}

  if (flag_capshocks==1)
  {
    //Shock to C-firms' capital stock
    for (j=1; j<=N2; j++)
    {
        K_cur(j)=K(j);
        loss=0;
        for (i=1; i<=N1; i++)
        {
            for (tt=t0; tt<=(t); tt++)
            {
                loss+=gtemp[tt-1][i-1][j-1];
            }
        }
        loss=ROUND(loss*shocks_capstock(j));
        K_loss(j)=loss*dim_mach;
        while(loss>0)
        {
            rani=int(ran1(p_seed)*N1*N2)%N1+1;
            rant=int(ran1(p_seed)*t0*(t))%((t)-t0+1)+t0;
            if (gtemp[rant-1][rani-1][j-1]>0)
            {
                if(loss<gtemp[rant-1][rani-1][j-1])
                {
                    if(rant==t)
                    {
                      deltaCapitalStock(1,j)-=loss*g_price[rant-1][rani-1][j-1];
                      Loss_Capital_mol(j)+=loss*g_price[rant-1][rani-1][j-1];
                      gtemp[rant-1][rani-1][j-1]-=loss;
                      if(SI(j)>0)
                      {
                        K(j)-=min(SI(j),loss*dim_mach);
                        loss=max(0.0,loss*dim_mach-SI(j))/dim_mach;
                      }
                      EI(1,j)-=loss*dim_mach;
                      EI(1,j)=max(EI(1,j),0.0);
                      loss=0;
                    }
                    else
                      {
                        Loss_Capital(j)+=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*loss);
                        Loss_Capital_mol(j)+=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*loss);
                        CapitalStock(1,j)-=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*loss);
                        K(j)-=loss*dim_mach;
                        K_cur(j)-=loss*dim_mach;
                        gtemp[rant-1][rani-1][j-1]-=loss;
                        g_c2[rant-1][rani-1][j-1]-=loss;
                        g_c3[rant-1][rani-1][j-1]-=loss;
                        loss=0;
                      }
                }
                else
                {
                    if(rant==t)
                    {
                      deltaCapitalStock(1,j)-=g_price[rant-1][rani-1][j-1]*gtemp[rant-1][rani-1][j-1];
                      Loss_Capital_mol(j)+=g_price[rant-1][rani-1][j-1]*gtemp[rant-1][rani-1][j-1];
                      if(SI(j)>0)
                      {
                        K(j)-=gtemp[rant-1][rani-1][j-1]*dim_mach-EI(1,j);
                      }
                      EI(1,j)=0;
                      loss-=gtemp[rant-1][rani-1][j-1];
                      gtemp[rant-1][rani-1][j-1]=0;
                    }
                    else
                    {
                      Loss_Capital(j)+=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*gtemp[rant-1][rani-1][j-1]);
                      Loss_Capital_mol(j)+=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*gtemp[rant-1][rani-1][j-1]);
                      CapitalStock(1,j)-=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*gtemp[rant-1][rani-1][j-1]);
                      K(j)-=gtemp[rant-1][rani-1][j-1]*dim_mach;
                      K_cur(j)-=gtemp[rant-1][rani-1][j-1]*dim_mach;
                      loss-=gtemp[rant-1][rani-1][j-1];
                      g_c2[rant-1][rani-1][j-1]-=gtemp[rant-1][rani-1][j-1];
                      g_c3[rant-1][rani-1][j-1]-=gtemp[rant-1][rani-1][j-1];
                      gtemp[rant-1][rani-1][j-1]=0;
                    }
                }
            }
        }
      if(shocks_capstock(j)>0)
      {
        Q2(j)=min(Q2(j),K_cur(j));
        ADJUSTEMISSENLAB();
      }
    }
  }

  if (flag_capshocks==2)
  {
    //Shock to C-firms' capital stock
    loss=0;
    K_temp_sum=0;
    K_cur=K;
    for (j=1; j<=N2; j++)
    {
      for (i=1; i<=N1; i++)
      {
          for (tt=t0; tt<=(t); tt++)
          {
              loss+=gtemp[tt-1][i-1][j-1];
              K_temp_sum+=gtemp[tt-1][i-1][j-1];
          }
      }
    }

    if(loss==N2)
    {
      loss=0;
    }
    else
    {
      loss=min(K_temp_sum-N2,ROUND(loss*shocks_capstock(1)));
    }

    while(loss>0)
    {
      ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      lossj=0;
      for (i=1; i<=N1; i++)
      {
          for (tt=t0; tt<=(t); tt++)
          {
              lossj+=gtemp[tt-1][i-1][ranj-1];
          }
      }

      lossj=min(loss,(lossj-1));
      K_loss(ranj)=lossj*dim_mach;
      while(lossj>0)
      {
          rani=int(ran1(p_seed)*N1*N2)%N1+1;
          rant=int(ran1(p_seed)*t0*(t))%((t)-t0+1)+t0;
          if (gtemp[rant-1][rani-1][ranj-1]>0)
          {
              if(lossj<gtemp[rant-1][rani-1][ranj-1])
              {
                  if(rant==t)
                  {
                      deltaCapitalStock(1,ranj)-=lossj*g_price[rant-1][rani-1][ranj-1];
                      Loss_Capital_mol(ranj)+=lossj*g_price[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]-=lossj;
                      loss-=lossj;
                      if(SI(ranj)>0)
                      {
                        K(ranj)-=min(SI(ranj),lossj*dim_mach);
                        lossj=max(0.0,lossj*dim_mach-SI(ranj))/dim_mach;
                      }
                      EI(1,ranj)-=lossj*dim_mach;
                      EI(1,ranj)=max(EI(1,ranj),0.0);
                      lossj=0;
                  }
                  else
                  {
                      Loss_Capital(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*lossj);
                      Loss_Capital_mol(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*lossj);
                      CapitalStock(1,ranj)-=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*lossj);
                      K(ranj)-=lossj*dim_mach;
                      K_cur(ranj)-=lossj*dim_mach;
                      gtemp[rant-1][rani-1][ranj-1]-=lossj;
                      g_c2[rant-1][rani-1][ranj-1]-=lossj;
                      g_c3[rant-1][rani-1][ranj-1]-=lossj;
                      loss-=lossj;
                      lossj=0;
                  }
              }
              else
              {
                  if(rant==t)
                  {
                      deltaCapitalStock(1,ranj)-=g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1];
                      Loss_Capital_mol(ranj)+=g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1];
                      if(SI(ranj)>0)
                      {
                        K(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach-EI(1,ranj);
                      }
                      EI(1,ranj)=0;
                      lossj-=gtemp[rant-1][rani-1][ranj-1];
                      loss-=gtemp[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]=0;
                  }
                  else
                  {
                      Loss_Capital(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1]);
                      Loss_Capital_mol(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1]);
                      CapitalStock(1,ranj)-=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1]);
                      K(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach;
                      K_cur(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach;
                      lossj-=gtemp[rant-1][rani-1][ranj-1];
                      loss-=gtemp[rant-1][rani-1][ranj-1];
                      g_c2[rant-1][rani-1][ranj-1]-=gtemp[rant-1][rani-1][ranj-1];
                      g_c3[rant-1][rani-1][ranj-1]-=gtemp[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]=0;
                  }
              }
          }
      }
    }

    if(K_loss.Sum()>0)
    {
      for (j=1; j<=N2; j++)
      {
        if(K_loss(j)>0)
        {
          Q2(j)=min(Q2(j),K_cur(j));
          ADJUSTEMISSENLAB();
        }
      }
    }
  }

  if (flag_capshocks==3)
  {
    //Shock to C-firms' capital stock
    loss=0;
    K_temp_sum=0;
    K_cur=K;
    for (j=1; j<=N2; j++)
    {
      for (i=1; i<=N1; i++)
      {
          for (tt=t0; tt<=(t); tt++)
          {
              loss+=gtemp[tt-1][i-1][j-1];
              K_temp_sum+=gtemp[tt-1][i-1][j-1];
          }
      }
    }

    if(loss==N2)
    {
      loss=0;
    }
    else
    {
      loss=min(K_temp_sum-N2,ROUND(loss*shocks_capstock(1)));
    }

    while(loss>0)
    {
      ranj=0;
      while(ranj==0)
      {
        ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      }

      lossj=0;
      for (i=1; i<=N1; i++)
      {
          for (tt=t0; tt<=(t); tt++)
          {
              lossj+=gtemp[tt-1][i-1][ranj-1];
          }
      }
      rnd=ran1(p_seed);
      lossj=min(lossj-1.0,ceil(rnd*lossj));
      lossj=min(loss,lossj);
      K_loss(ranj)=lossj*dim_mach;
      while(lossj>0)
      {
          rani=int(ran1(p_seed)*N1*N2)%N1+1;
          rant=int(ran1(p_seed)*t0*(t))%((t)-t0+1)+t0;
          if (gtemp[rant-1][rani-1][ranj-1]>0)
          {
              if(lossj<gtemp[rant-1][rani-1][ranj-1])
              {
                  if(rant==t)
                  {
                      deltaCapitalStock(1,ranj)-=lossj*g_price[rant-1][rani-1][ranj-1];
                      Loss_Capital_mol(ranj)+=lossj*g_price[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]-=lossj;
                      loss-=lossj;
                      if(SI(ranj)>0)
                      {
                        K(ranj)-=min(SI(ranj),lossj*dim_mach);
                        lossj=max(0.0,lossj*dim_mach-SI(ranj))/dim_mach;
                      }
                      EI(1,ranj)-=lossj*dim_mach;
                      EI(1,ranj)=max(EI(1,ranj),0.0);
                      lossj=0;
                  }
                  else
                  {
                      Loss_Capital(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*lossj);
                      Loss_Capital_mol(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*lossj);
                      CapitalStock(1,ranj)-=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*lossj);
                      K(ranj)-=lossj*dim_mach;
                      K_cur(ranj)-=lossj*dim_mach;
                      gtemp[rant-1][rani-1][ranj-1]-=lossj;
                      g_c2[rant-1][rani-1][ranj-1]-=lossj;
                      g_c3[rant-1][rani-1][ranj-1]-=lossj;
                      loss-=lossj;
                      lossj=0;
                  }
              }
              else
              {
                  if(rant==t)
                  {
                      deltaCapitalStock(1,ranj)-=g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1];
                      Loss_Capital_mol(ranj)+=g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1];
                      if(SI(ranj)>0)
                      {
                        K(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach-EI(1,ranj);
                      }
                      EI(1,ranj)=0;
                      lossj-=gtemp[rant-1][rani-1][ranj-1];
                      loss-=gtemp[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]=0;
                  }
                  else
                  {
                      Loss_Capital(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1]);
                      Loss_Capital_mol(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1]);
                      CapitalStock(1,ranj)-=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1]);
                      K(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach;
                      K_cur(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach;
                      lossj-=gtemp[rant-1][rani-1][ranj-1];
                      loss-=gtemp[rant-1][rani-1][ranj-1];
                      g_c2[rant-1][rani-1][ranj-1]-=gtemp[rant-1][rani-1][ranj-1];
                      g_c3[rant-1][rani-1][ranj-1]-=gtemp[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]=0;
                  }
              }
          }
      }
    }

    if(K_loss.Sum()>0)
    {
      for (j=1; j<=N2; j++)
      {
        if(K_loss(j)>0)
        {
          Q2(j)=min(Q2(j),K_cur(j));
          ADJUSTEMISSENLAB();
        }
      }
    }
  }

  //Calculate average environmental friendliness
  A1_ef_avg = 0;
  A2_ef_avg = 0;
  for (i=1; i <=N1; i++)
  {
    A1_ef_avg += A1p_ef(i);
  }
  for (j=1; j<=N2; j++)
  {
    A2_ef_avg += A2e_ef(j);
  }
  A1_ef_avg = A1_ef_avg/N1;
  A2_ef_avg = A2_ef_avg/N2;

  //Calculate average energy efficiency of C-firms
  A2_en_avg = 0;
  for (j=1; j<=N2; j++)
	{
		A2_en_avg += A2e_en(j);   //Can consider using not effective but A2_en (no shocks and no which machines used)
    if(D2_en_TOT>0){
		}
  }
  A2_en_avg /= N2;

  //Calculate effective carbon intensity of each C-firm and on average
  A2_ci_avg = 0;
  for (j=1; j<=N2; j++)
	{
    A2_ci(j) =  A2e_ef(j) / A2e_en(j);
    A2_ci_avg += A2_ci(j);
  }
  A2_ci_avg /= N2;

}

void ADJUSTEMISSENLAB(void)
{
  nmachprod=ceil(Q2(j)/dim_mach);
  nmp_temp=nmachprod;
  while (nmp_temp > 0)
  {
    cmin=10000000000000000;
    imin=0;
    jmin=0;
    tmin=0;

    for (i=1; i<=N1; i++)
    {
      for (tt=t0; tt<=t; tt++)
      {
        if (g_c2[tt-1][i-1][j-1] > 0 && (w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i)) < cmin)
        {
          cmin=w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i);
          imin=i;
          jmin=j;
          tmin=tt;
        }
      }
    }

    if (nmachprod>0)
    {
      if (g_c2[tmin-1][imin-1][jmin-1] >= nmp_temp)
      {
        A2e2(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*nmp_temp/nmachprod;
        A2e_en2(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*nmp_temp/nmachprod;
        A2e_ef2(j)+=A_ef(tmin,imin)*nmp_temp/nmachprod;
        g_c2[tmin-1][imin-1][jmin-1]-= nmp_temp;
        nmp_temp=0;
      }
      else
      {
        A2e2(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*g_c2[tmin-1][imin-1][jmin-1]/nmachprod;
        A2e_en2(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*g_c2[tmin-1][imin-1][jmin-1]/nmachprod;
        A2e_ef2(j)+=A_ef(tmin,imin)*g_c2[tmin-1][imin-1][jmin-1]/nmachprod;
        nmp_temp-=g_c2[tmin-1][imin-1][jmin-1];
        g_c2[tmin-1][imin-1][jmin-1]=0;
      }
    }
    else
    {
      std::cerr << "\n\n ERROR: nmachprod = 0!!!" << endl;
      Errors << "\n nmachprod = 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
  }

  if(A2e2(j)>0)
  {
    Ld2_control(j)=Q2(j)/A2e2(j);
  }
  else
  {
    if(Q2(j)>0)
    {
      std::cerr << "\n\n ERROR: Q2>0 and A2e2 == 0!!!" << endl;
      Errors << "\n  Q2>0 and A2e2 == 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      Ld2_control(j)=0;
    }
  }

  if(Ld2_control(j)>Ld2_wr(j))
  {
    A2e2(j)=0;
    A2e_en2(j)=0;
    A2e_ef2(j)=0;
    Ldtemp=Ld2_wr(j);
    nmp_temp=0;
    for (i=1; i<=N1; i++)
    {
      for (tt=t0; tt<=t; tt++)
      {
        if (g_c3[tt-1][i-1][j-1] > 0)
        {
          nmp_temp+=g_c3[tt-1][i-1][j-1];
        }
      }
    }
    nmachprod=0;
    while (Ldtemp > 0 && nmachprod<nmp_temp)
    {
      cmin=10000000000000000;
      imin=0;
      jmin=0;
      tmin=0;

      for (i=1; i<=N1; i++)
      {
        for (tt=t0; tt<=t; tt++)
        {
          if (g_c3[tt-1][i-1][j-1] > 0 && (w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i)) < cmin)
          {
            cmin=w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i);
            imin=i;
            jmin=j;
            tmin=tt;
          }
        }
      }

      if ((g_c3[tmin-1][imin-1][jmin-1]*dim_mach/A(tmin,imin)) > Ldtemp)
      {
        nmachprod+=ceil(Ldtemp*A(tmin,imin)/dim_mach);
        A2e2(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*ceil(Ldtemp*A(tmin,imin)/dim_mach);
        A2e_en2(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*ceil(Ldtemp*A(tmin,imin)/dim_mach);
        A2e_ef2(j)+=A_ef(tmin,imin)*ceil(Ldtemp*A(tmin,imin)/dim_mach);
        g_c3[tmin-1][imin-1][jmin-1]-=ceil(Ldtemp*A(tmin,imin)/dim_mach);
        Ldtemp=0;
      }
      else
      {
        Ldtemp-=(g_c3[tmin-1][imin-1][jmin-1]*dim_mach/A(tmin,imin));
        nmachprod+=g_c3[tmin-1][imin-1][jmin-1];
        A2e2(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*g_c3[tmin-1][imin-1][jmin-1];
        A2e_en2(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*g_c3[tmin-1][imin-1][jmin-1];
        A2e_ef2(j)+=A_ef(tmin,imin)*g_c3[tmin-1][imin-1][jmin-1];
        g_c3[tmin-1][imin-1][jmin-1]=0;
      }
    }
    A2e2(j)/=nmachprod;
    Q2(j)=A2e2(j)*Ld2_wr(j);
    A2e_en2(j)/=nmachprod;
    A2e_ef2(j)/=nmachprod;
  }

  if(A2e_en2(j)>0)
  {
    D2_en(j)=Q2(j)/A2e_en2(j);
  }
  else
  {
    if(Q2(j)>0)
    {
      std::cerr << "\n\n ERROR: Q2>0 and A2e_en2 == 0!!!" << endl;
      Errors << "\n  Q2>0 and A2e_en2 == 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      D2_en(j)=0;
    }
  }

  if(A2e_ef2(j)>0 && A2e_en2(j)>0)
  {
    Emiss2(j)= A2e_ef2(j)/A2e_en2(j)*Q2(j);
  }
  else
  {
    if(Q2(j)>0)
    {
      std::cerr << "\n\n ERROR: Q2>0 and A2e_ef2 or A2e_en2 == 0!!!" << endl;
      Errors << "\n  Q2>0 and A2e_ef2 or A2e_en2 == 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      Emiss2(j)=0;
    }
  }
}

void CANCMACH(void)
{
  //Based on actual replacement investment carried out, machines are scrapped
  scrapmax=(scrap_age(j)+SI(j))/dim_mach;
  indforn=int(fornit(j));

  scrap_n=0;
  if(scrapmax>0)
  {
    //First scrap machines which are too old, then ones with high production cost
    for (i=1; i<=N1; i++)
    {
      for (tt=t0; tt<=t; tt++)
      {
        if (scrapmax > 0 && g_pb[tt-1][i-1][j-1] > 0 && age[tt-1][i-1][j-1]>(agemax))
        {
          if (g_pb[tt-1][i-1][j-1] >= scrapmax)
          {
            gtemp[tt-1][i-1][j-1]-=scrapmax;
            g_pb[tt-1][i-1][j-1]-=scrapmax;
            scrap_n+=scrapmax*g_price[tt-1][i-1][j-1];
            scrapmax=0;
          }
          else
          {
            scrapmax-=g_pb[tt-1][i-1][j-1];
            gtemp[tt-1][i-1][j-1]-=g_pb[tt-1][i-1][j-1];
            scrap_n+=g_pb[tt-1][i-1][j-1]*g_price[tt-1][i-1][j-1];
            g_pb[tt-1][i-1][j-1]=0;
          }
        }
      }
    }
  }

  while (scrapmax > 0)
  {
    cmax=0;

    for (i=1; i<=N1; i++)
	  {
	 	  for (tt=t0; tt<=t; tt++)
		  {
        if (g_pb[tt-1][i-1][j-1] > 0 && C_pb[tt-1][i-1][j-1] > cmax)
        {
          ind_i=i;
          ind_tt=tt;
          cmax=C_pb[tt-1][i-1][j-1];
        }
      }
    }

    if (g_pb[ind_tt-1][ind_i-1][j-1] >= scrapmax)
    {
      gtemp[ind_tt-1][ind_i-1][j-1]-=scrapmax;
		  g_pb[ind_tt-1][ind_i-1][j-1]-=scrapmax;
      scrap_n+=scrapmax*g_price[ind_tt-1][ind_i-1][j-1];
		  scrapmax=0;
    }
	  else
	  {
		  scrapmax-=g_pb[ind_tt-1][ind_i-1][j-1];
      scrap_n+=g_pb[ind_tt-1][ind_i-1][j-1]*g_price[ind_tt-1][ind_i-1][j-1];
 		  gtemp[ind_tt-1][ind_i-1][j-1]-=g_pb[ind_tt-1][ind_i-1][j-1];
		  g_pb[ind_tt-1][ind_i-1][j-1]=0;
	  }
  }

  deltaCapitalStock(1,j)+=(SI_n(j)-scrap_n);
}

void EN_DEM_HOUSEHOLDS(void)
{
  c_en_h = c_en(2) * ratio_c_en_h_firms;

  // Calculate households energy demand based on previous period total expenditure
  Expenditure_en_h = 0;
  D_en_h = 0;
  for (const string& cl:classes_mh){
    double exp_en_from_income_mh_cl = energy_expenditure_sh_mh[cl] * Expenditure_tot_mh[cl];
    if (t==t_regime_shifts & flag_change_income_shares==1){ //If change in income, correct by income change to avoid having a bump in energy expenditure share
      Expenditure_en_mh[cl] = Expenditure_en_mh[cl] * (Income_gross_h(2) * income_sh_regime_shift_mh[cl] / Income_gross_mh[cl](2));
    }
    double dem_en_from_income_mh_cl = exp_en_from_income_mh_cl / c_en_h;

    //Persistency in energy demand (D_en_mh[cl] has still previous timestep value)
    if (t==1){
      D_en_mh[cl] = dem_en_from_income_mh_cl;
    }
    D_en_mh[cl] = demand_persistency_h * D_en_mh[cl] + (1 - demand_persistency_h) * dem_en_from_income_mh_cl;
    Expenditure_en_mh[cl] = D_en_mh[cl] * c_en_h;

    Expenditure_en_h += Expenditure_en_mh[cl];
    D_en_h += D_en_mh[cl];
  }

  //Calculate government energy demand
  energy_expenditure_sh_g = umap_sum_over_keys(energy_expenditure_sh_mh) / energy_expenditure_sh_mh.size(); // Unweighted average of households energy expenditure shares
  Exp_en_g = Exp_tot_g * energy_expenditure_sh_g;
  D_en_g = Exp_en_g / c_en_h;
}

void PAY_LAB_INV(void)
{
  //C-firms pay wages
  for (j=1; j<=N2; j++)
	{
    sendingBank=BankingSupplier_2(j);

    Wages_2_i(j)=0;
    for (const string& cl:classes_mh){
      Wages_2_i_mh[cl](j)=w_mh[cl](1)*Ld2_i_mh[cl](j); //Wages paid by firm j to each h class
      Wages_2_i(j)+=Wages_2_i_mh[cl](j); //Add to total wages paid by firm j
    }

    if(Deposits_2(1,j)>=Wages_2_i(j))
    {
      Deposits_2(1,j)-=Wages_2_i(j);
    }
    else
    {
      if(Wages_2_i(j)>0)
      {
        deviation=((Wages_2_i(j)-Deposits_2(1,j))/Wages_2_i(j)); //Should be negative
        if(deviation<=tolerance) //If can pay basically all wages (apart from small tolerance)
        {
          for (const string& cl:classes_mh){ //Scale proportionally wages to each hosehold class
            Wages_2_i_mh[cl](j)*=Deposits_2(1,j)/Wages_2_i(j);
          }
          Wages_2_i(j)=Deposits_2(1,j);
          Deposits_2(1,j)-=Wages_2_i(j);
        }
        else
        {
          std::cerr << "\n\n ERROR: C-Firm " << j << " cannot pay wages in period " << t << endl;
          Errors << "\n C-Firm " << j << " cannot pay wages in period " << t << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    Wages+=Wages_2_i(j);
    for (const string& cl:classes_mh){
      Wages_mh[cl]+=Wages_2_i_mh[cl](j);
    }
    Deposits_b(1,sendingBank)-=Wages_2_i(j);
    Outflows(sendingBank)+=Wages_2_i(j);

    //C-firms pay for investment
    indforn=int(fornit(j));
    Investment_2(j)=EI_n(j)+SI_n(j);
    receivingBank=BankingSupplier_1(indforn);
    if(Deposits_2(1,j)>=Investment_2(j))
    {
      Deposits_2(1,j)-=Investment_2(j);
      Deposits_b(1,sendingBank)-=Investment_2(j);
      Outflows(sendingBank)+=Investment_2(j);

      Deposits_1(1,indforn)+=Investment_2(j);
      Deposits_b(1,receivingBank)+=Investment_2(j);
      Inflows(receivingBank)+=Investment_2(j);
    }
    else
    {
      if(Investment_2(j)>0)
      {
        deviation=((Investment_2(j)-Deposits_2(1,j))/Investment_2(j));
        if(deviation<=tolerance) //Simply correct, matter of rounding
        {
          S1(indforn)-=Investment_2(j);
          S1_pre(indforn)-=Investment_2(j);
          S1_post(indforn)-=Investment_2(j);
          S1(indforn)+=Deposits_2(1,j);
          S1_pre(indforn)+=Deposits_2(1,j);
          S1_post(indforn)+=Deposits_2(1,j);
          Investment_2(j)=Deposits_2(1,j);
          Deposits_2(1,j)-=Investment_2(j);
          Deposits_b(1,sendingBank)-=Investment_2(j);
          Outflows(sendingBank)+=Investment_2(j);

          Deposits_1(1,indforn)+=Investment_2(j);
          Deposits_b(1,receivingBank)+=Investment_2(j);
          Inflows(receivingBank)+=Investment_2(j);
        }
        else
        {
          std::cerr << "\n\n ERROR: C-Firm " << j <<  " cannot pay for investment in period " << t << endl;
          Errors << "\n C-Firm " << j <<  " cannot pay for investment in period " << t << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  Wages_2_mh=umap_sum_each_key(Wages_2_i_mh); //Total wages to each household class
  //K-firms pay wages
  for (i=1; i<=N1; i++)
	{
    sendingBank=BankingSupplier_1(i);

    Wages_1_i(i)=0;
    for (const string& cl:classes_mh){
      Wages_1_i_mh[cl](i)=w_mh[cl](1)*Ld1_i_mh[cl](i); //Wages paid by firm i to each household class
      Wages_1_i(i)+=Wages_1_i_mh[cl](i); //Add to total wages paid by firm i
    }

    if(Deposits_1(1,i)>=Wages_1_i(i))
    {
      Deposits_1(1,i)-=Wages_1_i(i);
    }
    else
    {
      if(Wages_1_i(i)>0)
      {
        deviation=((Wages_1_i(i)-Deposits_1(1,i))/Wages_1_i(i));
        for (const string& cl:classes_mh){ //Scale proportionally wages to each hosehold class
          Wages_1_i_mh[cl](i)*=Deposits_1(1,i)/Wages_1_i(i);
        }
        Wages_1_i(i)=Deposits_1(1,i); //Pay wages in any case
        Deposits_1(1,i)-=Wages_1_i(i);
        if(deviation>tolerance) //If not able however exit
        {
          exiting_1(i)=1;
          exiting_1_payments(i)=1;
        }
      }
    }
    Wages+=Wages_1_i(i);
    for (const string& cl:classes_mh){
      Wages_mh[cl]+=Wages_1_i_mh[cl](i);
    }
    Deposits_b(1,sendingBank)-=Wages_1_i(i);
    Outflows(sendingBank)+=Wages_1_i(i);
  }

  Wages_1_mh=umap_sum_each_key(Wages_1_i_mh); //Total wages to each household class


  //Energy sector pays wages
  if(Deposits_e(1)>=Wages_en)
  {
    Deposits_e(1)-=Wages_en;
  }
  else
  {
    deviation=((Wages_en-Deposits_e(1))/Wages_en);
    if(deviation<=tolerance) //If can pay basically all wages (apart from small tolerance)
    {
      for (const string& cl:classes_mh){ //Scale proportionally wages to each hosehold class
      Wages_en_mh[cl]*=Deposits_e(1)/Wages_en;
    }
    Wages_en=Deposits_e(1);
      Deposits_e(1)-=Wages_en;
      }
    else
    {
      std::cerr << "\n\n ERROR: Energy sector cannot pay wages in period " << t << endl;
      Errors << "\n Energy sector cannot pay wages in period " << t << endl;
      exit(EXIT_FAILURE);
    }
  }
  Wages+=Wages_en;
  for (const string& cl:classes_mh){
    Wages_mh[cl]+=Wages_en_mh[cl];
  }

  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)-=Wages_en*DepositShare_e(i);
    Outflows(i)+=Wages_en*DepositShare_e(i);
    Deposits_eb(1,i)-=Wages_en*DepositShare_e(i);
  }

  for (const string& cl:classes_mh){
    Deposits_mh[cl](1)+=(Wages_mh[cl]+Benefits_mh[cl]+govTransfers_mh[cl]);
  }
  Deposits_h(1)+=(Wages+Benefits+govTransfers); //Here govTransfer still does not include Transfers from energy shock, added later

  //Households receive wages, Benefits & government transfers
  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)+=(Wages+Benefits+govTransfers)*DepositShare_h(i);
    Inflows(i)+=(Wages+Benefits+govTransfers)*DepositShare_h(i);
    Deposits_hb(1,i)+=(Wages+Benefits+govTransfers)*DepositShare_h(i);
  }
}

void COMPET2(void)
{
  //C-firms' market shares are updated using quasi-replicator dynamics. Firms with too low market share exit
  p2m(1)=0;
  l2m=0;
  mu2m=0;

  for (j=1; j<=N2; j++)
  {
    if (exiting_2(j) == 0)
    // Exclude exiting firms from the computation of average C-firm competitiveness.
    // They are removed from the average price and the average unsatisfied demand (see below).
    // They are de facto removed from ftot since f2 is set to zero when they exit
    // Their competitiveness is calculated but not included in the average C-firm competitiveness since their market share f2 is set to 0 when they exit
    {
      p2m(1)+=p2(j)/(N2r-exiting_2.Sum());
      l2m+=l2(j)/(N2r-exiting_2.Sum());
      mu2m+=mu2(1,j)/(N2r-exiting_2.Sum());
    }
  }

	ftot=0;
  double old_Em2 = 0; //
  RowVector old_E2; //
  RowVector old_f2; //
  old_E2.ReSize(N2); //
  old_f2.ReSize(N2); //
  double old_ftot = 0; //

  for (j=1; j<=N2; j++)
  {
    //Calculate competiveness
    E2(j) = -pow(p2(j)/p2m(1),omega1) - pow(l2(j)/l2m,omega2) - flag_env_friendl_in_market_shares * pow(A2_ci(j)/A2_ci_avg, omega4);
    old_E2(j) = -pow(p2(j)/p2m(1),omega1) - pow(l2(j)/l2m,omega2); //
    ftot(2) += f2(2,j);
	}

	if (ftot(2)>0)
	{
		for (j=1; j<=N2; j++)
		{
      //Calculate average competitiveness based on previous period market shares
			f2(2,j)/=ftot(2);
			Em2(1)+=E2(j)*f2(2,j);

      old_Em2 += old_E2(j)*f2(2,j); //
		}
	}
	else
	{
		std::cerr << "\n\n ERROR: f2 = 0 for all firms in period " << t << endl;
    Errors << "\n f2 = 0 for all firms in period " << t << endl;
    exit(EXIT_FAILURE);
	}

	ftot=0;
	if (Em2(1)==0)
  {
		std::cerr << "\n\n ERROR: Em2(1)=0 in period " << t << endl;
    Errors << "\n Em2(1)=0 in period " << t << endl;
    exit(EXIT_FAILURE);
  }


  // std::cout << "\n <t>imestep" << t; //
  
	for (j=1; j<=N2; j++)
	{
    old_f2(j) = f2(2,j)*((2*omega3)/(1+exp((-chi)*((old_E2(j)-old_Em2)/old_Em2)))+(1-omega3)); //
    f2(1,j) = f2(2,j)*((2*omega3)/(1+exp((-chi)*((E2(j)-Em2(1))/Em2(1))))+(1-omega3));
    
    //If market share too low, exit
    if (f2(1,j) <= (1/(N2r*500)))
    {
      f2(1,j)=0;
      f2(2,j)=0;
      f2(3,j)=0;
      if(exiting_2(j)==0 && exit_payments2(j)==0)
      {
        exit_marketshare2(j)=1;
      }
    }
    

		ftot(1)+=f2(1,j);
		ftot(2)+=f2(2,j);
		ftot(3)+=f2(3,j);

    if (old_f2(j) <= (1/(N2r*500))){ //
      old_f2(j) = 0; //
    } //
    old_ftot += old_f2(j); //
      
	}

	if (ftot(1)==0 || ftot(2)==0 || ftot(3)==0)
  {
		std::cerr << "\n\n ERROR: ftot=0 in period " << t << endl;
    Errors << "\n ftot=0 in period " << t << endl;
    exit(EXIT_FAILURE);
  }

	for (j=1; j<=N2; j++)
	{
    //Normalise to sum to 1
		f2(1,j)/=ftot(1);
		f2(2,j)/=ftot(2);
		f2(3,j)/=ftot(3);
    old_f2(j) /= old_ftot; //
    // std::cout << "\n Firm" << j-1 << " "  << A2_ci(j)/A2_ci_avg << " " << (f2(1,j)-f2(2,j))/(old_f2(j)-f2(2,j)) << " " << (old_E2(j)-old_Em2)/old_Em2 << " " << (E2(j)-Em2(1))/Em2(1);// << " " << pow(A2_ci(j)/A2_ci_avg, omega4); //
    // std::cout << "\n Firm" << j-1 << " "  << (-pow(p2(j)/p2m(1),omega1))/E2(j) << " " << (-pow(l2(j)/l2m,omega2))/E2(j) << " " << (-pow(A2_ci(j)/A2_ci_avg, omega4))/E2(j);// << " " << pow(A2_ci(j)/A2_ci_avg, omega4);
    // std::cout << "\n Firm" << j-1 << " "  << E2(j) << " " << (E2(j)-Em2(1))/Em2(1)<< " " << f2(1,j) << " " << old_E2(j) << " " << old_f2(j);// << " " << pow(A2_ci(j)/A2_ci_avg, omega4);

  }
}

void PROFIT(void)
{
  //Calculate K-firm profit
  for (i=1; i<=N1; i++)
	{
    kpi+=Q1(i)/Q1.Sum()*p1(i);
    EnergyPayments_1(i)=c_en(1)*D1_en(i);
    Pi1(1,i)=S1(i)+InterestDeposits_1(i)+Transfer_shock_f1(i)-Wages_1_i(i)-EnergyPayments_1(i)-t_CO2*Emiss1(i);
    sendingBank=BankingSupplier_1(i);

    //K-firms pay for energy; exit if they are unable
    if(Deposits_1(1,i)>=EnergyPayments_1(i))
    {
      Deposits_1(1,i)-=EnergyPayments_1(i);
      EnergyPayments+=EnergyPayments_1(i);
      Deposits_b(1,sendingBank)-=EnergyPayments_1(i);
      Outflows(sendingBank)+=EnergyPayments_1(i);
    }
    else if (Deposits_1(1,i)>=0)
    {
      deviation=((EnergyPayments_1(i)-Deposits_1(1,i))/EnergyPayments_1(i));
      EnergyPayments_1(i)=Deposits_1(1,i);
      Deposits_1(1,i)-=EnergyPayments_1(i);
      EnergyPayments+=EnergyPayments_1(i);
      Deposits_b(1,sendingBank)-=EnergyPayments_1(i);
      Outflows(sendingBank)+=EnergyPayments_1(i);
      if(deviation>tolerance)
      {
        exiting_1(i)=1;
        exiting_1_payments(i)=1;
      }
    }

    //K-firms pay CO2 tax
    if(exiting_1(i)==0)
    {
      Taxes_CO2_1(i)=t_CO2*Emiss1(i);
      if(Deposits_1(1,i)>=Taxes_CO2_1(i))
      {
        Taxes_CO2(1)+=Taxes_CO2_1(i);
        Deposits_1(1,i)-=Taxes_CO2_1(i);
        Deposits_b(1,sendingBank)-=Taxes_CO2_1(i);
        Outflows(sendingBank)+=Taxes_CO2_1(i);
      }
      else if(Deposits_1(1,i)>=0)
      {
        Taxes_CO2_1(i)=Deposits_1(1,i);
        Deposits_1(1,i)-=Taxes_CO2_1(i);
        Taxes_CO2(1)+=Taxes_CO2_1(i);
        Deposits_b(1,sendingBank)-=Taxes_CO2_1(i);
        Outflows(sendingBank)+=Taxes_CO2_1(i);
      }
    }
    else
    {
      Taxes_CO2_1(i)=0;
    }

    if (Pi1(1,i) > 0 && exiting_1(i)==0)
		{
			//If profit is positive, pay tax
      if(Deposits_1(1,i)>=aliq*Pi1(1,i))
      {
        Taxes_1(i)=aliq*Pi1(1,i);
        Deposits_1(1,i)-=Taxes_1(i);
        Deposits_b(1,sendingBank)-=Taxes_1(i);
        Outflows(sendingBank)+=Taxes_1(i);
        Taxes_g+=Taxes_1(i);
      }
      else if (Deposits_1(1,i)>=0)
      {
        Taxes_1(i)=Deposits_1(1,i);
        Deposits_1(1,i)-=Taxes_1(i);
        Deposits_b(1,sendingBank)-=Taxes_1(i);
        Outflows(sendingBank)+=Taxes_1(i);
        Taxes_g+=Taxes_1(i);
      }
      
      //Calculate bonuses supposed to be paid (if deposits sufficient)
      double Pi1_diff = Pi1(1,i) - Pi1(2,i);
      if(Pi1_diff > 0 & Pi1(2,i) > 0) //If profit increased and in previous period were positive, pay bonus as fraction of profits increase
      {
        Bonuses_1_i(i)=bonuses_share*Pi1_diff;
      }
      else if(Pi1_diff > 0 & Pi1(1,i) > 0)  //If profit increased and in previous period were negative but in current positive, pay bonus as fraction of new positive profits
      {
        Bonuses_1_i(i)=bonuses_share*Pi1(1,i);
      }
      else //If profit decreased, or are negative, no bonuses
      {
        Bonuses_1_i(i)=0;
      }

      //Pay bonuses and dividends, based on dividends availability
      if(Deposits_1(1,i)>=(Bonuses_1_i(i)+d1*(Pi1(1,i)-Bonuses_1_i(i))))
      {
        Dividends_1_i(i)=d1*(Pi1(1,i)-Bonuses_1_i(i));
      }
      else if (Deposits_1(1,i)>=Bonuses_1_i(i)) //If can afford to pay more than bonuses, pay bonuses completely and dividends as residuals
      {
        Dividends_1_i(i)=Deposits_1(1,i)-Bonuses_1_i(i);
      }
      else if (Deposits_1(1,i)>=0) //If cannot afford to pay bonuses entirely, pay only bonuses as much as possible and no dividends
      {
        Bonuses_1_i(i)=Deposits_1(1,i);
        Dividends_1_i(i)=0;
      }
      else //Deposits_1(1,i)<0
      {
        std::cerr << "\n\n ERROR: K-firm " << i << " has negative deposits in period " << t << endl;
        Errors << "\n K-firm " << i << " has negative deposits in period " << t << endl;
        exit(EXIT_FAILURE);
      }

      Deposits_1(1,i)-=(Dividends_1_i(i)+Bonuses_1_i(i));
      Deposits_b(1,sendingBank)-=(Dividends_1_i(i)+Bonuses_1_i(i));
      Outflows(sendingBank)+=(Dividends_1_i(i)+Bonuses_1_i(i));
			Dividends_1+=Dividends_1_i(i);
      Dividends_h(1)+=Dividends_1_i(i);
      Bonuses_1+=Bonuses_1_i(i);
      Bonuses_h+=Bonuses_1_i(i);
      Bonuses_1_mh["ma"]+=Bonuses_1_i(i);
      Bonuses_mh["ma"](1)+=Bonuses_1_i(i);
      for (const string& cl:classes_mh){ //Distribute to household classes
        Dividends_1_i_mh[cl](i)=Dividends_1_i(i) * Ownership_sh_1_i_mh[cl](i); //Dividends split based on firm ownership shares of each household class
        Dividends_1_mh[cl]+=Dividends_1_i_mh[cl](i);
        Dividends_mh[cl](1)+=Dividends_1_i_mh[cl](i);
      }
		}

		Pitot1+=Pi1(1,i);

    //Calculate mean deposits of non-failing K-firms for entry below
		if (exiting_1(i)==0)
		{
			mD1+=Deposits_1(1,i);
			ns1++;
		}

    if(nclient(i) < 1 && exiting_1(i)==0)
    {
      exiting_1(i)=1;
    }

    //Prepare failing K-firms for exit
    if(exiting_1(i)==1)
    {
      //Deposits recovered from failing K-firms overall and by each household class
      Deposits_recovered_1+=Deposits_1(1,i);
      for (const string& cl:classes_mh){
        Deposits_recovered_1_mh[cl] += Deposits_1(1,i) * Ownership_sh_1_i_mh[cl](i);
      }
      sendingBank=BankingSupplier_1(i);
      baddebt_1(i)=-Deposits_1(1,i);
      Deposits_b(1,sendingBank)-=Deposits_1(1,i);
      Outflows(sendingBank)+=Deposits_1(1,i);
      Deposits_1(1,i)=0;
    }
  }

  //Households receive K-firm dividends and bonuses and deposits recovered and pay tax on all wages and on K-firms dividends and bonuses
  Taxes_w_h = 0;
  Taxes_div_1_h = 0;
  Taxes_bon_1_h = 0;
  for (const string& cl : classes_mh){
    //Taxes on wages
    double taxes_w_cl=aliqw_mh[cl]*Wages_mh[cl];
    Taxes_w_mh[cl]=taxes_w_cl; //Each household class
    Taxes_w_h+=taxes_w_cl;     //Overall
    //Taxes on K-firms dividends
    double taxes_div_cl=aliqdiv*Dividends_1_mh[cl];
    Taxes_div_1_mh[cl]=taxes_div_cl; //Each household class
    Taxes_div_1_h+=taxes_div_cl;     //Overall
    //Taxes on K-firms bonuses
    double taxes_bon_cl=aliqw_mh[cl]*Bonuses_1_mh[cl];
    Taxes_bon_1_mh[cl]=taxes_bon_cl; //Each household class
    Taxes_bon_1_h+=taxes_bon_cl;     //Overall

    Deposits_mh[cl](1)+=(Deposits_recovered_1_mh[cl]+Bonuses_1_mh[cl]+Dividends_1_mh[cl]-Taxes_w_mh[cl]-Taxes_div_1_mh[cl]-Taxes_bon_1_mh[cl]);
  }

  double delta_Deposits_h = Deposits_recovered_1+Dividends_1+Bonuses_1-Taxes_w_h-Taxes_div_1_h-Taxes_bon_1_h;
  Deposits_h(1)+=delta_Deposits_h;
  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)+=delta_Deposits_h*DepositShare_h(i);
    Inflows(i)+=delta_Deposits_h*DepositShare_h(i);
    Deposits_hb(1,i)+=delta_Deposits_h*DepositShare_h(i);
  }

  if(flag_energyshocks>0 && t>t_regime_shifts){
    //Transfer to households
    if(flag_energyshocks==1 || flag_energyshocks==2){
      Transfer_shock=tre*max(0.0,(c_en(1)-c_en_preshock)*(D2_en.Sum()+D1_en.Sum())+ratio_c_en_h_firms*(c_en(1)-c_en_preshock)*(D_en_h+D_en_g));
    } else if(flag_energyshocks==3 || flag_energyshocks==4){
      //Transfer equal to overprofits of both fossil fuel and energy sectors
      Transfer_shock=tre*(max(0.0,FuelCost-(FuelCost/pf)*pf_preshock)+max(0.0,(max(0.0,EnergyPayments-(EnergyPayments/c_en(1))*(c_en_preshock)))-(max(0.0,FuelCost-(FuelCost/pf)*pf_preshock))));
    }

    for (const string& cl:classes_mh)
    {
      Transfer_shock_mh[cl]=Transfer_shock*Transfer_shock_sh_mh[cl];
      govTransfers_mh[cl]+=Transfer_shock_mh[cl];
      Deposits_mh[cl]+=Transfer_shock_mh[cl];
    }
    govTransfers+=Transfer_shock;
    Deposits_h(1)+=Transfer_shock;
    for(i=1; i<=NB; i++)
    {
      Deposits_b(1,i)+=(Transfer_shock)*DepositShare_h(i);
      Inflows(i)+=(Transfer_shock)*DepositShare_h(i);
      Deposits_hb(1,i)+=(Transfer_shock)*DepositShare_h(i);
    }
  }

  //Determine consumption demand
    //Households expenditure and C-goods consumption demand is determined; Check if consumption can be financed
  Exp_li=0;
  Exp_ki=0;
  Exp_dep=0;
  Exp_u_benefit=0;
  Exp_govTransf=0;
  Expenditure_tot_h=0;
  Cons_h=0;
  Consumption_h=0;
  for (const string& cl:classes_mh){
    Exp_li_mh[cl] = a1_mh[cl]*(Wages_mh[cl]-Taxes_w_mh[cl]+Bonuses_mh[cl](2)-Taxes_bon_mh[cl](2));
    Exp_ki_mh[cl] = a2_mh[cl]*(InterestDeposits_mh[cl]+Dividends_mh[cl](2)-Taxes_div_mh[cl](2));
    Exp_dep_mh[cl] = a3_mh[cl]*Deposits_mh[cl](2);
    Exp_govTransf_mh[cl] = a4_mh[cl]*govTransfers_mh[cl];
    Exp_u_benefit_mh[cl] = a5_mh[cl]*Benefits_mh[cl];
    Expenditure_tot_mh[cl] = (1-shock_cons)*(Exp_li_mh[cl]+Exp_ki_mh[cl]+Exp_dep_mh[cl]+Exp_govTransf_mh[cl]+Exp_u_benefit_mh[cl]);
    if(Deposits_mh[cl](1)<Expenditure_tot_mh[cl]) //Max target expenditure = households deposits
    {
      Expenditure_tot_mh[cl] = Deposits_mh[cl](1);
    }
    Exp_li += Exp_li_mh[cl];
    Exp_ki += Exp_ki_mh[cl];
    Exp_dep += Exp_dep_mh[cl];
    Exp_govTransf += Exp_govTransf_mh[cl];
    Exp_u_benefit += Exp_u_benefit_mh[cl];

    //Consumption_mh[cl] is the same as in previous period since not re-set to 0
    double real_consumption_based_on_income_cl = Expenditure_tot_mh[cl] * (1 - energy_expenditure_sh_mh[cl]) / p2m(1);
    if (t==1){ //Initialise variables that were not already at t==1
      Consumption_mh[cl] = real_consumption_based_on_income_cl;
      //p2m(2)=p2m(1);
    }
    double previous_period_real_cons = Consumption_mh[cl] / cpi(2);
    double cons_mh_real = demand_persistency_h * previous_period_real_cons + (1 - demand_persistency_h) * (real_consumption_based_on_income_cl);
    Cons_mh[cl] = cons_mh_real * p2m(1); //C-goods consumption
    Cons_h += Cons_mh[cl];

    //Need to recalculate the target total expenditure (then adjusted if Cons_mh is above what remains of the deposits)
    Expenditure_tot_mh[cl] = Cons_mh[cl] + Expenditure_en_mh[cl];
    Expenditure_tot_h += Expenditure_tot_mh[cl]; //Expenditure overall

    //Initialise actual consumption to desired one (so that we can calculate reduction based on Cons_mh)
    Consumption_mh[cl] = Cons_mh[cl];
    Consumption_h += Consumption_mh[cl];
  }
  
  //Households pay for energy and if deposits are not enough reduce consumption
  for (const string& cl:classes_mh){
    if ((Deposits_mh[cl](1)-Expenditure_en_mh[cl])/Deposits_mh[cl](1) < tolerance){ //Case in which all deposits are consumed, to avoid negativitity for subtracting first C-consumption
      Deposits_mh[cl](1) = 0;
    } else if (Deposits_mh[cl](1) >= Expenditure_en_mh[cl]){
      Deposits_mh[cl](1) -= Expenditure_en_mh[cl];
    } else{
      std::cerr << "\n\n ERROR: Household class " << cl <<  " cannot pay for energy in period " << t << endl;
      Errors << "\n Household class " << cl <<  " cannot pay for energy in period " << t << endl;
      exit(EXIT_FAILURE);
    }

    //Reduce actual consumption in case higher than residual deposits
    if(Deposits_mh[cl](1)<Consumption_mh[cl]){
      Consumption_h-=(Consumption_mh[cl]-Deposits_mh[cl](1));
      Consumption_mh[cl]=Deposits_mh[cl](1);
    }
  }
  Deposits_h(1) -= Expenditure_en_h;
  EnergyPayments += Expenditure_en_h;

  //Calculate C-expenditure shares (to updated later in case cannot be satisfied)
  for (const string& cl:classes_mh){
    Cons_sh_mh[cl] = Consumption_mh[cl] / Consumption_h;
  }

    //Goverment
  if (t==1){ //since GDP not yet calculated here in the first step
    Exp_tot_g = Consumption_h * (g_cons_GDP_ratio+0.3) / (1-g_cons_GDP_ratio);
  } else{
    Exp_tot_g = g_cons_GDP_ratio * GDP_n_new(2);
  }
      //Eventually add tax revenues to public spending
      //ISSUE TO SOLVE: if tax rates goes from a value !0 0 to 0 (else cases) we need to change to re-calculate tax rate
  if (t>=t_regime_shifts && flag_change_public_spending==2 && (flag_change_wage_tax_rate!=0 || flag_change_dividends_tax_rate!=0 || flag_change_wealth_tax_rate!=0)){
    double delta_div_taxes;
    if (aliqdiv!=0){
      delta_div_taxes = Taxes_div_h - Taxes_div_h * aliqdiv_prechange / aliqdiv;
    }else{
      delta_div_taxes = 0; 
    }
    double delta_wage_taxes = 0;
    double delta_wealth_tax_rate = 0;
    for (const string& cl:classes_mh){
      if (aliqw_mh[cl]!=0){
        delta_wage_taxes += Taxes_w_mh[cl] - Taxes_w_mh[cl] * aliqw_mh_prechange[cl] / aliqw_mh[cl];
      } else {
        delta_wage_taxes += 0;
      }
      if (aliqwealth_mh[cl]!=0){
        delta_wealth_tax_rate += Taxes_wealth_mh[cl] - Taxes_wealth_mh[cl] * aliqwealth_mh_prechange[cl] / aliqwealth_mh[cl];
      } else {
        delta_wealth_tax_rate += 0;
      }
    }

    Exp_tot_g += delta_div_taxes + delta_wage_taxes + delta_wealth_tax_rate;
  }
  Cons_g = Exp_tot_g - Exp_en_g;
      //Total
  Cons = Cons_g + Consumption_h;

  //Calculate CPI
  for (j=1; j<=N2; j++)	// Average C-firms price weighted on market shares
	{
    cpi(1)+=p2(j)*f2(1,j);
  }

  if (cpi(1) < 0.01)
	{
		std::cerr << "\n\n ERROR: CPI < 0.01 in period " << t << endl;
    Errors << "\n CPI < 0.01 in period " << t << endl;
    exit(EXIT_FAILURE);
	}

  //Consumption of C-goods takes place
  ALLOC();

  //Government determines and pays environmental subsidies for non-exiting C-firms. NB: EI and SI are expressed in terms of productive capacity, not number of machines, hence need to be divided by dim_mach
	for (j=1; j<=N2; j++)
	{
    if ((flag_environmental_subsidies_C_firms==1 || flag_environmental_subsidies_C_firms==2) && t>=t_regime_shifts && exiting_2(j)==0){
      EnvSubsidies_2_i(j) = EI(1,j)/dim_mach * EnvSubsidiesPerMachine_1_i(indforn) + SI(j)/dim_mach * EnvSubsidiesPerMachine_1_i(indforn); //Value of subsidy per machine is calculated in previous period. EnvSubsidies_2_i(j) set to 0 in SETVARS() in case no subsidies or exiting
      EnvSubsidies_2 += EnvSubsidies_2_i(j); //EnvSubsidies_2 set to 0 in SETVARS()

      if (env_subsidy_per_machine > 0)
      { //Subsidies
        Deposits_2(1,j)+=EnvSubsidies_2_i(j);
        Deposits_b(1,BankingSupplier_2(j))+=EnvSubsidies_2_i(j);
        Inflows(BankingSupplier_2(j))+=EnvSubsidies_2_i(j);
      } 
      else
      { //If taxes, need to check that firms can pay them out of deposits
        if(Deposits_2(1,j)>=fabs(EnvSubsidies_2_i(j)))
        {
          Deposits_2(1,j)+=EnvSubsidies_2_i(j);
          Deposits_b(1,BankingSupplier_2(j))+=EnvSubsidies_2_i(j);
          Outflows(BankingSupplier_2(j))-=EnvSubsidies_2_i(j); //- since it has a negative sign if it's a tax
        }
        else if(Deposits_2(1,j)>=0)
        {
          EnvSubsidies_2-=EnvSubsidies_2_i(j);
          EnvSubsidies_2_i(j)=-Deposits_2(1,j);
          Deposits_2(1,j)=0;
          EnvSubsidies_2+=EnvSubsidies_2_i(j);
          Deposits_b(1,BankingSupplier_2(j))+=EnvSubsidies_2_i(j);
          Outflows(BankingSupplier_2(j))-=EnvSubsidies_2_i(j);
        }
        else{
          std::cerr << "\n\n ERROR: C-firm " << j << " has negative deposits in period " << t << endl;
          Errors << "\n C-firm " << j << " has negative deposits in period " << t << endl;
          exit(EXIT_FAILURE);          
        }
      }
    }
  }

  //Calculate C-firm profits
	for (j=1; j<=N2; j++)
	{
		dN(j)=N(1,j)-N(2,j);
    Inventories(1,j)=N(1,j)*p2(j);
    dNm(j)=Inventories(1,j)-Inventories(2,j);
    dNtot+=dN(j);
    dNmtot+=dNm(j);
    EnergyPayments_2(j)=c_en(1)*D2_en(j);
    mol(j)=S2(1,j)-Wages_2_i(j)-EnergyPayments_2(j)-Loss_Capital_mol(j);
    if(exiting_2(j)==0)
    {
      LoanInterest_2(j)=r_deb_h(j)*Loans_2(1,j);
    }
    else
    {
      LoanInterest_2(j)=0;
    }

    //NB: Profit here also includes changes in nominal value of tangible assets
    Pi2(1,j)=S2(1,j)+InterestDeposits_2(j)+dNm(j)+deltaCapitalStock(1,j)+Transfer_shock_f2(j)+EnvSubsidies_2_i(j)-Investment_2(j)-Wages_2_i(j)-EnergyPayments_2(j)-LoanInterest_2(j)-t_CO2*Emiss2(j);

    //C-firm pays energy; if unable to do so it exits
    if(Deposits_2(1,j)>=EnergyPayments_2(j))
    {
      sendingBank=BankingSupplier_2(j);
      Deposits_2(1,j)-=EnergyPayments_2(j);
      EnergyPayments+=EnergyPayments_2(j);
      Deposits_b(1,sendingBank)-=EnergyPayments_2(j);
      Outflows(sendingBank)+=EnergyPayments_2(j);
    }
    else if(Deposits_2(1,j)>=0)
    {
      sendingBank=BankingSupplier_2(j);
      EnergyPayments_2(j)=Deposits_2(1,j);
      Deposits_2(1,j)=0;
      Deposits_b(1,sendingBank)-=EnergyPayments_2(j);
      Outflows(sendingBank)+=EnergyPayments_2(j);
      EnergyPayments+=EnergyPayments_2(j);
      Taxes_2(j)=0;
      LoanInterest_2(j)=0;
      DebtRemittances2(j)=0;
      exiting_2(j)=1;
      exit_payments2(j)=1;
      Pi2(1,j)=S2(1,j)+InterestDeposits_2(j)+dNm(j)+deltaCapitalStock(1,j)+Transfer_shock_f2(j)+EnvSubsidies_2_i(j)-Investment_2(j)-Wages_2_i(j)-EnergyPayments_2(j)-LoanInterest_2(j)-t_CO2*Emiss2(j);
    }

    //Non-exiting firms make principal & interest payments on loans; if unable to do so they exit
    if(exiting_2(j)==0)
    {
      DebtRemittances2(j)=repayment_share*Loans_2(1,j);
      if(Deposits_2(1,j)>=(LoanInterest_2(j)+DebtRemittances2(j)))
      {
        receivingBank=BankingSupplier_2(j);
        Deposits_2(1,j)-=(LoanInterest_2(j)+DebtRemittances2(j));
        Deposits_b(1,receivingBank)-=(LoanInterest_2(j)+DebtRemittances2(j));
        DebtService_2(1,j)=LoanInterest_2(j)+DebtRemittances2(j);
        Loans_2(1,j)-=DebtRemittances2(j);
        Loans_b(1,receivingBank)-=DebtRemittances2(j);
        LoanInterest(receivingBank)+=LoanInterest_2(j);
      }
      else
      {
        LoanInterest_2(j)=0;
        DebtRemittances2(j)=0;
        Taxes_2(j)=0;
        exiting_2(j)=1;
        exit_payments2(j)=1;
        Pi2(1,j)=S2(1,j)+InterestDeposits_2(j)+dNm(j)+deltaCapitalStock(1,j)+Transfer_shock_f2(j)+EnvSubsidies_2_i(j)-Investment_2(j)-Wages_2_i(j)-EnergyPayments_2(j)-LoanInterest_2(j)-t_CO2*Emiss2(j);
      }
    }

    //C-firms pay CO2 tax; if unable to do so they exit
    if(exiting_2(j)==0)
    {
      Taxes_CO2_2(j)=t_CO2*Emiss2(j);
      sendingBank=BankingSupplier_2(j);
      if(Deposits_2(1,j)>=Taxes_CO2_2(j))
      {
        Taxes_CO2(1)+=Taxes_CO2_2(j);
        Deposits_2(1,j)-=Taxes_CO2_2(j);
        Deposits_b(1,sendingBank)-=Taxes_CO2_2(j);
        Outflows(sendingBank)+=Taxes_CO2_2(j);
      }
      else if(Deposits_2(1,j)>=0)
      {
        Taxes_CO2_2(j)=Deposits_2(1,j);
        Deposits_2(1,j)-=Taxes_CO2_2(j);
        Taxes_CO2(1)+=Taxes_CO2_2(j);
        Deposits_b(1,sendingBank)-=Taxes_CO2_2(j);
        Outflows(sendingBank)+=Taxes_CO2_2(j);
      }
    }
    else
    {
      Taxes_CO2_2(j)=0;
    }

    //C-firm pays profit tax; if unable to do so they exit
    if((Pi2(1,j)-EnvSubsidies_2_i(j))>0 && exiting_2(j)==0)
    {
      Taxes_2(j)=aliq*(Pi2(1,j)-EnvSubsidies_2_i(j));
      sendingBank=BankingSupplier_2(j);
      if(Deposits_2(1,j)>=Taxes_2(j))
      {
        Deposits_2(1,j)-=Taxes_2(j);
        Taxes_g+=Taxes_2(j);
        Deposits_b(1,sendingBank)-=Taxes_2(j);
        Outflows(sendingBank)+=Taxes_2(j);
      }
      else if(Deposits_2(1,j)>=0)
      {
        Taxes_2(j)=Deposits_2(1,j);
        Deposits_2(1,j)=0;
        Taxes_g+=Taxes_2(j);
        Deposits_b(1,sendingBank)-=Taxes_2(j);
        Outflows(sendingBank)+=Taxes_2(j);
      }
    }

    if(exiting_2(j)==0){
      //Calculate bonuses supposed to be paid (if deposits sufficient)
      double Pi2_diff = Pi2(1,j) - Pi2(2,j);
      if(Pi2_diff > 0 & Pi2(2,j) > 0) //If profit increased and in previous period were positive, pay bonus as fraction of profits increase
      {
        Bonuses_2_i(j)=bonuses_share*Pi2_diff;
      }
      else if(Pi2_diff > 0 & Pi2(1,j) > 0)  //If profit increased and in previous period were negative but in current positive, pay bonus as fraction of new positive profits
      {
        Bonuses_2_i(j)=bonuses_share*Pi2(1,j);
      }
      else //If profit decreased, or are negative, no bonuses
      {
        Bonuses_2_i(j)=0;
      }

      //Non-exiting firms pay dividends
      if (Pi2(1,j) > 0)
      {
        sendingBank=BankingSupplier_2(j);
        if(Deposits_2(1,j)>=(Bonuses_2_i(j)+d2*(Pi2(1,j)-Bonuses_2_i(j))))
        {
          Dividends_2_i(j)=d2*(Pi2(1,j)-Bonuses_2_i(j));
        }
        else if (Deposits_2(1,j)>=Bonuses_2_i(j)) //If can afford to pay more than bonuses, pay bonuses completely and dividends as residuals
        {
          Dividends_2_i(j)=Deposits_2(1,j)-Bonuses_2_i(j);
        }
        else if (Deposits_2(1,j)>=0) //If cannot afford to pay bonuses entirely, pay only bonuses as much as possible and no dividends
        {
          Bonuses_2_i(j)=Deposits_2(1,j);
          Dividends_2_i(j)=0;
        }
        else //Deposits_2(1,j)<0
        {
          std::cerr << "\n\n ERROR: C-firm " << j << " has negative deposits in period " << t << endl;
          Errors << "\n C-firm " << j << " has negative deposits in period " << t << endl;
          exit(EXIT_FAILURE);
        }
        Deposits_2(1,j)-=(Dividends_2_i(j)+Bonuses_2_i(j));
        Deposits_b(1,sendingBank)-=(Dividends_2_i(j)+Bonuses_2_i(j));
        Outflows(sendingBank)+=(Dividends_2_i(j)+Bonuses_2_i(j));
        Dividends_h(1)+=Dividends_2_i(j);
        Dividends_2+=Dividends_2_i(j);
        Bonuses_2+=Bonuses_2_i(j);
        Bonuses_h+=Bonuses_2_i(j);
        Bonuses_2_mh["ma"]+=Bonuses_2_i(j);
        Bonuses_mh["ma"](1)+=Bonuses_2_i(j);
        for (const string& cl:classes_mh){ //Distribute to household classes
          Dividends_2_i_mh[cl](j)=Dividends_2_i(j) * Ownership_sh_2_i_mh[cl](j);
          Dividends_2_mh[cl]+=Dividends_2_i_mh[cl](j);
          Dividends_mh[cl](1)+=Dividends_2_i_mh[cl](j);
        }
      }
    }
    
		Pitot2+=Pi2(1,j);
	}

  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)-=Expenditure_en_h*DepositShare_h(i);
    Outflows(i)+=Expenditure_en_h*DepositShare_h(i);
    Deposits_hb(1,i)-=Expenditure_en_h*DepositShare_h(i);
  }

  //Government pays for energy
  EnergyPayments += Exp_en_g;

  //Energy sector receives revenue
  Deposits_e(1)+=EnergyPayments;

  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)+=EnergyPayments*DepositShare_e(i);
    Deposits_eb(1,i)+=EnergyPayments*DepositShare_e(i);
    Inflows(i)+=EnergyPayments*DepositShare_e(i);
  }

  //Profit of energy sector
  ProfitEnergy(1)=EnergyPayments+InterestDeposits_e-Wages_en+(CapitalStock_e(1)-CapitalStock_e(2))-t_CO2_en*Emiss_en-FuelCost;

  //Energy shock
  if((flag_energyshocks==1 && t>t_regime_shifts) || (flag_energyshocks==2 && t>t_regime_shifts))
  {
    Transfer_shock=tre*(EnergyPayments-(EnergyPayments/c_en(1))*(c_en_preshock));  //Update transfer based on actual revenues of energy sector for following timestep
    for (const string& cl:classes_mh){
      Transfer_shock_mh[cl]=Transfer_shock*Transfer_shock_sh_mh[cl];
    }

    Taxes_e_shock=aliqee*max(0.0, EnergyPayments-(EnergyPayments/c_en(1))*(c_en_preshock)); //Tax the additional revenues due to increased price (difference with revenues if price would have been the same)
    if(Deposits_e(1)>=Taxes_e_shock)
    {
      Deposits_e(1)-=Taxes_e_shock;
      for(i=1; i<=NB; i++)
      {
        Deposits_b(1,i)-=Taxes_e_shock*DepositShare_e(i);
        Outflows(i)+=Taxes_e_shock*DepositShare_e(i);
        Deposits_eb(1,i)-=Taxes_e_shock*DepositShare_e(i);
      }
    }
    else
    {
      std::cerr << "\n\n Energy sector cannot pay excess profit tax in period " << t << endl;
      Errors << "\n Energy sector cannot pay excess profit tax in period " << t << endl;
      exit(EXIT_FAILURE);
    }
  }

  //Energy sector pays fuel
  if(Deposits_e(1)>=FuelCost)
  {
    Deposits_e(1)-=FuelCost;
    for(i=1; i<=NB; i++)
    {
      Deposits_b(1,i)-=FuelCost*DepositShare_e(i);
      Outflows(i)+=FuelCost*DepositShare_e(i);
      Deposits_eb(1,i)-=FuelCost*DepositShare_e(i);
    }
  }
  else
  {
    std::cerr << "\n\n Energy sector cannot pay for fuel in period " << t << endl;
    Errors << "\n Energy sector cannot pay for fuel in period " << t << endl;
    exit(EXIT_FAILURE);
  }

  //Energy sector pays C02 tax
  if(Deposits_e(1)>=t_CO2_en*Emiss_en)
  {
    Taxes_CO2_e=t_CO2_en*Emiss_en;
    Taxes_CO2(1)+=Taxes_CO2_e;
  }
  else if(Deposits_e(1)>=0)
  {
    Taxes_CO2_e=Deposits_e(1);
    Taxes_CO2(1)+=Taxes_CO2_e;
  }
  else
  {
    std::cerr << "\n\n Energy sector has negative deposits in period " << t << endl;
    Errors << "\n Energy sector has negative deposits in period " << t << endl;
    exit(EXIT_FAILURE);
  }
  Deposits_e(1)-=Taxes_CO2_e;
  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)-=Taxes_CO2_e*DepositShare_e(i);
    Outflows(i)+=Taxes_CO2_e*DepositShare_e(i);
    Deposits_eb(1,i)-=Taxes_CO2_e*DepositShare_e(i);
  }

  //Energy sector pays taxes on profits
  if(Deposits_e(1) >= aliq * ProfitEnergy(1))
  {
    Taxes_e = aliq * ProfitEnergy(1);
    Taxes_g += Taxes_e;
  }
  else if(Deposits_e(1)>=0)
  {
    Taxes_e = Deposits_e(1);
    Taxes_g += Taxes_e;
  }
  else
  {
    std::cerr << "\n\n Energy sector has negative deposits in period " << t << endl;
    Errors << "\n Energy sector has negative deposits in period " << t << endl;
    exit(EXIT_FAILURE);
  }
  Deposits_e(1)-=Taxes_e;
  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)-=Taxes_e*DepositShare_e(i);
    Outflows(i)+=Taxes_e*DepositShare_e(i);
    Deposits_eb(1,i)-=Taxes_e*DepositShare_e(i);
  }

  //Fossil fuel agent pays taxes following an eventual energy price shock
  if((flag_energyshocks==3 && t>t_regime_shifts) || (flag_energyshocks==4 && t>t_regime_shifts))
  {
    //Energy sector pays a tax on excess profits following a shock on fossil fuel price
    Taxes_e_ff_shock=aliqee*max(0.0,(max(0.0,EnergyPayments-(EnergyPayments/c_en(1))*(c_en_preshock)))-(max(0.0,FuelCost-(FuelCost/pf)*pf_preshock)));
    if(Deposits_e(1)>=Taxes_e_ff_shock)
    {
      Deposits_e(1)-=Taxes_e_ff_shock;
      for(i=1; i<=NB; i++)
      {
        Deposits_b(1,i)-=Taxes_e_ff_shock*DepositShare_e(i);
        Outflows(i)+=Taxes_e_ff_shock*DepositShare_e(i);
        Deposits_eb(1,i)-=Taxes_e_ff_shock*DepositShare_e(i);
      }
    }
    else
    {
      std::cerr << "\n\n Energy sector cannot pay excess profit tax in period " << t << endl;
      Errors << "\n Energy sector cannot pay excess profit tax in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    //Fossil fuel agent pays a tax on excess profits following a shock on fossil fuel price
    Taxes_f_ff_shock=aliqef*max(0.0,FuelCost-(FuelCost/pf)*pf_preshock);
  }

  //Energy sector dividends and bonuses
  double NetProfEn = ProfitEnergy(1)-Taxes_e_shock-Taxes_e_ff_shock;

  if(NetProfEn>0)
  {
    double PEn_diff = ProfitEnergy(1) - ProfitEnergy(2);
    if(PEn_diff > 0 & ProfitEnergy(2) > 0) //If profit increased and in previous period were positive, pay bonus as fraction of profits increase
    {
      Bonuses_e=bonuses_share*PEn_diff;
    }
    else if(PEn_diff > 0 & ProfitEnergy(1) > 0)  //If profit increased and in previous period were negative but in current positive, pay bonus as fraction of new positive profits
    {
      Bonuses_e=bonuses_share*ProfitEnergy(1);
    }
    else //If profit decreased, or are negative, no bonuses
    {
      Bonuses_e=0;
    }

    if(Deposits_e(1)>=(Bonuses_e+de*(NetProfEn-Bonuses_e)))
    {
      Dividends_e=de*(NetProfEn-Bonuses_e);
    }
    else if(Deposits_e(1)>=Bonuses_e) //If can afford to pay more than bonuses, pay bonuses completely and dividends as residuals
    {
      Dividends_e=Deposits_e(1)-Bonuses_e;
    }
    else if(Deposits_e(1)>=0) //If cannot afford to pay bonuses entirely, pay only bonuses as much as possible and no dividends
    {
      Bonuses_e=Deposits_e(1); 
      Dividends_e=0; 
    }
    else
    {
      std::cerr << "\n\n Energy sector has negative deposits in period " << t << endl;
      Errors << "\n Energy sector has negative deposits in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    Dividends_h(1)+=Dividends_e;
    Bonuses_h+=Bonuses_e;
    Deposits_e(1)-=Dividends_e+Bonuses_e;
    for(i=1; i<=NB; i++)
    {
      Deposits_b(1,i)-=(Dividends_e+Bonuses_e)*DepositShare_e(i);
      Outflows(i)+=(Dividends_e+Bonuses_e)*DepositShare_e(i);
      Deposits_eb(1,i)-=(Dividends_e+Bonuses_e)*DepositShare_e(i);
    }
  }
  else
  {
    Dividends_e=0;
    Bonuses_e=0;
  }

  Bonuses_e_mh["ma"]=Bonuses_e;
  Bonuses_mh["ma"]+=Bonuses_e_mh["ma"];
  for (const string& cl:classes_mh){ //Split among household classes
    Dividends_e_mh[cl]=Dividends_e*Ownership_sh_e_mh[cl];
    Dividends_mh[cl](1)+=Dividends_e_mh[cl];
  }

  //Fossil fuel agent receives fuel payment and transfers a portion to households
  Deposits_fuel(1)+=(FuelCost-Taxes_f_ff_shock);
  Deposits_fuel_cb(1)+=(FuelCost-Taxes_f_ff_shock);
  TransferFuel=d_f*Deposits_fuel(1);
  for (const string& cl:classes_mh){
      TransferFuel_mh[cl]=TransferFuel*Ownership_sh_ff_mh[cl];
    }
  Deposits_fuel(1)-=TransferFuel;
  Deposits_fuel_cb(1)-=TransferFuel;

  //Households receive energy sector, C-firm dividend and small transfer from fossil fuel agent
  double Taxes_div_no1_h = 0; //All taxes on dividends apart from K-firms
  double Taxes_bon_no1_h = 0; //All taxes on bonuses apart from K-firms
  Taxes_wealth_h = 0;
  for (const string& cl:classes_mh){ //Split among household classes
    //Add taxes on C-firms and energy sector dividends and bonuses and fossil fuel sector transfers to taxes on dividends and bonuses
    double taxes_div_cl=aliqdiv*(Dividends_2_mh[cl]+Dividends_e_mh[cl]+TransferFuel_mh[cl]);
    double taxes_bon_cl=aliqw_mh[cl]*(Bonuses_2_mh[cl]+Bonuses_e_mh[cl]);
    Taxes_div_mh[cl](1)=Taxes_div_1_mh[cl] + taxes_div_cl; //Each household class
    Taxes_bon_mh[cl](1)=Taxes_bon_1_mh[cl] + taxes_bon_cl; //Each household class

    Taxes_mh[cl] = Taxes_w_mh[cl] + Taxes_div_mh[cl](1) + Taxes_bon_mh[cl](1); //Total taxes
    Deposits_mh[cl](1)+=(Dividends_2_mh[cl]+Bonuses_2_mh[cl]+Dividends_e_mh[cl]+Bonuses_e_mh[cl]+TransferFuel_mh[cl]-taxes_div_cl-taxes_bon_cl);

    Taxes_div_no1_h += taxes_div_cl;
    Taxes_bon_no1_h += taxes_bon_cl;

    //Huoseholds pay taxes on deposits 
    Taxes_wealth_mh[cl] = Deposits_mh[cl](1) * aliqwealth_mh[cl];
    Deposits_mh[cl](1) -= Taxes_wealth_mh[cl];
    Taxes_mh[cl] += Taxes_wealth_mh[cl];
    Taxes_wealth_h += Taxes_wealth_mh[cl];
  }

  Taxes_div_h = Taxes_div_1_h + Taxes_div_no1_h;        //Overall
  Taxes_bon_h = Taxes_bon_1_h + Taxes_bon_no1_h;        //Overall

  double delta_Deposits_h_2 = Dividends_2+Bonuses_2+Dividends_e+Bonuses_e+TransferFuel-Taxes_div_no1_h-Taxes_wealth_h-Taxes_bon_no1_h; //Overall change in deposits of households
  Deposits_h(1)+=delta_Deposits_h_2;
  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)+=delta_Deposits_h_2*DepositShare_h(i);
    Deposits_hb(1,i)+=delta_Deposits_h_2*DepositShare_h(i);
    Inflows(i)+=delta_Deposits_h_2*DepositShare_h(i);
  }

  //Total households taxes to government
  Taxes_h = Taxes_w_h + Taxes_div_h + Taxes_wealth_h + Taxes_bon_h; //Overall
  Taxes_g += Taxes_h;


  //C-firms which exit due to negative equity, inability to make payments or low market share are prepared for exit
  for(j=1; j<=N2; j++)
  {
    NW_2(1,j)=CapitalStock(1,j)+deltaCapitalStock(1,j)+Inventories(1,j)+Deposits_2(1,j)-Loans_2(1,j);
    if(NW_2(1,j)<0 && exit_payments2(j)==0 && exiting_2(j)==0 && exit_marketshare2(j)==0)
    {
      exit_equity2(j)=1;
      exiting_2(j)=1;
    }

    if(exit_payments2(j)==1 || exit_equity2(j)==1)
    {
      sendingBank=BankingSupplier_2(j);
      baddebt_2(j)=Loans_2(1,j)-Deposits_2(1,j);
      if(Loans_2(1,j)>Deposits_2(1,j))
      {
        baddebt_b(sendingBank)+=(Loans_2(1,j)-Deposits_2(1,j));
        Loans_b(1,sendingBank)=max(0.0,Loans_b(1,sendingBank)-Loans_2(1,j));
        Deposits_b(1,sendingBank)-=Deposits_2(1,j);
      }
      else
      {
        Deposits_recovered_2+=Deposits_2(1,j)-Loans_2(1,j);
        for (const string& cl:classes_mh){
          Deposits_recovered_2_mh[cl] += Deposits_2(1,j) * Ownership_sh_2_i_mh[cl](j);
        }
        Outflows(sendingBank)+=Deposits_2(1,j)-Loans_2(1,j);
        Loans_b(1,sendingBank)=max(0.0,Loans_b(1,sendingBank)-Loans_2(1,j));
        Deposits_b(1,sendingBank)-=Deposits_2(1,j);
      }
      Deposits_2(1,j)=0;
      Loans_2(1,j)=0;
    }

    if(exit_marketshare2(j)==1 && exit_payments2(j)==0 && exit_equity2(j)==0)
    {
      exiting_2(j)=1;
      sendingBank=BankingSupplier_2(j);
      baddebt_2(j)=Loans_2(1,j)-Deposits_2(1,j);
      if(Loans_2(1,j)>Deposits_2(1,j))
      {
        baddebt_b(sendingBank)+=(Loans_2(1,j)-Deposits_2(1,j));
        Loans_b(1,sendingBank)=max(0.0,Loans_b(1,sendingBank)-Loans_2(1,j));
        Deposits_b(1,sendingBank)-=Deposits_2(1,j);
      }
      else
      {
        Deposits_recovered_2+=Deposits_2(1,j)-Loans_2(1,j);
        for (const string& cl:classes_mh){
          Deposits_recovered_2_mh[cl] += Deposits_2(1,j) * Ownership_sh_2_i_mh[cl](j);
        }
        Outflows(sendingBank)+=Deposits_2(1,j)-Loans_2(1,j);
        Loans_b(1,sendingBank)=max(0.0,Loans_b(1,sendingBank)-Loans_2(1,j));
        Deposits_b(1,sendingBank)-=Deposits_2(1,j);
      }
      Deposits_2(1,j)=0;
      Loans_2(1,j)=0;
    }

    if(exiting_2(j)==1 && exit_payments2(j)==0 && exit_equity2(j)==0 && exit_marketshare2(j)==0){
      exit_payments2(j)=1;
    }
	}

  //Households receive deposits recovered
  for (const string& cl:classes_mh){
    Deposits_mh[cl]+=Deposits_recovered_2_mh[cl];
  }
  Deposits_h(1)+=Deposits_recovered_2;

  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)+=Deposits_recovered_2*DepositShare_h(i);
    Deposits_hb(1,i)+=Deposits_recovered_2*DepositShare_h(i);
    Inflows(i)+=Deposits_recovered_2*DepositShare_h(i);
  }
}

void ALLOC(void)
{
  n=1;
	ftot=0;
  Utilisation=Q2.Sum()/K.Sum();
	Cres=Cons;
  cpi_temp=cpi(1);

  //Output shocks
  if(flag_outputshocks==1 || flag_outputshocks==4)
  {
    for (j=1; j<=N2; j++)
	  {
      C_loss(j)=shocks_output2(j)*Q2(j);
      Q2(j)=(1-shocks_output2(j))*Q2(j);
    }
  }

  if(flag_outputshocks==2 || flag_outputshocks==5)
  {
    loss=Q2.Sum()*shocks_output2(1);
    while(loss>0)
    {
      ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      if(Q2(ranj)>=loss)
      {
        Q2(ranj)-=loss;
        C_loss(ranj)=loss;
        loss=0;
      }
      else
      {
        loss-=Q2(ranj);
        C_loss(ranj)=Q2(ranj);
        Q2(ranj)=0;
      }
    }
  }

  if(flag_outputshocks==3 || flag_outputshocks==6)
  {
    loss=Q2.Sum()*shocks_output2(1);
    while(loss>0)
    {
      ranj=0;
      while(ranj==0)
      {
        ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      }
      rnd=ran1(p_seed);
      lossj=rnd*Q2(ranj);
      if(lossj>=loss)
      {
        Q2(ranj)-=loss;
        C_loss(ranj)=loss;
        loss=0;
      }
      else
      {
        loss-=lossj;
        Q2(ranj)-=lossj;
        C_loss(ranj)=lossj;
      }
    }
  }

  //Investment shocks
  if(flag_inventshocks==1)
  {
    for (j=1; j<=N2; j++)
	  {
      if(N(2,j)>0)
      {
        ptemp=Inventories(2,j)/N(2,j);
        Loss_Inventories(j)+=shocks_invent(j)*N(2,j)*ptemp;
        N(2,j)=(1-shocks_invent(j))*N(2,j);
        Inventories(2,j)-=Loss_Inventories(j);
      }
    }
  }

  if(flag_inventshocks==2)
  {
    loss=N.Row(2).Sum()*shocks_invent(1);
    while(loss>0)
    {
      ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      if(N(2,ranj)>0)
      {
        ptemp=Inventories(2,ranj)/N(2,ranj);
        if(N(2,ranj)>=loss)
        {
          ptemp=Inventories(2,ranj)/N(2,ranj);
          Loss_Inventories(ranj)+=loss*ptemp;
          Inventories(2,ranj)-=Loss_Inventories(ranj);
          N(2,ranj)-=loss;
          loss=0;
        }
        else
        {
          Loss_Inventories(ranj)+=Inventories(2,ranj);
          Inventories(2,ranj)=0;
          loss-=N(2,ranj);
          N(2,ranj)=0;
        }
      }
    }
  }

  if(flag_inventshocks==3)
  {
    loss=N.Row(2).Sum()*shocks_invent(1);
    while(loss>0)
    {
      ranj=0;
      while(ranj==0)
      {
        ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      }
      rnd=ran1(p_seed);
      lossj=rnd*N(2,ranj);
      if(lossj>=loss && N(2,ranj)>0)
      {
        ptemp=Inventories(2,ranj)/N(2,ranj);
        Loss_Inventories(ranj)+=loss*ptemp;
        Inventories(2,ranj)-=Loss_Inventories(ranj);
        N(2,ranj)-=loss;
        loss=0;
      }
      else
      {
        ptemp=Inventories(2,ranj)/N(2,ranj);
        Loss_Inventories(ranj)+=lossj*ptemp;
        Inventories(2,ranj)-=lossj*ptemp;
        loss-=lossj;
        N(2,ranj)-=lossj;
      }
    }
  }

	for (j=1; j<=N2; j++)
	{
		f_temp2(j)=f2(1,j);
		ftot(1)+=f_temp2(j);
		Q2temp(j)=Q2(j)+N(2,j);
	}

  //Consumption demand is distributed among C-firms based on market shares
	while (Cres >= 1 && ftot(1) > 0)
	{
    Cresbis=Cres;
		for (j=1; j<=N2; j++)
		{
			if (f_temp2(j) > 0)
			{
				D_temp2(j)=Cres/cpi_temp*f_temp2(j);

				if (n==1)
				{
					D2(1,j)+=D_temp2(j);
				}

				if (D_temp2(j) <= Q2temp(j)) //If can satisfy demand
				{
					if (n > 1)
          {
						D2(1,j)+=D_temp2(j);
          }
          S2(1,j)+=p2(j)*D_temp2(j); //Update revenues
          Cresbis-=D_temp2(j)*p2(j); //Update residual consumption
          if(n==1)
          {
						l2(j)=1;                 //Unsatisfied consumption demand
          }
          Q2temp(j)-=D_temp2(j);		 //Residual capacity
				}
				else						           //If cannot satisfy demand
				{
					if (n > 1)
          {
						D2(1,j)+=Q2temp(j);
          }
          S2(1,j)+=p2(j)*Q2temp(j);
          Cresbis-=Q2temp(j)*p2(j);
          f_temp2(j)=0;
          if(n==1)
          {
            l2(j)=1+(D_temp2(j)-Q2temp(j));
          }
          Q2temp(j)=0;
				}
			}
		}
		ftot(1)=f_temp2.Sum();
	  f_temp2/=ftot(1);
		Cres=Cresbis;
    cpi_temp=0;
    for (j=1; j<=N2; j++)
    {
      cpi_temp+=p2(j)*f_temp2(j);
    }
		n++;
	}

  //Total nominal consumption is calculated (based on firms revenues)
  Consumption=S2.Row(1).Sum();
    //This is done to ensure that household deposits do not become negative due to consumption (may happen due to rounding issues when liquidity constraint is binding)
  while(Consumption>Cons)
  {
    for (j=1; j<=N2; j++)
	  {
      if(S2(1,j)>(S2(1,j)/S2.Row(1).Sum()*(Consumption-Cons)))
      {
        S2(1,j)-=(S2(1,j)/S2.Row(1).Sum()*(Consumption-Cons));
      }
    }
    Consumption=S2.Row(1).Sum();
  }

  //Calculate government final consumption (reduced proportionally) and update total expenditure
  Consumption_g = Cons_g * Consumption/Cons;
  Exp_tot_g = Exp_tot_g - Cons_g + Consumption_g;

  //Households pay for consumption
  Consumption_h = Consumption - Consumption_g;
  double cons_pre_adjustment = Consumption_h;
  for (const string& cl:classes_mh)
  {
    //Calculate each class final consumption based on shares previously calculated
    Consumption_mh[cl] = cons_pre_adjustment * Cons_sh_mh[cl];
    //Check to avoid deposits of each class being negative due to rounding issues (ensuring < tolerance, no issues for SFC with firms revenues)
    if(Deposits_mh[cl](1)<Consumption_mh[cl] && fabs(Deposits_mh[cl](1)-Consumption_mh[cl])<tolerance)
    {
      Consumption_h-=(Consumption_mh[cl]-Deposits_mh[cl](1));
      Consumption_mh[cl]=Deposits_mh[cl](1);
    }
    Deposits_mh[cl](1)-=Consumption_mh[cl];
    //Update total expenditure
    Expenditure_tot_mh[cl] = Expenditure_tot_mh[cl] - Cons_mh[cl] + Consumption_mh[cl];
  }
  Expenditure_tot_h = Expenditure_tot_h - Cons_h + Consumption_h;
  Deposits_h(1)-=Consumption_h;

  for(i=1; i<=NB; i++)
  {
    Deposits_b(1,i)-=Consumption_h*DepositShare_h(i);
    Outflows(i)+=Consumption_h*DepositShare_h(i);
    Deposits_hb(1,i)-=Consumption_h*DepositShare_h(i);
  }

  //Real consumption is calculated
  for(j=1; j<=N2; j++)
  {
    Consumption_r+=S2(1,j)/p2(j);
  }

  //C-firms receive revenue & we re-compute the CPI & we update inventories
  cpi(1)=S2.Row(1).Sum()/Consumption_r;
  for (j=1; j<=N2; j++)
	{
    receivingBank=BankingSupplier_2(j);
    Deposits_2(1,j)+=S2(1,j);
    Deposits_b(1,receivingBank)+=S2(1,j);
    Inflows(receivingBank)+=S2(1,j);
    N(1,j)=N(2,j)+flag_inventories*(Q2(j)-S2(1,j)/p2(j));
    if(N(1,j)<0)
    {
      if(fabs(N(1,j))/Q2.Sum()<tolerance)
      {
        N(1,j)=0;
      }
      else
      {
        std::cerr << "\n\n ERROR: Inventories of C-firm " << j << " are negative in period " << t << endl;
        Errors << "\n Inventories of C-firm " << j << " are negative in period " << t << endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  // Compute households' energy price
  c_en_h = c_en(1) * ratio_c_en_h_firms;
  S1_temp.Row(1)=S1;
  S2_temp.Row(1)=S2.Row(1);
}

void ENTRYEXIT(void)
{
  //Save sales as S1 & S2 are reset during entry

  Sales1=S1;
  Sales2=S2.Row(1);

  for (i=1; i<=N1; i++)
	{
    if(Deposits_1(1,i)<0)
    {
      std::cerr << "\n\n ERROR: K-firm " << i << " has negative deposits in period " << t << endl;
      Errors << "\n  K-firm " << i << " has negative deposits in period " << t << endl;
      exit(EXIT_FAILURE);
    }
	}

  if (ns1 > 0)
  {
		mD1/=ns1;
  }
	else
  {
    mD1=Deposits_1.Row(2).Sum()/N1r;
  }

  //Exiting K-firms lose customers; K-firm to be copied is chosen
  for (i=1; i<=N1; i++)
  {
    if (exiting_1(i)==1)
    {
      flag=0;
      for (j=1; j<=N2; j++)
      {
        Match(j,i)=0;
        if(fornit(j)==i)
        {
          fornit(j)=0;
        }
      }
      if(exiting_1.Sum()<N1r)
      {
        while (flag == 0)
        {
          rni=int(ran1(p_seed)*N1*N2)%N1+1;
          if (exiting_1(rni)==0)
          {
            ee1(i)=rni;
            flag=1;
          }
        }
      }
      else
      {
        ee1(i)=i;
      }
    }
  }

  for (i=1; i<=N1; i++)
  {
    if (exiting_1(i)==1)
    {
      //For entering K-firms, most variables are copied from a surviving K-firm as in the original DSK
      iii=int(ee1(i));
      f1(1,i)=0;
      f1(2,i)=0;
      for (const string& cl:classes_mh){
        Ld1_rd_i_mh[cl](i)=0;
      }
      A1(i)=A1(iii);
      shocks_labprod1(i)=shocks_labprod1(iii);
      shocks_eneff1(i)=shocks_eneff1(iii);
      p1(i)=p1(iii);
      A1(i)=A1(iii);
      A1p(i)=A1p(iii);
      A1_ef(i)=A1_ef(iii);
      A1p_ef(i)=A1p_ef(iii);
      A1_en(i)=A1_en(iii);
      A1p_en(i)=A1p_en(iii);

      EnvSubsidiesPerMachine_1_i(i)=EnvSubsidiesPerMachine_1_i(iii);

      //Entering K-firms' deposits, however, are received as a transfer from households
      receivingBank=BankingSupplier_1(i);
      multip_entry=ran1(p_seed);
      multip_entry=w1inf+multip_entry*(w1sup-w1inf);
      if(Deposits_h(1)>=multip_entry*mD1)
      {
        //If households have sufficient deposits, transfer is equal to mean deposits of surviving firms
        injection=multip_entry*mD1;
        Deposits_h(1)-=injection;
        for (const string& cl:classes_mh){ //Injection from each household class
          injection_mh[cl]=injection*Entry_financing_sh_mh[cl];
          Deposits_mh[cl](1)-=injection_mh[cl];
          if (Deposits_mh[cl](1)<0){
            std::cerr << "\n\n Household class " << cl << " cannot finance K-firm entry in period " << t << "CHANGE CODE" << endl;
            Errors << "\n Household class " << cl << " cannot finance K-firm entry in period " << t << "CHANGE CODE" << endl;
            //exit(EXIT_FAILURE);
          }
        }
        FirmTransfers+=injection;
        FirmTransfers_1+=injection;
        Injection_1(i)=injection;
        for(j=1; j<=NB; j++)
        {
          Deposits_hb(1,j)-=DepositShare_h(j)*injection;
          Deposits_b(1,j)-=DepositShare_h(j)*injection;
          Outflows(j)+=DepositShare_h(j)*injection;
        }
        Deposits_1(1,i)+=injection;
        Deposits_b(1,receivingBank)+=injection;
        Inflows(receivingBank)+=injection;
      }
      else
      {
        //If households cannot finance K-firm entry, this is either financed by gov. (flag_entry=0) or banks simply create the deposits (without corresponding loan) and book this as a loss (flag_entry=1)
        std::cerr << "\n\n Households cannot finance K-firm entry in period " << t << endl;
        Errors << "\n Households cannot finance K-firm entry in period " << t << endl;
        injection=0;  //In case Deposits_h not positive, this has to be 0 (all financed by either government or banks). But it should also give an error.
        if(Deposits_h(1)>0)  //Households finance as much as they can
        {
          injection=Deposits_h(1);
          Deposits_h(1)-=injection;
          for (const string& cl:classes_mh){ //Injection from each household class
            injection_mh[cl]=injection*Entry_financing_sh_mh[cl];
            Deposits_mh[cl](1)-=injection_mh[cl];
            if (Deposits_mh[cl](1)<0){
              std::cerr << "\n\n Household class " << cl << " cannot finance K-firm entry in period " << t << "CHANGE CODE" << endl;
              Errors << "\n Household class " << cl << " cannot finance K-firm entry in period " << t << "CHANGE CODE" << endl;
              //exit(EXIT_FAILURE);
            }
          }
          FirmTransfers+=injection;
          FirmTransfers_1+=injection;
          Injection_1(i)=injection;
          for(j=1; j<=NB; j++)
          {
            Deposits_b(1,j)-=Deposits_hb(1,j);
            Outflows(j)+=Deposits_hb(1,j);
            Deposits_hb(1,j)=0;
            DepositShare_h(j)=(NL_1(j)+NL_2(j))/(N1+N2);
          }
          Deposits_1(1,i)+=injection;
          Deposits_b(1,receivingBank)+=injection;
          Inflows(receivingBank)+=injection;
        }

        injection2=multip_entry*mD1-injection; //Injection left to be financed
        if(flag_entry==1) //Financed by banks
        {
          LossEntry_b(receivingBank)+=injection2;
        }
        else              //Financed by government
        {
          EntryCosts+=injection2;
          FirmTransfers_1+=injection2;
          Inflows(receivingBank)+=injection2;
        }
        Injection_1(i)+=injection2;
        Deposits_1(1,i)+=injection2;
        Deposits_b(1,receivingBank)+=injection2;
      }

      S1(i)=p1(i)*step;
      stepbis=step;
      while (stepbis > 0)
      {
        rni=int(ran1(p_seed)*N1*N2)%N2+1;
        if (Match(rni,i) == 0)
        {
          Match(rni,i)=1;
          stepbis--;
        }
      }

      //Update ownership based on contribution by each household class
      if (injection!=0){ //If no injection no need for updating (division by 0 would give an error)
        for (const string& cl:classes_mh){
          Ownership_sh_1_i_mh[cl](i)=injection_mh[cl]/injection; //Based on injection only of households (as if neither gov or banks contributed)
        }
      }

    }
  }


  // C-Firms
  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1)
    {
      f2_exit+=f2(1,j);
    }
  }

  if(exiting_2.Sum()>0){
    CurrentDemand=D2.Row(1).Sum();
    if(f2_exit>0){
      n_mach_needed=max(exiting_2.Sum(),ceil(f2_exit*CurrentDemand/dim_mach/u));
    }else{
      f2_exit=f2_entry_min*exiting_2.Sum();
      n_mach_needed=max(exiting_2.Sum(),ceil(f2_exit*CurrentDemand/dim_mach/u));
    }
  }

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1)
    {
      for (i=1; i<=N1; i++)
      {
        for (tt=t0; tt<=t; tt++)
        {
          if (gtemp[tt-1][i-1][j-1]>0)
          {
            n_mach_exit+=gtemp[tt-1][i-1][j-1];
          }
        }
      }
    }
  }

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1)
    {
      for (i=1; i<=N1; i++)
      {
        for (tt=t; tt>=t0; tt--)
        {
          if (gtemp[tt-1][i-1][j-1]>0)
          {
            C_secondhand(tt,i)=C(tt,i);
          }
        }
      }
    }
    else
    {
      //Calculate mean deposits of surviving C-firms
      mD2+=Deposits_2(1,j);
			ns2++;
    }
  }

  //The capital stock of exiting C-firms is transferred to their respective banks up to the value of bad debt; the value of second-hand capital is marked down depending on age of the machines
  n_mach_exit2=min(n_mach_needed,n_mach_exit);
  while(n_mach_exit2>0)
  {
    for(j=1; j<=N2; j++)
    {
      if(exiting_2(j)==1)
      {
        receivingBank=BankingSupplier_2(j);
        baddebt_2_temp=baddebt_2(j);
        for (i=1; i<=N1; i++)
        {
          for (tt=t; tt>=t0; tt--)
          {
            if (gtemp[tt-1][i-1][j-1]>0 && n_mach_exit2>0 && C(tt,i)<=C_secondhand.Minimum())
            {
              markdownCapital=max(0.0,(1-(double)age[tt-1][i-1][j-1]/(agemax)));
              g_secondhand[tt-1][i-1]+=min(n_mach_exit2,gtemp[tt-1][i-1][j-1]);
              age_secondhand[tt-1][i-1]=age[tt-1][i-1][j-1];
              g_secondhand_p[tt-1][i-1]=markdownCapital*g_price[tt-1][i-1][j-1];
              if(baddebt_2_temp>0)
              {
                capitalRecovered(receivingBank)+=min(n_mach_exit2,gtemp[tt-1][i-1][j-1])*g_secondhand_p[tt-1][i-1];
                baddebt_2_temp-=min(n_mach_exit2,gtemp[tt-1][i-1][j-1])*g_secondhand_p[tt-1][i-1];
              }
              n_mach_exit2-=min(n_mach_exit2,gtemp[tt-1][i-1][j-1]);
              C_secondhand(tt,i)=1000000;
            }
          }
        }
      }
    }
  }

	if (ns2 > 0)
	{
		mD2/=ns2;
	}
	else
  {
    std::cerr << "\n\n ERROR: All C-firms are exiting in period " << t << endl;
    Errors << "\n All C-firms are exiting in period " << t << endl;
    exit(EXIT_FAILURE);
  }

  //K-firms lose exiting C-firms as customers
  for (j=1; j<=N2; j++)
  {
    if (exiting_2(j)==1)
    {
      flag=0;
      indforn=int(fornit(j));
      if (indforn>=1)
      {
        Match(j,fornit(j))=0;
        fornit(j)=0;
      }
    }
  }

  //Determine the number of machines which each entering C-firm will have based on number of available second-hand machines
  n_exit2=exiting_2.Sum();
  n_mach_resid=min(n_mach_needed,n_mach_exit);
  if(n_exit2>n_mach_resid)
  {
    std::cerr << "\n\n ERROR: Not enough second hand capital in period " << t << endl;
    Errors << "\n Not enough second hand capital in period " << t << endl;
    exit(EXIT_FAILURE);
  }

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1)
    {
      n_mach_entry(j)=1;
      n_mach_resid--;
      k_entry(j)=ran1(p_seed);
    }
  }

  n_mach_resid2=n_mach_resid;

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1 && n_mach_resid>0)
    {
      n_mach_entry(j)+=floor(k_entry(j)/k_entry.Sum()*n_mach_resid2);
      n_mach_resid-=floor(k_entry(j)/k_entry.Sum()*n_mach_resid2);
    }
  }


  if(n_mach_resid<0 && n_mach_exit<n_mach_needed)
  {
    std::cerr << "\n\n ERROR: Remaining second hand machines are negative in period " << t << endl;
    Errors << "\n Remaining second hand machines are negative " << t << endl;
    exit(EXIT_FAILURE);
  }

  while(n_mach_resid>0)
  {
    rni=int(ran1(p_seed)*N1*N2)%N2+1;
    if(exiting_2(rni)==1 && n_mach_resid>0)
    {
      n_mach_entry(rni)++;
      n_mach_resid--;
    }
  }

  //Second-hand capital is purchased by households; below will be transferred to newly entering firms
  capitalRecoveredTot=capitalRecovered.Sum();
  if(Deposits_h(1)>=(capitalRecoveredTot)) //Households can finance entirely
  {
    Deposits_h(1)-=capitalRecoveredTot;
    for (const string& cl:classes_mh){
      capitalRecoveredTot_mh[cl]=capitalRecoveredTot*Entry_financing_sh_mh[cl];
      Deposits_mh[cl]-=capitalRecoveredTot_mh[cl];
      if (Deposits_mh[cl](1)<0){
        std::cerr << "\n\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
        Errors << "\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
        //exit(EXIT_FAILURE);
      }
    }
    FirmTransfers+=capitalRecoveredTot;
    FirmTransfers_2+=capitalRecoveredTot;
    for(i=1; i<=NB; i++)
    {
      Deposits_hb(1,i)-=DepositShare_h(i)*capitalRecoveredTot;
      Deposits_b(1,i)-=DepositShare_h(i)*capitalRecoveredTot;
      if(capitalRecovered(i)>=DepositShare_h(i)*capitalRecoveredTot)
      {
        Inflows(i)+=(capitalRecovered(i)-DepositShare_h(i)*capitalRecoveredTot);
      }
      else
      {
        Outflows(i)+=(DepositShare_h(i)*capitalRecoveredTot-capitalRecovered(i));
      }
    }
  }
  else if (Deposits_h(1)>=0)   //Households finance with all the remaining deposits
  {
    //Household deposits are insufficient to purchase second hand capital; households buy as much as they can, the rest is financed either by government (flag_entry=0) or booked as a loss by banks (flag_entry=1)
    std::cerr << "\n\n CODE NOT UNDERSTOOD Households cannot purchase second-hand capital in period " << t << endl;
    Errors << "\n CODE NOT UNDERSTOOD Households cannot purchase second-hand capital in period " << t << endl;
    FirmTransfers+=Deposits_h(1);
    FirmTransfers_2+=Deposits_h(1);
    capitalRecovered2=capitalRecovered;
    for(i=1; i<=NB; i++)
    {
      if(capitalRecovered2(i)>=Deposits_hb(1,i))
      {
        capitalRecovered2(i)-=Deposits_hb(1,i);
        Deposits_b(1,i)-=Deposits_hb(1,i);
        Deposits_h(1)-=Deposits_hb(1,i);
        for (const string& cl:classes_mh){
          Deposits_mh[cl]-=Deposits_hb(1,i)*Entry_financing_sh_mh[cl];
          if (Deposits_mh[cl](1)<0){
            std::cerr << "\n\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
            Errors << "\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
            //exit(EXIT_FAILURE);
          }
        }
        Deposits_hb(1,i)=0;
      }
      else
      {
        Deposits_b(1,i)-=Deposits_hb(1,i);
        Deposits_hb(1,i)-=capitalRecovered2(i);
        Deposits_h(1)-=capitalRecovered2(i);
        for (const string& cl:classes_mh){
          Deposits_mh[cl]-=capitalRecovered2(i)*Entry_financing_sh_mh[cl];
          if (Deposits_mh[cl](1)<0){
            std::cerr << "\n\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
            Errors << "\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
            //exit(EXIT_FAILURE);
          }
        }
        capitalRecovered2(i)=0;
        Outflows(i)+=Deposits_hb(1,i);
        Deposits_hb(1,i)=0;
      }
      DepositShare_h(i)=(NL_1(i)+NL_2(i))/(N1+N2);
    }

    for(i=1; i<=NB; i++) //Calculate how much left to be recovered from each bank
    {
      capitalRecoveredShare(i)=capitalRecovered2(i)/capitalRecovered2.Sum();
      capitalRecovered2(i)-=(capitalRecoveredShare(i)*Deposits_h(1));
      Inflows(i)+=(capitalRecoveredShare(i)*Deposits_h(1));
    }

    if(flag_entry==1) //Financed by banks
    {
      for(i=1; i<=NB; i++)
      {
        LossEntry_b(i)+=capitalRecovered2(i);
      }
    }
    else              //Financed by government
    {
      EntryCosts+=capitalRecovered2.Sum();
      BankTransfer+=capitalRecovered2.Sum();
      for(i=1; i<=NB; i++)
      {
        Inflows(i)+=capitalRecovered2(i);
      }
    }
    Deposits_h(1)=0;
    for (const string& cl:classes_mh){
      Deposits_mh[cl]=0;
    }
  }
  else
  {
    std::cerr << "\n\n ERROR: Household deposits are negative in period " << t << endl;
    Errors << "\n Household deposits are negative in period " << t << endl;
    exit(EXIT_FAILURE);
  }

  //Transfer of deposits to C-firms
  for(j=1; j<=N2; j++)
  {
    if (exiting_2(j)==1)
    {
      N(1,j)=0;
      N(2,j)=0;
      //Any inventories of exiting firms are destroyed and hence booked as a loss
      Injection_2(j)=-Inventories(1,j);
      Inventories(1,j)=0;
      Inventories(2,j)=0;
      //Households give entering C-firms a transfer of deposits
      receivingBank=BankingSupplier_2(j);
      multip_entry=ran1(p_seed);
      multip_entry=w2inf+multip_entry*(w2sup-w2inf);
      if(Deposits_h(1)>=multip_entry*mD2) //Households have sufficient deposits
      {
        injection=multip_entry*mD2;
        Deposits_h(1)-=injection;
        for (const string& cl:classes_mh){ //Injection from each household class
          injection_mh[cl]=injection*Entry_financing_sh_mh[cl];
          Deposits_mh[cl](1)-=injection_mh[cl];
          if (Deposits_mh[cl](1)<0){
            std::cerr << "\n\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
            Errors << "\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
            //exit(EXIT_FAILURE);
          }
        }
        FirmTransfers+=injection;
        FirmTransfers_2+=injection;
        Injection_2(j)+=injection;
        for(i=1; i<=NB; i++)
        {
          Deposits_hb(1,i)-=DepositShare_h(i)*injection;
          Deposits_b(1,i)-=DepositShare_h(i)*injection;
          Outflows(i)+=DepositShare_h(i)*injection;
        }
        Deposits_2(1,j)=injection;
        Deposits_b(1,receivingBank)+=injection;
        Inflows(receivingBank)+=injection;
      }
      else
      {
        //If households cannot finance C-firm entry, this is done by government or banks as in the case of K-firms
        std::cerr << "\n\n Households cannot finance C-firm " << j <<" entry in period " << t << endl;
        Errors << "\n Households cannot finance C-firm " << j <<" entry in period " << t << endl;
        injection=0;  //In case Deposits_h not positive, this has to be 0 (all financed by either government or banks)
        if(Deposits_h(1)>0) //Households finance as much as they can
        {
          injection=Deposits_h(1);
          Deposits_h(1)-=injection;
          for (const string& cl:classes_mh){ //Injection from each household class
            injection_mh[cl]=injection*Entry_financing_sh_mh[cl];
            Deposits_mh[cl](1)-=injection_mh[cl];
            if (Deposits_mh[cl](1)<0){
              std::cerr << "\n\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
              Errors << "\n Household class " << cl << " cannot finance C-firm entry in period " << t << "CHANGE CODE" << endl;
              //exit(EXIT_FAILURE);
            }
          }
          FirmTransfers+=injection;
          FirmTransfers_2+=injection;
          Injection_2(j)+=injection;
          for(i=1; i<=NB; i++)
          {
            Deposits_b(1,i)-=Deposits_hb(1,i);
            Outflows(i)+=Deposits_hb(1,i);
            Deposits_hb(1,i)=0;
            DepositShare_h(i)=(NL_1(i)+NL_2(i))/(N1+N2);
          }
          Deposits_2(1,j)=injection;
          Deposits_b(1,receivingBank)+=injection;
          Inflows(receivingBank)+=injection;
        }

        injection2=multip_entry*mD2-injection;  //Injection left to be financed
        if(flag_entry==1) //Financed by banks
        {
          LossEntry_b(receivingBank)+=injection2;
        }
        else //Financed by government
        {
          EntryCosts+=injection2;
          FirmTransfers_2+=injection2;
          Inflows(receivingBank)+=injection2;
        }
        Injection_2(j)+=injection2;
        Deposits_2(1,j)+=injection2;
        Deposits_b(1,receivingBank)+=injection2;
      }

      //Transfer of second-hand capital stock
      n_mach(j)=0;
      K(j)=0;
        //First subtract the capital stock previously held by the exiting firm
      Injection_2(j)-=(CapitalStock(1,j)+deltaCapitalStock(1,j));
      CapitalStock(1,j)=0;
        //Clear the exiting firms' entries in the frequency arrays
      n_mach_resid=n_mach_entry(j);
      for (i=1; i<=N1; i++)
      {
        for (tt=t0; tt<=t; tt++)
        {
          g[tt-1][i-1][j-1]=0;
          gtemp[tt-1][i-1][j-1]=0;
          g_c[tt-1][i-1][j-1]=0;
          g_c2[tt-1][i-1][j-1]=0;
          g_c3[tt-1][i-1][j-1]=0;
          age[tt-1][i-1][j-1]=0;
        }
      }
        //Give the entering firm an initial capital stock from the pool of second-hand capital
      while(n_mach_resid>0)
      {
        rni=int(ran1(p_seed)*N1*N2)%N1+1;
        for (tt=t0; tt<=t; tt++)
        {
          if(g_secondhand[tt-1][rni-1]>0)
          {
            if(g_secondhand[tt-1][rni-1]>=n_mach_resid)
            {
              g[tt-1][rni-1][j-1]+=n_mach_resid;
              gtemp[tt-1][rni-1][j-1]+=n_mach_resid;
              g_c[tt-1][rni-1][j-1]+=n_mach_resid;
              g_c2[tt-1][rni-1][j-1]+=n_mach_resid;
              g_c3[tt-1][rni-1][j-1]+=n_mach_resid;
              age[tt-1][rni-1][j-1]=age_secondhand[tt-1][rni-1];
              n_mach(j)+=n_mach_resid;
              CapitalStock(1,j)+=n_mach_resid*g_secondhand_p[tt-1][rni-1];
              g_price[tt-1][rni-1][j-1]=g_secondhand_p[tt-1][rni-1];
              K(j)+=n_mach_resid*dim_mach;
              g_secondhand[tt-1][rni-1]-=n_mach_resid;
              n_mach_resid=0;
            }
            else
            {
              g[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              gtemp[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              g_c[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              g_c2[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              g_c3[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              age[tt-1][rni-1][j-1]=age_secondhand[tt-1][rni-1];
              n_mach(j)+=g_secondhand[tt-1][rni-1];
              CapitalStock(1,j)+=g_secondhand[tt-1][rni-1]*g_secondhand_p[tt-1][rni-1];
              g_price[tt-1][rni-1][j-1]=g_secondhand_p[tt-1][rni-1];
              K(j)+=g_secondhand[tt-1][rni-1]*dim_mach;
              n_mach_resid-=g_secondhand[tt-1][rni-1];
              g_secondhand[tt-1][rni-1]=0;
            }
          }
        }
      }
        //The recovered capital stock is added to the net worth injection
      Injection_2(j)+=CapitalStock(1,j);
      rni=int(ran1(p_seed)*N1*N2)%N1+1;
      fornit(j)=rni;
      Match(j,rni)=1;
      EI(1,j)=0;
      scrap_age(j)=0;
      deltaCapitalStock(1,j)=0;

      //Set the newly entering firm's cost, mark-up and price
      c2(j)=0;
      c2p(j)=0;
      for (i=1; i<=N1; i++)
      {
        for (tt=t0; tt<=t; tt++)
        {
          c2(j)+=(w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/((1-shocks_eneff2(j))*A_en(tt,i))+t_CO2*A_ef(tt,i)/((1-shocks_eneff2(j))*A_en(tt,i)))*g[tt-1][i-1][j-1]/n_mach(j);
          if(pass_2(j)==1)
          {
            c2p(j)+=(w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/((1-shocks_eneff2(j))*A_en(tt,i))+t_CO2*A_ef(tt,i)/((1-shocks_eneff2(j))*A_en(tt,i)))*g[tt-1][i-1][j-1]/n_mach(j);
          }
          else
          {
            c2p(j)+=(w_tot_for_1_wr_mh(2)/((1-shocks_labprod2(j))*A(tt,i))+(c_en_preshock+pass_2(j)*(c_en(2)-c_en_preshock))/((1-shocks_eneff2(j))*A_en(tt,i))+t_CO2*A_ef(tt,i)/((1-shocks_eneff2(j))*A_en(tt,i)))*g[tt-1][i-1][j-1]/n_mach(j);
          }
        }
      }
      mu2(1,j)=mi2;
      p2(j)=(1+mu2(1,j))*c2p(j);
      p2_entry+=p2(j);
      DebtService_2(1,j)=0;

      //Update ownership based on contribution by each household class weighted on capital stock and injection
      for (const string& cl:classes_mh){
        double ShareCapitalContributionClass;
        if (capitalRecoveredTot!=0){
          ShareCapitalContributionClass = capitalRecoveredTot_mh[cl] / capitalRecoveredTot;
        } else{
          ShareCapitalContributionClass = Ownership_sh_2_mh[cl];
        }
        double CapitalContributionClass = CapitalStock(1,j) * ShareCapitalContributionClass;      //How much capital passed from the households class to the K-firm j (assuming no second-hand capital from banks or government)
        Ownership_sh_2_i_mh[cl](j)=(injection_mh[cl]+CapitalContributionClass)/(injection+CapitalStock(1,j)); //Based on injection only of households (as if neither gov or banks contributed)
      }
    }
  }

  if(exiting_2.Sum()>0)
  {
    CurrentDemand=D2.Row(1).Sum();
    if(f2_exit>0){
      p2_entry/=exiting_2.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        {
          CompEntry(j)=-p2(j)/p2_entry;
        }
      }
      CompEntry_m=CompEntry.Sum()/exiting_2.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        {
          EntryShare(j)=(1/exiting_2.Sum())*((2*omega3)/(1+exp((-chi)*((CompEntry(j)-CompEntry_m)/CompEntry_m)))+(1-omega3));
        }
      }
      EntryShare=EntryShare/EntryShare.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        {
          f2(1,j)=max(f2_entry_min,EntryShare(j)*f2_exit);
          f2(2,j)=f2(1,j);
          f2(3,j)=f2(1,j);
        }
      }
      // Market shares are normalised after entering C-firms get a market share
      ftot(1)=f2.Row(1).Sum();
      ftot(2)=f2.Row(2).Sum();
      ftot(3)=f2.Row(3).Sum();

      for(j=1; j<=N2; j++)
      {
        f2(1,j)/=ftot(1);
        f2(2,j)/=ftot(2);
        f2(3,j)/=ftot(3);
        if (exiting_2(j)==1)
        {
          D2(1,j)=min(K(j),f2(1,j)*CurrentDemand);
          l2(j)=1+(f2(1,j)*CurrentDemand-D2(1,j));
          De(j)=D2(1,j);
          S2(1,j)=p2(j)*D2(1,j);
          mol(j)=S2(1,j)-D2(1,j)*c2(j);
          if(mol(j)<0)
          {
            mol(j)=0;
          }
        }
      }
    }else{
      f2_exit=f2_entry_min*exiting_2.Sum();
      p2_entry/=exiting_2.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        {
          CompEntry(j)=-p2(j)/p2_entry;
        }
      }
      CompEntry_m=CompEntry.Sum()/exiting_2.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        {
          EntryShare(j)=(1/exiting_2.Sum())*((2*omega3)/(1+exp((-chi)*((CompEntry(j)-CompEntry_m)/CompEntry_m)))+(1-omega3));
        }
      }
      EntryShare=EntryShare/EntryShare.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        {
          f2(1,j)=max(f2_entry_min,EntryShare(j)*f2_exit);
          f2(2,j)=f2(1,j);
          f2(3,j)=f2(1,j);
        }
      }
      // Market shares are normalised after entering C-firms get a market share
      ftot(1)=f2.Row(1).Sum();
      ftot(2)=f2.Row(2).Sum();
      ftot(3)=f2.Row(3).Sum();

      for(j=1; j<=N2; j++)
      {
        f2(1,j)/=ftot(1);
        f2(2,j)/=ftot(2);
        f2(3,j)/=ftot(3);
        if (exiting_2(j)==1)
        {
          D2(1,j)=min(K(j),f2(1,j)*CurrentDemand);
          l2(j)=1+(f2(1,j)*CurrentDemand-D2(1,j));
          De(j)=D2(1,j);
          S2(1,j)=p2(j)*D2(1,j);
          mol(j)=S2(1,j)-D2(1,j)*c2(j);
          if(mol(j)<0)
          {
            mol(j)=0;
          }
        }
      }
    }
  }

  //Update C-firm K-firm network
	nclient=0;
	for (i=1; i<=N1; i++)
	{
		for (j=1; j<=N2; j++)
		{
			nclient(i)+=Match(j,i);
		}
	}
	for (i=1; i<=N1; i++)
	{
		if (nclient(i) == 0)
		{
			stepbis=step;
			while (stepbis > 0)
			{
				rni=int(ran1(p_seed)*N1*N2)%N2+1;
				if (Match(rni,i) == 0)
				{
					Match(rni,i)=1;
					stepbis--;
				}
			}
		}
	}

  nclient=0;
	for (i=1; i<=N1; i++)
	{
		for (j=1; j<=N2; j++)
		{
			nclient(i)+=Match(j,i);
		}
	}
}

void TECHANGEND(void)
{
  //Endogenous technological change
  Inn=0;
  Imm=0;

  A1inn=0.00001;
  A1pinn=0.00001;
  A1imm=0.00001;
  A1pimm=0.00001;

  EE_inn=0.00001;
  EEp_inn=0.00001;
  EE_imm=0.00001;
  EEp_imm=0.00001;

  EF_inn=100000;
  EFp_inn=100000;
  EF_imm=100000;
  EFp_imm=100000;

  A1_en_max=0;
  A1p_en_max=0;
  A1_ef_max=100000;
  A1p_ef_max=100000;


  //Calculate average price K-firms
  p1_avg = 0;
  for (i=1; i<=N1; i++){
      p1_avg += (p1(i) * S1(i) / S1.Sum());
  }
  
  //Calculate what needed to determine subsidies
  if (flag_environmental_subsidies_C_firms==1 && t>=(t_regime_shifts-1)){
    //Update average carbon intensity over K-firms and on average and average price across K-firms
    A1_ci_avg = 0;
    A1_ci_min = 1000000;
    for (i=1; i<=N1; i++)
    {
      A1_ci(i) = A1_ef(i)/A1_en(i);
      if (A1_ci(i) < A1_ci_min){ //Store minimum carbon intensity for calculation of subsidy
        A1_ci_min = A1_ci(i);
      }
      A1_ci_avg += A1_ci(i);
      
    }
    A1_ci_avg /= N1;

    //Calculate max subsidy
    EnvSubsidiesPerMachines_Max = p1_avg * env_subsidy_per_machine;
  } else if (flag_environmental_subsidies_C_firms==2 && t>=(t_regime_shifts-1)){
    //Update average carbon intensity over K-firms and on average and average price across K-firms
    A1_ci_avg = 0;
    A1_ci_max = 0;
    for (i=1; i<=N1; i++)
    {
      A1_ci(i) = A1_ef(i)/A1_en(i);
      if (A1_ci(i) > A1_ci_max && A1_ci(i) < 100000){ //Store maximmum carbon intensity for calculation of subsidy. < 100000 to avoid choosing machines not innovated/imitated with super high carbon intensity assinged by default (even if here should not be a problem since ones from previous period)
        A1_ci_max = A1_ci(i);
      }
      A1_ci_avg += A1_ci(i);
    }
    A1_ci_avg /= N1;

    //Calculate max subsidy
    EnvSubsidiesPerMachines_Max = p1_avg * env_subsidy_per_machine;
  }

  for (i=1; i<=N1; i++)
  {
    //K-firms determine R&D spending
    if(S1_pre(i)>S1_post(i) && shocks_output1(i)>0)
    {
      if(exiting_1(i)==0)
      {
        RD(1,i)=max(0.0,nu*(S1_post(i)-(S1_pre(i)-S1_post(i))));
      }
      else
      {
        RD(1,i)=max(0.0,nu*min(S1(i),(S1_post(i)-(S1_pre(i)-S1_post(i)))));
      }
    }
    else
    {
      RD(1,i)=nu*S1(i);
      if (S1(i)==0)
      {
        RD(1,i)=RD(2,i);

        if (nclient(i) < 1)
        {
          std::cerr << "\n\n ERROR: nclient < 1 for K-firm " << i << " in period " << t << endl;
          Errors << "\n nclient < 1 for K-firm " << i << " in period " << t << endl;
          exit(EXIT_FAILURE);
        }
      }
    }

    //K-firms determine labour demand for R&D
    Ld1_rd_i_mh["wr"](i)=RD(1,i)/w_tot_for_1_wr_mh(1); //Demand based on entire wages to pay for each Worker
    for (const string& cl:classes_mh){
      Ld1_rd_i_mh[cl](i)=Ld1_rd_i_mh["wr"](i)*ld_ratios_mh[cl];
    }
    LD1_rd_mh=umap_sum_each_key(Ld1_rd_i_mh); //Total L demand for each household class

    //Divide between innovation and imitation
    RDin(i)=Ld1_rd_i_mh["wr"](i)*xi;
    RDim(i)=Ld1_rd_i_mh["wr"](i)*(1-xi);
    //If a shock to resources devoted to R&D applies, reduce them correspondingly
    if(flag_RDshocks==4)
    {
      RDin(i)=RDin(i)*(1-shocks_rd(i));
      RDim(i)=RDim(i)*(1-shocks_rd(i));
    }

    //Determine whether firm innovates and/or imitates (also taking into account possible shocks to R&D effectiveness)
    parber=1-exp(-o1*RDin(i));
    if(flag_RDshocks==3)
    {
      parber=(1-shocks_rd(i))*parber;
    }
    Inn(i)=bnldev(parber,1,p_seed);

    parber=1-exp(-o2*RDim(i));
    if(flag_RDshocks==3)
    {
      parber=(1-shocks_rd(i))*parber;
    }
    Imm(i)=bnldev(parber,1,p_seed);

    if (Inn(i) == 1)
    {
      //If firm innovates, determine characteristics of new technology, possibly taking into account shocks

      //Energy efficiency
      if(flag_RDshocks==1)
      {
        b_a2_shock=(1-shocks_rd(i))*b_a2;
        if(b_a2_shock<=0)
        {
            b_a2_shock=0.0001;
        }
        rnd=betadev(b_a2_shock,b_b2,p_seed);
      }
      else
      {
        rnd=betadev(b_a2,b_b2,p_seed);
      }
      if(flag_RDshocks==2)
      {
        rnd=(uu31*(1+shocks_rd(i)))+rnd*((uu41*(1-shocks_rd(i)))-(uu31*(1+shocks_rd(i))));
      }
      else
      {
        rnd=uu31+rnd*(uu41-uu31);
      }
      EE_inn(i)=A1_en(i)*(1+rnd);

      if(flag_RDshocks==1)
      {
        b_a2_shock=(1-shocks_rd(i))*b_a2;
        if(b_a2_shock<=0)
        {
            b_a2_shock=0.0001;
        }
        rnd=betadev(b_a2_shock,b_b2,p_seed);
      }
      else
      {
        rnd=betadev(b_a2,b_b2,p_seed);
      }
      if(flag_RDshocks==2)
      {
        rnd=(uu32*(1+shocks_rd(i)))+rnd*((uu42*(1-shocks_rd(i)))-(uu32*(1+shocks_rd(i))));
      }
      else
      {
        rnd=uu32+rnd*(uu42-uu32);
      }
      EEp_inn(i)=A1p_en(i)*(1+rnd);

      //Environmental friendliness
      if(flag_RDshocks==1)
      {
        b_a3_shock=(1-shocks_rd(i))*b_a3;
        if(b_a3_shock<=0)
        {
            b_a3_shock=0.0001;
        }
        rnd=betadev(b_a3_shock,b_b3,p_seed);
      }
      else
      {
        rnd=betadev(b_a3,b_b3,p_seed);
      }

      if(flag_RDshocks==2)
      {
        rnd=(uu51*(1+shocks_rd(i)))+rnd*((uu61*(1-shocks_rd(i)))-(uu51*(1+shocks_rd(i))));
      }
      else
      {
        rnd=uu51+rnd*(uu61-uu51);
      }
      EF_inn(i)=A1_ef(i)*(1-rnd);

      if(flag_RDshocks==1)
      {
        b_a3_shock=(1-shocks_rd(i))*b_a3;
        if(b_a3_shock<=0)
        {
            b_a3_shock=0.0001;
        }
        rnd=betadev(b_a3_shock,b_b3,p_seed);
      }
      else
      {
        rnd=betadev(b_a3,b_b3,p_seed);
      }
      if(flag_RDshocks==2)
      {
        rnd=(uu52*(1+shocks_rd(i)))+rnd*((uu62*(1-shocks_rd(i)))-(uu52*(1+shocks_rd(i))));
      }
      else
      {
        rnd=uu52+rnd*(uu62-uu52);
      }
      EFp_inn(i)=A1p_ef(i)*(1-rnd);

  
      //Labour productivity
      if(flag_RDshocks==1)
      {
        b_a1_shock=(1-shocks_rd(i))*b_a11;
        if(b_a1_shock<=0)
        {
            b_a1_shock=0.0001;
        }
        rnd=betadev(b_a1_shock,b_b11,p_seed);
      }
      else
      {
        rnd=betadev(b_a11,b_b11,p_seed);
      }

      if(flag_RDshocks==2)
      {
        rnd=(uu11*(1+shocks_rd(i)))+rnd*((uu21*(1-shocks_rd(i)))-(uu11*(1+shocks_rd(i))));
      }
      else
      {
        rnd=uu11+rnd*(uu21-uu11);
      }

      if (flag_correlate_prod_and_green_tech==1){ //Case in which we add an influence of carbon intensity on labour productivity
        double CI_old = A1_ef(i)/A1_en(i);        //Carbon intensity of machine produced
        double CI_inn = EF_inn(i)/EE_inn(i);      //Carbon intensity of innovated machine
        double CI_g = (CI_inn - CI_old)/CI_old;   //Relative difference in carbon intensity
        rnd -= correlation_prod_and_green_tech * CI_g;
      } else if (flag_correlate_prod_and_green_tech==2 & t >= t_regime_shifts){
        double CI_old = A1_ef(i)/A1_en(i);        //Carbon intensity of machine produced
        double CI_inn = EF_inn(i)/EE_inn(i);      //Carbon intensity of innovated machine
        double CI_g = (CI_inn - CI_old)/CI_old;   //Relative difference in carbon intensity
        rnd -= correlation_prod_and_green_tech * CI_g;
      }
      A1inn(i)=A1(i)*(1+rnd);

      if(flag_RDshocks==1)
      {
        b_a1_shock=(1-shocks_rd(i))*b_a12;
        if(b_a1_shock<=0)
        {
            b_a1_shock=0.0001;
        }
        rnd=betadev(b_a1_shock,b_b12,p_seed);
      }
      else
      {
        rnd=betadev(b_a12,b_b12,p_seed);
      }
      if(flag_RDshocks==2)
      {
        rnd=(uu12*(1+shocks_rd(i)))+rnd*((uu22*(1-shocks_rd(i)))-(uu12*(1+shocks_rd(i))));
      }
      else
      {
        rnd=uu12+rnd*(uu22-uu12);
      }
      A1pinn(i)=A1p(i)*(1+rnd);

      if ((flag_correlate_prod_and_green_tech==3 || flag_correlate_prod_and_green_tech==4) & t >= t_regime_shifts){ // In case, override the previously calculated changes
        //Labour productivity
        rnd=betadev(b_a11,b_b11,p_seed);
        rnd=uu11+rnd*(uu21-uu11);

        double rel_change_lp;
        double change_ee;
        double change_ef;
        if (rnd >= 0){
          rel_change_lp = (rnd - 0)/(uu21 - 0); //Calculate change in LP compared to max change in the relative direction
          change_ee = rel_change_lp * uu41; //Energy efficiency
          change_ef = rel_change_lp * uu61; //Environmental friendliness (opposite direction, but we take same bound since then changing by opposite direction)
        }else{
          rel_change_lp = (0 - rnd)/(0 - uu11);
          change_ee = rel_change_lp * uu31; //Energy efficiency
          change_ef = rel_change_lp * uu51; //Environmental friendliness (opposite direction, but we take same bound since then changing by opposite direction))
        }
        
        if (flag_correlate_prod_and_green_tech==4){
          rnd*=rs_uu_lp;
          change_ee*=rs_uu_ee;
          change_ef*=rs_uu_ef;
        }
        
        A1inn(i)=A1(i)*(1+rnd);
        EE_inn(i)=A1_en(i)*(1+change_ee);
        EF_inn(i)=A1_ef(i)*(1-change_ef);
      }

      if (A1pinn(i)==0 || A1inn(i)==0 || A1p(i)==0 || A1(i)==0)
      {
        std::cerr << "\n\n ERROR: A1=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
      if (EEp_inn(i)==0 || EE_inn(i)==0 || A1p_en(i)==0 || A1_en(i)==0)
      {
        std::cerr << "\n\n ERROR: A1_en=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1_en=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
    }


    if (Imm(i) == 1)
    {
      //If K-firm imitates, determine which other firm's technology it will imitate
      Tdtot=0;
      for (ii=1; ii<=N1; ii++)
      {
          Td.element(ii)=sqrt(((A1(ii)-A1(i))*(A1(ii)-A1(i))) + ((A1p(ii)-A1p(i))*(A1p(ii)-A1p(i))) + ((A1_en(ii)-A1_en(i))*(A1_en(ii)-A1_en(i))) + ((A1_ef(ii)-A1_ef(i))*(A1_ef(ii)-A1_ef(i))) + ((A1p_en(ii)-A1p_en(i))*(A1p_en(ii)-A1p_en(i))) + ((A1p_ef(ii)-A1p_ef(i))*(A1p_ef(ii)-A1p_ef(i))));
          if (Td.element(ii)>0)
          {
            Td.element(ii)=1/Td.element(ii);
          }
          else
          {
            Td.element(ii)=0;
          }
          Tdtot+=Td.element(ii);
      }
      for (ii=1; ii<=N1; ii++)
      {
        Td.element(ii)/=Tdtot;
        Td.element(ii)+=Td.element(ii-1);
      }
      rnd=ran1(p_seed);
      for (ii=1; ii<=N1; ii++)
      {
        if (rnd <= Td.element(ii) && rnd > Td.element(ii-1))
        {
          A1imm(i)=A1(ii);
          A1pimm(i)=A1p(ii);
          EE_imm(i)=A1_en(ii);
          EEp_imm(i)=A1p_en(ii);
          EF_imm(i)=A1_ef(ii);
          EFp_imm(i)=A1p_ef(ii);
        }
      }

      if (A1pimm(i)==0 || A1imm(i)==0 || A1p(i)==0 || A1(i)==0)
      {
        std::cerr << "\n\n ERROR: A1=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
      if (EEp_imm(i)==0 || EE_imm(i)==0 || A1p_en(i)==0 || A1_en(i)==0)
      {
        std::cerr << "\n\n ERROR: A1_en=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1_en=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
      if (EFp_imm(i)==0 || EF_imm(i)==0 || A1p_ef(i)==0 || A1_ef(i)==0)
      {
        std::cerr << "\n\n ERROR: A1_ef=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1_ef=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
    }

    //Decision on technology adoption (and eventual subsidy determination)
    double EnvSubsidiesPerMachine_imm;
    double EnvSubsidiesPerMachine_inn;
    if (flag_environmental_subsidies_C_firms==1 && t>=(t_regime_shifts-1)){
      //Calculate carbon intensity of innovated and imitated vintages
      double CI_imm = EF_imm(i)/EE_imm(i);
      double CI_inn = EF_inn(i)/EE_inn(i);
      
      //Calculate subsidy for current, innovated and imitated vintages
      if (A1_ci(i) < A1_ci_avg){
        EnvSubsidiesPerMachine_1_i(i) = EnvSubsidiesPerMachines_Max * (A1_ci_avg - A1_ci(i))/(A1_ci_avg - A1_ci_min);
      } else{
        EnvSubsidiesPerMachine_1_i(i) = 0;
      }

      if (CI_imm < A1_ci_avg){
        EnvSubsidiesPerMachine_imm = EnvSubsidiesPerMachines_Max * (A1_ci_avg - CI_imm)/(A1_ci_avg - A1_ci_min);
      } else{
        EnvSubsidiesPerMachine_imm = 0;
      }

      if (CI_inn < A1_ci_avg){
        EnvSubsidiesPerMachine_inn = EnvSubsidiesPerMachines_Max * (A1_ci_avg - CI_inn)/(A1_ci_avg - A1_ci_min);
      } else{
        EnvSubsidiesPerMachine_inn = 0;
      }
    } else if(flag_environmental_subsidies_C_firms==2 && t>=(t_regime_shifts-1)){ //Subsidies for dirty machines
      //Calculate carbon intensity of innovated and imitated vintages
      double CI_imm = EF_imm(i)/EE_imm(i);
      double CI_inn = EF_inn(i)/EE_inn(i);
      
      //Calculate subsidy for current, innovated and imitated vintages
      if (A1_ci(i) > A1_ci_avg && A1_ci(i) < 100000){
        EnvSubsidiesPerMachine_1_i(i) = EnvSubsidiesPerMachines_Max * (A1_ci(i) - A1_ci_avg)/(A1_ci_max - A1_ci_avg);
      } else{
        EnvSubsidiesPerMachine_1_i(i) = 0;
      }

      if (CI_imm > A1_ci_avg && CI_imm < 100000){
        EnvSubsidiesPerMachine_imm = EnvSubsidiesPerMachines_Max * (CI_imm - A1_ci_avg)/(A1_ci_max - A1_ci_avg);
      } else{
        EnvSubsidiesPerMachine_imm = 0;
      }

      if (CI_inn > A1_ci_avg && CI_inn < 100000){
        EnvSubsidiesPerMachine_inn = EnvSubsidiesPerMachines_Max * (CI_inn - A1_ci_avg)/(A1_ci_max - A1_ci_avg);
      } else{
        EnvSubsidiesPerMachine_inn = 0;
      }
    } else{
      EnvSubsidiesPerMachine_imm = 0;
      EnvSubsidiesPerMachine_inn = 0;
    }

      //Calculate attractiveness
    double Attr_imm = (1+mi1)*(w_tot_for_1_wr_mh(1)/(A1pimm(i))+c_en(1)/EEp_imm(i)+t_CO2*EFp_imm(i)/EEp_imm(i))-EnvSubsidiesPerMachine_imm + (w_tot_for_1_wr_mh(1)/A1imm(i)+c_en(1)/EE_imm(i)+t_CO2*EF_imm(i)/EE_imm(i))*b;
    double Attr_inn = (1+mi1)*(w_tot_for_1_wr_mh(1)/(A1pinn(i))+c_en(1)/EEp_inn(i)+t_CO2*EFp_inn(i)/EEp_inn(i))-EnvSubsidiesPerMachine_inn + (w_tot_for_1_wr_mh(1)/A1inn(i)+c_en(1)/EE_inn(i)+t_CO2*EF_inn(i)/EE_inn(i))*b;
    double Attr_cur = (1+mi1)*(w_tot_for_1_wr_mh(1)/(A1p(i))+c_en(1)/A1p_en(i)+t_CO2*A1p_ef(i)/A1p_en(i))-EnvSubsidiesPerMachine_1_i(i)    + (w_tot_for_1_wr_mh(1)/A1(i)+c_en(1)/A1_en(i)+t_CO2*A1_ef(i)/A1_en(i))*b;

      //If the imitated technology is superior, adopt it
    if ( Attr_imm < Attr_cur)
    {
      A1(i)=A1imm(i);
      A1p(i)=A1pimm(i);
      A1_en(i)=EE_imm(i);
      A1p_en(i)=EEp_imm(i);
      A1_ef(i)=EF_imm(i);
      A1p_ef(i)=EFp_imm(i);

      if (EnvSubsidiesPerMachine_imm > EnvSubsidiesPerMachines_Max){
        EnvSubsidiesPerMachine_1_i(i) = EnvSubsidiesPerMachines_Max;
      } else{
        EnvSubsidiesPerMachine_1_i(i) = EnvSubsidiesPerMachine_imm;
      }
    }

      //If the innovated technology is superior, adopt it
    if (Attr_inn < Attr_cur)
    {
      A1(i)=A1inn(i);
      A1p(i)=A1pinn(i);
      A1_en(i)=EE_inn(i);
      A1p_en(i)=EEp_inn(i);
      A1_ef(i)=EF_inn(i);
      A1p_ef(i)=EFp_inn(i);

      if (EnvSubsidiesPerMachine_inn > EnvSubsidiesPerMachines_Max){
        EnvSubsidiesPerMachine_1_i(i) = EnvSubsidiesPerMachines_Max;
      } else{
        EnvSubsidiesPerMachine_1_i(i) = EnvSubsidiesPerMachine_inn;
      }
    }
    
    //Update carbon intensity (just for eventual reporting)
    if ((flag_environmental_subsidies_C_firms==1 || flag_environmental_subsidies_C_firms==2) && t>=t_regime_shifts){
      A1_ci(i) = A1_ef(i)/A1_en(i);
    }

    //Update the productivity matrices
    if (t < T )
    {
      A(t+1,i)=A1(i);
      A_en(t+1,i)=A1_en(i);
      A_ef(t+1,i)=A1_ef(i);
    }
  }

  //Determine the best technologies in the system post-R&D
  A1top=A1(1);
  A1ptop=A1p(1);
  A1_en_top=A1_en(1);
  A1_ef_top=A1_ef(1);
  A1p_en_top=A1p_en(1);
  A1p_ef_top=A1p_ef(1);

  for (i=1; i<=N1; i++)
  {
    if (A1(i) > A1top)
    {
      A1top=A1(i);
    }
    if (A1p(i) > A1ptop)
    {
      A1ptop=A1p(i);
    }
    if (A1_en(i) > A1_en_top)
    {
      A1_en_top=A1_en(i);
    }
    if (A1p_en(i) > A1p_en_top)
    {
      A1p_en_top=A1p_en(i);
    }
    if (A1_ef(i) < A1_ef_top)
    {
      A1_ef_top=A1_ef(i);
    }
    if (A1p_ef(i) < A1p_ef_top)
    {
      A1p_ef_top=A1p_ef(i);
    }
  }


  for (i=1; i<=N1; i++)
  {
    if (A1p(i)>A1pmax)
    {
      A1pmax=A1p(i);
    }
    if (A1_en(i)>A1_en_max)
    {
      A1_en_max=A1_en(i);
    }
    if (A1p_en(i)>A1p_en_max)
    {
      A1p_en_max=A1p_en(i);
    }
    if (A1_ef(i)<A1_ef_max)
    {
      A1_ef_max=A1_ef(i);
    }
    if (A1p_ef(i)<A1p_ef_max)
    {
      A1p_ef_max=A1p_ef(i);
    }
  }
}

void ENERGYINEQUALITY(void)
{
  cons_dir_en_fp_h = 0;
  cons_indir_en_fp_h = 0;
  invest_en_fp_h = 0;
  publ_cons_en_fp_h = 0;
  en_fp_h = 0;
  for (const string& cl:classes_mh){
    //Direct energy footprint
    cons_dir_en_fp_mh[cl] = D_en_mh[cl];
    cons_dir_en_fp_h += cons_dir_en_fp_mh[cl];

    //Indirect energy footprint from consumption
    //cons_indir_en_fp_mh[cl] = Consumption_mh[cl] / Am_en_2; // ALSO HERE CAN MATCH SINGLE FIRMS WITH INVESTMENT, OTHERWISE CONTROL WOULD NOT WORK
    cons_indir_en_fp_mh[cl] = D2_en_TOT * Consumption_mh[cl]/Consumption_h;
    cons_indir_en_fp_h += cons_indir_en_fp_mh[cl];

    //Investment energy footprint
    // invest_en_fp_mh[cl] = Ownership_sh_2_mh[cl] *  / Am_en_1; //VARAIBLE FOR INVESTMENT + CAN BE DONE AT THE SINGLE FIRM LEVEL! MATCHING PROD WITH AMOUNT OF INV, OTHERWISE CONTROL WOULD NOT WORK
    invest_en_fp_mh[cl] = D1_en_TOT * Ownership_sh_2_mh[cl];
    invest_en_fp_h += invest_en_fp_mh[cl];

    //Public spending energy footprint
    publ_cons_en_fp_mh[cl] = D_en_g * LS_sh_mh[cl];
    publ_cons_en_fp_h += publ_cons_en_fp_mh[cl];

    //Total energy footprint
    en_fp_mh[cl] = cons_dir_en_fp_mh[cl] + cons_indir_en_fp_mh[cl] + invest_en_fp_mh[cl] + publ_cons_en_fp_mh[cl];
    en_fp_h += en_fp_mh[cl];
  }

  // Check that sum of energy footprints is equal to total energy demand
  if (t!=1){
    deviation=fabs((en_fp_h-D_en_TOT(1))/D_en_TOT(1));
    if(deviation>1.e-03) //Not using tolerance because often values of 1.e-05, but it's fine
    {
      std::cerr<<"Total households energy footprint differs from total energy demand in period " << t << endl;
      Errors << "\n Total households energy footprint differs from total energy demand in period " << t << endl;
    }
  }
}

void CARBONINEQUALITY(void)
{
  cons_dir_carb_fp_h = 0;
  cons_indir_carb_fp_h = 0;
  invest_carb_fp_h = 0;
  publ_cons_carb_fp_h = 0;
  carb_fp_h = 0;
  double Emiss1_TOT_and_en = Emiss1_TOT + Emiss_en * D1_en_TOT / D_en_TOT(1); //Total emissions of K-sector including also the ones embedded in its energy demand
  double Emiss2_TOT_and_en = Emiss2_TOT + Emiss_en * D2_en_TOT / D_en_TOT(1); //Total emissions of C-sector including also the ones embedded in its energy demand
  double Emiss_g_TOT = Emiss_en * D_en_g / D_en_TOT(1) + Emiss2_TOT_and_en * Consumption_g/Consumption; //Total emissions of government
  for (const string& cl:classes_mh){
    //Direct carbon footprint
    cons_dir_carb_fp_mh[cl] = Emiss_en * D_en_mh[cl] / D_en_TOT(1);
    cons_dir_carb_fp_h += cons_dir_carb_fp_mh[cl];

    //Indirect carbon footprint from consumption (accounting also portion from energy sector proportionally to sector energy demand share)
    //cons_indir_carb_fp_mh[cl] = Consumption_mh[cl] / Am_carb_2; // ALSO HERE CAN MATCH SINGLE FIRMS WITH INVESTMENT, OTHERWISE CONTROL WOULD NOT WORK
    cons_indir_carb_fp_mh[cl] = Emiss2_TOT_and_en * Consumption_h/Consumption * Consumption_mh[cl]/Consumption_h;
    cons_indir_carb_fp_h += cons_indir_carb_fp_mh[cl];

    //Investment carbon footprint
    // invest_carb_fp_mh[cl] = Ownership_sh_2_mh[cl] *  / Am_carb_1; //VARIABLE FOR INVESTMENT + CAN BE DONE AT THE SINGLE FIRM LEVEL! MATCHING PROD WITH AMOUNT OF INV, OTHERWISE CONTROL WOULD NOT WORK
    invest_carb_fp_mh[cl] = Emiss1_TOT_and_en * Ownership_sh_2_mh[cl];
    invest_carb_fp_h += invest_carb_fp_mh[cl];

    //Public spending carbon footprint
    publ_cons_carb_fp_mh[cl] = Emiss_g_TOT * LS_sh_mh[cl];
    publ_cons_carb_fp_h += publ_cons_carb_fp_mh[cl];

    //Total carbon footprint
    carb_fp_mh[cl] = cons_dir_carb_fp_mh[cl] + cons_indir_carb_fp_mh[cl] + invest_carb_fp_mh[cl] + publ_cons_carb_fp_mh[cl];
    carb_fp_h += carb_fp_mh[cl];
  }

  // Check that sum of carbon footprints is equal to total emissions
  if (t!=1){
    deviation=fabs((carb_fp_h-Emiss_TOT(1))/Emiss_TOT(1));
    if(deviation>1.e-03) //Not using tolerance becuse often values of 1.e-06, but it's fine
    {
      std::cerr<<"Total households carbon footprint differs from total carbon demand in period " << t << endl;
      Errors << "\n Total households carbon footprint differs from total carbon demand in period " << t << endl;
    }
  }
}

void DEPOSITCHECK(void)
{
  for(i=1; i<=NB; i++)
  {
    deviation=fabs(DepositShare_e(i)-Deposits_eb(1,i)/Deposits_eb.Row(1).Sum());
    if(deviation>tolerance && Deposits_eb(1,i)>tolerance)
    {
      std::cerr<<"Share error Deposits_eb for bank " << i << " in period " << t << endl;
      Errors << "\n Share error Deposits_eb for bank " << i << " in period " << t << endl;
    }
    deviation=fabs(DepositShare_h(i)-Deposits_hb(1,i)/Deposits_hb.Row(1).Sum());
    if(deviation>tolerance)
    {
      std::cerr<<"Share error Deposits_hb for bank " << i << " in period " << t << endl;
      Errors << "\n Share error Deposits_hb for bank " << i << " in period " << t << endl;
    }
    // Check that for each bank total deposits - deposits held for households - deposits held for banks = deposits held for C- & K-firms
    DepositsCheck_1=Deposits_b(1,i)-Deposits_hb(1,i)-Deposits_eb(1,i);
    DepositsCheck_2=0;
    for(j=1; j<=N1; j++)
    {
      if(BankMatch_1(j,i)==1)
      {
        DepositsCheck_2+=Deposits_1(1,j);
      }
    }
    for(j=1; j<=N2; j++)
    {
      if(BankMatch_2(j,i)==1)
      {
        DepositsCheck_2+=Deposits_2(1,j);
      }
    }
    deviation=fabs((DepositsCheck_1-DepositsCheck_2)/Deposits_b(1,i));
    if(deviation>tolerance)
    {
      std::cerr<<"Share error firm deposits for bank " << i << " in period " << t << endl;
      Errors << "\n Share error firm deposits for bank " << i << " in period " << t << endl;
    }
  }
}

void NEGATIVITYCHECK(void)
{
  for (j=1; j<=N2; ++j)
  {
    if(Loans_2(1,j)<0)
    {
      std::cerr<<"Error loans for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error loans for C-firm " << j << " in period " << t << endl;
    }
    if(Deposits_2(1,j)<0)
    {
      std::cerr<<"Error deposits for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error deposits for C-firm " << j << " in period " << t << endl;
    }
    if(CapitalStock(1,j)<0)
    {
      std::cerr<<"Error Capital for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Capital for C-firm " << j << " in period " << t << endl;
    }
    if(Inventories(1,j)<0)
    {
      std::cerr<<"Error Inventories for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Inventories for C-firm " << j << " in period " << t << endl;
    }
    if(Investment_2(j)<0)
    {
      std::cerr<<"Error Investment for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Investment for C-firm " << j << " in period " << t << endl;
    }
    if(Taxes_2(j)<0)
    {
      std::cerr<<"Error Taxes for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Taxes for C-firm " << j << " in period " << t << endl;
    }
    if(Wages_2_i(j)<0)
    {
      std::cerr<<"Error Wages for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Wages for C-firm " << j << " in period " << t << endl;
    }
    if(EnergyPayments_2(j)<0)
    {
      std::cerr<<"Error EnergyPayments for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error EnergyPayments for C-firm " << j << " in period " << t << endl;
    }
    if(Dividends_2_i(j)<0)
    {
      std::cerr<<"Error Dividends for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Dividends for C-firm " << j << " in period " << t << endl;
    }
    if(Bonuses_2_i(1)<0)
    {
      std::cerr<<"Error Bonuses for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Bonuses for C-firm " << j << " in period " << t << endl;
    }
    if(LoanInterest_2(j)<0)
    {
      std::cerr<<"Error LoanInterest for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error LoanInterest for C-firm " << j << " in period " << t << endl;
    }
    if(InterestDeposits_2(j)<0)
    {
      std::cerr<<"Error InterestDeposits for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error InterestDeposits for C-firm " << j << " in period " << t << endl;
    }
    if(S2(1,j)<0)
    {
      std::cerr<<"Error S2 for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error S2 for C-firm " << j << " in period " << t << endl;
    }
    if(DebtRemittances2(j)<0)
    {
      std::cerr<<"Error DebtRemittances for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error DebtRemittances for C-firm " << j << " in period " << t << endl;
    }
    if(Taxes_CO2_2(j)<0)
    {
      std::cerr<<"Error Taxes_CO2 for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Taxes_CO2 for C-firm " << j << " in period " << t << endl;
    }
  }

  for (i=1; i<=N1; i++)
  {
    if(Deposits_1(1,i)<0)
    {
      std::cerr<<"Error Deposits for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Deposits for K-firm " << i << " in period " << t << endl;
    }
    if(S1(i)<0)
    {
      std::cerr<<"Error S1 for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error S1 for K-firm " << i << " in period " << t << endl;
    }
    if(Taxes_1(i)<0)
    {
      std::cerr<<"Error Taxes for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Taxes for K-firm " << i << " in period " << t << endl;
    }
    if(Wages_1_i(i)<0)
    {
      std::cerr<<"Error Wages for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Wages for K-firm " << i << " in period " << t << endl;
    }
    if(EnergyPayments_1(i)<0)
    {
      std::cerr<<"Error EnergyPayments for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error EnergyPayments for K-firm " << i << " in period " << t << endl;
    }
    if(Dividends_1_i(i)<0)
    {
      std::cerr<<"Error Dividends for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Dividends for K-firm " << i << " in period " << t << endl;
    }
    if(Bonuses_1_i(1)<0)
    {
      std::cerr<<"Error Bonuses for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Bonuses for K-firm " << i << " in period " << t << endl;
    }
    if(InterestDeposits_1(i)<0)
    {
      std::cerr<<"Error InterestDeposits for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error InterestDeposits for K-firm " << i << " in period " << t << endl;
    }
    if(Taxes_CO2_1(i)<0)
    {
      std::cerr<<"Error Taxes_CO2 for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Taxes_CO2 for K-firm " << i << " in period " << t << endl;
    }
  }

  for (i=1; i<=NB; i++)
  {
    if(Loans_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Loans for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Loans for Bank " << i << " in period " << t << endl;
    }
    if(Deposits_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Deposits for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Deposits for Bank " << i << " in period " << t << endl;
    }
    if(Deposits_hb(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Deposits_hb for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Deposits_hb  for Bank " << i << " in period " << t << endl;
    }
    if(Deposits_eb(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Deposits_eb for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Deposits_eb for Bank " << i << " in period " << t << endl;
    }
    if(GB_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error GB for Bank " << i << " in period " << t << endl;
      Errors << "\n Error GB for Bank " << i << " in period " << t << endl;
    }
    if(Reserves_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Reserves for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Reserves for Bank " << i << " in period " << t << endl;
    }
    if(Advances_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Advances for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Advances for Bank " << i << " in period " << t << endl;
    }
    if(Taxes_b(i)<0)
    {
      std::cerr<<"Error Taxes for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Taxes for Bank " << i << " in period " << t << endl;
    }
    if(Dividends_b_i(i)<0)
    {
      std::cerr<<"Error Dividends for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Dividends for Bank " << i << " in period " << t << endl;
    }
    if(LoanInterest(i)<0)
    {
      std::cerr<<"Error LoanInterest for Bank " << i << " in period " << t << endl;
      Errors << "\n Error LoanInterest for Bank " << i << " in period " << t << endl;
    }
    if(InterestDeposits(i)<0)
    {
      std::cerr<<"Error InterestDeposits for Bank " << i << " in period " << t << endl;
      Errors << "\n Error InterestDeposits for Bank " << i << " in period " << t << endl;
    }
    if(InterestBonds_b(i)<0)
    {
      std::cerr<<"Error InterestBonds for Bank " << i << " in period " << t << endl;
      Errors << "\n Error InterestBonds for Bank " << i << " in period " << t << endl;
    }
    if(BondRepayments_b(i)<0)
    {
      std::cerr<<"Error BondRepayments for Bank " << i << " in period " << t << endl;
      Errors << "\n Error BondRepayments for Bank " << i << " in period " << t << endl;
    }
    if(InterestReserves_b(i)<0)
    {
      std::cerr<<"Error InterestReserves for Bank " << i << " in period " << t << endl;
      Errors << "\n Error InterestReserves for Bank " << i << " in period " << t << endl;
    }
    if(InterestAdvances_b(i)<0)
    {
      std::cerr<<"Error InterestAdvances for Bank " << i << " in period " << t << endl;
      Errors << "\n Error InterestAdvances for Bank " << i << " in period " << t << endl;
    }
    if(Bailout_b(i)<0)
    {
      std::cerr<<"Error Bailout for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Bailout for Bank " << i << " in period " << t << endl;
    }
  }

  //Households
  if(Deposits_h(1)<0)
  {
    std::cerr<<"Error Deposits_h in period " << t << endl;
    Errors << "\n Error Deposits_h in period " << t << endl;
  }
  if(Consumption_h<0)
  {
    std::cerr<<"Error Consumption_h in period " << t << endl;
    Errors << "\n Error Consumption_h in period " << t << endl;
  }
  if(Consumption_g<0)
  {
    std::cerr<<"Error Consumption_g in period " << t << endl;
    Errors << "\n Error Consumption_g in period " << t << endl;
  }
  if(Benefits<0)
  {
    std::cerr<<"Error Benefits in period " << t << endl;
    Errors << "\n Error Benefits in period " << t << endl;
  }
  if(Taxes_h<0)
  {
    std::cerr<<"Error Taxes_h in period " << t << endl;
    Errors << "\n Error Taxes_h in period " << t << endl;
  }
  if(Wages<0)
  {
    std::cerr<<"Error Wages in period " << t << endl;
    Errors << "\n Error Wages in period " << t << endl;
  }
  if(Dividends_h(1)<0)
  {
    std::cerr<<"Error Dividends in period " << t << endl;
    Errors << "\n Error Dividends in period " << t << endl;
  }
  if(Bonuses_h<0)
  {
    std::cerr<<"Error Bonuses_h in period " << t << endl;
    Errors << "\n Error Bonuses_h in period " << t << endl;
  }
  if(InterestDeposits_h<0)
  {
    std::cerr<<"Error InterestDeposits_h in period " << t << endl;
    Errors << "\n Error InterestDeposits_h in period " << t << endl;
  }

  //Government
  if(G<0)
  {
    std::cerr<<"Error G in period " << t << endl;
    Errors << "\n Error G in period " << t << endl;
  }
  if(Taxes_g<0)
  {
    std::cerr<<"Error Taxes in period " << t << endl;
    Errors << "\n Error Taxes in period " << t << endl;
  }
  if(Bailout<0)
  {
    std::cerr<<"Error Bailout in period " << t << endl;
    Errors << "\n Error Bailout in period " << t << endl;
  }
  if(InterestReserves<0)
  {
    std::cerr<<"Error InterestReserves in period " << t << endl;
    Errors << "\n Error InterestReserves in period " << t << endl;
  }
  if(InterestAdvances<0)
  {
    std::cerr<<"Error InterestAdvances in period " << t << endl;
    Errors << "\n Error InterestAdvances in period " << t << endl;
  }

  //Energy sector
  if(EnergyPayments<0)
  {
    std::cerr<<"Error EnergyPayments in period " << t << endl;
    Errors << "\n Error EnergyPayments in period " << t << endl;
  }
  if(Wages_en<0)
  {
    std::cerr<<"Error Wages_en in period " << t << endl;
    Errors << "\n Error Wages_en in period " << t << endl;
  }
  if(Dividends_e<0)
  {
    std::cerr<<"Error Dividends_e in period " << t << endl;
    Errors << "\n Error Dividends_e in period " << t << endl;
  }
  if(Bonuses_e<0)
  {
    std::cerr<<"Error Bonuses_e in period " << t << endl;
    Errors << "\n Error Bonuses_e in period " << t << endl;
  }
  if(InterestDeposits_e<0)
  {
    std::cerr<<"Error InterestDeposits_e in period " << t << endl;
    Errors << "\n Error InterestDeposits_e in period " << t << endl;
  }

  //Household classes
  for (const string& cl:classes_mh){
    if(Deposits_mh[cl](1)<0)
    {
      std::cerr<<"Error Deposits_mh of class " << cl <<" in period " << t << endl;
      Errors << "\n Error Deposits_mh of class " << cl <<" in period " << t << endl;
    }
    if(Consumption_mh[cl]<0)
    {
      std::cerr<<"Error Consumption_mh of class " << cl <<" in period " << t << endl;
      Errors << "\n Error Consumption_mh of class " << cl <<" in period " << t << endl;
    }
    if(Benefits_mh[cl]<0)
    {
      std::cerr<<"Error Benefits_mh of class " << cl <<" in period " << t << endl;
      Errors << "\n Error Benefits_mh of class " << cl <<" in period " << t << endl;
    }
    if(Taxes_mh[cl]<0)
    {
      std::cerr<<"Error Taxes_mh of class " << cl <<" in period " << t << endl;
      Errors << "\n Error Taxes_mh of class " << cl <<" in period " << t << endl;
    }
    if(Wages_mh[cl]<0)
    {
      std::cerr<<"Error Wages_mh of class " << cl <<" in period " << t << endl;
      Errors << "\n Error Wages_mh of class " << cl <<" in period " << t << endl;
    }
    if(Dividends_mh[cl](1)<0)
    {
      std::cerr<<"Error Dividends_mh of class " << cl <<" in period " << t << endl;
      Errors << "\n Error Dividends_mh of class " << cl <<" in period " << t << endl;
    }
    if(Bonuses_mh[cl](1)<0)
    {
      std::cerr<<"Error Bonuses_mh of class " << cl <<" in period " << t << endl;
      Errors << "\n Error Bonuses_mh of class " << cl <<" in period " << t << endl;
    }
    if(InterestDeposits_mh[cl]<0)
    {
      std::cerr<<"Error InterestDeposits_mh of class " << cl <<" in period " << t << endl;
      Errors << "\n Error InterestDeposits_mh of class " << cl <<" in period " << t << endl;
    }
  }
}

void CHECKSUMS(void)
{

  //Check that sum of household classes variables sum to households' overall
  double Deposits_mh_tot=0;
  double Consumption_mh_tot=0;
  double Wages_mh_tot=0;
  double Dividends_mh_tot=0;
  double Bonuses_mh_tot=0;
  double Taxes_w_mh_tot=0;
  double Taxes_div_mh_tot=0;
  double Taxes_bon_mh_tot=0;
  for (const string& cl:classes_mh){
    Deposits_mh_tot+=Deposits_mh[cl](1);
    Consumption_mh_tot+=Consumption_mh[cl];
    Wages_mh_tot+=Wages_mh[cl];
    Dividends_mh_tot+=Dividends_mh[cl](1);
    Bonuses_mh_tot+=Bonuses_mh[cl](1);
    Taxes_w_mh_tot+=Taxes_w_mh[cl];
    Taxes_div_mh_tot+=Taxes_div_mh[cl](1);
    Taxes_bon_mh_tot+=Taxes_bon_mh[cl](1);
  }
  deviation=fabs((Deposits_mh_tot-Deposits_h(1)))/Deposits_mh_tot;
  if(deviation>tolerance && Deposits_mh_tot>0 && Deposits_h(1)>0)
  {
    std::cerr << "\n\n ERROR: Deposits of household classes do not add up in period " << t << endl;
    Errors << "\n Deposits of household classes do not add up in period " << t << endl;
  }
  deviation=fabs((Consumption_mh_tot-Consumption_h))/Consumption_mh_tot;
  if(deviation>tolerance && Consumption_mh_tot>0 && Consumption_h>0)
  {
    std::cerr << "\n\n ERROR: Consumption of household classes do not add up in period " << t << endl;
    Errors << "\n Consumption of household classes do not add up in period " << t << endl;
  }
  deviation=fabs((Wages_mh_tot-Wages))/Wages_mh_tot;
  if(deviation>tolerance && Wages_mh_tot>0 && Wages>0)
  {
    std::cerr << "\n\n ERROR: Wages of household classes do not add up in period " << t << endl;
    Errors << "\n Wages of household classes do not add up in period " << t << endl;
  }
  deviation=fabs((Dividends_mh_tot-Dividends_h(1)))/Dividends_mh_tot;
  if(deviation>tolerance && Dividends_mh_tot>0 && Dividends_h(1)>0)
  {
    std::cerr << "\n\n ERROR: Dividends of household classes do not add up in period " << t << endl;
    Errors << "\n Dividends of household classes do not add up in period " << t << endl;
  }
  deviation=fabs((Bonuses_mh_tot-Bonuses_h))/Bonuses_mh_tot;
  if(deviation>tolerance && Bonuses_mh_tot>0 && Bonuses_h>0)
  {
    std::cerr << "\n\n ERROR: Bonuses of household classes do not add up in period " << t << endl;
    Errors << "\n Bonuses of household classes do not add up in period " << t << endl;
  }
  deviation=fabs((Taxes_w_mh_tot-Taxes_w_h))/Taxes_w_mh_tot;
  if(deviation>tolerance && Taxes_w_mh_tot>0 && Taxes_w_h>0)
  {
    std::cerr << "\n\n ERROR: Taxes on wages of household classes do not add up in period " << t << endl;
    Errors << "\n Taxes on wages of household classes do not add up in period " << t << endl;
  }
  deviation=fabs((Taxes_div_mh_tot-Taxes_div_h))/Taxes_div_mh_tot;
  if(deviation>tolerance && Taxes_div_mh_tot>0 && Taxes_div_h>0)
  {
    std::cerr << "\n\n ERROR: Taxes on dividends of household classes do not add up in period " << t << endl;
    Errors << "\n Taxes on dividends of household classes do not add up in period " << t << endl;
  }
  deviation=fabs((Taxes_bon_mh_tot-Taxes_bon_h))/Taxes_bon_mh_tot;
  if(deviation>tolerance && Taxes_bon_mh_tot>0 && Taxes_bon_h>0)
  {
    std::cerr << "\n\n ERROR: Taxes on bonuses of household classes do not add up in period " << t << endl;
    Errors << "\n Taxes on bonuses of household classes do not add up in period " << t << endl;
  }


  deviation=fabs((Deposits_h(1)-Deposits_hb.Row(1).Sum())/Deposits_hb.Row(1).Sum());
  if(deviation>tolerance && Deposits_hb.Row(1).Sum()>tolerance && Deposits_h(1)>tolerance)
  {
    std::cerr<<"\n ERROR: Sum error Deposits_h in period " << t << endl;
    Errors << "\n Sum error Deposits_h in period " << t << endl;
  }
  deviation=fabs((Deposits_e(1)-Deposits_eb.Row(1).Sum())/Deposits_eb.Row(1).Sum());
  if(deviation>tolerance && Deposits_eb.Row(1).Sum()>tolerance && Deposits_e(1)>tolerance)
  {
    std::cerr<<"Sum error Deposits_e in period " << t << endl;
    Errors << "\n Sum error Deposits_e in period " << t << endl;
  }
  deviation=fabs((GB_cb(1)+GB_b.Row(1).Sum()-GB(1))/GB(1));
  if(deviation>tolerance)
  {
    std::cerr<<"Sum error GB in period " << t << endl;
    Errors << "\n Sum error GB in period " << t << endl;
  }
  deviation=fabs((Deposits_1.Row(1).Sum()+Deposits_2.Row(1).Sum()+Deposits_hb.Row(1).Sum()+Deposits_eb.Row(1).Sum()-Deposits_b.Row(1).Sum())/Deposits_b.Row(1).Sum());
  if(deviation>tolerance && Deposits_b.Row(1).Sum()>tolerance)
  {
    std::cerr<<"Sum error Deposits in period " << t << endl;
    Errors << "\n Sum error Deposits in period " << t << endl;
  }
  deviation=fabs((Reserves(1)-Reserves_b.Row(1).Sum())/Reserves_b.Row(1).Sum());
  if(deviation>tolerance && Reserves_b.Row(1).Sum()>tolerance && Reserves(1)>tolerance)
  {
    std::cerr<<"Sum error Reserves in period " << t << endl;
    Errors << "\n Sum error Reserves in period " << t << endl;
  }
  deviation=fabs((Advances(1)-Advances_b.Row(1).Sum())/Advances_b.Row(1).Sum());
  if(deviation>tolerance && Advances_b.Row(1).Sum()>tolerance && Advances(1)>tolerance)
  {
    std::cerr<<"Sum error Advances in period " << t << endl;
    Errors << "\n Sum error Advances in period " << t << endl;
  }
  deviation=fabs((Loans_2.Row(1).Sum()+Loans_e(1)-Loans_b.Row(1).Sum())/Loans_b.Row(1).Sum());
  if(deviation>tolerance && Loans_b.Row(1).Sum()>tolerance && Loans_2.Row(1).Sum()>tolerance)
  {
    std::cerr<<"Sum error Loans in period " << t << endl;
    Errors << "\n Sum error Loans in period " << t << endl;
  }

}

void ADJUSTSTOCKS(void)
{
  deviation=0;
  for (i=1; i<=NB; i++)
  {
    prior(i)=Loans_b(1,i)+Reserves_b(1,i)+GB_b(1,i)-Deposits_b(1,i)-Advances_b(1,i);
  }
  prior_cb=Advances(1)+GB_cb(1)-Reserves(1);

  if(Advances(1)<=0 || Advances_b.Row(1).Sum()<=0)
  {
    Advances(1)=0;
    for(i=1; i<=NB; i++)
    {
      Advances_b(1,i)=0;
    }
  }

  if(Reserves(1)<=0 || Reserves_b.Row(1).Sum()<=0)
  {
    Reserves(1)=0;
    for(i=1; i<=NB; i++)
    {
      Reserves_b(1,i)=0;
    }
  }

  if(GB(1)<0)
  {
    GB_cb(1)=GB(1);
    for (i=1; i<=NB; i++)
    {
      GB_b(1,i)=0;
    }
  }

  if(GB(1)>0 && fabs(GB_cb(1)-GB(1))/GB(1)<tolerance && GB_b.Row(1).Sum()<tolerance)
  {
    GB_cb(1)=GB(1);
    for (i=1; i<=NB; i++)
    {
      GB_b(1,i)=0;
    }
  }

  if(GB(1)>0 && fabs(GB_cb(1)-GB(1))/GB(1)>tolerance && GB_b.Row(1).Sum()>tolerance)
  {
    if(GB_b.Row(1).Sum()>0)
    {
      for (i=1; i<=NB; i++)
      {
        ShareBonds(i)=GB_b(1,i)/GB_b.Row(1).Sum();
      }
    }
    else
    {
      for (i=1; i<=NB; i++)
      {
        ShareBonds(i)=(NL_1(i)+NL_2(i))/(N1+N2);
      }
    }
  }

  for(i=1; i<=NB; i++)
  {
    if(Deposits_h(1)>0 && Deposits_hb.Row(1).Sum()>0)
    {
      DepositShare_h(i)=Deposits_hb(1,i)/Deposits_hb.Row(1).Sum();
    }
    else
    {
      DepositShare_h(i)=(NL_2(i)+NL_1(i))/(N1+N2);
    }

    if(Deposits_e(1)>0 && Deposits_eb.Row(1).Sum()>0)
    {
      DepositShare_e(i)=Deposits_eb(1,i)/Deposits_eb.Row(1).Sum();
    }
    else
    {
      DepositShare_e(i)=(NL_2(i)+NL_1(i))/(N1+N2);
    }

    if(Reserves(1)>0)
    {
      ShareReserves(i)=Reserves_b(1,i)/Reserves_b.Row(1).Sum();
    }

    if(Advances(1)>0)
    {
      ShareAdvances(i)=Advances_b(1,i)/Advances_b.Row(1).Sum();
    }
  }

  for (i=1; i<=NB; i++)
  {
    Loans_b(1,i)=0;
    Deposits_b(1,i)=0;

    for(j=1; j<=N2; j++)
    {
      if(BankMatch_2(j,i)==1)
      {
        Loans_b(1,i)+=Loans_2(1,j);
        Deposits_b(1,i)+=Deposits_2(1,j);
      }
    }
    for(j=1; j<=N1; j++)
    {
      if(BankMatch_1(j,i)==1)
      {
        Deposits_b(1,i)+=Deposits_1(1,j);
      }
    }

    if(Deposits_h(1)>0)
    {
      Deposits_hb(1,i)=DepositShare_h(i)*Deposits_h(1);
      Deposits_b(1,i)+=(DepositShare_h(i)*Deposits_h(1));
    }
    else
    {
      Deposits_hb(1,i)=0;
    }

    if(Deposits_e(1)>0)
    {
      Deposits_eb(1,i)=DepositShare_e(i)*Deposits_e(1);
      Deposits_b(1,i)+=(DepositShare_e(i)*Deposits_e(1));
    }
    else
    {
      Deposits_eb(1,i)=0;
    }

    if(Loans_e(1)>0)
    {
      Loans_b(1,i)+=(DepositShare_e(i)*Loans_e(1));
    }

    if(Reserves(1)>0)
    {
      Reserves_b(1,i)=ShareReserves(i)*Reserves(1);
    }
    else
    {
      Reserves_b(1,i)=0;
    }

    if(Advances(1)>0)
    {
      Advances_b(1,i)=ShareAdvances(i)*Advances(1);
    }
    else
    {
      Advances_b(1,i)=0;
    }

    if(GB(1)>0 && fabs(GB_cb(1)-GB(1))/GB(1)>tolerance)
    {
      GB_b(1,i)=ShareBonds(i)*(GB(1)-GB_cb(1));
    }

    post=Loans_b(1,i)+Reserves_b(1,i)+GB_b(1,i)-Deposits_b(1,i)-Advances_b(1,i);
    Adjustment(i)=post-prior(i);
    deviation+=fabs(Adjustment(i));
  }

  post_cb=Advances(1)+GB_cb(1)-Reserves(1);
  Adjustment_cb=post_cb-prior_cb;
  deviation+=fabs(Adjustment_cb);

  if(deviation/GDP_n(1)>tolerance)
  {
    std::cerr << "\n\n ERROR: Adjustment in stocks exceeds tolerance in period " << t << endl;
    Errors << "\n Adjustment in stocks exceeds tolerance in period " << t << endl;
  }
}

void SFC_CHECK(void)
{
  //Calculate the sectoral balances
  Balance_h=Wages+Benefits+InterestDeposits_h+Dividends_h(1)+Bonuses_h+TransferFuel-Taxes_h-Consumption_h-Expenditure_en_h-FirmTransfers+govTransfers;
  Balance_1=Sales1.Sum()+InterestDeposits_1.Sum()+FirmTransfers_1+Transfer_shock_f1.Sum()-Wages_1_i.Sum()-EnergyPayments_1.Sum()-Dividends_1_i.Sum()-Bonuses_1_i.Sum()-Taxes_1.Sum()-Taxes_CO2_1.Sum();
  Balance_2=Sales2.Sum()+InterestDeposits_2.Sum()+FirmTransfers_2+Transfer_shock_f2.Sum()+EnvSubsidies_2-Wages_2_i.Sum()-Investment_2.Sum()-LoanInterest_2.Sum()-EnergyPayments_2.Sum()-Dividends_2_i.Sum()-Bonuses_2_i.Sum()-Taxes_2.Sum()-Taxes_CO2_2.Sum();
  Balance_b=LoanInterest.Sum()+InterestBonds_b.Sum()+InterestReserves_b.Sum()+Bailout_b.Row(1).Sum()+BankTransfer-InterestDeposits.Sum()-Taxes_b.Sum()-InterestAdvances_b.Sum()-Dividends_b_i.Sum()-Bonuses_b_i.Sum();
  Balance_e=EnergyPayments+InterestDeposits_e-Wages_en-Dividends_e-Bonuses_e-Loan_interest_e-Taxes_e-Taxes_CO2_e-FuelCost-Taxes_e_shock-Taxes_e_ff_shock;
  Balance_cb=InterestBonds_cb+InterestAdvances-InterestReserves-TransferCB;
  Balance_g=Taxes_g+TransferCB+Taxes_CO2(1)+Taxes_e_shock+Taxes_f_ff_shock+Taxes_e_ff_shock-InterestBonds-Bailout-EntryCosts-G-govTransfers-Transfer_shock_f-Exp_tot_g-EnvSubsidies_2;
  Balance_f=FuelCost-TransferFuel-Taxes_f_ff_shock;

  //Sectoral balances should sum to zero
  BalanceSum=Balance_h+Balance_1+Balance_2+Balance_b+Balance_e+Balance_cb+Balance_g+Balance_f;
  //Deviation needs to be scaled somehow since model variables (and hence possibly deviations due to rounding) will grow over time
  deviation=fabs(BalanceSum)/(fabs(Balance_h)+fabs(Balance_1)+fabs(Balance_2)+fabs(Balance_b)+fabs(Balance_e)+fabs(Balance_cb)+fabs(Balance_g)+fabs(Balance_f));
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Sectoral balances do not sum to zero in period " << t << endl;
    Errors << "\n Sectoral balances do not sum to zero in period " << t << endl;
  }

  //Household classes checks
    //Check that household classes balances sum to total households balance
  double Balance_mh_tot=0;
  for (const string& cl:classes_mh){
    Balance_mh[cl]=Wages_mh[cl]+Benefits_mh[cl]+InterestDeposits_mh[cl]+Bonuses_mh[cl](1)+Dividends_mh[cl](1)+TransferFuel_mh[cl]-Taxes_mh[cl]-Consumption_mh[cl]-Expenditure_en_mh[cl]-FirmTransfers*Entry_financing_sh_mh[cl]+govTransfers_mh[cl];
    Balance_mh_tot+=Balance_mh[cl];
  }
  deviation=fabs(Balance_mh_tot-Balance_h)/Balance_h;
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Household classes balances do not sum to total household class balance in period " << t << endl;
    Errors << "\n Household classes balances do not sum to total household class balance in period " << t << endl;
  }

  //Compare stock and flow measures of bank net worth
  for(i=1; i<=NB; i++)
  {
    NW_b(1,i)+=Adjustment(i);
    if(NW_b(1,i)<=0 && Bank_active(i)==1)
    {
      std::cerr << "\n\n ERROR: NW of active bank " << i << " is negative in period " << t << endl;
      Errors << "\n NW of active bank " << i << " is negative in period " << t << endl;
    }
    NW_b_c(i)=Loans_b(1,i)+GB_b(1,i)+Reserves_b(1,i)-Deposits_b(1,i)-Advances_b(1,i);
  }
  deviation=fabs((NW_b_c.Sum()-NW_b.Row(1).Sum())/NW_b_c.Sum());
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for BANKS are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for BANKS are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of K-firm net worth
  for(i=1; i<=N1; i++)
  {
    Balances_1(i)=Sales1(i)+InterestDeposits_1(i)+Transfer_shock_f1(i)-Wages_1_i(i)-EnergyPayments_1(i)-Dividends_1_i(i)-Bonuses_1_i(i)-Taxes_1(i)-Taxes_CO2_1(i);
    NW_1(1,i)=Deposits_1(1,i);
    NW_1_c(i)=NW_1(2,i)+Balances_1(i)+baddebt_1(i)+Injection_1(i);
  }
  deviation=fabs((NW_1_c.Sum()-NW_1.Row(1).Sum())/NW_1_c.Sum());
  if(deviation>tolerance && NW_1_c.Sum()>tolerance && NW_1.Row(1).Sum()>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for K-FIRMS are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for K-FIRMS are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of C-firm net worth
  for(i=1; i<=N2; i++)
  {
    NW_2(1,i)=CapitalStock(1,i)+deltaCapitalStock(1,i)+Inventories(1,i)+Deposits_2(1,i)-Loans_2(1,i);
    NW_2_c(i)=NW_2(2,i)+Pi2(1,i)+(t_CO2*Emiss2(i)-Taxes_CO2_2(i))+baddebt_2(i)+Injection_2(i)-Dividends_2_i(i)-Bonuses_2_i(i)-Taxes_2(i)-Loss_Capital(i)-Loss_Inventories(i); //Needed a correction for CO2 tax, which calculated in Pi2 as t_CO2*Emiss2(i)

    deviation=fabs((NW_2_c(i)-NW_2(1,i))/NW_2_c(i));
    if(deviation>tolerance)
    {
      std::cerr << "\n\n ERROR: Stock and flow measures of net worth for C-FIRM " << i << " are not consistent in period " << t << endl;
      Errors << "\n Stock and flow measures of net worth for C-FIRM " << i << " are not consistent in period " << t << endl;
    }
  }

  double dev = fabs(NW_2_c.Sum()-NW_2.Row(1).Sum());
  deviation=fabs(dev/NW_2_c.Sum());
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for C-FIRMS are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for C-FIRMS are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of Energy sector net worth
  NW_e(1)=Deposits_e(1)+CapitalStock_e(1);
  NW_e_c=NW_e(2)+Balance_e+CapitalStock_e(1)-CapitalStock_e(2);
  deviation=fabs((NW_e(1)-NW_e_c)/NW_e_c);
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the ENERGY SECTOR are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for the ENERGY SECTOR are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of Households net worth
    //Update ownership of each household class
  for (const string& cl:classes_mh){
    Ownership_1_i_mh[cl]=SP(Ownership_sh_1_i_mh[cl], NW_1.Row(1)); //For each K-firm
    Ownership_2_i_mh[cl]=SP(Ownership_sh_2_i_mh[cl], NW_2.Row(1)); //For each C-firm
  }
  Ownership_1_mh=umap_sum_each_key(Ownership_1_i_mh);
  Ownership_2_mh=umap_sum_each_key(Ownership_2_i_mh);
    //Calculate total ownership value of each household class
  Ownership_tot_h = 0;
  for (const string& cl:classes_mh){
    Ownership_e_mh[cl] = Ownership_sh_e_mh[cl] * NW_e(1);
    Ownership_b_mh[cl] = Ownership_sh_b_mh[cl] * NW_b.Row(1).Sum();
    Ownership_tot_mh[cl] = Ownership_1_mh[cl] + Ownership_2_mh[cl] + Ownership_e_mh[cl] + Ownership_b_mh[cl];

    Ownership_tot_h += Ownership_tot_mh[cl];
  }
    //Households overall
  NW_h(1)=Deposits_h(1);
  NW_h_c=NW_h(2)+Balance_h+Deposits_recovered_1+Deposits_recovered_2;
  double abs_dev = fabs(NW_h(1)- NW_h_c);
  deviation=fabs(abs_dev/NW_h_c);
  if(deviation>tolerance && NW_h(1)>0 && NW_h_c>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for HOUSEHOLDS are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for HOUSEHOLDS are not consistent in period " << t << endl;
  }
    //Household classes
  for (const string& cl:classes_mh)
  {
    NW_mh[cl](1)=Deposits_mh[cl](1);
    NW_mh_c[cl]=NW_mh[cl](2)+Balance_mh[cl]+Deposits_recovered_1_mh[cl]+Deposits_recovered_2_mh[cl];
    double abs_dev = fabs(NW_mh[cl](1)- NW_mh_c[cl]);
    deviation=fabs(abs_dev/NW_mh_c[cl]); 
    if(deviation>tolerance && NW_mh[cl](1)>0 && NW_mh_c[cl]>tolerance)
    {
      std::cerr << "\n\n ERROR: Stock and flow measures of net worth for HOUSEHOLD CLASS "<< cl <<" are not consistent in period " << t << endl;
      Errors << "\n Stock and flow measures of net worth for HOUSEHOLD CLASS "<< cl <<" are not consistent in period " << t << endl;
    }
  }

  //Compare stock and flow measures of CB net worth
  NW_cb(1)=GB_cb(1)+Advances(1)-Reserves(1)-Deposits_fuel_cb(1);
  NW_cb_c=NW_cb(2)+Balance_cb+Adjustment_cb;

  deviation=fabs((NW_cb(1)-NW_cb_c)/NW_cb_c);
  if(deviation>tolerance && fabs(NW_cb(1)/GDP_n(1))>tolerance && fabs(NW_cb_c/GDP_n(1))>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the CENTRAL BANK are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for the CENTRAL BANK are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of government net worth
  NW_gov(1)=-GB(1);
  NW_gov_c=NW_gov(2)+Balance_g;
  deviation=fabs((NW_gov(1)-NW_gov_c)/NW_gov_c);
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the GOVERNMENT are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for the GOVERNMENT are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of Fossil fuel net worth
  NW_f(1)=Deposits_fuel(1);
  NW_f_c=Deposits_fuel(2)+FuelCost-TransferFuel-Taxes_f_ff_shock;
  deviation=fabs((NW_f(1)-NW_f_c)/NW_f_c);
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the FOSSIL FUEL SECTOR are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for the FOSSIL FUEL SECTOR are not consistent in period " << t << endl;
  }

  //Sum of all sectoral net worths should be equal to nominal value of tangible assets in the economy
  NWSum=NW_h(1)+NW_1.Row(1).Sum()+NW_2.Row(1).Sum()+NW_b.Row(1).Sum()+NW_e(1)+NW_cb(1)+NW_gov(1)+NW_f(1);
  RealAssets=CapitalStock.Row(1).Sum()+deltaCapitalStock.Row(1).Sum()+Inventories.Row(1).Sum()+CapitalStock_e(1);
  deviation=fabs((NWSum-RealAssets)/RealAssets);
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Aggregate net worth not equal to tangible assets in period " << t << endl;
    Errors << "\n Aggregate net worth not equal to tangible assets in period " << t << endl;
  }
}

void OVERBOOST(void)
{
  //Reset t0 to shorten time taken to iterate over technology arrays
  t00=t0;
  flag=0;
  for (tt=t00; tt<=t && flag == 0; tt++)
  {
    for (i=1; i<=N1 && flag == 0; i++)
    {
      for (j=1; j<=N2 && flag == 0; j++)
      {
        if (g[tt-1][i-1][j-1] > 0 || gtemp[tt-1][i-1][j-1] > 0)
          flag=1;
      }
    }
    if (flag == 0)
    {
      t0++;
    }
  }
}

void UPDATE(void)
{
  if (cpi(2)<=0)
	{
		std::cerr << "\n\n ERROR: cpi(t-1)<=0 in period " << t << endl;
    Errors << "\n cpi(t-1)<=0 in period " << t << endl;
    exit(EXIT_FAILURE);
	}
	if (Am(2)<=0)
	{
		std::cerr << "\n\n ERROR: Am(t-1)=0 in period " << t << endl;
    Errors << "\n Am(t-1)=0 in period " << t << endl;
    exit(EXIT_FAILURE);
	}

  //Update mark-up in the energy sector
  if(t>1)
  {
    dw2=kappa*dw2+(1-kappa)*w_tot_for_1_wr_mh(1)/w_tot_for_1_wr_mh(2);
    d_cpi=kappa*d_cpi+(1-kappa)*cpi(1)/cpi(2);
    mi_en*=dw2; //markup energy sector increasing with nominal wage
    CF_ge*=dw2; //Green energy plant cost increasing with nominal wage
    mi_en_preshock*=dw2;
    mi_en_shock*=dw2;
    pf*=dw2;
    pf_preshock*=dw2;
    pf_shock*=dw2;
    c_en_preshock*=dw2;

    if(flag_energyshocks==1 && t>(t_regime_shifts+3))
    {
      mi_en=mi_en+0.175*(mi_en_preshock-mi_en); //Since mi_en_preshock < mi_en, second term is negative and brings fossil fuel price back to value pre-shock
    }

    if(flag_energyshocks==3 && t>(t_regime_shifts+3))
    {
      pf=pf+0.175*(pf_preshock-pf); //Since pf_preshock < pf, second term is negative and brings fossil fuel price back to value pre-shock
    }
  }

    //Check that Ownership is equal to total Net Worth of firms
  double Ownership_1_tot = umap_sum_over_keys(Ownership_1_mh);
  deviation=fabs((Ownership_1_tot-NW_1.Row(1).Sum()))/Ownership_1_tot;
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Total K-firms ownership of household classes does not match their total net worth in period " << t << endl;
    Errors << "\n Total K-firms ownership of household classes does not match their total net worth in period " << t << endl;
  }
  double Ownership_2_tot = umap_sum_over_keys(Ownership_2_mh);
  deviation=fabs((Ownership_2_tot-NW_2.Row(1).Sum()))/Ownership_2_tot;
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Total C-firms ownership of household classes does not match their total net worth in period " << t << endl;
    Errors << "\n Total C-firms ownership of household classes does not match their total net worth in period " << t << endl;
  }
  //Calculate ownerhship shares (for reporting)
  for (const string& cl:classes_mh){
    Ownership_sh_1_mh[cl]= Ownership_1_mh[cl]/Ownership_1_tot;
    Ownership_sh_2_mh[cl]= Ownership_2_mh[cl]/Ownership_2_tot;
  }

  //Update lagged variables needed in next period
  Taxes_CO2(2)=Taxes_CO2(1);
  Deposits_h(2)=Deposits_h(1);
  Dividends_h(2)=Dividends_h(1);
  NW_h(2)=NW_h(1);
  w_tot_for_1_wr_mh(2)=w_tot_for_1_wr_mh(1);
  Income_gross_h(2)=Income_gross_h(1);
  for (const string& cl:classes_mh){
    Deposits_mh[cl](2)=Deposits_mh[cl](1);
    Dividends_mh[cl](2)=Dividends_mh[cl](1);
    Bonuses_mh[cl](2)=Bonuses_mh[cl](1);
    Taxes_div_mh[cl](2)=Taxes_div_mh[cl](1);
    Taxes_bon_mh[cl](2)=Taxes_bon_mh[cl](1);
    NW_mh[cl](2)=NW_mh[cl](1);
    w_mh[cl](2)=w_mh[cl](1);
    Income_gross_mh[cl](2)=Income_gross_mh[cl](1);
  }
  Deposits_e(2)=Deposits_e(1);
  ProfitEnergy(2)=ProfitEnergy(1);
  Deposits_fuel(2)=Deposits_fuel(1);
  Deposits_fuel_cb(2)=Deposits_fuel_cb(1);
  NW_f(2)=NW_f(1);
  GB_cb(2)=GB_cb(1);
  GB(2)=GB(1);
  Advances(2)=Advances(1);
  Reserves(2)=Reserves(1);
  CapitalStock_e(2)=CapitalStock_e(1);
  Loans_e(2)=Loans_e(1);
  NW_e(2)=NW_e(1);
  NW_gov(2)=NW_gov(1);
  NW_cb(2)=NW_cb(1);
  U(2)=U(1);
  Em2(2)=Em2(1);
  ProfitCB(2)=ProfitCB(1);
  c_en(2)=c_en(1);
  GDP_r(2)=GDP_r(1);
  GDP_n(2)=GDP_n(1);
  GDP_r_new(2)=GDP_r_new(1);
  GDP_n_new(2)=GDP_n_new(1);
  p2m(2)=p2m(1);
  cpi(5)=cpi(4);
  cpi(4)=cpi(3);
  cpi(3)=cpi(2);
  cpi(2)=cpi(1);

  cpi_new(5)=cpi_new(4);
  cpi_new(4)=cpi_new(3);
  cpi_new(3)=cpi_new(2);
  cpi_new(2)=cpi_new(1);
  if(t==2) //Save values for climate module
  {
    cpi_init=cpi(1);
    GDP_init=GDP_n(1);
  }
  Am(2)=Am(1);
  Am_en(2)=Am_en(1);

  for (i=1; i<=N1; i++)
  {
    Deposits_1(2,i)=Deposits_1(1,i);
    Pi1(2,i)=Pi1(1,i);
    RD(2,i)=RD(1,i);
    f1(2,i)=f1(1,i);
    NW_1(2,i)=NW_1(1,i);
    S1_temp(2,i)=S1_temp(1,i);
  }

  for (j=1; j<=N2; j++)
  {
    Deposits_2(2,j)=Deposits_2(1,j);
    Pi2(2,i)=Pi2(1,i);
    Loans_2(2,j)=Loans_2(1,j);
    DebtService_2(2,j)=DebtService_2(1,j);
    f2(3,j)=f2(2,j);
    f2(2,j)=f2(1,j);
    D2(2,j)=D2(1,j);
    N(2,j)=N(1,j);
    Inventories(2,j)=Inventories(1,j);
    EI(2,j)=EI(1,j);
    deltaCapitalStock(2,j)=deltaCapitalStock(1,j);
    S2(2,j)=S2(1,j);
    S2_temp(2,j)=S2_temp(1,j);
    mu2(2,j)=mu2(1,j);
    CapitalStock(2,j)=CapitalStock(1,j);
    NW_2(2,j)=NW_2(1,j);

    for (i=1; i<=N1; i++)
    {
      for (tt=t0; tt<=t; tt++)
      {
        if (gtemp[tt-1][i-1][j-1] > 0){
          age[tt-1][i-1][j-1]=age[tt-1][i-1][j-1]+1;
        }
      }
    }
  }

	for (i=1; i<=NB; i++)
  {
    fB(2,i)=fB(1,i);
    Deposits_b(2,i)=Deposits_b(1,i);
    BankProfits(2,i)=BankProfits(1,i);
    Deposits_hb(2,i)=Deposits_hb(1,i);
    Deposits_eb(2,i)=Deposits_eb(1,i);
    GB_b(2,i)=GB_b(1,i);
    Loans_b(2,i)=Loans_b(1,i);
    Advances_b(2,i)=Advances_b(1,i);
    Reserves_b(2,i)=Reserves_b(1,i);
    NW_b(2,i)=NW_b(1,i);
  }
}

///////////WRITE OUTPUT/////////////////////////////

void INITOUTPUTCOLLECTION(void)
{
  if(reducedoutput==1){
    //General
    output_pairs.push_back(std::make_pair("t", [&](){ return t;}));
    output_pairs.push_back(std::make_pair("C_consumption", [&](){ return Consumption_r;}));
    output_pairs.push_back(std::make_pair("unemployment", [&](){ return U(1);}));
    output_pairs.push_back(std::make_pair("public_debt_to_GDP", [&](){ return GB(1)/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("failing_C", [&](){ return n_exit2;}));
    output_pairs.push_back(std::make_pair("failing_banks", [&](){ return counter_bankfailure;}));
    //GDP accounting
    output_pairs.push_back(std::make_pair("GDP_r", [&](){ return GDP_r(1);}));
    output_pairs.push_back(std::make_pair("GDP_r_new", [&](){ return GDP_r_new(1);}));
    //R&D
    output_pairs.push_back(std::make_pair("RandD_num_wr", [&](){ return LD1_rd_mh["wr"];}));   // Determines probability of innovation and imitation in K-sector (thus productivity)
    //Labour
    output_pairs.push_back(std::make_pair("wage_tot_1wr", [&](){ return w_tot_for_1_wr_mh(1)/cpi(1);}));
    //Income inequality
    output_pairs.push_back(std::make_pair("income_net_wr", [&](){ return (Income_gross_mh["wr"](1)-Taxes_mh["wr"]+Taxes_wealth_mh["wr"])/cpi(1);}));
    output_pairs.push_back(std::make_pair("ma_wr_net_income_ratio", [&](){ return ((Income_gross_mh["ma"](1)-Taxes_mh["ma"]+Taxes_wealth_mh["ma"])/LS_sh_mh["ma"])/((Income_gross_mh["wr"](1)-Taxes_mh["wr"]+Taxes_wealth_mh["wr"])/LS_sh_mh["wr"]);})); 

    for (const string& cl:classes_mh){
      //Consumer price indexes
      output_pairs.push_back(std::make_pair("cpi_"+cl, [&](){ return cpi_new_mh[cl];}));
      //Real expenditure properly calculated
      output_pairs.push_back(std::make_pair("real_exp_"+cl, [&](){ return Expenditure_tot_mh[cl]/cpi_new_mh[cl];}));
      //Real expenditure in energy properly calculated
      output_pairs.push_back(std::make_pair("exp_en_"+cl, [&](){ return Expenditure_en_mh[cl]/c_en_h;}));
      //Real consumption properly calculated
       output_pairs.push_back(std::make_pair("C_consumption_"+cl, [&](){ return Consumption_mh[cl]/cpi(1);}));
    }
    //Prices and inflation
    output_pairs.push_back(std::make_pair("CPI_core_abs", [&](){ return  cpi(1);}));
    output_pairs.push_back(std::make_pair("CPI_headline", [&](){ return cpi_new(1);}));
    output_pairs.push_back(std::make_pair("CPI_goods", [&](){ return cpi_goods;})); 
    output_pairs.push_back(std::make_pair("CPI_energy", [&](){ return cpi_energy;}));
    output_pairs.push_back(std::make_pair("energy_cpi_over_C_cpi", [&](){ return cpi_energy/cpi_goods;}));

    output_pairs.push_back(std::make_pair("energy_price", [&](){ return c_en(1)/cpi_new(1);}));
    output_pairs.push_back(std::make_pair("inframarginal_energy_price", [&](){ return c_infra/cpi_new(1);}));
    output_pairs.push_back(std::make_pair("energy_price_over_C_price", [&](){ return c_en(1)/p2m(1);}));
    //Climate
    output_pairs.push_back(std::make_pair("emissions_yearly_EU", [&](){ return Emiss_TOT(1);}));
    
    output_pairs.push_back(std::make_pair("emissions_industry_direct", [&](){ return Emiss1_TOT+Emiss2_TOT;}));
    output_pairs.push_back(std::make_pair("emissions_industry_indirect", [&](){ return Emiss_en * (D1_en_TOT + D2_en_TOT) / D_en_TOT(1);})); 
    output_pairs.push_back(std::make_pair("emissions_residential_indirect", [&](){ return Emiss_en * D_en_h / D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("emissions_government_indirect", [&](){ return Emiss_en * D_en_g / D_en_TOT(1);})); 

    output_pairs.push_back(std::make_pair("en_prod_mean_C", [&](){ return Am_en_2;}));
    output_pairs.push_back(std::make_pair("mean_prod_C", [&](){ return Am2;}));
    output_pairs.push_back(std::make_pair("mean_env_friendl_C", [&](){ return A2_ef_avg;}));
    output_pairs.push_back(std::make_pair("mean_carbon_intensity_C", [&](){ return A2_ci_avg;}));

    //Transition
    output_pairs.push_back(std::make_pair("green_sh_en", [&](){ return Q_ge/D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("green_sh_capacity_en", [&](){ return K_ge/(K_de+K_ge);}));

    //Initialise output stream string with headers
    INITOUTPUTSTRING(output_string, output_pairs);
  } 
  else  //General output
  {
    //General
    output_pairs.push_back(std::make_pair("t", [&](){ return t;}));
    output_pairs.push_back(std::make_pair("C_consumption", [&](){ return Consumption_r;}));
    output_pairs.push_back(std::make_pair("unemployment", [&](){ return U(1);}));
    output_pairs.push_back(std::make_pair("mean_prod_CK", [&](){ return Am(1);}));
    output_pairs.push_back(std::make_pair("mean_prod_C", [&](){ return Am2;}));
    output_pairs.push_back(std::make_pair("mean_prod_K", [&](){ return Am1;}));
    output_pairs.push_back(std::make_pair("mean_prod_CK_growth", [&](){ return log(Am(1))-log(Am(2));}));
    //GDP accounting
    output_pairs.push_back(std::make_pair("GDP_r", [&](){ return GDP_r(1);}));
    output_pairs.push_back(std::make_pair("GDP_growth", [&](){ return GDP_rg;}));
    output_pairs.push_back(std::make_pair("GDP_n", [&](){ return GDP_n(1);}));
    output_pairs.push_back(std::make_pair("GDP_r_new", [&](){ return GDP_r_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_growth_new", [&](){ return GDP_rg_new;}));
    output_pairs.push_back(std::make_pair("GDP_n_new", [&](){ return GDP_n_new(1);}));
      //Expenditure approach
    output_pairs.push_back(std::make_pair("GDP_exp_sh_priv_cons", [&](){ return Expenditure_tot_h/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_exp_sh_pub_cons", [&](){ return Exp_tot_g/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_exp_sh_invest", [&](){ return Investment_n/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_exp_sh_exp", [&](){ return Exports/GDP_n_new(1);}));
      //Production approach
    output_pairs.push_back(std::make_pair("GDP_prod_sh_K", [&](){ return VA_1/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_prod_sh_C", [&](){ return VA_2/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_prod_sh_en", [&](){ return VA_en/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_prod_sh_b", [&](){ return VA_b/GDP_n_new(1);}));
      //Income approach
    output_pairs.push_back(std::make_pair("GDP_inc_sh_labour", [&](){ return PaymentsToLabour/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_inc_sh_gov", [&](){ return PaymentsToGovernment/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("GDP_inc_sh_cap", [&](){ return PaymentsToCapital/GDP_n_new(1);}));
    //C-firms
    output_pairs.push_back(std::make_pair("capacity_utilisation_C", [&](){ return Utilisation;}));
    output_pairs.push_back(std::make_pair("capacity_productive_C", [&](){ return K.Sum();}));
    output_pairs.push_back(std::make_pair("capital_replacement_rate", [&](){ return SI.Sum()/K.Sum();}));
    output_pairs.push_back(std::make_pair("mean_price_C", [&](){ return p2m(1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("mean_markup_C", [&](){ return mu2m;}));
    output_pairs.push_back(std::make_pair("mean_unsatisfied_demand_C", [&](){ return l2m/cpi(1);}));
    output_pairs.push_back(std::make_pair("profits_tot_2", [&](){ return Pitot2/cpi(1);}));
    output_pairs.push_back(std::make_pair("unsold_output_perc_2", [&](){ return Q2temp.Sum()/(Q2.Sum()+N.Row(2).Sum());}));
    output_pairs.push_back(std::make_pair("inventories_2", [&](){ return N.Row(1).Sum();}));
    output_pairs.push_back(std::make_pair("num_2_with_unsold_output", [&](){ return with_unsold_Q_2;}));
    output_pairs.push_back(std::make_pair("num_2_with_negative_profits", [&](){ return with_neg_profits_2;}));    
      //Investment
    output_pairs.push_back(std::make_pair("investment_tot", [&](){ return Investment_r;}));
    output_pairs.push_back(std::make_pair("investment_repl", [&](){ return ReplacementInvestment_r;}));
    output_pairs.push_back(std::make_pair("investment_expan", [&](){ return ExpansionInvestment_r;}));
    output_pairs.push_back(std::make_pair("env_subs_over_invest", [&](){ return EnvSubsidies_2/Investment_n;}));
    ///K-firms
    output_pairs.push_back(std::make_pair("profits_tot_K", [&](){ return Pitot1/cpi(1);}));
    output_pairs.push_back(std::make_pair("sales_tot_K", [&](){ return S1.Sum()/cpi(1);}));
    output_pairs.push_back(std::make_pair("mean_price_K", [&](){ return p1_avg/cpi(1);}));
    output_pairs.push_back(std::make_pair("tot_investment_prod_capacity", [&](){ return I.Sum();}));
    output_pairs.push_back(std::make_pair("mean_price_C_K_ratio", [&](){ return p2m(1)/p1_avg;}));
    //Entry-exit
    output_pairs.push_back(std::make_pair("entry_tot_h", [&](){ return FirmTransfers/cpi(1);}));
    output_pairs.push_back(std::make_pair("entry_tot", [&](){ return (FirmTransfers_1+FirmTransfers_2)/cpi(1);}));
    output_pairs.push_back(std::make_pair("entry_sh_h", [&](){ return FirmTransfers/(FirmTransfers_1+FirmTransfers_2);}));
    output_pairs.push_back(std::make_pair("failing_C", [&](){ return n_exit2;}));
    output_pairs.push_back(std::make_pair("failing_C_payments", [&](){ return exit_payments2.Sum();}));
    output_pairs.push_back(std::make_pair("failing_C_neg_equity", [&](){ return exit_equity2.Sum();}));
    output_pairs.push_back(std::make_pair("failing_C_marketshare", [&](){ return exit_marketshare2.Sum();}));
    //Government
    output_pairs.push_back(std::make_pair("r_bonds", [&](){ return r_bonds;}));
    output_pairs.push_back(std::make_pair("deficit_to_GDP", [&](){ return Deficit/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("public_debt_to_GDP", [&](){ return GB(1)/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("cons_sh_from_g", [&](){ return Consumption_g/Consumption;}));
    output_pairs.push_back(std::make_pair("u_benefit", [&](){ return Benefits/cpi(1);}));
      //Balance
    output_pairs.push_back(std::make_pair("u_benefit_to_GDP", [&](){ return Benefits/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("govtransf_to_GDP", [&](){ return govTransfers/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("transfer_CB_to_GDP", [&](){ return TransferCB/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("interests_bonds_to_GDP", [&](){ return InterestBonds/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("bailouts_to_GDP", [&](){ return Bailout/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("public_cons_to_GDP", [&](){ return Exp_tot_g/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("entry_cost_to_GDP", [&](){ return EntryCosts/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("env_subsidies_to_GDP", [&](){ return EnvSubsidies_2/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_to_GDP", [&](){ return Taxes_tot_g/GDP_n_new(1);}));
        //Tax sources
    output_pairs.push_back(std::make_pair("taxes_h_w_to_GDP", [&](){ return Taxes_w_h/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_h_div_to_GDP", [&](){ return Taxes_div_h/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_h_bon_to_GDP", [&](){ return Taxes_bon_h/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_h_wealth_to_GDP", [&](){ return Taxes_wealth_h/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_C_to_GDP", [&](){ return Taxes_2.Sum()/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_K_to_GDP", [&](){ return Taxes_1.Sum()/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_en_to_GDP", [&](){ return Taxes_e/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_CO2_to_GDP", [&](){ return Taxes_CO2(1)/GDP_n_new(1);}));
    output_pairs.push_back(std::make_pair("taxes_en_shock_to_GDP", [&](){ return (Taxes_f_ff_shock+Taxes_e_ff_shock)/GDP_n_new(1);}));    
      //Taxation
    output_pairs.push_back(std::make_pair("taxes_sh_C02", [&](){ return Taxes_CO2(1)/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("taxes_sh_en_shock", [&](){ return (Taxes_e_shock+Taxes_f_ff_shock+Taxes_e_ff_shock)/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("taxes_sh_h_w", [&](){ return Taxes_w_h/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("taxes_sh_h_div", [&](){ return Taxes_div_h/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("taxes_sh_h_bon", [&](){ return Taxes_bon_h/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("taxes_sh_h_wealth", [&](){ return Taxes_wealth_h/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("taxes_sh_K", [&](){ return Taxes_1.Sum()/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("taxes_sh_C", [&](){ return Taxes_2.Sum()/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("taxes_sh_b", [&](){ return Taxes_b.Sum()/Taxes_tot_g;}));
    output_pairs.push_back(std::make_pair("wage_tax_rate_wr", [&](){ return aliqw_mh["wr"];}));
    output_pairs.push_back(std::make_pair("wage_tax_rate_pr", [&](){ return aliqw_mh["pr"];}));
    output_pairs.push_back(std::make_pair("wage_tax_rate_ma", [&](){ return aliqw_mh["ma"];}));
    output_pairs.push_back(std::make_pair("divid_tax_rate", [&](){ return aliqdiv;}));
    output_pairs.push_back(std::make_pair("wealth_tax_rate_ma", [&](){ return aliqwealth_mh["ma"];}));
    //R&D
    output_pairs.push_back(std::make_pair("RandD_expenditure", [&](){ return (RD_en_de+RD_en_ge+RD.Row(1).Sum())/cpi(1);}));
    output_pairs.push_back(std::make_pair("RandD_expenditure_K", [&](){ return RD.Row(1).Sum()/cpi(1);}));
    output_pairs.push_back(std::make_pair("RandD_sh_K", [&](){ return RD.Row(1).Sum()/(RD_en_de+RD_en_ge+RD.Row(1).Sum());}));
    output_pairs.push_back(std::make_pair("RandD_num_wr", [&](){ return LD1_rd_mh["wr"];}));   // Determines probability of innovation and imitation in K-sector (thus productivity)
    output_pairs.push_back(std::make_pair("num_K_firms_innovating", [&](){ return Inn.Sum();})); 
    output_pairs.push_back(std::make_pair("num_K_firms_imitating", [&](){ return Imm.Sum();}));
    //Labour
    output_pairs.push_back(std::make_pair("wage_wr", [&](){ return w_mh["wr"](1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("wage_tot_1wr", [&](){ return w_tot_for_1_wr_mh(1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("wage_growth", [&](){ return (w_mh["wr"](1)-w_mh["wr"](2))/w_mh["wr"](2);}));
    output_pairs.push_back(std::make_pair("wages_tot_h", [&](){ return Wages/cpi(1);}));
    output_pairs.push_back(std::make_pair("dividends_h", [&](){ return Dividends_h(1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("bonuses_h", [&](){ return Bonuses_h/cpi(1);}));
    output_pairs.push_back(std::make_pair("labour_income_h", [&](){ return (Bonuses_h+Wages)/cpi(1);}));
    output_pairs.push_back(std::make_pair("deposits_h", [&](){ return Deposits_h(1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("labour_sh_C", [&](){ return LD2_wr/LD_mh["wr"];}));
    output_pairs.push_back(std::make_pair("labour_sh_RandD", [&](){ return (LDen_rd_de_mh["wr"]+LDen_rd_ge_mh["wr"]+LD1_rd_mh["wr"])/LD_mh["wr"];}));
    output_pairs.push_back(std::make_pair("labour_sh_en", [&](){ return LDen_tot_mh["wr"]/LD_mh["wr"];}));
    output_pairs.push_back(std::make_pair("wealth_sh_dep_h", [&](){ return Deposits_h(1)/(NW_h(1)+Ownership_tot_h);}));
    output_pairs.push_back(std::make_pair("taxes_h", [&](){ return Taxes_h/cpi(1);}));
    output_pairs.push_back(std::make_pair("firms_entry_transf_h", [&](){ return FirmTransfers/cpi(1);}));
      //Expenditure
    output_pairs.push_back(std::make_pair("expenditure_tot", [&](){ return (Expenditure_tot_h+Exp_tot_g)/cpi(1);}));
    output_pairs.push_back(std::make_pair("expenditure_h", [&](){ return Expenditure_tot_h/cpi(1);}));
    output_pairs.push_back(std::make_pair("expenditure_g", [&](){ return Exp_tot_g/cpi(1);}));
    output_pairs.push_back(std::make_pair("expenditure_C_h", [&](){ return Consumption_h/cpi(1);}));
    output_pairs.push_back(std::make_pair("expenditure_en_h", [&](){ return Expenditure_en_h/cpi(1);}));
    output_pairs.push_back(std::make_pair("expenditure_en_g", [&](){ return Exp_en_g/cpi(1);}));    
        //Expenditure (desired) out of income and wealth
    output_pairs.push_back(std::make_pair("exp_from_li", [&](){ return Exp_li/cpi(1);}));
    output_pairs.push_back(std::make_pair("exp_from_ki", [&](){ return Exp_ki/cpi(1);}));
    output_pairs.push_back(std::make_pair("exp_from_dep", [&](){ return Exp_dep/cpi(1);}));
    output_pairs.push_back(std::make_pair("exp_from_u_benefits", [&](){ return Exp_u_benefit/cpi(1);}));
    output_pairs.push_back(std::make_pair("exp_from_govtransf", [&](){ return Exp_govTransf/cpi(1);}));
    output_pairs.push_back(std::make_pair("exp_sh_from_li", [&](){ return Exp_li/Expenditure_tot_h;}));
    output_pairs.push_back(std::make_pair("exp_sh_from_ki", [&](){ return Exp_ki/Expenditure_tot_h;}));
    output_pairs.push_back(std::make_pair("exp_sh_from_dep", [&](){ return Exp_dep/Expenditure_tot_h;}));
    output_pairs.push_back(std::make_pair("exp_sh_from_u_benefits", [&](){ return Exp_u_benefit/Expenditure_tot_h;}));
    output_pairs.push_back(std::make_pair("exp_sh_from_govtransf", [&](){ return Exp_govTransf/Expenditure_tot_h;}));
      //Income sources
    output_pairs.push_back(std::make_pair("income_gross_h", [&](){ return Income_gross_h(1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("taxes_sh_of_inc_h", [&](){ return Taxes_h/Income_gross_h(1);}));
    output_pairs.push_back(std::make_pair("wage_sh_of_inc_h", [&](){ return Wages/Income_gross_h(1);}));
    output_pairs.push_back(std::make_pair("dole_sh_of_inc_h", [&](){ return Benefits/Income_gross_h(1);}));
    output_pairs.push_back(std::make_pair("interests_sh_of_inc_h", [&](){ return InterestDeposits_h/Income_gross_h(1);}));
    output_pairs.push_back(std::make_pair("divid_sh_of_inc_h", [&](){ return Dividends_h(1)/Income_gross_h(1);}));
    output_pairs.push_back(std::make_pair("bonuses_sh_of_inc_h", [&](){ return Bonuses_h/Income_gross_h(1);}));
        //Dividends sources
    output_pairs.push_back(std::make_pair("divid_sh_C_h", [&](){ return Dividends_2/Dividends_h(1);}));
    output_pairs.push_back(std::make_pair("divid_sh_K_h", [&](){ return Dividends_1/Dividends_h(1);}));
    output_pairs.push_back(std::make_pair("divid_sh_en_h", [&](){ return Dividends_e/Dividends_h(1);}));
    output_pairs.push_back(std::make_pair("divid_sh_b_h", [&](){ return Dividends_b/Dividends_h(1);}));
        //Bonuses sources
    output_pairs.push_back(std::make_pair("bonuses_sh_C_h", [&](){ return Bonuses_2/Bonuses_h;}));
    output_pairs.push_back(std::make_pair("bonuses_sh_K_h", [&](){ return Bonuses_1/Bonuses_h;}));
    output_pairs.push_back(std::make_pair("bonuses_sh_en_h", [&](){ return Bonuses_e/Bonuses_h;}));
    output_pairs.push_back(std::make_pair("bonuses_sh_b_h", [&](){ return Bonuses_b/Bonuses_h;}));
      //Inequality
        //Piketty indicators
    /* output_pairs.push_back(std::make_pair("r_wealth_wr", [&](){ return r_wealth_mh["wr"];}));
    output_pairs.push_back(std::make_pair("r_wealth_pr", [&](){ return r_wealth_mh["pr"];}));
    output_pairs.push_back(std::make_pair("r_wealth_ma", [&](){ return r_wealth_mh["ma"];}));
    output_pairs.push_back(std::make_pair("r_wealth_h", [&](){ return r_wealth_h;}));
    output_pairs.push_back(std::make_pair("income_growth_wr", [&](){ return income_gross_growth_mh["wr"];}));
    output_pairs.push_back(std::make_pair("income_growth_pr", [&](){ return income_gross_growth_mh["pr"];}));
    output_pairs.push_back(std::make_pair("income_growth_ma", [&](){ return income_gross_growth_mh["ma"];}));
    output_pairs.push_back(std::make_pair("income_growth_h", [&](){ return income_gross_growth_h;}));
    output_pairs.push_back(std::make_pair("r_g_wr", [&](){ return r_wealth_mh["wr"]-income_gross_growth_mh["wr"];}));
    output_pairs.push_back(std::make_pair("r_g_pr", [&](){ return r_wealth_mh["pr"]-income_gross_growth_mh["pr"];}));
    output_pairs.push_back(std::make_pair("r_g_ma", [&](){ return r_wealth_mh["ma"]-income_gross_growth_mh["ma"];}));
    output_pairs.push_back(std::make_pair("r_g_h", [&](){ return r_wealth_h-income_gross_growth_h;}));
    output_pairs.push_back(std::make_pair("wealth_income_ratio_wr", [&](){ return (NW_mh["wr"](1)+Ownership_tot_mh["wr"])/Income_gross_mh["wr"](1);}));
    output_pairs.push_back(std::make_pair("wealth_income_ratio_pr", [&](){ return (NW_mh["pr"](1)+Ownership_tot_mh["pr"])/Income_gross_mh["pr"](1);}));
    output_pairs.push_back(std::make_pair("wealth_income_ratio_ma", [&](){ return (NW_mh["ma"](1)+Ownership_tot_mh["ma"])/Income_gross_mh["ma"](1);}));
    output_pairs.push_back(std::make_pair("wealth_income_ratio_h", [&](){ return (NW_h(1)+Ownership_tot_h)/Income_gross_h(1);})); */
        //Income
    output_pairs.push_back(std::make_pair("wage_ratio_pr_wr", [&](){ return w_ratios_mh["pr"];}));
    output_pairs.push_back(std::make_pair("wage_ratio_ma_wr", [&](){ return w_ratios_mh["ma"];}));
    output_pairs.push_back(std::make_pair("wage_sh_wr", [&](){ return Wages_mh["wr"]/Wages;}));
    output_pairs.push_back(std::make_pair("wage_sh_pr", [&](){ return Wages_mh["pr"]/Wages;}));
    output_pairs.push_back(std::make_pair("wage_sh_ma", [&](){ return Wages_mh["ma"]/Wages;}));
    output_pairs.push_back(std::make_pair("lab_inc_gross_sh_wr", [&](){ return (Wages_mh["wr"]+Bonuses_mh["wr"](1))/PaymentsToLabour;}));
    output_pairs.push_back(std::make_pair("lab_inc_gross_sh_pr", [&](){ return (Wages_mh["pr"]+Bonuses_mh["pr"](1))/PaymentsToLabour;}));
    output_pairs.push_back(std::make_pair("lab_inc_gross_sh_ma", [&](){ return (Wages_mh["ma"]+Bonuses_mh["ma"](1))/PaymentsToLabour;}));
    output_pairs.push_back(std::make_pair("income_sh_gross_wr", [&](){ return Income_gross_mh["wr"](1)/Income_gross_h(1);}));
    output_pairs.push_back(std::make_pair("income_sh_gross_pr", [&](){ return Income_gross_mh["pr"](1)/Income_gross_h(1);}));
    output_pairs.push_back(std::make_pair("income_sh_gross_ma", [&](){ return Income_gross_mh["ma"](1)/Income_gross_h(1);}));
    output_pairs.push_back(std::make_pair("income_sh_net_wr", [&](){ return (Income_gross_mh["wr"](1)-Taxes_mh["wr"]+Taxes_wealth_mh["wr"])/(Income_gross_h(1)-Taxes_h+Taxes_wealth_h);}));
    output_pairs.push_back(std::make_pair("income_sh_net_pr", [&](){ return (Income_gross_mh["pr"](1)-Taxes_mh["pr"]+Taxes_wealth_mh["pr"])/(Income_gross_h(1)-Taxes_h+Taxes_wealth_h);}));
    output_pairs.push_back(std::make_pair("income_sh_net_ma", [&](){ return (Income_gross_mh["ma"](1)-Taxes_mh["ma"]+Taxes_wealth_mh["ma"])/(Income_gross_h(1)-Taxes_h+Taxes_wealth_h);}));
    output_pairs.push_back(std::make_pair("income_net_wr", [&](){ return (Income_gross_mh["wr"](1)-Taxes_mh["wr"]+Taxes_wealth_mh["wr"])/cpi(1);}));
    output_pairs.push_back(std::make_pair("income_net_pr", [&](){ return (Income_gross_mh["pr"](1)-Taxes_mh["pr"]+Taxes_wealth_mh["pr"])/cpi(1);}));
    output_pairs.push_back(std::make_pair("income_net_ma", [&](){ return (Income_gross_mh["ma"](1)-Taxes_mh["ma"]+Taxes_wealth_mh["ma"])/cpi(1);}));
    output_pairs.push_back(std::make_pair("ma_wr_gross_income_ratio", [&](){ return (Income_gross_mh["ma"](1)/LS_sh_mh["ma"])/(Income_gross_mh["wr"](1)/LS_sh_mh["wr"]);}));
    output_pairs.push_back(std::make_pair("ma_wr_net_income_ratio", [&](){ return ((Income_gross_mh["ma"](1)-Taxes_mh["ma"]+Taxes_wealth_mh["ma"])/LS_sh_mh["ma"])/((Income_gross_mh["wr"](1)-Taxes_mh["wr"]+Taxes_wealth_mh["wr"])/LS_sh_mh["wr"]);})); 
    output_pairs.push_back(std::make_pair("ma_wr_lab_income_gross_ratio", [&](){ return ((Wages_mh["ma"]+Bonuses_mh["ma"](1))/LS_sh_mh["ma"])/((Wages_mh["wr"]+Bonuses_mh["wr"](1))/LS_sh_mh["wr"]);}));
        //Taxation
    output_pairs.push_back(std::make_pair("taxes_sh_of_inc_wr", [&](){ return Taxes_mh["wr"]/Income_gross_mh["wr"](1);}));
    output_pairs.push_back(std::make_pair("taxes_sh_of_inc_pr", [&](){ return Taxes_mh["pr"]/Income_gross_mh["pr"](1);}));
    output_pairs.push_back(std::make_pair("taxes_sh_of_inc_ma", [&](){ return Taxes_mh["ma"]/Income_gross_mh["ma"](1);}));
          //Carbon tax incidence
    for (const string& cl:classes_mh){
      output_pairs.push_back(std::make_pair("carbon_tax_incidence_"+cl, [&](){ return tax_CO2_incidence_on_income_mh[cl]*100;}));
    }
        //Wealth
    output_pairs.push_back(std::make_pair("deposits_wr", [&](){ return Deposits_mh["wr"](1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("deposits_pr", [&](){ return Deposits_mh["pr"](1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("deposits_ma", [&](){ return Deposits_mh["ma"](1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("deposits_sh_wr", [&](){ return Deposits_sh_mh["wr"];}));
    output_pairs.push_back(std::make_pair("deposits_sh_pr", [&](){ return Deposits_sh_mh["pr"];}));
    output_pairs.push_back(std::make_pair("deposits_sh_ma", [&](){ return Deposits_sh_mh["ma"];}));
    /* output_pairs.push_back(std::make_pair("ownership_K_wr", [&](){ return Ownership_1_mh["wr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("ownership_K_pr", [&](){ return Ownership_1_mh["pr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("ownership_K_ma", [&](){ return Ownership_1_mh["ma"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("ownership_sh_K_wr", [&](){ return Ownership_sh_1_mh["wr"];}));
    output_pairs.push_back(std::make_pair("ownership_sh_K_pr", [&](){ return Ownership_sh_1_mh["pr"];}));
    output_pairs.push_back(std::make_pair("ownership_sh_K_ma", [&](){ return Ownership_sh_1_mh["ma"];}));
    output_pairs.push_back(std::make_pair("ownership_C_wr", [&](){ return Ownership_2_mh["wr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("ownership_C_pr", [&](){ return Ownership_2_mh["pr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("ownership_C_ma", [&](){ return Ownership_2_mh["ma"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("ownership_sh_C_wr", [&](){ return Ownership_sh_1_mh["wr"];}));
    output_pairs.push_back(std::make_pair("ownership_sh_C_pr", [&](){ return Ownership_sh_2_mh["pr"];}));
    output_pairs.push_back(std::make_pair("ownership_sh_C_ma", [&](){ return Ownership_sh_2_mh["ma"];})); */
    output_pairs.push_back(std::make_pair("NW_ratio_wr", [&](){ return (NW_mh["wr"](1)+Ownership_tot_mh["wr"])/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_pr", [&](){ return (NW_mh["pr"](1)+Ownership_tot_mh["pr"])/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_ma", [&](){ return (NW_mh["ma"](1)+Ownership_tot_mh["ma"])/(GDP_n(1)*4);}));
        //Expenditure
          //Consumer price indexes
    for (const string& cl:classes_mh){
      output_pairs.push_back(std::make_pair("cpi_"+cl, [&](){ return cpi_new_mh[cl];}));
    }
          //Expenditure properly calculated
    for (const string& cl:classes_mh){
      output_pairs.push_back(std::make_pair("real_exp_"+cl, [&](){ return Expenditure_tot_mh[cl]/cpi_new_mh[cl];}));
      output_pairs.push_back(std::make_pair("real_exp_2_"+cl, [&](){ return Consumption_mh[cl]/cpi_goods + Expenditure_en_mh[cl]/cpi_energy;}));
    }              
          //Each class share of total expenditure
    output_pairs.push_back(std::make_pair("expenditure_wr", [&](){ return Expenditure_tot_mh["wr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("expenditure_pr", [&](){ return Expenditure_tot_mh["pr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("expenditure_ma", [&](){ return Expenditure_tot_mh["ma"]/cpi(1);}));   
    output_pairs.push_back(std::make_pair("exp_sh_wr", [&](){ return Expenditure_tot_mh["wr"]/Expenditure_tot_h;}));
    output_pairs.push_back(std::make_pair("exp_sh_pr", [&](){ return Expenditure_tot_mh["pr"]/Expenditure_tot_h;}));
    output_pairs.push_back(std::make_pair("exp_sh_ma", [&](){ return Expenditure_tot_mh["ma"]/Expenditure_tot_h;}));
          //Consumption of each household class
    output_pairs.push_back(std::make_pair("C_consumption_wr", [&](){ return Consumption_mh["wr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("C_consumption_pr", [&](){ return Consumption_mh["pr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("C_consumption_ma", [&](){ return Consumption_mh["ma"]/cpi(1);}));;
    output_pairs.push_back(std::make_pair("C_consumption_sh_wr", [&](){ return Cons_sh_mh["wr"];}));
    output_pairs.push_back(std::make_pair("C_consumption_sh_pr", [&](){ return Cons_sh_mh["pr"];}));
    output_pairs.push_back(std::make_pair("C_consumption_sh_ma", [&](){ return Cons_sh_mh["ma"];}));
          //Consumption reduction due to lack of deposits and/or of lack of supply
    output_pairs.push_back(std::make_pair("C_cons_reduction_wr", [&](){ return (Cons_mh["wr"] - Consumption_mh["wr"])/Cons_mh["wr"];}));
    output_pairs.push_back(std::make_pair("C_cons_reduction_pr", [&](){ return (Cons_mh["pr"] - Consumption_mh["pr"])/Cons_mh["pr"];}));
    output_pairs.push_back(std::make_pair("C_cons_reduction_ma", [&](){ return (Cons_mh["ma"] - Consumption_mh["ma"])/Cons_mh["ma"];}));
          //Desired consumption demand based on income and APCs
    output_pairs.push_back(std::make_pair("C_consumption_des_wr", [&](){ return Cons_mh["wr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("C_consumption_des_pr", [&](){ return Cons_mh["pr"]/cpi(1);}));
    output_pairs.push_back(std::make_pair("C_consumption_des_ma", [&](){ return Cons_mh["ma"]/cpi(1);}));
          // Effective average propensities to consume of housheolds overall
    output_pairs.push_back(std::make_pair("APC_eff_net_income_wr", [&](){ return Expenditure_tot_mh["wr"]/(Income_gross_mh["wr"](1)-Taxes_mh["wr"]+Taxes_wealth_mh["wr"]);}));
    output_pairs.push_back(std::make_pair("APC_eff_net_income_pr", [&](){ return Expenditure_tot_mh["pr"]/(Income_gross_mh["pr"](1)-Taxes_mh["pr"]+Taxes_wealth_mh["pr"]);}));
    output_pairs.push_back(std::make_pair("APC_eff_net_income_ma", [&](){ return Expenditure_tot_mh["ma"]/(Income_gross_mh["ma"](1)-Taxes_mh["ma"]+Taxes_wealth_mh["ma"]);}));
    output_pairs.push_back(std::make_pair("APC_eff_net_income_h", [&](){ return Expenditure_tot_h/(Income_gross_h(1)-Taxes_h+Taxes_wealth_h);}));
    output_pairs.push_back(std::make_pair("APC_out_of_li_h", [&](){ return Exp_li/(Wages+Benefits+Bonuses_h-Taxes_h);}));
    output_pairs.push_back(std::make_pair("APC_out_of_ki_h", [&](){ return Exp_ki/(InterestDeposits_h+Dividends_h(1));}));
    output_pairs.push_back(std::make_pair("APC_out_of_dep_h", [&](){ return Exp_dep/Deposits_h(1);}));
          //Expenditure in energy
    output_pairs.push_back(std::make_pair("exp_sh_en_h", [&](){ return effective_en_exp_sh_h;}));
    for (const string& cl:classes_mh){
      output_pairs.push_back(std::make_pair("exp_sh_en_"+cl, [&](){ return effective_en_exp_sh_mh[cl];}));
      output_pairs.push_back(std::make_pair("exp_en_"+cl, [&](){ return Expenditure_en_mh[cl]/c_en_h;}));
      output_pairs.push_back(std::make_pair("exp_en_2_"+cl, [&](){ return Expenditure_en_mh[cl]/cpi(1);}));
    }
        //Energy
    /* output_pairs.push_back(std::make_pair("en_fp_wr", [&](){ return en_fp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("en_fp_pr", [&](){ return en_fp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("en_fp_ma", [&](){ return en_fp_mh["ma"];}));
    output_pairs.push_back(std::make_pair("en_fp_sh_wr", [&](){ return en_fp_mh["wr"]/en_fp_h;}));
    output_pairs.push_back(std::make_pair("en_fp_sh_pr", [&](){ return en_fp_mh["pr"]/en_fp_h;}));
    output_pairs.push_back(std::make_pair("en_fp_sh_ma", [&](){ return en_fp_mh["ma"]/en_fp_h;}));
    output_pairs.push_back(std::make_pair("cons_dir_en_fp_sh_h", [&](){ return cons_dir_en_fp_h/en_fp_h;}));
    output_pairs.push_back(std::make_pair("cons_indir_en_fp_sh_h", [&](){ return cons_indir_en_fp_h/en_fp_h;}));
    output_pairs.push_back(std::make_pair("inv_en_fp_sh_h", [&](){ return invest_en_fp_h/en_fp_h;}));
    output_pairs.push_back(std::make_pair("publ_exp_en_fp_sh_h", [&](){ return publ_cons_en_fp_h/en_fp_h;}));
    output_pairs.push_back(std::make_pair("cons_dir_en_fp_sh_wr", [&](){ return cons_dir_en_fp_mh["wr"]/en_fp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("cons_indir_en_fp_sh_wr", [&](){ return cons_indir_en_fp_mh["wr"]/en_fp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("inv_en_fp_sh_wr", [&](){ return invest_en_fp_mh["wr"]/en_fp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("cons_dir_en_fp_sh_pr", [&](){ return cons_dir_en_fp_mh["pr"]/en_fp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("cons_indir_en_fp_sh_pr", [&](){ return cons_indir_en_fp_mh["pr"]/en_fp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("inv_en_fp_sh_pr", [&](){ return invest_en_fp_mh["pr"]/en_fp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("cons_dir_en_fp_sh_ma", [&](){ return cons_dir_en_fp_mh["ma"]/en_fp_mh["ma"];}));
    output_pairs.push_back(std::make_pair("cons_indir_en_fp_sh_ma", [&](){ return cons_indir_en_fp_mh["ma"]/en_fp_mh["ma"];}));
    output_pairs.push_back(std::make_pair("inv_en_fp_sh_ma", [&](){ return invest_en_fp_mh["ma"]/en_fp_mh["ma"];}));
    output_pairs.push_back(std::make_pair("en_price_over_wage_wr", [&](){ return c_en(1)/w_mh["wr"](1);}));
    output_pairs.push_back(std::make_pair("en_price_over_wage_pr", [&](){ return c_en(1)/w_mh["pr"](1);}));
    output_pairs.push_back(std::make_pair("en_price_over_wage_ma", [&](){ return c_en(1)/w_mh["ma"](1);}));
    output_pairs.push_back(std::make_pair("en_sh_prod_costs_C", [&](){ return EnergyPayments_2.Sum()/(Wages_2_i.Sum()+Investment_2.Sum()+LoanInterest_2.Sum()+EnergyPayments_2.Sum());}));
    output_pairs.push_back(std::make_pair("en_sh_prod_costs_K", [&](){ return EnergyPayments_1.Sum()/(Wages_1_i.Sum()+EnergyPayments_1.Sum());})); */
      // Carbon
    /* output_pairs.push_back(std::make_pair("carb_fp_wr", [&](){ return carb_fp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("carb_fp_pr", [&](){ return carb_fp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("carb_fp_ma", [&](){ return carb_fp_mh["ma"];}));
    output_pairs.push_back(std::make_pair("ind_carb_fp_wr", [&](){ return carb_fp_mh["wr"]/(LS_sh_mh["wr"]*LS);}));
    output_pairs.push_back(std::make_pair("ind_carb_fp_pr", [&](){ return carb_fp_mh["pr"]/(LS_sh_mh["pr"]*LS);}));
    output_pairs.push_back(std::make_pair("ind_carb_fp_ma", [&](){ return carb_fp_mh["ma"]/(LS_sh_mh["ma"]*LS);}));
    output_pairs.push_back(std::make_pair("ind_carb_fp_pr_to_wr", [&](){ return carb_fp_mh["pr"]/LS_sh_mh["pr"]/(carb_fp_mh["wr"]/LS_sh_mh["wr"]);}));
    output_pairs.push_back(std::make_pair("ind_carb_fp_ma_to_wr", [&](){ return carb_fp_mh["ma"]/LS_sh_mh["ma"]/(carb_fp_mh["wr"]/LS_sh_mh["wr"]);}));
    output_pairs.push_back(std::make_pair("carb_fp_sh_wr", [&](){ return carb_fp_mh["wr"]/carb_fp_h;}));
    output_pairs.push_back(std::make_pair("carb_fp_sh_pr", [&](){ return carb_fp_mh["pr"]/carb_fp_h;}));
    output_pairs.push_back(std::make_pair("carb_fp_sh_ma", [&](){ return carb_fp_mh["ma"]/carb_fp_h;}));
    output_pairs.push_back(std::make_pair("cons_dir_carb_fp_sh_h", [&](){ return cons_dir_carb_fp_h/carb_fp_h;}));
    output_pairs.push_back(std::make_pair("cons_indir_carb_fp_sh_h", [&](){ return cons_indir_carb_fp_h/carb_fp_h;}));
    output_pairs.push_back(std::make_pair("inv_carb_fp_sh_h", [&](){ return invest_carb_fp_h/carb_fp_h;}));
    output_pairs.push_back(std::make_pair("publ_exp_carb_fp_sh_h", [&](){ return publ_cons_carb_fp_h/carb_fp_h;}));
    output_pairs.push_back(std::make_pair("cons_dir_carb_fp_sh_wr", [&](){ return cons_dir_carb_fp_mh["wr"]/carb_fp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("cons_indir_carb_fp_sh_wr", [&](){ return cons_indir_carb_fp_mh["wr"]/carb_fp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("inv_carb_fp_sh_wr", [&](){ return invest_carb_fp_mh["wr"]/carb_fp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("cons_dir_carb_fp_sh_pr", [&](){ return cons_dir_carb_fp_mh["pr"]/carb_fp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("cons_indir_carb_fp_sh_pr", [&](){ return cons_indir_carb_fp_mh["pr"]/carb_fp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("inv_carb_fp_sh_pr", [&](){ return invest_carb_fp_mh["pr"]/carb_fp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("cons_dir_carb_fp_sh_ma", [&](){ return cons_dir_carb_fp_mh["ma"]/carb_fp_mh["ma"];}));
    output_pairs.push_back(std::make_pair("cons_indir_carb_fp_sh_ma", [&](){ return cons_indir_carb_fp_mh["ma"]/carb_fp_mh["ma"];}));
    output_pairs.push_back(std::make_pair("inv_carb_fp_sh_ma", [&](){ return invest_carb_fp_mh["ma"]/carb_fp_mh["ma"];})); */
    //Capital stocks
    output_pairs.push_back(std::make_pair("capital_C", [&](){ return CapitalStock.Row(1).Sum()/cpi(1);}));
    output_pairs.push_back(std::make_pair("capital_en", [&](){ return CapitalStock_e(1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("inventories", [&](){ return N.Row(1).Sum();}));
    output_pairs.push_back(std::make_pair("GDP_over_tot_capital", [&](){ return GDP_n(1)/(CapitalStock.Row(1).Sum()+CapitalStock_e(1));}));
    output_pairs.push_back(std::make_pair("GDP_over_C_capital", [&](){ return GDP_n(1)/CapitalStock.Row(1).Sum();}));
    //Checks
      //Stock-flow ratios (net worth over GDP)
    output_pairs.push_back(std::make_pair("NW_ratio_h", [&](){ return NW_h(1)/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_K", [&](){ return NW_1.Row(1).Sum()/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_C", [&](){ return NW_2.Row(1).Sum()/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_en", [&](){ return NW_e(1)/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_fossil", [&](){ return NW_f(1)/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_gov", [&](){ return NW_gov(1)/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_banks", [&](){ return NW_b.Row(1).Sum()/(GDP_n(1)*4);}));
    output_pairs.push_back(std::make_pair("NW_ratio_cb", [&](){ return NW_cb(1)/(GDP_n(1)*4);}));
      //Sectors relative variables
    output_pairs.push_back(std::make_pair("prod_ratio_CK", [&](){ return Am2/Am1;}));
    output_pairs.push_back(std::make_pair("en_prod_ratio_CK", [&](){ return Am_en_2/Am_en_1;}));
    output_pairs.push_back(std::make_pair("price_ratio_CK", [&](){ return p2m(1)/(p1.Sum()/N1r);}));
    output_pairs.push_back(std::make_pair("env_friendl_ratio_CK", [&](){ return A2_ef_avg/A1_ef_avg;}));
    //Central Bank
    output_pairs.push_back(std::make_pair("monetary_policy_rate", [&](){ return r;}));
    output_pairs.push_back(std::make_pair("inflation_target", [&](){ return d_cpi_target_a;}));
    //Prices and inflation
    output_pairs.push_back(std::make_pair("CPI_core_abs", [&](){ return  cpi(1);}));
    output_pairs.push_back(std::make_pair("CPI_headline", [&](){ return cpi_new(1);}));
    output_pairs.push_back(std::make_pair("CPI_goods", [&](){ return cpi_goods;})); 
    output_pairs.push_back(std::make_pair("CPI_energy", [&](){ return cpi_energy;}));
    output_pairs.push_back(std::make_pair("energy_cpi_over_C_cpi", [&](){ return cpi_energy/cpi_goods;}));

    output_pairs.push_back(std::make_pair("mean_unit_cost_L_C", [&](){ return uc_L_avg_2;}));
    output_pairs.push_back(std::make_pair("mean_unit_cost_en_C", [&](){ return uc_en_avg_2;}));
    output_pairs.push_back(std::make_pair("mean_unit_cost_CO2_C", [&](){ return uc_CO2_avg_2;}));
    output_pairs.push_back(std::make_pair("mean_unit_cost_sh_L_C", [&](){ return uc_L_avg_2/uc_avg_2;}));
    output_pairs.push_back(std::make_pair("mean_unit_cost_sh_en_C", [&](){ return uc_en_avg_2/uc_avg_2;}));
    output_pairs.push_back(std::make_pair("mean_unit_cost_sh_CO2_C", [&](){ return uc_CO2_avg_2/uc_avg_2;}));
    output_pairs.push_back(std::make_pair("inflation_rate_core", [&](){ return inflation_a;}));
    output_pairs.push_back(std::make_pair("inflation_rate_headline", [&](){ return inflation_a_new;}));
    //Banks
    output_pairs.push_back(std::make_pair("banks_bad_debt_over_GDP", [&](){ return baddebt_b.Sum()/GDP_n(1);}));
    output_pairs.push_back(std::make_pair("banks_mean_CAP", [&](){ return capitalAdequacyRatio.Sum()/NB;}));
    output_pairs.push_back(std::make_pair("failing_banks", [&](){ return counter_bankfailure;}));
    output_pairs.push_back(std::make_pair("banks_deposits", [&](){ return Deposits_b.Row(1).Sum()/cpi(1);}));
    output_pairs.push_back(std::make_pair("banks_credit_supply", [&](){ return CreditSupply_all/cpi(1);}));
    output_pairs.push_back(std::make_pair("credit_demand", [&](){ return CreditDemand_all/cpi(1);}));
    output_pairs.push_back(std::make_pair("unmet_credit_demand_perc", [&](){ return (CreditDemand_all-CreditSupply_all)/CreditDemand_all;}));
    output_pairs.push_back(std::make_pair("unmet_credit_demand", [&](){ return (CreditDemand_all-CreditSupply_all)/cpi(1);}));    
    output_pairs.push_back(std::make_pair("loans_C", [&](){ return Loans_2.Row(1).Sum()/cpi(1);}));
    output_pairs.push_back(std::make_pair("banks_profits_over_GDP", [&](){ return BankProfits.Sum()/GDP_n(1);}));
    //Energy
      //Productivity
    output_pairs.push_back(std::make_pair("en_prod_mean_CK", [&](){ return Am_en(1);}));
    output_pairs.push_back(std::make_pair("en_prod_mean_C", [&](){ return Am_en_2;}));
    output_pairs.push_back(std::make_pair("en_prod_mean_K", [&](){ return Am_en_1;}));
    output_pairs.push_back(std::make_pair("en_prod_max_C", [&](){ return A1_en_top;}));
    output_pairs.push_back(std::make_pair("en_prod_mean_growth", [&](){ return (Am_en(1)-Am_en(2))/Am_en(2);}));
      //Energy demand
    output_pairs.push_back(std::make_pair("GDP_r_over_en_dem_f_tot", [&](){ return GDP_r(1)/D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("en_dem_f_tot", [&](){ return D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("en_dem_f_sh_h", [&](){ return D_en_h/D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("en_dem_f_sh_C", [&](){ return D2_en_TOT/D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("en_dem_f_sh_K", [&](){ return D1_en_TOT/D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("en_dem_f_sh_g", [&](){ return D_en_g/D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("en_dem_growth_tot", [&](){ return (D_en_TOT(1)-D_en_TOT(2))/D_en_TOT(2);}));
    output_pairs.push_back(std::make_pair("en_dem_f_h", [&](){ return D_en_h;}));
    output_pairs.push_back(std::make_pair("en_dem_f_wr", [&](){ return D_en_mh["wr"];}));
    output_pairs.push_back(std::make_pair("en_dem_f_pr", [&](){ return D_en_mh["pr"];})); 
    output_pairs.push_back(std::make_pair("en_dem_f_ma", [&](){ return D_en_mh["ma"];}));    
      //Energy sector
    output_pairs.push_back(std::make_pair("en_bill_over_GDP", [&](){ return EnergyPayments/GDP_n_new(1);}));
        //Energy market
    output_pairs.push_back(std::make_pair("energy_price", [&](){ return c_en(1)/cpi_new(1);}));
    output_pairs.push_back(std::make_pair("inframarginal_energy_price", [&](){ return c_infra/cpi_new(1);}));
    output_pairs.push_back(std::make_pair("energy_price_over_C_price", [&](){ return c_en(1)/p2m(1);}));
    output_pairs.push_back(std::make_pair("markup_en", [&](){ return mi_en/cpi_new(1);}));
    output_pairs.push_back(std::make_pair("fossil_fuel_price", [&](){ return pf/cpi(1);}));
    output_pairs.push_back(std::make_pair("profits_en", [&](){ return ProfitEnergy(1)/cpi(1);}));
    output_pairs.push_back(std::make_pair("en_price_over_wage_equiv", [&](){ return c_en(1)/w_tot_for_1_wr_mh(1);}));
        //Transition
    output_pairs.push_back(std::make_pair("cost_green_plants_en", [&](){ return CF_ge(t)/cpi(1);}));
    output_pairs.push_back(std::make_pair("green_sh_en", [&](){ return Q_ge/D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("green_sh_capacity_en", [&](){ return K_ge/(K_de+K_ge);}));
        //Energy shock
    output_pairs.push_back(std::make_pair("transf_en_shock_to_h", [&](){ return Transfer_shock/cpi(1);}));
    //Climate
    output_pairs.push_back(std::make_pair("emissions_sh_en", [&](){ return Emiss_en/Emiss_TOT(1);}));
    output_pairs.push_back(std::make_pair("emissions_sh_K", [&](){ return Emiss1_TOT/Emiss_TOT(1);}));
    output_pairs.push_back(std::make_pair("emissions_sh_C", [&](){ return Emiss2_TOT/Emiss_TOT(1);}));
    output_pairs.push_back(std::make_pair("GDP_r_over_yearly_emissions_EU", [&](){ return GDP_r(1)/Emiss_TOT(1);}));
        output_pairs.push_back(std::make_pair("emissions_en", [&](){ return Emiss_en;}));
    output_pairs.push_back(std::make_pair("emissions_K", [&](){ return Emiss1_TOT;}));
    output_pairs.push_back(std::make_pair("emissions_C", [&](){ return Emiss2_TOT;}));

    output_pairs.push_back(std::make_pair("emissions_industry_direct", [&](){ return Emiss1_TOT+Emiss2_TOT;}));
    output_pairs.push_back(std::make_pair("emissions_industry_indirect", [&](){ return Emiss_en * (D1_en_TOT + D2_en_TOT) / D_en_TOT(1);})); 
    output_pairs.push_back(std::make_pair("emissions_residential_indirect", [&](){ return Emiss_en * D_en_h / D_en_TOT(1);}));
    output_pairs.push_back(std::make_pair("emissions_government_indirect", [&](){ return Emiss_en * D_en_g / D_en_TOT(1);})); 

    output_pairs.push_back(std::make_pair("emissions_industry_direct_sh", [&](){ return (Emiss1_TOT+Emiss2_TOT) /Emiss_TOT(1);}));
    output_pairs.push_back(std::make_pair("emissions_industry_indirect_sh", [&](){ return (Emiss_en * (D1_en_TOT + D2_en_TOT) / D_en_TOT(1)) /Emiss_TOT(1);})); 
    output_pairs.push_back(std::make_pair("emissions_residential_indirect_sh", [&](){ return (Emiss_en * D_en_h / D_en_TOT(1)) /Emiss_TOT(1);}));
    output_pairs.push_back(std::make_pair("emissions_government_indirect_sh", [&](){ return (Emiss_en * D_en_g / D_en_TOT(1)) /Emiss_TOT(1);})); 
    
    output_pairs.push_back(std::make_pair("mean_env_friendl_K", [&](){ return A1_ef_avg;}));
    output_pairs.push_back(std::make_pair("mean_env_friendl_C", [&](){ return A2_ef_avg;}));
      //Energy intensities
    /* output_pairs.push_back(std::make_pair("en_int_energy_exp", [&](){ return en_int_energy_exp;}));
    output_pairs.push_back(std::make_pair("en_int_goods_exp", [&](){ return en_int_goods_exp;}));
    output_pairs.push_back(std::make_pair("en_int_energy_goods_exp_ratio", [&](){ return en_int_energy_exp/en_int_goods_exp;}));
    output_pairs.push_back(std::make_pair("en_int_exp_wr", [&](){ return en_int_exp_mh["wr"];}));
    output_pairs.push_back(std::make_pair("en_int_exp_pr", [&](){ return en_int_exp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("en_int_exp_ma", [&](){ return en_int_exp_mh["ma"];}));
    output_pairs.push_back(std::make_pair("en_int_exp_h", [&](){ return en_int_exp_h;}));   
    output_pairs.push_back(std::make_pair("en_int_exp_ratio_wr_to_pr", [&](){ return en_int_exp_mh["wr"] / en_int_exp_mh["pr"];}));
    output_pairs.push_back(std::make_pair("en_int_exp_ratio_wr_to_ma", [&](){ return en_int_exp_mh["wr"] / en_int_exp_mh["ma"];})); */
      //Carbon intensities
    output_pairs.push_back(std::make_pair("mean_carbon_intensity_C", [&](){ return A2_ci_avg;}));
    output_pairs.push_back(std::make_pair("mean_carbon_intensity_K", [&](){ return A1_ci_avg;}));
       //EU
    output_pairs.push_back(std::make_pair("emissions_yearly_EU", [&](){ return Emiss_TOT(1);}));
    output_pairs.push_back(std::make_pair("emissions_growth_EU", [&](){ return log(Emiss_TOT(1))-log(Emiss_TOT(2));}));
      //Global
    output_pairs.push_back(std::make_pair("emissions_yearly_gl", [&](){ return Emiss_yearly_calib(1);}));
    output_pairs.push_back(std::make_pair("emissions_growth_gl", [&](){ return (Emiss_yearly_calib(1)-Emiss_yearly_calib(2))/Emiss_yearly_calib(2);}));
    output_pairs.push_back(std::make_pair("emissions_in_atm_gl", [&](){ return Cat(1);}));
    output_pairs.push_back(std::make_pair( "temperature",[&](){ return Tmixed(1);}));
    //Policies
    output_pairs.push_back(std::make_pair("carbon_tax_en", [&](){ return t_CO2_en/cpi(1);})); 
    output_pairs.push_back(std::make_pair("carbon_tax", [&](){ return t_CO2/cpi(1);}));
    
    //Initialise output stream string with headers
    INITOUTPUTSTRING(output_string, output_pairs);
  }
}

void INITOUTPUTSTRING(std::ostringstream& out_string, std::vector<std::pair<std::string, std::function<double()>>>& out_vector)
{
  for(auto it=out_vector.begin();it!=out_vector.end();it++) {
      out_string << it->first;

      // Check if it's not the last element
      if (it != out_vector.end() - 1) {
        out_string << ",";
    }
  }
  out_string << std::endl;
}

void ADDTVALUESTOOUTPUTSTRING(std::ostringstream& out_string, std::vector<std::pair<std::string, std::function<double()>>>& out_vector)
{
  for(auto it=out_vector.begin();it!=out_vector.end();it++) {
    out_string << it->second();

    // Check if it's not the last element
    if (it != out_vector.end() - 1) {
        out_string << ",";
    }
  }
  out_string << std::endl;
}

void SAVE(void)
/* LEGEND
- Agent for the variable reported after variable name: _C, _K, _gov, _h, _cb, _banks, _en, _fossil
- Other abbreviations: sh - share; NW - net worth; en - energy; eff - efficiency; dem - demand; gl - global; inc - income; divid - dividends
- End of the name: _n for nominal value; _r for real value, nothing specified: nominal (or a ratio)
*/
{
  //Calculate additional variables only for reporting
  for (const string& cl:classes_mh){
        Income_gross_mh[cl](1) = Wages_mh[cl]+Benefits_mh[cl]+InterestDeposits_mh[cl]+Dividends_mh[cl](1)+Bonuses_mh[cl](1)+govTransfers_mh[cl];
        r_wealth_mh[cl] = (InterestDeposits_mh[cl]+Dividends_mh[cl](1)-Taxes_div_mh[cl](1))/(NW_mh[cl](1)+Ownership_tot_mh[cl]);
        income_gross_growth_mh[cl] = log(Income_gross_mh[cl](1)) - log(Income_gross_mh[cl](2));
  }
  Income_gross_h(1) = Wages+Benefits+InterestDeposits_h+Dividends_h(1)+Bonuses_h+govTransfers;
  r_wealth_h = (InterestDeposits_h+Dividends_h(1))/(NW_h(1)+Ownership_tot_h);
  income_gross_growth_h = log(Income_gross_h(1)) - log(Income_gross_h(2));

  Taxes_tot_g = Taxes_g + Taxes_CO2(1)+Taxes_e_shock+Taxes_f_ff_shock+Taxes_e_ff_shock;
  
  uc_L_avg_2 = w_tot_for_1_wr_mh(1)/Am2;
  uc_en_avg_2 = c_en(1)/Am_en_2;
  uc_CO2_avg_2 = t_CO2 * A2_ci_avg;
  uc_avg_2 = uc_L_avg_2 + uc_en_avg_2 + uc_CO2_avg_2;

  en_int_energy_exp = D_en_h/(Expenditure_en_h/cpi(1));
  en_int_goods_exp = (D1_en_TOT+D2_en_TOT)/(Consumption_h/cpi(1)); 
  effective_en_exp_sh_h = Expenditure_en_h/Expenditure_tot_h;

  for (const string& cl:classes_mh){
    effective_en_exp_sh_mh[cl] = Expenditure_en_mh[cl]/Expenditure_tot_mh[cl];
    en_int_exp_mh[cl] = en_int_energy_exp * effective_en_exp_sh_mh[cl] + en_int_goods_exp * (1-effective_en_exp_sh_mh[cl]);
  }
  en_int_exp_h = en_int_energy_exp * effective_en_exp_sh_h + en_int_goods_exp * (1-effective_en_exp_sh_h);

  with_unsold_Q_2 = 0;
  with_neg_profits_2 = 0;
  for(j=1; j<=N2; j++)
  {
    if (Q2temp(j) > 0){
      with_unsold_Q_2 += 1;
    }
    if (Pi2(1,j) < 0){
      with_neg_profits_2 += 1;
    }
  }

    //Carbon tax incidence
  for (const string& cl:classes_mh){
    tax_CO2_incidence_C_mh[cl] =  Taxes_CO2_2.Sum() * Consumption_h/Consumption * Cons_sh_mh[cl];
    tax_CO2_incidence_en_mh[cl] = Taxes_CO2_e * ((D_en_h + D2_en_TOT)/D_en_TOT(1)) * effective_en_exp_sh_mh[cl];

    tax_CO2_incidence_on_income_mh[cl] = (tax_CO2_incidence_C_mh[cl] + tax_CO2_incidence_en_mh[cl]) / Income_gross_mh[cl](1);
  }

  // Relative energy price index
  cpi_energy = c_en(1) / mi_en0;

  // Relative goods price index
  if (t==1 | t==2){
    cpi_goods = 1;
  }else{
    cpi_goods = cpi(1) / cpi_init;
  }

  //Consumer price index for different households classes
  for (const string& cl:classes_mh){
    cpi_new_mh[cl] = (cpi_goods * Consumption_mh[cl] + cpi_energy * Expenditure_en_mh[cl]) / Expenditure_tot_mh[cl];
  }

  //Consumer price index and inflation considering energy (headline) 
  cpi_new(1) = (cpi_goods * Consumption_h + cpi_energy * Expenditure_en_h) / Expenditure_tot_h;
  inflation_a_new=cpi_new(1)/cpi_new(5)-1;

  ADDTVALUESTOOUTPUTSTRING(output_string, output_pairs);

  //Output file is written at the end of all timesteps
}

/// @brief Write in output files header and variable values (called in SAVE())
/// @param file_name Name of the output file
/// @param out_string Output stream string with values to be printed
void WRITEFILE(char *file_name, std::ostringstream& out_string)
{
    std::ofstream inv;

    inv.open(file_name, ofstream::trunc);
    inv.setf(ios::fixed);
    inv.precision(4);
    inv.setf(ios::right);

    //Print variables
    inv << out_string.str();

    //Close file
    inv << endl;
    inv.close();


}

///////////GENERATE OUTPUT FOLDERS, FILES & NAMES/////////////////////
//These functions generate the directories for saving output and the names of the .txt files in which model output is saved
int make_directory(const char* name)
{
  #ifdef __linux__
    return mkdir(name,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  #elif __APPLE__
    return mkdir(name,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  #else
    return mkdir(name);
  #endif
}

void FOLDERS(char *path, char* name_run)
{
  //Create a folder called "output" in the same directory as the executable
  //Also create a subdirectory of "output" called "errors"
  std::string outstr(path);

  for(j=outstr.length(); j>0; j--) //Remove executable name
  {
    if(outstr[j-1]=='/')
    {
      break;
    }
    else
    {
      outstr.pop_back();
    }
  }
  outstr.pop_back(); //Remove "/" between build fodler and executable names
    //When running (not when debugging) it adds another / in front of the executable name that we need to remove before "build"
  if (outstr[outstr.length()-1]=='/'){
    outstr.pop_back();
  }
  for(j=outstr.length(); j>0; j--) //Remove name of build folder
  {
    if(outstr[j-1]=='/')
    {
      break;
    }
    else
    {
      outstr.pop_back();
    }
  }
  outstr+="output";
  std::string gen_outstr = outstr;
  gen_outstr+="/";
  gen_outstr+=name_run;
  std::string errstr=gen_outstr + "/errors";
  char gen_out_dir[gen_outstr.length()+1];
  char out_dir[outstr.length()+1];
  char err_dir[errstr.length()+1];
  strcpy(gen_out_dir,gen_outstr.c_str());
  strcpy(out_dir,outstr.c_str());
  strcpy(err_dir,errstr.c_str());
	const int out_fol=make_directory(out_dir);
  const int gen_out=make_directory(gen_out_dir);
  const int err_fol=make_directory(err_dir);
}

void GENFILEOUTPUT1(char *path, const char *s1, char runname[], char const* seednumber)
{
	strcpy(nomefile1,path);
	char* name1=strcat(nomefile1,s1);
  name1=strcat(nomefile1,"_");
	name1=strcat(nomefile1,runname);
	name1=strcat(nomefile1,"_");
	name1=strcat(nomefile1,seednumber);
	strcat(nomefile1,".csv");
}

void GENFILEYMC(char *path, const char *s2,char runname[], char const* seednumber)
{
  //Standard file with selected aggregate variables
	strcpy(nomefile2,path);
	char* name2=strcat(nomefile2,s2);
  name2=strcat(nomefile2,"_");
	name2=strcat(nomefile2,runname);
	name2=strcat(nomefile2,"_");
	name2=strcat(nomefile2,seednumber);
	strcat(nomefile2,".csv");
}

void GENFILESHOCKEXP(char *path, const char *s3, char runname[], char const* seednumber)
{
	//File with aggregate variables to save when running climate shock experiment
  strcpy(nomefile3,path);
  char* name3=strcat(nomefile3,s3);
  name3=strcat(nomefile3,"_");
	name3=strcat(nomefile3,runname);
	name3=strcat(nomefile3,"_");
	name3=strcat(nomefile3,seednumber);
	strcat(nomefile3,".txt");
}

void GENFILEPROD1(char *path, const char *s4, char runname[], char const* seednumber)
{
	//File to save log deviation of K-Firm productivity from mean
  strcpy(nomefile4,path);
	char* name5=strcat(nomefile4,s4);
  name5=strcat(nomefile4,"_");
	name5=strcat(nomefile4,runname);
	name5=strcat(nomefile4,"_");
	name5=strcat(nomefile4,seednumber);
	strcat(nomefile4,".csv");
}

void GENFILEPROD2(char *path, const char *s5, char runname[], char const* seednumber)
{
	//File to save log deviation of C-Firm productivity from mean
  strcpy(nomefile5,path);
	char* name4=strcat(nomefile5,s5);
  name4=strcat(nomefile5,"_");
	name4=strcat(nomefile5,runname);
	name4=strcat(nomefile5,"_");
	name4=strcat(nomefile5,seednumber);
	strcat(nomefile5,".txt");
}

void GENFILEPRODALL1(char *path, const char *s6, char runname[], char const* seednumber)
{
	//File to save untransformed K-Firm productivity
  strcpy(nomefile6,path);
	char* name6=strcat(nomefile6,s6);
  name6=strcat(nomefile6,"_");
	name6=strcat(nomefile6,runname);
	name6=strcat(nomefile6,"_");
	name6=strcat(nomefile6,seednumber);
	strcat(nomefile6,".csv");
}

void GENFILEPRODALL2(char *path, const char *s7, char runname[], char const* seednumber)
{
	//File to save untransformed C-Firm productivity
  strcpy(nomefile7,path);
	char* name7=strcat(nomefile7,s7);
  name7=strcat(nomefile7,"_");
	name7=strcat(nomefile7,runname);
	name7=strcat(nomefile7,"_");
	name7=strcat(nomefile7,seednumber);
	strcat(nomefile7,".csv");
}

void GENFILEPRODALL1_en(char *path, const char *s8, char runname[], char const* seednumber)
{
  //File to save energy efficiency of K-Firms
  strcpy(nomefile8,path);
  const char* name8=strcat(nomefile8,s8);
  name8=strcat(nomefile8,"_");
  name8=strcat(nomefile8,runname);
  name8=strcat(nomefile8,"_");
  name8=strcat(nomefile8,seednumber);
  strcat(nomefile8,".csv");
}

void GENFILEPRODALL2_en(char *path, const char *s9, char runname[], char const* seednumber)
{
  //File to save energy efficiency of C-Firms
  strcpy(nomefile9,path);
  const char* name9=strcat(nomefile9,s9);
  name9=strcat(nomefile9,"_");
  name9=strcat(nomefile9,runname);
  name9=strcat(nomefile9,"_");
  name9=strcat(nomefile9,seednumber);
  strcat(nomefile9,".csv");
}


void GENFILEPRODALL1_ef(char *path, const char *s10, char runname[], char const* seednumber)
{
  //File to save environmental friendliness of K-Firms
  strcpy(nomefile10,path);
  const char* name10=strcat(nomefile10,s10);
  name10=strcat(nomefile10,"_");
  name10=strcat(nomefile10,runname);
  name10=strcat(nomefile10,"_");
  name10=strcat(nomefile10,seednumber);
  strcat(nomefile10,".csv");
}

void GENFILEPRODALL2_ef(char *path, const char *s11, char runname[], char const* seednumber)
{
  //File to save environmental friendliness of C-Firms
  strcpy(nomefile11,path);
  const char* name11=strcat(nomefile11,s11);
  name11=strcat(nomefile11,"_");
  name11=strcat(nomefile11,runname);
  name11=strcat(nomefile11,"_");
  name11=strcat(nomefile11,seednumber);
  strcat(nomefile11,".csv");
}

void GENFILENWALL1(char *path, const char *s12, char runname[], char const* seednumber)
{
	strcpy(nomefile12,path);
	char*name12=strcat(nomefile12,s12);
  name12=strcat(nomefile12,"_");
  name12=strcat(nomefile12,runname);
  name12=strcat(nomefile12,"_");
  name12=strcat(nomefile12,seednumber);
  strcat(nomefile12,".csv");
}

void GENFILENWALL2(char *path, const char *s13, char runname[], char const* seednumber)
{
	strcpy(nomefile13,path);
	char*name13=strcat(nomefile13,s13);
  name13=strcat(nomefile13,"_");
  name13=strcat(nomefile13,runname);
  name13=strcat(nomefile13,"_");
  name13=strcat(nomefile13,seednumber);
  strcat(nomefile13,".csv");
}

void GENFILENWALL3(char *path, const char *s14, char runname[], char const* seednumber)
{
	strcpy(nomefile14,path);
	char*name14=strcat(nomefile14,s14);
  name14=strcat(nomefile14,"_");
  name14=strcat(nomefile14,runname);
  name14=strcat(nomefile14,"_");
  name14=strcat(nomefile14,seednumber);
  strcat(nomefile14,".csv");
}

void GENFILEDEBALL2(char *path, const char *s15, char runname[], char const* seednumber)
{
	strcpy(nomefile15,path);
	char*name15=strcat(nomefile15,s15);
  name15=strcat(nomefile15,"_");
  name15=strcat(nomefile15,runname);
  name15=strcat(nomefile15,"_");
  name15=strcat(nomefile15,seednumber);
  strcat(nomefile15,".csv");
}

void GENFILEVALIDATION1(char *path, const char *s16, char const* seednumber)
{
	//File to save aggregate variables needed for validation
  strcpy(nomefile16,path);
	char* name16=strcat(nomefile16,s16);
  name16=strcat(nomefile16,"_");
  name16=strcat(nomefile16,seednumber);
  strcat(nomefile16,".csv");
}

void GENFILEVALIDATION2(char *path, const char *s17, char const* seednumber)
{
	//File to save C-Firms' sales for validation
  strcpy(nomefile17,path);
	char* name17=strcat(nomefile17,s17);
  name17=strcat(nomefile17,"_");
  name17=strcat(nomefile17,seednumber);
  strcat(nomefile17,".csv");
}

void GENFILEVALIDATION3(char *path, const char *s18, char const* seednumber)
{
	//File to save K-Firms' sales for validation
  strcpy(nomefile18,path);
	char* name18=strcat(nomefile18,s18);
  name18=strcat(nomefile18,"_");
  name18=strcat(nomefile18,seednumber);
  strcat(nomefile18,".csv");
}

void GENFILEVALIDATION4(char *path, const char *s19, char const* seednumber)
{
	//File to save log deviation of C-Firms' productivity from mean for validation
  strcpy(nomefile19,path);
	char* name19=strcat(nomefile19,s19);
  name19=strcat(nomefile19,"_");
  name19=strcat(nomefile19,seednumber);
  strcat(nomefile19,".csv");
}

void GENFILEVALIDATION5(char *path, const char *s20, char const* seednumber)
{
	//File to save log deviation of K-Firms' productivity from mean for validation
  strcpy(nomefile20,path);
	char* name20=strcat(nomefile20,s20);
  name20=strcat(nomefile20,"_");
  name20=strcat(nomefile20,seednumber);
  strcat(nomefile20,".csv");
}

void GENFILEVALIDATION6(char *path, const char *s21, char const* seednumber)
{
	//File to save log deviation of C-Firms' energy efficiency from mean for validation
  strcpy(nomefile21,path);
	char* name21=strcat(nomefile21,s21);
  name21=strcat(nomefile21,"_");
  name21=strcat(nomefile21,seednumber);
  strcat(nomefile21,".csv");
}

void GENFILEVALIDATION7(char *path, const char *s22, char const* seednumber)
{
	//File to save log deviation of C-Firms' environmental friendliness from mean for validation
  strcpy(nomefile22,path);
	char* name22=strcat(nomefile22,s22);
  name22=strcat(nomefile22,"_");
  name22=strcat(nomefile22,seednumber);
  strcat(nomefile22,".csv");
}

void GENFILEVALIDATION8(char *path, const char *s23, char const* seednumber)
{
	//File to save log deviation of K-Firms' energy efficiency from mean for validation
  strcpy(nomefile23,path);
	char* name23=strcat(nomefile23,s23);
  name23=strcat(nomefile23,"_");
  name23=strcat(nomefile23,seednumber);
  strcat(nomefile23,".csv");
}

void GENFILEVALIDATION9(char *path, const char *s24, char const* seednumber)
{
	//File to save log deviation of K-Firms' environmental friendliness from mean for validation
  strcpy(nomefile24,path);
	char* name24=strcat(nomefile24,s24);
  name24=strcat(nomefile24,"_");
  name24=strcat(nomefile24,seednumber);
  strcat(nomefile24,".csv");
}

void GENFILEVALIDATION10(char *path, const char *s25, char const* seednumber)
{
	//File to save C-Firms' investment for validation
  strcpy(nomefile25,path);
	char* name25=strcat(nomefile25,s25);
  name25=strcat(nomefile25,"_");
  name25=strcat(nomefile25,seednumber);
  strcat(nomefile25,".csv");
}

void GENFILEVALIDATION11(char *path, const char *s26, char const* seednumber)
{
	//File to save growth rate of individual C-Firm sales for validation
  strcpy(nomefile26,path);
	char* name26=strcat(nomefile26,s26);
  name26=strcat(nomefile26,"_");
  name26=strcat(nomefile26,seednumber);
  strcat(nomefile26,".csv");
}

void GENFILEVALIDATION12(char *path, const char *s27, char const* seednumber)
{
  //File to save growth rate of individual K-Firm sales for validation
	strcpy(nomefile27,path);
	char* name27=strcat(nomefile27,s27);
  name27=strcat(nomefile27,"_");
  name27=strcat(nomefile27,seednumber);
  strcat(nomefile27,".csv");
}

void GENFILESHOCKPARS(char *path, const char *s28,char runname[], char const* seednumber)
{
	strcpy(nomefile28,path);
  char* name3=strcat(nomefile28,s28);
  name3=strcat(nomefile28,"_");
	name3=strcat(nomefile28,runname);
	name3=strcat(nomefile28,"_");
	name3=strcat(nomefile28,seednumber);
	strcat(nomefile28,".txt");
}

//This function generates the actual output files
void INTFILE(void)
{
  ofstream Errors(errorfilename);
  //Standard aggregate output file
  ofstream inv_ymc(nomefile2);
}

///////////AUXILIARY/////////////////////

void catchAlarm(int sig) {
    //If the time taken to perform the run exceeds the threshold set above, abort the simulation

    std::cerr << "\n\n Run timed out!" << endl;
    Errors << "\n Run timed out! " << endl;
    exit(EXIT_FAILURE);
}

double ROUND(double x)
{
  //Rounds a double to the closest integer
  double x_floor=floor(x);
  double resto=x-x_floor;
  if (resto > 0.5) x=x_floor+1;
  else x=x_floor;
  return x;
}

void ALLOCATEBANKCUSTOMERS(void)
{
  //Used during initialisation to assign C-Firms and K-Firms to banks
  //Initialise number of C-firm customers of each bank to 0
  NL_2=0;
  double sum_NL_2;
  sum_NL_2=0;

  //re-perform the random drawing of customer numbers until the sum of the random values drawn is equal to the number of C-Firms
  while (sum_NL_2!=N2){
       sum_NL_2=0;
       for (i=1; i<=NB; i++)
        {
          //Draw from truncated pareto. pareto_a=shape parameter, pareto_k=lower bound, pareto_p=upper bound
          pareto_rv = bpareto(pareto_a, pareto_k, pareto_p);
          NL_2(i)=pareto_rv;
          sum_NL_2+=pareto_rv;
        }
  }

  //Initialise number of C-firm customers of each bank to 0
  NL_1=0;
  double sum_NL_1;
  sum_NL_1=0;
  //re-perform the random drawing of customer numbers until the sum of the random values drawn is equal to the number of K-Firms
  while (sum_NL_1!=N1){
       sum_NL_1=0;
       for (i=1; i<=NB; i++)
        {
          //Draw from truncated pareto. pareto_a=shape parameter, pareto_k=lower bound, pareto_p=upper bound
          //Adjust upper and lower bounds to reflect the number of K-Firms relative to C-Firms
          pareto_rv = bpareto(pareto_a, min(pareto_k*N1r/N2r,0.9), ceil(pareto_p*N1r/N2r));
          NL_1(i)=pareto_rv;
          sum_NL_1+=pareto_rv;
        }
  }
}

double bpareto(double par_a, double par_k, double par_p)
{

  double z;     // Uniform random number from 0 to 1
  double rv;    // RV to be returned

  // Pull a uniform RV (0 < z < 1)
  do
  {
    z=double(ran1(p_seed));
  }
  while ((z == 0) || (z == 1));

  // Generate the bounded Pareto rv using the inversion method
  rv = pow((pow(par_k, par_a) / (z*pow((par_k/par_p), par_a) - z + 1)), (1.0/par_a));
  // make the variable an integer
  rv=ceil(rv);

  return(rv);
}