#ifndef PARS_H
#define PARS_H

int         N1;			                    // Number of K-Firms
int         N2; 			                // Number of C-firms
int         NB;                             // Numbe of Banks
int         T;			                    // Simulation periods

double      varphi; 		                // Maximum bonds to Loans ratio of banks
double      nu;                             // Fraction of revenue devoted to R&D
double      xi;                             // Divides R&D expenditure between innovation and imitation
double      o1;			                    // Effectiveness of R&D expenditure on innovation
double      o2;			                    // Effectiveness of R&D expenditure on imitation
double      uu11;			                // In DSK Case: Lower bound on machine productivity changes due to R&D; in K+S case: Suppport of uniform distribution for exogenous process innovation
double      uu21;			                // In DSK Case: Upper bound on machine productivity changes due to R&D; in K+S case: Suppport of uniform distribution for exogenous product innovation
double      uu12;			                // In DSK Case: Lower bound on own productivity changes due to R&D
double      uu22;			                // In DSK Case: Upper bound on own productivity changes due to R&D
double      uu31;			                // Lower bound on energy efficiency changes due to R&D; C-firms
double      uu41;			                // Upper bound on energy efficiency changes due to R&D; C-firms
double      uu32;			                // Lower bound on energy efficiency changes due to R&D; K-firms
double      uu42;			                // Upper bound on energy efficiency changes due to R&D; K-firms
double      uu51; 		                    // Lower bound on environmental friendliness changes due to R&D; C-firms
double      uu61; 		                    // Upper bound on environmental friendliness changes due to R&D; C-firms
double      uu52; 		                    // Lower bound on environmental friendliness changes due to R&D; K-firms
double      uu62; 		                    // Upper bound on environmental friendliness changes due to R&D; K-firms
double      uinf;		                    // Lower bound on exogenous changes in technological frontier
double      usup;		                    // Upper bound on exogenous changes in technological frontier
double      b_a11;                          // Parameter alpha for Beta distribution governing machine productivity innovation
double      b_a12;                          // Parameter alpha for Beta distribution governing own productivity innovation
double      b_a1_shock;                     // Parameter alpha incorporating climate shock to R&D
double      b_b11;                          // Parameter beta for Beta distribution governing machine productivity innovation
double      b_b12;                          // Parameter beta for Beta distribution governing own productivity innovation
double      b_a2;		                    // Parameter alpha for Beta distribution governing energy efficiency innovation
double      b_a2_shock;                     // Parameter alpha incorporating climate shock to R&D
double      b_b2;		                    // Parameter beta for Beta distribution governing energy efficiency innovation
double      b_a3;		                    // in DSK: Parameter alpha for Beta distribution governing environmental friendliness innovation; in K+S: Parameter alpha for Beta distribution governing change in technological frontier
double      b_a3_shock;                     // Parameter alpha incorporating climate shock to R&D
double      b_b3;		                    // in DSK: Parameter beta for Beta distribution governing environmental friendliness innovation; in K+S: Parameter beta for Beta distribution governing change in technological frontier
double      mi1;                            // Mark-up K-Firms
double      mi2;                            // Initial Mark-up C-firms
double      Gamma;                          // Determining number of potential new clients contacted by K-firms
double      chi;    	                    // Governing replicator dynamics of C-firm market share
double      omega1;		                    // Weight of relative price in C-firm competitiveness
double      omega2;		                    // Weight of unsatisfied demand in C-firm competitiveness
double      omega4;		                    // Weight of environmental friendliness in C-firm competitiveness
double      psi1;                           // Wage sensitivity to inflation
double      psi2;                           // Wage sensitivity to productivity
double      psi3;                           // Wage sensitivity to unemployment
double      deltami2;                       // Sensitivity of C-firm mark-up to change in market share
double      w_min;	                        // Subsistence wage
double      pmin;		                    // lower bound on prices
double      theta;		                    // Probability of being able to re-set price (K-firms and C-firms)
double      u;		                        // C-firms' desired capacity utilisation
double      alfa;                           // Persistence of C-firms' adaptive expectations
double      b;                              // Pay-back period for new technologies
double      dim_mach;                       // Output producible by 1 machine
double      agemax;                         // Maximum lifespan of machine tools
double      credit_multiplier;              // Governing banks' credit supply
double      beta_basel;                     // Sensitivity of credit supply to defaults
double      floor_default_probability;      // Floor level of the probabiity of default of credit demanders (according to the internal risk-based approach of Basel Accords)
double      upsilon;                        // Sensitivity of default probability to leverage ratio
double      lambdaB1;                       // Weight of the credit risk in the interest rate
double      lambdaB2;                       // Weight of the capital adequacy ratio gap in the interest rate
double      riskWeightLoans;                // Risk weight associated to loans with no colateral	
double      riskWeightGovBonds;             // Risk weight associated to government bonds	
double      capitalAdequacyRatioTarget;     // Regulatory Capital Adequacy Ratio target
double      bankmarkdown;                   // Banks' markdown on deposit interest rate
double      centralbankmarkdown;            // CB's markdown on reserve interest rate
double      d1;                             // K-firm dividend payout rate
double      d2;                             // C-firm dividend payout rate
double      db;                             // Bank dividend payout rate
double      repayment_share;                // Loan repayment share
double      bonds_share;                    // Gov. bond repayment share
double      pareto_a;                       // Alpha parameter for pareto distribution function
double      pareto_k;                       // K parameter for pareto distribution function
double      pareto_p;                       // P parameter for pareto distribution function
double      d_cpi_target;                   // Inflation target
double      ustar;		                    // Target unemployment
double      w1sup;                  		// Upper bound on distribution governing deposits transferred to newly entering K-firms
double      w1inf;                  		// Lower bound on distribution governing deposits transferred to newly entering K-firms
double      w2sup;                          // Upper bound on distribution governing deposits transferred to newly entering C-firms
double      w2inf;                          // Lower bound on distribution governing deposits transferred to newly entering C-firms
double      k_const;                        // Governing interest rates charged to individual bank borrowers 
double      k_const2;                       // Governing interest rates charged to individual bank borrowers (additive markup case)
unordered_map<string, double> aliqw_mh;     // Tax rate on wages
unordered_map<string, double> aliqw_mh_prechange;        // Tax rate on wages before tax change
unordered_map<string, double> aliqwealth_mh;             // Tax rate on wealth
unordered_map<string, double> aliqwealth_mh_prechange;   // Tax rate on wealth before tax change
double      aliqw_progressivity_regime_shift;            // Parameter determining the progressivity of labour income taxation when flag_change_wage_tax_rate = 4
double      aliq_households_increase;       // Tax rate increase when flag_highly_progressive_taxation > 0
double      aliqdiv;                        // Tax rate on households dividends
double      aliqdiv_prechange;              // Tax rate on households dividendsr before tax change
double      taylor1;                        // Taylor rule inflation sensitivity
double      taylor2;                        // Taylor rule unemployment sensitivity
double      bondsmarkdown;                  // Markdown on gov. bond interest rate
double      mdw;                            // Maximum variation in factors governing changes in wage
double      phi2;                           // Governing maximum amount of loans demanded by C-firms           
double      b1sup;                          // Upper bound on distribution of net worth of bailed-out banks
double      b1inf;                          // Lower bound on distribution of net worth of bailed-out banks
double      b2sup;                          // Upper bound on distribution of net worth of bailed-out banks (for case when all banks have failed)
double      b2inf;                          // Lower bound on distribution of net worth of bailed-out banks (for case when all banks have failed)
double      aliq;                           // Tax rate on firm profits
double      aliqb;                          // Tax rate on bank profits
double      wu;                             // Unemployment benefit rate (as a share of wage)
double      r_base;                         // Taylor rule intercept
double      de;                             // Dividend payout rate energy sector
unordered_map<string, double> a1_mh;        // Propensity to consume out of wage income
unordered_map<string, double> a2_mh;        // Propensity to consume out of dividend & interest income
unordered_map<string, double> a3_mh;        // Propensity to consume out of wealth (deposits)
unordered_map<string, double> a4_mh;        // Propensity to consume out gov transfers (excluded unemployment benefits)
unordered_map<string, double> a5_mh;        // Propensity to consume out unemployment benefits
double      f2_entry_min;                   // Minimum market share for entering C-firms
double      kappa;                          // Persistence of long-run average mean productivity change
double      taylor;                         // Persistence parameter for CB rate adjustment
double      omicron;                        // Desired inventory to output ratio for C-firms
double      I_max;                          // Maximum desired expansion investment as % of current capital stock
double      persistence;                    // Persistence of one-off climate shocks
double      omega3;                         // Minimum market share (as % of previous share) remaining post-replicator
double      d_f;                            // Share of fossil fuel cost transfered to households
unordered_map<string, double> d_f_sh_mh;    // Share of fossil fuel total transfer that goes to each houseold class
double      g_ls;                           // Growth rate of labour force
double      aliqee;                         // Rate of tax on excess energy profit during energy price shock
double      aliqef;                         // Rate of tax on excess fossil fuel profit during energy price shock
double      tre;                            // Household transfer payment rate during energy price shock
unordered_map<string, double> Transfer_shock_sh_mh; //Shares of transfer following an energy price shock that goes to each households class
double      tref;                           // Firm transfer payment rate during energy price shock
double      passthrough;                    // Probability that a firm will pass through the energy price shock
vector<string>                classes_mh;   // Keys for households classes
unordered_map<string, double> LS_sh_mh;     // Household classes population shares
unordered_map<string, double> w_ratios_mh;  // Household classes wage ratios relative to Workers
unordered_map<string, double> ld_ratios_mh; // Ratios of labour demand for production based on number of Workers
unordered_map<string, double> Dividends_sh_mh; // Share of dividends to each household class
unordered_map<string, double> Ownership_sh_e_mh; // Share of dividends from energy setor to each household class
unordered_map<string, double> Ownership_sh_ff_mh; // Share of dividends from fossil fuel setor to each household class
unordered_map<string, double> Ownership_sh_b_mh; // Share of dividends from banks to each household class
unordered_map<string, double> Ownership_sh_1_mh0; // Inital share of ownership of each K-firm of each household class
unordered_map<string, double> Ownership_sh_2_mh0; // Inital share of ownership of each C-firm of each household class
unordered_map<string, double> energy_expenditure_sh_mh; // Share of expenditure spent in energy by each households class     
double      redistribute_co2TaxRev;             // Share of carbon tax revenue redistributed to households overall
unordered_map<string, double> co2TaxRev_sh_mh;  // Share of carbon tax revenue redistributed to each household class
double      profit_share_energy_inv;            // Profit share of nominal energy investment 
double      c_price_increase_perc;              // Percentual increase in C-firms markup (only if flag_c_price_increase!=0)
double      aliqdiv_regime_shift;            // If flag_change_aliqdiv==1, new tax rate on dividends after t_regime_shift; If flag_change_aliqdiv==2, share of personal income taxation shifting from labour to capital income taxation;  
unordered_map<string, double>  aliqw_regime_shift_mh;  // New wages tax rates (if flag_change_wage_tax_rate==1) or new ratios of wages tax rate compared to Workers (if flag_change_wage_tax_rate==2) after t_regime_shift 
unordered_map<string, double>  aliqwealth_regime_shift_mh;  // New wealth tax rates (if flag_change_wealth_tax_rate==1)
double      demand_persistency_h;              // Persistency in households demand of energy and goods

//Government
double      g_cons_GDP_ratio;               // Desired public consumption/GDP ratio
double      g_cons_GDP_ratio_regime_shift;  // If flag_change_public_spending==1 desired public consumption/GDP ratio after regime shift; If flag_change_public_spending==2 share of additional tax income to spend

//Energy sector
double      share_RD_en;                    // Share of energy sector revenues devoted to R&D
double      share_de_0;                     // Share of energy R&D devoted to dirty energy
int         payback_en;                     // Payback period for green energy plants
int         life_plant;                     // Maximum lifespan of green energy plants
double      exp_quota;                      // Maximum expansion of green energy capacity
double      o1_en;                          // Effectiveness of R&D expenditure in energy
double      uu1_en;		                    // Lower bound for innovation gain in energy sector
double      uu2_en;		                    // Upper bound for innovation gain in energy sector
double      ratio_c_en_h_firms;             // Ratio of energy price paid by households and firms
double      ratio_mi_en_shock;              //Ratio of new markup over old one when energy price shock through energy sector markup hits (flag_energyshocks==1 or 2)

//Climate
double      exp_quota_param;                // Governs the evolution of the constraint on green energy investment when flag_endogenous_exp_quota=1
RowVector   a_0;             		        // Governing influence of temperature anomaly on disaster generating function
RowVector   b_0;                            // Governing influence of temperature variance on disaster generating function
int         nshocks;                        // Number of different climate shock channels
int         t_start_climbox;  	            // Period in which the climate box is activated
double      T_pre;            	            // Pre-industrial global mean surface temperature
double      intercept_temp;		            // Intercept of temperature anomaly function
double      slope_temp;			            // Slope of temperature anomaly function
double      tc1;                            // Scaling parameter for linearly-growing carbon tax
double      tc2;                            // Growth rate of carbon tax for exponential case
int         ndep;                           // Number of ocean layers
RowVector   laydep;                         // Depth of ocean layers (in meters); beginning with top layer
double      fertil;                         // Carbon fertilisation effect of atmospheric carbon on NPP 
double      heatstress;                     // Heat stress effect of warming on NPP
double      humtime;                        // Decaying time of humus (in years)
double      biotime;                        // Decay time of carbon in biosphere (in years)
double      humfrac;                        // Fraction of decaying biomass carbon that ends up in humus
double      eddydif;                        // Eddy diffusion coefficient in m^2/year
double      Conref;                         // Reference amount of carbon in upper ocean layer (= pre-industrial) for atmospheric carbon exchange
double      ConrefT;                        // Influence of temperature on Conref
double      rev0;                           // Standard Revelle factor 
double      revC;                           // Impact of carbon content on revelle factor 
int         niterclim;                      // Number of iterations in the carbon exchange between ocean and atmosphere
double      forCO2;                         // Radiative forcing from e-folding CO2 (W/m^2)
double      otherforcefac;                  // Factor by which CO2-induced forcing is multiplied to account for non-CO2 forcing; if this is changed, initial conditions must also be adjusted
double      outrad;                         // W/m^2/K; outgoing radiation per 1K of surface warming 
double      seasurf;                        // Fraction of planetary surface covered by sea
double      heatcap;                        // J/m^3; heat capacity of water
double      secyr;                          // Seconds per year 
int         freqclim;                       // Climate box is called every time after freqclim time-steps have passed in the economic model (set to 4 if econ model is quarterly)
double      g_emiss_global;                 // Exogenous growth rate of global emissions (exluding endogenous EU emissions)
double      emiss_share;                    // Initial share of EU in global emissions
RowVector   shockexponent1;                 // Exponents used for endogenous change in shape of "disaster generating function"
RowVector   shockexponent2;                 // Exponents used for endogenous change in shape of "disaster generating function"
double      a2_nord;                        // Parameter for DICE-type damage function
double      sd_nord;                        // Governing standard deviation of DICE-based shocks at agent level

//Experiments & policies
int         t_regime_shifts;                // Timestep to trigger policy and/or regime shifts
unordered_map<string, double> income_sh_regime_shift_mh; // New income shares in case of income changes (only if flag_change_income_shares==1)
unordered_map<string, double> a1_mh_regime_shift; // New average propensities to consume out of wage & benefit income in case of their change (only if flag_change_apc==1)
unordered_map<string, double> a2_mh_regime_shift; // New average propensities to consume out of dividend & interest income in case of their change (only if flag_change_apc==1)
unordered_map<string, double> a3_mh_regime_shift; // New average propensities to consume out of wealth (deposits) in case of their change (only if flag_change_apc==1)
unordered_map<string, double> en_exp_sh_mh_regime_shift; // New share of expenditure spent in energy in case of their change (only if flag_change_en_exp_shares==1)
double      env_subsidy_per_machine;        // Defines environmnetal subsidy based on value of flag_environmental_subsidies_C_firms;
double correlation_prod_and_green_tech;     // Defines how much a change in carbon intensity transfers into a change in labour productivity (only if flag_correlate_prod_and_green_tech == 1 or 2)
double rs_uu_lp;                            // When flag_correlate_prod_and_green_tech==4, re-scaling of labour productivity change intervals for discovered technology
double rs_uu_ee;                            // When flag_correlate_prod_and_green_tech==4, re-scaling of energy efficiency change intervals for discovered technology
double rs_uu_ef;                            // When flag_correlate_prod_and_green_tech==4, re-scaling of environmental friendliness change intervals for discovered technology
double K_ge_END_perc;                       // When flag_energy_exp = 3 or 4, final green energy share after t_length_energy_transition timesteps
int t_length_energy_transition;             // When flag_energy_exp = 3 or 4, number of timesteps through which green energy share increases from K_ge0_perc to K_ge_END_perc
double renew_impact_on_p_e;                 // WHen flag_energy_exp = 4, reduction of energy price for each additional point of green energy share
double bonuses_share;                       // Fraction of additional profits that is paid to Top 10% households as bonuses by firms, banks and the energy sector

#endif
