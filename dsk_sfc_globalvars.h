// Define global variables
#ifndef GLOBAL_VAR_H
#define GLOBAL_VAR_H

// #include "dsk_sfc_parameters.h"
// #include "dsk_sfc_inits.h"

#include <map>

int              reducedoutput;                                 // Dummy indicating whether full output is saved
long int         seed;				                         // Seed for random number generation
long int         *p_seed;			                         // Pointer to seed
RowVector        shocks_kfirms;                              // Vector of climate shocks to K-firms
RowVector        shocks_cfirms;                              // Vector of climate shocks to C-firms
double           shock_scalar;                               // Scalar climate shock
ofstream         Errors;                                     // File for saving error messages                 
int              cerr_enabled;                               // Determines whether error messages to console are printed to console
int              verbose;                                    // Determines whether simulation progress updates are printed to console

int              i;	                                         // Index
int              ii;	                                     // Index
int              iii;					                     // Index
int              j;                                          // Index               
int              jjj;                                        // Index
int              t;                                          // Current simulation period
int              tt;                                         // Index
int              rni;                                        // Random number
int              t0;                                         // Index 
int              t00;                                        // Index
int              n;                                          // Counter
int              iterations;                                 // Counter
double           pareto_rv;                                  // Pareto random number
double           tolerance;                                  // Tolerance level for deviations from accounting consistency
double           deviation;                                  // Deviation from accounting consistency
double           parber;                                     // Input for draw from Bernoulli
double           rnd;                                        // Random number
double           N1r;						                 // Number of K-firms as double
double           N2r;                                        // Number of C-firms as double
int              step;                                       // Counter
int              stepbis;                                    // Counter
int              cont;                                       // Counter               
double           age0;                                       // Age of initial machines
double           Amax;                                       // Maximum labour productivity C-firms
double           A1pmax;                                     // Maximum labour productivity K-firms
double           A1_en_max;                                  // Maximum energy efficiency C-firms
double           A1_ef_max;                                  // Maximum environmental friendliness C-firms
double           A1p_en_max;                                 // Maximum energy efficiency K-firms
double           A1p_ef_max;                                 // Maximum environmental friendliness K-firms
double           D20;                                        // Initial demand for consumption goods of households overall
double           Dtot0;                                      // Initial expenditure of households overall
unordered_map<string, double> Den0_mh;                       // Initial demand for consumption goods of each household class
double           Den0;                                       // Initial demand for energy
int              DS2_min_index;                              // Index of minimum debt service to sales ratio
int              newbroch;                                   // Brochures sent to potential customers
int              indforn;                                    // C-firm's machine supplier
int              flag;                                       // Flag
double           payback;                                    // Payback variable
int              jmax;                                       // Index
int              tmax;                                       // Index
int              imax;                                       // Index
double           nmachprod;                                  // Number of machines used
double           nmp_temp;                                   // Number of machines used (temporary value)
double           cmin;                                       // Index
int              imin;                                       // Index
int              jmin;                                       // Index
int              tmin;                                       // Index
double           MaxFunds;                                   // Maximum liquid funding C-Firm expects to have
double           prestmax;                                   // Maximum loans demanded by C-firm
double           p1prova;                                    // Temporary storage for K-firm price         
int              rated_firm_2;                               // Index
double           Qpast;                                      // Temporary storage for K-firm output
double           Ipast;                                      // Temporary storage for investment
double           scrapmax;                                   // Maximum number of machines to be scrapped
double           cmax;                                       // Maximum production cost
int              ind_i;                                      // Index
int              ind_tt;                                     // Index
double           scrap_n;                                    // Nominal value of scrapped machines
int              sendingBank;                                // Index
int              receivingBank;                              // Index
double           c_de_min;                                   // Minimum cost of dirty energy plant
double           cf_min_ge;                                  // Minimum cost of green energy capacity expansion
RowVector        G_de_temp;                                  // Temporary storage for dirty energy capacity
double           Q_de_temp;                                  // Temporary storage for dirty energy produced
int              idmin;                                      // Index 
double           parber_en_de;                               // Input for draw from Bernoulli
double           parber_en_ge;                               // Input for draw from Bernoulli
double           cpi_temp;                                   // Temporary cpi
RowVector        Bond_share;                                 // Share of gov. bonds owned by each bank
int              maxbank;                                    // Index
double           max_equity;                                 // Maximum bank equity
double           multip_bailout;                             // Multiplier for bailout
double           min_equity;                                 // Minimum bank equity
double           multip_entry;                               // Multiplier for firm entry
double           injection;                                  // Liquidity injection for entering firms by hosueholds overall
double           injection2;                                 // Alternative storage for liquidity injection for entering firms
unordered_map<string, double> injection_mh;                  // Liquidity injection for entering firm by each household class
double           n_mach_exit;                                // Number of machines of exiting firms
double           n_mach_exit2;                               // Temporary storage for n_mach_exit
double           n_mach_needed;                              // Number of machines needed for entering firms
double           n_mach_resid;                               // Number of remaining machines for entering firms
double           n_mach_resid2;                              // Temporary storage for remaining machines for entering firms
double           n_exit2;                                    // Number of exiting C-firms
RowVector        k_entry;                                    // Share of second-hand capital allocated to each firm
double           cpi_init;                                   // Initial value of cpi needed for climate policy
double           GDP_init;                                   // Initial value of nominal GDP needed for climate policy
double           cpi_t_regime_shifts;                        // Value of cpi at regime shifts, required for carbon tax increment
double           baddebt_2_temp;                             // Temporary storage for C-firms' bad debt
double           markdownCapital;                            // Markdown applied to capital goods sold on second hand market
RowVector        prior;                                      // Used to determine deviations in bank balance sheets
double           post;                                       // Used to determine deviations in bank balance sheets
double           prior_cb;                                   // Used to determine deviations in CB balance sheet
double           post_cb;                                    // Used to determine deviations in CB balance sheet
double           DepositsCheck_1;                            // Used to detect errors in distribution of firm deposits
double           DepositsCheck_2;                            // Used to detect errors in distribution of firm deposits
double           p2_entry;                                   // Price of newly entering C-firms
double           f2_exit;                                    // Sum of market shares of exiting C-firms
double           CurrentDemand;                              // Sum of demand experienced by C-firms in current period
RowVector        EntryShare;                                 // Share of available market captured by entering C-firms
RowVector        CompEntry;                                  // Pseudo-competitiveness of entering C-firms
double           CompEntry_m;                                // Mean pseudo-competitiveness of entering C-firms
double           K_gap;                                      // Gap between desired and actual capital
RowVector        K_temp;                                     // Temporary storage for capital stock
double           K_top;                                      // Upper limit for expansion investment
int              loss;                                       // Temporary storage for capital stock lost due to climate shocks
int              lossj;                                      // Temporary storage for capital stock lost at firm level
RowVector        K_loss;                                     // Temporary storage for capital stock lost at firm level
RowVector        C_loss;                                     // Temporary storage for output lost of C-Firms
RowVector        I_loss;                                     // Temporary storage for output lost of K-Firms
int              rani;                                       // Random integer
int              rant;                                       // Random integer
int              ranj;                                       // Random integer
double           reduction;                                  // Temporary storage for reduction in production due to insufficient labour supply
RowVector        marker_age;                                 // Indicates whether a firm has only 1 unit of capital left which is also older than agemax
double           K_temp_sum;                                 // Temporary storage for overall capital stock   
double           mi_en_preshock;                             // Temporary storage for pre-shock energy markup
double           pf_preshock;                                // Temporary storage for pre-shock fossil fuel price
double           mi_en_shock;                                // Temporary storage for shocked energy markup
double           c_en_preshock;                              // Temporary storage for pre-shock energy price
double           pf_shock;                                   // Temporary storage for shocked fossil fuel price
double           c_infra_t;                                  // Target inframarginal energy cost to achieve desired energy price shock
double           ptemp;                                      // Temporary storage for C-firm price
RowVector        pass_1;                                     // Indicator for K-firm passthrough of energy price shock
RowVector        pass_2;                                     // Indicator for C-firm passthrough of energy price shock
double           Ldtemp;                                     // Temporary storage for C-Firm labour
string           tech;                                       // Storage for name of energy technology
double           old_capitalStock;                           // Storage for existing energy capital stock
double           sumDepEn;                                   // Used for SFC check of energy sector deposits

// Balance sheet items
RowVector        Deposits_h(2);                              // Household total deposits
unordered_map<string, RowVector> Deposits_mh;                // Household classes deposits
RowVector        Deposits_e(2);                              // Energy Sector deposits
RowVector        Loans_e(2);                                 // Aggregate loans to the energy sector
Matrix           Deposits_1;                                 // K-firm deposits
Matrix           Deposits_2;                                 // C-firm deposits
Matrix           Deposits_b;                                 // Total deposits from banks' side
Matrix           Deposits_hb;                                // Household deposits from banks' side
Matrix           Deposits_eb;                                // Energy sector deposits from banks' side
Matrix           GB_b;                                       // Government bonds held by banks 
RowVector        GB_cb(2);                                   // Government bonds held by CB
RowVector        GB(2);                                      // Government bonds
RowVector        Deposits_fuel(2);                           // Deposits of fossil fuel sector with CB
RowVector        Deposits_fuel_cb(2);                        // Deposits of fossil fuel sector from CB side
Matrix           Loans_2;                                    // Loans of C-firms
Matrix           Loans_b;                                    // Loans to C-firms and energy sector from banks' side
Matrix           Advances_b;                                 // CB advances to banks
RowVector        Advances(2);                                // CB advances
Matrix           Reserves_b;                                 // CB reserves held by banks
RowVector        Reserves(2);                                // CB Reserves
Matrix           CapitalStock;                               // Nominal value of C-firms' capital stock (machines)
Matrix           deltaCapitalStock;                          // Change in nominal value of C-firms' capital stock
RowVector        CapitalStock_e(2);                          // Nominal value of Energy sectors' capital stock (plants)
Matrix           Inventories;                                // Nominal value of C-firms' inventories
RowVector        NW_h(2);                                    // Net worth of households
unordered_map<string, RowVector> NW_mh;                      // Net worth of each household class
Matrix           NW_1;                                       // Net worth of K-firms
Matrix           NW_2;                                       // Net worth of C-firms
Matrix           NW_b;                                       // Net worth of banks
RowVector        NW_gov(2);                                  // Net worth of Government
RowVector        NW_cb(2);                                   // Net worth of CB
RowVector        NW_e(2);                                    // Net worth of Energy sector
RowVector        NW_f(2);                                    // Net worth of fossil fuel sector

double           NW_h_c;                                     // Net worth of households (control; for SFC-check)
unordered_map<string, double> NW_mh_c;                       // Net worth of household classes (control; for SFC-check)
RowVector        NW_1_c;                                     // Net worth of K-firms (control; for SFC-check)
RowVector        NW_2_c;                                     // Net worth of C-firms (control; for SFC-check)
RowVector        NW_b_c;                                     // Net worth of banks (control; for SFC-check)
double           NW_gov_c;                                   // Net worth of Government (control; for SFC-check)
double           NW_cb_c;                                    // Net worth of CB (control; for SFC-check)
double           NW_e_c;                                     // Net worth of energy sector (control; for SFC-check)
double           NW_f_c;                                     // Net worth of fossil fuel sector (control; for SFC-check)
double           NWSum;                                      // Sum of net worths
double           RealAssets;                                 // Nominal value of all real assets in the economy

//Additional TFM Items
RowVector        Wages_2_i;                                  // Wages paid by each C-firm to households overall
unordered_map<string, double> Wages_2_mh;                    // Wages paid by all C-firms to each household class (only for reporting)
unordered_map<string, RowVector> Wages_2_i_mh;               // Wages paid by each C-firm to each household class
RowVector        Wages_1_i;                                  // Wages paid by each K-firm to households overall
unordered_map<string, double> Wages_1_mh;                    // Wages paid by all K-firm to each household class (only for reporting)
unordered_map<string, RowVector> Wages_1_i_mh;               // Wages paid by each K-firm to each household class
double           Wages_en;                                   // Wages paid by Energy sector to households overall
unordered_map<string, double> Wages_en_mh;                   // Wages paid by Energy sector to each household class
double           Wages;                                      // Wages received by households (total)
unordered_map<string, double> Wages_mh;                      // Total wages received by each household class
double           Benefits;                                   // Unemployment benefits received by households
unordered_map<string, double> Benefits_mh;                   // Unemployment benefits received by each household class
RowVector        Dividends_h(2);                             // Dividend payments received by households overall (1-lag)
unordered_map<string, RowVector> Dividends_mh;               // Dividend payments received by each household class
double           Bonuses_h;                                  // Bonuses received by households overall
unordered_map<string, RowVector> Bonuses_mh;                 // Bonuses received by each household class
RowVector        Dividends_1_i;                              // Dividends paid by each K-firm overall
unordered_map<string, RowVector> Dividends_1_i_mh;           // Dividends paid by each K-firm to each household class
double           Dividends_1;                                // Total K-firm dividends paid to households overall
unordered_map<string, double> Dividends_1_mh;                // Total K-firm dividends paid to each household class
RowVector        Bonuses_1_i;                                // Bonuses paid by each K-firm to managers
double           Bonuses_1;                                  // Bonuses paid by all K-firms to managers
unordered_map<string, double> Bonuses_1_mh;               // Bonuses paid by all K-firms to each household class
RowVector        Dividends_2_i;                              // Dividends paid by each C-firm overall
unordered_map<string, RowVector> Dividends_2_i_mh;           // Dividends paid by each C-firm to each household class
double           Dividends_2;                                // Total C-firm dividends paid to households overall
unordered_map<string, double> Dividends_2_mh;                // Total C-firm dividends paid to each household class
RowVector        Bonuses_2_i;                                // Bonuses paid by each C-firm to managers
double           Bonuses_2;                                  // Bonuses paid by all C-firms to managers
unordered_map<string, double> Bonuses_2_mh;               // Bonuses paid by all C-firms to each household class
RowVector        Dividends_b_i;                              // Dividends paid by each bank to households overall
double           Dividends_b;                                // Dividends paid by all banks to households overall
unordered_map<string, double> Dividends_b_mh;                // Dividends paid by all banks to each household class
RowVector        Bonuses_b_i;                                // Dividends paid by each bank to managers
double           Bonuses_b;                                  // Dividends paid by all banks to managers
unordered_map<string, double> Bonuses_b_mh;                // Dividends paid by all banks to each household class
double           Dividends_e;                                // Dividends paid by energy sector overall
unordered_map<string, double> Dividends_e_mh;                // Dividends paid by energy sector to each household class
double           Bonuses_e;                                  // Bonuses paid by energy sector to managers
unordered_map<string, double> Bonuses_e_mh;                  // Bonuses paid by energy sector to each household class
RowVector        Investment_2;                               // Investment expenditures of C-firms
RowVector        EnergyPayments_1;                           // Payments of K-firms for energy
RowVector        EnergyPayments_2;                           // Payments of C-firms for energy
double           EnergyPayments;                             // Energy payments received by energy sector
RowVector        Taxes_1;                                    // Taxes on profits paid by K-firms
RowVector        Taxes_2;                                    // Taxes on profits paid by C-firms
double           Taxes_e;                                    // Taxes on profits paid by the energy sector
double           Taxes_h;                                    // Total taxes paid by households overall
unordered_map<string, double> Taxes_mh;                      // Total taxes paid by each household class
double           Taxes_w_h;                                  // Taxes on wage paid by households overall
unordered_map<string, double> Taxes_w_mh;                    // Taxes on wage paid by each household class
double           Taxes_div_h;                                // Taxes on dividends paid by households overall
unordered_map<string, RowVector> Taxes_div_mh;               // Taxes on all dividends paid by each household class (1-lag for consumption function)
double           Taxes_div_1_h;                                // Taxes on K-firms dividends paid by households overall
unordered_map<string, double> Taxes_div_1_mh;                // Taxes on K-firms dividends paid by each household class
double           Taxes_bon_h;                                // Taxes on bonuses paid by households overall
unordered_map<string, RowVector> Taxes_bon_mh;               // Taxes on bonuses paid by each household class (1-lag for consumption function)ù
double           Taxes_bon_1_h;                              // Taxes on K-firms bonuses paid by households overall
unordered_map<string, double> Taxes_bon_1_mh;                // Taxes on K
double           Taxes_wealth_h;                             // Taxes on wage paid by households overall
unordered_map<string, double> Taxes_wealth_mh;               // Taxes on wage paid by each household class
RowVector        Taxes_b;                                    // Taxes on profits paid by banks
RowVector        Taxes_CO2_1;                                // Taxes on CO2 paid by K-firms
RowVector        Taxes_CO2_2;                                // Taxes on CO2 paid by C-firms
double           Taxes_CO2_e;                                // Taxes on CO2 paid by energy sector
RowVector        Taxes_CO2(2);                               // Total CO2 taxes collected by the government
RowVector        InterestDeposits_1;                         // Deposit interest paid to K-firms
RowVector        InterestDeposits_2;                         // Deposit interest paid to C-firms
RowVector        InterestDeposits;                           // Deposit interest paid by banks
double           InterestDeposits_h;                         // Deposit interest paid to households overall
unordered_map<string, double> InterestDeposits_mh;           // Deposit interest paid to each household class
double           InterestDeposits_e;                         // Deposit interest paid to energy sector
double           InterestBonds;                              // Bond interest paid by government
double           InterestBonds_cb;                           // Bond interest paid to CB
RowVector        InterestBonds_b;                            // Bond interest paid to banks
double           BondRepayments_cb;                          // Bond repayments made to CB
RowVector        BondRepayments_b;                           // Bond repayments made to banks
double           Taxes_g;                                    // Taxes collected by government
double           Taxes_tot_g;                                // Taxes collected by government (including CO2 and energy shock)           
RowVector        LoanInterest;                               // Loan interest paid to banks
RowVector        LoanInterest_2;                             // Loan interest paid by C-firms
double           FirmTransfers;                              // Transfer payments for firm entry made by households _mh we need to differentiate classes here
double           FirmTransfers_1;                            // Transfer payments received by entering K-firms
double           FirmTransfers_2;                            // Transfer payments received by entering C-firms
RowVector        Injection_1;                                // Total injection of net worth for entering K-firms
RowVector        Injection_2;                                // Total injection of net worth for entering C-firms
double           InterestReserves;                           // Interest on reserves paid by CB 
double           InterestAdvances;                           // Interest on advances received by CB
RowVector        InterestReserves_b;                         // Interest on reserves received by banks
RowVector        InterestAdvances_b;                         // Interest on advances paid by banks
double           TransferCB;                                 // Transfers of CB profits to government 
double           FuelCost;                                   // Cost of fossil fuels for Energy sector
double           TransferFuel;                               // Transfer from fossil fuel sector to households
unordered_map<string, double> TransferFuel_mh;               // Transfer from fossil fuel sector to each household class
double           Taxes_e_shock;                              // Excess profit tax on energy sector following an energy price shock through increased markup
double           Taxes_f_ff_shock;                           // Excess profit tax on fossil fuel sector following an energy price shock through fossil price
double           Taxes_e_ff_shock;                           // Excess profit tax on energy sector following an energy price shock through fossil price
double           Transfer_shock;                             // Transfer to households to combat energy price shock
unordered_map<string, double> Transfer_shock_mh;            // Transfer to each household class to combat energy price shock
double           Transfer_shock_f;                           // Transfer to firms to combat energy price shock
RowVector        Transfer_shock_f1;                          // Transfer to K-Firms
RowVector        Transfer_shock_f2;                          // Transfer to C-Firms

double           Balance_h;                                  // Sectoral balance households
unordered_map<string, double> Balance_mh;                    // Sectoral balance household classes
double           Balance_1;                                  // Sectoral balance K-firms
double           Balance_2;                                  // Sectoral balance C-firms
double           Balance_e;                                  // Sectoral balance energy sector
double           Balance_b;                                  // Sectoral balance banks
double           Balance_g;                                  // Sectoral balance government
double           Balance_cb;                                 // Sectoral balance CB
double           Balance_f;                                  // Sectoral balance fossil fuels
double           BalanceSum;                                 // Sum of sectoral balances
RowVector        Balances_1;                                 // Individual balances K-firms

//Households
unordered_map<string, RowVector> w_mh;                       // Household classes LD1_wr rate (1-lag)
RowVector        w_tot_for_1_wr_mh(2);                       // Wages to pay in total for production for each worker employed (total cost of unitary productive labour) (1-lag)
double           LS;                                         // Total labour supply
unordered_map<string, double> LS_mh;                         // Household classes labour supply
RowVector        U(2);                                       // Unemployment rate overall
unordered_map<string, double> U_mh;                          // Unemployment rate of each househols class
double           Exp_li;                                     // Overall households' total expenditure demand out of labour income
double           Exp_ki;                                     // Overall households' total expenditure demand out of capital income
double           Exp_dep;                                    // Overall households' total expenditure demand out of deposits
double           Exp_govTransf;                              // Overall households' total consumption demand out of government transfers from carbon tax revenues
double           Exp_u_benefit;                              // Overall households' total consumption demand out of unemployment benefit
unordered_map<string, double> Exp_li_mh;                    // Household classes expenditure out of labour income
unordered_map<string, double> Exp_ki_mh;                    // Household classes expenditure demand out of capital income
unordered_map<string, double> Exp_dep_mh;                   // Household classes expenditure demand out of deposits
unordered_map<string, double> Exp_govTransf_mh;             // Household classes consumption demand out of government transfers from carbon tax revenues
unordered_map<string, double> Exp_u_benefit_mh;             // Household classes consumption demand out of unemployment benefit
double           Cons_h;                                     // Overall households' total desired C-goods demand (before constrains on depostis and supply)
unordered_map<string, double> Cons_mh;                       // Household classes total C-goods consumption demand
unordered_map<string, double> Cons_sh_mh;                    // Household classes C-goods consumption shares
double           Cons;                                       // Total C-goods consumption demand (households + government)
double           Cres;                                       // Residual C-goods consumption
double           Cresbis;                                    // Temporary storage for residual C-goods consumption
double           Deposits_recovered_1;                       // Liquidity recovered from failing K-firms from households overall
double           Deposits_recovered_2;                       // Liquidity recovered from failing C-firms from households overall
unordered_map<string, double> Deposits_recovered_1_mh;       // Liquidity recovered from failing K-firms from each household class
unordered_map<string, double> Deposits_recovered_2_mh;       // Liquidity recovered from failing C-firms from each household class
RowVector        Income_gross_h(2);                          // Households' total income before taxes (only for output)
unordered_map<string, RowVector> Income_gross_mh;            // Household classes total income before taxes (only for output)
unordered_map<string, double> r_wealth_mh;                   // Household classes return on wealth (only for output)
double           r_wealth_h;                                 // Households overall return on wealth (only for output)
unordered_map<string, double> income_gross_growth_mh;        // Household classes gross income growth rate (only for output)
double           income_gross_growth_h;                      // Households overall gross income growth rate (only for output)
int              n_classes_mh;                               // Number of households' classes
double           Consumption;                                // C-goods final total consumption expenditure overall (government + households)
double           Consumption_g;                              // C-goods final consumption expenditure of government
double           Consumption_h;                              // C-goods final/actual consumption expenditure of households overall (after constrains on supply and deposits)
unordered_map<string, double> Consumption_mh;                // C-goods consumption expenditure of each household class
double           Consumption_r;                              // Total real C-goods consumption
unordered_map<string, RowVector> Ownership_sh_1_i_mh;        // Share of ownership of each K-firm of each household class
unordered_map<string, RowVector> Ownership_1_i_mh;           // Ownership of each K-firm of each household class
unordered_map<string, double> Ownership_1_mh;                // Ownership of all K-firms of each household class
unordered_map<string, double> Ownership_sh_1_mh;             // Share of ownership of all K-firms of each household class
unordered_map<string, RowVector> Ownership_sh_2_i_mh;        // Share of ownership of each C-firm of each household class
unordered_map<string, RowVector> Ownership_2_i_mh;           // Ownership of each C-firm of each household class
unordered_map<string, double> Ownership_2_mh;                // Ownership of all C-firms of each household class
unordered_map<string, double> Ownership_sh_2_mh;             // Share of ownership of all C-firms of each household class
unordered_map<string, double> Ownership_e_mh;                // Ownership of energy sector of each household class
unordered_map<string, double> Ownership_b_mh;                // Ownership of all banks of each household class
double           Ownership_tot_h;                            // Ownership of all sectors of each household class
unordered_map<string, double> Ownership_tot_mh;              // Ownership of all sectors of households overall
unordered_map<string, double> Deposits_sh_mh;                // Household classes deposit shares
unordered_map<string, double> Entry_financing_sh_mh;         // Household classes shares of contribution to entry of C-/K- firms 
unordered_map<string, double> Expenditure_en_mh;             // Expenditure in energy of each household class
double           Expenditure_en_h;                           // Expenditure in energy of all households 
unordered_map<string, double> D_en_mh;                       // Energy demand of each household class
double           D_en_h;                                     // Energy demand of all households 
unordered_map<string, double> Expenditure_tot_mh;            // Total expenditure of each household class
double           Expenditure_tot_h;                          // Total expenditure of all households
unordered_map<string, double> cons_dir_en_fp_mh;             // Energy footprint from direct energy demand of each household class
unordered_map<string, double> cons_indir_en_fp_mh;           // Energy footprint from C-goods consumption of each household class
unordered_map<string, double> invest_en_fp_mh;               // Energy footprint from investment of each household class
unordered_map<string, double> publ_cons_en_fp_mh;            // Energy footprint from public spending of each household class
unordered_map<string, double> en_fp_mh;                      // Total energy footprint of each household class
double           cons_dir_en_fp_h;                           // Energy footprint from direct energy demand of huseholds overall
double           cons_indir_en_fp_h;                         // Energy footprint from C-goods consumption of huseholds overall
double           invest_en_fp_h;                             // Energy footprint from investment of huseholds overall
double           publ_cons_en_fp_h;                          // Energy footprint from public spending of huseholds overall
double           en_fp_h;                                    // Total energy footprint of household overall
unordered_map<string, double> cons_dir_carb_fp_mh;           // Carbon footprint from direct energy demand of each household class
unordered_map<string, double> cons_indir_carb_fp_mh;         // Carbon footprint from C-goods consumption of each household class
unordered_map<string, double> invest_carb_fp_mh;             // Carbon footprint from investment of each household class
unordered_map<string, double> publ_cons_carb_fp_mh;          // Carbon footprint from public spending of each household class
unordered_map<string, double> carb_fp_mh;                    // Total carbon footprint of each household class
double           cons_dir_carb_fp_h;                         // Carbon footprint from direct energy demand of huseholds overall
double           cons_indir_carb_fp_h;                       // Carbon footprint from C-goods consumption of huseholds overall
double           invest_carb_fp_h;                           // Carbon footprint from investment of huseholds overall
double           publ_cons_carb_fp_h;                        // Carbon footprint from public spending of huseholds overall
double           carb_fp_h;                                  // Total carbon footprint of household overall
double           effective_en_exp_sh_h;                      // Effective energy expenditure share of households overall (only for reporting)
unordered_map<string, double> effective_en_exp_sh_mh;        // Effective energy expenditure share of each household class (only for reporting)
double           en_int_energy_exp;                          // Energy intensity of households' energy expenditure
double           en_int_goods_exp;                           // Energy intensity of households' goods expenditure
unordered_map<string, double> en_int_exp_mh;                 // Energy intensity of each household class expenditure (only for reporting)
double           en_int_exp_h;                               // Energy intensity of households' expenditure overall
unordered_map<string, double> tax_CO2_incidence_C_mh;        // Incidence of carbon tax paid by C-firms on each household class
unordered_map<string, double> tax_CO2_incidence_en_mh;       // Incidence of carbon tax paid by energy sector on each household class
unordered_map<string, double> tax_CO2_incidence_on_income_mh;  // Incidence of carbon tax passed to each household class over relative to its income


//C-Firms
int              ns2;                                        // Number of active C-firms
double           mD2;                                        // Mean deposits C-firms
RowVector        p2;                                         // C-firms' prices
RowVector        BankingSupplier_2;                          // C-firms' suppliers of banking services
Matrix           BankMatch_2;                                // Matrix matching C-firms to banks
RowVector        A2;                                         // Labour productivity C-firms
RowVector        A2_mprod;                                   // Labour productivity C-firms for computing mean productivity
RowVector        A2_en;			                             // Energy efficiency C-firms
RowVector        A2_ef;		                                 // Environmental friendliness C-firms
RowVector        c2;                                         // Production cost C-firms
RowVector        c2p;                                        // Production cost C-firms use to set prices (for energy price shock experiment)
RowVector        c1p;                                        // Production cost K-firms use to set prices (for energy price shock experiment)
RowVector        l2;                                         // Unsatisfied consumption demand by each firm
double           l2m;                                        // Mean unsatisfied demand
RowVector        p2m(2);                                     // Mean C-firm price (1-lag)
double           mu2m;                                       // Mean C-firm markup
RowVector        Q2temp;                                     // Temporary storage for C-firm output (produced + inventories)
RowVector        f_temp2;                                    // Temporary storage for C-firm market shares
RowVector        D_temp2;                                    // Temporary storage for consumption demand
RowVector        DebtRemittances2;                           // Loan repayments C-firms
RowVector        baddebt_2;                                  // Bad debt of exiting C-firms
Matrix           S2;                                         // C-firm revenues
RowVector        Sales2;                                     // Temporary storage for C-firm revenues
RowVector        K;                                          // C-firms' productive capacity
RowVector        K_cur;                                      // Needed for shocks to capital stock
Matrix           f2;				                         // C-firms' market share
RowVector        E2;                                         // C-firms' competitiveness
RowVector        CreditDemand;                               // C-firms' credit demand
RowVector        I;					                         // Total investment of C-firms in terms of productive capacity
Matrix           EI;	                                     // Expansion investment of C-firms in terms of productive capacity
RowVector        EI_n;			                             // Nominal value of C-firms' expansion investment
RowVector        SI;                                         // Substitution investment of C-firms in terms of productive capacity
RowVector        SI_n;                                       // Nominal value of C-firms' substitution investment
RowVector        Q2;                                         // Quantity produced by C-firms
Matrix           mu2;                                        // Mark-up of C-firms 
RowVector        fornit;                                     // C-firms' supplier of machine tools
RowVector        n_mach;                                     // C-firms' number of machines
Matrix           D2;				                         // Demand for consumption goods
RowVector        De;                                         // Expected demand
Matrix           N;                                          // Inventories (real)
RowVector        Ne;                                         // Desired inventories
RowVector        DebtServiceToSales2;                        // C-firms' debt service to sales ratio
Matrix           DebtService_2;                              // C-firms' debt service
Real             DS2_min;                                    // Minimum debt service to sales ratio
RowVector        DebtServiceToSales2_temp;                   // Temporary storage for C-firms' debt service to sales ratio
RowVector        k;                                          // Ranking of C-firms' debt service to sales ratio
RowVector        r_deb_h;                                    // Borrowing rate charged to individual C-firms
RowVector        FirmDefaultProbability;                     // Probability of default of individual C-firms estimated by banks
RowVector        EId;                                        // Desired expansion investment
RowVector        SId;                                        // Desired replacement investment
RowVector        EIp;                                        // Expansion investment post determination of maximum acceptable borrowing
RowVector        SIp;                                        // Replacement investment post determination of maximum acceptable borrowing
RowVector        Ip;                                         // Overall investment post determination of maximum acceptable borrowing
RowVector        Cmach;			                             // Cost of overall investment
RowVector        CmachEI;                                    // Cost of expansion investment
RowVector        CmachSI;                                    // Cost of replacement investment
RowVector        Qd;                                         // Quantity demanded from C-firms
RowVector        Kd;                                         // Desired capital stock of C-firms
RowVector        Ktrig;                                      // Current capital stock determining expansion investment  
RowVector        A2e;                                        // Effective labour productivity  
RowVector        c2e;                                        // Effective unit cost
RowVector        A2e_en;                                     // Effective energy efficiency (considering shocks and machines used)
RowVector        A2e_ef;                                     // Effective environmental friendliness
RowVector        A2e2;                                       // Needed for capital stock shocks
RowVector        A2e_en2;                                    // Needed for capital stock shocks
RowVector        A2e_ef2;                                    // Needed for capital stock shocks
RowVector        Ld2_wr;                                     // C-firms' labour demand for Workers
RowVector        Ld2_control;                                // Needed for capital stock shocks
unordered_map<string, RowVector> Ld2_i_mh;                   // Demand for each household class (key) of each C-firm (RowVector index)
unordered_map<string, double> LD2_mh;                        // Demand for each household class (key) of all C-firms
RowVector        mol;                                        // C-firms gross profit
RowVector        exiting_2;                                  // Indicating whether firm is exiting
RowVector        exit_payments2;                             // Indicating whether firm is exiting due to inability to make a payment
RowVector        exit_equity2;                               // Indicating whether firm is exiting due to negative equity
RowVector        exit_marketshare2;                          // Indicating whether firm is exiting due to loss of market share
RowVector        D2_en;		                                 // C-firms' energy demand	
RowVector        Emiss2;                                     // C-firms' emissions
RowVector        dN;                                         // Change in inventories
RowVector        dNm;                                        // Change in nominal value of inventories
Matrix           Pi2;                                        // C-firms' profit
RowVector        n_mach_entry;                               // Number of machines of entering firms
RowVector        scrap_age;                                  // Number of machines scrapped due to age
RowVector        EnvSubsidies_2_i;                           // Environmental subsidies to be paid to each C-firm
double           EnvSubsidies_2;                             // Environmental subsidies to be paid to C-firms overall
double           uc_L_avg_2;                                 // Average unit cost of labour for C-firms
double           uc_en_avg_2;                                // Average unit cost of energy for C-firms
double           uc_CO2_avg_2;                               // Average unit cost of emissions (due to CO2 tax) for C-firms
double           uc_avg_2;                                   // Average unit cost for C-firms
int              with_unsold_Q_2;                            // Number of C-firms with unsold output left (only for reporting)
int              with_neg_profits_2;                         // Number of C-firms with negative profits (only for reporting)

std::vector<std::vector<std::vector<int>>> age;              // Age of exiting machines
std::vector<std::vector<std::vector<double>>> g_c;           // Frequency of machines for cost calculation
std::vector<std::vector<std::vector<double>>> g_c2;          // Needed for shocks to capital stock
std::vector<std::vector<std::vector<double>>> g_c3;          // Needed for shocks to capital stock
std::vector<std::vector<std::vector<double>>> g;             // Frequency of machines
std::vector<std::vector<std::vector<double>>> g_price;       // Array containing original purchase prices of machines
std::vector<std::vector<std::vector<double>>> gtemp;         // Temporary storage for frequency of machines
std::vector<std::vector<std::vector<double>>> C_pb;          // Array containing production costs of machines to be scrapped
std::vector<std::vector<std::vector<double>>> g_pb;          // Array containing machines to be scrapped
std::vector<std::vector<double>> g_secondhand;               // Machines to be sold on second-hand market
std::vector<std::vector<double>> g_secondhand_p;             // Prices of machines to be sold on second-hand market
std::vector<std::vector<int>> age_secondhand;                // Age of machines to be sold on second-hand market

//K-firms
int              ns1;                                        // Number of active K-firms
double           mD1;                                        // Mean deposits K-firms
RowVector        p1;                                         // Prices of K-firms
double           p1_avg;                                     // Average price across K-firms weighted on revenues shares
RowVector        BankingSupplier_1;                          // K-firms' suppliers of banking services
Matrix           BankMatch_1;                                // Matrix matching K-firms to banks
RowVector        baddebt_1;                                  // Bad debt of exiting K-firms
RowVector        A1;                                         // Productivity of machines produced by K-firm
RowVector        A1p;                                        // Productivity of K-firm production technique
RowVector        A1_en;                                      // Energy Efficiency of machines produced by K-firm
RowVector        A1_ef;                                      // Environmental Friendliness of machines produced by K-firm
RowVector        A1p_en;                                     // Energy Efficiency of K-firm production technique
RowVector        A1p_ef;                                     // Environmental Friendliness of K-firm production technique
Matrix           A_en;                                       // Matrix containing energy efficiencies of existing machine tools
Matrix           A_ef;                                       // Matrix containing environmental friendliness of existing machine tools			
Matrix           A;					                         // Matrix containing productivities of existing machine tools
Matrix           C;					                         // Matrix containing cost of existing machine tools
Matrix           C_secondhand;                               // Matrix containing cost of machine tools available on secondhand market
RowVector        c1;				                         // Production cost of K-firms
Matrix           f1;                                         // Market share of K-firms
RowVector        Q1;                                         // Quantity produced by K-firms
RowVector        Td;                                         // Technological distance 
Matrix           Match;                                      // Matrix matching K-firms to customers
RowVector        S1;                                         // Revenues of K-firms
RowVector        S1_pre;                                     // Revenues of K-firms pre-shock
RowVector        S1_post;                                    // Revenues of K-firms post-shock
RowVector        Sales1;                                     // Temporary storage for revenues of K-firms
double           A1top;			                             // Maximum productivity of machines produced
double           A1ptop;			                         // Maximum productivity of K-firm production process
RowVector        A1f;				                         // Productivity of machines produced by foreign firms
RowVector        A1pf;                                       // Productivity of production process used by foreign firms
double           A1_en_top;                                  // Maximum energy efficiency of machines produced
double           A1p_en_top;                                 // Maximum energy efficiency of K-firm production process
double           A1_ef_top;                                  // Maximum environmental friendliness of machines produced
double           A1p_ef_top;                                 // Maximum environmental friendliness of K-firm production process
double           A2_ef_avg;                                  // Average environmental friendliness of C-firms machines
double           A1_ef_avg;                                  // Average environmental friendliness of K-firm production process
RowVector        A1_ci;                                      // Carbon intensity of machine produced by each K-firm
double           A1_ci_avg;                                  // Average carbon intensity of machines produced by K-firms (not weighted)
Matrix           RD;                                         // K-firms' R&D expenditures
RowVector        Inn;				                         // Indicates whether K-firm innovates
RowVector        Imm;                                        // Indicates whether K-firm imitates
RowVector        A1inn;			                             // Productivity of innovated machine
RowVector        A1pinn;			                         // Productivity of innovated production process
RowVector        A1imm;			                             // Productivity of imitated machine
RowVector        A1pimm;                                     // Productivity of imitated production process
RowVector        EE_inn;                                     // Energy efficiency of innovated machine
RowVector        EEp_inn;                                    // Energy efficiency of innovated production process
RowVector        EE_imm;                                     // Energy efficiency of imitated machine
RowVector        EEp_imm;                                    // Energy efficiency of imitated production process
RowVector        EF_inn;                                     // Environmental friendliness of innovated machine
RowVector        EFp_inn;                                    // Environmental friendliness of innovated production process
RowVector        EF_imm;                                     // Environmental friendliness of imitated machine
RowVector        EFp_imm;                                    // Environmental friendliness of imitated production process
RowVector        nclient;                                    // K-firms' number of clients
RowVector        RDin;				                         // Part of R&D expenditure devoted to innovation
RowVector        RDim;	                                     // Part of R&D expenditure devoted to imitation
RowVector        exiting_1;	                                 // Indicates whether K-firm is exiting
RowVector        exiting_1_payments;                         // Indicates whether K-firm is exiting due to inability to make a payment
RowVector        D1_en;		                                 // Energy demand of K-firms
RowVector        D1;                                         // Amount of machines demanded from K-firms
RowVector        Ld1_prod_wr;                                // Labour demand of Workers of each K-firm (only for production)
unordered_map<string, RowVector> Ld1_prod_i_mh;              // Demand for each household class (key) of each K-firm (RowVector index) (only for production)
unordered_map<string, RowVector> Ld1_i_mh;                   // Demand for each household class (key) of each K-firm (RowVector index) (production + R&D)
unordered_map<string, double> LD1_mh;                        // Demand for each household class (key) of all K-firms (production + R&D)
unordered_map<string, RowVector> Ld1_rd_i_mh;                // Labour demand for R&D of each K-firm for each household class
unordered_map<string, double> LD1_rd_mh;                     // Labour demand for R&D of each household class (key) of all K-firms
RowVector        Emiss1;                                     // Emissions of K-firms
Matrix           Pi1;	                                     // Profit of K-firms
RowVector        Anew;                                       // Productivity of new machines after exogenous technological change
RowVector        ee1;                                        // K-firm from which newly entering one will copy
RowVector        Attractivess_1_i;                           // Attractiveness of machines sold by each K-firm
RowVector        EnvSubsidiesPerMachine_1_i;                 // Subsidies paid for each machine sold, for each K-firm
double           EnvSubsidiesPerMachines_Max;                // Maximum environemtal subsidy paid (if flag_environmental_subsidies_C_firms==1)
double           A1_ci_min;                                  // Minimum carbon intensity of machines produced by K-firms (pre-innovation)  (if flag_environmental_subsidies_C_firms==1)
double           A1_ci_max;                                  // Maximum carbon intensity of machines produced by K-firms (pre-innovation)  (if flag_environmental_subsidies_C_firms==2)



//Banks
RowVector        NbClient_1;                                 // Banks' number of K-firm clients
RowVector        NbClient_2;                                 // Banks' number of C-firm clients
RowVector        NL_1;                                       // Banks' initial number of K-firm clients
RowVector        NL_2;                                       // Banks' initial number of C-firm clients
Matrix           fB;                                         // Market Share of banks
double           r_depo;                                     // Interest rate on bank deposits
RowVector        r_deb;                                      // Base loan rate of banks
RowVector        capitalAdequacyRatio;                       // Banks' capital adequacy ratio	
RowVector        riskWeightedAssets;                         // Banks' risk weighted assets
RowVector        Bank_active;                                // Indicates whether bank is active or has been bought
RowVector        bankmarkup;                                 // Mark-up of banks over monetary policy rate
Matrix           BankProfits;                                // Profit of banks
RowVector        BankProfits_temp;                           // Temporary storage for banks' profit
RowVector        BankCredit;                                 // Remaining credit supply 
RowVector        BaselBankCredit;                            // Credit supply given by regulatory ratio
RowVector        buffer;                                     // Used for calculating maximum credit according to regulatory ratio                  
RowVector        bonds_dem;                                  // Banks' demand for government bonds
double           bonds_dem_tot;                              // Total demand for government bonds from banks
Matrix           DebtServiceToSales2_bank;                   // Matrix storing debt service to sales ratio for clients of each bank
Matrix           DS2_rating;                                 // Matrix ranking clients of each bank by debt service to sales ratio
RowVector        Outflows;                                   // Transactions implying outflows of reserves
RowVector        Inflows;                                    // Transactions implying inflows of reserves
RowVector        DepositShare_e;                             // Banks' share in deposit "market"
RowVector        DepositShare_h;                             // Banks' share in deposit "market"
RowVector        LoanShare_e;                                // Banks' share in energy loan "market"
RowVector        baddebt_b;                                  // Bad debt held by banks
RowVector        bonds_purchased;                            // Gov. bonds purchases in the current period
RowVector        BankEquity_temp;                            // Temporary storage for bank net worth      
RowVector        Bailout_b;                                  // Government bailouts received by banks
RowVector        LossAbsorbed;                               // Losses absorbed when banks buy failing banks
RowVector        capitalRecovered;                           // Nominal value of second-hand capital recovered from each failing C-firms by households overall      
RowVector        capitalRecovered2;                          // Temporary storage for nominal value of second-hand capital recovered from failing firms
double           capitalRecoveredTot;                        // Nominal value of second-hand capital recovered from all failing C-firms by households overall
unordered_map<string, double> capitalRecoveredTot_mh;        // Nominal value of second-hand capital recovered from all failing C-firms by each household class 
RowVector        capitalRecoveredShare;                      // Share of nominal value of second-hand capital recovered from failing firms
RowVector        ReserveBalance;                             // Balance between inflows and outflows of reserves
RowVector        LossEntry_b;                                // Losses from firm entry taken by banks
RowVector        ShareBonds;                                 // Share of gov. bonds held by each bank
RowVector        ShareReserves;                              // Share of reserves held by each bank
RowVector        ShareAdvances;                              // Share of advances held by each bank
RowVector        Adjustment;                                 // Adjustment term to preserve stock-flow consistency

//Government
double           r_bonds;                                    // Interest rate on gov. bonds
double           Bailout;                                    // Cost of bank bailouts
double           G;                                          // Government spending on unemployment benefits in total 
double           Deficit;                                    // Public deficit
double           PSBR;                                       // Public sector borrowing requirement (includes bond repayments to be financed)
double           NewBonds;                                   // Newly issued bonds
double           EntryCosts;                                 // Firm entry costs paid by government
double           BankTransfer;                               // Transfers from gov. to banks for firm entry
double           govTransfers;                               // Government transfers to households overall
unordered_map<string, double> govTransfers_mh;               // Government transfers to each household class
double           Exp_tot_g;                                  // Government total expenditure
double           Cons_g;                                     // Government C-goods consumption demand
double           Exp_en_g;                                   // Government expenditure in energy
double           D_en_g;                                     // Government energy demand
double           energy_expenditure_sh_g;                    // Energy expenditure share of government

//Central Bank
double           r_cbreserves;                               // Interest rate on CB reserves
double           r_a;                                        // Annual policy rate (used in quarterly calibration)
double           r;                                          // Monetary policy rate
RowVector        ProfitCB(2);                                // Profits of the central bank
double           Adjustment_cb;                              // Adjustment term to preserve stock-flow consistency
double           d_cpi_target_a;                             // Annualised inflation target
double           inflation_a;                                // Annualised core inflation rate (excluding energy price)
double           inflation_a_new;                            // Annualised headline inflation rate (including energy price)

//Energy
RowVector        A_de;                                       // Thermal efficiency of dirty energy plants
RowVector        EM_de;                                      // Emissions per unit of energy produced of dirty plants
double           pf;                                         // Price of fossil fuel
double           mi_en;                                      // Mark-up of energy producer
RowVector        G_de;                                       // Number of dirty energy plants of different vintages
RowVector        G_ge;                                       // number of green energy plants of different vintages
RowVector        G_ge_n;                                     // Book value of green energy plants of different vintages
RowVector        CF_ge;                                      // Cost of expanding green energy capacity (values at all t)
RowVector        c_en(2);                                    // Price of energy
double           c_en_h;                                     // Price of energy for households
double           D1_en_TOT;                                  // Total energy demand from K-firms
double           D2_en_TOT;                                  // Total energy demand from C-firms
RowVector        D_en_TOT(2);                                // Total energy demand (inlcuding also households')
double           K_ge;                                       // Total capacity of green energy
double           K_de;                                       // Total capacity of dirty energy
double           K_ge_target_perc;                           // Target green energy share in supply when flag_energy_exp = 3
double           K_gelag;                                    // Lagged value of green capacity 
double           K_delag;                                    // Lagged value of dirty capacity
double           Q_ge;                                       // Quantity of green energy produced
double           Q_de;                                       // Quantity of dirty energy produced
RowVector        C_de;                                       // Unit production cost of dirty energy plants
double           EI_en;                                      // Investment in productive capacity of energy sector
double           EI_en_de;                                   // Investment in dirty energy capacity
double           EI_en_ge;                                   // Investment in green energy capacity
double           IC_en;                                      // Total cost of green energy investment
RowVector        IC_en_quota;                                // Investment cost for green energy paid in each period (ammortisation)
double           LDexp_en;                                   // Labour demand for green energy investment
double           LDprod_en;                                  // Labour demand for energy production
double           LDmaint_en;                                 // Labour demand for maintenance of energy capacity
double           PC_en;                                      // Production cost of energy sector
double           c_infra;                                    // Unit cost of the inframarginal dirty plant used
double           share_de;                                   // Share of dirty energy in total energy capacity
double           Rev_en;                                     // Revenue of energy sector
double           RD_en_de;                                   // R&D expenditure on dirty technology
double           RD_en_ge;                                   // R&D expenditure on green technology
unordered_map<string, double> LDen_rd_de_mh;                 // Labour demand for R&D in dirty technology
unordered_map<string, double> LDen_rd_ge_mh;                 // Labour demand for R&D in green technology
unordered_map<string, double> LDen_exp_mh;                   // Labour demand for green energy investment for each household class
double           Inn_en_ge;                                  // Indicates whether innovation has taken place in green energy
double           Inn_en_de;                                  // Indicates whether innovation has taken place in dirty energy
double           A_de_inn;                                   // Thermal efficiency of new dirty technology
double           EM_de_inn;                                  // Emission intensity of new dirty technology
double           CF_ge_inn;                                  // Installation cost of new green technology
RowVector        ProfitEnergy(2);                               // Profit of the energy sector
double           G_de_0;                                     // Initial capacity of dirty energy
double           G_ge_0;                                     // Initial capacity of green energy
double           G_ge_n_0;                                   // Initial nominal value of green capacity
double           Loan_interest_e(2);                         // Interest payments charged to the energy sector 
RowVector        CreditDemand_e(2);                          // Aggregate credit demand of the energy sector 
double           DebtService_e;                              // Aggregate debt service to be paid by the energy sector
double           DeafaultedDebtRecovered_e(2);               // Agggregate recovered bad debt of the energy sector 
double           DebtRemittances_e;                          // Loan repayments made by energy sector
double           DebtWrittenOff_e;                           // Energy sector loans written off due to default
double           DefaultedDeposits_e;                        // Remaining deposits of energy sector after default
double           ShareDefaultedDebt_e;                       // Share of energy sector loans in default
double           BadDebt_e;                                  // Total energy sector debt in default
double           Loans_preDefault_e;                         // Stock of energy sector loans prior to defaults

//Climate
RowVector        Tmixed(2);                                  // Temperature in the mixed layer
RowVector        Emiss_yearly_calib(2);                      // Emissions calibrated w.r.t. year 2010
double           g_rate_em_y;                                // Annual growth rate of emissions in EU (compared to initial value)
RowVector        Emiss_TOT;                                  // Total emissions (in EU)
RowVector        Emiss_yearly(2);                            // Annual emissions (EU + global)
double           Emiss1_TOT;                                 // Total emissions K-firms
double           Emiss2_TOT;                                 // Total emissions C-firms
double           Emiss_en;                                   // Total emissions Energy sector
double           NPP;                                        // Net primary production
double		     Cum_emissions;	                             // Cumulative emissions
RowVector        shocks_machprod;                            // Shocks to productivity of machines
RowVector        shocks_techprod;                            // Shocks to productivity of K-firm production processes
RowVector        shocks_machprod_max;                        // Max size of shocks to productivity of machines
RowVector        shocks_techprod_max;                        // Max size of shocks to productivity of K-firm production processes
RowVector        shocks_encapstock_de;                       // Shocks to capacity of dirty energy
RowVector        shocks_encapstock_ge;                       // Shocks to capacity of green energy
RowVector        shocks_capstock;                            // Shocks to stock of machines
RowVector        shocks_invent;                              // Shocks to inventories
RowVector        shocks_rd;                                  // Shocks to R&D effectiveness
RowVector        shocks_labprod1;                            // Shocks to K-firm labour productivity
RowVector        shocks_labprod2;                            // Shocks to C-firm labour productivity
RowVector        shocks_eneff1;                              // Shocks to K-firm energy efficiency
RowVector        shocks_eneff2;                              // Shocks to C-firm energy efficiency
RowVector        shocks_output1;                             // Shocks to K-firm output
RowVector        shocks_output2;                             // Shocks to C-firm output
double           shock_pop;                                  // Shock to population of households _mh need to manage to individual classes
double           shock_cons;                                 // Shock to aggregate consumption demand
double           shock_nord;                                 // Shock value given by DICE-type damage function
RowVector        X_a;                                        // Location parameters of beta distribution for climate shocks
RowVector        X_b;                                        // Scale parameters of beta distribution for climate shocks
RowVector        Loss_Capital;                               // Nominal value of capital stock destroyed by climate shock (used for SFC check)
RowVector        Loss_Capital_mol;                           // Nominal value of capital stock destroyed by climate shock (used to adjust gross profit margin)
RowVector        Loss_Inventories;                           // Nominal value of inventories destroyed by climate shock
double           t_CO2;                                      // Carbon tax on firms
double           t_CO2_en;                                   // Carbon tax on energy sector
double           Emiss_gauge;                                // Reference value for emissions
RowVector        Cat(2);                                     // Atmospheric carbon in GtC
double           humrelease;                                 // Carbon released from decaying humus 
RowVector        hum(2);                                     // Carbon stored in humus (GtC)
double           biorelease;                                 // Carbon released from decaying biomass
RowVector        biom(2);                                    // Carbon stored in biomass (GtC)
double           Cat1;                                       // Temporary guess for Cat after biosphere and emissions
double           dCat1;                                      // Carbon added to atmosphere and ocean 
RowVector        fluxC;                                      // Carbon flux between ocean layers
Matrix           Con;                                        // Carbon content in the ndep ocean layers in GtC
Matrix           Hon;                                        // Heat content of the ocean layers  (Joule per m^2 of ocean surface)
Matrix           Ton;                                        // Temperature of the ocean layers 
double           Con1;                                       // Temporary guess for carbon in upper ocean layer (GtC)
double           Ctot1;                                      // Temporary guess for carbon in upper ocean layer + atmosphere (GtC)
RowVector        Cax;                                        // Cat used during iteration for equlibrating atmospheric and mixed-layer carbon
RowVector        Caxx;                                       // Cat used during iteration for equlibrating atmospheric and mixed-layer carbon
RowVector        Cay;                                        // Residual used during iteration
RowVector        Cayy;                                       // Residual used during iteration 
RowVector        Caa;                                        // Estimated slope duing the iteration 
double           FCO2;                                       // Radiative forcing from CO2
double           Fin;                                        // Radiative forcing input due to all greenhouse gases
double           Fout;                                       // Additional out-radiation due to global warming
RowVector        fluxH;                                      // Heat flux between ocean layers
double           Emiss_global;                               // Global GHG emissions out of EU (used for EU calibration)

//Others; Aggregate; Reporting etc.
RowVector        Am(2);                 	                 // Mean productivity across K and C-firms	
double           Am_a;                     	                 // Mean productivity across K and C-firms (alternative definition)
double           Am2;                                        // Mean productivity across C-firms
double           Am1;                                        // Mean productivity across K-firms
RowVector        ftot(3);                                    // Sum of C-firm market shares
RowVector        Em2(2);                                     // Mean competitiveness of C-firms
RowVector        cpi(5);                                     // Consumer price index (level) excluding energy price (for core inflation)
RowVector        cpi_new(5);                                 // Consumer price index including energy price (for headline inflation)
double           cpi_energy;                                 // Consumer price index (relative) for energy
double           cpi_goods;                                  // Consumer price index (relative) for goods            
unordered_map<string, double> cpi_new_mh;                    // Consumer price index for each households class based on nominal expenditure shares of goods and energy
double           kpi;                                        // Price index of capital goods
RowVector        Am_en(2);                                   // Mean effective energy efficiency across K and C-firms weighted on firms' energy demand
double           Am_en_2;                                    // Mean effective energy efficiency across C-firms weighted on firms' energy demand
double           Am_en_1;                                    // Mean effective energy efficiency across K-firms weighted on firms' energy demand
double           A2_en_avg;                                  // Mean effective energy efficiency across C-firms (not weighted)
RowVector        A2_ci;                                      // Carbon intensity of each C-firms (effective)
double           A2_ci_avg;                                  // Mean carbon intensity across C-firms (effective)
double           Tdtot;                                      // Sum of technological distances
unordered_map<string, double> LDen_tot_mh;                   // Total labour demand of energy sector of each household class
double           LD1_wr;					                 // Total labour demand of Workers for production from K-firms
double           LD2_wr;                                     // Total labour demand of Workers for production from C-firms
double           LSe_wr;                                     // Remaining labour supply of Workers
double           LD;                                         // Total labour demand overall
unordered_map<string, double> LD_mh;                         // Total labour demand of each houseold class
double           LD2;                                        // Total labour demand for production from firms
double           Pitot1;				                     // Total profit of K-firms
double           Pitot2;                                     // Total profit of C-firms
double           dNtot;                                      // Total change in real value of inventories
double           dNmtot;                                     // Total change in nominal value of inventories
double           ExpansionInvestment_r;                      // Total real expansion investment
double           ExpansionInvestment_n;                      // Total nominal expansion investment
double           ReplacementInvestment_r;                    // Total real replacement investment
double           ReplacementInvestment_n;                    // Total nominal replacement investment
double           Investment_r;                               // Total real investment
double           Investment_n;                               // Total nominal investment
Real             CreditDemand_all;                           // Total demand for bank credit
Real             CreditSupply_all;                           // Total supply of bank credit
double           Q2tot;                                      // Total output of C-firms
double           Q1tot;                                      // Total output of K-firms
double	         Q2dtot;                                     // Total desired output of C-firms
double           D2tot;                                      // Total demand for consumption goods
double           A_mi;                                       // Mean of log C-firm productivity
double           A1_mi;                                      // Mean of log K-firm productivity
double           A2_en_mi;                                   // Mean of log C-firm energy efficiency
double           A2_ef_mi;                                   // Mean of log C-firm environmental friendliness
double           A1_en_mi;                                   // Mean of log log K-firm energy efficiency
double           A1_ef_mi;                                   // Mean of log K-firm environmental friendliness
double           A_sd;                                       // Standard deviation of C-firm productivity
double           H1;                                         // Herfindahl index K-firms
double           H2;                                         // Herfindahl index C-firms
double           HB;                                         // Herfindahl index banks                                        
double           d_U;                                        // Change in unemployment rate
double           d_cpi;                                      // Change in consumer price index
double           d_Am;                                       // Change in mean productivity
double           dw;                                         // Change in LD1_wr
double           dw2;                                        // Change in wage
double           A2scr;                                      // Log deviation of C-firm productivity from mean
double           A1scr;                                      // Log deviation of K-firm productivity from mean
Matrix           S1_temp;                                    // Temporary storage for K-firm sales
Matrix           S2_temp;                                    // Temporary storage for C-firm sales
double           Utilisation;                                // C-firms' aggregate capacity utilisation
double           counter_bankfailure;                        // Number of failing banks
    // GDP
double           GDP_rg;                                     // Growth rate real GDP
double           GDP_ng;                                     // Growth rate nominal GDP
RowVector        GDP_r(2);                                   // Real GDP
RowVector        GDP_n(2);                                   // Nominal GDP
double           GDP_n_exp;                                  // Nominal GDP calculated through expenditure approach     
double           GDP_n_inc;                                  // Nominal GDP calculated through income approach
double           GDP_n_prod;                                 // Nominal GDP calculated through production approach
RowVector        GDP_r_new(2);                                   // Real GDP with new approach
RowVector        GDP_n_new(2);                                   // Nominal GDP with new approach
double           GDP_rg_new;                                     // Growth rate real GDP with new approach
double           GDP_ng_new;                                     // Growth rate nominal GDP with new approach
double           Exports;                                   //Total exports for expenditure-approach GDP
double           Imports;                                   //Total imports for expenditure-approach GDP
double           VA_1;                                      //Value added of K-sector for production-approach GDP
double           VA_2;                                      //Value added of C-sector for production-approach GDP
double           VA_en;                                      //Value added of energy sector for production-approach GDP
double           VA_b;                                      //Value added of banks for production-approach GDP
double           PaymentsToLabour;                          //Total payments to labour for income-approach GDP
double           PaymentsToGovernment;                      //Total payments to government for income-approach GDP
double           PaymentsToCapital;                         //Total payments to capital for income-approach GDP

//Ouput
std::vector<std::pair<std::string, std::function<double()>>> output_pairs; //Vector with variable_name - function pairs for printing general output
std::ostringstream output_string;                           //Output string stream with variables to print for general output
    //Filenames
char nomefile1[64];                                          // File "out" (inv_output1)
char nomefile2[64];                                          // File "ymc" (inv_ymc)
char nomefile3[64];                                          // File "resultsexp" (inv_res)
char nomefile4[64];                                          // File "A1" (inv_prod1)
char nomefile5[64];                                          // File "A2" (inv_prod2)
char nomefile6[64];                                          // File "A1_all" (inv_prodall1)  
char nomefile7[64];                                          // File "A2_all" (inv_prodall2)
char nomefile8[64];                                          // File "A1all_en" (inv_prodall1_en)
char nomefile9[64];                                          // File "A2all_en" (inv_prodall2_en)
char nomefile10[64];                                         // File "A1all_ef" (inv_prodall1_ef)
char nomefile11[64];                                         // File "A2all_ef" (inv_prodall2_ef)
char nomefile12[64];                                         // File "NW1all" (inv_nwall1)
char nomefile13[64];                                         // File "NW2all" (inv_nwall2)
char nomefile14[64];                                         // File "NWBall" (inv_nwall3)
char nomefile15[64];                                         // File "Deb2all" (inv_deball2)
char nomefile16[64];                                         // File "validation1" (inv_val1)
char nomefile17[64];                                         // File "validation2" (inv_val2)
char nomefile18[64];                                         // File "validation3" (inv_val3)
char nomefile19[64];                                         // File "validation4" (inv_val4)
char nomefile20[64];                                         // File "validation5" (inv_val5)
char nomefile21[64];                                         // File "validation6" (inv_val6)
char nomefile22[64];                                         // File "validation7" (inv_val7)
char nomefile23[64];                                         // File "validation8" (inv_val8)
char nomefile24[64];                                         // File "validation9" (inv_val9)
char nomefile25[64];                                         // File "validation10" (inv_val10)
char nomefile26[64];                                         // File "validation10" (inv_val10)
char nomefile27[64];                                         // File "validation10" (inv_val10)
char nomefile28[64];                                         // File "shockpars" (inv_shockpars)
char errorfilename[150];                                     // Name of error file

#endif