#ifndef MODULE_MACRO_H
#define MODULE_MACRO_H

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <fenv.h>

//#include <string>
#include <string.h>
#include <sstream>

#include <unordered_map>

// Include Newmat Libraries and random number generators
#include "../newmat10/include.h"
#include "../newmat10/newmat.h"
#include "../newmat10/newmatio.h"
#include "../auxiliary/ran1.h"

// Include functions from other modules
#include "../dsk_sfc_functions.h"


// -- Functions -- //
void LABOR(void);                                                     // Allocates labour supply; scales down production of firms if labour supply is insufficient
void MACRO(void);                                                     // Calculates macroeconomic aggregates
void GDP(void);                                                       // Calculates GDP with different accounting and check consistency
void WAGE(void);                                                      // Determines change in wage rate
void GOV_BUDGET(void);                                                // Determines unemployment benefits & government deficit; implements bond market
void TAYLOR(void);                                                    // Determines change in monetary policy rate

//-- Pars --//
extern int              N1;
extern int              N2;
extern double           N1r;
extern double           N2r;
extern int              NB;
extern double           dim_mach;
extern double           wu;
extern double           bonds_share;
extern double           aliqb;
extern double           taylor1;
extern double           taylor2;
extern double           bankmarkdown;
extern double           bankmarkup_init;
extern double           lambdaB;
extern double           CapitalAdequacyRatioTarget;
extern double           centralbankmarkdown;
extern double           bondsmarkdown;
extern double           ustar;
extern double           d_cpi_target;
extern double           mdw;
extern double           psi1;
extern double           psi2;
extern double           psi3;
extern double           w_min;
extern int              reducedoutput;
extern double           kappa;
extern double           taylor;
extern double           varphi;
extern double           g_ls;
extern long             *p_seed;
extern unordered_map<string, double> ld_ratios_mh;
extern unordered_map<string, double> w_ratios_mh;
extern vector<string>   classes_mh;
extern double           tolerance;
extern int              t_regime_shifts;

// -- Vars -- //
extern ofstream         Errors;
extern int              i;
extern int              j;
extern int              t;
extern double           LD1_wr;					
extern double           LD2_wr;
extern unordered_map<string, double> LD1_rd_mh;
extern unordered_map<string, double> LDen_tot_mh;
extern RowVector        Ld1_prod_wr;
extern RowVector        Ld2_wr;
extern unordered_map<string, double> LD1_mh;
extern unordered_map<string, double> LD2_mh;
extern double           LSe_wr;
extern double           LS;
extern unordered_map<string, double>  LS_mh;
extern double           LD;
extern unordered_map<string, double> LD_mh; 
extern double           LD2;
extern RowVector        Q2;
extern RowVector        Q1;
extern double           Qpast;
extern Matrix           Match;
extern RowVector        I;					
extern Matrix           EI;		
extern RowVector        SI;
extern RowVector        GB;
extern RowVector        GB_cb;
extern double           Deficit;
extern double           Taxes_g;
extern double           Bailout;
extern double           r_bonds;
extern double           r_cbreserves;
extern double           PSBR;
extern double           NewBonds;
extern Matrix           GB_b;
extern RowVector        Bond_share;
extern RowVector        Inflows;
extern double           InterestBonds;
extern double           InterestBonds_cb;
extern RowVector        InterestBonds_b;
extern double           BondRepayments_cb;
extern RowVector        BondRepayments_b;
extern RowVector        bonds_dem;
extern double           bonds_dem_tot;
extern Matrix           BankProfits;
extern RowVector        BankProfits_temp;
extern RowVector        DepositShare_h;
extern Matrix           Deposits_b;
extern Matrix           Deposits_hb;
extern double           r;
extern double           r_depo;
extern RowVector        r_deb;
extern double           r_base;
extern RowVector        bankmarkup;
extern RowVector        U;
extern unordered_map<string, double> U_mh;
extern double           d_cpi;
extern RowVector        Outflows;
extern RowVector        bonds_purchased;
extern RowVector        Am;
extern double           Am_a;
extern double           Am1;
extern double           Am2;
extern double           ExpansionInvestment_r;
extern double           ExpansionInvestment_n;
extern double           ReplacementInvestment_r;
extern double           ReplacementInvestment_n;
extern double           Investment_r;
extern double           Investment_n;
extern double           Consumption_r;
extern Real             CreditDemand_all;
extern Real             CreditSupply_all;
extern double           Q2tot;
extern double           Q1tot;
extern double	        Q2dtot;
extern double           D2tot;
extern double           A_mi;
extern double           A1_mi;
extern double           A2_en_mi;
extern double           A2_ef_mi;
extern double           A1_en_mi;
extern double           A1_ef_mi;
extern double           A_sd;
extern double           H1;
extern double           H2;
extern RowVector        Qd;
extern Matrix           D2;
extern RowVector        EI_n;
extern RowVector        SI_n;
extern RowVector        A2e;
extern Matrix           f2;
extern Matrix           f1;
extern RowVector        A2e_en;
extern RowVector        A2e_ef;
extern double           D2_en_TOT;
extern RowVector        Emiss2;
extern double           Emiss2_TOT;
extern double           D1_en_TOT;
extern RowVector        Emiss1;
extern double           Emiss1_TOT;
extern RowVector        A2;
extern RowVector        A2_mprod;
extern RowVector        A2_en;
extern RowVector        A2_ef;
extern RowVector        A1p;
extern RowVector        A1p_en;
extern RowVector        A1p_ef;
extern RowVector        BaselBankCredit;
extern RowVector        CreditDemand;
extern double           HB;
extern Matrix           fB;
extern double           Pitot1;
extern double           Pitot2;
extern double           d_U;
extern RowVector        cpi;
extern double           d_Am;
extern double           dw;
extern RowVector        ProfitCB;
extern double           TransferCB;
extern double           G;
extern unordered_map<string, double> Benefits_mh;
extern double           Benefits;
extern double           EntryCosts;
extern RowVector        Taxes_CO2;
extern RowVector        p1;
extern RowVector        p2;
extern double           d_cpi_target_a;
extern double           inflation_a;
extern RowVector        D2_en;
extern RowVector        D1_en;
extern RowVector        Am_en;
extern Matrix           Loans_b;
extern int              ranj;
extern RowVector        shocks_labprod1;
extern double           reduction;
extern double           Ipast;
extern RowVector        K;
extern RowVector        nclient;
extern double           r_a;
extern double           Taxes_e_shock;
extern double           Taxes_f_ff_shock;
extern double           Taxes_e_ff_shock;
extern double           Transfer_shock_f;
extern double           govTransfers;
extern unordered_map<string, double> govTransfers_mh;
extern double           redistribute_co2TaxRev; 
extern unordered_map<string, double> co2TaxRev_sh_mh;
extern unordered_map<string, double>  Transfer_shock_mh;
extern double           Transfer_shock;
extern unordered_map<string, RowVector> w_mh;
extern RowVector        w_tot_for_1_wr_mh;
extern double           D_en_h;
extern double           Am_en_2;
extern double           Am_en_1;
extern double           EnvSubsidies_2;
extern RowVector        EnvSubsidies_2_i;
    //GDP
extern RowVector        GDP_r;
extern RowVector        GDP_n;
extern RowVector        GDP_r_new;
extern RowVector        GDP_n_new;
extern double           GDP_n_exp;
extern double           GDP_n_inc;
extern double           GDP_n_prod;
extern double           GDP_rg;
extern double           GDP_ng;
extern double           GDP_rg_new;
extern double           GDP_ng_new;
extern double           Consumption;
extern double           Consumption_g;
extern double           Exp_tot_g;
extern double           FuelCost;
extern double           Wages;
extern RowVector        Wages_1_i;
extern RowVector        Wages_2_i;
extern double           Wages_en;
extern double           InterestDeposits_h;
extern RowVector        Dividends_h;
extern double           Bonuses_h;
extern RowVector        Bonuses_1_i;
extern RowVector        Bonuses_2_i;
extern double           Bonuses_e;
extern RowVector        S1;
extern Matrix           S2;
extern double           EnergyPayments;
extern RowVector        EnergyPayments_1;
extern RowVector        EnergyPayments_2;
extern Matrix           Inventories;
extern RowVector        Taxes_1;
extern RowVector        Taxes_2;            
extern RowVector        Taxes_CO2_1;
extern RowVector        Taxes_CO2_2;
extern double           Taxes_CO2_e;
extern double           Taxes_e_shock;
extern RowVector        Taxes_b;
extern RowVector        LoanInterest;
extern RowVector        LoanInterest_2;
extern double           Exports;
extern double           Imports;
extern double           VA_1;
extern double           VA_2;
extern double           VA_en;
extern double           VA_b;
extern double           PaymentsToLabour;  
extern double           PaymentsToGovernment;
extern double           PaymentsToCapital;
extern double           Expenditure_tot_h;


//-- Flags --//
extern int              flag_rate_setting_markup;
extern int              flag_energyshocks;
extern int              flag_energyshocks_MP;


#endif