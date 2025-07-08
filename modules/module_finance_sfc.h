#ifndef MODULE_FINANCE_H
#define MODULE_FINANCE_H

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
void TOTCREDIT(void);                                                       // Determines maximum amount of credit banks will extend
void LOANRATES(void);                                                       // Determines loan rates charged to individual borrowers
void BANKING(void);                                                         // Determines bank profits; Banks receive second-hand capital from failing firms
void BAILOUT(void);                                                         // Failing banks are bailed out by government or bought by other banks depending on scenario
void SETTLEMENT(void);                                                      // Settlement of interbank transactions; granting/repayment of CB Advances

//-- Flags --//
extern int              flagbailout;
extern int              flag_rate_setting_markup;
extern int              flag_rate_setting_loans;

//-- Pars --//
extern double           credit_multiplier;
extern double           beta_basel;
extern double           floor_default_probability;
extern double           upsilon;
extern double           lambdaB1;
extern double           lambdaB2;
extern double           riskWeightLoans;	
extern double           riskWeightGovBonds;	
extern double           capitalAdequacyRatioTarget;
extern int              NB;
extern int              N2;
extern int              N1;
extern double           varphi;
extern double           k_const;
extern double           k_const2;
extern double           db;
extern double           aliqb;
extern long             *p_seed;
extern double           b1sup;
extern double           b1inf;
extern double           b2sup;
extern double           b2inf;
extern double           agemax;
extern double           FirmDefaultProbability_init;
extern vector<string>   classes_mh;
extern unordered_map<string, double> Ownership_sh_b_mh;
extern double           bonuses_share;

// -- Vars -- //
extern int              i;
extern int              j;
extern int              t;
extern int              tt;
extern int              t0;
extern double           tolerance;
extern RowVector        BankCredit;
extern RowVector        BaselBankCredit;
extern Matrix           NW_b;
extern Matrix           NW_2;
extern RowVector        buffer;
extern RowVector        capitalAdequacyRatio;	
extern RowVector        riskWeightedAssets;
extern RowVector        bonds_dem;
extern double           bonds_dem_tot;
extern Matrix           S2;
extern Matrix           Deposits_2;
extern Matrix           Loans_2;
extern RowVector        CreditDemand;
extern RowVector        DebtServiceToSales2; 
extern Real             DS2_min;   
extern Matrix           DebtServiceToSales2_bank;
extern Matrix           BankMatch_2;
extern Matrix           BankMatch_1;
extern RowVector        DebtServiceToSales2_temp; 
extern int              DS2_min_index; 
extern Matrix           DS2_rating;
extern Matrix           DebtService_2;
extern RowVector        k;
extern RowVector        NL_2;
extern RowVector        NL_1;
extern RowVector        r_deb; 
extern RowVector        r_deb_h;
extern RowVector        FirmDefaultProbability;
extern RowVector        LoanInterest;
extern RowVector        InterestDeposits;
extern RowVector        baddebt_b;  
extern double           r_bonds;
extern Matrix           GB_b;
extern double           r;
extern double           r_cbreserves;
extern Matrix           Reserves_b;
extern Matrix           Advances_b;
extern RowVector        Dividends_b_i;
extern double           Dividends_b;
extern unordered_map<string, double> Dividends_b_mh;
extern unordered_map<string, RowVector> Dividends_mh;
extern RowVector        Bonuses_b_i;
extern double           Bonuses_b;
extern unordered_map<string, RowVector> Bonuses_mh;
extern unordered_map<string, double> Bonuses_b_mh;
extern double           Bonuses_h;
extern double           Taxes_g;
extern RowVector        Taxes_b;
extern Matrix           Deposits_b;
extern Matrix           Deposits_hb;
extern Matrix           Deposits_eb;
extern RowVector        Deposits_h;
extern unordered_map<string, RowVector> Deposits_mh;
extern RowVector        Outflows;
extern RowVector        Inflows;
extern Matrix           BankProfits;
extern RowVector        BankProfits_temp;
extern RowVector        Bank_active;
extern double           Bailout;
extern RowVector        BankEquity_temp; 
extern int              maxbank;
extern double           max_equity;
extern double           multip_bailout;
extern RowVector        Bailout_b;
extern RowVector        LossAbsorbed;
extern Matrix           Loans_b;
extern double           min_equity;
extern RowVector        BankingSupplier_2;
extern RowVector        BankingSupplier_1;
extern RowVector        exiting_2;
extern int              receivingBank;
extern RowVector        capitalRecovered;
extern double           n_mach_exit;
extern RowVector        ReserveBalance;
extern RowVector        Advances;
extern RowVector        Reserves;
extern RowVector        Dividends_h;
extern double           Bonuses_h;
extern unordered_map<string, RowVector> Bonuses_mh;
extern double           InterestReserves;
extern double           InterestAdvances;
extern RowVector        InterestReserves_b;
extern RowVector        InterestAdvances_b;
extern RowVector        ProfitCB;
extern double           InterestBonds_cb;
extern RowVector        LossEntry_b;
extern RowVector        DepositShare_h;
extern RowVector        DepositShare_e;
extern double           counter_bankfailure;

#endif