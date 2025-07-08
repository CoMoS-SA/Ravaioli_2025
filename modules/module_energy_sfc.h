#ifndef MODULE_ENERGY_H
#define MODULE_ENERGY_H

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

#include <string.h>
#include <sstream>
#include <map>
#include <vector>
#include <unordered_map>

// Include Newmat Libraries and random number generators
#include "../newmat10/include.h"
#include "../newmat10/newmat.h"
#include "../newmat10/newmatio.h"
#include "../auxiliary/ran1.h"
#include "../auxiliary/betadev.h"
#include "../auxiliary/bnldev.h"

// Include functions from other modules
#include "../dsk_sfc_functions.h"

// -- Functions -- //
void EN_DEM_TOT(void);                                  // Aggregates energy demand coming from firms
void ENERGY_INV_PROD(void);                         // Energy production and investment in productive capacity
void ENERGY_RandD(void);                            // Energy sector R&D
void ENERGY_LABOUR(void);                           // Calculates energy sector labour demand and wages to pay                             
void EMISS_IND(void);                               // Calculates aggregate emissions from firms

//-- Flags --//
extern int              flag_energy_exp;
extern int              flag_share_END;
extern int              flag_energyshocks;
extern int              flag_endogenous_exp_quota;

//-- Pars --//
extern int              N1;
extern int              N2;
extern int              payback_en;
extern double           exp_quota;
extern double           t_CO2_en;
extern double           share_de_0;
extern double           share_RD_en;
extern double           o1_en;
extern long             *p_seed;
extern double           b_a11;
extern double           b_b11;
extern double           uu1_en;
extern double           uu2_en;
extern int              life_plant;
extern int              T;
extern double           exp_quota_param;
extern double           profit_share_energy_inv; 
extern double           rescale_ene_price; 
extern unordered_map<string, double> ld_ratios_mh;
extern vector<string>   classes_mh;
extern unordered_map<string, double> Transfer_shock_sh_mh;
extern int              t_regime_shifts;
extern double           ratio_mi_en_shock;
extern double           K_ge_END_perc;
extern int              t_length_energy_transition;
extern int              t_regime_shift;
extern double           renew_impact_on_p_e;

//-- Inits --//
extern double           K_ge0_perc;

// -- Vars -- //
extern RowVector        D1_en;		
extern RowVector        D2_en;			
extern double           D1_en_TOT;              
extern double           D2_en_TOT;              
extern RowVector        D_en_TOT;
extern int              i;
extern int              j;
extern double           tolerance;
extern RowVector        A1p_ef;
extern RowVector        A1p_en;
extern RowVector        Q1;
extern RowVector        Q2;
extern RowVector        A2e_ef;
extern RowVector        A2e_en;
extern RowVector        Emiss1;
extern RowVector        Emiss2;
extern double           Emiss1_TOT;            
extern double           Emiss2_TOT;
extern int              t;
extern int              tt;
extern double           K_ge;                   
extern double           K_de;                  
extern double           Q_ge;                   
extern double           Q_de; 
extern RowVector        C_de;                   
extern RowVector        G_de;  
extern RowVector        G_ge;
extern RowVector        G_ge_n;
extern double           pf;
extern double           mi_en;
extern RowVector        c_en;
extern RowVector        A_de;
extern double           EI_en;
extern double           c_de_min;
extern double           cf_min_ge;
extern RowVector        CF_ge;
extern double           EI_en_de;
extern double           EI_en_ge;
extern double           IC_en;
extern RowVector        IC_en_quota;
extern unordered_map<string, double> LDen_exp_mh;
extern double           LDprod_en; 
extern double           LDmaint_en;
extern RowVector        w_tot_for_1_wr_mh;
extern unordered_map<string, RowVector> w_mh;
extern double           PC_en;
extern double           Emiss_en;
extern RowVector        EM_de;
extern RowVector        G_de_temp;
extern double           Q_de_temp;
extern int              idmin;
extern double           c_infra;
extern double           share_de;
extern double           Rev_en;
extern double           RD_en_de;
extern double           RD_en_ge;
extern unordered_map<string, double> LDen_rd_de_mh;
extern unordered_map<string, double> LDen_rd_ge_mh;
extern unordered_map<string, double> LDen_tot_mh;
extern double           parber_en_de;
extern double           parber_en_ge;
extern double           rnd;
extern double           Inn_en_ge;
extern double           Inn_en_de;
extern double           A_de_inn;
extern double           EM_de_inn;
extern double           CF_ge_inn;
extern RowVector        Emiss_TOT;
extern RowVector        CapitalStock_e;
extern double           Wages_en;                                         
extern RowVector        cpi;
extern unordered_map<string, double> Wages_en_mh;
extern double           FuelCost;
extern double           mi_en_preshock;
extern double           pf_preshock;
extern double           mi_en_shock;
extern double           c_en_preshock;
extern double           pf_shock;
extern double           c_infra_t;
extern double           K_gelag;
extern double           K_delag;
extern double           K_ge_target_perc;
extern double           G_de_0; 
extern double           G_ge_0;
extern double           G_ge_n_0;
extern double           D_en_h;
extern string           tech;
extern double           old_capitalStock;
extern double           LDexp_en;
extern double           D_en_g;


#endif