#ifndef FLAGS_H
#define FLAGS_H

int flag_exogenousshocks;               // Setting = 1 activates experiment with exogenous climate shocks

int flag_clim_tech;                     // Activates key DSK (as opposed to KS) features including climate box, shocks, energy tech., etc.
                                        // 0 = KS
                                        // 1 = [BASELINE] DSK

int flag_cum_emissions;                 // Switches between C-Roads climate box and simple cumulative emissions one
                                        // = 0 [BASELINE] C-ROADS
                                        // = 1  simple cumulative emission linear relation with temp

int flag_tax_CO2;                       // Activates C02 tax 
                                        // = 0 [BASELINE] off
                                        // = 1 on and increasing with inflation
                                        // = 2 on and increasing linearly with time
                                        // = 3 on and increasing exponentially with time + inflation correction
                                        // = 4 on and increasing with nominal GDP
                                        // = 5 as 1, but introduced from regime shift + initial tax on energy set at same value as for industry (to have only 1 policy parameter)
                                        // = 6 as 3, but introduced from regime shift + initial tax on energy set at same value as for industry (to have only 1 policy parameter)

int flag_encapshocks;                   // Shocks to energy sector's productive capacity
                                        // = 0 no shock
                                        // = 1 Energy sector loses some percentage of both brown and green capacity

int flag_popshocks;                     // Shocks to the population/labour force
                                        // = 0 no shock
                                        // = 1 reduce labour force by some percentage

int flag_demandshocks;                  // Shocks to aggregate consumption demand
                                        // = 0 no shock
                                        // = 1 aggregate consumption demand reduced by x%

int flag_capshocks;                     // Shocks to C-firms' capital stocks
                                        // = 0 no shock
                                        // = 1 x% shocks to capital stocks of all firms 
                                        // = 2 x% shock to aggregate capital stock, affecting firms with uniform probability
                                        // = 3 x% shock to aggregate capital stock, with some firms having higher prob. of being affected

int flag_outputshocks;                  // Shocks to current output of C- and K-firms
                                        // = 0 no shock
                                        // = 1 x% shocks to current output of all C- and K-firms 
                                        // = 2 x% shock to aggregate output of both cons. and cap. goods, affecting firms with uniform probability
                                        // = 3 x% shock to aggregate output of both cons. and cap. goods, with some firms having higher prob. of being affected

int flag_inventshocks;                  // Shocks to inventories of C-firms
                                        // = 0 no shock
                                        // = 1 x% shocks to inventories of all C-firms
                                        // = 2 x% shock to aggregate inventories of C-firms, affecting firms with uniform probability
                                        // = 3 x% shock to aggregate inventories of C-firms, with some firms having higher prob. of being affected

int flag_RDshocks;                      // Shocks to the R&D process of K-firms
                                        // = 0 no shocks
                                        // = 1 on shape of beta distribution (alpha)
                                        // = 2 on support of beta distribution
                                        // = 3 on Bernoulli determining whether firm innovates/imitates at all
                                        // = 4 on amount of resources devoted to R&D

int flag_prodshocks1;                   // Shocks to productivity affecting the characteristics of capital vintages (roughly similar to TFP shocks in conventional models)
                                        // = 0 no shock
                                        // = 1 On labour productivity of current vintages
                                        // = 2 On energy efficiency of current vintages
                                        // = 3 Both labour producitivity and energy efficiency of current vintages
                                        // = 4 On labour productivity of all existing vintages
                                        // = 5 On energy efficiency of all existing vintages
                                        // = 6 Both labour producitivity and energy efficiency of all existing vintages

int flag_prodshocks2;                   // Shocks to productivity (not affecting characteristics of capital vintages)
                                        // = 0 no shock
                                        // = 1 On labour productivity of C-firms and K-firms
                                        // = 2 On energy efficiency of C-firms and K-firms
                                        // = 3 Both labour producitivity and energy efficiency

int flag_energyshocks;                  // Shock to energy price
                                        // = 0 no shock
                                        // = 1 transitory shock to energy price through mark-up
                                        // = 2 permanent shock to energy price through mark-up
                                        // = 3 transitory shock to energy price through fossil fuel price
                                        // = 4 permanent shock to energy price through fossil fuel price

int flag_energyshocks_MP;               // Indicates whether monetary policy reacts to inflation during energy price shock
                                        // = 0 No
                                        // = 1 (Default) Yes 

int flag_share_END;                     // Switches on endogenous share of R&D expenditures in dirty vs green energy
                                        // = 0 exogenous, given by share_de_0
                                        // = 1 endogenous, based on share in total energy sector productive capacity
                                        // = 2 endogenous, based on share in total energy sector supply

int flag_energy_exp;                    //Determines whether maximum expansion of green energy capacity per period is constrained
                                        // = 0 not constrained
                                        // = 1 [BASELINE] constrained
                                        // = 2 green energy capacity is expanded in order to keep the green share equal to the initial one
                                        // = 3 green energy capacity is expanded to reach a green share equal to K_ge_END_perc, increasing it linearly starting from t_regime_shifts to reach the target after t_length_energy_transition timesteps
                                        // = 4 as for 3, but energy price markup's indexation on wages is reduced proportionally to the increase in energy share, with proportionality defined by renew_impact_on_p_e NOT IMPLEMENTED FOR ENERGY PRICE SHOCKS
                                        // > 4 expansion is constrained but minimum investment in green to keep share equal to initial one

int flag_endogenous_exp_quota;          //Determines whether the constraint on green investment is exogenous or endogenous
                                        // = 0 exogenous (fixed value of exp_quota)
                                        // = 1 endogenous (determined by exp_quota_param)

int flag_pop_growth;                    // Switches on population growth at exogenous rate g_ls from period 200
                                        // = 0 off
                                        // = 1 [BASELINE] on

int flagbailout;                        // Switches between bailout rules for banks
                                        // = 0 [BASELINE] Banks are always bailed out by government
                                        // = 1 Failing banks are purchased by largest bank; government steps in as last resort

int flag_entry;                         // Determines what happens if households cannot finance K-firm and/or C-firm entry
                                        // = 0 The remaining entry costs are paid by government 
                                        // = 1 The remaining entry costs are booked as a loss for the banks

int flag_nonCO2_force;                  // Determines whether or not to consider non-CO2 radiative forcing in the C-Roads climate box
                                        // = 0 Non-CO2 forcing not included
                                        // = 1 [BASELINE] Non-CO2 forcing included

int flag_inventories;                   // Switches C-firm planned inventories on or off
                                        // = 0 Off (unsold output is scrapped)
                                        // = 1 On (MODEL NOT STABLE)

int flag_uniformshocks;                 // Determines whether one-off micro-level climate shocks should be uniform or heterogeneous across agents
                                        // = 0 heterogeneous shocks
                                        // = 1 uniform shocks (i.e. same size shock for all agents)

int flag_rate_setting_markup;           // Determines how the interest rates are set (additive VS multiplicative)
                                        // = 0 depends on an additive markup over the central bank interest rate
                                        // = 1 depends on a multiplicative markup over the central bank interest rate (old version)

int flag_rate_setting_loans;            // Determines how the interest rate on loans is set
                                        // = 0 depends on an additive markup over the central bank interest rate, and the ranking of demanders according to their debt service-to-sales ratio
                                        // = 1 depends on an additive markup over the central bank interest rate, and the default probability of the demander
                                        // = 2 depends on an additive markup over the central bank interest rate, and the CAR of the supplier
                                        // = 3 depends on an additive markup over the central bank interest rate, the default probability of the demander, and the CAR of the supplier

int flag_c_price_increase;              // Determines whether firms increase their markup at t_regime_shifts
                                        // = 0 no increase
                                        // = 1 set value to c_price_increase
                                        // = 2 increase value by c_price_increase %

int flag_change_income_shares;          // Determines whether there is a change in income shares during the simulation at t_regimes_shifts
                                        // = 0 no change
                                        // = 1 change income shares to ones defined in income_sh_regime_shift_mh
                                        // = 2 change wage ratio to values defined in income_sh_regime_shift_mh

int flag_change_apc;                    // Determines whether there is a change in average propensities to consume during the simulation at t_regimes_shifts
                                        // = 0 no change
                                        // = 1 change      

int flag_change_en_exp_shares;          // Determines whether there is a change in households energy expenditure shares during the simulation at t_regimes_shifts
                                        // = 0 no change
                                        // = 1 change to values specified by a.._mh_simul                           

int flag_change_dividends_tax_rate;     // Determines whether there is a change in dividends tax rate during the simulation at t_regimes_shifts
                                        // = 0 no change
                                        // = 1 change to aliqdiv_regime_shift
                                        // = 2 shift a share of personal income tax equal to aliqdiv_regime_shift from labour income to capital income, keeping taxes from labour income shares unaltered

int flag_change_wage_tax_rate;          // Determines whether there is a change in wages tax rates during during the simulation at t_regimes_shifts
                                        // = 0 no change
                                        // = 1 change to aliqw_regime_shift_mh
                                        // = 2 change to maintain same total tax on wages using rates ratios defined in aliqw_regime_shift_mh
                                        // = 3 change both tax rates on wages (using rates ratios defined in aliqw_regime_shift_mh) and on dividends (to aliqdiv_regime_shift) to maintain same total taxes on households
                                        // = 4 change tax rates based only on one parameter determining the progressivity of the system, aliqw_progressivity_regime_shift

int flag_change_wealth_tax_rate;        // Determines whether there is a change in wealth tax rate during the simulation at t_regimes_shifts
                                        // = 0 no change
                                        // = 1 change to aliqwealth_regime_shift

int flag_highly_progressive_taxation;   // Determines whether there is a change in labour income or dividends rate during the simulation at t_regimes_shifts
                                        // = 0 no change
                                        // = 1 increase tax rate on Managers labour income by aliq_households_increase
                                        // = 2 increase dividends' tax rate by aliq_households_increase
                                        // = 3 increase tax rate on Workers' labour income by aliq_households_increase

int flag_change_public_spending;        // Determines wether there is a change in public spending target during the simulation at t_regimes_shifts
                                        // = 0 no change
                                        // = 1 change to g_cons_GDP_ratio_regime_shift
                                        // = 2 change by g_cons_GDP_ratio_regime_shift * the change in tax revenues due to tax rate changes (if flag_change_dividends_tax_rate or flag_change_wage_tax_rate != 0)

int flag_env_friendl_in_market_shares;  // Determine wether environmental friendliness influences the market shares of the firms 
                                        // = 0 no
                                        // = 1 yes

int flag_environmental_subsidies_C_firms;   // Determine wether the government subsidises C-firms to purchase more environmentally-friendly machines 
                                            // = 0 no
                                            // = 1 yes, based on carbon intensity difference from average of machines supplied and as a percentage of average K-firms price (with max set by env_subsidy_per_machine)                                
                                            // = 2 as 1, but government subsidises C-firms to purchase dirtier machines (or taxes them if env_subsidy_per_machine has a negative value)

int flag_correlate_prod_and_green_tech;     // Determine whethere labour productivity growth is influenced by change in carbon intensity (NOT COMPATIBLE YET WITH PRODUCTIVITY SHOCKS)
                                            // = 0 no
                                            // = 1 yes, for C_firms growth rate of  productivity is influenced by growth rate of carbon intensity baed on the parameter correlation_prod_and_green_tech;         
                                            // = 2 yes, as 1 but only from t_regime_shift       
                                            // = 3 yes, perfect correlation from t_regime_shift. Energy intensity and environmental friendliness are determined from change in L productivity.
                                            // = 4 yes, as 3 but in addition re-scale the technology jump intervals according to rs_uu_lp, rs_uu_ee, rs_uu_ef

                                            #endif