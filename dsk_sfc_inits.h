#ifndef INITS_H
#define INITS_H

double           A0;                                                 // Initial labour productivity
double           LS0_wr;                                                // Initial labour supply of workers
double           W10;                                                // Initial deposits K-firms
double           W20;                                                // Initial deposits C-firms
double           L0;                                                 // Initial loans C-firms
double           w0_wr;                                                 // Initial wage rate of Workers
double           K0;                                                 // Initial capital stock C-firms
double           bankmarkup_init;                                    // Initial bank mark-up
double           FirmDefaultProbability_init;                        // Initial C-firms' default probability            
double           A0_en;                                              // Initial energy efficiency
double           A0_ef;                                              // Initial environmental friendliness of machines (use) 
double           A_ef_ratio0;                                        // Ratio of initial environmental friendliness of production technique over of machine (use)       
double           K_ge0_perc;                                         // Initial capacity of green energy
double           pf0;                                                // Initial price of fossil fuel
double           mi_en0;                                             // Initial mark-up of energy sector
double           A_de0;                                              // Initial thermal efficiency of dirty energy plants
double           EM0;                                                // Initial emissions per unit of energy produced
double           CF_ge0;                                             // Initial cost of expanding green energy capacity
double           t_CO2_0;                                            // Initial carbon tax on firms
double           t_CO2_en_0;                                         // Initial carbon tax on energy sector
double           d_Am_init;                                          // Initial change in mean productivity
double           D_h0;                                               // Initial deposits households
RowVector        D_h0_sh_mh;                                         // Initial deposit shares household classes
double           NW_b0;                                              // Initial net worth banking sector
double           pm;                                                 // Initial productivity ratio between K and C sectors
double           D_e0;                                               // Initial energy sector deposits


double           Emiss_yearly_0;                                     // Initial annual emissions
double           Cum_emissions_0;                                    // Initial cumulative emissions
double           T_0_cumemiss;                                       // Initial temperature anomaly using cumulative emissions climate box
double           NPP0;                                               // Initial net primary production of plants in GtC per year
double           Con00;                                              // Pre-industrial ocean carbon content per meter of depth
double           Cat0;                                               // Pre-industrial atmospheric carbon
double           Catinit0;                                           // Initial atmospheric carbon in GtC
RowVector        Coninit0;                                           // Oceanic carbon in each layer in GtC
RowVector        Honinit0;                                           // Ocean heat content in each layer (per m^2 of ocean surface) in J/m^2
double           Tmixedinit0;                                        // Initial surface temperature (above pre-industrial value) 
double           biominit0;                                          // Initial carbon stock of biosphere
double           huminit0;                                           // Initial carbon stock of humus 
double           Catinit1;                                           // Initial atmospheric carbon in GtC
RowVector        Coninit1;                                           // Oceanic carbon in each layer in GtC
RowVector        Honinit1;                                           // Ocean heat content in each layer (per m^2 of ocean surface) in J/m^2
double           Tmixedinit1;                                        // Initial surface temperature (above pre-industrial value) 
double           biominit1;                                          // Initial carbon stock of biosphere
double           huminit1;                                           // Initial carbon stock of humus 

#endif
