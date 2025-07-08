#include "module_energy_sfc.h"

void EN_DEM_TOT(void)
{
    //Sum up energy demands coming from households, C-firms and K-firms
    D1_en_TOT=D1_en.Sum();
    D2_en_TOT=D2_en.Sum();
    
    D_en_TOT(1)=ROUND(D1_en_TOT+D2_en_TOT+D_en_h+D_en_g);
} 


void ENERGY_INV_PROD(void)
{
    //INVESTMENT-----------------------------------------------------------------------------
    //Determine existing productive capacity for energy
    K_ge=G_ge.Sum(); 
    K_de=G_de.Sum();

    //If existing capacity is insufficient to satisfy demand, expansion investment takes place
    if ((K_ge+K_de)<=D_en_TOT(1))
    {
        EI_en=D_en_TOT(1)-K_ge-K_de;   
        
        c_de_min=C_de(1)*100000;
        cf_min_ge=CF_ge(1)*100000;

        //Determine best dirty & green technology
        for (tt=1; tt<=t; tt++)
        {
            C_de(tt)=pf/A_de(tt)+t_CO2_en*EM_de(tt);
            if (C_de(tt)<c_de_min)
            {
                c_de_min=C_de(tt);
            }
            
            if (CF_ge(tt)<cf_min_ge) 
            {
                cf_min_ge=CF_ge(tt);
            }
        }

        if(flag_endogenous_exp_quota==1)
        {
            exp_quota=max(0.0, tanh(((c_de_min*payback_en-cf_min_ge)/cf_min_ge)/exp_quota_param));
        }
        
        //If green investment is not constrained, make all investment green if green is superior
        if (flag_energy_exp==0)
        {
            if (c_de_min*payback_en<cf_min_ge)
            {
                EI_en_de=EI_en;
                G_de(t)+=EI_en_de; 
            }
            else
            {
                EI_en_ge=EI_en;
                G_ge(t)+=EI_en_ge;
                IC_en_quota(t)=EI_en*cf_min_ge/payback_en;
                G_ge_n(t)+=EI_en_ge*cf_min_ge;
            }
        }
        //Otherwise, may be forced to also invest in dirty even when green is better
        else if(flag_energy_exp==1)
        {
            if (c_de_min*payback_en>cf_min_ge)
            {    
                if (EI_en<(exp_quota*K_gelag+(K_gelag-K_ge)))
                {
                    EI_en_ge=EI_en;
                    G_ge(t)+=EI_en_ge;
                    IC_en_quota(t)=EI_en*cf_min_ge/payback_en;
                    G_ge_n(t)+=EI_en_ge*cf_min_ge;
                }
                else
                {
                    EI_en_ge=(exp_quota*K_gelag+(K_gelag-K_ge));
                    G_ge(t)+=EI_en_ge;
                    IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
                    G_ge_n(t)+=EI_en_ge*cf_min_ge;
                    EI_en_de=EI_en-EI_en_ge;
                    G_de(t)+=EI_en_de; 
                }
            }
            else 
            {
                EI_en_de=EI_en;
                G_de(t)+=EI_en_de; 
            }
        }
        //Case in which share of investment in green is fixed exogenously
        else if (flag_energy_exp==2)
        {
            EI_en_ge=K_ge0_perc*EI_en;
            G_ge(t)+=EI_en_ge;
            IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
            G_ge_n(t)+=EI_en_ge*cf_min_ge;
            EI_en_de=EI_en-EI_en_ge;
            G_de(t)+=EI_en_de; 
        } 
        //Case in which investment in green changes to increase linearly green share until specified timestep and value
        else if (flag_energy_exp==3 || flag_energy_exp==4)
        {
            if (t > (t_regime_shifts + t_length_energy_transition)){ //Transition completed
                EI_en_ge=K_ge_END_perc*EI_en;
            } else if (t > t_regime_shifts){ //Transition ongoing
                // Target green share in current timestep
                K_ge_target_perc = (K_ge_END_perc - K_ge0_perc) / (t_length_energy_transition) * (t - t_regime_shifts) + K_ge0_perc;
                // Target expansion of green energy
                double K_ge_target = D_en_TOT(1) * K_ge_target_perc;

                EI_en_ge=K_ge_target-K_ge;
                if (EI_en_ge<0){
                    EI_en_ge=0;
                }
            } else { //Transition not started, as for when flag equal to 2
                EI_en_ge=K_ge0_perc*EI_en;
            }
            G_ge(t)+=EI_en_ge;
            IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
            G_ge_n(t)+=EI_en_ge*cf_min_ge;
            // Residual is invested in brown plants. If we invest in more green plants than what required for increasing productive capacity we still do it.
            if (EI_en_ge > EI_en){
                EI_en_de = 0;
            } else{
                EI_en_de=EI_en-EI_en_ge;
            }
            G_de(t)+=EI_en_de;
        }
        else //flag_energy_exp>4
        {
            if (c_de_min*payback_en>cf_min_ge)
            {    
                if (EI_en<(exp_quota*K_gelag+(K_gelag-K_ge)) || (K_delag/(K_gelag+K_delag))<=(exp_quota))
                {
                    EI_en_ge=EI_en;
                    G_ge(t)+=EI_en_ge;
                    IC_en_quota(t)=EI_en*cf_min_ge/payback_en;
                    G_ge_n(t)+=EI_en_ge*cf_min_ge;
                }
                else
                {
                    EI_en_ge=exp_quota*K_gelag+(K_gelag-K_ge);
                    G_ge(t)+=EI_en_ge;
                    IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
                    G_ge_n(t)+=EI_en_ge*cf_min_ge;
                    EI_en_de=EI_en-EI_en_ge;
                    G_de(t)+=EI_en_de; 
                }
            }
            else 
            {
                EI_en_ge=K_ge0_perc*EI_en;
                G_ge(t)+=EI_en_ge;
                IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
                G_ge_n(t)+=EI_en_ge*cf_min_ge;
                EI_en_de=EI_en-EI_en_ge;
                G_de(t)+=EI_en_de; 
            }
        }
        // New capacity
        K_ge=G_ge.Sum(); 
        K_de=G_de.Sum();
    } else if(flag_energy_exp==3 & (t > t_regime_shifts) & ((Q_ge/D_en_TOT(1)) < K_ge_END_perc)){  //To still invest to meet green energy share during the transition, even if this means to overinvest in energy capacity since there would not be the need
        if (t > (t_regime_shifts + t_length_energy_transition)){ //Transition completed
            K_ge_target_perc = K_ge_END_perc;
        } else{
            // Target green share in current timestep
            K_ge_target_perc = (K_ge_END_perc - K_ge0_perc) / (t_length_energy_transition) * (t - t_regime_shifts) + K_ge0_perc;
        }
        // Target expansion of green energy
        double K_ge_target = D_en_TOT(1) * K_ge_target_perc;

        EI_en_ge=K_ge_target-K_ge;

        if (EI_en_ge<0){
            EI_en_ge=0;
        }
        G_ge(t)+=EI_en_ge;
        IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
        G_ge_n(t)+=EI_en_ge*cf_min_ge;

        EI_en_de = 0;
        G_de(t)+=EI_en_de;
    }
    
    //Calculate ammortised investment cost from expansion of green energy
    for (tt=1; tt<=t; tt++)
    {
        if (t-tt<=payback_en && G_ge(tt)>0)
        {
            IC_en+=IC_en_quota(tt);
        } 
    }

    // PRODUCTION-----------------------------------------
    //Quantity of dirty energy is the residual
    Q_de=D_en_TOT(1)-K_ge; 

    //Produce dirty energy only if green is insufficient
    if (Q_de<0) //Only green energy
    {
        Q_de=0;
        Q_ge=D_en_TOT(1);
        //constant marginal cost for green energy
        c_infra=0;
        //Energy mark-up shock
        if((flag_energyshocks==1 && t==t_regime_shifts) || (flag_energyshocks==2 && t==t_regime_shifts))
        {
            mi_en_preshock=mi_en;
            mi_en=ratio_mi_en_shock*mi_en;
            mi_en_shock=mi_en;
        }
        
        c_en(1)=mi_en;
        Emiss_en=0;
    }
    else //Green energy NOT sufficient
    {
        Q_ge=K_ge;
        for (tt=1; tt<=t; tt++)
        {   
            C_de(tt)=pf/A_de(tt)+t_CO2_en*EM_de(tt);
        }
        G_de_temp=G_de;
        Q_de_temp=Q_de;
        
        //If dirty energy is needed, successively activate dirty plants, starting with the most efficient
        while (Q_de_temp>0)
        {
            c_de_min=C_de(1)*10;
            idmin=1;
            for (tt=1; tt<=t; tt++)
            {
                if (G_de_temp(tt)>0)
                {
                    if (C_de(tt)<=c_de_min)
                    {
                    idmin=tt;
                    c_de_min=C_de(idmin);
                    }
                }
            }
            
            if (Q_de_temp>G_de_temp(idmin))
            {
                PC_en+=G_de(idmin)*C_de(idmin);
                Emiss_en+=G_de(idmin)*EM_de(idmin);
                Q_de_temp-=G_de_temp(idmin);
                FuelCost+=G_de(idmin)*pf/A_de(idmin);
            }
            else
            {
                PC_en+=Q_de_temp*C_de(idmin);
                Emiss_en+=Q_de_temp*EM_de(idmin);
                FuelCost+=Q_de_temp*pf/A_de(idmin);
                Q_de_temp=0;
            }
            G_de_temp(idmin)=0;
            if(Q_de_temp<tolerance)
            {
                Q_de_temp=0;
            }
        }

        //Increasing marginal cost in dirty energy
        c_infra=c_de_min;

        //Determine energy price
        //Fossil fuel price shock
        if((flag_energyshocks==3 && t==t_regime_shifts) || (flag_energyshocks==4 && t==t_regime_shifts))
        {
            c_en_preshock=c_infra+mi_en;
            pf_preshock=pf;
            c_infra_t=2*c_en(2)-mi_en;
            i=0;
            while(c_infra<c_infra_t)
            {
                i=i+1;
                pf=pf_preshock*pow((1.01),i);
                for (tt=1; tt<=t; tt++)
                {   
                    C_de(tt)=pf/A_de(tt)+t_CO2_en*EM_de(tt);
                }
                G_de_temp=G_de;
                Q_de_temp=Q_de;
                PC_en=0;
                Emiss_en=0;
                FuelCost=0;
                
                while (Q_de_temp>0)
                {
                    c_de_min=C_de(1)*10;
                    idmin=1;
                    for (tt=1; tt<=t; tt++)
                    {
                        if (G_de_temp(tt)>0)
                        {
                            if (C_de(tt)<=c_de_min)
                            {
                            idmin=tt;
                            c_de_min=C_de(idmin);
                            }
                        }
                    }
                    
                    if (Q_de_temp>G_de_temp(idmin))
                    {
                        PC_en+=G_de(idmin)*C_de(idmin);
                        Emiss_en+=G_de(idmin)*EM_de(idmin);
                        Q_de_temp-=G_de_temp(idmin);
                        FuelCost+=G_de(idmin)*pf/A_de(idmin);
                    }
                    else
                    {
                        PC_en+=Q_de_temp*C_de(idmin);
                        Emiss_en+=Q_de_temp*EM_de(idmin);
                        FuelCost+=Q_de_temp*pf/A_de(idmin);
                        Q_de_temp=0;
                    }
                    G_de_temp(idmin)=0;
                    if(Q_de_temp<tolerance)
                    {
                        Q_de_temp=0;
                    }
                }
                c_infra=c_de_min;
            }
            pf_shock=pf;
        }

        //Energy mark-up shock
        if((flag_energyshocks==1 && t==t_regime_shifts) || (flag_energyshocks==2 && t==t_regime_shifts))
        {
            c_en_preshock=c_infra+mi_en;
            mi_en_preshock=mi_en;
            mi_en=ratio_mi_en_shock*(mi_en+c_infra)-c_infra;
            mi_en_shock=mi_en;
        }
        c_en(1)=c_infra+mi_en;
    }

    if (flag_energy_exp==4  && t>=t_regime_shifts){ //Reduce energy price relatively to green enery share increase
        c_en(1) *= (1 - renew_impact_on_p_e*(Q_ge/D_en_TOT(1)-K_ge0_perc));
    }
    
    if (flag_share_END==1)
    {
        share_de=K_de/(K_de+K_ge);
    }
    else if (flag_share_END==2)
    {
        share_de=1-Q_ge/D_en_TOT(1);
    }
    else
    {
        share_de=share_de_0;
    }
}

void ENERGY_RandD()
{
    //Determine energy sector revenues
    Rev_en=c_en(1)*D_en_TOT(1);

    //Determine R&D expenditure and split between clean and dirty
        //Check that desired expenditure in R&D is less than current availability (Revenues - green investments - cost of fossil production)
    if (Rev_en*share_RD_en<Rev_en-IC_en-PC_en) 
    {
        RD_en_de=max(0.0,share_RD_en*share_de*Rev_en);
        RD_en_ge=max(0.0,share_RD_en*(1-share_de)*Rev_en);
    }
    else
    {
        RD_en_de=max(0.0,share_de*(Rev_en-PC_en-IC_en));
        RD_en_ge=max(0.0,(1-share_de)*(Rev_en-PC_en-IC_en));
    }
    
    // Innovation in dirty energy sector
    parber_en_de=1-exp(-o1_en*RD_en_de);     
    Inn_en_de=bnldev(parber_en_de,1,p_seed);      
   
    if (Inn_en_de == 1) //If innovation takes place, update technologies
    {
        rnd=betadev(b_a11,b_b11,p_seed);          
        rnd=uu1_en+rnd*(uu2_en-uu1_en);                 
        A_de_inn=A_de(t)*(1+rnd);
        
        if (A_de_inn>1)
        {
            A_de_inn=1;
        }
        
        rnd=betadev(b_a11,b_b11,p_seed);          
        rnd=uu1_en+rnd*(uu2_en-uu1_en);                
        EM_de_inn=EM_de(t)*(1-rnd);
        
        if (EM_de_inn<0)
        {
            EM_de_inn=0;
        }

    }
    else
    {
        A_de_inn=A_de(t);
        EM_de_inn=EM_de(t);
    }
    
    //Innovation in green energy sector
    parber_en_ge=1-exp(-o1_en*RD_en_ge);    
    Inn_en_ge=bnldev(parber_en_ge,1,p_seed);

    if (Inn_en_ge == 1 && CF_ge(t)>0)
    {
        rnd=betadev(b_a11,b_b11,p_seed);         
        rnd=uu1_en+rnd*(uu2_en-uu1_en);
        CF_ge_inn=CF_ge(t)*(1-rnd);
        
        if (CF_ge_inn<0)
        {
            CF_ge_inn=0;
        }
    }
    else
    {
        CF_ge_inn=CF_ge(t);
    }
    
    K_gelag=max(K_ge,(K_ge+K_de)*K_ge0_perc);
    K_delag=K_de;

    if(t<=life_plant)
    {
        G_de(1)-=G_de_0/life_plant;
        G_ge(1)-=G_ge_0/life_plant;
        G_ge_n(1)-=G_ge_n_0/life_plant;
    }

    //Update technology vectors
    if (t < T)
    {
        if (pf/A_de_inn+t_CO2_en*EM_de_inn < pf/A_de(t)+t_CO2_en*EM_de(t))
        {
            A_de(t+1)=A_de_inn;        
            EM_de(t+1)=EM_de_inn;
            C_de(t+1)= pf/A_de_inn+t_CO2_en*EM_de_inn;
        }
        else
        {
            A_de(t+1)=A_de(t);         
            EM_de(t+1)=EM_de(t);
            C_de(t+1)= pf/A_de(t)+t_CO2_en*EM_de(t);
        }

        if (CF_ge_inn < CF_ge(t))
        {
            CF_ge(t+1)=CF_ge_inn;
        }
        else
        {
            CF_ge(t+1)=CF_ge(t);
        }
    
        for (tt=1; tt<=t; tt++)
        {
            if (t-tt>life_plant)
            {
                G_de(tt)=0;
                G_ge(tt)=0;
                G_ge_n(tt)=0;
            }
        }
    }

    //Update nominal value of Energy capital stock; Energy capital valued at production cost --> 0 for dirty
    CapitalStock_e(1)=G_ge_n.Sum();

    //Calculate total emissions
    Emiss_TOT(1)= Emiss_en+Emiss2_TOT+Emiss1_TOT;

    ENERGY_LABOUR();
}

void ENERGY_LABOUR(void)
{  
    //Compute labour demand associated to R&D
    LDen_rd_de_mh["wr"]=RD_en_de/w_tot_for_1_wr_mh(2);     
    LDen_rd_ge_mh["wr"]=RD_en_ge/w_tot_for_1_wr_mh(2);
    for (const string& cl:classes_mh){
      LDen_rd_de_mh[cl]=LDen_rd_de_mh["wr"]*ld_ratios_mh[cl];
      LDen_rd_ge_mh[cl]=LDen_rd_ge_mh["wr"]*ld_ratios_mh[cl];
    }

    //Compute labour demand associated to investment and total
    double LDexp_en_wr=IC_en/w_tot_for_1_wr_mh(2); //Number of workers for expansion investment
    for (const string& cl:classes_mh){
        double num_h = LDexp_en_wr*ld_ratios_mh[cl];
        LDen_exp_mh[cl]=num_h;
        LDen_tot_mh[cl]=num_h+LDen_rd_de_mh[cl]+LDen_rd_ge_mh[cl];  //Total number (including R&D)
    }

    //Calculate wages to be paid by energy sector
    Wages_en=0;
    for (const string& cl:classes_mh){
      Wages_en_mh[cl]=w_mh[cl](2)*LDen_tot_mh[cl]; //Wages paid by firm i to each h class
      Wages_en+=Wages_en_mh[cl]; //Add to total wages paid
    }
}

void EMISS_IND(void)
{
    //Calculate total emissions from C-firms and K-firms
    for (i=1; i<=N1; i++)
    {
        Emiss1(i)=A1p_ef(i)/A1p_en(i)*Q1(i);
    }
    
    for (j=1; j<=N2; j++)
    {
        Emiss2(j)= A2e_ef(j)/A2e_en(j)*Q2(j);
    }
    
    Emiss1_TOT=Emiss1.Sum();
    Emiss2_TOT=Emiss2.Sum();
}