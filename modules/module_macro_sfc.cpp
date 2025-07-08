#include "module_macro_sfc.h"

void LABOR(void)
{		
	//Update labour supply (if switched on)
	if (t>200)
	{
		LS=0;
		for (const string& cl : classes_mh){
			LS_mh[cl]*=g_ls;
			LS+=LS_mh[cl];
		}
	}

	//Check on L demand from last period.
	if((LD1_rd_mh["wr"]+LDen_tot_mh["wr"]) > LS_mh["wr"])
	{
		std::cerr << "\n\n ERROR: Remaining labour supply of Workers is negative in period " << t << endl;
		Errors << "\n Remaining labour supply of Workers is negative in period " << t << endl;
		exit(EXIT_FAILURE);
	}

	//Check on Labour supply and demand of Workers and in case scale back production
		//Calculate total labour demand for Workers
	LD1_wr=Ld1_prod_wr.Sum();
	LD2_wr=Ld2_wr.Sum();

	LSe_wr=LS_mh["wr"]; //Remaining labour supply of Workers
	LSe_wr-=(LD1_rd_mh["wr"]+LDen_tot_mh["wr"]);

		//Check
	if (LD2_wr + LD1_wr <= LSe_wr)
  	{
		LSe_wr=LSe_wr-LD1_wr-LD2_wr;							
  	}
	else //If total labour demand exceeds supply, production is scaled back
	{
		for (j=1; j<=N2; j++)
		{
            //Q2(j)=Q2(j)*LSe/(LD1tot+LD2tot); 
			Ld2_wr(j)=Ld2_wr(j)*LSe_wr/(LD1_wr+LD2_wr);
			Q2(j)=Ld2_wr(j)*A2e(j);
		}

		for (i=1; i<=N1; i++)
		{
    		Qpast=Q1(i);
			
      		if (Qpast > 0)
			{
       			Ld1_prod_wr(i)=Ld1_prod_wr(i)*LSe_wr/(LD1_wr+LD2_wr);
				Q1(i)=floor(Ld1_prod_wr(i)*((1-shocks_labprod1(i))*A1p(i)));
				reduction=Qpast-Q1(i);
				while(reduction>0)
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
						reduction-=1;
					}
				}
			}
		}
		// Recalculate total reduced labour demand
		LD1_wr=Ld1_prod_wr.Sum();
		LD2_wr=Ld2_wr.Sum();
		for (const string& cl:classes_mh){
			LD1_mh[cl]=LD1_wr*ld_ratios_mh[cl]+LD1_rd_mh[cl];
			LD2_mh[cl]=LD2_wr*ld_ratios_mh[cl];
		}
		for(j=1; j<=N2; j++)
		{
      		LD2_wr+=Ld2_wr(j);
		}
	}
	
    LD = 0;
	Benefits=0; 
	for (const string& cl:classes_mh){ 
		//Total labour demand
		double LD_cl = LD1_mh[cl]+LD2_mh[cl]+LDen_tot_mh[cl];
		LD_mh[cl]=LD_cl;	//Of the household class
		LD+=LD_cl;			//Overall

		//Unemployment rate of the household class
		U_mh[cl]=(LS_mh[cl]-LD_mh[cl])/LS_mh[cl];

		//Unemployment benefit payments by government
		double G_cl; //Transfer to the class
		if (LS_mh[cl]>LD_mh[cl]){
			G_cl = U_mh[cl]*LS_mh[cl]*w_mh[cl](2)*wu; //Based on corresponding wage
		}
		else {		//No unemployed in the class
			G_cl = 0;
		}
		// G_cl += Transfer_shock_mh[cl]; //Add transfer following energy price shock (in case)
		
		Benefits_mh[cl]=G_cl;		 //To the household class
		Benefits+=G_cl;				 //Overall
	}
	G=Benefits;

	//Overall unemployment rate
	U(1)=(LS-LD)/LS; //Overall
	//Determine carbon tax revenue transfer
	double co2TaxRev = Taxes_CO2(2)*redistribute_co2TaxRev;
		//Add to this energy price shock transfer in case to calculate total government transfer
	govTransfers=0;
	for (const string& cl:classes_mh){
		govTransfers_mh[cl] = co2TaxRev * co2TaxRev_sh_mh[cl] + Transfer_shock_mh[cl];
		govTransfers+=govTransfers_mh[cl];
	}
}

void MACRO(void)
{
	//Calculate macroeconomic aggregates, mean values etc
	ExpansionInvestment_r=EI.Row(1).Sum();
	ExpansionInvestment_n=EI_n.Sum();
	ReplacementInvestment_r=SI.Sum();
	ReplacementInvestment_n=SI_n.Sum();
	Investment_r=ExpansionInvestment_r+ReplacementInvestment_r;
	Investment_n=ExpansionInvestment_n+ReplacementInvestment_n;
	Q2tot=Q2.Sum();
	Q2dtot=Qd.Sum();
	D2tot=D2.Row(1).Sum();
  	Q1tot=Q1.Sum();
	LD2=LD1_wr+LD2_wr;

	//Labour and energy productivity
	Am_en_2 = 0;
	for (j=1; j<=N2; j++)
	{
		if (LD2>0)
		{
			Am_a+=Ld2_wr(j)/LD2*A2e(j);
			Am2+=Ld2_wr(j)/LD2_wr*A2e(j);
		}
		if((D2_en_TOT+D1_en_TOT)>0){
			Am_en(1)+=D2_en(j)/(D2_en_TOT+D1_en_TOT)*A2e_en(j);
		}
		if(D2_en_TOT>0){
			Am_en_2 += D2_en(j)/D2_en_TOT*A2e_en(j);
		}
		Am(1)+=A2_mprod(j);
		A_mi+=log(A2(j));
		A2_en_mi+=log(A2_en(j));
		A2_ef_mi+=log(A2_ef(j));
		H2+=f2(1,j)*f2(1,j);
  	}

	if(D1_en_TOT>0){
		for (i=1; i<=N1; i++){
			Am_en_1 += D1_en(i)/D1_en_TOT*A1p_en(i);
		}
	}

	A_mi/=N2r;
	A2_en_mi/=N2r;
	A2_ef_mi/=N2r;
	H2=(H2-1/N2r)/(1-1/N2r);

	for (j=1; j<=N2; j++)
	{
		A_sd+=(log(A2(j))-A_mi)*(log(A2(j))-A_mi);
	}

  	A_sd=sqrt(A_sd/N2);              

	for (i=1; i <=N1; i++)
	{
		if (Q1tot>0)
		{
			f1(1,i)=Q1(i)/Q1tot;
		}
		else
		{ 
			f1(1,i)=f1(2,i);
		}

		H1+=f1(1,i)*f1(1,i);                    
		A1_mi+=log(A1p(i));
		A1_en_mi+=log(A1p_en(i));
		A1_ef_mi+=log(A1p_ef(i));

		if (LD2>0)
		{
			Am_a+=Ld1_prod_wr(i)/LD2*A1p(i);
			Am1+=Ld1_prod_wr(i)/LD1_wr*A1p(i);
		}
		Am(1)+=A1p(i);
		if((D2_en_TOT+D1_en_TOT)>0){
			Am_en(1)+=D1_en(i)/(D2_en_TOT+D1_en_TOT)*A1p_en(i);
		}
	}

	Am(1)/=(N1r+N2r);
	A1_mi/=N1r;
	A1_en_mi/=N1r;
	A1_ef_mi/=N1r;
  	H1=(H1-1/N1r)/(1-1/N1r);                 

	CreditSupply_all=BaselBankCredit.Sum();
	CreditDemand_all=CreditDemand.Sum();

	for (i=1; i <=NB; i++)
	{
		if (CreditSupply_all > 0)
		{
			fB(1,i)=(BaselBankCredit(i)/CreditSupply_all);
		}
		else
		{
			fB(1,i)=fB(2,i);
		}

		HB+=fB(1,i)*fB(1,i);
	}
	
	//Calculate GDP
	GDP();

	//Update wage rate
  	WAGE();
}

void GDP(void){
	//OLD APPROACH---------------------
	GDP_r(1)=Q1tot*dim_mach+Q2tot;
	GDP_n(1)=0;
	double GDP_n_1 = 0; //Contribution to GDP of K-sector
	double GDP_n_2 = 0;	//Contribution to GDP of C-sector
	for (i=1; i <=N1; i++)
	{
		GDP_n_1+=Q1(i)*p1(i);
	}
	for (i=1; i <=N2; i++)			//Goods consumption
	{
		GDP_n_2+=Q2(i)*p2(i);
	}
	GDP_n(1)+=GDP_n_1;
	GDP_n(1)+=GDP_n_2;
	
	// Calculate growth rates
	if(t>1)
	{
		GDP_rg= log(GDP_r(1))-log(GDP_r(2));
		GDP_ng= log(GDP_n(1))-log(GDP_n(2));
	}

	//NEW APPROACH---------------------
	//Expenditure approach
	RowVector InventoriesChange_i = Inventories.Row(1) - Inventories.Row(2);
	Exports = FuelCost;
	Imports = 0;
	GDP_n_exp = Expenditure_tot_h + Exp_tot_g + InventoriesChange_i.Sum() + Investment_n - Exports + Imports;

	//Production approach
	RowVector VA_1_i(N1);						//Value added by each K-firm
	VA_1_i = S1 - EnergyPayments_1;
	RowVector VA_2_i(N2);						//Value added by each C-firm
	VA_2_i = S2.Row(1) + InventoriesChange_i - EnergyPayments_2 - LoanInterest_2;
	VA_1 = VA_1_i.Sum(); 				//Value added by K-sector
	VA_2 = VA_2_i.Sum(); 				//Value added by C-sector
	VA_en = EnergyPayments - FuelCost; 	//Value added by energy sector
	VA_b = LoanInterest.Sum(); //+InterestBonds_b.Sum()+InterestReserves_b.Sum();   //Value added by banks
	double TaxesOnProduction = 0;
	double SubsidiesOnProduction = 0;
	GDP_n_prod = VA_1 + VA_2 + VA_en + VA_b + TaxesOnProduction - SubsidiesOnProduction;

	//Income approach
	PaymentsToLabour = Wages+Bonuses_h;
	double TotalTaxes = Taxes_1.Sum() + Taxes_CO2_1.Sum() + Taxes_2.Sum() + Taxes_CO2_2.Sum() + Taxes_CO2_e + Taxes_e_shock + Taxes_e_ff_shock + Taxes_b.Sum(); //+ Taxes_f_ff_shock; //All apart from households
	double TotalSubsidies = EnvSubsidies_2;
	PaymentsToGovernment = TotalTaxes - TotalSubsidies;
		//GOS = VA - payments to labour and government
	double GOS_1; 						//Gross operating surplus K-sector
	GOS_1 = VA_1 - Wages_1_i.Sum() - Bonuses_1_i.Sum() - Taxes_1.Sum() - Taxes_CO2_1.Sum();
	double GOS_2; 						//Gross operating surplus C-sector
	GOS_2 = VA_2 - Wages_2_i.Sum() - Bonuses_2_i.Sum() - Taxes_2.Sum() - Taxes_CO2_2.Sum() + EnvSubsidies_2;
	double GOS_en;  					//Gross operating surplus energy sector
	GOS_en = VA_en - Wages_en - Bonuses_e - Taxes_CO2_e + Taxes_e_shock + Taxes_e_ff_shock;
	double GOS_b; 						//Gross operating surplus banks
	GOS_b = VA_b - Taxes_b.Sum();
	PaymentsToCapital = GOS_1 + GOS_2 + GOS_en + GOS_b;
	GDP_n_inc = PaymentsToLabour + PaymentsToGovernment + PaymentsToCapital;
	//double GDP_n_inc_check = Wages + InterestDeposits_h + Dividends_h(1) + undistr profits + Taxes_1.Sum() + Taxes_2.Sum() + transfers; //Disaggregating GOS

		// Real GDP
	GDP_n_new(1) = GDP_n_exp;
	GDP_r_new(1) = GDP_n_exp/cpi(1);

		// Checks on GDP calculation with different approaches
	if (GDP_n_exp < 0 | GDP_n_prod < 0 | GDP_n_inc < 0){
		std::cerr << "\n\n ERROR: GDP negative in period " << t << endl;
		Errors << "\n GDP negative in period "<< t << endl;
	}
	double dev=fabs((GDP_n_exp-GDP_n_prod))/GDP_n_exp;
	if(dev>tolerance)
	{
		std::cerr << "\n\n ERROR: GDP calculated through expenditure and production approach differ by "<< dev*100 << "% in period " << t << endl;
		Errors << "\n GDP calculated through expenditure and production approach differ by "<< dev*100 << "% in period "<< t << endl;
	}
	dev=fabs((GDP_n_exp-GDP_n_inc))/GDP_n_exp;
	if(dev>tolerance)
	{
		std::cerr << "\n\n ERROR: GDP calculated through expenditure and income approach differ by "<< dev*100 << "% in period " << t << endl;
		Errors << "\n GDP calculated through expenditure and income approach differ by "<< dev*100 << "% in period "<< t << endl;
	}
	dev=fabs((GDP_n_inc-GDP_n_prod))/GDP_n_inc;
	if(dev>tolerance)
	{
		std::cerr << "\n\n ERROR: GDP calculated through income and production approach differ by "<< dev*100 << "% in period " << t << endl;
		Errors << "\n GDP calculated through income and production approach differ by "<< dev*100 << "% in period "<< t << endl;
	}

	// Growth rates
	if(t>1)
	{
		GDP_rg_new = log(GDP_r_new(1))-log(GDP_r_new(2));
		GDP_ng_new = log(GDP_n_new(1))-log(GDP_n_new(2));
	}
		//Check with previous way of calculating it
	/* dev=(GDP_ng-GDP_ng_new)/GDP_ng;
	if(fabs(dev)>tolerance)
	{
		std::cerr << "\n\n ERROR: Nominal GDP growth calculated through expenditure and OLD approach differ by "<< dev*100 << "% in period " << t << endl;
		Errors << "\n Nominal GDP growth calculated through expenditure and OLD approach differ by "<< dev*100 << "% in period "<< t << endl;
	}
	dev=(GDP_rg-GDP_rg_new)/GDP_rg;
	if(fabs(dev)>tolerance)
	{
		std::cerr << "\n\n ERROR: Real GDP growth calculated through expenditure and OLD approach differ by "<< dev*100 << "% in period " << t << endl;
		Errors << "\n Real GDP growth calculated through expenditure and OLD approach differ by "<< dev*100 << "% in period "<< t << endl;
	} */
}

void WAGE(void)
{
	d_U=(U(1)-U(2));

	d_cpi=cpi(1)/cpi(5)-1;

	d_Am=kappa*d_Am+(1-kappa)*((Am(1)-Am(2))/Am(2));

	dw=d_cpi_target + psi1*(pow(1+d_cpi-d_cpi_target_a,0.25)-1) + psi2*d_Am - psi3*d_U;

	if(dw>mdw) //Limit to wage change
	{
		dw=mdw;
	}
	if(dw<(-mdw))
	{
		dw=(-mdw);
	}

	//Wage update    
	w_mh["wr"](1)=w_mh["wr"](2)*(1+dw);
	if (w_mh["wr"](1) < w_min-0.001) //Minimum wage check
	{
		w_mh["wr"](1)=w_min;
	}
	w_tot_for_1_wr_mh(1)=0;
	for (const string& cl : classes_mh){
		w_mh[cl](1) = w_mh["wr"](1) * w_ratios_mh[cl];	 //Household classes wage update
		w_tot_for_1_wr_mh(1) += (w_mh[cl](1) * ld_ratios_mh[cl]);  //Update total wage to pay for each worker 
	}
} 

void GOV_BUDGET(void)
{
	//If outstanding government debt is greater than 0, need to take bond repayments & interest into account when calculating borrowing requirement
	if(GB(2)>0)
	{
		TransferCB=ProfitCB(2);
		InterestBonds=r_bonds*GB(2);
		InterestBonds_cb=r_bonds*GB_cb(2);
		Deficit=G+govTransfers+r_bonds*GB(2)+Bailout+EntryCosts+Transfer_shock_f+Exp_tot_g+EnvSubsidies_2-Taxes_g-TransferCB-Taxes_CO2(1)-Taxes_e_shock-Taxes_f_ff_shock-Taxes_e_ff_shock;
		if((-Deficit)>GB(2))
		{
			PSBR=Deficit;
			BondRepayments_cb=0;
			for(i=1; i<=NB;i++)
			{
				InterestBonds_b(i)=r_bonds*GB_b(2,i);
				BondRepayments_b(i)=0;
			}
		}
		else if((-Deficit)>GB_b.Row(2).Sum())
		{
			PSBR=Deficit+bonds_share*GB_cb(2);
			BondRepayments_cb=bonds_share*GB_cb(2);
			for(i=1; i<=NB;i++)
			{
				InterestBonds_b(i)=r_bonds*GB_b(2,i);
				BondRepayments_b(i)=0;
			}
		}
		else
		{
			PSBR=Deficit+bonds_share*GB(2);
			BondRepayments_cb=bonds_share*GB_cb(2);
			for(i=1; i<=NB;i++)
			{
				InterestBonds_b(i)=r_bonds*GB_b(2,i);
				BondRepayments_b(i)=bonds_share*GB_b(2,i);
			}
		}
	}
	else
	{
    	//If gov. debt is negative, government earns reserve rate on deposits with CB
		InterestBonds=-r_cbreserves*GB(2);
		InterestBonds_cb=-r_cbreserves*GB_cb(2);
		TransferCB=ProfitCB(2);
		BondRepayments_cb=0;
		for(i=1; i<=NB;i++)
		{
      		InterestBonds_b(i)=0;
      		BondRepayments_b(i)=0;
    	}
		Deficit=G+r_bonds*GB(2)+Bailout+EntryCosts+Transfer_shock_f+Exp_tot_g-Taxes_g-TransferCB-Taxes_CO2(1)-Taxes_e_shock-Taxes_f_ff_shock-Taxes_e_ff_shock;
		PSBR=Deficit;
	}

	//If government debt is smaller than 0 it is treated as a government deposit at the CB. This is first run down before new borrowing happens
	if(PSBR>0 && GB(2)<0)
	{
		if((-GB(2))>=PSBR)
		{
			GB(1)+=PSBR;
			GB_cb(1)+=PSBR;
			PSBR=0;
		}
		else
		{
			PSBR=PSBR+GB(2);
			GB_cb(1)=0;
			GB(1)=0;
		}
	}

	//Government needs to borrow
	if(PSBR>=0)
	{
		//Determine supply of new bonds and possibly banks' demand for bonds
		NewBonds=PSBR;
		for(i=1; i<=NB;i++)
		{
			if(BankProfits(1,i)>0)
			{
				BankProfits_temp(i)=(1-aliqb)*BankProfits(1,i);
			}
			else
			{
				BankProfits_temp(i)=0;
			}
    	}

		for(i=1; i<=NB; i++)
		{
			bonds_dem(i) = max(0.0,varphi * Loans_b(1,i)-GB_b(1,i));
		}
		
		bonds_dem_tot=bonds_dem.Sum();

		for(i=1; i<=NB;i++)
		{	
			//If there is excess demand for bonds, banks buy minimum between their demand and a share determined by their relative profits
			if (bonds_dem_tot >= PSBR & bonds_dem(i) >= (BankProfits_temp(i)/BankProfits_temp.Sum())*PSBR)
			{
				bonds_purchased(i) = (BankProfits_temp(i)/BankProfits_temp.Sum())*PSBR;
				GB_b(1,i) += bonds_purchased(i);
				GB(1)+=bonds_purchased(i);
				Outflows(i) +=bonds_purchased(i);
				NewBonds -= bonds_purchased(i);
			}
			else if (bonds_dem_tot >= PSBR & bonds_dem(i) < (BankProfits_temp(i)/BankProfits_temp.Sum())*PSBR)
			{
				bonds_purchased(i) = bonds_dem(i);
				GB_b(1,i) += bonds_purchased(i);
				GB(1)+=bonds_purchased(i);
				Outflows(i) +=bonds_purchased(i);
				NewBonds -= bonds_purchased(i);
			}
			//If there is excess supply of bonds, demand is fully satisfied
			else if (bonds_dem_tot < PSBR)
			{
				bonds_purchased(i) = bonds_dem(i);
				GB_b(1,i) += bonds_purchased(i);
				GB(1)+=bonds_purchased(i);
				Outflows(i) +=bonds_purchased(i);
				NewBonds -= bonds_purchased(i);
			}
		} 
		//Central bank buys remaining bonds
		GB_cb(1)+=max(0.0,NewBonds);
		GB(1)+=max(0.0,NewBonds);
	}
	//Government is running a surplus
	else
	{
		//If surplus is sufficient to repay all outstanding bonds held by banks, repay them and then repay the CB (possibly making GB_cb negative)
		if((-PSBR)>=GB_b.Row(2).Sum())
		{
			for(i=1; i<=NB;i++)
			{
				Inflows(i)+=GB_b(2,i);
				GB(1)-=GB_b(2,i);
				PSBR+=GB_b(2,i);
				GB_b(1,i)=0;
			}
			GB_cb(1)+=PSBR;
			GB(1)+=PSBR;
		}
		//Otherwise repay on bonds held by banks
		else
		{
			Bond_share=GB_b.Row(2)/GB_b.Row(2).Sum();
			for(i=1; i<=NB;i++)
			{
				Inflows(i)-=(PSBR*Bond_share(i));
				GB_b(1,i)+=(PSBR*Bond_share(i));
				GB(1)+=(PSBR*Bond_share(i));
			}
		}
	}

	//Make interest and principal payments on bonds
	for(i=1; i<=NB; i++)
	{
		if(GB_b(1,i)>0)
		{
			Inflows(i)+=InterestBonds_b(i)+BondRepayments_b(i);
			GB_b(1,i)-=BondRepayments_b(i);
			GB(1)-=BondRepayments_b(i);
		}
		else
		{
			Inflows(i)+=InterestBonds_b(i);
		}
	}

	if(GB_cb(1)>0)
	{
		GB_cb(1)-=BondRepayments_cb;
		GB(1)-=BondRepayments_cb;
	}

}

void TAYLOR(void)
{ 
  //Update monetary policy rate & all other rates linked to it
  if(flag_energyshocks>0 && flag_energyshocks_MP==0 && t>t_regime_shifts && t<(t_regime_shifts+9))
  {
	inflation_a=d_cpi_target_a;
  }
  else
  {
  	inflation_a=cpi(1)/cpi(5)-1;
  }
  
  r_a=(r_base+ taylor1*(inflation_a-d_cpi_target_a)+taylor2*(ustar-U(1)));
  r=taylor*r+(1-taylor)*(pow((1+r_a),0.25)-1);

  if(r<=0)
  {
    r = 0.000001;
  }


  //INTEREST RATE SETTING

  if(flag_rate_setting_markup==0 || flag_rate_setting_markup==2)
  {
	r_deb=r+bankmarkup;
      
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

  if(flag_rate_setting_markup==1)
  {
    r_depo=r*(1-bankmarkdown);
    r_cbreserves=r*(1-centralbankmarkdown);
    r_bonds=r*(1-bondsmarkdown);
  }
  
} 
