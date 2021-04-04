/******************************************************************************
PART 1. Replication of the main findings of:

The War on Poverty's Experiment in Public Medicine: 
The Impact of Community Health Centers on the Mortality of Older Americans
by Martha Bailey and Andrew Goodman-Bacon
******************************************************************************/
* Set up the working directories (PLEASE MAKE SURE YOU CHANGED DIRECTORIES BEFORE RUNNING THE CODE!!!)

*dofile directory (where this file is stored)
global dofile "C:\Users\btpta\Desktop\Metrics 2\Replication"

*data directory (where posted datasets are stored)
global data "C:\Users\btpta\Desktop\Metrics 2\Replication\20120070_data\aer_data"

*output directory (where regression output, figures, and logs are saved)
global output "C:\Users\btpta\Desktop\Metrics 2\Replication\Output files"


cd "$output"


/***********************************************
FIGURE 1. AMR by Age Group, 1959-1988 
***********************************************/
clear
clear matrix
clear mata

*collapse mean annual mortality rates by group
use fips stfips cofips year amr* imr nnmr pnmr copop* births if year<1989 using "$data\aer_data",clear

*make manual weights for each age group's mortality rate
cap drop testo
egen testo = total(copop), by(year)
replace amr  = amr*(copop/testo)  // AMR all ages

cap drop testo
egen testo = total(births), by(year)
replace imr  = imr*(births/testo) 		// infant MR
replace nnmr  = nnmr*(births/testo)		//  Neonatal Infant Mortality Rate
replace pnmr  = pnmr*(births/testo)		//Postneonatal Infant Mortality Rate

cap drop testo
egen testo = total(copop_ch), by(year) 		//Population, 1-19
replace amr_ch  = amr_ch*(copop_ch/testo)	

cap drop testo
egen testo = total(copop_ad), by(year)		// Population, 20-49
replace amr_ad  = amr_ad*(copop_ad/testo)

cap drop testo
egen testo = total(copop_eld), by(year)		//Population, 50+
replace amr_eld  = amr_eld*(copop_eld/testo)

collapse (sum) amr imr nnmr pnmr amr_ch amr_ad amr_eld, by(year)

graph set window fontface "Times New Roman"							

*AMR
#delimit ;
scatter amr year,
connect(l l l l l)	
msymbol(i i)								
ytitle("Deaths per 100,000 Residents")	
xtitle("")		
title("A. Age-Adjusted Mortality", size(medium) color(black))
saving(panel1, replace);
#delimit cr;	

*Infant Mortality
#delimit ;
scatter imr nnmr pnmr year,
connect(l l l l l)		
msymbol(p Th Sh)	
mcolor(navy  red green)	
lcolor(navy  red green)
legend(off)
text(26 1962 "Total", size(medium) place(e) yaxis(1)) 
text(20 1959 "Neonatal", size(medium) place(e) yaxis(1)) 
text(8.5 1959 "Post-Neonatal", size(medium) place(e) yaxis(1)) 
xtitle("")
ytitle("Deaths per 1,000 Live Births")
title("B. Infants", size(medium) color(black))
saving(panel2, replace);
#delimit cr;	

*Child and Adult Mortality
#delimit ;
scatter amr_ch amr_ad year,
connect(l l l l l)
msymbol (p Th )	
mcolor(navy red)
lcolor(navy red)	
legend(off)
xtitle("Year")
ytitle("Deaths per 100,000 Residents")
title("C. Children (1-19) and Adults (20-49)", size(medmsall) color(black))
text(250 1980 "Adults", size(medium))
text(100 1965 "Children", size(medium))			
saving(panel3, replace);
#delimit cr;	

*Older Adult Mortality
#delimit ;
scatter amr_eld year,
connect(l l l l l)	
msymbol(i i)	
legend(off)
xtitle("Year")
ytitle("Deaths per 100,000 Residents")
title("D. Older Adults (50+)", size(medmsall) color(black))
saving(panel4, replace);
#delimit cr;	


graph combine panel1.gph panel2.gph panel3.gph panel4.gph, col(2)  xcommon xsize(8.5) ysize(5.5) 
graph display, xsize(8.5) ysize(5.5)

graph export "$output/figure2.wmf", replace

/******************************************************
TABLE 1: Mean 1960 Characteristics 
******************************************************/
clear all

use "$data/aer_data.dta", clear
egen amr65 = total(amr*(year==1965)),by(fips)
cap drop _59medfaminc 
gen dms = _tot_med_stud>0
drop _tot_med_stud
recode stfips (9 23 25 33 44 50 34 36 42 = 1) ///
			  (18 17 26 39 55 19 20 27 29 31 38 46 = 2) ///
			  (10 11 12 13 24 37 45 51 54 1 21 28 47 5 22 40 48 = 4) ///
			  (4 8 16 35 30 49 32 56 6 41 53 = 5), gen(region)
xi i.region, noomit
drop if fips==6037|fips==17031|fips==36061|fips==3011
xi i.exp1, pref(_T)
qui xtreg amr D* R* H* _Texp* [aw=popwt], fe 		/*to get paper estimation sample*/
keep if e(sample)
keep if year==1960
keep *fips* chc_year_exp _* amr dms copop year region
drop _60medschlmt24
replace _tot_act_md = _tot_act_md/copop*1000
drop _Texp*

merge 1:1 stfips cofips using "$output/pscore_temp"
keep if _merge==3
drop _merge

*rescale DFL weights
egen totdfl = total(dflwgt1), by(treat1)
replace dflwgt1 = dflwgt1/totdfl
count if treat1
replace dflwgt1 = 1/r(N) if treat1

*by timing
preserve
	egen chc = cut(chc_year_exp), at(1965, 1968, 1971, 1975, 1981)
	keep if chc<.
	collapse (mean) _* copop amr dms (count) fips, by(chc)
	xpose, varname clear
	ren v1 _65_67
	ren v2 _68_70
	ren v3 _71_74
	ren v4 _75_80
	save "$output/table1", replace
restore

*by treatment
preserve
	gen chc = chc_year_exp<=1974
	collapse (mean) _* copop amr dms (count) fips, by(chc)
	xpose, varname clear
	ren v1 _other
	ren v2 _pre74
	merge 1:1 _varname using "$output/table1"
	drop _merge
	save "$output/table1", replace
restore

preserve
	gen chc = chc_year_exp<=1974
	gen p = .
	gen varname = ""
	local i = 2
	foreach var of varlist _* copop dms amr{
		reg `var' chc
		testparm chc
		replace p = r(p) in `i'
		replace varname = "`var'" in `i'
		local i = `i'+1
	}
	ren varname _varname
	keep if p<.
	keep p _varname
	merge 1:1 _varname using "$output/table1"
	drop _merge
	save "$output/table1", replace
restore

*reweighted means
preserve
	gen chc = chc_year_exp<=1974
	xi i.region, noomit

	collapse (mean) _* dms amr copop (count) fips [aw=dflwgt1], by(chc)
	drop if chc

	xpose, clear varname
	ren v1 dfl_
	merge 1:1 _varname using "$output/table1"
	drop _merge
	save "$output/table1", replace
restore 

*do the DFL-weighted t-tests
preserve
	*resample pscores 1000 times and calculate a t-stat for EACH covariate on EACH rep
	cap gen chc = chc_year_exp<=1974
	xi i.region, noomit
	gen _varname = ""
	gen pdfl = .

	set seed 12345
	forval i = 1/1000{
		cap drop b_*
		*draw new error terms and apply them to the index and feed through normal CDF
		gen b_pscore = normal(index1 + rnormal())
		*gen dfl weights
		sum chc
		local ED 					= r(mean)
		gen b_dflwgt				= (b_pscore/`ED')*((1-`ED')/(1-b_pscore))*(1-chc)			/*see DiNardo (2002) about only applying these weights to the control group for TOT*/
		egen b_totdfl 				= total(b_dflwgt), by(chc)
		replace b_dflwgt			= b_dflwgt/b_totdfl
		count if chc
		replace b_dflwgt 			= 1/r(N) if chc
		*do 1000 t-test of differences in dfl weighted means to get a bootstrap t-distribution
		foreach var of varlist copop _60pcturban _60pctrurf _I* _60pct04years _60pctmt64years _60pctnonwhit _60pctlt4schl	_60pctmt12schl _pct59inclt3k _pct59incmt10k _tot_act_md dms amr{
			reg `var' chc [aw=b_dflwgt]
			testparm chc
			cap gen pct_`var' = sqrt(r(F)) in `i'
			replace pct_`var' = sqrt(r(F)) in `i'
		}
	}
	gen ind = _n
	sum pct*
	keep ind pct*
	keep if pct_copop<.
	save "$data/pct_t_bootstrap", replace
restore

*do t-test and compare to the bootstrap t-distributions (not tabulated t-dist)
gen chc = chc_year_exp<=1974
xi i.region, noomit

gen _varname = ""
gen pdfl = .
gen ind = .
local i = 2
foreach var of varlist copop _60pcturban _60pctrurf _Iregion_1 _Iregion_2 _Iregion_4 _Iregion_5 _60pct04years _60pctmt64years _60pctnonwhit _60pctlt4schl	_60pctmt12schl _pct59inclt3k _pct59incmt10k _tot_act_md dms amr{
	reg `var' chc [aw=dflwgt1]
	testparm chc
	local t = sqrt(r(F))
	preserve
		use "$data/pct_t_bootstrap", replace
		count if pct_`var'>=abs(`t')
		local p = r(N)/1000
	restore
	replace pdfl = `p' in `i'
	replace _varname = "`var'" in `i'
	replace ind = `i' if _varname=="`var'"
	local i = `i'+1
}
keep if pdfl<.
keep pdfl _varname ind
merge 1:1 _varname using "$output/table1"
drop _merge
order _varname _65 _68 _71 _75 _pre _other p dfl_ pdfl
replace ind = 1 if _varname=="fips"
replace ind = 0 if _varname=="chc"
sort ind
save "$output/table1", replace

export excel using "$output/table1.xls", replace firstr(var)

erase "$output/table1.dta"
erase "$data/pct_t_bootstrap.dta"


*********************************************************************
*******************Calculate weights based on pscore ****************
/************************************
 ESTIMATE PROPENSITY SCORES EQUATION
 *We construct propensity scores by estimating a probit with the binary dependent variable equal to 1 if a county received a CHC from 1965 to 1974 using the following covariates
*************************************/

use "$data/aer_pscore_data", clear

*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61					
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	

drop if stfips==2|stfips==15
ren _copop x_copop					//Total population over 50 in 1960
ren _copop2 x_copop2				//Total population over 50 in 1960 Squared

*1. 
local early "NO POPULATION, LINEARLY, WEIGHTED, EARLY CENTERS"
probit treat1 _* [pw=x_copop]				  // models using the early CHCs 1965-1974 on the county sample 
estimates store ps1							 //observed 1959- 1988 

*2. 
local all "NO POPULATION, LINEARLY, WEIGHTED, ALL CENTERS"
probit treat2 _* [pw=x_copop]				// // models using all CHCs
estimates store ps2

*give missing values the state mean in order to predict pscores...imputation not used in the estimation...this may be wrong, but its only for 36 small counties
foreach var of varlist _*{
	egen testo = mean(`var'),by(stfips)
	replace `var' = testo if `var'==.
	drop testo
}

*PREDICT PROPENSITY SCORES
forval i = 1/2{
	estimates restore ps`i'
	predict index`i', xb
	gen pscore`i' = normal(index`i')
}
label var pscore1 early
label var pscore2 late

keep *fips treat? pscore? ind*

*CREATE DFL WEIGHTS (the method proposed by DiNardo,Fortin and Lemieux (1996) (henceforth DFL) to analyze the effect of covariates)

sum treat1
local ED 					= r(mean)					// share if treated counties early 
gen dflwgt1 			= (pscore1/`ED')*((1-`ED')/(1-pscore1))*(1-treat1)	

sum treat2
local ED 					= r(mean)					// share if treated counties
gen dflwgt2 			= (pscore2/`ED')*((1-`ED')/(1-pscore2))*(1-treat2)	

sort stfips cofips
save "$output/pscore_temp", replace


/*********************************************************************
FIGURE 2: Predicting the timing of CHC grant with 1960 Characteristics
*********************************************************************/
clear all

use *fips* year chc_year_exp amr _* copop using "$data/aer_data" if chc_year_exp<=1974, clear
drop if fips==6037|fips==17031|fips==36061

*add the 1965 AMR value to all obs and then keep year 1960
egen amr65 = total(amr*(year==1965)), by(fips)
keep if year==1960
gen damr = amr65-amr


************
* 1965 AMR *
************
*univariate regression
reg amr65 chc_year_exp [aw = copop], robust
predict amrhat
local bu = floor(_b[chc_year_exp]*10)/10
local seu = floor(_se[chc_year_exp]*10)/10

*multivariate
#delimit ;
reg amr65	
_60pcturban _60pctrurf _60pct04years 
_60pctmt64years _60pctnonwhit _60pctmt12schl 
_60pctlt4schl _pct59inclt3k _pct59incmt10k 
_tot_act_md  chc_year_exp [aw = copop] , robust;
#delimit cr;
local ba = floor(_b[chc_year_exp]*10)/10
local sea = floor(_se[chc_year_exp]*10)/10

*get best-fit line for average x's and variable chc_year
preserve
foreach var of varlist _60pcturban _60pctrurf _60pct04years _60pctmt64years _60pctnonwhit _60pctmt12schl _60pctlt4schl _pct59inclt3k _pct59incmt10k _tot_act_md{
	egen testo = mean(`var')
	replace `var' = testo
	drop testo
}
predict amrhatadj
keep stfips cofips  amrhatadj
save "amrhatadj", replace
restore
merge 1:1 stfips cofips  using "amrhatadj"
drop _merge

******
*DAMR*
******
*univariate regression
reg damr chc_year_exp [aw = copop] 
predict damrhat
bys chc_year_exp: replace damrhat=. if _n>1
local bdu = floor(_b[chc_year_exp]*10)/10
local sedu = floor(_se[chc_year_exp]*10)/10

*multivariate
#delimit ;
reg damr	
_60pcturban _60pctrurf _60pct04years 
_60pctmt64years _60pctnonwhit _60pctmt12schl 
_60pctlt4schl _pct59inclt3k _pct59incmt10k 
_tot_act_md  chc_year_exp [aw = copop] , robust;
#delimit cr;
local bda = floor(_b[chc_year_exp]*10)/10
local seda = floor(_se[chc_year_exp]*10)/10

*get best-fit line for average x's and variable chc_year
preserve
foreach var of varlist _60pcturban _60pctrurf _60pct04years _60pctmt64years _60pctnonwhit _60pctmt12schl _60pctlt4schl _pct59inclt3k _pct59incmt10k _tot_act_md{
	egen testo = mean(`var')
	replace `var' = testo
	drop testo
}
predict damrhatadj
bys chc_year_exp: replace damrhatadj=. if _n>1
keep stfips cofips  damrhatadj
save "damrhatadj", replace
restore
merge 1:1 stfips cofips  using "damrhatadj"
drop _merge 


*MAKE GRAPHS
#delimit ;
twoway scatter amr chc_year_exp [aw=copop] ,
connect(n)										
msymbol(oh)								
legend(off)										
xtitle("")										
legend(off)
yaxis(1)
||
scatter amrhat amrhatadj chc_year_exp ,
connect(l l)								
msymbol(i T)									
mcolor( green red)				
lpattern( solid solid)							
lwidth( medthick medthick)			
lcolor(green red)							
legend(off)										
xtitle("")										
ytitle("Deaths per 100,000 Residents", size(large))
title("A. 1965 AMR", size(vlarge) color(black))
legend(off)
yaxis(1) saving("$output/panel_a.gph", replace);	

#delimit ;
twoway scatter damr chc_year_exp [aw=copop] ,
connect(n)										
cmissing(n) 							
msymbol(Oh)										
xtitle("")					
ytitle("Change in Deaths per 100,000 Residents", size(large))
title("B. 1960-1965 Change in AMR", size(vlarge) color(black))	
legend(off)
yaxis(1)
||
scatter damrhat damrhatadj chc_year_exp ,
connect(l l)																		
msymbol(i T)									
mcolor( green red)							
lpattern( solid solid)												
lcolor(green red)														
xtitle("")									
title("B. 1960-1965 Change in AMR", size(vlarge) color(black))
legend(order(- "Fitted Values: " 2 3) rows(1) label(2 "Univariate") label(3 "Multivariate") size(medium) region(style(none)))
yaxis(1) saving("$output/panel_b", replace);
#delimit cr;	

cd "$dofile"		
grc1leg "$output/panel_a.gph" "$output/panel_b.gph", legendfrom("$output/panel_b.gph") xsize(5) ysize(2.75) col(2) imargin(medium) 
graph display, xsize(5) ysize(2.75)
	
erase "$output/amrhatadj.dta"
erase "$output/damrhatadj.dta"
erase "$output/panel_a.gph"
erase "$output/panel_b.gph"	



/******************************************************
Main Event-Study Results for AMR
******************************************************/
clear all
*create a file to store estimates later on
set obs 33
gen time = _n - 8
save "$output/amr_chc_es_results", replace emptyok
use "$data/aer_data", clear
xtset 
*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	

*Make fixed effects
*urban categories by which to generate year-FE
cap drop _urb
cap drop Durb
egen _urb = total(D_60pcturban*(year==1960)), by(fips) /*D_60pcturban is 1960: Urban x Trend form raw data*/
egen Durb = cut(_urb), at(0, 1, 25, 50, 75, 110)		/*quarters with a zero*/

*make year dummies by urban category
xi i.year*i.Durb i.year*i.stfips 						/*this makes state-by-year FE and urban-by-year FE*/
cap drop _IDurb*
cap drop _Istfi*

*preferred specification: year FE, urban-by-year FE, state-by-year effects, 1960 char trends, REIS vars, AHA varls
local X "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _Texp* [aw=popwt]"	


*************************
*Early CHCs, All Counties
*************************
char exp1[omit] -1
xi i.exp1, pref(_T)
xtset fips year
xtreg amr `X' if year<=1988, cluster(fips) fe
testparm _Texp1_2-_Texp1_6
local ppre = r(p)
testparm _Texp1_8-_Texp1_22
local ppost = r(p)
*STORE RESULTS IN A STATA FILE
preserve
	use "$output/amr_chc_es_results", clear
	quietly{
		gen b_X_early_all				= .
		gen se_X_early_all				= .	
		forval h = 1/23{
				if `h'==7{
					replace b_X_early_all = 0 in `h'
				}
				else{
					replace b_X_early_all = _b[_Texp1_`h'] in `h'
					replace se_X_early_all = _se[_Texp1_`h'] in `h'
				}
		}
	}
	save "$output/amr_chc_es_results", replace
restore

*************************
*All CHCs, All Counties
*************************
local X "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _Texp* [aw=popwt]"	
char exp2[omit] -1
xi i.exp2, pref(_T)

xtreg amr `X' if year<=1988, cluster(fips) fe
testparm _Texp2_2-_Texp2_6
local ppre = r(p)
testparm _Texp2_8-_Texp2_16
local ppost = r(p)
*STORE RESULTS IN A STATA FILE
preserve
	use "$output/amr_chc_es_results", clear
	quietly{
		gen b_X_all_all				= .
		gen se_X_all_all				= .	
		forval h = 1/17{
				if `h'==7{
					replace b_X_all_all = 0 in `h'
				}
				else{
					replace b_X_all_all = _b[_Texp2_`h'] in `h'
					replace se_X_all_all = _se[_Texp2_`h'] in `h'
				}
		}
	}
	save "$output/amr_chc_es_results", replace
restore

**********************************
*Early CHCs, All Years (1959-1998)
**********************************
*Baltimore and Arlington can't be combined with counties after 1988
drop if fips==24510
drop if fips==51013

local X "_Iyear* _IyeaXDu* _IyeaXst* D_* _Texp* [aw=popwt]"	

char exp1_1998[omit] -1
xi i.exp1_1998, pref(_T)

xtreg amr `X' if samp8998, cluster(fips) fe
testparm _Texp1_1998_2-_Texp1_1998_6
local ppre = r(p)
testparm _Texp1_1998_8-_Texp1_1998_32
local ppost = r(p)
*STORE RESULTS IN A STATA FILE
preserve
	use "$output/amr_chc_es_results", clear
	quietly{
		gen b_X_early_5998				= .
		gen se_X_early_5998				= .	
		forval h = 1/33{
				if `h'==7{
					replace b_X_early_5998 = 0 in `h'
				}
				else{
					replace b_X_early_5998 = _b[_Texp1_1998_`h'] in `h'
					replace se_X_early_5998 = _se[_Texp1_1998_`h'] in `h'
				}
		}
	}
	save "$output/amr_chc_es_results", replace
restore

*Figure 3 Main results 

use "$output/amr_chc_es_results", clear


#delimit ;					
for var *all_all: replace X = . if time>8;
for var *early_all:  replace X = . if time>14;
for var *early_5998:  replace X = . if time>23;

*confidence intervals for early centers
cap drop ub* lb*;
gen ub 			= b_X_early_all + 1.96*se_X_early_all ;
gen lb 			= b_X_early_all - 1.96*se_X_early_all ;		

twoway (scatter b_X_all_all b_X_early_all b_X_early_5998 ub lb
		time if time>=-6 & time<=24,							
		xline(-1, lcolor(black)) 
		yline(0, lcolor(black)) 								
		connect(l l l l l l l l l)							
		msymbol(i i i i i i i)									
		mcolor(green navy red navy navy blue blue)	
		lpattern( solid solid solid dash dash dot dot)		
		lwidth( medthick medthick medthick thin thin medium medium medium medium)		
		lcolor(green  navy red navy navy blue blue)		
		legend(off)			
		xlabel(-6(3)24, labsize(medium))  					
		ylabel(,  labsize(medium))								
		xtitle("Years Since CHC Establishment", size(medium))
		ytitle("Deaths per 100,000 Residents" " ", size(medium))
		title("", size(medium) color(black))	
		text(7 9 "Year Before CHCs Began Operating", j(left)))
		(pcarrowi -6 6 -9 5 , lcolor(black) mcolor(black) lwidth(medthick) mlwidth(medthick)
		text(-5 8 "All CHCs, 1959-1988", j(left)))		
		(pcarrowi -12 14 -17 11 , lcolor(black) mcolor(black) lwidth(medthick) mlwidth(medthick)
		text(-10 18 "Early CHCs (funded 1965-1974);" "3,044 counties observed 1959-88", j(left)))		
		(pcarrowi -23 22 -21 21 , lcolor(black) mcolor(black) lwidth(medthick) mlwidth(medthick)
		text(-25 18 "Early CHCs (funded 1965-1974);" "388 counties observed 1959-98", j(left)))
		; 
		#delimit cr;
		
graph export "$output/figure5.wmf", replace

/******************************************************
PART 2. Heterogeneous treatment effects (we aplly estimator proposed by Sun, Liyang, and Sarah Abraham. "Estimating dynamic treatment effects in event studies with heterogeneous treatment effects." Journal of Econometrics (2020)
*******************************************************/

*---------------------------------------------------------------------
*---------------------------- this is for weights and it plots 3 graphs but we can build it for every cohort
/*Our cohort categorical variable based on when the first CHC was established in the country is 'chc_year_exp'
Our relative time varibales is 
'exp1' for the early centers (event-time, 1965-1974 CHCs);
'exp2' for the event-time all centers all counties;  
'exp1_1998' for Early CHCs, All Years (1959-1998)

We consider the time period from 1965 to 1998 for the relative time periods start from -7 to +14.
The specification includes leads =-7,-6,-5,-4,-3,-2, we exclude the lead=-1 following  Bailey & Goodman-Bacon (2014), and lags=0,1,2 to estimate the dynamic effect of CHC status on mortality.  We first generate these relative time indicators*/

*---------------------------------------------Weights for early CHCs -------------------------------
clear all
use "$data/aer_data", clear
xtset
cd "$output"



*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61					
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	

drop if year > 1988

*We will need to drop obs in -7, -1 and 15 to make the data balanced 
drop if exp1==-7 
drop if exp1==-1 
drop if exp1==15


**Weights for early CHCs 1965-74

gen g_6 = exp1 == -6
gen g_5 = exp1 == -5
gen g_4 = exp1 == -4
gen g_3 = exp1 == -3
gen g_2 = exp1 == -2
gen g0 = exp1 == 0
gen g1 = exp1 == 1
gen g2 = exp1 == 2
gen g3 = exp1 == 3
gen g4 = exp1 == 4 
gen g5 = exp1 == 5
gen g6 = exp1 == 6
gen g7= exp1 == 7
gen g8 = exp1 == 8
gen g9 = exp1 == 9
gen g10 = exp1 == 10
gen g11 = exp1 == 11
gen g12 = exp1 == 12
gen g13 = exp1 == 13
gen g14 = exp1 == 14


*For the coefficient associate with each of the above relative time indicators in a two-way fixed effects regression, we estimate the weights and export to a spreadsheet "weights.xlsx".

eventstudyweights  g_6 g_5 g_4 g_3 g_2 g0 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14  , controls(i.fips i.year) cohort(chc_year_exp) rel_time(exp1) saveweights("weights_early_65_74.xlsx") 
cap drop  g_7 g_6 g_5 g_4 g_3 g_2 g0 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15

*plot the weights 

import excel "$output/weights_early_65_74.xlsx", clear firstrow
preserve 
keep g_2 chc_year_exp exp1
reshape wide g_2, i(exp1) j(chc_year_exp)
egen w_sum = rowtotal(g_2*) /* this is just to check the properties of weights described in the paper*/
graph twoway line g_2* exp1, ytitle("Weights") xtitle("Relative wave") title("Estimated weights underlying {&mu} {subscript:-2}", size(medium) color(black)) legend(off) saving(g_2_early_65_74, replace)
	restore

	preserve 
keep g_5 chc_year_exp exp1
reshape wide g_5, i(exp1) j(chc_year_exp)
egen w_sum = rowtotal(g_5*) /* this is just to check the properties of weights described in the paper*/
graph twoway line g_5* exp1, ytitle("Weights") xtitle("Relative wave") title("Estimated weights underlying {&mu} {subscript:-5}", size(medium) color(black)) legend(off) saving(g_5_early_65_74, replace)
	restore


	preserve 
keep g4 chc_year_exp exp1
reshape wide g4, i(exp1) j(chc_year_exp)
egen w_sum = rowtotal(g4*) /* this is just to check the properties of weights described in the paper*/
graph twoway line g4* exp1, ytitle("Weights") xtitle("Relative wave") title("Estimated weights underlying {&mu} {subscript:4}", size(medium) color(black)) legend(off) saving(g4_early_65_74, replace)
	restore

		preserve 
keep g10 chc_year_exp exp1
reshape wide g10, i(exp1) j(chc_year_exp)
egen w_sum = rowtotal(g10*) /* this is just to check the properties of weights described in the paper*/
graph twoway line g10* exp1, ytitle("Weights") xtitle("Relative wave") title("Estimated weights underlying {&mu} {subscript:10}", size(medium) color(black)) legend(size(*0.5) cols(5) label(1 "1965") label(2	"1967") label (3 "1968") label (4 "1969") label (5 "1967") label (6 "1971") label (7 "1972") label (8 "1973") label (9 "1974")) saving(g10_early_65_74, replace) 
	restore
				
grc1leg "$output/g_2_early_65_74.gph" "$output/g_5_early_65_74.gph" "$output/g4_early_65_74.gph" "$output/g10_early_65_74.gph", legendfrom("$output/g10_early_65_74.gph") xsize(3) ysize(2) col(2) imargin(small) 

	
graph export "$output/weights_yearly_64_75.png", replace


*--------------------------------------Weights for all counties and all years -----------------------------

clear all
use "$data/aer_data", clear
xtset
cd "$output"

*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61					
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	

drop if year > 1988

drop if exp2==-1 
drop if exp2==-7
drop if exp2==9

**Weights for All CHCs all counties

gen g_6 = exp2 == -6
gen g_5 = exp2 == -5
gen g_4 = exp2 == -4
gen g_3 = exp2 == -3
gen g_2 = exp2 == -2
gen g0 = exp2 == 0
gen g1 = exp2 == 1
gen g2 = exp2 == 2
gen g3 = exp2 == 3
gen g4 = exp2 == 4 
gen g5 = exp2 == 5
gen g6 = exp2 == 6
gen g7= exp2 == 7
gen g8 = exp2 == 8



eventstudyweights  g_6 g_5 g_4 g_3 g_2 g0 g1 g2 g3 g4 g5 g6 g7 g8 , controls(i.fips i.year) cohort(chc_year_exp) rel_time(exp2) saveweights("weights_all_counties") 
drop  g_6 g_5 g_4 g_3 g_2 g0 g1 g2 g3 g4 g5 g6 g7 g8 

**-------------------------------------Weights for Early CHCs, All Years (1959-1998)-------------------------
clear all
use "$data/aer_data", clear
xtset
cd "$output"
*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61					
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	
*Baltimore and Arlington can't be combined with counties after 1988
drop if fips==24510
drop if fips==51013

drop if exp1_1998==-1 
drop if exp1_1998==-7
drop if exp1_1998==25 


gen g_6 = exp1_1998 == -6
gen g_5 = exp1_1998 == -5
gen g_4 = exp1_1998 == -4
gen g_3 = exp1_1998 == -3
gen g_2 = exp1_1998 == -2
gen g0 = exp1_1998 == 0
gen g1 = exp1_1998 == 1
gen g2 = exp1_1998 == 2
gen g3 = exp1_1998 == 3
gen g4 = exp1_1998 == 4 
gen g5 = exp1_1998 == 5
gen g6 = exp1_1998 == 6
gen g7= exp1_1998 == 7
gen g8 = exp1_1998 == 8
gen g9 = exp1_1998 == 9
gen g10 = exp1_1998 == 10
gen g11 = exp1_1998 == 11
gen g12 = exp1_1998 == 12
gen g13 = exp1_1998 == 13
gen g14 = exp1_1998 == 14
gen g15 = exp1_1998 == 15
gen g16 = exp1_1998 == 16 
gen g17 = exp1_1998 == 17
gen g18 = exp1_1998 == 18
gen g19 = exp1_1998 == 19
gen g20 = exp1_1998 == 20
gen g21 = exp1_1998 == 21
gen g22 = exp1_1998 == 22
gen g23 = exp1_1998 == 23
gen g24 = exp1_1998 == 24



	eventstudyweights  g_6 g_5 g_4 g_3 g_2 g0 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15 g16 g17 g18 g19 g20 g21 g22 g23 g24 , controls(i.fips i.year) cohort(chc_year_exp) rel_time(exp1_1998) saveweights("weights_early_all_years")
	drop   g_6 g_5 g_4 g_3 g_2 g0 g1 g2 g3 g4 g5 g6 g7 g8 g9 g10 g11 g12 g13 g14 g15 g16 g17 g18 g19 g20 g21 g22 g23 g24 

	
	

*------------------------------------------NEW ESTIMATOR---- Early CHCs, All Counties---------------------

clear all

set obs 180
gen time = _n
save "$output/CATT Early CHCs All Countries", replace emptyok

use "$data/aer_data", clear
xtset
cd "$output"
*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61					
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	


drop if exp1==-7 
drop if exp1==-1 
drop if exp1==15


* Step 1 following Sun & Abraham (2020): estimate the the CATT using the interacted specification (relative time with treatment dummy)

*preferred specification: year FE, urban-by-year FE, state-by-year effects, 1960 char trends, REIS vars, AHA varls
local X "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _Ttreat* [aw=popwt]"	

*Make fixed effects
*urban categories by which to generate year-FE
cap drop _urb
cap drop Durb
egen _urb = total(D_60pcturban*(year==1960)), by(fips)
egen Durb = cut(_urb), at(0, 1, 25, 50, 75, 110)	

*make year dummies by urban category
xi i.year*i.Durb i.year*i.stfips 				
cap drop _IDurb*
cap drop _Istfi*

*Generate interacted CATT treatment variables:
gen treat = exp1*chc_year_exp
*We will need to drop obs in -7, -1 and 15 to make the data balanced 



*CATT Early CHCs, All Counties
xi i.treat, pref(_T)
xtreg amr `X' if year<=1988, cluster(fips) fe

*STORE RESULTS IN A STATA FILE
preserve
	use "$output/CATT Early CHCs All Countries", clear
	quietly{
		gen b_X_early_all				=0.
	forval h = 2/172{
					replace b_X_early_all = _b[_Ttreat_`h'] in `h'
					}
		}
	 save "$output/Early CHCs All Countries", replace
restore


****Match the coefficients with weights


import excel "$output\weights_early_65_74.xlsx", clear firstrow
g time=_n
save weights_early_65_74.dta, replace

use "$output/Early CHCs All Countries", clear
merge 1:1  time using  weights_early_65_74
drop _merge

******Calculate the updated estimates 

foreach var of varlist g* {
	g new_b_`var'=`var'*b_X_early_all
} 

collapse (sum) new_b* 
xpose, clear

rename v1  b_X_early_all
g newb_X_early_all=b_X_early_all/_N
drop b_X_early_all
g time=_n-7 if _n<6
replace time=_n-6 if time==. 
save newbeta_early_65_74.dta, replace


*----------------------------------------------------NEW ESTIMATOR ---- All CHCs, All Counties----------------


clear all
set obs 196
gen time = _n 
save "$output/CATT all all", replace emptyok

use "$data/aer_data", clear
xtset
cd "$output"
*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61					
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	

drop if exp2==-1 
drop if exp2==-7
drop if exp2==9

*preferred specification: year FE, urban-by-year FE, state-by-year effects, 1960 char trends, REIS vars, AHA varls
local X "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _Ttreat* [aw=popwt]"	

*Make fixed effects
*urban categories by which to generate year-FE
cap drop _urb
cap drop Durb
egen _urb = total(D_60pcturban*(year==1960)), by(fips)
egen Durb = cut(_urb), at(0, 1, 25, 50, 75, 110)	

*make year dummies by urban category
xi i.year*i.Durb i.year*i.stfips 				
cap drop _IDurb*
cap drop _Istfi*

*Generate interacted CATT treatment variables:
gen treat = exp2*chc_year_exp
*We will need to drop obs in -7, -1 and 9 to make the data balanced 



*CATT Early CHCs, All Counties

xi i.treat, pref(_T)
xtreg amr `X' if year<=1988, cluster(fips) fe


*STORE RESULTS IN A STATA FILE
preserve
	use "$output/CATT all all", clear
	quietly{
		gen b_X_all_all				= 0.
		forval h = 2/196{
			replace b_X_all_all = _b[_Ttreat_`h'] in `h'
					}
		}
	save "$output/CATT all all", replace
restore


****Match the coefficients with weights


import excel "$output\weights_all_counties.xlsx", clear firstrow
gene time=_n
save weights_early_all_years.dta, replace

use "$output/CATT all all", clear
merge 1:1  time using  weights_early_all_years
drop _merge


******Calculate the updated estimates 

foreach var of varlist g* {
	g new_b_`var'=`var'*b_X_all_all
} 

collapse (sum) new_b* 
xpose, clear
rename v1  b_X_all_all
g newb_X_all_all=b_X_all_all/_N
drop b_X_all_all
g time=_n-7 if _n<6
replace time=_n-6 if time==. 
save newbeta_all_all.dta, replace




*--------------------------------------------NEW Estimator ----------Early CHCs, All Years (1959-1998)

clear all
set obs 262
gen time = _n
save "$output/CATT Early CHCs all years", replace 


use "$data/aer_data", clear
xtset
cd "$output"


*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61					
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	

*Baltimore and Arlington can't be combined with counties after 1988
drop if fips==24510
drop if fips==51013

drop if exp1_1998==-1 
drop if exp1_1998==-7
drop if exp1_1998==25 

*preferred specification: year FE, urban-by-year FE, state-by-year effects, 1960 char trends, REIS vars, AHA varls
local X "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _Ttreat* [aw=popwt]"	

*Make fixed effects
*urban categories by which to generate year-FE
cap drop _urb
cap drop Durb
egen _urb = total(D_60pcturban*(year==1960)), by(fips)
egen Durb = cut(_urb), at(0, 1, 25, 50, 75, 110)	

*make year dummies by urban category
xi i.year*i.Durb i.year*i.stfips 				
cap drop _IDurb*
cap drop _Istfi*

*Generate interacted CATT treatment variables:
gen treat = exp1_1998*chc_year_exp
*We will need to drop obs in -7, -1 and 15 to make the data balanced 



*CATT Early CHCs, All Counties

xi i.treat, pref(_T)
xtreg amr `X' if samp8998, cluster(fips) fe



*STORE RESULTS IN A STATA FILE
preserve
	use "$output/CATT Early CHCs all years", clear
	quietly{
		gen b_X_early_5998				= 0.
		forval h = 2/262{
			replace b_X_early_5998 = _b[_Ttreat_`h'] in `h'
				}
		}
	
	save "$output/CATT Early CHCs all years", replace
restore


****Match the coefficients with weights


import excel "$output\weights_early_all_years.xlsx", clear firstrow
gene time=_n
save weights_early_all_years.dta, replace

use "$output/CATT Early CHCs all years", clear
merge 1:1  time using  weights_early_all_years
drop _merge


******Calculate the updated estimates 

foreach var of varlist g* {
	g new_b_`var'=`var'*b_X_early_5998 

} 

collapse (sum) new_b* 
xpose, clear
rename v1  b_X_early_5998
g newb_X_early_5998=b_X_early_5998/_N
drop b_X_early_5998
g time=_n-7 if _n<6
replace time=_n-6 if time==. 
save newbeta_early_5998.dta, replace

*---------------------------------make a table to compare coefficients----------------------

use "$output/amr_chc_es_results", clear
drop se*
merge 1:1  time using newbeta_early_65_74
drop _merge
merge 1:1  time using newbeta_all_all
drop _merge
merge 1:1  time using newbeta_early_5998
drop _merge
order time b_X_early_all  newb_X_early_all  b_X_all_all  newb_X_all_all  b_X_early_5998  newb_X_early_5998 

