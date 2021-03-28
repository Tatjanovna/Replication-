/******************************************************************************
Replication for:

The War on Poverty's Experiment in Public Medicine: 
The Impact of Community Health Centers on the Mortality of Older Americans
by Martha Bailey and Andrew Goodman-Bacon

Date: 10/03/2021
******************************************************************************/
* Set up the working directories

*dofile directory (where this file is stored)
global dofile "C:\Users\Corrine\Documents\博士课件\PhD Econometrics 2\Replication--mainReplication\Replication"

*data directory (where posted datasets are stored)
global data "C:\Users\Corrine\Documents\博士课件\PhD Econometrics 2\112871-V1\aer_data"

*output directory (where regression output, figures, and logs are saved)
global output "C:\Users\Corrine\Documents\博士课件\PhD Econometrics 2\Replication--main\Output files"


cd "$output"


/***********************************************
FIGURE 2. AMR by Age Group, 1959-1988
***********************************************/


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
msymbol(i i)							/*add lines*/ 
ylabel(#4, labsize(small))
xlabel(1959(5)1984 1988, labsize(small))
xtitle("Year")
ytitle("Deaths per 100,000 Residents")			
title("{it: A. Age-Adjusted Mortality}", size(medmsall) color(black))
graphregion(fcolor(white) color(white) icolor(white)) saving(panel1, replace);
#delimit cr;	

*Infant Mortality
#delimit ;
scatter imr nnmr pnmr year,
connect(l l l l l)
msymbol(i Oh S)		
ylabel(#4, labsize(small))
xlabel(1959(5)1984 1988, labsize(small))
legend(off)
text(27 1964 "Total", size(small) place(e) yaxis(1)) 
text(14 1959 "Neonatal", size(small) place(e) yaxis(1)) 
text(3 1959 "Post-Neonatal", size(small) place(e) yaxis(1)) 
xtitle("Year")
ytitle("Deaths per 1,000 Live Births")
title("{it: B. Infants}", size(medmsall) color(black))
graphregion(fcolor(white) color(white) icolor(white)) saving(panel2, replace);
#delimit cr;	

*Child and Adult Mortality
#delimit ;
scatter amr_ch amr_ad year,
connect(l l l l l)	
msymbol(i i)
ylabel(#4, labsize(small))
xlabel(1959(5)1984 1988, labsize(small))
legend(off)
xtitle("Year")
ytitle("Deaths per 100,000 Residents")
title("{it: C. Children (1-19) and Adults (20-49)}", size(medmsall) color(black))
text(250 1980 "Adults", size(medsmall))
text(100 1965 "Children", size(medsmall))			
graphregion(fcolor(white) color(white) icolor(white)) saving(panel3, replace);
#delimit cr;	

*Older Adult Mortality
#delimit ;
scatter amr_eld year,
connect(l l l l l)	
msymbol(i i)
ylabel(#4, labsize(small))
xlabel(1959(5)1984 1988, labsize(small))
legend(off)
xtitle("Year")
ytitle("Deaths per 100,000 Residents")
title("{it: D. Older Adults (50+)}", size(medmsall) color(black))
graphregion(fcolor(white) color(white) icolor(white)) saving(panel4, replace);
#delimit cr;	


graph combine panel1.gph panel2.gph panel3.gph panel4.gph, col(2) imargin(tiny) xcommon xsize(8.5) ysize(5.5) graphregion(fcolor(white) color(white) icolor(white))
graph display, xsize(8.5) ysize(5.5)

graph export "$output/figure2.wmf", replace


/*********************************************************************
FIGURE 4: Predicting the timing of CHC grant with 1960 Characteristics
*********************************************************************/

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
connect(n)										/*LINES: add lines to the fitted values, but not to the ASMR scatters*/
cmissing(n) 							/*MISSING: do not connect lines between missing points (ie. confidence intervals)*/
msymbol(Oh)									/*MARKER SYMBOL: O is closed circle, T is closed triangle, Oh and Th are open circle and triangle*/
xlabel(1965(3)1974, labsize(large))				/*XLABEL: this makes the axis tight [it doesn't work with xscale(range()), only xlabel]*/
ylabel(200(400)1400, labsize(large) axis(1))						/*YLABEL: make there be 4 ticks so they don't overlap and make the tick-mark text small*/							
graphregion(fcolor(white) color(white) icolor(white)) 
yaxis(1)
||
scatter amrhat amrhatadj chc_year_exp ,
connect(l l)										/*LINES: add lines to the fitted values, but not to the ASMR scatters*/
cmissing(y y)								/*MISSING: do not connect lines between missing points (ie. confidence intervals)*/
msymbol(i T)									/*MARKER SYMBOL: O is closed circle, T is closed triangle, Oh and Th are open circle and triangle*/
msize( . medlarge)
mcolor( forest_green maroon)								/*MARKER COLOR:*/ 
lpattern( solid solid)							/*LINE PATTERN:*/
lwidth( medthick medthick)						/*LINE WIDTHS:*/
lcolor(forest_green maroon)						
xlabel(1965(3)1974, labsize(large))				/*XLABEL: this makes the axis tight [it doesn't work with xscale(range()), only xlabel]*/
ylabel(200(400)1400, labsize(large) axis(1))						/*YLABEL: make there be 4 ticks so they don't overlap and make the tick-mark text small*/
xtitle("")										/*XTITLE: no title on this one, only on the bottom panels to save vertical space*/
ytitle("Deaths per 100,000 Residents", size(large))
title("{it: A. 1965 AMR}", size(vlarge) color(black))	/*TITLE: panel titles are defined in the locals above*/
legend(order(- "Fitted Values: " 2 3) rows(1) label(2 "Univariate") label(3 "Multivariate") size(medium) region(style(none)))
graphregion(fcolor(white) color(white) icolor(white)) 
yaxis(1) saving("$output/panel_a.gph", replace);


#delimit ;
twoway scatter damr chc_year_exp [aw=copop] ,
connect(n)										/*LINES: add lines to the fitted values, but not to the ASMR scatters*/
cmissing(n) 							/*MISSING: do not connect lines between missing points (ie. confidence intervals)*/
msymbol(Oh)									/*MARKER SYMBOL: O is closed circle, T is closed triangle, Oh and Th are open circle and triangle*/
							
xlabel(1965(3)1974, labsize(large))				/*XLABEL: this makes the axis tight [it doesn't work with xscale(range()), only xlabel]*/
ylabel(-400(200)200, labsize(large) axis(1))						/*YLABEL: make there be 4 ticks so they don't overlap and make the tick-mark text small*/
xtitle("")										/*XTITLE: no title on this one, only on the bottom panels to save vertical space*/
ytitle("Change in Deaths per 100,000 Residents", size(large))
title("{it: B. 1960-1965 Change in AMR}", size(vlarge) color(black))	/*TITLE: panel titles are defined in the locals above*/

graphregion(fcolor(white) color(white) icolor(white)) 
yaxis(1)
||
scatter damrhat damrhatadj chc_year_exp ,
connect(l l)										/*LINES: add lines to the fitted values, but not to the ASMR scatters*/
cmissing(y y)								/*MISSING: do not connect lines between missing points (ie. confidence intervals)*/
msymbol(i T)									/*MARKER SYMBOL: O is closed circle, T is closed triangle, Oh and Th are open circle and triangle*/
msize( . medlarge)
mcolor( forest_green maroon)								/*MARKER COLOR:*/ 
lpattern( solid solid)							/*LINE PATTERN:*/
lwidth( medthick medthick)						/*LINE WIDTHS:*/
lcolor(forest_green maroon)									/*LINE COLOR:*/ 
xlabel(1965(3)1974, labsize(large))				/*XLABEL: this makes the axis tight [it doesn't work with xscale(range()), only xlabel]*/
ylabel(-400(200)200, labsize(large) axis(1))						/*YLABEL: make there be 4 ticks so they don't overlap and make the tick-mark text small*/
xtitle("")										/*XTITLE: no title on this one, only on the bottom panels to save vertical space*/
title("{it: B. 1960-1965 Change in AMR}", size(vlarge) color(black))	/*TITLE: panel titles are defined in the locals above*/
legend(order(- "Fitted Values: " 2 3) rows(1) label(2 "Univariate") label(3 "Multivariate") size(medium) region(style(none)))
graphregion(fcolor(white) color(white) icolor(white)) 
yaxis(1) saving("$output/panel_b", replace);
#delimit cr;	

		
graph combine panel_a.gph panel_b.gph
graph display, xsize(8.5) ysize(4)
graph export "$output/figure4.wmf", replace

	



/******************************************************
FIGURE 5: Main Event-Study Results for AMR
******************************************************/

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
egen _urb = total(D_60pcturban*(year==1960)), by(fips)
egen Durb = cut(_urb), at(0, 1, 25, 50, 75, 110)		/*quarters with a zero*/

*make year dummies by urban category
xi i.year*i.Durb i.year*i.stfips 						/*this makes state-by-year FE and urban-by-year FE*/
cap drop _IDurb*
cap drop _Istfi*

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

**********************************
*Early CHCs, All Years (1959-1998)
**********************************
*Baltimore and Arlington can't be combined with counties after 1988
drop if fips==24510
drop if fips==51013

*can't include AHA or REIS covariates in the later years
local X "_Iyear* _IyeaXDu* _IyeaXst* D_* _Texp* [aw=popwt]"	

char exp1_1998[omit] -1
xi i.exp1_1998, pref(_T)

xtreg amr `X' if samp8998, cluster(fips) fe
testparm _Texp1_1998_2-_Texp1_1998_6
local ppre = r(p)
testparm _Texp1_1998_8-_Texp1_1998_32
local ppost = r(p)

		
#delimit ;					
for var *all_all: replace X = . if time>8;
for var *early_all:  replace X = . if time>14;
for var *early_5998:  replace X = . if time>23;

cap drop ub* lb*;
gen ub 			= b_X_early_all + 1.96*se_X_early_all ;
gen lb 			= b_X_early_all - 1.96*se_X_early_all ;		

twoway (scatter b_X_all_all b_X_early_all b_X_early_5998 ub lb
		time if time>=-6 & time<=24,							/*SCATTER: allows markers, and is actually equivalent to "line" with the connect option*/
		xline(-1, lcolor(black)) 								/*XLINE: refers to a vertical line that crosses the x-axis. Put one at the omitted category, -1, so all "pre" periods have NO CHC*/
		yline(0, lcolor(black)) 								/*YLINE: horizontal black line at 0*/
		connect(l l l l l l l l l)								/*LINES: add lines (of course)*/
		cmissing(n n n n n n n n n)								/*MISSING: do not connect lines between missing points (ie. confidence intervals)*/
		msymbol(O i Th i i i i)									/*MARKER SYMBOL: Th is open triangle, X are x's, and make the preferred specification be a solid line*/
		msize(medium . medium medium medium . . . .)			/*MARKER SIZE: only matters when a marker is specified*/
		mcolor(maroon navy forest_green navy navy blue blue)	/*MARKER COLOR: */
		lpattern( solid solid solid dash dash dot dot)			/*LINE PATTERN: estimates solid, CI dashed*/
		lwidth( medthick vthick medthick medthick medthick medium medium medium medium)			/*LINE WIDTHS: make preferred specification thick*/
		lcolor(maroon navy forest_green navy navy blue blue)	/*LINE COLOR: main results are always navy*/ 		
		legend(off)				
		xlabel(-6(3)24, labsize(medium))  						/*XLABEL: this makes the axis tight [it doesn't work with xscale(range()), only xlabel]*/
		ylabel(,  labsize(medium))								/*YLABEL*/
		xtitle("Years Since CHC Establishment", size(medium))	/*XTITLE*/
		ytitle("Deaths per 100,000 Residents" " ", size(medium))
		title("", size(medium) color(black))	
		graphregion(fcolor(white) color(white) icolor(white) margin(small))  /*BACKGROUNDS: this code takes away all the borders and color from the stata graph background and maximizes graph size within its region*/
		plotregion(margin(vsmall)))
		(pcarrowi 7 2 7 -.5, lcolor(black) mcolor(black) lwidth(medthick) mlwidth(medthick)
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
TABLE 1: Mean 1960 Characteristics 
******************************************************/

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

*RESCALING HAS A HUGE EFFECT ON THE T-STAT!  IN MY ORIGINAL CODE I DIDN'T DO A GOOD JOB RESCALING THINGS
	
	
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


/******************************************************
TABLE 2: DD Estimates by Specification
******************************************************/
clear
clear matrix
set mat 5000
clear mata
set maxvar 10000
set matsize 10000
set more off, perm
pause on
capture log close
log using "$output/log_table2", replace text

use "$data/aer_data" if year<=1988, clear
*Drop NY/LA/Chicago
drop if stfips==36 & cofips==61
drop if stfips==6  & cofips==37
drop if stfips==17 & cofips==31	

*Make fixed effects
*urban categories by which to generate year-FE
cap drop _urb
cap drop Durb
egen _urb = total(D_60pcturban*(year==1960)), by(fips)
egen Durb = cut(_urb), at(0, 1, 25, 50, 75, 110)		/*quarters with a zero*/

*make year dummies by urban category, state FE and county trends
xi	i.year*i.Durb i.fips*year i.year*i.stfips			/*this makes county trends, state-by-year FE and urban-by-year FE*/
cap drop _Ifips*
cap drop _IDurb*
cap drop _Istfi*	

*TREATMENT VARS
char did1[omit] -1
xi i.did1, pref(_DD)
		
		
*PANEL A		
*define the specifications
local X1 "_Iyear* _IyeaXDu* _DD* [aw=popwt]"										
local X2 "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _DD* [aw=popwt]"	
local X3 "_Iyear* _IyeaXDu* _IyeaXst* _IfipXy* R_* H_* _DD* [aw=popwt]"								
local X4 "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _DD* [aw=dflpopwgt1]"								

*for display in the log file
local DisX1 "FE + UxY FE"	
local DisX2 "FE + UxY FE + X + STxY FE"
local DisX3 "FE + UxY FE + X + STxY FE + County Trends"	
local DisX4 "FE + UxY FE + X + STxY FE, DFL Weights"

local r replace
foreach i of numlist 1 2 3 4{		
	di "amr: `DisX`i''"
	xtreg amr  `X`i'', cluster(fips) fe
			
	*OUTREG RESULTS
	outreg2 using "$output/table2.xls", `r' keep(_DDdid1_2 _DDdid1_4 _DDdid1_5 _DDdid1_6) noparen noaster ctitle("AMR: `e(cmdline)'") 
	local r append
}


*PANEL B
*define the specifications
local X1 "_Iyear* _IyeaXDu* _DD* [aw=popwt_eld]"										
local X2 "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _DD* [aw=popwt_eld]"	
local X3 "_Iyear* _IyeaXDu* _IyeaXst* _IfipXy* R_* H_* _DD* [aw=popwt_eld]"								
local X4 "_Iyear* _IyeaXDu* _IyeaXst* D_* R_* H_* _DD* [aw=dflpopwgt1_eld]"								

*for display in the log file
local DisX1 "FE + UxY FE"	
local DisX2 "FE + UxY FE + X + STxY FE"
local DisX3 "FE + UxY FE + X + STxY FE + County Trends"	
local DisX4 "FE + UxY FE + X + STxY FE, DFL Weights"

*adjust hospital variables to be per-50+ year old
replace H_bpc = H_bpc*(copop/copop_eld)
replace H_hpc = H_hpc*(copop/copop_eld)

foreach i of numlist 1 2 3 4{		
	di "amr_eld: `DisX`i''"
	xtreg amr_eld `X`i'', cluster(fips) fe
			
	*OUTREG RESULTS
	outreg2 using "$output/table2.xls", `r' keep(_DDdid1_2 _DDdid1_4 _DDdid1_5 _DDdid1_6) noparen noaster ctitle("AMR 50+: `e(cmdline)'") 
	local r append
}

log close
