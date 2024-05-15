******************************************************************;
*
* PEV_FRAILTY
* ===========
*
* SAS-macro for computation of Proportion of Explained Variation
* (PEV) for models with frailty term
*
* Gleiss, A. & Schemper, M. Explained variation in shared frailty 
* models,
* Statistics in Medicine 2018, 37:1482 - 1490
*
* Author:  Andreas Gleiss
* Version: 1.7 (maddala option) 
* Date:    14 Feb 2019
*
* Macro parameters:
* =================
* 
* data		name of SAS data set
* time		survival time variable 
* status	survival status variable
* varlist	list of independent variables for Cox model
* frailty	variable to be used as random effect (e.g. center)
* modopt	options for model statement in model with random 
*			statement (e.g., %str(/ maxiter=100));
* censval	value of status indicating censored observations 
*			(default: 0)
* dist		distribution of random effect:
*			lognormal (default) or gamma
* odssel	control output from phreg (default: none)
* print 	=1 to print results (default), =0 (e.g., for simulations)
* vmax		=1 to limit PEV estimates to non-negative values
* by		by processing for simulations (default: _all)
*			Caution: dist=gamma with by statement gives wrong 
*					 results for EV and CovParms
* margonly	=1 to calculate only marginal PEV (default: 0)
* maddala	=1 to calculate likelihood-based Maddala-R-squared
*
******************************************************************;

%macro pev_frailty(data,time,status,varlist,frailty,modopt,censval=0,
					dist=lognormal,odssel=none,print=1,vmax=1,by=_all,
					margonly=0, maddala=0);

%let nvar=0;
%do %while(%scan(&varlist,&nvar+1)~=);
	%let nvar=%eval(&nvar+1);
	%let var&nvar=%scan(&varlist,&nvar);
	%end;
data _data;
	set &data;
	_all=1;
	run;

*** Calculate PEV for model without frailty;

%if &margonly=0 %then %do;
	ods select &odssel;
	proc phreg data=_data ev;
		*class &class;
		model &time*&status(&censval)=&varlist;
		ods output NObs=_nobs0 GlobalTests=_likel0 ExplainedVariation=_ev0;
		by &by;
		run;
	data _r2m0;
		merge _likel0(where=(test="Likelihood Ratio"))
			  _nobs0(keep=nobsused);
		*by &by;
		R2M=100*(1-(exp(-chisq))**(1/nobsused));
		run; * s. Schemper (StatMed, 1993, S. 2378);
	%end;


*** Calculate marginal PEV for frailty;

* output random effect estimates;
ods select &odssel;
proc phreg data=_data; * ev option not possible with random statement;
	class &frailty;
	model &time*&status(&censval)= &modopt;
	random &frailty / solution dist=&dist;
	ods output NObs=_nobs1 GlobalTests=_likel1 SolutionR=est_frailty0 CovParms=_gest1; 
	by &by;
	run;
ods select all;
* merge with original data;
proc sort data=_data;
	by &by &frailty;
	run;
data offset_mod0;
	merge _data 
		  est_frailty0(keep=&frailty estimate &by
					   rename=(estimate=frailty_est));
		by &by &frailty;
	_dummy=1;
	run;
proc standard data=offset_mod0 out=offset_mod0_std mean=0;
	var frailty_est;
	by &by;
	run; * necessary to standardize offset "by hand" before EV calculation;

* input frailty estimates as offset => estimate marginal PEV;
ods select &odssel;
proc phreg data=offset_mod0_std ev;
	model &time*&status(&censval)=_dummy /* SAS needs a dummy in addition to the offset */
		/ offset=frailty_est;
	ods output ExplainedVariation=_ev1;
	by &by;
	run; 
data _r2m1;
	merge _likel1(where=(test="Likelihood Ratio"))
		  _nobs1(keep=nobsused);
	*by &by;
	R2M=100*(1-(exp(-chisq))**(1/nobsused));
	run; * s. Schemper (StatMed, 1993, S. 2378);


*** Calculate PEV for model including frailty;

%if &margonly=0 %then %do;
	* output random effect estimates;
	ods select &odssel;
	proc phreg data=_data outest=est; * ev option not possible with random statement;
		class /*&class*/ &frailty;
		model &time*&status(&censval)=&varlist &modopt;
		random &frailty / solution dist=&dist;
		ods output NObs=_nobs2 GlobalTests=_likel2 SolutionR=est_frailty CovParms=_gest2;
		by &by;
		run;
	ods select all;
	* merge with original data;
	proc sort data=_data;
		by &by &frailty;
		run;
	data offset_mod;
		length _type_ $8.;
		merge _data 
			  est_frailty(keep=&frailty estimate &by
						  rename=(estimate=frailty_est));
			by &by &frailty;
		_dummy=1;
		_type_="PARMS"; * for next merge;
		run;
	data offset_mod;
		merge offset_mod est(keep=_type_ &varlist &by
							 rename=(
								%do iv=1 %to &nvar;
									&&var&iv=beta_&&var&iv
									%end; ));
		by &by _type_;
		xbeta=frailty_est
			%do iv=1 %to &nvar;
				+beta_&&var&iv*&&var&iv
				%end;;
		run;
	proc standard data=offset_mod out=offset_mod_std mean=0;
		var frailty_est xbeta;
		by &by;
		run; * necessary to standardize offset "by hand" before EV calculation;

	* input linear predictor estimates as offset => estimate PEV;
	ods select &odssel;
	proc phreg data=offset_mod_std ev;
		model &time*&status(&censval)=_dummy /* SAS needs a dummy in addition to the offset */
			/ offset=xbeta;
		ods output ExplainedVariation=_ev2;
		by &by;
		run; 
	data _r2m2;
		merge _likel2(where=(test="Likelihood Ratio"))
			  _nobs2(keep=nobsused);
		*by &by;
		R2M=100*(1-(exp(-chisq))**(1/nobsused));
		run; * s. Schemper (StatMed, 1993, S. 2378);

	ods select all;
	data _ev2;
		merge _ev2 _gest2(keep=estimate rename=(estimate=g_hat));
		run;
	%end;

data _ev1;
	merge _ev1 _gest1(keep=estimate rename=(estimate=g_hat));
	run;

%if &margonly %then %do;
	data _ev_all;
		set _ev1;
		model="Random effect";
		%if &vmax=1 %then %do;
			v=max(0,v);
			%end;
		run;
	data _r2m_all;
		set _r2m1;
		model="Random effect";
		run;
	%end;
%else %do;
	data _ev_all;
		set _ev0(in=zero) _ev1(in=one) _ev2(in=two);
		if zero then model="Fixed effects";
		if one  then model="Random effect";
		if two  then model="Full model";
		%if &vmax=1 %then %do;
			v=max(0,v);
			%end;
		run;
	data _r2m_all;
		set _r2m0(in=zero) _r2m1(in=one) _r2m2(in=two);
		if zero then model="Fixed effects";
		if one  then model="Random effect";
		if two  then model="Full model";
		run;
	%end;
%if &print=1 %then %do;
	title "Marginal Proportions of Explained Variation";
	proc print data=_ev_all noobs label;
		var model 
			%if &by~=_all %then %do;
				&by
				%end;
			d dx v;
		label model="Model";
		run;
	%if &maddala=1 %then %do;
		title "Likelihood-based R-square by Maddala";
		proc print data=_r2m_all noobs label;
			var model 
				%if &by~=_all %then %do;
					&by
					%end;
				R2M;
			label model="Model" R2M="R-squared (Maddala)";
			format R2M f5.2;
			run;
		%end;
	title;
	%end;
%mend;

*%pev_frailty(data=fpev.osdata_ct, time=surt,status=surs,censval=0,
			 varlist=lev inf, frailty=ctnum_, odssel=all);*, dist=gamma);

*%pev_frailty(data=fpev.osdata_ct, time=surt,status=surs,censval=0,
			 varlist=age sex tumor pos_lk_3 grad tu_lo2 tu_lo3 lev inf,
			 frailty=ctnum_);

