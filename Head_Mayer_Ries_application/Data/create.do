# delimit ;
clear all;
cap log close;
cap file close _all;
cap program drop _all;
set scrollbufsize 1000000;
set more off;
log using create.log, replace;


*******************************************;
* Head, Mayer & Ries (2010) data;
* Montiel Olea & Plagborg-Moller 2017-03-25;
*******************************************;

*** SETTINGS ***;

local do_export = 1; // =1: export to CSV for use in Matlab;

*** END SETTINGS ***;


* Form data set;

use col_regfile09;

egen dyad = group(iso_o iso_d);
xtset dyad year;
order dyad year iso_o iso_d;


* Create variables;

foreach v of varlist pop_? gdpcap_? distw {;
	gen `v'_log = log(`v');
};

gen flow_roundlog = log(round(flow, 0.01));
gen byte gatt_both = gatt_o*gatt_d;
gen byte col_a = missing(indepdate) & (col_h==1);

gen yrsindep = year - indepdate;
replace yrsindep = 0 if missing(yrsindep) | yrsindep<0;
gen yrsindepcap = cond(yrsindep>=60, 60, yrsindep);


* Regressions;

local lhs = "flow_roundlog";
local rhs = "pop_?_log gdpcap_?_log distw_log contig comlang_off comleg col_hist col_a rta gatt_both comcur acp_to_eu";
local fe = "yrsindepcap year";

local fe_i = "";
foreach v of varlist `fe' {;
	local fe_i = "`fe_i'i.`v' ";
};

reg `lhs' `rhs' `fe_i', vce(cluster dyad);
gen reg_sample = e(sample);


* Export;

if `do_export'==1 {;

	egen dyad_id = group(dyad) if reg_sample==1;

	export delim dyad_id `lhs' `fe' `rhs' if reg_sample==1 using gravity.csv, replace;

};


log close;

