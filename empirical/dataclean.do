clear all

import delimited "[Your directory]\filename.csv"
// import delimited "" // Put your CSV file directory and name

// 1. Data from 2010 to 2019

keep if sic<1000 | (sic>=1000&sic<=1500) | (sic>=1500&sic<=1799) | (sic>=2000&sic<=3999) | (sic>=5000&sic<=5199) | (sic>=5200&sic<=5999) | (sic>=7000&sic<=8999)

keep if (fyear>=2009 & fyear<=2019)
keep if fic=="USA"

drop if gvkey==66552 // Duplicates, anyway data of this firm is missing
drop if sale <=0 | capx <= 0 | ppent <= 0

keep gvkey fyear sale at csho prcc_f ceq dlc dltt oibdp capx ppent dp ib che sic
foreach v of var * { 
	drop if missing(`v') // Drop any observations with missing value
}

xtset gvkey fyear
bysort gvkey: drop if _N<11 // making balanced panel



// Adjust for inflation (CPI)

matrix cpi = (100, 101.63657, 104.82760, 107.00088, 108.56948, 110.32338, 110.45703, 111.85691, 114.24108, 117.02742, 119.14901)

local x = 1
local y = 1
local years 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019

foreach i of local years{
	replace sale = sale*cpi[1,`x']/100 if fyear==`i'
	replace at = at*cpi[1,`x']/100 if fyear==`i'
	replace ib = ib*cpi[1,`x']/100 if fyear==`i'
	replace dp = dp*cpi[1,`x']/100 if fyear==`i'
	replace csho = csho*cpi[1,`x']/100 if fyear==`i'
	replace prcc_f= prcc_f*cpi[1,`x']/100 if fyear==`i'
	replace ceq = ceq*cpi[1,`x']/100 if fyear==`i'
	replace ppent = ppent*cpi[1,`x']/100 if fyear==`i'
	replace dlc = dlc*cpi[1,`x']/100 if fyear==`i'
	replace dltt = dltt*cpi[1,`x']/100 if fyear==`i'
	replace oibdp = oibdp*cpi[1,`x']/100 if fyear==`i'
	replace che = che*cpi[1,`x']/100 if fyear==`i'
	replace capx = capx*cpi[1,`x']/100 if fyear==`i'
	local x =`x'+1
}


generate sg = (sale)/(l.sale)-1
generate inv= (capx)/l.at
generate ta = at
generate fsize = log(at)
generate ppe = ppent/at
generate roa = ib/at
generate tq = (at+(csho*prcc_f)-ceq)/at
generate lev = (dltt+dlc)/at
generate cf = (oibdp)/l.at
generate liq = che/at
generate y = sale/at
replace sale = log(sale/1000)

keep gvkey fyear sale sg inv ta fsize tq lev cf ppe y roa liq
tab fyear

export excel using "[Your directory]\data1.xls", firstrow(variables)