program nlCesFe

     /*
    at[.,.] should be a row vector whose columns correspond to the
    value of beta, lambda, alphaF, alphaD...in that order.
    */
    version 18
    syntax varlist if, at(name)
    local lhs: word 1 of `varlist'
    local F: word 2 of `varlist'
    local D: word 3 of `varlist'
    local K: word 4 of `varlist'

    tempname beta lambda alphaF alphaD alpha delta
    scalar `beta' = `at'[1,1]
    scalar `lambda' = `at'[1,2]
    scalar `alphaF' = `at'[1,3]
    scalar `alphaD' = `at'[1,4]
    scalar `alpha' = `at'[1,5]
    scalar `delta' = `at'[1,6]

    * Compute the CES part
    replace `lhs' = `delta' + `alpha' * log(`K') + ///
    (1-`alpha')/`beta' * log( ///
    `lambda'^(1-`beta') * (`alphaF' * `F')^(`beta') + ///
    (1-`lambda')^(1 - `beta') * (`alphaD' * `D')^(`beta')) `if'

    * Add in the year fixed effects
    loc TfePos = $M + 1            // There are the 5 params (M=5) in the above CES part -> position 6 is where FEs start
    loc StartYear = $FirstYear + 1
    forvalues y = `StartYear'/$LastYear { // Calculate the time fixed effects 

        replace `lhs' = `lhs' + `at'[1,`TfePos']*(year == `y') `if'
        loc TfePos = `TfePos' + 1
    }
    
    * Add in the State fixed effects
    loc SfePos = ($M + 1) + ($LastYear - $FirstYear) 
    forvalues c = 2/$NRegions { // Calculate county fixed effects
        replace `lhs' = `lhs' + ///
        `at'[1,`SfePos'] * (state == `c') `if'
        local SfePos = `SfePos' + 1
    }
    
end