program nlCesFe

     /*
    at[.,.] should be a row vector whose columns correspond to the
    value of beta, lambda, alphaF, alphaD...in that order.
    */
    version 18
    syntax varlist(min=3 max=3) if, at(name)
    local lhs: word 1 of `varlist'
    local F: word 2 of `varlist'
    local D: word 3 of `varlist'

    tempname beta lambda alphaF alphaD delta
    scalar `beta' = `at'[1,1]
    scalar `lambda' = `at'[1,2]
    scalar `alphaF' = `at'[1,3]
    scalar `alphaD' = `at'[1,4]
    scalar `delta' = `at'[1,5]

    * Compute the CES part
    replace `lhs' = `delta' + (1/`beta') * log( ///
    `lambda'^(1-`beta') * (`alphaF'*`F')^(`beta') + ///
    (1-`lambda')^(1 - `beta') * (`alphaD' * `D')^(`beta')) `if'

    * Add in the year fixed effects
    loc TfePos = 6            // There are the four CES parameters -> position 5 is where TFEs start
    forvalues y = 2006/2023 { // Calculate the time fixed effects 

        replace `lhs' = `lhs' + `at'[1,`TfePos']*(year == `y') `if'
        loc TfePos = `TfePos' + 1
    }
    
    * Add in the County fixed effects
    loc CfePos = 6 + 18   // We have 19 years in the data (2005-2023), but leave one out to avoid collinearity
    forvalues c = 2/508 { // Calculate county fixed effects
        replace `lhs' = `lhs' + ///
        `at'[1,`CfePos'] * (fipcodenum == `c') `if'
        local CfePos = `CfePos' + 1
    }
    
end