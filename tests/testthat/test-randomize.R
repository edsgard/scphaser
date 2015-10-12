
library('scphaser')
context('Randomization of acset data-structure')

test_that('racset swaps allele counts', {

    ##create alt and ref matrixes with no shared values
    nels = 100
    nrows = 5
    ncols = nels / nrows
    refcount = matrix(sample(1:nels, nels), nrow = nrows, dimnames = list(1:nrows, 1:ncols))
    altcount = matrix(sample((nels+1):(nels*2), nels), nrow = nrows, dimnames = list(1:nrows, 1:ncols))

    acset = new_acset(1:rownames(refcount), refcount, altcount)

    ##randomize
    rac = racset(acset)

    
    ##*###
    ##check that total number of elements hasn't changed
    ##*###
    
    ##expected
    origsize = length(refcount)

    ##observed
    refsize = length(rac[['refcount']])
    altsize = length(rac[['altcount']])

    ##test
    expect_identical(refsize, origsize)
    expect_identical(altsize, origsize)

    
    ##*##
    ##check that half of the elements are identical to the original
    ##*##

    ##expected
    nswapped = as.integer(floor(origsize / 2))
    nkept = origsize - nswapped

    ##observed
    nref_kept = length(which(rac[['refcount']] == refcount))
    nref_swapped = length(which(rac[['refcount']] == altcount))
    nalt_kept = length(which(rac[['altcount']] == altcount))
    nalt_swapped = length(which(rac[['altcount']] == refcount))

    ##test
    expect_identical(nref_swapped, nswapped)
    expect_identical(nalt_swapped, nswapped)
    expect_identical(nref_kept, nkept)
    expect_identical(nalt_kept, nkept)    
    
})
