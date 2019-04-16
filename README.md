Introduction
------------

This `R` package proposes different functions to make inference on a
random changepoint framework for longitudinal data. The work is still in
progress and the package is not yet fully functional. Already
implemented or being implemented are : \* testRCPMM: a test for the
existence of a random changepoint for longitudinal data \* rcpme: an
estimation algorithm for random changepoint mixed models \* bircpme: an
estimation algorithm for bivariate random changepoint mixed models
taking into account an eventual correlation between two markers

In the following, we will see how to use these functions and how to
manipulate them with a toy dataset.

The dataset: PAQUID cohort
--------------------------

The PAQUID cohort is an epidemiologic study on cogntive ageing launched
in South-West France in 1988. In the `rcpm` package, a sample of one
hundred randomly selected demented subjects has been extracted and
provided as a toy dataset in the package. Let us import the dataset into
the working environment and let us look at the data !

    library(rcpm)
    data(paquid)
    head(paquid)

    ##   id ist      delay el
    ## 1  1  25 -1.9700001  1
    ## 2  1  25 -1.8621301  1
    ## 3  1  28 -1.6521401  1
    ## 4  1  23 -1.1809501  1
    ## 5  1  16 -0.1474333  1
    ## 6  1  15  0.1474333  1

The dataset contains one line per visit. On each line, we can read `id`
the patient identifier , `ist` its current Isaacs Set Test valuewhich
measures memory impairment, `delay` the current delay to patient
diagnostic of dementia divided by ten and `el` a binary variable for
educational level (0 for non primary school certificate patients, 1 for
other).

It is interesting to look at the all longitudinal trajectory according
to educational level

![](Readme_files/figure-markdown_strict/unnamed-chunk-2-1.png)

Testing the existence of the random changepoint
-----------------------------------------------

A first step is to assess if there is a random changepoint in these
trajectories for both groups. This will be done with the `rcpm` function
`testRCPMM`. But first let us create two datasets for each educational
level.

    paquid_el0 <- paquid[paquid$el == 0, ]
    paquid_el1 <- paquid[paquid$el == 1, ]

To apply the test, we need to specify a formula indicating the score
variable, the time variable and the grouping variable considered. This
formula must be written this way `score ~ time | group` so that in our
example, it should look like `ist ~ delay | id`. This first two
arguments are essential for the function to run. Next argument is
`covariate` but the underlying functionnality has not been developed yet
so I will not give details about it. The `gamma` argument is used to
smooth the trajectory at the date of the changepoint, it has to be fixed
at a small value regarding the timescale, here `delay`. Be careful, by
default, its value is 0.1 so if your timescale is very reduced, let us
say from 0 to 0.1, you need to chose a smaller value for `gamma`. The
argument `nbnodes` indicates the number of nodes used for the
pseudo-adaptive gaussian quadrature which is used tu numerically compute
the integrals. I would advise not to touch this value. Finally, the
`nbpert` is important because it fixed the size of the empirical null
test statistics which will be used to compute the empirical pvalue.
Default is 500 which might induced quite a long computing time.

    test_el0 <- testRCPMM(longdata = paquid_el0, formu = ist ~ delay | id)
    test_el0

    ## $`empirical p-value`
    ## [1] 0.092
    ## 
    ## $`obs. stat`
    ## [1] -5.513509

    test_el1 <- testRCPMM(longdata = paquid_el1, formu = ist ~ delay | id)
    test_el1

    ## $`empirical p-value`
    ## [1] 0
    ## 
    ## $`obs. stat`
    ## [1] -23.24749

Estimating the random changepoint model
---------------------------------------

Now that we tested for the existence of a random changepoint and did not
rejected its absence. We can estimate the random changepoint model.

    esti <- rcpme(paquid, ist ~ delay | id, nbnodes = 10, model = "bw", link = "linear")
    esti

    ## $call
    ## $call[[1]]
    ## rcpme
    ## 
    ## $call$longdata
    ## paquid
    ## 
    ## $call$formu
    ## ist ~ delay | id
    ## 
    ## $call$nbnodes
    ## [1] 10
    ## 
    ## $call$model
    ## [1] "bw"
    ## 
    ## $call$link
    ## [1] "linear"
    ## 
    ## 
    ## $Loglik
    ## [1] -1896.379
    ## 
    ## $formula
    ## ist ~ delay | id
    ## 
    ## $fixed
    ##           par se(par)   ICinf  ICsup Wald stat. pvalue
    ## beta0  26.007   0.809  24.422 27.592     32.153      0
    ## beta1 -11.248   1.464 -14.118 -8.378     -7.681      0
    ## beta2  -9.152   1.475 -12.042 -6.261         NA     NA
    ## mutau  -0.182   0.079  -0.337 -0.027         NA     NA
    ## 
    ## $sdres
    ## [1] 3.138221
    ## 
    ## $VarEA
    ##           [,1]      [,2]      [,3]       [,4]
    ## [1,]  32.61864 -25.47839 -31.48569 0.00000000
    ## [2,] -25.47839  35.50965  42.82904 0.00000000
    ## [3,] -31.48569  42.82904  53.71886 0.00000000
    ## [4,]   0.00000   0.00000   0.00000 0.01468409
    ## 
    ## $optpar
    ##  [1]  26.0069299 -11.2480216  -9.1516442  -0.1821298   3.1382205
    ##  [6]   5.7112729  -4.4610698  -5.5129022   3.9507599   4.6157204
    ## [11]  -1.4219320  -0.1211779
    ## 
    ## $covariate
    ## [1] "NULL"
    ## 
    ## $REadjust
    ## [1] "no"
    ## 
    ## $invhessian
    ##               [,1]         [,2]         [,3]          [,4]         [,5]
    ##  [1,]  0.654257291  0.091924284 -0.180450652 -0.0312181127  0.003524816
    ##  [2,]  0.091924284  2.144498379  2.041761036 -0.0945021345  0.021690701
    ##  [3,] -0.180450652  2.041761036  2.174515993 -0.0787160980  0.020998281
    ##  [4,] -0.031218113 -0.094502134 -0.078716098  0.0062475238 -0.001113132
    ##  [5,]  0.003524816  0.021690701  0.020998281 -0.0011131317  0.011341151
    ##  [6,] -0.029866978 -0.026504131 -0.004039567  0.0019591489 -0.005801985
    ##  [7,]  0.293922532  1.002138630  0.850571568 -0.0604840426  0.016489529
    ##  [8,]  0.295706586  0.968659647  0.803175664 -0.0578148513  0.022345045
    ##  [9,] -0.197322332 -0.747552424 -0.653813342  0.0414511813 -0.011697451
    ## [10,] -0.195667632 -0.804812803 -0.699929822  0.0431733014 -0.010901295
    ## [11,] -0.030139575 -0.076785880 -0.053525128  0.0050868850  0.009890698
    ## [12,]  0.003881245  0.006592269  0.004588295 -0.0007719516  0.000252938
    ##                [,6]        [,7]         [,8]        [,9]        [,10]
    ##  [1,] -0.0298669776  0.29392253  0.295706586 -0.19732233 -0.195667632
    ##  [2,] -0.0265041305  1.00213863  0.968659647 -0.74755242 -0.804812803
    ##  [3,] -0.0040395665  0.85057157  0.803175664 -0.65381334 -0.699929822
    ##  [4,]  0.0019591489 -0.06048404 -0.057814851  0.04145118  0.043173301
    ##  [5,] -0.0058019851  0.01648953  0.022345045 -0.01169745 -0.010901295
    ##  [6,]  0.3307126946 -0.26140054 -0.431186527  0.01749367  0.011827161
    ##  [7,] -0.2614005359  1.30700037  1.393113575 -0.36937547 -0.502638089
    ##  [8,] -0.4311865265  1.39311358  1.694822899 -0.35430316 -0.487110816
    ##  [9,]  0.0174936723 -0.36937547 -0.354303159  0.77400882  0.656925998
    ## [10,]  0.0118271613 -0.50263809 -0.487110816  0.65692600  0.825453349
    ## [11,] -0.0227756993 -0.09639364 -0.017113718 -0.01912353  0.084162050
    ## [12,] -0.0001203908  0.01042521  0.009907703  0.01114598  0.005329226
    ##              [,11]         [,12]
    ##  [1,] -0.030139575  0.0038812452
    ##  [2,] -0.076785880  0.0065922691
    ##  [3,] -0.053525128  0.0045882948
    ##  [4,]  0.005086885 -0.0007719516
    ##  [5,]  0.009890698  0.0002529380
    ##  [6,] -0.022775699 -0.0001203908
    ##  [7,] -0.096393639  0.0104252055
    ##  [8,] -0.017113718  0.0099077027
    ##  [9,] -0.019123526  0.0111459752
    ## [10,]  0.084162050  0.0053292261
    ## [11,]  0.272768253 -0.0029124696
    ## [12,] -0.002912470  0.0016784404
    ## 
    ## $conv
    ## [1] 1
    ## 
    ## $init
    ##  [1] 20.5376137 -6.3719949 -0.5000000 -0.5504476  3.9294680  4.0595700
    ##  [7]  0.0000000  0.0000000  1.0000000  0.0000000  1.0000000  1.0000000
    ## 
    ## $model
    ## [1] "bw"
    ## 
    ## $gamma
    ## [1] 0.1
    ## 
    ## $link
    ## [1] "linear"

    pred <- unlist(IndPred(esti))
    datapred <- cbind(paquid, pred)
    datapred[1:25,]

    ##    id ist      delay el     pred
    ## 1   1  25 -1.9700001  1 25.97140
    ## 2   1  25 -1.8621301  1 25.67891
    ## 3   1  28 -1.6521401  1 25.10431
    ## 4   1  23 -1.1809501  1 23.77074
    ## 5   1  16 -0.1474333  1 18.72309
    ## 6   1  15  0.1474333  1 14.85978
    ## 7   1  NA  0.3637235  1 11.61085
    ## 8   2  24 -2.0594318  1 23.45999
    ## 9   2  21 -1.9567618  1 23.14135
    ## 10  2  25 -1.7566218  1 22.51987
    ## 11  2  21 -1.5406018  1 21.84842
    ## 12  2  18 -1.2722918  1 21.01288
    ## 13  2  20 -1.0423202  1 20.29435
    ## 14  2  17 -0.7674400  1 19.42936
    ## 15  2  21 -0.5683983  1 18.79335
    ## 16  2  NA -0.2684463  1 17.78659
    ## 17  2  17  0.2684463  1 15.67042
    ## 18  3  31 -1.1492129  1 32.59069
    ## 19  3  34 -1.0522929  1 32.39302
    ## 20  3  37 -0.8507829  1 31.86614
    ## 21  3  31 -0.6405229  1 30.98521
    ## 22  3  27 -0.3796029  1 28.58188
    ## 23  3  22 -0.1400411  1 23.56903
    ## 24  3  14  0.1400411  1 14.90966
    ## 25  4  25 -1.1649555  1 25.69294

    sampID <- sample(datapred$id, size = 12)
    datapred$el <- as.factor(datapred$el)
    datapred <- transform(datapred, EL = ifelse(el==0, "LEL", "HEL"))
    #ggplot(data = datapred[datapred$id %in% sampID,], aes(x = delay, y = ist)) +  geom_point() + geom_line(data=datapred[datapred$ID %in% sampID,], aes(x=delay, y=pred)) + facet_wrap(~ id + el, labeller = label_wrap_gen(multi_line=FALSE))
    ggplot(data = datapred[datapred$id %in% sampID,], aes(x = delay, y = ist, colour = el)) +  geom_point() + geom_line(data=datapred[datapred$id %in% sampID,], aes(x=delay, y=pred)) + facet_wrap(~ id + EL, labeller = label_wrap_gen(multi_line=FALSE))

    ## Warning: Removed 7 rows containing missing values (geom_point).

![](Readme_files/figure-markdown_strict/unnamed-chunk-6-1.png)
