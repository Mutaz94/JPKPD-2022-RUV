Wed Sep 29 21:49:05 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa1/A1/dat9.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m9.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1921.46145527466        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1233E+02 -1.7471E+01 -5.3706E+00 -1.7433E+00  8.5026E+01  2.1623E+01  1.5753E+00 -7.4031E-01  8.4617E+00 -1.3076E+01
            -2.4407E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1958.25243379279        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.3028E-01  9.4905E-01  8.8349E-01  1.0475E+00  8.6218E-01  9.8590E-01  9.6020E-01  9.7336E-01  9.1507E-01  9.3135E-01
             1.2881E+00
 PARAMETER:  2.7727E-02  4.7706E-02 -2.3881E-02  1.4645E-01 -4.8293E-02  8.5798E-02  5.9389E-02  7.3003E-02  1.1250E-02  2.8879E-02
             3.5316E-01
 GRADIENT:   1.9985E+02  3.8585E-01  1.6240E+00  1.2725E+01  8.1791E+00  1.9810E+01 -3.5965E+00  6.4192E+00 -6.3235E+00  7.8534E+00
             2.0057E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1962.50661982196        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      264
 NPARAMETR:  9.2757E-01  7.2381E-01  8.1093E-01  1.2270E+00  7.1795E-01  9.8537E-01  1.2526E+00  7.2811E-01  8.9142E-01  8.1673E-01
             1.2560E+00
 PARAMETER:  2.4814E-02 -2.2322E-01 -1.0958E-01  3.0461E-01 -2.3136E-01  8.5259E-02  3.2520E-01 -2.1730E-01 -1.4936E-02 -1.0245E-01
             3.2792E-01
 GRADIENT:  -1.9490E+01  2.8898E+01  9.7912E+00  5.5039E+01 -1.9147E+01  5.5400E-01 -1.4705E+00  1.2187E+00  1.0198E+01  1.4007E+00
             3.6573E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1964.67568950227        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  9.3714E-01  6.9021E-01  6.9431E-01  1.2007E+00  6.4889E-01  9.8231E-01  1.4506E+00  5.1422E-01  8.2428E-01  7.4604E-01
             1.2456E+00
 PARAMETER:  3.5072E-02 -2.7076E-01 -2.6484E-01  2.8294E-01 -3.3249E-01  8.2153E-02  4.7197E-01 -5.6511E-01 -9.3244E-02 -1.9298E-01
             3.1959E-01
 GRADIENT:   2.7685E+00  1.1683E+01  1.0021E+01  8.7470E+00 -1.4154E+01 -4.2147E-01  3.1607E-01 -9.3094E-01  2.5389E-03 -8.9521E-01
            -1.6329E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1965.34482275505        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  9.3080E-01  4.7802E-01  8.3000E-01  1.3344E+00  6.6074E-01  9.7607E-01  1.7116E+00  7.2538E-01  7.8594E-01  8.1208E-01
             1.2486E+00
 PARAMETER:  2.8286E-02 -6.3811E-01 -8.6329E-02  3.8851E-01 -3.1439E-01  7.5779E-02  6.3744E-01 -2.2106E-01 -1.4088E-01 -1.0816E-01
             3.2203E-01
 GRADIENT:  -2.8734E+00  2.8651E+00  6.9839E-01  4.5842E-01 -3.7418E+00 -1.1766E+00  6.5331E-02  1.2195E+00 -1.5037E+00  1.8692E+00
             3.0027E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1965.82992776384        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  9.2840E-01  2.6496E-01  9.3394E-01  1.4691E+00  6.6064E-01  9.7731E-01  2.1761E+00  8.2807E-01  7.6112E-01  8.3098E-01
             1.2384E+00
 PARAMETER:  2.5703E-02 -1.2282E+00  3.1656E-02  4.8463E-01 -3.1454E-01  7.7052E-02  8.7751E-01 -8.8661E-02 -1.7297E-01 -8.5148E-02
             3.1383E-01
 GRADIENT:   3.2450E+00  1.3389E+00  1.1115E+00  5.1824E+00  5.3482E-01  9.4946E-01 -4.7263E-01 -1.0431E+00  5.2667E-01 -2.0683E+00
            -4.6552E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1965.98585764403        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      972
 NPARAMETR:  9.2364E-01  1.4677E-01  9.8158E-01  1.5439E+00  6.5423E-01  9.7137E-01  2.8476E+00  9.1170E-01  7.3854E-01  8.5371E-01
             1.2462E+00
 PARAMETER:  2.0565E-02 -1.8189E+00  8.1409E-02  5.3429E-01 -3.2430E-01  7.0953E-02  1.1465E+00  7.5538E-03 -2.0308E-01 -5.8162E-02
             3.2011E-01
 GRADIENT:  -1.9718E+00  1.5004E+00  1.9054E+00  1.2630E+01 -6.1188E+00 -4.5124E-01 -3.0528E-01  4.8967E-01 -6.7603E-01  1.2248E+00
             1.3322E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1966.06968507281        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  9.2220E-01  8.5813E-02  1.0078E+00  1.5791E+00  6.5315E-01  9.7033E-01  3.7549E+00  9.4974E-01  7.2817E-01  8.5460E-01
             1.2471E+00
 PARAMETER:  1.9007E-02 -2.3556E+00  1.0777E-01  5.5683E-01 -3.2594E-01  6.9885E-02  1.4231E+00  4.8430E-02 -2.1723E-01 -5.7119E-02
             3.2084E-01
 GRADIENT:  -1.9712E+00  7.4465E-01  1.4528E+00  7.6492E+00 -4.3386E+00 -3.9130E-01 -8.5006E-02  3.2632E-01 -6.0544E-01  8.7872E-01
             1.9171E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1966.10527265935        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1323
 NPARAMETR:  9.2188E-01  4.9201E-02  1.0262E+00  1.6010E+00  6.5377E-01  9.7039E-01  5.0300E+00  9.7423E-01  7.2244E-01  8.5330E-01
             1.2449E+00
 PARAMETER:  1.8661E-02 -2.9118E+00  1.2591E-01  5.7065E-01 -3.2500E-01  6.9946E-02  1.7154E+00  7.3890E-02 -2.2512E-01 -5.8644E-02
             3.1907E-01
 GRADIENT:  -5.4799E-01  4.2540E-01  1.1296E+00  5.5294E+00 -2.6299E+00 -1.0994E-01  3.8235E-02  2.1932E-02 -4.1772E-01  2.4576E-01
             3.6490E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1966.12079779226        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1503
 NPARAMETR:  9.2199E-01  2.7715E-02  1.0324E+00  1.6054E+00  6.5457E-01  9.7047E-01  6.1207E+00  9.8615E-01  7.2078E-01  8.5288E-01
             1.2444E+00
 PARAMETER:  1.8780E-02 -3.4858E+00  1.3190E-01  5.7340E-01 -3.2377E-01  7.0029E-02  1.9117E+00  8.6054E-02 -2.2742E-01 -5.9137E-02
             3.1864E-01
 GRADIENT:   1.2256E+00 -8.0200E-02 -1.9477E+00 -1.6014E+01  5.1006E+00  1.2903E-01 -9.7315E-02  6.9442E-02  2.3295E-01 -3.0102E-01
             2.6591E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1966.13949784533        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1640
 NPARAMETR:  9.2190E-01  2.8150E-02  1.0340E+00  1.6087E+00  6.5368E-01  9.7038E-01  6.7285E+00  9.8599E-01  7.2052E-01  8.5405E-01
             1.2442E+00
 PARAMETER:  1.8677E-02 -3.4702E+00  1.3344E-01  5.7540E-01 -3.2514E-01  6.9928E-02  2.0064E+00  8.5886E-02 -2.2778E-01 -5.7770E-02
             3.1850E-01
 GRADIENT:   8.7770E-01  1.3811E-01  2.6326E-01 -7.4266E+00  5.9292E-01  6.6351E-02  9.5548E-02  3.1658E-03  3.0952E-01  5.1812E-02
             1.1083E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1966.14057409398        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1815
 NPARAMETR:  9.1975E-01  2.0620E-02  1.0327E+00  1.6120E+00  6.5146E-01  9.6925E-01  7.8109E+00  9.8715E-01  7.1839E-01  8.5236E-01
             1.2438E+00
 PARAMETER:  1.6342E-02 -3.7815E+00  1.3219E-01  5.7747E-01 -3.2854E-01  6.8768E-02  2.1555E+00  8.7070E-02 -2.3074E-01 -5.9744E-02
             3.1816E-01
 GRADIENT:  -4.0861E+00  9.3657E-02  2.9176E-01 -8.7749E+00  4.1326E-01 -3.2550E-01  5.2654E-02 -4.9165E-02 -2.6351E-01 -5.9854E-02
            -1.1754E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1966.14070609858        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1991
 NPARAMETR:  9.1886E-01  1.2817E-02  1.0296E+00  1.6154E+00  6.4831E-01  9.6874E-01  9.8192E+00  9.8752E-01  7.1703E-01  8.5047E-01
             1.2436E+00
 PARAMETER:  1.5383E-02 -4.2569E+00  1.2914E-01  5.7957E-01 -3.3339E-01  6.8244E-02  2.3843E+00  8.7441E-02 -2.3263E-01 -6.1961E-02
             3.1799E-01
 GRADIENT:  -5.8190E+00  5.6121E-02  3.2674E-01 -9.6826E+00 -1.2797E-01 -4.6488E-01  1.9937E-02 -4.8388E-02 -4.7176E-01 -7.3408E-02
            -1.7838E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1966.14754840412        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     2154
 NPARAMETR:  9.2266E-01  1.0000E-02  1.0276E+00  1.6157E+00  6.4713E-01  9.7028E-01  1.0919E+01  9.8830E-01  7.1791E-01  8.5034E-01
             1.2438E+00
 PARAMETER:  1.9508E-02 -4.5345E+00  1.2719E-01  5.7977E-01 -3.3520E-01  6.9825E-02  2.4905E+00  8.8232E-02 -2.3141E-01 -6.2119E-02
             3.1818E-01
 GRADIENT:   4.0596E+00  5.7300E-02 -2.2047E-01 -1.2203E+01  5.2800E-01  2.1695E-01 -1.3296E-03  9.5711E-02  1.2567E-01  5.4729E-02
             2.6972E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1966.15057742447        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2343
 NPARAMETR:  9.2086E-01  1.0340E-02  1.0276E+00  1.6178E+00  6.4705E-01  9.6984E-01  1.0721E+01  9.8691E-01  7.1782E-01  8.5008E-01
             1.2438E+00
 PARAMETER:  1.7556E-02 -4.4717E+00  1.2721E-01  5.8104E-01 -3.3533E-01  6.9378E-02  2.4722E+00  8.6829E-02 -2.3153E-01 -6.2420E-02
             3.1818E-01
 GRADIENT:  -5.6444E-01  2.4617E-01 -1.6437E-01 -6.4738E+00 -1.9283E-01  1.9308E-02 -3.7481E-03 -2.8981E-02  5.2817E-02 -4.9896E-03
            -4.9473E-02

0ITERATION NO.:   71    OBJECTIVE VALUE:  -1966.15057742447        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     2367
 NPARAMETR:  9.2086E-01  1.0340E-02  1.0276E+00  1.6178E+00  6.4705E-01  9.6984E-01  1.0721E+01  9.8691E-01  7.1782E-01  8.5008E-01
             1.2438E+00
 PARAMETER:  1.7556E-02 -4.4717E+00  1.2721E-01  5.8104E-01 -3.3533E-01  6.9378E-02  2.4722E+00  8.6829E-02 -2.3153E-01 -6.2420E-02
             3.1818E-01
 GRADIENT:  -7.0954E-01  4.4286E-03 -1.6200E-01  2.1588E+00 -9.5350E-02 -2.1835E-02 -9.5433E-03 -2.9634E-02 -1.5897E-02 -5.2578E-03
            -4.5082E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2367
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0654E-03  1.2481E-03 -1.8521E-02 -5.9215E-03 -2.2933E-02
 SE:             2.9777E-02  2.0230E-03  1.7845E-02  2.9133E-02  2.1063E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7146E-01  5.3728E-01  2.9932E-01  8.3893E-01  2.7623E-01

 ETASHRINKSD(%)  2.4416E-01  9.3223E+01  4.0217E+01  2.4017E+00  2.9438E+01
 ETASHRINKVR(%)  4.8773E-01  9.9541E+01  6.4260E+01  4.7458E+00  5.0210E+01
 EBVSHRINKSD(%)  5.5196E-01  9.3726E+01  4.1694E+01  2.8398E+00  2.8499E+01
 EBVSHRINKVR(%)  1.1009E+00  9.9606E+01  6.6004E+01  5.5990E+00  4.8876E+01
 RELATIVEINF(%)  9.4511E+01  2.1796E-02  5.2103E+00  6.6511E+00  6.6855E+00
 EPSSHRINKSD(%)  3.3147E+01
 EPSSHRINKVR(%)  5.5307E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1966.1505774244717     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1047.2120442197991     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.80
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1966.151       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.21E-01  1.03E-02  1.03E+00  1.62E+00  6.47E-01  9.70E-01  1.07E+01  9.87E-01  7.18E-01  8.50E-01  1.24E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.38E+03
 
 TH 2
+       -3.63E+01  1.14E+03
 
 TH 3
+       -6.94E+00  1.00E+02  3.71E+02
 
 TH 4
+       -8.85E+00  3.86E+02 -6.69E+01  7.78E+02
 
 TH 5
+        8.64E+00 -3.85E+02 -7.35E+02 -1.21E+02  1.88E+03
 
 TH 6
+        3.34E+00 -5.27E+00 -2.27E-01 -3.36E+00 -2.28E+00  2.06E+02
 
 TH 7
+        1.63E-03  1.97E+00  1.82E-03 -9.12E-03  6.18E-03  3.12E-04  4.28E-03
 
 TH 8
+        5.11E-02  7.82E+00 -3.54E+01 -3.42E+00  1.25E+00 -3.86E-02  2.22E-03  3.69E+01
 
 TH 9
+        3.98E+00 -2.73E+01  1.38E+01 -2.36E+00  4.20E+00 -9.09E-01  4.24E-02 -2.66E+00  3.44E+02
 
 TH10
+        1.38E+00  2.24E+01 -3.49E+00  4.75E+00 -9.13E+01  5.76E-01  5.58E-03  2.70E+01  3.87E+00  7.67E+01
 
 TH11
+       -1.15E+01 -2.62E+00 -1.15E+01 -1.28E+01 -7.67E+00  2.76E+00  6.42E-03  1.22E+01  1.21E+01  1.76E+01  2.57E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       42.940
Stop Time:
Wed Sep 29 21:49:50 CDT 2021
