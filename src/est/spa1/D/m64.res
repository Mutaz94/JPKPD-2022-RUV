Thu Sep 30 03:23:37 CDT 2021
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
$DATA ../../../../data/spa1/D/dat64.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14295.4700794663        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8174E+02  4.0382E+02 -5.1635E+01  1.6703E+02  5.6602E+02 -2.8504E+03 -9.1416E+02 -1.1801E+02 -1.7652E+03 -1.1021E+03
            -2.5830E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -665.610515142928        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1276E+00  8.8544E-01  8.4835E-01  2.6445E+00  1.1682E+00  2.7694E+00  1.4096E+00  1.0370E+00  2.8132E+00  1.5272E+00
             1.2098E+01
 PARAMETER:  2.2009E-01 -2.1672E-02 -6.4468E-02  1.0725E+00  2.5546E-01  1.1186E+00  4.4328E-01  1.3633E-01  1.1343E+00  5.2343E-01
             2.5930E+00
 GRADIENT:  -6.6206E+01  4.3521E+01 -3.9349E+01  1.1877E+02 -1.0847E+01  2.3802E+00  1.8083E+00  7.9187E+00 -3.2584E+01  1.0734E+01
             1.7443E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -714.059875807748        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0127E+00  9.9341E-01  3.2102E+00  3.1380E+00  2.2377E+00  3.3776E+00  5.4507E+00  8.6629E-01  4.1636E+00  3.8091E+00
             9.2281E+00
 PARAMETER:  1.1263E-01  9.3384E-02  1.2663E+00  1.2436E+00  9.0546E-01  1.3172E+00  1.7957E+00 -4.3534E-02  1.5264E+00  1.4374E+00
             2.3222E+00
 GRADIENT:  -4.7258E+01  2.0095E+01  4.1752E+00  8.9249E+01 -4.6240E+01  9.2730E+01  1.2640E+01 -2.1437E-01  5.3892E+01  2.9911E+01
             9.8867E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -721.825837682687        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      324
 NPARAMETR:  1.0299E+00  8.9074E-01  3.2181E+00  3.1879E+00  2.1998E+00  3.2980E+00  6.1623E+00  8.6739E-01  3.8924E+00  3.9260E+00
             9.1875E+00
 PARAMETER:  1.2945E-01 -1.5699E-02  1.2688E+00  1.2594E+00  8.8838E-01  1.2933E+00  1.9184E+00 -4.2270E-02  1.4590E+00  1.4676E+00
             2.3178E+00
 GRADIENT:  -5.2803E+01  1.9023E+01  4.6355E+00  6.1223E+01 -5.0379E+01  4.1928E+00  6.2822E+00 -2.6703E-01 -2.6069E+00  3.3028E+01
             7.4322E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -773.615022304130        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      425
 NPARAMETR:  1.0296E+00  1.1156E-01  3.2240E+00  3.1378E+00  1.0442E+01  3.1766E+00  6.1243E+00  3.6905E+00  3.8569E+00  3.4604E+00
             8.6814E+00
 PARAMETER:  1.2917E-01 -2.0932E+00  1.2706E+00  1.2435E+00  2.4458E+00  1.2558E+00  1.9123E+00  1.4058E+00  1.4499E+00  1.3414E+00
             2.2612E+00
 GRADIENT:  -3.5625E+01 -7.5682E-03  1.7729E+01  1.0477E+02  1.5993E-01  1.0013E+02  1.9123E+00 -6.9935E+00  1.2495E+02  8.2400E-03
             1.0539E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -780.193661437760        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0290E+00  1.0000E-02  3.2265E+00  2.9514E+00  6.7738E+01  2.8610E+00  6.0481E+00  3.8240E+00  3.6937E+00  3.6556E+00
             7.3087E+00
 PARAMETER:  1.2863E-01 -4.7091E+00  1.2714E+00  1.1823E+00  4.3157E+00  1.1512E+00  1.8998E+00  1.4413E+00  1.4066E+00  1.3962E+00
             2.0891E+00
 GRADIENT:  -2.6506E+01  0.0000E+00  2.2324E+01  1.2163E+02  3.3314E-02  7.4385E+01  3.6894E-02 -4.3211E+00  1.3832E+02  1.3375E-04
            -4.8461E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -797.850486530666        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  1.0292E+00  1.0000E-02  3.2114E+00  2.8110E+00  1.3294E+02  2.7731E+00  6.0446E+00  7.6923E+00  3.5582E+00  3.7157E+00
             7.2664E+00
 PARAMETER:  1.2874E-01 -6.1146E+00  1.2667E+00  1.1335E+00  4.9899E+00  1.1200E+00  1.8992E+00  2.1402E+00  1.3693E+00  1.4126E+00
             2.0833E+00
 GRADIENT:  -2.9804E+01  0.0000E+00  2.4834E+00  6.9070E+01  1.9030E-02  1.9533E+00  2.7425E-02  1.8723E+01  7.6964E+01  2.5563E-06
            -7.9733E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -855.761556062240        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      834
 NPARAMETR:  1.0316E+00  1.0000E-02  2.9993E+00  1.2435E+00  1.7321E+06  1.8465E+00  6.0188E+00  9.4188E+00  1.8838E+00  2.1552E+00
             8.4714E+00
 PARAMETER:  1.3111E-01 -1.7126E+01  1.1984E+00  3.1795E-01  1.4465E+01  7.1331E-01  1.8949E+00  2.3427E+00  7.3328E-01  8.6787E-01
             2.2367E+00
 GRADIENT:   1.4027E+01  0.0000E+00  1.8128E+01 -5.6975E+01 -2.9140E-06 -4.2966E+01  3.2195E-02  5.3886E+01  5.8481E+01  0.0000E+00
             7.0867E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -884.630921008720        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1019             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0421E+00  1.0000E-02  5.3964E-01  1.5107E+00  7.1591E+05  2.2430E+00  3.3171E+00  8.6350E+00  1.0970E+00  2.3414E+00
             7.9012E+00
 PARAMETER:  1.4124E-01 -1.5755E+01 -5.1685E-01  5.1256E-01  1.3581E+01  9.0781E-01  1.2991E+00  2.2558E+00  1.9258E-01  9.5075E-01
             2.1670E+00
 GRADIENT:   2.8464E+01  0.0000E+00 -1.0978E+01  1.4431E+02  1.0600E-05  4.1008E+01  1.0300E-03  1.4175E+02 -6.2451E+00  0.0000E+00
            -1.9350E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -951.377568507858        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1115
 NPARAMETR:  7.2937E-01  1.0000E-02  6.4718E-02  5.6818E-01  7.5535E+05  1.9468E+00  7.1201E-01  2.8391E+00  2.3786E-01  2.4450E+00
             8.0525E+00
 PARAMETER: -2.1558E-01 -1.5755E+01 -2.6377E+00 -4.6532E-01  1.3635E+01  7.6620E-01 -2.3967E-01  1.1435E+00 -1.3361E+00  9.9406E-01
             2.1860E+00
 GRADIENT:   6.1561E+01  0.0000E+00 -3.0815E+01 -1.0072E+01  1.9456E-06 -3.0503E+01  3.9268E-04  6.2020E+01  1.9669E+00  0.0000E+00
             2.8281E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -970.878784951145        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1291
 NPARAMETR:  5.9849E-01  1.0000E-02  5.4162E-02  4.7670E-01  6.3333E+05  2.1003E+00  4.6649E-01  1.8557E+00  1.4640E-01  2.4759E+00
             7.7673E+00
 PARAMETER: -4.1335E-01 -1.5755E+01 -2.8158E+00 -6.4087E-01  1.3459E+01  8.4210E-01 -6.6251E-01  7.1828E-01 -1.8214E+00  1.0066E+00
             2.1499E+00
 GRADIENT:  -6.1742E-01  0.0000E+00  3.7743E+00 -5.2083E-01  2.6411E-06 -1.2700E+00  1.6651E-04 -1.4573E+00  3.5990E-01  0.0000E+00
             4.0865E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -971.706613227433        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1467            RESET HESSIAN, TYPE II
 NPARAMETR:  5.6941E-01  1.0000E-02  4.5501E-02  4.2828E-01  3.2450E+05  2.1126E+00  3.9940E-01  1.7719E+00  4.7950E-02  2.4828E+00
             7.7080E+00
 PARAMETER: -4.6316E-01 -1.5755E+01 -2.9900E+00 -7.4797E-01  1.2790E+01  8.4791E-01 -8.1779E-01  6.7204E-01 -2.9376E+00  1.0094E+00
             2.1423E+00
 GRADIENT:   5.6825E+01  0.0000E+00  8.0083E+01  3.6790E+01  4.9381E-06  4.2395E+01  1.7078E-04  5.1305E+00  6.8322E-02  0.0000E+00
             2.0424E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -971.713140750666        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1644
 NPARAMETR:  5.6771E-01  1.0000E-02  4.5707E-02  4.2762E-01  3.2450E+05  2.1052E+00  3.9930E-01  1.7764E+00  4.2347E-02  2.4828E+00
             7.7106E+00
 PARAMETER: -4.6615E-01 -1.5755E+01 -2.9855E+00 -7.4952E-01  1.2790E+01  8.4443E-01 -8.1804E-01  6.7460E-01 -3.0618E+00  1.0094E+00
             2.1426E+00
 GRADIENT:   1.9456E-01  0.0000E+00 -1.6195E-01 -3.8957E-01  4.5033E-06  1.8942E-01  1.2955E-04 -1.0606E-01  3.0108E-02  0.0000E+00
             1.0983E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -971.731263554398        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1827
 NPARAMETR:  5.6828E-01  1.0000E-02  4.5539E-02  4.2783E-01  1.2808E+05  2.1098E+00  3.5029E-01  1.7788E+00  1.3421E-02  2.4825E+00
             7.7083E+00
 PARAMETER: -4.6514E-01 -1.5755E+01 -2.9892E+00 -7.4904E-01  1.1860E+01  8.4661E-01 -9.4901E-01  6.7593E-01 -4.2109E+00  1.0093E+00
             2.1423E+00
 GRADIENT:   7.5157E-01  0.0000E+00 -2.6268E+00  2.8630E+00  1.2314E-05  8.9724E-01  9.9827E-05  1.5137E-01  3.0783E-03  1.1264E-11
            -4.9602E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -971.798098517291        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     2022             RESET HESSIAN, TYPE I
 NPARAMETR:  5.6612E-01  1.0000E-02  4.5099E-02  4.2385E-01  3.7579E+01  2.1095E+00  1.3533E-01  1.7715E+00  1.0000E-02  2.6096E+00
             7.7091E+00
 PARAMETER: -4.6895E-01 -1.5755E+01 -2.9989E+00 -7.5838E-01  3.7265E+00  8.4645E-01 -1.9001E+00  6.7184E-01 -6.2872E+00  1.0592E+00
             2.1424E+00
 GRADIENT:   5.6520E+01  0.0000E+00  8.4933E+01  3.0942E+01  3.7924E-02  4.2050E+01  2.6320E-05  6.3060E+00  0.0000E+00  2.2079E-04
             2.1852E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -971.879614431353        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     2100
 NPARAMETR:  5.6428E-01  1.0000E-02  4.4897E-02  4.2350E-01  1.0844E+01  2.1028E+00  1.3504E-01  1.7700E+00  1.0000E-02  2.4494E+00
             7.7046E+00
 PARAMETER: -4.7221E-01 -1.5755E+01 -3.0034E+00 -7.5920E-01  2.4836E+00  8.4326E-01 -1.9022E+00  6.7100E-01 -6.2872E+00  9.9585E-01
             2.1418E+00
 GRADIENT:   5.4163E+01  0.0000E+00  8.5487E+01  3.2201E+01 -2.1362E-02  4.0596E+01  2.6779E-05  7.0020E+00  0.0000E+00  7.0126E-03
             2.1793E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -971.893256740572        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2286             RESET HESSIAN, TYPE I
 NPARAMETR:  5.6853E-01  1.0000E-02  4.4824E-02  4.2454E-01  1.1053E+01  2.1115E+00  1.2380E-01  1.7632E+00  1.0000E-02  1.1066E+00
             7.7031E+00
 PARAMETER: -4.6471E-01 -1.5755E+01 -3.0050E+00 -7.5675E-01  2.5027E+00  8.4738E-01 -1.9891E+00  6.6711E-01 -6.2872E+00  2.0129E-01
             2.1416E+00
 GRADIENT:   5.6489E+01  0.0000E+00  8.2590E+01  3.5383E+01  6.0687E-02  4.2065E+01  2.2555E-05  5.9219E+00  0.0000E+00  9.0955E-04
             2.0494E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -971.914294908100        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     2432
 NPARAMETR:  5.6461E-01  1.0000E-02  4.4189E-02  4.2104E-01  1.0623E+01  2.1056E+00  1.0000E-02  1.7582E+00  1.0000E-02  9.4070E-01
             7.7052E+00
 PARAMETER: -4.7163E-01 -1.5755E+01 -3.0193E+00 -7.6502E-01  2.4630E+00  8.4460E-01 -5.2846E+00  6.6427E-01 -6.2872E+00  3.8870E-02
             2.1419E+00
 GRADIENT:   5.5894E+01  0.0000E+00  8.2416E+01  3.7413E+01 -5.0879E-03  4.1012E+01  0.0000E+00  6.1643E+00  0.0000E+00  1.5362E-03
             2.1066E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -971.929010018262        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:     2582
 NPARAMETR:  5.6183E-01  1.0000E-02  4.3917E-02  4.1667E-01  1.0930E+01  2.0985E+00  1.0000E-02  1.7516E+00  1.0000E-02  9.1501E-01
             7.7008E+00
 PARAMETER: -4.7655E-01 -1.5755E+01 -3.0255E+00 -7.7545E-01  2.4915E+00  8.4121E-01 -5.2846E+00  6.6053E-01 -6.2872E+00  1.1181E-02
             2.1413E+00
 GRADIENT:  -6.0244E-01  0.0000E+00  6.2169E-01 -2.7526E+00 -1.5537E-02 -1.0963E+00  0.0000E+00  4.2627E-01  0.0000E+00  8.0492E-04
            -6.2983E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -971.943300490263        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2779
 NPARAMETR:  5.6367E-01  1.0000E-02  4.3813E-02  4.1709E-01  1.0875E+01  2.1112E+00  1.0000E-02  1.7493E+00  1.0000E-02  8.9766E-01
             7.7046E+00
 PARAMETER: -4.7329E-01 -1.5755E+01 -3.0278E+00 -7.7446E-01  2.4865E+00  8.4726E-01 -5.2846E+00  6.5920E-01 -6.2872E+00 -7.9674E-03
             2.1418E+00
 GRADIENT:   8.5211E-01  0.0000E+00 -1.4206E+00 -3.9898E-01 -6.5811E-03  1.0329E+00  0.0000E+00  2.2727E-01  0.0000E+00  8.5875E-04
             4.3910E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -971.949712130399        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     2978
 NPARAMETR:  5.6303E-01  1.0000E-02  4.3649E-02  4.1586E-01  1.0869E+01  2.1115E+00  1.0000E-02  1.7473E+00  1.0000E-02  8.7208E-01
             7.7047E+00
 PARAMETER: -4.7443E-01 -1.5755E+01 -3.0316E+00 -7.7741E-01  2.4859E+00  8.4741E-01 -5.2846E+00  6.5806E-01 -6.2872E+00 -3.6870E-02
             2.1418E+00
 GRADIENT:   9.7658E-01  0.0000E+00 -1.1114E+00 -1.0539E+00 -9.6610E-03  1.1253E+00  0.0000E+00  2.7516E-01  0.0000E+00  8.3292E-04
             1.7926E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -971.959670823585        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3154
 NPARAMETR:  5.5914E-01  1.0000E-02  4.3017E-02  4.1342E-01  1.0860E+01  2.1024E+00  1.0000E-02  1.7434E+00  1.0000E-02  8.4808E-01
             7.6965E+00
 PARAMETER: -4.8135E-01 -1.5755E+01 -3.0462E+00 -7.8330E-01  2.4851E+00  8.4310E-01 -5.2846E+00  6.5582E-01 -6.2872E+00 -6.4776E-02
             2.1408E+00
 GRADIENT:   5.7051E+01  0.0000E+00  8.4010E+01  3.8080E+01  2.3920E-02  4.0656E+01  0.0000E+00  5.9756E+00  0.0000E+00  8.3525E-04
             2.0430E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -971.966813866025        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     3298
 NPARAMETR:  5.6166E-01  1.0000E-02  4.3130E-02  4.1297E-01  1.0970E+01  2.1158E+00  1.0000E-02  1.7414E+00  1.0000E-02  8.3338E-01
             7.7038E+00
 PARAMETER: -4.7686E-01 -1.5755E+01 -3.0435E+00 -7.8437E-01  2.4951E+00  8.4945E-01 -5.2846E+00  6.5469E-01 -6.2872E+00 -8.2271E-02
             2.1417E+00
 GRADIENT:   5.8774E+01  0.0000E+00  8.6176E+01  3.3542E+01  3.9309E-02  4.3307E+01  0.0000E+00  5.8930E+00  0.0000E+00  6.6427E-04
             2.1514E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -971.972200647837        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3495             RESET HESSIAN, TYPE I
 NPARAMETR:  5.6036E-01  1.0000E-02  4.3042E-02  4.1182E-01  1.0908E+01  2.1110E+00  1.0000E-02  1.7397E+00  1.0000E-02  8.2858E-01
             7.7029E+00
 PARAMETER: -4.7917E-01 -1.5755E+01 -3.0456E+00 -7.8716E-01  2.4895E+00  8.4716E-01 -5.2846E+00  6.5372E-01 -6.2872E+00 -8.8043E-02
             2.1416E+00
 GRADIENT:   5.8561E+01  0.0000E+00  8.7471E+01  3.2196E+01  2.3113E-02  4.2445E+01  0.0000E+00  5.8833E+00  0.0000E+00  7.4486E-04
             2.1631E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -971.977976002003        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     3691             RESET HESSIAN, TYPE I
 NPARAMETR:  5.5996E-01  1.0000E-02  4.2873E-02  4.1102E-01  1.0912E+01  2.1122E+00  1.0000E-02  1.7383E+00  1.0000E-02  7.6018E-01
             7.7021E+00
 PARAMETER: -4.7989E-01 -1.5755E+01 -3.0495E+00 -7.8911E-01  2.4899E+00  8.4771E-01 -5.2846E+00  6.5292E-01 -6.2872E+00 -1.7420E-01
             2.1415E+00
 GRADIENT:   5.8828E+01  0.0000E+00  8.7099E+01  3.3023E+01  2.6479E-02  4.2681E+01  0.0000E+00  5.9195E+00  0.0000E+00  6.4276E-04
             2.1505E+01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -971.982260064690        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3888             RESET HESSIAN, TYPE I
 NPARAMETR:  5.5922E-01  1.0000E-02  4.2748E-02  4.0999E-01  1.0904E+01  2.1114E+00  1.0000E-02  1.7363E+00  1.0000E-02  7.2913E-01
             7.7020E+00
 PARAMETER: -4.8121E-01 -1.5755E+01 -3.0524E+00 -7.9162E-01  2.4892E+00  8.4736E-01 -5.2846E+00  6.5177E-01 -6.2872E+00 -2.1591E-01
             2.1415E+00
 GRADIENT:   5.8941E+01  0.0000E+00  8.7720E+01  3.2473E+01  2.2588E-02  4.2584E+01  0.0000E+00  5.8856E+00  0.0000E+00  6.1321E-04
             2.1606E+01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -971.986873656558        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     4085             RESET HESSIAN, TYPE I
 NPARAMETR:  5.5873E-01  1.0000E-02  4.2607E-02  4.0917E-01  1.0910E+01  2.1118E+00  1.0000E-02  1.7349E+00  1.0000E-02  6.7011E-01
             7.7014E+00
 PARAMETER: -4.8209E-01 -1.5755E+01 -3.0557E+00 -7.9364E-01  2.4897E+00  8.4756E-01 -5.2846E+00  6.5092E-01 -6.2872E+00 -3.0032E-01
             2.1414E+00
 GRADIENT:   5.9142E+01  0.0000E+00  8.7718E+01  3.2765E+01  2.4081E-02  4.2686E+01  0.0000E+00  5.8976E+00  0.0000E+00  5.2926E-04
             2.1566E+01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -971.990522591442        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     4281             RESET HESSIAN, TYPE I
 NPARAMETR:  5.5804E-01  1.0000E-02  4.2511E-02  4.0872E-01  1.0918E+01  2.1104E+00  1.0000E-02  1.7329E+00  1.0000E-02  6.3148E-01
             7.7005E+00
 PARAMETER: -4.8332E-01 -1.5755E+01 -3.0580E+00 -7.9472E-01  2.4905E+00  8.4686E-01 -5.2846E+00  6.4982E-01 -6.2872E+00 -3.5969E-01
             2.1413E+00
 GRADIENT:   5.8984E+01  0.0000E+00  8.7491E+01  3.3517E+01  2.5516E-02  4.2406E+01  0.0000E+00  5.7350E+00  0.0000E+00  4.7135E-04
             2.1422E+01

0ITERATION NO.:  140    OBJECTIVE VALUE:  -971.994044689916        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     4477             RESET HESSIAN, TYPE I
 NPARAMETR:  5.5758E-01  1.0000E-02  4.2393E-02  4.0776E-01  1.0899E+01  2.1108E+00  1.0000E-02  1.7316E+00  1.0000E-02  5.2711E-01
             7.7006E+00
 PARAMETER: -4.8414E-01 -1.5755E+01 -3.0608E+00 -7.9709E-01  2.4887E+00  8.4708E-01 -5.2846E+00  6.4904E-01 -6.2872E+00 -5.4034E-01
             2.1413E+00
 GRADIENT:   5.9248E+01  0.0000E+00  8.8053E+01  3.2907E+01  2.1695E-02  4.2535E+01  0.0000E+00  5.7890E+00  0.0000E+00  3.6005E-04
             2.1549E+01

0ITERATION NO.:  145    OBJECTIVE VALUE:  -971.998071305119        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     4673             RESET HESSIAN, TYPE I
 NPARAMETR:  5.5721E-01  1.0000E-02  4.2257E-02  4.0706E-01  1.0917E+01  2.1115E+00  1.0000E-02  1.7305E+00  1.0000E-02  4.5497E-01
             7.6999E+00
 PARAMETER: -4.8481E-01 -1.5755E+01 -3.0640E+00 -7.9879E-01  2.4903E+00  8.4740E-01 -5.2846E+00  6.4840E-01 -6.2872E+00 -6.8753E-01
             2.1412E+00
 GRADIENT:   5.9464E+01  0.0000E+00  8.7824E+01  3.3485E+01  2.5702E-02  4.2679E+01  0.0000E+00  5.8228E+00  0.0000E+00  2.7173E-04
             2.1456E+01

0ITERATION NO.:  150    OBJECTIVE VALUE:  -972.000785463458        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     4870
 NPARAMETR:  5.5656E-01  1.0000E-02  4.2179E-02  4.0652E-01  1.0912E+01  2.1101E+00  1.0000E-02  1.7287E+00  1.0000E-02  4.2999E-01
             7.6995E+00
 PARAMETER: -4.8598E-01 -1.5755E+01 -3.0658E+00 -8.0012E-01  2.4899E+00  8.4676E-01 -5.2846E+00  6.4736E-01 -6.2872E+00 -7.4399E-01
             2.1412E+00
 GRADIENT:   8.2406E-01  0.0000E+00 -1.8494E+00 -7.3476E-01 -5.2766E-03  1.0294E+00  0.0000E+00  1.8276E-01  0.0000E+00  2.0809E-04
             3.1234E-02

0ITERATION NO.:  155    OBJECTIVE VALUE:  -972.003082953788        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     5053
 NPARAMETR:  5.5623E-01  1.0000E-02  4.2096E-02  4.0624E-01  1.0928E+01  2.1099E+00  1.0000E-02  1.7277E+00  1.0000E-02  1.0000E-02
             7.6987E+00
 PARAMETER: -4.8657E-01 -1.5755E+01 -3.0678E+00 -8.0082E-01  2.4913E+00  8.4666E-01 -5.2846E+00  6.4679E-01 -6.2872E+00 -5.6767E+00
             2.1411E+00
 GRADIENT:   7.3931E-01  0.0000E+00 -2.4082E+00  1.1424E-01  1.0763E-03  9.7448E-01  0.0000E+00  1.3447E-01  0.0000E+00  0.0000E+00
            -1.2412E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     5053
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0678E-03 -1.8547E-06  6.1465E-03 -4.1465E-04 -9.8082E-07
 SE:             2.9264E-02  7.7160E-07  2.6183E-02  2.3281E-04  1.1959E-06
 N:                     100         100         100         100         100

 P VAL.:         8.3574E-01  1.6227E-02  8.1440E-01  7.4899E-02  4.1211E-01

 ETASHRINKSD(%)  1.9630E+00  9.9997E+01  1.2284E+01  9.9220E+01  9.9996E+01
 ETASHRINKVR(%)  3.8874E+00  1.0000E+02  2.3058E+01  9.9994E+01  1.0000E+02
 EBVSHRINKSD(%)  1.7773E+00  9.9997E+01  1.1326E+01  9.9208E+01  9.9996E+01
 EBVSHRINKVR(%)  3.5230E+00  1.0000E+02  2.1370E+01  9.9994E+01  1.0000E+02
 RELATIVEINF(%)  2.0590E+01  1.4277E-08  1.6197E+00  8.5792E-05  4.9237E-08
 EPSSHRINKSD(%)  1.2839E+01
 EPSSHRINKVR(%)  2.4029E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -972.00308295378784     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -53.064549749115145     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   105.37
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     9.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -972.003       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.56E-01  1.00E-02  4.21E-02  4.06E-01  1.09E+01  2.11E+00  1.00E-02  1.73E+00  1.00E-02  1.00E-02  7.70E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.21E+00
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -5.90E+02  0.00E+00  1.58E+05
 
 TH 4
+        8.67E+01  0.00E+00 -2.32E+04  3.40E+03
 
 TH 5
+        9.34E-03  0.00E+00 -2.49E+00  3.66E-01  3.95E-05
 
 TH 6
+       -3.40E-01  0.00E+00  9.07E+01 -1.33E+01 -1.44E-03  5.22E-02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        4.79E-01  0.00E+00 -1.28E+02  1.88E+01  2.02E-03 -7.37E-02  0.00E+00  1.04E-01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.15E-01  0.00E+00  1.11E+02 -1.63E+01 -1.75E-03  6.37E-02  0.00E+00 -8.99E-02  0.00E+00  0.00E+00  7.78E-02
 
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
+        7.71E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.64E+02  0.00E+00  1.61E+05
 
 TH 4
+       -4.18E+02  0.00E+00 -2.36E+04  4.04E+03
 
 TH 5
+        1.59E-01  0.00E+00 -2.55E+00  3.74E-01  6.75E-03
 
 TH 6
+        4.78E+00  0.00E+00  9.15E+01 -2.49E+01  5.76E-03  3.99E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        3.77E+00  0.00E+00 -1.41E+02 -5.63E+01 -1.99E-02  1.60E+00  0.00E+00  3.90E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -9.73E+00  0.00E+00  1.13E+02 -1.90E+01 -3.90E-03  1.25E+00  0.00E+00  2.19E+00  0.00E+00  0.00E+00  8.57E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        7.79E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.29E+03  0.00E+00  1.65E+05
 
 TH 4
+       -2.77E+02  0.00E+00 -2.29E+04  3.84E+03
 
 TH 5
+        1.64E-01  0.00E+00 -2.93E+00  3.64E-01  1.02E-04
 
 TH 6
+        1.07E+02  0.00E+00 -6.78E+02 -1.58E+00  3.32E-02  6.57E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.78E+01  0.00E+00  9.22E+01 -7.65E+01 -1.76E-02  9.84E-01  0.00E+00  2.60E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.05E+02  0.00E+00  1.14E+03 -1.50E+02 -3.89E-02  8.88E+00  0.00E+00 -2.80E-01  0.00E+00  0.00E+00  1.95E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,      115.004
Stop Time:
Thu Sep 30 03:25:34 CDT 2021
