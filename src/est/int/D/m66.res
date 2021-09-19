Sat Sep 18 07:20:36 CDT 2021
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
$DATA ../../../../data/int/D/dat66.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   33862.2953514065        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.4607E+02  4.1249E+02 -1.0770E+01  2.8074E+02  1.8932E+02 -2.9923E+03 -1.3609E+03 -9.2402E+01 -2.0378E+03 -1.1043E+03
            -6.7636E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -868.732892163377        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  2.3419E+00  2.1723E+00  9.4007E-01  1.3522E+00  9.0691E-01  5.0658E+00  7.9237E+00  9.9553E-01  2.6962E+00  1.8504E+00
             1.1547E+01
 PARAMETER:  9.5096E-01  8.7580E-01  3.8202E-02  4.0175E-01  2.2905E-03  1.7225E+00  2.1699E+00  9.5521E-02  1.0918E+00  7.1540E-01
             2.5465E+00
 GRADIENT:   4.8981E+01  1.7170E+01 -3.5456E+01  2.0308E+01 -5.3517E+01  1.0129E+02  8.9688E+01  4.6694E+00  8.0017E+01  3.8687E+01
             3.0943E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -937.175804317006        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.2998E+00  4.0094E+00  1.8452E+01  2.5414E+00  1.9203E+00  4.0369E+00  1.5892E+01  9.6635E-01  2.9222E+00  1.9185E+00
             1.1530E+01
 PARAMETER:  3.6224E-01  1.4886E+00  3.0152E+00  1.0327E+00  7.5248E-01  1.4955E+00  2.8658E+00  6.5769E-02  1.1723E+00  7.5155E-01
             2.5450E+00
 GRADIENT:   1.0637E+01  2.4938E+01 -1.0448E+00  2.1264E+01 -3.7655E+01  1.1214E+02  3.9124E+01  8.6763E-02  4.0317E+01  4.6454E+01
             3.1467E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1052.46028424147        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  1.7325E+00  7.0519E-01  1.8671E+01  1.6246E+00  2.6695E+00  3.0229E+00  7.6957E+00  1.9231E+00  1.6221E+00  1.0528E+00
             1.1043E+01
 PARAMETER:  6.4954E-01 -2.4929E-01  3.0270E+00  5.8528E-01  1.0819E+00  1.2062E+00  2.1407E+00  7.5395E-01  5.8374E-01  1.5144E-01
             2.5018E+00
 GRADIENT:   8.4241E+01 -4.5506E+00 -5.3495E+00  2.1188E+00  2.8406E+01  2.5756E+01  1.2790E+01  3.1354E-01  1.3903E+01  1.3584E+01
             2.5164E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1110.79420820075        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0815E+00  6.2512E-01  2.0433E+01  1.4071E+00  2.3849E+00  2.4638E+00  7.0891E+00  1.2014E+00  1.5780E+00  6.5096E-01
             9.4206E+00
 PARAMETER:  1.7835E-01 -3.6981E-01  3.1172E+00  4.4154E-01  9.6918E-01  1.0017E+00  2.0586E+00  2.8351E-01  5.5615E-01 -3.2931E-01
             2.3429E+00
 GRADIENT:  -1.7516E+00 -1.1444E+01 -7.5520E-01 -1.2025E+01 -3.7142E+00  1.2106E+01 -1.5456E-01  4.5894E-02  8.8715E+00  4.6878E+00
             6.3709E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1114.77088978078        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0662E+00  1.4300E+00  8.1600E+00  1.0772E+00  2.3114E+00  2.3389E+00  5.8111E+00  9.3157E-01  1.1158E+00  4.8261E-01
             9.1422E+00
 PARAMETER:  1.6406E-01  4.5770E-01  2.1992E+00  1.7437E-01  9.3783E-01  9.4967E-01  1.8598E+00  2.9115E-02  2.0959E-01 -6.2855E-01
             2.3129E+00
 GRADIENT:  -3.6968E+00  1.6645E+00 -6.7525E-01  7.1009E+00 -6.2186E-01 -4.9824E+00  2.4496E-01  4.9128E-02  1.1029E+00  2.5257E+00
            -9.8588E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1115.03497481809        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0809E+00  1.6563E+00  5.4163E+00  9.2570E-01  2.2414E+00  2.3832E+00  5.4029E+00  7.4620E-01  8.3342E-01  3.7359E-01
             9.1854E+00
 PARAMETER:  1.7776E-01  6.0458E-01  1.7894E+00  2.2798E-02  9.0709E-01  9.6846E-01  1.7869E+00 -1.9276E-01 -8.2223E-02 -8.8460E-01
             2.3176E+00
 GRADIENT:   2.2232E-01  2.9930E-01  2.7208E-01 -1.2629E+00 -3.8504E-01  1.8131E-01 -4.5298E-01  6.7400E-02  3.3410E-01  1.6174E+00
            -4.2028E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1115.13003608910        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.0816E+00  1.7605E+00  3.7548E+00  8.5786E-01  2.1463E+00  2.3946E+00  5.2245E+00  5.3323E-01  6.4840E-01  2.5122E-01
             9.1988E+00
 PARAMETER:  1.7843E-01  6.6563E-01  1.4230E+00 -5.3312E-02  8.6374E-01  9.7321E-01  1.7534E+00 -5.2880E-01 -3.3325E-01 -1.2814E+00
             2.3191E+00
 GRADIENT:   5.9834E-02 -5.7749E-01 -3.6673E-01 -2.1578E+00  1.1587E+00  1.3120E+00  2.5728E-01  7.3710E-02  1.6356E-01  7.7685E-01
             1.0318E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1115.27097170381        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.0807E+00  1.7981E+00  3.2350E+00  8.4136E-01  2.0907E+00  2.3732E+00  5.1920E+00  3.5049E-01  4.9810E-01  1.2331E-01
             9.2188E+00
 PARAMETER:  1.7759E-01  6.8674E-01  1.2740E+00 -7.2740E-02  8.3748E-01  9.6424E-01  1.7471E+00 -9.4841E-01 -5.9696E-01 -1.9931E+00
             2.3212E+00
 GRADIENT:  -4.4719E-01  1.0035E-01  2.6258E-02  2.1693E+00 -1.2294E+00 -1.3414E+00 -3.1563E-03  3.5424E-02 -5.8889E-01  1.7932E-01
            -1.8275E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1116.01820170288        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      770
 NPARAMETR:  1.0859E+00  1.4307E+00  8.4644E+00  1.0565E+00  2.3884E+00  2.4309E+00  6.1960E+00  7.9637E-01  9.5420E-01  2.5793E-01
             9.2428E+00
 PARAMETER:  1.8237E-01  4.5819E-01  2.2359E+00  1.5495E-01  9.7063E-01  9.8825E-01  1.9239E+00 -1.2769E-01  5.3115E-02 -1.2551E+00
             2.3238E+00
 GRADIENT:   5.4600E-01  1.2635E+00 -9.5828E-01 -4.5412E-01  8.0191E+00  4.9374E+00  4.2661E+00  2.8406E-02 -9.4217E-01  6.5797E-01
             3.2703E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1116.65317091679        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      946
 NPARAMETR:  1.0817E+00  1.0787E+00  1.9679E+01  1.2441E+00  2.5050E+00  2.3776E+00  6.7541E+00  6.6985E-01  1.2662E+00  8.2893E-02
             9.2165E+00
 PARAMETER:  1.7852E-01  1.7578E-01  3.0795E+00  3.1844E-01  1.0183E+00  9.6610E-01  2.0102E+00 -3.0071E-01  3.3605E-01 -2.3902E+00
             2.3210E+00
 GRADIENT:   6.9069E-01 -2.9664E-01 -5.9125E-02  1.8643E-01  2.5129E-01 -5.4897E-01 -1.0956E-01 -2.2505E-03 -4.7722E-01  5.0744E-02
            -1.0844E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1116.68694447166        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1121
 NPARAMETR:  1.0787E+00  1.0148E+00  2.9468E+01  1.2942E+00  2.5599E+00  2.3792E+00  6.9206E+00  5.2857E-01  1.3435E+00  3.1741E-02
             9.2183E+00
 PARAMETER:  1.7579E-01  1.1464E-01  3.4833E+00  3.5788E-01  1.0400E+00  9.6678E-01  2.0345E+00 -5.3757E-01  3.9526E-01 -3.3502E+00
             2.3212E+00
 GRADIENT:  -2.0702E-02 -1.9187E-02  1.7971E-02 -2.4972E-03 -4.0930E-02  7.4127E-02  1.9602E-02 -6.7757E-04 -2.7704E-02  7.3053E-03
            -4.2188E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1116.68698113789        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1297
 NPARAMETR:  1.0788E+00  1.0179E+00  2.9695E+01  1.2931E+00  2.5611E+00  2.3783E+00  6.9129E+00  5.2130E-01  1.3429E+00  3.0224E-02
             9.2189E+00
 PARAMETER:  1.7589E-01  1.1773E-01  3.4910E+00  3.5705E-01  1.0404E+00  9.6640E-01  2.0334E+00 -5.5143E-01  3.9481E-01 -3.3991E+00
             2.3213E+00
 GRADIENT:   1.8522E-02  1.1791E-02  2.6745E-02  3.9444E-03 -6.3741E-02 -1.5985E-02 -3.6599E-02 -6.5687E-04 -2.5684E-02  6.6270E-03
             2.2920E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1116.68704480131        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1473
 NPARAMETR:  1.0789E+00  1.0268E+00  2.8904E+01  1.2891E+00  2.5581E+00  2.3782E+00  6.8935E+00  5.1377E-01  1.3395E+00  2.9059E-02
             9.2189E+00
 PARAMETER:  1.7595E-01  1.2643E-01  3.4640E+00  3.5391E-01  1.0393E+00  9.6636E-01  2.0306E+00 -5.6598E-01  3.9226E-01 -3.4384E+00
             2.3213E+00
 GRADIENT:   2.6703E-02  8.6976E-02  2.8143E-02  7.5239E-02 -7.8967E-02 -5.7305E-02 -5.8999E-02 -6.7670E-04  6.9142E-03  6.1345E-03
            -3.7435E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1116.68801059053        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1652
 NPARAMETR:  1.0790E+00  1.0525E+00  2.5295E+01  1.2737E+00  2.5418E+00  2.3794E+00  6.8332E+00  4.4259E-01  1.3223E+00  1.8885E-02
             9.2188E+00
 PARAMETER:  1.7607E-01  1.5119E-01  3.3306E+00  3.4190E-01  1.0329E+00  9.6686E-01  2.0218E+00 -7.1512E-01  3.7941E-01 -3.8694E+00
             2.3212E+00
 GRADIENT:  -1.0452E-02  1.6032E-01  6.7222E-03  1.6226E-01 -3.2065E-03 -6.2080E-02 -7.6236E-02 -6.4052E-04  6.8974E-02  2.6019E-03
            -6.5097E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1116.69016180062        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1831
 NPARAMETR:  1.0790E+00  1.0385E+00  2.5865E+01  1.2788E+00  2.5438E+00  2.3798E+00  6.8652E+00  3.5307E-01  1.3255E+00  1.0000E-02
             9.2188E+00
 PARAMETER:  1.7605E-01  1.3777E-01  3.3529E+00  3.4589E-01  1.0337E+00  9.6700E-01  2.0265E+00 -9.4108E-01  3.8176E-01 -4.6035E+00
             2.3212E+00
 GRADIENT:   1.3866E-03  1.2996E-02  2.5702E-03 -8.1796E-02 -2.3258E-02  6.5982E-03  1.0762E-01 -3.7681E-04  1.9740E-02  0.0000E+00
             3.7461E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1116.69019837809        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2005
 NPARAMETR:  1.0790E+00  1.0379E+00  2.5860E+01  1.2792E+00  2.5439E+00  2.3798E+00  6.8676E+00  3.5371E-01  1.3255E+00  1.0000E-02
             9.2187E+00
 PARAMETER:  1.7607E-01  1.3721E-01  3.3527E+00  3.4625E-01  1.0337E+00  9.6700E-01  2.0268E+00 -9.3927E-01  3.8180E-01 -4.6035E+00
             2.3212E+00
 GRADIENT:   5.8897E-03  2.1620E-02 -1.4404E-04 -3.3662E-02 -1.0745E-02  6.8675E-03  1.2257E-01 -3.7680E-04  2.4977E-03  0.0000E+00
            -1.3735E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1116.69022420472        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2182            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0790E+00  1.0369E+00  2.5971E+01  1.2800E+00  2.5444E+00  2.3797E+00  6.8678E+00  3.5492E-01  1.3268E+00  1.0000E-02
             9.2188E+00
 PARAMETER:  1.7604E-01  1.3625E-01  3.3570E+00  3.4689E-01  1.0339E+00  9.6696E-01  2.0268E+00 -9.3588E-01  3.8274E-01 -4.6035E+00
             2.3212E+00
 GRADIENT:   9.8702E-01  1.8324E-01  4.2452E-03  1.0608E+00  3.2016E-01  2.2994E+00  1.6386E+01 -3.6612E-04  1.2439E-01  0.0000E+00
             4.1791E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1116.69023210998        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2361
 NPARAMETR:  1.0790E+00  1.0359E+00  2.5982E+01  1.2800E+00  2.5445E+00  2.3797E+00  6.8685E+00  3.5486E-01  1.3267E+00  1.0000E-02
             9.2187E+00
 PARAMETER:  1.7605E-01  1.3527E-01  3.3574E+00  3.4689E-01  1.0339E+00  9.6697E-01  2.0269E+00 -9.3603E-01  3.8267E-01 -4.6035E+00
             2.3212E+00
 GRADIENT:   2.0922E-03 -1.7267E-02 -1.3340E-03 -3.5680E-02  5.9096E-04  4.8159E-03  3.2553E-02 -3.7244E-04  6.3830E-03  0.0000E+00
             2.4549E-02

0ITERATION NO.:   93    OBJECTIVE VALUE:  -1116.69023759147        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     2455
 NPARAMETR:  1.0790E+00  1.0362E+00  2.5996E+01  1.2801E+00  2.5445E+00  2.3797E+00  6.8684E+00  3.5503E-01  1.3266E+00  1.0000E-02
             9.2188E+00
 PARAMETER:  1.7605E-01  1.3559E-01  3.3579E+00  3.4694E-01  1.0339E+00  9.6698E-01  2.0269E+00 -9.3556E-01  3.8261E-01 -4.6035E+00
             2.3212E+00
 GRADIENT:   2.5925E-03 -3.1378E-03 -5.7263E-04  2.7302E-04 -9.9237E-03  7.8354E-03  3.3337E-02 -3.7431E-04 -3.6091E-03  0.0000E+00
             7.3306E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2455
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4214E-02  3.7081E-02 -1.1456E-04 -7.9327E-02 -1.9744E-06
 SE:             2.8828E-02  2.4165E-02  2.1001E-04  1.4170E-02  1.2594E-04
 N:                     100         100         100         100         100

 P VAL.:         6.2198E-01  1.2491E-01  5.8540E-01  2.1720E-08  9.8749E-01

 ETASHRINKSD(%)  3.4235E+00  1.9045E+01  9.9296E+01  5.2528E+01  9.9578E+01
 ETASHRINKVR(%)  6.7298E+00  3.4463E+01  9.9995E+01  7.7464E+01  9.9998E+01
 EBVSHRINKSD(%)  5.3972E+00  1.3042E+01  9.9182E+01  5.6099E+01  9.9543E+01
 EBVSHRINKVR(%)  1.0503E+01  2.4382E+01  9.9993E+01  8.0727E+01  9.9998E+01
 RELATIVEINF(%)  8.9424E+01  4.0242E+01  1.5022E-03  1.0345E+01  4.6590E-04
 EPSSHRINKSD(%)  6.5536E+00
 EPSSHRINKVR(%)  1.2678E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1116.6902375914706     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       537.39912217694018     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    73.38
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1116.690       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  1.04E+00  2.60E+01  1.28E+00  2.54E+00  2.38E+00  6.87E+00  3.55E-01  1.33E+00  1.00E-02  9.22E+00
 


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
+        1.46E+02
 
 TH 2
+       -4.37E+00  1.58E+01
 
 TH 3
+       -1.11E-02 -4.58E-02  1.53E-03
 
 TH 4
+       -2.22E+00  2.22E+01 -2.94E-02  1.04E+02
 
 TH 5
+       -8.23E-01 -2.69E+00 -1.83E-01 -7.08E+00  3.27E+01
 
 TH 6
+       -6.37E+00  8.01E-01 -2.83E-03  1.07E+00  1.44E+00  2.24E+01
 
 TH 7
+        4.36E-01  2.10E+00 -3.25E-03 -7.28E+00  6.11E-01  8.81E-02  2.41E+00
 
 TH 8
+       -2.08E+00  3.65E+00  2.35E-03  1.82E+00  3.87E-01  3.69E-01 -4.91E-02 -1.40E-01
 
 TH 9
+       -1.36E+00 -3.12E+00 -4.05E-02 -2.66E+01  3.72E+00 -7.69E-01  1.55E+00 -1.85E-01  1.75E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.33E+00 -1.71E+00  1.82E-03 -7.05E+00 -1.96E-01  2.23E+00  2.10E-01 -1.43E-02  2.88E+00  0.00E+00  1.11E+01
 
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
 #CPUT: Total CPU Time in Seconds,       90.558
Stop Time:
Sat Sep 18 07:22:08 CDT 2021
