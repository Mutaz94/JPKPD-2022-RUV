Sat Sep 25 14:20:53 CDT 2021
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
$DATA ../../../../data/spa/D/dat45.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF DATA RECS IN DATA SET:      500
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

 TOT. NO. OF OBS RECS:      400
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12600.0435022734        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4184E+02  2.4856E+02 -6.2210E+01  2.1001E+02  1.5845E+02 -1.7719E+03 -7.9214E+02 -4.4293E+01 -1.2405E+03 -3.3647E+02
            -2.4156E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -559.016143033992        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2574E+00  1.0518E+00  9.9113E-01  1.4604E+00  1.3065E+00  2.0024E+00  1.3401E+00  9.6778E-01  1.4092E+00  1.0483E+00
             1.4279E+01
 PARAMETER:  3.2907E-01  1.5053E-01  9.1086E-02  4.7871E-01  3.6734E-01  7.9434E-01  3.9274E-01  6.7249E-02  4.4301E-01  1.4716E-01
             2.7588E+00
 GRADIENT:   7.4605E+00  3.0725E+00 -6.9110E+00  6.7491E-01 -5.0897E+00  4.9193E+01  4.5594E-01  3.5535E+00  1.0398E+01  2.2938E+00
             1.4683E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -567.501735922994        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2603E+00  7.6882E-01  1.2684E+00  1.8567E+00  1.8365E+00  1.8600E+00  3.8180E+00  6.6567E-01  1.6126E+00  1.8645E+00
             1.3042E+01
 PARAMETER:  3.3134E-01 -1.6290E-01  3.3773E-01  7.1878E-01  7.0787E-01  7.2059E-01  1.4397E+00 -3.0696E-01  5.7783E-01  7.2297E-01
             2.6682E+00
 GRADIENT:   1.7244E+01  1.5220E+01 -6.5345E+00  3.2458E+01 -3.2077E+00  4.2441E+00  1.0847E+01  1.1617E+00  2.6782E+01  1.4346E+00
             9.8979E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -587.594705104583        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0402E+00  6.2013E-01  1.6308E+00  1.4918E+00  2.7368E+00  1.7243E+00  2.4812E+00  3.5944E-01  8.6832E-01  7.0220E+00
             1.2086E+01
 PARAMETER:  1.3938E-01 -3.7782E-01  5.8909E-01  5.0000E-01  1.1068E+00  6.4482E-01  1.0088E+00 -9.2321E-01 -4.1201E-02  2.0490E+00
             2.5920E+00
 GRADIENT:  -3.5302E+01  5.0031E+00  8.6554E+00 -1.1766E+01 -1.5064E+01  2.5031E+01  2.7794E+00  1.1118E-02  6.0698E-01  6.5122E+00
             7.9471E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -602.220984296564        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0340E+00  3.1969E-01  8.7770E-01  1.5042E+00  3.8182E+00  1.5019E+00  3.6465E+00  3.6388E-02  6.9770E-01  7.3545E+00
             1.0330E+01
 PARAMETER:  1.3342E-01 -1.0404E+00 -3.0455E-02  5.0827E-01  1.4398E+00  5.0670E-01  1.3938E+00 -3.2135E+00 -2.5997E-01  2.0953E+00
             2.4350E+00
 GRADIENT:   1.8042E+01  8.7106E+00  1.6113E+00  1.4822E+01 -5.7153E+00 -1.4682E+01  1.5138E+00 -9.3605E-04 -8.2360E+00 -1.6171E+00
            -5.0791E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -626.744979977309        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  8.7561E-01  1.1615E-01  3.9461E-01  1.4715E+00  2.7766E+01  1.6948E+00  5.3374E+00  1.0000E-02  8.4845E-01  1.2321E+01
             1.0394E+01
 PARAMETER: -3.2839E-02 -2.0528E+00 -8.2986E-01  4.8629E-01  3.4238E+00  6.2758E-01  1.7747E+00 -7.9091E+00 -6.4340E-02  2.6113E+00
             2.4412E+00
 GRADIENT:  -4.8455E+01  5.4051E+00  3.9993E+00  1.1600E+02  1.0029E+00 -1.7173E+01  3.8406E-01  0.0000E+00 -5.6930E+00  4.0066E+00
            -1.2183E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -654.649534892006        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  6.5208E-01  2.9286E-02  8.3906E-02  6.7691E-01  1.8631E+03  1.6816E+00  3.3784E+00  1.0000E-02  5.2939E-01  1.6368E+01
             1.0153E+01
 PARAMETER: -3.2758E-01 -3.4307E+00 -2.3781E+00 -2.9022E-01  7.6300E+00  6.1972E-01  1.3174E+00 -1.8942E+01 -5.3604E-01  2.8954E+00
             2.4178E+00
 GRADIENT:   8.7179E+00  1.8489E-01 -1.9322E+01  4.3608E+01  1.4248E-03  1.0884E+01  2.0796E-02  0.0000E+00 -1.5161E+01 -1.1264E-05
             1.2616E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -667.516565223091        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  4.3776E-01  1.3085E-02  2.9808E-02  3.2506E-01  3.4253E+04  1.4404E+00  1.9053E+00  1.0000E-02  1.0253E+00  1.2127E+01
             8.8944E+00
 PARAMETER: -7.2609E-01 -4.2363E+00 -3.4130E+00 -1.0237E+00  1.0542E+01  4.6491E-01  7.4465E-01 -2.6098E+01  1.2495E-01  2.5954E+00
             2.2854E+00
 GRADIENT:  -2.1825E+01  2.6597E-01  2.6441E+01 -1.4554E+01 -7.3300E-05 -2.4764E+01  6.7835E-03  0.0000E+00 -1.9879E-01 -1.2703E-08
            -1.3757E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -668.336913860460        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  4.6761E-01  1.5223E-02  3.4628E-02  3.5763E-01  2.3893E+04  1.5136E+00  1.9942E+00  1.0000E-02  9.8646E-01  1.2259E+01
             9.1672E+00
 PARAMETER: -6.6012E-01 -4.0849E+00 -3.2631E+00 -9.2825E-01  1.0181E+01  5.1447E-01  7.9027E-01 -2.5159E+01  8.6362E-02  2.6062E+00
             2.3156E+00
 GRADIENT:  -2.0161E+01  2.7158E-01  3.2796E+01 -3.9394E+01 -9.1243E-05 -9.6429E+00  1.1419E-02  0.0000E+00 -3.7851E+00 -2.2159E-08
             1.6941E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -670.092234176351        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  5.3450E-01  1.6867E-02  4.7828E-02  4.6383E-01  7.8064E+03  1.6025E+00  2.4243E+00  1.0000E-02  1.0138E+00  1.1123E+01
             9.0442E+00
 PARAMETER: -5.2643E-01 -3.9824E+00 -2.9401E+00 -6.6824E-01  9.0627E+00  5.7157E-01  9.8555E-01 -2.2797E+01  1.1372E-01  2.5090E+00
             2.3021E+00
 GRADIENT:   6.2734E-02 -8.0987E-02  3.0350E-02  2.7801E-02  6.5090E-05 -1.0530E-01  5.9660E-03  0.0000E+00  4.8446E-02 -2.0776E-07
             2.8725E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -670.099480202359        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1042
 NPARAMETR:  5.3546E-01  1.9799E-02  4.7909E-02  4.6426E-01  6.5535E+03  1.6029E+00  2.2898E+00  1.0000E-02  1.0141E+00  1.1140E+01
             9.0388E+00
 PARAMETER: -5.2464E-01 -3.8221E+00 -2.9385E+00 -6.6731E-01  8.8878E+00  5.7179E-01  9.2846E-01 -2.2508E+01  1.1402E-01  2.5105E+00
             2.3015E+00
 GRADIENT:  -1.2334E-02  4.7051E-02 -1.6740E-01  1.7424E-01  6.1100E-05  2.0467E-02  1.4676E-02  0.0000E+00 -2.4320E-02 -2.7719E-07
            -1.2401E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -670.101396059602        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1181
 NPARAMETR:  5.3517E-01  1.8801E-02  4.7908E-02  4.6419E-01  6.7155E+03  1.6031E+00  1.6909E+00  1.0000E-02  1.0143E+00  1.1146E+01
             9.0410E+00
 PARAMETER: -5.2516E-01 -3.8739E+00 -2.9385E+00 -6.6746E-01  8.9122E+00  5.7195E-01  6.2529E-01 -2.2550E+01  1.1420E-01  2.5111E+00
             2.3018E+00
 GRADIENT:   3.9240E+00 -4.0966E-02  4.1325E+00  1.6779E+00  6.9031E-05  1.0396E+00  4.5174E-03  0.0000E+00  1.1424E-01 -2.6186E-07
             1.9432E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -670.105633219548        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1280
 NPARAMETR:  5.3367E-01  2.0040E-02  4.7532E-02  4.6169E-01  6.5330E+03  1.6016E+00  1.0000E-02  1.0000E-02  1.0151E+00  1.1320E+01
             9.0332E+00
 PARAMETER: -5.2797E-01 -3.8100E+00 -2.9463E+00 -6.7286E-01  8.8846E+00  5.7102E-01 -5.0135E+00 -2.2550E+01  1.1501E-01  2.5266E+00
             2.3009E+00
 GRADIENT:  -5.7162E-01 -6.7690E-03  6.4093E-02  5.8201E-02  6.0046E-05 -5.5381E-02  0.0000E+00  0.0000E+00  8.2742E-02 -2.8640E-07
            -2.8344E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -670.106348434975        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1455
 NPARAMETR:  5.3487E-01  2.0284E-02  4.7714E-02  4.6291E-01  6.6242E+03  1.6022E+00  1.2072E-01  1.0000E-02  1.0144E+00  1.1227E+01
             9.0387E+00
 PARAMETER: -5.2572E-01 -3.7979E+00 -2.9425E+00 -6.7022E-01  8.8985E+00  5.7135E-01 -2.0143E+00 -2.2550E+01  1.1432E-01  2.5183E+00
             2.3015E+00
 GRADIENT:   1.0256E-02 -3.9118E-03 -1.0953E-01  6.2898E-02  6.3874E-05 -6.6847E-03  4.4744E-05  0.0000E+00 -3.4909E-03 -2.7330E-07
            -5.5565E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -670.106366603970        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1632
 NPARAMETR:  5.3498E-01  2.0360E-02  4.7744E-02  4.6309E-01  6.6295E+03  1.6022E+00  1.7676E-01  1.0000E-02  1.0144E+00  1.1215E+01
             9.0391E+00
 PARAMETER: -5.2552E-01 -3.7942E+00 -2.9419E+00 -6.6984E-01  8.8993E+00  5.7138E-01 -1.6330E+00 -2.2550E+01  1.1430E-01  2.5173E+00
             2.3016E+00
 GRADIENT:  -3.2182E-02  8.4365E-05 -4.3058E-02 -8.8862E-03  6.2868E-05 -9.0928E-03  1.0009E-04  0.0000E+00 -1.3949E-03 -2.7156E-07
            -1.9391E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -670.106431308747        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1812
 NPARAMETR:  5.3560E-01  2.0427E-02  4.7888E-02  4.6410E-01  6.1694E+03  1.6026E+00  2.3638E-01  1.0000E-02  1.0143E+00  1.1203E+01
             9.0394E+00
 PARAMETER: -5.2437E-01 -3.7909E+00 -2.9389E+00 -6.6766E-01  8.8274E+00  5.7162E-01 -1.3423E+00 -2.2550E+01  1.1423E-01  2.5162E+00
             2.3016E+00
 GRADIENT:   5.2206E-02 -2.9525E-04 -1.3029E-01  1.2943E-01  6.9636E-05 -8.5790E-03  1.8075E-04  0.0000E+00 -1.1772E-02 -3.1397E-07
            -6.7188E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -670.106552188259        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1987
 NPARAMETR:  5.3603E-01  2.0508E-02  4.8007E-02  4.6483E-01  5.2201E+03  1.6029E+00  2.6993E-02  1.0000E-02  1.0143E+00  1.1265E+01
             9.0405E+00
 PARAMETER: -5.2356E-01 -3.7870E+00 -2.9364E+00 -6.6608E-01  8.6603E+00  5.7182E-01 -3.5122E+00 -2.2550E+01  1.1420E-01  2.5217E+00
             2.3017E+00
 GRADIENT:  -6.2402E-03  3.1994E-04  1.3610E-02 -1.2766E-02  8.1767E-05  1.1634E-03  2.4077E-06  0.0000E+00  1.2456E-03 -4.4349E-07
             7.2837E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -670.121659094193        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2171
 NPARAMETR:  5.3586E-01  2.0227E-02  4.7905E-02  4.6441E-01  3.8519E+01  1.6027E+00  1.0000E-02  1.0000E-02  1.0142E+00  1.4890E+01
             9.0378E+00
 PARAMETER: -5.2388E-01 -3.8007E+00 -2.9385E+00 -6.6699E-01  3.7512E+00  5.7168E-01 -1.0965E+02 -2.2550E+01  1.1411E-01  2.8007E+00
             2.3014E+00
 GRADIENT:  -6.1254E-01  3.7367E-02  2.0858E-01 -1.2718E-01 -5.1564E-03 -1.3658E-01  0.0000E+00  0.0000E+00  8.4844E-02  3.4701E-02
             2.2697E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -670.134368548913        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2350
 NPARAMETR:  5.3934E-01  1.9480E-02  4.8617E-02  4.6964E-01  2.0304E+01  1.6058E+00  1.0000E-02  1.0000E-02  1.0133E+00  1.0218E+01
             9.0395E+00
 PARAMETER: -5.1742E-01 -3.8384E+00 -2.9238E+00 -6.5578E-01  3.1108E+00  5.7361E-01 -1.1350E+02 -2.2550E+01  1.1322E-01  2.4241E+00
             2.3016E+00
 GRADIENT:  -7.8595E-02  2.7525E-02 -2.1359E-01  4.3822E-01 -1.8550E-02 -5.3476E-02  0.0000E+00  0.0000E+00  6.9073E-02  5.1481E-02
            -1.5152E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -670.140357945176        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2525
 NPARAMETR:  5.3914E-01  1.7786E-02  4.8497E-02  4.6896E-01  1.4361E+01  1.6066E+00  1.0000E-02  1.0000E-02  1.0117E+00  7.5274E+00
             9.0412E+00
 PARAMETER: -5.1778E-01 -3.9293E+00 -2.9263E+00 -6.5723E-01  2.7645E+00  5.7410E-01 -1.2159E+02 -2.2550E+01  1.1161E-01  2.1185E+00
             2.3018E+00
 GRADIENT:   8.0895E-02 -5.1759E-03  3.0940E-02 -5.7766E-02 -1.1090E-02  1.3545E-02  0.0000E+00  0.0000E+00 -2.5615E-02  1.3888E-02
             8.8120E-03

0ITERATION NO.:   99    OBJECTIVE VALUE:  -670.140499851388        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     2652
 NPARAMETR:  5.3893E-01  1.7862E-02  4.8444E-02  4.6864E-01  1.4202E+01  1.6064E+00  1.0000E-02  1.0000E-02  1.0119E+00  7.3278E+00
             9.0406E+00
 PARAMETER: -5.1817E-01 -3.9251E+00 -2.9274E+00 -6.5791E-01  2.7534E+00  5.7399E-01 -1.2244E+02 -2.2550E+01  1.1179E-01  2.0917E+00
             2.3017E+00
 GRADIENT:  -1.8671E-04  9.8378E-05 -8.3743E-03  1.0749E-02  1.4925E-04 -9.3119E-04  0.0000E+00  0.0000E+00  5.1545E-04 -9.6418E-05
            -3.6137E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2652
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9504E-03 -3.7477E-07  7.0317E-05 -2.3075E-02 -8.1832E-04
 SE:             2.8692E-02  1.5078E-06  1.7955E-04  2.4012E-02  8.0981E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4580E-01  8.0370E-01  6.9534E-01  3.3657E-01  3.1225E-01

 ETASHRINKSD(%)  3.8793E+00  9.9995E+01  9.9398E+01  1.9555E+01  9.7287E+01
 ETASHRINKVR(%)  7.6081E+00  1.0000E+02  9.9996E+01  3.5286E+01  9.9926E+01
 EBVSHRINKSD(%)  3.3701E+00  9.9994E+01  9.9423E+01  1.9676E+01  9.7774E+01
 EBVSHRINKVR(%)  6.6267E+00  1.0000E+02  9.9997E+01  3.5480E+01  9.9950E+01
 RELATIVEINF(%)  5.8434E+00  6.1322E-08  5.5508E-05  9.9266E-01  5.6851E-03
 EPSSHRINKSD(%)  1.4060E+01
 EPSSHRINKVR(%)  2.6143E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -670.14049985138831     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       65.010326712349865     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.84
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -670.140       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.39E-01  1.79E-02  4.84E-02  4.69E-01  1.42E+01  1.61E+00  1.00E-02  1.00E-02  1.01E+00  7.33E+00  9.04E+00
 


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
+        1.37E+03
 
 TH 2
+       -3.71E+02  1.63E+03
 
 TH 3
+       -4.40E+03 -1.94E+02  1.49E+05
 
 TH 4
+       -1.01E+02  1.39E+02 -1.95E+04  3.01E+03
 
 TH 5
+        1.74E-01 -9.64E-01 -1.45E+00  1.67E-01  2.79E-03
 
 TH 6
+        1.90E+00  1.25E+01 -4.26E+01 -2.12E+01  4.88E-03  6.41E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.23E-01 -1.75E+01  5.76E+02 -8.19E+01 -7.14E-03 -3.21E-01  0.00E+00  0.00E+00  8.41E+01
 
 TH10
+        1.11E-02  6.81E-01 -2.93E-01  2.37E-02 -3.54E-03  3.91E-03  0.00E+00  0.00E+00 -2.36E-03  6.05E-03
 
 TH11
+       -1.76E+01  4.37E+00  1.23E+02 -1.14E+01 -2.98E-03  1.54E+00  0.00E+00  0.00E+00  2.90E+00 -1.30E-03  4.58E+00
 
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
 #CPUT: Total CPU Time in Seconds,       45.766
Stop Time:
Sat Sep 25 14:21:46 CDT 2021
