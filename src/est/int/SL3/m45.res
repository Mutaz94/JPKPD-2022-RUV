Wed Sep 29 04:19:55 CDT 2021
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
$DATA ../../../../data/int/SL3/dat45.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      983
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      883
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1116.98309075205        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0746E+02  1.3546E+02  1.2636E+02  9.4818E+01  2.4688E+02 -1.2332E+00 -1.0269E+02 -6.6978E+02 -2.3791E+02 -2.6910E+01
            -4.2166E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2839.22450104056        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.4530E-01  1.0863E+00  1.1152E+00  9.6083E-01  1.0177E+00  1.2201E+00  9.3068E-01  1.1092E+00  1.0886E+00  8.7416E-01
             2.1952E+00
 PARAMETER:  4.3744E-02  1.8276E-01  2.0900E-01  6.0043E-02  1.1756E-01  2.9892E-01  2.8160E-02  2.0360E-01  1.8487E-01 -3.4495E-02
             8.8627E-01
 GRADIENT:   4.1403E+01  1.0702E+01 -7.3549E+00  2.1048E+01  1.3071E+01  7.8056E+01  1.5735E+00 -5.6133E+00 -3.1571E+00 -1.0846E+01
            -9.6985E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2845.80946494148        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      179
 NPARAMETR:  9.6204E-01  1.2931E+00  1.6879E+00  8.8200E-01  1.3483E+00  1.1806E+00  5.9999E-01  1.4816E+00  1.0285E+00  1.3872E+00
             2.2336E+00
 PARAMETER:  6.1297E-02  3.5707E-01  6.2346E-01 -2.5565E-02  3.9882E-01  2.6606E-01 -4.1084E-01  4.9312E-01  1.2811E-01  4.2727E-01
             9.0363E-01
 GRADIENT:  -3.5927E+01 -2.0547E+00 -2.1207E+01  4.5677E+01  2.7735E+01  5.5764E+00 -9.7951E+00 -3.5427E+00 -4.0329E+01  3.0633E+01
            -7.9498E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2858.48849500552        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      355
 NPARAMETR:  9.8558E-01  1.3691E+00  2.1739E+00  8.2195E-01  1.4949E+00  1.1656E+00  8.3020E-01  2.5432E+00  1.1061E+00  1.2251E+00
             2.2553E+00
 PARAMETER:  8.5472E-02  4.1417E-01  8.7653E-01 -9.6075E-02  5.0203E-01  2.5321E-01 -8.6091E-02  1.0334E+00  2.0083E-01  3.0299E-01
             9.1326E-01
 GRADIENT:   2.3414E+00 -2.3275E+01 -2.9837E+01  3.3740E+00  3.9862E+01  7.9939E-01  1.0781E+01 -3.2687E-01 -6.7117E+00 -1.0576E+00
            -2.0842E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2863.86918638166        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.8892E-01  1.4584E+00  2.6026E+00  7.5925E-01  1.5406E+00  1.1571E+00  6.5785E-01  2.9725E+00  1.2861E+00  1.2083E+00
             2.2827E+00
 PARAMETER:  8.8860E-02  4.7736E-01  1.0565E+00 -1.7542E-01  5.3218E-01  2.4592E-01 -3.1877E-01  1.1894E+00  3.5162E-01  2.8925E-01
             9.2536E-01
 GRADIENT:   7.4588E+00 -5.8661E+00 -1.0123E+01 -9.8478E+00  1.1690E+01 -2.0299E+00  5.1648E+00  1.3325E+00  4.5162E-01 -1.1111E+01
             1.2641E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2864.42007555779        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  9.8903E-01  1.4546E+00  2.6512E+00  7.5817E-01  1.5466E+00  1.1322E+00  5.0252E-01  2.9830E+00  1.3049E+00  1.3083E+00
             2.2918E+00
 PARAMETER:  8.8972E-02  4.7471E-01  1.0750E+00 -1.7684E-01  5.3607E-01  2.2417E-01 -5.8812E-01  1.1929E+00  3.6609E-01  3.6874E-01
             9.2932E-01
 GRADIENT:   1.0035E+02  1.2684E+02 -6.2866E+00  1.2205E+01  4.0324E+01  2.8759E+01  1.2902E+00  3.6455E+00 -4.7071E+00  7.9975E+00
             3.8395E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2864.90391095629        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      846
 NPARAMETR:  9.8903E-01  1.4546E+00  2.6515E+00  7.5818E-01  1.5466E+00  1.1630E+00  4.9124E-01  2.9831E+00  1.3785E+00  1.2863E+00
             2.2914E+00
 PARAMETER:  8.8969E-02  4.7474E-01  1.0751E+00 -1.7683E-01  5.3604E-01  2.5098E-01 -6.1083E-01  1.1930E+00  4.2099E-01  3.5175E-01
             9.2915E-01
 GRADIENT:   7.4655E+00 -3.0886E-01 -1.0381E+01 -4.0033E+00  8.2871E+00  8.7414E-02 -9.2115E-02  2.4125E-01 -2.5762E-01  3.0027E-01
             2.3161E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2866.09844247897        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1028             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8457E-01  1.4463E+00  3.1612E+00  7.7095E-01  1.5584E+00  1.1652E+00  4.8092E-01  3.1488E+00  1.3699E+00  1.2946E+00
             2.2692E+00
 PARAMETER:  8.4445E-02  4.6898E-01  1.2510E+00 -1.6013E-01  5.4365E-01  2.5288E-01 -6.3205E-01  1.2470E+00  4.1476E-01  3.5817E-01
             9.1943E-01
 GRADIENT:   9.5263E+01  1.4669E+02  4.6157E+00  1.9712E+01  2.9193E+01  4.6418E+01  3.1903E+00  2.5577E+00  7.5627E+00  3.3232E+00
             2.0813E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2866.24549305482        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1186
 NPARAMETR:  9.8137E-01  1.4328E+00  3.1965E+00  7.7054E-01  1.5654E+00  1.1564E+00  4.8406E-01  3.1819E+00  1.3663E+00  1.2968E+00
             2.2658E+00
 PARAMETER:  8.1194E-02  4.5962E-01  1.2620E+00 -1.6066E-01  5.4815E-01  2.4532E-01 -6.2555E-01  1.2575E+00  4.1210E-01  3.5988E-01
             9.1791E-01
 GRADIENT:   9.0557E+01  1.2433E+02  3.8222E+00  9.2675E+00  3.6517E+01  4.2620E+01  4.0068E+00  3.5113E+00  8.4876E+00  3.8466E+00
             1.9489E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2866.26673223009        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1329
 NPARAMETR:  9.7889E-01  1.4287E+00  3.1950E+00  7.7185E-01  1.5592E+00  1.1637E+00  4.6950E-01  3.1612E+00  1.3645E+00  1.2955E+00
             2.2697E+00
 PARAMETER:  7.8664E-02  4.5677E-01  1.2616E+00 -1.5897E-01  5.4416E-01  2.5161E-01 -6.5610E-01  1.2509E+00  4.1077E-01  3.5888E-01
             9.1964E-01
 GRADIENT:  -8.5749E+00  5.6510E+00 -8.8391E-01 -1.1746E+01 -4.8745E-01  2.0956E-01  1.0776E-02 -1.2230E+00 -9.5916E-01  1.4705E-01
             7.1779E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2866.38633307323        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1463
 NPARAMETR:  9.7781E-01  1.4114E+00  3.2182E+00  7.8292E-01  1.5502E+00  1.1505E+00  4.4775E-01  3.1397E+00  1.3684E+00  1.2908E+00
             2.2809E+00
 PARAMETER:  7.7560E-02  4.4461E-01  1.2688E+00 -1.4472E-01  5.3839E-01  2.4016E-01 -7.0352E-01  1.2441E+00  4.1366E-01  3.5526E-01
             9.2459E-01
 GRADIENT:   8.3003E+01  1.1729E+02  4.3479E+00  7.5201E+00  3.1678E+01  3.9486E+01  3.1927E+00  2.5947E+00  8.8497E+00  4.1999E+00
             3.4962E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2866.40630340654        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:     1566
 NPARAMETR:  9.7206E-01  1.3907E+00  3.1963E+00  7.9349E-01  1.5341E+00  1.1399E+00  4.2909E-01  3.1375E+00  1.3572E+00  1.2795E+00
             2.2821E+00
 PARAMETER:  7.1664E-02  4.2980E-01  1.2620E+00 -1.3132E-01  5.2792E-01  2.3090E-01 -7.4608E-01  1.2434E+00  4.0543E-01  3.4651E-01
             9.2508E-01
 GRADIENT:  -2.0931E+01  2.7416E+00 -1.6177E+00 -1.3408E+01 -4.3098E+00 -8.1020E+00 -2.5404E-01 -1.2790E+00 -1.3339E-01  3.7116E-01
             2.1013E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2869.51275998389        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1746
 NPARAMETR:  9.8211E-01  1.0318E+00  4.1346E+00  1.0556E+00  1.4280E+00  1.1722E+00  2.9877E-01  3.1132E+00  1.0996E+00  1.1544E+00
             2.2923E+00
 PARAMETER:  8.1944E-02  1.3126E-01  1.5194E+00  1.5411E-01  4.5628E-01  2.5890E-01 -1.1081E+00  1.2357E+00  1.9490E-01  2.4357E-01
             9.2956E-01
 GRADIENT:  -3.3530E+00  7.4997E+00  6.1462E-01  8.7265E+00  2.1295E-01  3.3767E+00 -2.7830E-02 -6.3582E+00 -6.5152E-01  6.4357E-01
             4.3489E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2869.70706883966        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1924
 NPARAMETR:  9.8552E-01  1.0066E+00  4.0605E+00  1.0670E+00  1.4044E+00  1.1799E+00  2.4853E-01  3.1151E+00  1.0926E+00  1.1266E+00
             2.2816E+00
 PARAMETER:  8.5417E-02  1.0660E-01  1.5013E+00  1.6484E-01  4.3961E-01  2.6542E-01 -1.2922E+00  1.2363E+00  1.8854E-01  2.1917E-01
             9.2490E-01
 GRADIENT:   2.5597E+00  3.5752E+00 -2.1706E-01  2.7857E+00 -3.9240E+00  5.9203E+00 -3.6213E-02 -6.5167E+00 -5.5152E-01 -1.0896E+00
             3.4609E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2870.33406320639        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     2080
 NPARAMETR:  9.8276E-01  9.8891E-01  4.1427E+00  1.0745E+00  1.4063E+00  1.1609E+00  1.8399E-01  3.2969E+00  1.0976E+00  1.1374E+00
             2.2435E+00
 PARAMETER:  8.2609E-02  8.8850E-02  1.5213E+00  1.7185E-01  4.4097E-01  2.4923E-01 -1.5929E+00  1.2930E+00  1.9312E-01  2.2879E-01
             9.0802E-01
 GRADIENT:   9.5814E+01  1.0050E+01  8.9071E+00  3.9118E+01  3.1957E+01  4.5730E+01  7.0217E-01  7.1575E+00  6.1240E+00  3.3611E+00
             1.8921E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2870.40598090275        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     2214
 NPARAMETR:  9.8260E-01  9.8997E-01  4.3032E+00  1.0738E+00  1.4061E+00  1.1594E+00  1.9328E-01  3.3773E+00  1.0947E+00  1.1195E+00
             2.2367E+00
 PARAMETER:  8.2452E-02  8.9920E-02  1.5594E+00  1.7121E-01  4.4081E-01  2.4788E-01 -1.5436E+00  1.3171E+00  1.9048E-01  2.1287E-01
             9.0498E-01
 GRADIENT:  -1.2875E+00  2.8135E-01  7.2158E-01 -2.7752E+00 -2.6441E+00 -9.4568E-01  2.3783E-02 -1.4985E+00  4.2587E-01 -1.3971E+00
            -3.1609E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2870.44023656676        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2390
 NPARAMETR:  9.8372E-01  9.7969E-01  4.3595E+00  1.0820E+00  1.4093E+00  1.1628E+00  1.3700E-01  3.4333E+00  1.0935E+00  1.1264E+00
             2.2386E+00
 PARAMETER:  8.3587E-02  7.9478E-02  1.5723E+00  1.7882E-01  4.4308E-01  2.5081E-01 -1.8878E+00  1.3335E+00  1.8936E-01  2.1904E-01
             9.0586E-01
 GRADIENT:   5.1497E-01 -4.9820E-01 -5.7731E-02 -6.4712E-01  6.2426E-01  1.9372E-01  1.2332E-02  1.4866E-01  5.7059E-01  1.3587E-03
             5.2774E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2870.44212046517        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     2549
 NPARAMETR:  9.8353E-01  9.7243E-01  4.3684E+00  1.0868E+00  1.4056E+00  1.1633E+00  1.1420E-01  3.4371E+00  1.0897E+00  1.1235E+00
             2.2382E+00
 PARAMETER:  8.3397E-02  7.2041E-02  1.5744E+00  1.8327E-01  4.4045E-01  2.5127E-01 -2.0698E+00  1.3346E+00  1.8592E-01  2.1645E-01
             9.0569E-01
 GRADIENT:   2.3413E-01 -6.3871E-01 -8.9500E-02 -5.6015E-01 -2.1001E-02  3.8084E-01  8.5163E-03  2.0065E-01  4.5210E-01  7.6542E-02
             6.3070E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2870.44253780645        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2726
 NPARAMETR:  9.8366E-01  9.7272E-01  4.3697E+00  1.0870E+00  1.4056E+00  1.1635E+00  1.1616E-01  3.4365E+00  1.0888E+00  1.1232E+00
             2.2382E+00
 PARAMETER:  8.3521E-02  7.2337E-02  1.5747E+00  1.8342E-01  4.4047E-01  2.5145E-01 -2.0528E+00  1.3344E+00  1.8504E-01  2.1617E-01
             9.0569E-01
 GRADIENT:   4.2263E-01 -3.2718E-02 -3.7987E-02 -4.5771E-02 -1.5103E-01  4.4886E-01  5.6975E-03  1.5973E-01  1.5267E-01  1.6975E-02
            -4.3041E-02

0ITERATION NO.:   94    OBJECTIVE VALUE:  -2870.44267787166        NO. OF FUNC. EVALS.: 110
 CUMULATIVE NO. OF FUNC. EVALS.:     2836
 NPARAMETR:  9.8363E-01  9.7276E-01  4.3688E+00  1.0870E+00  1.4056E+00  1.1636E+00  1.2722E-01  3.4366E+00  1.0880E+00  1.1230E+00
             2.2382E+00
 PARAMETER:  8.3490E-02  7.2379E-02  1.5745E+00  1.8340E-01  4.4046E-01  2.5151E-01 -1.9619E+00  1.3345E+00  1.8430E-01  2.1604E-01
             9.0569E-01
 GRADIENT:   3.7395E-01 -2.0468E-01 -6.0571E-02 -2.1880E-01 -1.0640E-01  4.7012E-01  8.1673E-03  1.9177E-01  1.5996E-01  2.4295E-02
            -1.7647E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2836
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0647E-03 -9.9520E-03 -3.8361E-02 -5.1742E-04 -3.1395E-02
 SE:             2.9670E-02  2.9343E-03  2.0243E-02  2.8574E-02  2.2469E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4452E-01  6.9479E-04  5.8089E-02  9.8555E-01  1.6233E-01

 ETASHRINKSD(%)  6.0208E-01  9.0170E+01  3.2183E+01  4.2730E+00  2.4725E+01
 ETASHRINKVR(%)  1.2005E+00  9.9034E+01  5.4009E+01  8.3634E+00  4.3337E+01
 EBVSHRINKSD(%)  9.3937E-01  9.1720E+01  3.3919E+01  4.3526E+00  2.2428E+01
 EBVSHRINKVR(%)  1.8699E+00  9.9314E+01  5.6333E+01  8.5157E+00  3.9825E+01
 RELATIVEINF(%)  9.8096E+01  1.1754E-01  2.3369E+01  1.7630E+01  3.4461E+01
 EPSSHRINKSD(%)  1.7673E+01
 EPSSHRINKVR(%)  3.2223E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          883
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1622.8454496394520     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2870.4426778716634     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1247.5972282322114     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    79.03
 Elapsed covariance  time in seconds:    13.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2870.443       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  9.73E-01  4.37E+00  1.09E+00  1.41E+00  1.16E+00  1.27E-01  3.44E+00  1.09E+00  1.12E+00  2.24E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.47E-02  2.35E-01  8.74E-01  1.60E-01  1.27E-01  7.18E-02  1.57E+00  4.53E-01  1.86E-01  1.81E-01  2.28E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.20E-03
 
 TH 2
+        7.42E-04  5.53E-02
 
 TH 3
+       -6.84E-04 -5.64E-02  7.64E-01
 
 TH 4
+       -4.25E-04 -3.63E-02  3.64E-02  2.57E-02
 
 TH 5
+        2.07E-04  2.57E-02 -2.84E-03 -1.67E-02  1.61E-02
 
 TH 6
+        5.29E-04 -2.06E-03  8.65E-03  1.14E-03 -1.09E-03  5.16E-03
 
 TH 7
+        7.44E-03  1.28E-01 -1.92E-01 -5.46E-02  7.02E-02 -4.48E-03  2.45E+00
 
 TH 8
+        1.31E-03 -5.57E-03  2.77E-01  1.39E-03  1.90E-04  5.25E-03 -1.49E-01  2.05E-01
 
 TH 9
+       -5.07E-04  2.46E-02 -2.11E-02 -1.96E-02  1.04E-02 -5.59E-04 -1.37E-01  1.20E-02  3.46E-02
 
 TH10
+        8.45E-05  2.56E-02 -1.20E-02 -1.90E-02  1.29E-02 -2.49E-03 -2.09E-03  2.05E-03  1.92E-02  3.26E-02
 
 TH11
+        4.40E-04  3.65E-02 -5.92E-02 -2.34E-02  1.30E-02 -1.06E-03  8.89E-02 -8.09E-03  1.59E-02  8.11E-03  5.20E-02
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.47E-02
 
 TH 2
+        9.09E-02  2.35E-01
 
 TH 3
+       -2.26E-02 -2.75E-01  8.74E-01
 
 TH 4
+       -7.64E-02 -9.64E-01  2.60E-01  1.60E-01
 
 TH 5
+        4.69E-02  8.61E-01 -2.56E-02 -8.19E-01  1.27E-01
 
 TH 6
+        2.12E-01 -1.22E-01  1.38E-01  9.89E-02 -1.19E-01  7.18E-02
 
 TH 7
+        1.37E-01  3.49E-01 -1.40E-01 -2.17E-01  3.53E-01 -3.98E-02  1.57E+00
 
 TH 8
+        8.36E-02 -5.23E-02  7.01E-01  1.92E-02  3.31E-03  1.61E-01 -2.10E-01  4.53E-01
 
 TH 9
+       -7.86E-02  5.63E-01 -1.30E-01 -6.58E-01  4.40E-01 -4.18E-02 -4.72E-01  1.43E-01  1.86E-01
 
 TH10
+        1.35E-02  6.04E-01 -7.64E-02 -6.55E-01  5.62E-01 -1.92E-01 -7.38E-03  2.51E-02  5.73E-01  1.81E-01
 
 TH11
+        5.56E-02  6.80E-01 -2.97E-01 -6.41E-01  4.47E-01 -6.49E-02  2.49E-01 -7.84E-02  3.74E-01  1.97E-01  2.28E-01
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.41E+02
 
 TH 2
+       -4.66E+01  5.60E+02
 
 TH 3
+        6.35E+00  1.17E+01  4.18E+00
 
 TH 4
+        1.24E+01  5.02E+02 -7.44E-01  8.20E+02
 
 TH 5
+        2.10E+01 -2.36E+02 -1.94E+01 -2.56E+01  3.80E+02
 
 TH 6
+       -1.05E+02  2.99E+01 -2.16E+00  4.23E+01  1.19E+01  2.24E+02
 
 TH 7
+        2.45E+00 -1.13E+01  1.61E-01 -2.62E+00 -1.72E+00 -1.99E+00  1.85E+00
 
 TH 8
+       -1.54E+01 -1.22E+01 -5.44E+00  7.13E-01  2.04E+01 -1.75E+00 -2.52E-02  1.25E+01
 
 TH 9
+        7.11E+01 -4.42E+01  4.01E+00  8.01E+01 -1.16E+01 -2.39E+01  1.51E+01 -7.47E+00  1.93E+02
 
 TH10
+       -1.44E+01 -6.83E+00 -2.35E+00  4.85E+01  8.12E+00  2.59E+01 -7.84E-01  3.09E+00 -2.48E+01  7.14E+01
 
 TH11
+        4.16E+00 -6.14E+01 -1.01E+00 -3.98E+00  4.54E+01  4.49E+00 -3.67E-01  1.46E+00 -8.38E+00  2.08E+01  4.83E+01
 
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
 #CPUT: Total CPU Time in Seconds,       92.842
Stop Time:
Wed Sep 29 04:21:30 CDT 2021
