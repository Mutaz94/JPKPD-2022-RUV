Wed Sep 29 01:46:44 CDT 2021
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
$DATA ../../../../data/int/SL1/dat16.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3339.51235828344        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7476E+02  1.1432E+01 -2.0689E+01  1.1822E+02  1.8428E+02  8.0743E+01 -6.4095E+01 -5.7099E+01 -1.1341E+01 -1.0564E+01
            -1.0547E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3525.85820794656        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      108
 NPARAMETR:  1.0308E+00  1.0183E+00  9.9560E-01  9.9749E-01  9.3677E-01  8.4172E-01  1.1752E+00  1.0251E+00  9.5716E-01  9.9737E-01
             1.3932E+00
 PARAMETER:  1.3035E-01  1.1818E-01  9.5590E-02  9.7483E-02  3.4683E-02 -7.2304E-02  2.6145E-01  1.2481E-01  5.6218E-02  9.7364E-02
             4.3157E-01
 GRADIENT:  -1.7934E+01 -2.4363E+01 -1.8577E+01  5.2554E-02  9.5407E+00 -1.8135E+01 -1.7075E+01 -4.0979E+00 -2.6666E+00  6.3840E+00
            -6.4931E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3527.73079062063        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      286
 NPARAMETR:  1.0274E+00  1.0876E+00  1.1153E+00  9.7195E-01  1.0285E+00  8.6113E-01  1.2248E+00  1.2031E+00  1.0135E+00  9.5393E-01
             1.3927E+00
 PARAMETER:  1.2704E-01  1.8393E-01  2.0914E-01  7.1551E-02  1.2810E-01 -4.9514E-02  3.0275E-01  2.8493E-01  1.1338E-01  5.2836E-02
             4.3124E-01
 GRADIENT:  -2.7431E+01 -1.3603E+01 -1.2772E+01  1.6214E+01  4.5730E+01 -9.0708E+00 -5.8443E+00 -4.7229E+00  7.4199E+00 -9.8961E+00
            -5.7946E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3531.15947778043        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  1.0370E+00  1.1865E+00  1.2904E+00  9.2895E-01  1.1167E+00  8.8015E-01  1.1293E+00  1.7333E+00  8.9762E-01  1.1579E+00
             1.4327E+00
 PARAMETER:  1.3630E-01  2.7102E-01  3.5493E-01  2.6305E-02  2.1035E-01 -2.7666E-02  2.2160E-01  6.5003E-01 -8.0108E-03  2.4665E-01
             4.5956E-01
 GRADIENT:  -7.7703E-01 -3.2635E+00 -5.6937E+00 -3.7392E+00 -1.8542E+00  1.8553E-01 -3.8861E-02  1.3687E+00  1.0361E+00 -4.4862E-01
             8.7336E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3532.75263131700        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      631
 NPARAMETR:  1.0356E+00  1.2624E+00  1.5662E+00  8.9598E-01  1.2313E+00  8.7933E-01  1.1017E+00  2.4077E+00  7.9709E-01  1.2651E+00
             1.4232E+00
 PARAMETER:  1.3494E-01  3.3303E-01  5.4865E-01 -9.8347E-03  3.0803E-01 -2.8591E-02  1.9685E-01  9.7866E-01 -1.2679E-01  3.3511E-01
             4.5290E-01
 GRADIENT:   2.9878E+02  1.6460E+02  4.3844E+00  2.1903E+01  6.6075E+01  1.7474E+01  8.2536E+00  8.8278E+00  6.0322E+00  1.0390E+01
             1.7010E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3533.45275179862        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      769
 NPARAMETR:  1.0295E+00  1.2164E+00  1.5663E+00  9.2621E-01  1.2149E+00  8.8325E-01  1.1367E+00  2.3898E+00  7.5725E-01  1.2420E+00
             1.4159E+00
 PARAMETER:  1.2903E-01  2.9590E-01  5.4869E-01  2.3349E-02  2.9466E-01 -2.4143E-02  2.2809E-01  9.7121E-01 -1.7806E-01  3.1672E-01
             4.4774E-01
 GRADIENT:  -2.1129E+01  2.5944E+00 -1.2495E+01 -5.8044E+00  8.0010E+00  1.2291E+00 -2.7022E+00 -2.8459E+00 -1.3653E-01  1.0081E+00
             3.4973E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3533.64360595655        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  1.0376E+00  1.2012E+00  1.5720E+00  9.3478E-01  1.2153E+00  8.7373E-01  1.1657E+00  2.3841E+00  7.5329E-01  1.2384E+00
             1.4188E+00
 PARAMETER:  1.3695E-01  2.8328E-01  5.5237E-01  3.2555E-02  2.9501E-01 -3.4984E-02  2.5333E-01  9.6884E-01 -1.8330E-01  3.1383E-01
             4.4978E-01
 GRADIENT:   1.2987E+00 -8.8403E-01 -1.4124E+01 -2.9186E+00  1.5049E+01 -2.8125E+00  4.8601E-01 -3.8244E+00  8.1089E-01  2.1151E+00
             9.5003E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3533.93819693139        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1085
 NPARAMETR:  1.0345E+00  1.1720E+00  1.6066E+00  9.5758E-01  1.2115E+00  8.9421E-01  1.1985E+00  2.4057E+00  7.1739E-01  1.2369E+00
             1.4129E+00
 PARAMETER:  1.3394E-01  2.5868E-01  5.7411E-01  5.6655E-02  2.9186E-01 -1.1810E-02  2.8104E-01  9.7784E-01 -2.3214E-01  3.1265E-01
             4.4567E-01
 GRADIENT:  -6.8693E+00  6.8210E-01 -1.5763E+01  1.1244E+01  2.2409E+01  6.0968E+00  4.2301E-01 -7.4790E+00 -1.3440E+00  4.1887E+00
             5.6532E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3534.32354814274        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     1207
 NPARAMETR:  1.0344E+00  1.1409E+00  1.6096E+00  9.6246E-01  1.1718E+00  8.7683E-01  1.1818E+00  2.3818E+00  7.2829E-01  1.1895E+00
             1.4093E+00
 PARAMETER:  1.3381E-01  2.3179E-01  5.7600E-01  6.1742E-02  2.5856E-01 -3.1444E-02  2.6708E-01  9.6787E-01 -2.1706E-01  2.7352E-01
             4.4310E-01
 GRADIENT:   2.9943E+02  8.4166E+01  2.8352E+00  2.7375E+01  6.0244E+01  1.6739E+01  9.5764E+00  2.9953E+00  1.8360E+00  9.8976E+00
             1.0068E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3534.59249732960        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1388
 NPARAMETR:  1.0351E+00  1.1440E+00  1.6095E+00  9.6640E-01  1.1699E+00  8.7740E-01  1.1932E+00  2.4263E+00  7.3501E-01  1.1806E+00
             1.4083E+00
 PARAMETER:  1.3449E-01  2.3452E-01  5.7592E-01  6.5818E-02  2.5695E-01 -3.0796E-02  2.7662E-01  9.8635E-01 -2.0787E-01  2.6600E-01
             4.4236E-01
 GRADIENT:  -5.1904E+00 -3.9591E+00 -1.2957E+01 -1.1264E+01 -1.0012E-01 -1.1727E+00 -2.6646E+00 -4.6977E+00 -1.3150E+00  1.4035E+00
             2.5879E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3534.61936393234        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1575
 NPARAMETR:  1.0352E+00  1.1444E+00  1.6093E+00  9.6717E-01  1.1702E+00  8.7854E-01  1.2059E+00  2.4244E+00  7.3646E-01  1.1790E+00
             1.4088E+00
 PARAMETER:  1.3461E-01  2.3492E-01  5.7580E-01  6.6624E-02  2.5719E-01 -2.9492E-02  2.8725E-01  9.8558E-01 -2.0590E-01  2.6464E-01
             4.4274E-01
 GRADIENT:  -4.8342E+00 -2.0572E+00 -1.3182E+01 -9.5498E+00  6.7046E-01 -6.7019E-01 -6.2614E-01 -4.7570E+00 -7.1273E-01  1.4858E+00
             3.5637E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -3534.62455680186        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:     1660
 NPARAMETR:  1.0353E+00  1.1442E+00  1.6092E+00  9.6709E-01  1.1705E+00  8.8087E-01  1.2113E+00  2.4229E+00  7.3799E-01  1.1787E+00
             1.4093E+00
 PARAMETER:  1.3460E-01  2.3494E-01  5.7591E-01  6.6639E-02  2.5716E-01 -2.6448E-02  2.9193E-01  9.8594E-01 -2.0588E-01  2.6465E-01
             4.4266E-01
 GRADIENT:  -4.0698E+02  4.5762E+02  1.6113E+02  1.0736E+03 -4.1819E+02  3.7740E-01  1.8865E-01  1.0707E+02 -4.8123E-01  2.0574E+02
            -2.4360E+02
 NUMSIGDIG:         2.3         2.3         2.8         2.3         2.3         1.7         2.5         2.3         1.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1660
 NO. OF SIG. DIGITS IN FINAL EST.:  1.3
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1309E-03 -1.4922E-02 -1.6695E-02  1.6215E-02 -3.3214E-02
 SE:             2.9763E-02  2.4319E-02  2.2244E-02  2.2479E-02  2.3188E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1622E-01  5.3949E-01  4.5291E-01  4.7071E-01  1.5204E-01

 ETASHRINKSD(%)  2.8978E-01  1.8528E+01  2.5481E+01  2.4692E+01  2.2317E+01
 ETASHRINKVR(%)  5.7872E-01  3.3624E+01  4.4469E+01  4.3286E+01  3.9653E+01
 EBVSHRINKSD(%)  5.8715E-01  1.8718E+01  2.9172E+01  2.8514E+01  1.8345E+01
 EBVSHRINKVR(%)  1.1709E+00  3.3932E+01  4.9833E+01  4.8898E+01  3.3324E+01
 RELATIVEINF(%)  9.8819E+01  2.8021E+01  3.8675E+01  2.1460E+01  4.0558E+01
 EPSSHRINKSD(%)  2.0691E+01
 EPSSHRINKVR(%)  3.7101E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3534.6245568018562     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1880.5351970334455     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    52.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3534.625       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.14E+00  1.61E+00  9.67E-01  1.17E+00  8.81E-01  1.21E+00  2.43E+00  7.36E-01  1.18E+00  1.41E+00
 


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
+        1.41E+05
 
 TH 2
+       -7.20E+04  3.76E+04
 
 TH 3
+        1.10E+01  2.88E+01  9.05E+03
 
 TH 4
+       -2.01E+05  1.04E+05 -6.83E+01  2.91E+05
 
 TH 5
+        4.60E+00  4.27E+00 -2.99E+02 -2.25E+02  3.00E+04
 
 TH 6
+        1.52E+01 -4.77E+00  2.71E+00 -1.18E+01  2.29E+00  2.51E+02
 
 TH 7
+       -4.44E+01  5.07E+01 -2.47E+01  4.82E+01 -1.84E+01 -2.57E-01  6.32E+01
 
 TH 8
+       -9.51E-01 -1.22E+01 -3.89E+01  3.66E+01  1.93E+01 -6.02E-01  2.24E+00  5.00E+02
 
 TH 9
+        5.04E+02  6.60E+04  5.38E+01 -7.01E+02  2.58E+02 -2.91E-01  2.16E+01 -2.09E+01  1.17E+05
 
 TH10
+       -4.09E+01  3.21E+04  9.66E+01  9.22E+01 -4.46E+01 -2.68E+00  2.69E+01 -1.21E+00  5.70E+04  2.78E+04
 
 TH11
+       -6.81E+00  2.53E+01 -1.88E+02 -1.70E+02 -1.16E+02  3.72E+00 -3.36E+00 -2.67E+01  1.28E+02  7.67E+00  7.63E+03
 
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
 #CPUT: Total CPU Time in Seconds,       68.219
Stop Time:
Wed Sep 29 01:47:54 CDT 2021
