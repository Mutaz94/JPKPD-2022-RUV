Wed Sep 29 18:26:22 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat74.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1692.90779519264        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3104E+02 -2.2390E+01 -1.4208E+01 -9.5376E+00 -1.1899E+01  5.6978E+01 -7.2215E+00  1.4487E+01  1.0091E+00  1.7147E+01
             2.3459E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1699.25027099131        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9327E-01  1.0589E+00  1.0852E+00  1.0203E+00  1.0639E+00  9.3688E-01  1.0445E+00  9.4108E-01  1.0122E+00  9.4859E-01
             9.3599E-01
 PARAMETER:  9.3251E-02  1.5723E-01  1.8179E-01  1.2011E-01  1.6193E-01  3.4795E-02  1.4352E-01  3.9270E-02  1.1214E-01  4.7225E-02
             3.3848E-02
 GRADIENT:  -1.2857E-01  1.0492E+00  1.8780E-02 -3.4516E-01 -7.9110E+00 -4.8884E+00 -6.5924E+00  6.2800E+00 -1.6207E+00 -4.2053E+00
            -9.9917E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1700.64765637926        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.9306E-01  1.0195E+00  1.1538E+00  1.0485E+00  1.1038E+00  9.5562E-01  1.1457E+00  7.3451E-01  9.7969E-01  1.0480E+00
             9.6429E-01
 PARAMETER:  9.3038E-02  1.1935E-01  2.4303E-01  1.4733E-01  1.9876E-01  5.4602E-02  2.3605E-01 -2.0855E-01  7.9481E-02  1.4686E-01
             6.3634E-02
 GRADIENT:   9.8905E-01  4.3562E-01  8.9127E-01 -1.3159E-01  7.7218E+00  3.5977E+00 -1.3548E+00  7.8140E-02 -2.2173E+00  3.3981E-01
             1.9848E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1701.46286557720        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  9.9223E-01  8.3726E-01  1.0356E+00  1.1582E+00  9.4821E-01  9.4499E-01  1.3794E+00  5.1698E-01  9.0104E-01  9.2296E-01
             9.5457E-01
 PARAMETER:  9.2204E-02 -7.7624E-02  1.3495E-01  2.4684E-01  4.6819E-02  4.3416E-02  4.2167E-01 -5.5975E-01 -4.2023E-03  1.9830E-02
             5.3503E-02
 GRADIENT:  -3.7196E-02  6.8305E+00  4.1294E+00  6.3542E+00 -5.3213E+00 -5.8346E-01  1.0856E-01 -3.1548E-01 -8.6499E-01 -7.9515E-01
            -1.3800E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1701.58109120883        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.9016E-01  6.4293E-01  1.0660E+00  1.2754E+00  8.8987E-01  9.4446E-01  1.6492E+00  4.8191E-01  8.4821E-01  9.1506E-01
             9.5961E-01
 PARAMETER:  9.0115E-02 -3.4172E-01  1.6388E-01  3.4323E-01 -1.6682E-02  4.2861E-02  6.0028E-01 -6.3000E-01 -6.4623E-02  1.1235E-02
             5.8767E-02
 GRADIENT:  -9.0616E-02  4.0149E+00  3.1192E+00  5.7893E+00 -4.0074E+00  1.9446E-01 -2.6640E-01 -9.9480E-01 -5.9587E-01 -1.4767E-01
             3.0226E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1701.59356295510        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.8898E-01  5.7129E-01  1.1034E+00  1.3252E+00  8.8051E-01  9.4354E-01  1.7809E+00  5.2622E-01  8.3149E-01  9.1877E-01
             9.6073E-01
 PARAMETER:  8.8917E-02 -4.5986E-01  1.9837E-01  3.8153E-01 -2.7250E-02  4.1883E-02  6.7714E-01 -5.4203E-01 -8.4532E-02  1.5280E-02
             5.9936E-02
 GRADIENT:  -2.2048E-01  6.6438E+00  6.0154E+00  1.2871E+01 -8.7878E+00  2.9519E-01 -3.1455E-01 -1.1550E+00 -1.0155E+00 -2.0445E-01
             5.5815E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1701.60634231000        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  9.8813E-01  5.2814E-01  1.1391E+00  1.3562E+00  8.8138E-01  9.4287E-01  1.8717E+00  5.8167E-01  8.2257E-01  9.2362E-01
             9.6086E-01
 PARAMETER:  8.8061E-02 -5.3839E-01  2.3028E-01  4.0466E-01 -2.6263E-02  4.1177E-02  7.2686E-01 -4.4186E-01 -9.5322E-02  2.0549E-02
             6.0076E-02
 GRADIENT:  -3.3300E-01  7.9056E+00  7.4923E+00  1.7379E+01 -1.1827E+01  3.3873E-01 -2.9147E-01 -1.0062E+00 -1.1677E+00 -1.2260E-01
             7.9356E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.62862582684        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1231
 NPARAMETR:  9.8683E-01  4.6473E-01  1.2074E+00  1.4005E+00  8.9105E-01  9.4154E-01  2.0260E+00  6.8839E-01  8.1046E-01  9.3365E-01
             9.5953E-01
 PARAMETER:  8.6740E-02 -6.6631E-01  2.8849E-01  4.3680E-01 -1.5354E-02  3.9760E-02  8.0604E-01 -2.7340E-01 -1.1015E-01  3.1349E-02
             5.8690E-02
 GRADIENT:  -4.0409E-01  8.1044E+00  7.7545E+00  2.0227E+01 -1.3560E+01  2.7343E-01 -1.1036E-01 -3.0373E-01 -1.1457E+00 -2.3144E-02
             7.4400E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.70077290530        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1406
 NPARAMETR:  9.8505E-01  3.6973E-01  1.2952E+00  1.4611E+00  9.0150E-01  9.3939E-01  2.3187E+00  7.9784E-01  7.9371E-01  9.4747E-01
             9.5744E-01
 PARAMETER:  8.4938E-02 -8.9500E-01  3.5863E-01  4.7918E-01 -3.6937E-03  3.7480E-02  9.4102E-01 -1.2585E-01 -1.3103E-01  4.6041E-02
             5.6506E-02
 GRADIENT:  -3.1579E-01  6.1402E+00  6.3530E+00  1.6899E+01 -1.1035E+01  8.5202E-02  1.3826E-01  2.3298E-01 -9.6462E-01 -2.2261E-01
             1.0741E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1701.85852723186        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1581
 NPARAMETR:  9.8252E-01  2.3060E-01  1.3977E+00  1.5487E+00  9.0478E-01  9.3697E-01  3.0021E+00  9.0923E-01  7.7295E-01  9.6067E-01
             9.5736E-01
 PARAMETER:  8.2366E-02 -1.3670E+00  4.3482E-01  5.3740E-01 -6.6419E-05  3.4898E-02  1.1993E+00  4.8479E-03 -1.5754E-01  5.9879E-02
             5.6429E-02
 GRADIENT:  -3.7576E-01  3.8969E+00  4.6282E+00  1.3106E+01 -8.0330E+00  4.0769E-02  2.1191E-01  3.5283E-01 -4.9731E-01 -1.4708E-01
             1.6342E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1702.24547269638        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1759
 NPARAMETR:  9.8017E-01  9.9207E-02  1.5063E+00  1.6339E+00  9.1033E-01  9.3484E-01  4.5758E+00  1.0208E+00  7.5184E-01  9.6905E-01
             9.5767E-01
 PARAMETER:  7.9970E-02 -2.2105E+00  5.0966E-01  5.9097E-01  6.0536E-03  3.2617E-02  1.6208E+00  1.2054E-01 -1.8523E-01  6.8561E-02
             5.6748E-02
 GRADIENT:  -5.6091E-01  3.3260E+00  2.7943E+00  1.0434E+01 -5.2831E+00 -6.8339E-02  3.6978E+00 -7.3004E-02 -2.3709E+00 -9.0751E-01
            -1.1402E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1702.85062462620        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1937
 NPARAMETR:  9.7959E-01  4.8220E-02  1.5499E+00  1.6579E+00  9.1771E-01  9.3415E-01  6.0456E+00  1.0631E+00  7.3818E-01  9.7559E-01
             9.5762E-01
 PARAMETER:  7.9379E-02 -2.9320E+00  5.3817E-01  6.0553E-01  1.4129E-02  3.1882E-02  1.8993E+00  1.6124E-01 -2.0356E-01  7.5288E-02
             5.6691E-02
 GRADIENT:   1.1051E-01  2.4481E+00 -9.6893E-01 -8.9634E+00  2.4309E+00 -3.1827E-03  5.7690E+00 -6.8830E-01 -2.5569E+00 -2.7551E+00
            -4.7193E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1702.98714339091        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     2102
 NPARAMETR:  9.7959E-01  4.1577E-02  1.4493E+00  1.6511E+00  8.8276E-01  9.3409E-01  6.2389E+00  9.7476E-01  7.4163E-01  9.5691E-01
             9.5659E-01
 PARAMETER:  7.9382E-02 -3.0802E+00  4.7105E-01  6.0147E-01 -2.4702E-02  3.1820E-02  1.9308E+00  7.4438E-02 -1.9890E-01  5.5957E-02
             5.5622E-02
 GRADIENT:   7.0751E-01 -1.6378E+00 -1.8571E+00 -3.3428E+00  1.4391E+00  5.7089E-02 -2.7739E+00  7.3687E-01  1.5519E+00  1.4081E+00
             2.7961E-01

0ITERATION NO.:   64    OBJECTIVE VALUE:  -1703.00440539593        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     2235
 NPARAMETR:  9.7930E-01  4.1599E-02  1.4542E+00  1.6522E+00  8.8219E-01  9.3390E-01  6.2396E+00  9.7293E-01  7.4177E-01  9.5410E-01
             9.5651E-01
 PARAMETER:  7.8625E-02 -3.0763E+00  4.7440E-01  6.0164E-01 -2.5415E-02  3.1487E-02  1.9325E+00  7.1552E-02 -1.9852E-01  5.2789E-02
             5.5027E-02
 GRADIENT:  -2.1947E+00  1.1392E+01 -1.9057E-01 -5.2201E+01 -4.5804E-01 -1.9444E-01  1.4658E+01 -1.7179E-01  4.7764E-01 -2.3814E-01
            -4.3215E-01
 NUMSIGDIG:         1.8         2.4         3.4         2.5         2.6         2.3         2.5         1.4         2.4         2.1
                    1.7

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2235
 NO. OF SIG. DIGITS IN FINAL EST.:  1.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3478E-03  1.3038E-02 -2.4787E-02 -1.5669E-02 -3.1154E-02
 SE:             2.9847E-02  9.1784E-03  1.5490E-02  2.8585E-02  2.1311E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6398E-01  1.5545E-01  1.0955E-01  5.8358E-01  1.4378E-01

 ETASHRINKSD(%)  1.0112E-02  6.9251E+01  4.8108E+01  4.2379E+00  2.8605E+01
 ETASHRINKVR(%)  2.0222E-02  9.0545E+01  7.3072E+01  8.2962E+00  4.9028E+01
 EBVSHRINKSD(%)  4.3422E-01  8.0062E+01  5.0086E+01  3.7519E+00  2.5006E+01
 EBVSHRINKVR(%)  8.6654E-01  9.6025E+01  7.5086E+01  7.3630E+00  4.3759E+01
 RELATIVEINF(%)  9.8965E+01  1.9388E+00  4.6062E+00  4.7382E+01  1.0177E+01
 EPSSHRINKSD(%)  4.3985E+01
 EPSSHRINKVR(%)  6.8623E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1703.0044053959305     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -967.85357883219228     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.35
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1703.004       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  4.17E-02  1.45E+00  1.65E+00  8.82E-01  9.34E-01  6.25E+00  9.72E-01  7.42E-01  9.54E-01  9.56E-01
 


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
+        1.32E+03
 
 TH 2
+       -2.62E+01  6.38E+04
 
 TH 3
+       -8.12E+00 -7.33E+01  2.40E+03
 
 TH 4
+       -8.89E+00 -1.94E+02 -2.15E+01  1.82E+03
 
 TH 5
+        7.63E+00  9.05E+01 -1.79E+04 -7.06E+01  1.21E+03
 
 TH 6
+        3.41E+00 -6.99E+00  5.54E+00 -2.31E+00 -5.18E+00  2.25E+02
 
 TH 7
+        3.89E-02 -1.84E+00 -9.91E-01 -4.56E+00  2.43E+00 -1.39E-02  8.20E+00
 
 TH 8
+       -1.25E+00 -6.51E+01  8.08E+01  1.96E+00 -8.31E+01  1.57E+00 -4.64E-01  4.83E+01
 
 TH 9
+       -2.85E+00 -6.27E+01  1.05E+04 -2.62E+01 -3.21E+02  4.44E+00  3.37E-02  8.19E+01  6.69E+02
 
 TH10
+       -4.49E+00 -2.49E+02  1.63E+04  1.46E+01 -3.54E+02  4.89E+00 -1.89E+00  9.01E+01  2.82E+02  3.18E+02
 
 TH11
+       -1.07E+01 -2.62E+01  4.05E+01 -1.00E+01 -5.02E+01  2.61E+00 -9.43E-02  2.73E+01  6.08E+01  6.29E+01  2.40E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.220
Stop Time:
Wed Sep 29 18:27:01 CDT 2021
