Sat Sep 25 02:54:59 CDT 2021
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
$DATA ../../../../data/int/SL3/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -34.7213541798119        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3231E+01 -1.1850E+02  1.3114E+02  1.5483E+02  1.8460E+02  3.6610E+01 -1.4742E+02 -2.3240E+02 -1.3488E+02 -4.9505E+01
            -7.0576E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2354.83887008476        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0786E+00  1.4006E+00  8.7771E-01  9.2345E-01  9.6811E-01  8.3773E-01  1.1840E+00  1.0699E+00  9.9092E-01  1.0859E+00
             5.1655E+00
 PARAMETER:  1.7563E-01  4.3692E-01 -3.0441E-02  2.0358E-02  6.7587E-02 -7.7061E-02  2.6888E-01  1.6760E-01  9.0875E-02  1.8242E-01
             1.7420E+00
 GRADIENT:   3.1531E+01  6.8790E+01  1.1559E+01  2.0676E+01 -6.8975E+01 -8.4382E+00  1.9662E+01  4.4300E+00  1.7740E+01  3.9306E+00
             7.9554E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2412.09045980897        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.0473E+00  2.0423E+00  5.5782E-01  6.2591E-01  1.1143E+00  9.2526E-01  9.1041E-01  4.5587E+00  2.4885E+00  1.5385E+00
             4.2658E+00
 PARAMETER:  1.4625E-01  8.1407E-01 -4.8372E-01 -3.6856E-01  2.0825E-01  2.2317E-02  6.1427E-03  1.6170E+00  1.0117E+00  5.3082E-01
             1.5506E+00
 GRADIENT:   9.3766E+00  1.2438E+02 -7.3789E+00  7.8344E+01 -1.7330E+02  1.9309E+01  3.1403E+01  2.9737E+01  6.3556E+01 -1.4232E+01
             6.1995E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2610.41778827914        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9818E-01  1.8345E+00  8.9360E-01  5.5567E-01  1.4382E+00  8.6575E-01  6.9468E-01  7.8703E+00  9.4439E-01  1.3439E+00
             2.7722E+00
 PARAMETER:  9.8176E-02  7.0678E-01 -1.2500E-02 -4.8759E-01  4.6342E-01 -4.4163E-02 -2.6430E-01  2.1631E+00  4.2781E-02  3.9557E-01
             1.1197E+00
 GRADIENT:  -5.7072E+01 -7.8315E-01 -1.1001E+01  4.2815E+01 -3.4558E+01 -1.1899E+01 -4.0508E+01  6.7603E+01 -1.0752E+01 -5.3815E+01
            -1.4271E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2665.99872048260        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0126E+00  1.8304E+00  1.5228E+00  5.4360E-01  1.5699E+00  8.8642E-01  8.1604E-01  3.7459E+00  1.1801E+00  1.5737E+00
             2.7622E+00
 PARAMETER:  1.1255E-01  7.0453E-01  5.2052E-01 -5.0953E-01  5.5102E-01 -2.0564E-02 -1.0329E-01  1.4207E+00  2.6560E-01  5.5340E-01
             1.1160E+00
 GRADIENT:  -1.1190E+01 -2.1776E+00 -7.9775E+00  6.9684E+00  5.2829E+00 -2.3789E+00 -8.5094E-01 -2.1913E+00  1.2001E+00 -1.1597E+00
            -4.8166E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2670.52104640627        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0175E+00  1.7855E+00  2.9650E+00  5.7326E-01  1.6741E+00  8.9031E-01  8.6324E-01  4.9347E+00  9.6169E-01  1.6508E+00
             2.7659E+00
 PARAMETER:  1.1739E-01  6.7971E-01  1.1869E+00 -4.5641E-01  6.1528E-01 -1.6184E-02 -4.7067E-02  1.6963E+00  6.0942E-02  6.0126E-01
             1.1174E+00
 GRADIENT:   2.8002E+00 -8.9237E+00 -1.4668E+00 -8.0283E+00  6.3515E+00 -2.4646E-01  8.6788E-01  1.2348E+00 -4.1987E-01  1.4412E+00
             3.8145E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2670.55471137968        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:      507
 NPARAMETR:  1.0174E+00  1.7853E+00  2.9755E+00  5.7464E-01  1.6725E+00  8.9078E-01  8.5709E-01  4.9417E+00  9.8046E-01  1.6491E+00
             2.7690E+00
 PARAMETER:  1.1723E-01  6.7961E-01  1.1904E+00 -4.5402E-01  6.1433E-01 -1.5662E-02 -5.4211E-02  1.6977E+00  8.0264E-02  6.0025E-01
             1.1185E+00
 GRADIENT:   2.1795E+00 -7.6540E+00 -1.4142E+00 -6.9233E+00  6.0593E+00 -5.5984E-02  2.1059E-01  1.5924E+00 -2.3250E-02  1.6382E+00
             6.3251E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2670.55493233406        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      580
 NPARAMETR:  1.0174E+00  1.7854E+00  2.9755E+00  5.7464E-01  1.6725E+00  8.9086E-01  8.5574E-01  4.9417E+00  9.8332E-01  1.6491E+00
             2.7689E+00
 PARAMETER:  1.1723E-01  6.7962E-01  1.1904E+00 -4.5402E-01  6.1433E-01 -1.5570E-02 -5.5783E-02  1.6977E+00  8.3182E-02  6.0025E-01
             1.1185E+00
 GRADIENT:   2.1231E+00 -7.7876E+00 -1.4105E+00 -6.8872E+00  6.1015E+00 -4.0963E-02  2.0107E-02  1.6204E+00  1.3410E-02  1.6751E+00
             6.3224E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2670.55559065886        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  1.0174E+00  1.7854E+00  2.9755E+00  5.7464E-01  1.6725E+00  8.9077E-01  8.5519E-01  4.9416E+00  9.8440E-01  1.6491E+00
             2.7689E+00
 PARAMETER:  1.1723E-01  6.7964E-01  1.1904E+00 -4.5401E-01  6.1431E-01 -1.5665E-02 -5.6429E-02  1.6977E+00  8.4277E-02  6.0025E-01
             1.1184E+00
 GRADIENT:   2.1762E+00 -7.7639E+00 -1.4078E+00 -6.8429E+00  6.1091E+00 -5.7362E-02 -5.8640E-02  1.6291E+00  2.2718E-02  1.6891E+00
             6.2666E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2670.64565220570        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  1.0173E+00  1.7932E+00  2.9783E+00  5.7566E-01  1.6677E+00  8.9125E-01  8.5062E-01  4.9297E+00  9.9832E-01  1.6478E+00
             2.7566E+00
 PARAMETER:  1.1719E-01  6.8400E-01  1.1914E+00 -4.5224E-01  6.1147E-01 -1.5134E-02 -6.1788E-02  1.6953E+00  9.8316E-02  5.9942E-01
             1.1140E+00
 GRADIENT:   2.2017E+00  4.0309E+00 -1.1021E+00 -1.9638E+00  4.0692E+00 -3.5320E-02 -4.2131E-01  1.3789E+00  1.1564E-01  1.7969E+00
            -4.3321E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2670.65289738668        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      912
 NPARAMETR:  1.0173E+00  1.7956E+00  2.9792E+00  5.7596E-01  1.6663E+00  8.9275E-01  8.5411E-01  4.9262E+00  9.9295E-01  1.6474E+00
             2.7530E+00
 PARAMETER:  1.1718E-01  6.8531E-01  1.1916E+00 -4.5171E-01  6.1062E-01 -1.3447E-02 -5.7691E-02  1.6946E+00  9.2928E-02  5.9917E-01
             1.1127E+00
 GRADIENT:  -4.7087E+00 -7.1420E+00 -1.1944E+00 -2.7641E+00  8.0181E-01 -4.2443E-02  2.0640E-02  7.4194E-01 -1.4816E-02  1.0157E+00
            -9.4652E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2670.68602035703        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1095
 NPARAMETR:  1.0172E+00  1.8014E+00  2.9783E+00  5.7651E-01  1.6646E+00  8.9267E-01  8.5307E-01  4.9281E+00  9.9421E-01  1.6460E+00
             2.7650E+00
 PARAMETER:  1.1710E-01  6.8857E-01  1.1913E+00 -4.5076E-01  6.0959E-01 -1.3544E-02 -5.8919E-02  1.6949E+00  9.4193E-02  5.9834E-01
             1.1170E+00
 GRADIENT:  -5.6124E+00  2.2273E-01 -1.1249E+00  4.6724E-01 -2.5146E-01 -6.7102E-02  2.3708E-02  7.3875E-01 -4.5118E-02  1.0135E+00
            -4.8407E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2670.71387989123        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1283
 NPARAMETR:  1.0193E+00  1.8016E+00  3.0524E+00  5.7631E-01  1.6647E+00  8.9279E-01  8.5272E-01  4.9071E+00  9.9668E-01  1.6419E+00
             2.7656E+00
 PARAMETER:  1.1915E-01  6.8869E-01  1.2159E+00 -4.5112E-01  6.0962E-01 -1.3402E-02 -5.9321E-02  1.6907E+00  9.6673E-02  5.9585E-01
             1.1173E+00
 GRADIENT:  -1.1279E-01  5.5445E-01 -4.3312E-01 -7.4261E-01 -1.5680E+00  2.5642E-02  3.8530E-02 -1.1473E-01 -1.6048E-02  9.5767E-02
             1.5903E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2670.71426300130        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1438
 NPARAMETR:  1.0192E+00  1.8005E+00  3.0505E+00  5.7702E-01  1.6642E+00  8.9264E-01  8.5238E-01  4.9113E+00  9.9489E-01  1.6408E+00
             2.7653E+00
 PARAMETER:  1.1905E-01  6.8806E-01  1.2153E+00 -4.4987E-01  6.0931E-01 -1.3576E-02 -5.9717E-02  1.6915E+00  9.4880E-02  5.9518E-01
             1.1171E+00
 GRADIENT:   6.7572E+00  1.5445E+01 -3.3060E-01  1.4179E+00  9.6172E-01  5.5840E-01 -1.9510E-02  5.0681E-01 -4.5476E-02  7.0688E-01
             1.8598E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2670.71441419224        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1512
 NPARAMETR:  1.0188E+00  1.7979E+00  3.0505E+00  5.7772E-01  1.6641E+00  8.9241E-01  8.5266E-01  4.9113E+00  9.9457E-01  1.6397E+00
             2.7648E+00
 PARAMETER:  1.1866E-01  6.8665E-01  1.2153E+00 -4.4867E-01  6.0931E-01 -1.3829E-02 -5.9393E-02  1.6915E+00  9.4556E-02  5.9451E-01
             1.1170E+00
 GRADIENT:   5.7423E+00  1.3565E+01 -3.8622E-01  7.2740E-01  1.3248E+00  4.7151E-01 -4.7660E-02  5.9785E-01 -1.9799E-02  6.3867E-01
             1.6265E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2670.71451280680        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1587
 NPARAMETR:  1.0186E+00  1.7955E+00  3.0505E+00  5.7887E-01  1.6641E+00  8.9230E-01  8.5344E-01  4.9113E+00  9.9284E-01  1.6386E+00
             2.7646E+00
 PARAMETER:  1.1846E-01  6.8531E-01  1.2153E+00 -4.4667E-01  6.0931E-01 -1.3948E-02 -5.8476E-02  1.6915E+00  9.2813E-02  5.9381E-01
             1.1169E+00
 GRADIENT:   5.2369E+00  1.2585E+01 -4.6730E-01  4.7617E-01  1.7111E+00  4.2988E-01 -3.5152E-02  6.9825E-01 -2.0731E-02  5.9078E-01
             1.4966E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2670.71453324858        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1663
 NPARAMETR:  1.0186E+00  1.7943E+00  3.0505E+00  5.7956E-01  1.6641E+00  8.9227E-01  8.5383E-01  4.9113E+00  9.9197E-01  1.6379E+00
             2.7645E+00
 PARAMETER:  1.1839E-01  6.8461E-01  1.2153E+00 -4.4548E-01  6.0931E-01 -1.3990E-02 -5.8023E-02  1.6915E+00  9.1936E-02  5.9344E-01
             1.1169E+00
 GRADIENT:   5.0571E+00  1.2203E+01 -5.1416E-01  4.1929E-01  1.9249E+00  4.1484E-01 -3.5726E-02  7.5549E-01 -2.0656E-02  5.7321E-01
             1.4482E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2670.71725106091        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:     1787
 NPARAMETR:  1.0192E+00  1.7933E+00  3.0505E+00  5.8112E-01  1.6641E+00  8.9275E-01  8.5670E-01  4.9113E+00  9.9265E-01  1.6365E+00
             2.7652E+00
 PARAMETER:  1.1898E-01  6.8407E-01  1.2153E+00 -4.4279E-01  6.0931E-01 -1.3444E-02 -5.4664E-02  1.6915E+00  9.2618E-02  5.9258E-01
             1.1171E+00
 GRADIENT:  -5.1702E-01 -1.3304E+00 -7.6385E-01 -9.1432E-01 -2.2646E-01  1.2593E-02  4.3331E-01  3.9139E-01  1.2205E-01 -1.0915E-01
             2.3897E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2670.72766795053        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1944
 NPARAMETR:  1.0194E+00  1.7924E+00  3.0881E+00  5.8262E-01  1.6641E+00  8.9274E-01  8.5555E-01  4.9029E+00  9.9043E-01  1.6365E+00
             2.7653E+00
 PARAMETER:  1.1919E-01  6.8355E-01  1.2275E+00 -4.4022E-01  6.0929E-01 -1.3461E-02 -5.6016E-02  1.6898E+00  9.0382E-02  5.9256E-01
             1.1172E+00
 GRADIENT:   7.1106E+00  1.5004E+01 -3.8261E-01  1.5746E+00  1.7307E+00  5.9537E-01  1.2663E-01  5.6408E-01  3.5095E-02  4.6919E-01
             1.9492E+00

0ITERATION NO.:   93    OBJECTIVE VALUE:  -2670.72841411848        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:     2026
 NPARAMETR:  1.0194E+00  1.7927E+00  3.0907E+00  5.8273E-01  1.6639E+00  8.9274E-01  8.5564E-01  4.9052E+00  9.9052E-01  1.6363E+00
             2.7662E+00
 PARAMETER:  1.1917E-01  6.8352E-01  1.2288E+00 -4.4017E-01  6.0934E-01 -1.3469E-02 -5.6009E-02  1.6898E+00  9.0427E-02  5.9260E-01
             1.1171E+00
 GRADIENT:  -3.9373E-02 -1.9132E+04  1.0642E+04 -1.4852E+04  2.1455E+04 -4.3798E-03 -1.6789E-02 -7.7574E+03 -1.2919E-03  2.2061E+04
            -1.1705E+04
 NUMSIGDIG:         3.9         3.3         3.3         3.3         3.3         3.9         2.8         3.3         3.1         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2026
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6994E-03 -1.6423E-02 -3.9594E-02  1.4133E-02 -3.0735E-02
 SE:             2.9237E-02  2.4629E-02  1.5614E-02  1.6300E-02  2.4542E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5365E-01  5.0487E-01  1.1218E-02  3.8590E-01  2.1044E-01

 ETASHRINKSD(%)  2.0515E+00  1.7491E+01  4.7691E+01  4.5394E+01  1.7782E+01
 ETASHRINKVR(%)  4.0610E+00  3.1923E+01  7.2638E+01  7.0181E+01  3.2402E+01
 EBVSHRINKSD(%)  2.1164E+00  1.7694E+01  5.4020E+01  5.0939E+01  1.3227E+01
 EBVSHRINKVR(%)  4.1880E+00  3.2257E+01  7.8859E+01  7.5930E+01  2.4705E+01
 RELATIVEINF(%)  9.5717E+01  5.7613E+00  9.7102E+00  1.9444E+00  4.6779E+01
 EPSSHRINKSD(%)  1.7218E+01
 EPSSHRINKVR(%)  3.1472E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          883
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1622.8454496394520     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2670.7284141184782     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1047.8829644790262     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    51.22
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2670.728       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.79E+00  3.09E+00  5.83E-01  1.66E+00  8.93E-01  8.56E-01  4.90E+00  9.90E-01  1.64E+00  2.77E+00
 


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
+        1.29E+03
 
 TH 2
+       -1.36E+02  2.18E+06
 
 TH 3
+        3.85E+01 -3.73E+02  2.27E+05
 
 TH 4
+       -5.96E+02  1.04E+07 -4.81E+02  4.97E+07
 
 TH 5
+        1.38E+02 -1.43E+03  5.22E+02 -1.73E+03  3.18E+06
 
 TH 6
+        1.25E+01 -2.30E+02  7.32E+01 -1.09E+03  2.73E+02  2.30E+02
 
 TH 7
+        4.65E+00  4.27E+02 -1.35E+02  1.97E+03 -5.03E+02 -9.20E+00  1.21E+02
 
 TH 8
+       -1.76E+01  1.70E+02 -6.56E+01  2.24E+02 -3.90E+05 -3.34E+01  6.17E+01  4.79E+04
 
 TH 9
+       -7.89E+00  1.50E+02 -5.00E+01  7.45E+02 -1.85E+02 -2.19E+00  2.41E+01  2.43E+01  1.89E+01
 
 TH10
+        1.49E+02 -1.47E+03  6.49E+00 -1.82E+03  1.73E+01  2.86E+02 -5.25E+02 -3.24E+00 -1.90E+02  3.48E+06
 
 TH11
+       -6.56E+01  8.64E+05 -1.08E+02  4.13E+06 -4.04E+02 -8.62E+01  1.70E+02  4.95E+01  6.35E+01 -1.09E+06  3.43E+05
 
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
 #CPUT: Total CPU Time in Seconds,       66.252
Stop Time:
Sat Sep 25 02:56:07 CDT 2021
