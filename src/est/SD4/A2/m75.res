Sun Oct 24 02:18:57 CDT 2021
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
$DATA ../../../../data/SD4/A2/dat75.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m75.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -809.944892655393        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0866E+02  1.2284E+01  6.6025E+01 -8.6444E+01  1.5508E+02  4.5111E+01 -5.5000E+01 -2.0036E+01 -1.2546E+02 -1.1117E+02
            -1.3640E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1312.61338848187        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0495E+00  8.4897E-01  8.4347E-01  1.1275E+00  7.3909E-01  1.0043E+00  1.1268E+00  9.4727E-01  1.3103E+00  1.1106E+00
             2.0193E+00
 PARAMETER:  1.4831E-01 -6.3737E-02 -7.0230E-02  2.1998E-01 -2.0233E-01  1.0426E-01  2.1934E-01  4.5826E-02  3.7023E-01  2.0494E-01
             8.0273E-01
 GRADIENT:   2.4099E+02  2.5308E+01  4.8365E+01 -1.1960E+01 -3.6550E+01  2.1485E+01 -1.7283E+00  5.5354E+00  1.5596E+01  9.0747E+00
            -1.6012E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1331.66924000010        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      206
 NPARAMETR:  1.0489E+00  5.6652E-01  4.1165E-01  1.3067E+00  3.7698E-01  9.8378E-01  1.2943E+00  5.8619E-01  1.1682E+00  6.5634E-01
             1.9252E+00
 PARAMETER:  1.4772E-01 -4.6825E-01 -7.8759E-01  3.6754E-01 -8.7557E-01  8.3649E-02  3.5794E-01 -4.3411E-01  2.5544E-01 -3.2108E-01
             7.5503E-01
 GRADIENT:   8.4038E+01  1.1704E+02  1.5992E+02  6.7863E+01 -2.2301E+02  2.9792E+00 -8.7253E+00 -3.0765E+00 -4.0705E+00  3.2156E+00
            -1.7179E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1369.22291625420        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      383
 NPARAMETR:  1.0320E+00  4.7897E-01  2.2340E-01  1.2532E+00  2.8025E-01  1.0007E+00  8.9384E-01  1.0756E-01  1.3215E+00  4.5207E-01
             2.1638E+00
 PARAMETER:  1.3153E-01 -6.3611E-01 -1.3988E+00  3.2568E-01 -1.1721E+00  1.0074E-01 -1.2231E-02 -2.1297E+00  3.7878E-01 -6.9391E-01
             8.7188E-01
 GRADIENT:   4.8821E+01  3.3398E+01  1.1385E+00  1.2152E+02 -3.2401E+01  5.1496E+00 -1.1566E+01 -4.9201E-01  2.8538E+00 -1.1721E+01
            -7.2532E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1382.84032590882        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      560
 NPARAMETR:  1.0142E+00  3.8292E-01  2.4885E-01  1.1895E+00  2.7678E-01  9.8011E-01  1.3911E+00  9.5661E-02  1.1866E+00  4.4340E-01
             2.4104E+00
 PARAMETER:  1.1406E-01 -8.5993E-01 -1.2909E+00  2.7357E-01 -1.1845E+00  7.9914E-02  4.3012E-01 -2.2469E+00  2.7107E-01 -7.1328E-01
             9.7979E-01
 GRADIENT:  -9.1258E-01  7.2606E+00  1.4399E+01  1.4798E+00 -2.2660E+01  2.9826E+00  8.0906E-02 -2.9744E-01  1.9331E+00 -3.6803E+00
            -5.5927E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1383.46001609437        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  1.0142E+00  2.9437E-01  2.6526E-01  1.2361E+00  2.7045E-01  9.6799E-01  1.5839E+00  8.9195E-02  1.1473E+00  5.4142E-01
             2.4398E+00
 PARAMETER:  1.1406E-01 -1.1229E+00 -1.2270E+00  3.1200E-01 -1.2077E+00  6.7466E-02  5.5991E-01 -2.3169E+00  2.3740E-01 -5.1357E-01
             9.9190E-01
 GRADIENT:  -3.9397E+00  3.2283E+00  6.0539E+00  1.9956E+00 -1.2451E+01  3.0612E-01 -8.4567E-01 -1.6319E-01 -1.6924E+00  4.4904E-01
             2.2067E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1385.30843355171        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      916
 NPARAMETR:  1.0034E+00  1.4615E-01  3.1299E-01  1.3245E+00  2.8600E-01  9.5783E-01  3.1751E+00  5.9553E-01  1.0760E+00  5.1220E-01
             2.4186E+00
 PARAMETER:  1.0337E-01 -1.8231E+00 -1.0616E+00  3.8101E-01 -1.1518E+00  5.6912E-02  1.2553E+00 -4.1830E-01  1.7324E-01 -5.6905E-01
             9.8319E-01
 GRADIENT:  -9.1775E+00 -7.3703E+00  1.0014E+00 -9.0747E+00  1.3647E+01  9.0186E-01 -7.9094E+00  3.8662E+00  2.0190E+00  1.2135E+01
             1.0511E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1387.01588176162        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1093
 NPARAMETR:  1.0119E+00  1.4124E-01  2.9989E-01  1.3227E+00  2.7534E-01  9.5369E-01  3.1719E+00  7.5057E-01  1.0975E+00  3.6456E-01
             2.3969E+00
 PARAMETER:  1.1184E-01 -1.8573E+00 -1.1044E+00  3.7966E-01 -1.1897E+00  5.2581E-02  1.2543E+00 -1.8693E-01  1.9306E-01 -9.0907E-01
             9.7418E-01
 GRADIENT:   9.8895E+00 -8.2600E+00  1.0941E+01  6.8435E+00 -1.0253E+01 -1.2968E+00 -1.0200E+01  4.7678E+00  5.5191E+00  4.4216E+00
             9.7436E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1388.20381944193        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1268
 NPARAMETR:  1.0115E+00  2.2303E-01  2.9231E-01  1.2878E+00  2.7784E-01  9.6220E-01  2.2696E+00  8.0010E-01  1.1149E+00  2.3773E-01
             2.3570E+00
 PARAMETER:  1.1139E-01 -1.4004E+00 -1.1299E+00  3.5291E-01 -1.1807E+00  6.1470E-02  9.1962E-01 -1.2302E-01  2.0877E-01 -1.3366E+00
             9.5739E-01
 GRADIENT:  -8.2091E-01 -2.4623E+00  8.8255E+00  3.6274E+00 -1.0875E+01  1.4845E-01 -5.4465E+00 -2.4802E-01  1.3660E-01  5.5403E-01
            -3.4533E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1388.38329853970        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1444
 NPARAMETR:  1.0095E+00  1.9723E-01  2.9245E-01  1.2941E+00  2.7601E-01  9.6072E-01  2.6172E+00  8.1071E-01  1.1064E+00  2.5403E-01
             2.3520E+00
 PARAMETER:  1.0943E-01 -1.5234E+00 -1.1295E+00  3.5782E-01 -1.1873E+00  5.9924E-02  1.0621E+00 -1.0985E-01  2.0109E-01 -1.2703E+00
             9.5525E-01
 GRADIENT:  -1.6867E+00  6.3617E-01 -4.9608E-02 -4.1948E-01 -1.0219E+00  2.0060E-01 -2.9890E-01  8.6251E-01 -4.5963E-01  1.0295E+00
            -1.3016E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1388.64188353316        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1621
 NPARAMETR:  1.0082E+00  1.5112E-01  2.9718E-01  1.3124E+00  2.7540E-01  9.5717E-01  3.4329E+00  8.1779E-01  1.0958E+00  1.9933E-01
             2.3497E+00
 PARAMETER:  1.0815E-01 -1.7897E+00 -1.1134E+00  3.7186E-01 -1.1895E+00  5.6225E-02  1.3334E+00 -1.0115E-01  1.9152E-01 -1.5128E+00
             9.5430E-01
 GRADIENT:   3.0870E+00  2.3932E-01  7.6461E+00 -4.1572E-01 -1.2643E+01 -4.9647E-01  5.2167E-01 -1.0046E+00  3.9268E+00  1.8552E-01
            -2.1354E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1388.76040327631        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1796
 NPARAMETR:  1.0065E+00  1.5207E-01  3.0197E-01  1.3185E+00  2.7946E-01  9.5763E-01  3.4166E+00  8.4342E-01  1.0797E+00  1.3248E-01
             2.3652E+00
 PARAMETER:  1.0645E-01 -1.7834E+00 -1.0974E+00  3.7651E-01 -1.1749E+00  5.6701E-02  1.3287E+00 -7.0293E-02  1.7669E-01 -1.9213E+00
             9.6087E-01
 GRADIENT:  -9.9073E-01 -2.0316E+00  3.1144E+00  1.6984E+00 -2.6959E+00  4.2975E-03 -2.2385E+00  7.7812E-01  1.2652E+00  8.1221E-02
             9.6166E-01

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1388.76571311917        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1951
 NPARAMETR:  1.0068E+00  1.5107E-01  3.0196E-01  1.3174E+00  2.7941E-01  9.5779E-01  3.4236E+00  8.4549E-01  1.0782E+00  1.1408E-01
             2.3618E+00
 PARAMETER:  1.0681E-01 -1.7883E+00 -1.0983E+00  3.7580E-01 -1.1751E+00  5.6908E-02  1.3320E+00 -6.6971E-02  1.7656E-01 -2.0504E+00
             9.5941E-01
 GRADIENT:   2.2519E-02  4.4897E+01 -3.5013E+01  1.0093E+00 -8.5262E-01  7.0437E-02  6.3047E+01  3.1103E-01  1.0257E+00  1.9544E-02
            -1.5609E-02
 NUMSIGDIG:         4.0         2.3         2.4         2.8         3.8         2.7         2.3         1.4         1.5         1.3
                    4.5

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1951
 NO. OF SIG. DIGITS IN FINAL EST.:  1.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3732E-03  3.3144E-02 -1.3111E-02 -1.6678E-02  4.9279E-03
 SE:             2.9159E-02  1.3897E-02  1.9894E-02  2.7027E-02  4.2349E-03
 N:                     100         100         100         100         100

 P VAL.:         9.6244E-01  1.7082E-02  5.0987E-01  5.3717E-01  2.4456E-01

 ETASHRINKSD(%)  2.3122E+00  5.3442E+01  3.3352E+01  9.4570E+00  8.5813E+01
 ETASHRINKVR(%)  4.5710E+00  7.8324E+01  5.5581E+01  1.8020E+01  9.7987E+01
 EBVSHRINKSD(%)  2.2841E+00  6.3586E+01  2.9818E+01  8.1954E+00  8.5813E+01
 EBVSHRINKVR(%)  4.5159E+00  8.6740E+01  5.0744E+01  1.5719E+01  9.7987E+01
 RELATIVEINF(%)  9.3729E+01  6.5559E+00  2.7023E+00  3.4829E+01  1.0253E-01
 EPSSHRINKSD(%)  3.6222E+01
 EPSSHRINKVR(%)  5.9324E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1388.7657131191690     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -653.61488655543087     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1388.766       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.51E-01  3.02E-01  1.32E+00  2.79E-01  9.58E-01  3.43E+00  8.46E-01  1.08E+00  1.16E-01  2.36E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       61.161
Stop Time:
Sun Oct 24 02:19:09 CDT 2021
