Sat Sep 18 15:08:13 CDT 2021
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
$DATA ../../../../data/spa/D/dat14.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1124.33599626464        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.7907E+02 -2.1506E+02 -1.2025E+02 -2.6262E+02  2.1301E+02 -4.1758E+02 -2.3591E+02  2.3253E-01 -4.2065E+02 -6.7386E+01
            -6.9537E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1452.14707645747        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0040E+00  1.0205E+00  2.5523E+00  8.6709E-01  1.4802E+00  1.4500E+00  2.9995E+00  1.0807E+00  2.2092E+00  1.1796E+00
             1.0713E+00
 PARAMETER:  1.0401E-01  1.2032E-01  1.0370E+00 -4.2611E-02  4.9221E-01  4.7156E-01  1.1984E+00  1.7765E-01  8.9264E-01  2.6517E-01
             1.6889E-01
 GRADIENT:  -6.0547E+01 -4.9404E+01  8.6341E+00 -6.1948E+01  1.1267E+02 -5.2534E+01  6.3114E+01 -1.2214E+00  3.6286E+01 -1.7873E+01
            -2.4006E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1473.24332834840        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.3323E-01  4.3138E-01  3.4604E+00  1.6912E+00  1.2200E+00  1.7544E+00  3.8121E+00  1.3543E+00  1.4284E+00  1.2194E+00
             9.6875E-01
 PARAMETER:  3.0898E-02 -7.4076E-01  1.3414E+00  6.2546E-01  2.9883E-01  6.6210E-01  1.4382E+00  4.0331E-01  4.5657E-01  2.9835E-01
             6.8247E-02
 GRADIENT:  -7.1305E+01  1.0974E+01  2.2639E+01  1.6912E+02 -7.8311E+00  7.9004E+01  2.3053E+00 -9.4420E+00 -9.8154E-01  6.3424E+00
            -7.8716E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1495.08625996103        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0634E+00  1.1839E+00  1.4165E+00  1.0099E+00  1.1103E+00  1.5391E+00  2.5236E+00  1.0041E+00  1.4215E+00  8.9213E-01
             1.1084E+00
 PARAMETER:  1.6143E-01  2.6877E-01  4.4821E-01  1.0981E-01  2.0461E-01  5.3120E-01  1.0257E+00  1.0412E-01  4.5169E-01 -1.4139E-02
             2.0290E-01
 GRADIENT:   2.0631E+01  1.3088E+00  9.9609E+00 -1.5789E+01 -1.1927E+01 -1.5034E+01  1.9815E+00  2.7371E+00  8.3165E-01  1.2380E-01
             4.8833E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1497.41336668026        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0471E+00  1.1460E+00  1.1591E+00  1.0333E+00  1.0296E+00  1.5710E+00  2.5462E+00  5.3307E-01  1.3395E+00  8.2908E-01
             1.0986E+00
 PARAMETER:  1.4600E-01  2.3626E-01  2.4763E-01  1.3272E-01  1.2921E-01  5.5171E-01  1.0346E+00 -5.2910E-01  3.9226E-01 -8.7443E-02
             1.9401E-01
 GRADIENT:   2.7453E+00  9.2102E-01  1.5122E+00 -1.6243E+00 -3.5405E+00 -2.5245E+00 -1.2498E+00  1.0069E+00 -1.8819E+00  1.6103E-01
             2.8600E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1497.72295683823        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0432E+00  1.1310E+00  1.0911E+00  1.0310E+00  1.0062E+00  1.5745E+00  2.5614E+00  2.2217E-01  1.3410E+00  8.0606E-01
             1.0899E+00
 PARAMETER:  1.4234E-01  2.2312E-01  1.8723E-01  1.3058E-01  1.0613E-01  5.5394E-01  1.0405E+00 -1.4043E+00  3.9343E-01 -1.1559E-01
             1.8606E-01
 GRADIENT:  -1.0767E+00 -1.7282E+00 -3.3699E-01 -1.8922E+00  3.8340E-01 -6.6681E-01  5.4043E-01  1.5132E-01 -9.8213E-01 -7.5873E-01
            -5.7969E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1497.81271418731        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.0436E+00  1.1441E+00  1.1019E+00  1.0272E+00  1.0144E+00  1.5768E+00  2.5421E+00  6.3817E-02  1.3586E+00  8.2237E-01
             1.0899E+00
 PARAMETER:  1.4270E-01  2.3466E-01  1.9705E-01  1.2685E-01  1.1426E-01  5.5539E-01  1.0330E+00 -2.6517E+00  4.0648E-01 -9.5559E-02
             1.8608E-01
 GRADIENT:  -4.1906E-01 -4.5048E-01  1.1397E+00 -9.1912E-01 -8.9856E-01  2.3174E-01 -3.4757E-01  9.7438E-03 -4.0714E-03 -1.4376E-01
            -6.8360E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1504.89712381200        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.0990E+00  9.7612E-01  1.1856E+00  1.1310E+00  9.9472E-01  1.7703E+00  3.1062E+00  1.0000E-02  1.2285E+00  8.4597E-01
             1.1086E+00
 PARAMETER:  1.9443E-01  7.5831E-02  2.7022E-01  2.2308E-01  9.4709E-02  6.7115E-01  1.2334E+00 -9.0221E+00  3.0582E-01 -6.7266E-02
             2.0314E-01
 GRADIENT:  -3.8285E+00 -1.0253E+00 -1.9837E+00 -1.3531E+01  9.9235E+00  1.4964E+00 -2.6076E+00  0.0000E+00 -3.2679E+00 -2.2830E+00
             5.4136E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1505.17750458338        NO. OF FUNC. EVALS.: 204
 CUMULATIVE NO. OF FUNC. EVALS.:      792
 NPARAMETR:  1.1011E+00  9.5380E-01  1.1922E+00  1.1467E+00  9.8838E-01  1.7532E+00  3.1606E+00  1.0000E-02  1.2188E+00  8.8034E-01
             1.1058E+00
 PARAMETER:  1.9629E-01  5.2697E-02  2.7581E-01  2.3691E-01  8.8315E-02  6.6143E-01  1.2508E+00 -9.6756E+00  2.9787E-01 -2.7450E-02
             2.0057E-01
 GRADIENT:  -2.4335E+00 -2.6038E-01 -2.3922E+00 -1.2974E+01  3.9568E+00 -2.3996E+00 -2.7005E+00  0.0000E+00 -1.4491E+00  1.6604E+00
             6.1115E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1505.32902617412        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  1.1013E+00  9.5375E-01  1.1953E+00  1.1544E+00  9.8634E-01  1.7620E+00  3.1701E+00  1.0000E-02  1.2276E+00  8.6652E-01
             1.0991E+00
 PARAMETER:  1.9645E-01  5.2648E-02  2.7836E-01  2.4354E-01  8.6247E-02  6.6643E-01  1.2538E+00 -9.6756E+00  3.0508E-01 -4.3268E-02
             1.9449E-01
 GRADIENT:   6.7030E+01  4.1361E+00 -1.7196E+00  1.1261E+01  5.6953E+00  5.6659E+01  3.2255E+01  0.0000E+00  1.2731E+00  4.4489E-02
             2.9969E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1505.33241468875        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  1.1011E+00  9.5370E-01  1.1951E+00  1.1546E+00  9.8627E-01  1.7643E+00  3.1722E+00  1.0000E-02  1.2274E+00  8.6655E-01
             1.0992E+00
 PARAMETER:  1.9633E-01  5.2598E-02  2.7825E-01  2.4375E-01  8.6171E-02  6.6773E-01  1.2544E+00 -9.6756E+00  3.0493E-01 -4.3232E-02
             1.9456E-01
 GRADIENT:  -2.3117E+00  1.4046E+00 -2.1208E+00 -8.3760E+00  5.2271E+00  9.3220E-02 -2.4692E+00  0.0000E+00 -1.0154E+00  3.1615E-03
             2.8851E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1505.34372464431        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1274
 NPARAMETR:  1.1012E+00  9.5346E-01  1.1957E+00  1.1556E+00  9.8617E-01  1.7641E+00  3.1730E+00  1.0000E-02  1.2277E+00  8.6655E-01
             1.0985E+00
 PARAMETER:  1.9638E-01  5.2339E-02  2.7872E-01  2.4459E-01  8.6070E-02  6.6766E-01  1.2547E+00 -9.6756E+00  3.0512E-01 -4.3234E-02
             1.9396E-01
 GRADIENT:  -2.2750E+00  1.5512E+00 -2.1283E+00 -7.8758E+00  5.2756E+00  6.6358E-02 -2.5395E+00  0.0000E+00 -9.8243E-01 -3.8974E-02
             2.6056E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1505.38379928293        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  1.1009E+00  9.5184E-01  1.1958E+00  1.1584E+00  9.8505E-01  1.7639E+00  3.1917E+00  1.0000E-02  1.2274E+00  8.6690E-01
             1.0981E+00
 PARAMETER:  1.9613E-01  5.0639E-02  2.7879E-01  2.4700E-01  8.4935E-02  6.6755E-01  1.2605E+00 -9.6756E+00  3.0492E-01 -4.2830E-02
             1.9354E-01
 GRADIENT:   6.6819E+01  4.9603E+00 -2.1286E+00  1.2945E+01  6.0346E+00  5.7304E+01  3.3421E+01  0.0000E+00  1.6547E+00  4.8636E-02
             2.6599E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1505.39404105520        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1558
 NPARAMETR:  1.1003E+00  9.5134E-01  1.1953E+00  1.1594E+00  9.8467E-01  1.7643E+00  3.2012E+00  1.0000E-02  1.2267E+00  8.6708E-01
             1.0984E+00
 PARAMETER:  1.9563E-01  5.0112E-02  2.7841E-01  2.4788E-01  8.4556E-02  6.6775E-01  1.2635E+00 -9.6756E+00  3.0429E-01 -4.2627E-02
             1.9382E-01
 GRADIENT:   6.6245E+01  5.2199E+00 -2.4487E+00  1.3317E+01  6.4286E+00  5.7398E+01  3.4022E+01  0.0000E+00  1.7555E+00  6.6497E-02
             2.8273E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1505.41346580061        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1705
 NPARAMETR:  1.1007E+00  9.4979E-01  1.1963E+00  1.1605E+00  9.8417E-01  1.7638E+00  3.2025E+00  1.0000E-02  1.2272E+00  8.6705E-01
             1.0975E+00
 PARAMETER:  1.9590E-01  4.8485E-02  2.7927E-01  2.4889E-01  8.4043E-02  6.6746E-01  1.2639E+00 -9.6756E+00  3.0477E-01 -4.2653E-02
             1.9308E-01
 GRADIENT:  -2.5786E+00  2.4829E+00 -2.5296E+00 -6.7869E+00  5.3805E+00 -3.0262E-02 -1.2766E+00  0.0000E+00 -4.8148E-01  2.3510E-02
             2.3334E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1505.48035956756        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1872
 NPARAMETR:  1.1005E+00  9.4984E-01  1.2156E+00  1.1762E+00  9.8412E-01  1.7638E+00  3.2045E+00  1.0000E-02  1.2271E+00  8.7381E-01
             1.0977E+00
 PARAMETER:  1.9581E-01  4.8535E-02  2.9524E-01  2.6228E-01  8.3993E-02  6.6746E-01  1.2646E+00 -9.6756E+00  3.0462E-01 -3.4887E-02
             1.9317E-01
 GRADIENT:  -2.5763E+00  5.5955E+00 -1.8881E-01  4.4157E-02  9.3277E-01 -1.3121E-02 -3.2678E+00  0.0000E+00 -6.5218E-01  1.8003E-01
             1.8468E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1505.49498035494        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     2017
 NPARAMETR:  1.1007E+00  9.4782E-01  1.2161E+00  1.1762E+00  9.8409E-01  1.7637E+00  3.2064E+00  1.0000E-02  1.2273E+00  8.7240E-01
             1.0974E+00
 PARAMETER:  1.9593E-01  4.6408E-02  2.9567E-01  2.6232E-01  8.3957E-02  6.6742E-01  1.2652E+00 -9.6756E+00  3.0485E-01 -3.6512E-02
             1.9298E-01
 GRADIENT:   6.6785E+01  8.0458E+00  1.8266E-01  2.1423E+01  1.7801E+00  5.7288E+01  3.1343E+01  0.0000E+00  1.7512E+00  4.7076E-02
             1.8439E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1505.50195413980        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2191
 NPARAMETR:  1.0999E+00  9.4815E-01  1.2149E+00  1.1773E+00  9.8374E-01  1.7641E+00  3.2207E+00  1.0000E-02  1.2260E+00  8.7244E-01
             1.0982E+00
 PARAMETER:  1.9524E-01  4.6759E-02  2.9463E-01  2.6325E-01  8.3606E-02  6.6763E-01  1.2696E+00 -9.6756E+00  3.0378E-01 -3.6465E-02
             1.9366E-01
 GRADIENT:  -2.9715E+00  5.7715E+00 -7.5096E-01  1.0483E-01  2.0533E+00  4.8656E-02 -2.3482E+00  0.0000E+00 -5.0696E-01  4.2018E-03
             2.0564E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1505.52377018750        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     2347
 NPARAMETR:  1.0998E+00  9.4565E-01  1.2143E+00  1.1780E+00  9.8342E-01  1.7643E+00  3.2318E+00  1.0000E-02  1.2256E+00  8.7208E-01
             1.0983E+00
 PARAMETER:  1.9510E-01  4.4117E-02  2.9414E-01  2.6379E-01  8.3286E-02  6.6775E-01  1.2731E+00 -9.6756E+00  3.0339E-01 -3.6872E-02
             1.9373E-01
 GRADIENT:   6.5788E+01  8.3480E+00 -8.7367E-01  2.1796E+01  3.4587E+00  5.7405E+01  3.3074E+01  0.0000E+00  2.0375E+00  6.6724E-03
             2.2941E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1505.53595394467        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     2493
 NPARAMETR:  1.0999E+00  9.4389E-01  1.2144E+00  1.1780E+00  9.8334E-01  1.7641E+00  3.2328E+00  1.0000E-02  1.2256E+00  8.7235E-01
             1.0980E+00
 PARAMETER:  1.9525E-01  4.2253E-02  2.9428E-01  2.6381E-01  8.3196E-02  6.6763E-01  1.2733E+00 -9.6756E+00  3.0347E-01 -3.6570E-02
             1.9353E-01
 GRADIENT:   6.5982E+01  8.0943E+00 -9.1402E-01  2.1454E+01  3.5311E+00  5.7354E+01  3.3101E+01  0.0000E+00  2.1048E+00  4.1589E-02
             2.2340E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1505.54422769713        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2671            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0997E+00  9.4307E-01  1.2139E+00  1.1786E+00  9.8307E-01  1.7641E+00  3.2416E+00  1.0000E-02  1.2250E+00  8.7252E-01
             1.0983E+00
 PARAMETER:  1.9501E-01  4.1381E-02  2.9382E-01  2.6435E-01  8.2921E-02  6.6765E-01  1.2761E+00 -9.6756E+00  3.0290E-01 -3.6374E-02
             1.9375E-01
 GRADIENT:   6.5685E+01  8.2023E+00 -1.2427E+00  2.1554E+01  3.9983E+00  5.7351E+01  3.3682E+01  0.0000E+00  2.2032E+00  5.9803E-02
             2.3834E+00

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1505.54907845220        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     2799
 NPARAMETR:  1.0995E+00  9.4247E-01  1.2136E+00  1.1790E+00  9.8290E-01  1.7649E+00  3.2460E+00  1.0000E-02  1.2246E+00  8.7260E-01
             1.0984E+00
 PARAMETER:  1.9503E-01  4.0644E-02  2.9391E-01  2.6440E-01  8.2855E-02  6.6758E-01  1.2762E+00 -9.6756E+00  3.0289E-01 -3.6377E-02
             1.9365E-01
 GRADIENT:   5.7635E+04 -1.1240E+05  3.8235E+04 -2.1261E+04  5.6205E+04 -9.7731E-02 -4.4060E+03  0.0000E+00  3.7104E+04 -1.1242E+05
            -5.8054E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2799
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7056E-03  1.5137E-02 -3.9612E-04 -2.6602E-02 -5.4045E-03
 SE:             2.9935E-02  2.4224E-02  1.3003E-04  2.2896E-02  2.0281E-02
 N:                     100         100         100         100         100

 P VAL.:         9.2798E-01  5.3204E-01  2.3168E-03  2.4530E-01  7.8986E-01

 ETASHRINKSD(%)  1.0000E-10  1.8847E+01  9.9564E+01  2.3296E+01  3.2058E+01
 ETASHRINKVR(%)  1.0000E-10  3.4141E+01  9.9998E+01  4.1165E+01  5.3838E+01
 EBVSHRINKSD(%)  1.7026E-01  1.5989E+01  9.9604E+01  2.6394E+01  3.0089E+01
 EBVSHRINKVR(%)  3.4022E-01  2.9422E+01  9.9998E+01  4.5822E+01  5.1125E+01
 RELATIVEINF(%)  9.9571E+01  1.9592E+01  2.9356E-04  1.2718E+01  9.3896E+00
 EPSSHRINKSD(%)  4.2579E+01
 EPSSHRINKVR(%)  6.7028E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1505.5490784521951     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -770.39825188845691     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1505.549       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.10E+00  9.42E-01  1.21E+00  1.18E+00  9.83E-01  1.76E+00  3.24E+00  1.00E-02  1.22E+00  8.73E-01  1.10E+00
 


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
+        6.11E+07
 
 TH 2
+        6.83E+03  3.16E+08
 
 TH 3
+       -2.05E+04  4.68E+04  2.21E+07
 
 TH 4
+       -1.11E+04  2.54E+04 -2.53E+07  2.89E+07
 
 TH 5
+       -2.23E+04 -3.03E+08 -4.51E+04 -2.41E+04  2.91E+08
 
 TH 6
+       -1.06E+01  2.74E+01 -6.54E+00  8.02E+00 -2.68E+01  6.39E+01
 
 TH 7
+       -5.03E+02  1.15E+03  1.06E+03  2.18E+06 -6.91E+06  5.06E-01  1.64E+05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.53E+07  3.94E+03 -1.19E+04 -6.42E+03 -1.29E+04 -7.23E+00 -2.86E+02  0.00E+00  2.04E+07
 
 TH10
+       -2.44E+04  5.55E+04 -9.03E+07  1.03E+08 -5.32E+04  2.93E+01  1.26E+03  0.00E+00 -1.41E+04  3.69E+08
 
 TH11
+       -1.07E+04  2.44E+04 -3.70E+07  4.24E+07 -2.34E+04  1.33E+01  3.19E+06  0.00E+00 -6.20E+03  1.51E+08  6.22E+07
 
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
 #CPUT: Total CPU Time in Seconds,       48.734
Stop Time:
Sat Sep 18 15:09:04 CDT 2021
