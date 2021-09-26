Fri Sep 24 22:31:05 CDT 2021
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
$DATA ../../../../data/int/A3/dat59.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -679.208523268732        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.3077E+01  1.9479E+02  3.3935E+02 -7.0228E+01  2.5482E+02 -1.3278E+01 -1.5479E+02 -3.3359E+02 -7.5660E+01 -2.2867E+02
            -5.6330E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2785.22514925296        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0746E+00  7.8548E-01  7.0891E-01  1.1828E+00  7.1392E-01  1.0498E+00  7.8348E-01  8.4838E-01  8.6346E-01  9.3922E-01
             2.6971E+00
 PARAMETER:  1.7196E-01 -1.4146E-01 -2.4402E-01  2.6791E-01 -2.3699E-01  1.4855E-01 -1.4401E-01 -6.4428E-02 -4.6810E-02  3.7298E-02
             1.0922E+00
 GRADIENT:   3.0355E+01  4.6990E+01 -1.6669E+00  5.0493E+01  1.1091E+01  1.1151E+01 -4.4931E+00  1.4419E+01 -3.4977E+01 -5.2896E-01
             1.8023E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2797.27788730945        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0863E+00  5.6122E-01  4.8405E-01  1.3368E+00  4.9486E-01  1.0076E+00  7.1168E-01  2.4781E-01  9.3228E-01  8.9414E-01
             2.6533E+00
 PARAMETER:  1.8277E-01 -4.7765E-01 -6.2556E-01  3.9028E-01 -6.0349E-01  1.0757E-01 -2.4013E-01 -1.2951E+00  2.9873E-02 -1.1889E-02
             1.0758E+00
 GRADIENT:   5.5439E+01  7.5212E+01 -3.7076E+01  2.3055E+02  1.2946E+01 -7.9904E+00 -1.1003E+01  1.3691E+00 -4.0952E+01  4.9599E+00
             1.8797E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2823.97852720321        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0351E+00  3.3710E-01  2.5443E-01  1.1725E+00  2.8930E-01  1.0674E+00  1.1328E+00  1.4714E-01  1.1050E+00  6.9774E-01
             2.3827E+00
 PARAMETER:  1.3449E-01 -9.8737E-01 -1.2687E+00  2.5912E-01 -1.1403E+00  1.6527E-01  2.2465E-01 -1.8163E+00  1.9987E-01 -2.5992E-01
             9.6825E-01
 GRADIENT:  -3.4811E+01  4.6875E+01 -3.5128E+01  5.8474E+01  2.8575E+01  9.2637E+00 -2.1151E+01 -1.4153E+00 -1.5739E+01 -7.3082E+00
            -9.2566E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2825.66246682418        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0460E+00  2.7797E-01  2.0158E-01  1.1062E+00  2.4056E-01  1.0586E+00  1.2952E+00  1.1658E-01  1.2101E+00  6.9702E-01
             2.3798E+00
 PARAMETER:  1.4497E-01 -1.1802E+00 -1.5016E+00  2.0095E-01 -1.3248E+00  1.5696E-01  3.5863E-01 -2.0492E+00  2.9068E-01 -2.6093E-01
             9.6702E-01
 GRADIENT:  -1.6416E+01  1.9992E+01 -9.4147E-02  4.6348E+01 -1.3112E+01  5.5647E+00 -9.8667E+00 -1.1202E+00 -8.9244E+00 -9.8106E+00
            -5.6520E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2829.82818337586        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      459
 NPARAMETR:  1.0612E+00  3.1903E-01  2.5393E-01  1.1362E+00  2.8508E-01  1.0395E+00  1.2852E+00  1.5487E-01  1.1468E+00  7.1841E-01
             2.4704E+00
 PARAMETER:  1.5939E-01 -1.0425E+00 -1.2707E+00  2.2768E-01 -1.1550E+00  1.3873E-01  3.5091E-01 -1.7652E+00  2.3694E-01 -2.3072E-01
             1.0044E+00
 GRADIENT:   2.1232E+00 -1.4212E+00  2.4577E+00  1.9685E+00  1.7459E+00 -1.2304E-02  8.8435E-01 -1.0365E+00 -2.5096E-01  7.7339E-01
             1.9683E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2842.83697611619        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.0608E+00  2.8483E-01  2.0449E-01  1.1054E+00  2.4378E-01  1.0519E+00  1.2844E+00  1.2831E+00  1.1876E+00  6.2382E-01
             2.3265E+00
 PARAMETER:  1.5898E-01 -1.1559E+00 -1.4872E+00  2.0022E-01 -1.3115E+00  1.5062E-01  3.5027E-01  3.4927E-01  2.7191E-01 -3.7189E-01
             9.4436E-01
 GRADIENT:  -1.9130E+00  3.2390E+01  2.2971E+01  3.2721E+01 -5.7216E+01  4.8364E-01 -4.7764E+00  7.0566E+00 -2.1221E+01  1.8066E+01
             2.4186E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2846.09760719882        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  1.0609E+00  2.6697E-01  1.8599E-01  1.0778E+00  2.2965E-01  1.0520E+00  1.3181E+00  1.3661E+00  1.2472E+00  5.8908E-01
             2.2890E+00
 PARAMETER:  1.5908E-01 -1.2206E+00 -1.5821E+00  1.7488E-01 -1.3712E+00  1.5069E-01  3.7622E-01  4.1199E-01  3.2091E-01 -4.2920E-01
             9.2812E-01
 GRADIENT:  -1.0996E+00  2.6357E+01  2.0631E+01  2.9525E+01 -5.9491E+01  5.9642E-02 -3.5615E+00  3.8378E+00 -1.7532E+01  1.4169E+01
             1.6681E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2846.99192104581        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      923
 NPARAMETR:  1.0608E+00  2.6705E-01  1.8591E-01  1.0514E+00  2.2973E-01  1.0494E+00  1.3180E+00  1.3096E+00  1.3155E+00  5.8914E-01
             2.2885E+00
 PARAMETER:  1.5904E-01 -1.2203E+00 -1.5825E+00  1.5014E-01 -1.3708E+00  1.4825E-01  3.7613E-01  3.6975E-01  3.7420E-01 -4.2909E-01
             9.2788E-01
 GRADIENT:   1.0117E+01  3.3209E+01  3.2712E+01  4.5517E-01 -4.6883E+00  3.4354E-01 -3.4180E+00 -1.4212E-01  5.6538E-01  1.4364E+01
             1.6354E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2846.99326636089        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      993
 NPARAMETR:  1.0608E+00  2.6703E-01  1.8587E-01  1.0506E+00  2.2974E-01  1.0475E+00  1.3180E+00  1.3138E+00  1.3104E+00  5.8914E-01
             2.2884E+00
 PARAMETER:  1.5904E-01 -1.2204E+00 -1.5827E+00  1.4933E-01 -1.3708E+00  1.4641E-01  3.7613E-01  3.7294E-01  3.7031E-01 -4.2910E-01
             9.2784E-01
 GRADIENT:   1.0052E+01  3.3253E+01  3.2439E+01 -4.5653E-01 -4.2861E+00 -3.7670E-01 -3.4136E+00  1.5387E-01 -6.2119E-01  1.4415E+01
             1.6499E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2847.01519167131        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  1.0608E+00  2.6678E-01  1.8531E-01  1.0471E+00  2.2980E-01  1.0410E+00  1.3180E+00  1.3319E+00  1.2917E+00  5.8910E-01
             2.2872E+00
 PARAMETER:  1.5904E-01 -1.2213E+00 -1.5858E+00  1.4601E-01 -1.3706E+00  1.4022E-01  3.7615E-01  3.8662E-01  3.5598E-01 -4.2916E-01
             9.2734E-01
 GRADIENT:   9.8396E+00  3.2799E+01  2.9787E+01 -4.1008E+00  5.1051E-01 -2.8156E+00 -3.4220E+00  1.3349E+00 -5.1953E+00  1.4693E+01
             1.6814E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2847.51439383874        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1144
 NPARAMETR:  1.0608E+00  2.6397E-01  1.7937E-01  1.0422E+00  2.3043E-01  1.0321E+00  1.3184E+00  1.3704E+00  1.2831E+00  5.8869E-01
             2.2748E+00
 PARAMETER:  1.5905E-01 -1.2319E+00 -1.6183E+00  1.4137E-01 -1.3678E+00  1.3155E-01  3.7638E-01  4.1510E-01  3.4926E-01 -4.2986E-01
             9.2191E-01
 GRADIENT:   9.8493E+00  2.3773E+01  6.0350E+00 -6.8849E+00  4.7090E+01 -6.3927E+00 -3.8507E+00  2.1349E+00 -9.3502E+00  1.5662E+01
             1.3278E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2847.65736828822        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1272             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0608E+00  2.6361E-01  1.7862E-01  1.0437E+00  2.3051E-01  1.0372E+00  1.3184E+00  1.3702E+00  1.2951E+00  5.8863E-01
             2.2733E+00
 PARAMETER:  1.5905E-01 -1.2333E+00 -1.6225E+00  1.4274E-01 -1.3674E+00  1.3655E-01  3.7642E-01  4.1499E-01  3.5860E-01 -4.2996E-01
             9.2122E-01
 GRADIENT:   1.0013E+01  2.2381E+01  3.2790E+00 -4.8032E+00  5.3091E+01 -4.4204E+00 -3.9144E+00  1.7190E+00 -6.7536E+00  1.5756E+01
             1.2986E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2847.71927064498        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1439
 NPARAMETR:  1.0609E+00  2.6352E-01  1.7869E-01  1.0506E+00  2.3059E-01  1.0514E+00  1.3185E+00  1.3701E+00  1.2950E+00  5.8857E-01
             2.2738E+00
 PARAMETER:  1.5909E-01 -1.2336E+00 -1.6221E+00  1.4937E-01 -1.3671E+00  1.5014E-01  3.7651E-01  4.1488E-01  3.5851E-01 -4.3006E-01
             9.2145E-01
 GRADIENT:  -2.9056E-01  1.4324E+01 -5.4013E+00  4.8102E-01  7.9046E+00 -1.3834E-01 -4.2232E+00  1.4289E+00 -7.7373E+00  1.5369E+01
             1.1734E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2849.60875979399        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1618             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0609E+00  2.6335E-01  1.7884E-01  1.0501E+00  2.2978E-01  1.0518E+00  1.3484E+00  1.3484E+00  1.3282E+00  4.8595E-01
             2.2599E+00
 PARAMETER:  1.5914E-01 -1.2343E+00 -1.6213E+00  1.4887E-01 -1.3706E+00  1.5049E-01  3.9895E-01  3.9894E-01  3.8379E-01 -6.2165E-01
             9.1530E-01
 GRADIENT:   1.0838E+01  2.3355E+01  8.3515E+00  4.0591E+00  3.6749E+01  8.6514E-01 -2.9065E+00 -6.5779E+00 -1.5948E+00  3.0312E+00
             2.3268E+00

0ITERATION NO.:   72    OBJECTIVE VALUE:  -2849.60875979399        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1684
 NPARAMETR:  1.0609E+00  2.6335E-01  1.7884E-01  1.0501E+00  2.2978E-01  1.0518E+00  1.3484E+00  1.3484E+00  1.3282E+00  4.8595E-01
             2.2599E+00
 PARAMETER:  1.5914E-01 -1.2343E+00 -1.6213E+00  1.4887E-01 -1.3706E+00  1.5049E-01  3.9895E-01  3.9894E-01  3.8379E-01 -6.2165E-01
             9.1530E-01
 GRADIENT:  -1.9957E+04  1.3028E+03 -1.9521E+03  1.0664E+04 -1.1601E+03 -1.0552E+04  7.9575E+03  7.9344E+03  8.2709E+03 -5.1031E+03
            -3.4586E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1684
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.9843E-04  6.9523E-04  4.5539E-03 -2.8660E-03  1.8346E-02
 SE:             2.9596E-02  2.5504E-02  2.4374E-02  2.8267E-02  1.6367E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7848E-01  9.7825E-01  8.5179E-01  9.1924E-01  2.6234E-01

 ETASHRINKSD(%)  8.4993E-01  1.4560E+01  1.8343E+01  5.3020E+00  4.5168E+01
 ETASHRINKVR(%)  1.6926E+00  2.7000E+01  3.3322E+01  1.0323E+01  6.9935E+01
 EBVSHRINKSD(%)  1.1351E+00  1.4524E+01  1.8705E+01  5.1358E+00  4.5730E+01
 EBVSHRINKVR(%)  2.2572E+00  2.6939E+01  3.3912E+01  1.0008E+01  7.0547E+01
 RELATIVEINF(%)  9.7729E+01  2.0227E+01  1.1919E+01  5.2891E+01  3.5267E+00
 EPSSHRINKSD(%)  2.1205E+01
 EPSSHRINKVR(%)  3.7913E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2849.6087597939895     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1195.5194000255788     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    45.30
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2849.609       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  2.63E-01  1.79E-01  1.05E+00  2.30E-01  1.05E+00  1.35E+00  1.35E+00  1.33E+00  4.86E-01  2.26E+00
 


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
+        2.79E+07
 
 TH 2
+        7.39E+02  7.53E+06
 
 TH 3
+       -1.63E+07  5.74E+04  9.39E+06
 
 TH 4
+        3.01E+07 -1.32E+04  1.45E+04  3.25E+07
 
 TH 5
+       -7.56E+02 -1.22E+04 -7.20E+04  1.28E+04  7.96E+06
 
 TH 6
+       -1.32E+03  6.78E+02 -7.58E+02  1.42E+03 -7.01E+02  3.17E+07
 
 TH 7
+        8.74E+06 -3.27E+02  3.99E+02 -7.93E+03  3.63E+02  4.12E+02  2.74E+06
 
 TH 8
+        8.76E+06 -2.27E+04  5.12E+06 -7.93E+03  2.35E+04  4.13E+02 -2.07E+02  2.73E+06
 
 TH 9
+        4.12E+03 -2.15E+03  2.49E+03 -4.45E+03  2.36E+03  4.36E+02 -1.29E+03 -1.29E+03  6.11E+06
 
 TH10
+       -1.56E+07  1.04E+04 -1.17E+04  1.41E+04 -1.06E+04 -7.35E+02  3.83E+02  4.90E+06  2.31E+03  8.69E+06
 
 TH11
+       -1.29E+02  7.22E+03 -9.14E+03  2.06E+03 -7.45E+03 -1.06E+02  6.30E+01  3.59E+03  3.45E+02 -1.64E+03  1.85E+05
 
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
 #CPUT: Total CPU Time in Seconds,       60.832
Stop Time:
Fri Sep 24 22:32:11 CDT 2021
