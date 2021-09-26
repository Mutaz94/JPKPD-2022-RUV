Sat Sep 25 12:45:11 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1742.96732009260        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -6.1840E+01 -3.9802E+01 -1.0687E+01 -6.9917E+01 -4.7095E+01  3.4177E+00  1.5392E-01  1.6976E+01  2.1479E+00  2.0434E+01
             3.5903E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1752.86061218768        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0533E+00  1.0506E+00  1.1957E+00  1.0032E+00  1.1524E+00  9.6964E-01  1.0091E+00  8.7080E-01  9.8015E-01  9.3710E-01
             8.9804E-01
 PARAMETER:  1.5195E-01  1.4939E-01  2.7872E-01  1.0319E-01  2.4184E-01  6.9172E-02  1.0904E-01 -3.8347E-02  7.9949E-02  3.5030E-02
            -7.5411E-03
 GRADIENT:   9.8228E+01 -1.0410E+01  1.0370E+01 -2.8681E+01  1.2541E+01 -3.5488E+00 -1.6754E+00 -2.9442E-01 -6.3933E+00 -1.8491E+01
            -1.9074E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1754.21303606995        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0421E+00  9.9239E-01  1.0037E+00  1.0384E+00  1.0299E+00  1.0031E+00  9.8308E-01  4.4917E-01  9.4576E-01  9.7175E-01
             8.9691E-01
 PARAMETER:  1.4128E-01  9.2358E-02  1.0367E-01  1.3763E-01  1.2942E-01  1.0305E-01  8.2935E-02 -7.0035E-01  4.4237E-02  7.1339E-02
            -8.7978E-03
 GRADIENT:   6.5151E+01 -1.1223E+01  5.3337E-01 -1.2143E+01  1.5807E-01  9.6684E+00 -9.1040E+00 -2.0854E-01 -8.6785E+00 -6.0987E-01
            -1.5077E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1754.71593654504        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0295E+00  9.9145E-01  9.6934E-01  1.0428E+00  1.0166E+00  9.9378E-01  1.1387E+00  4.3915E-01  9.0651E-01  9.5512E-01
             9.0614E-01
 PARAMETER:  1.2907E-01  9.1416E-02  6.8864E-02  1.4193E-01  1.1647E-01  9.3756E-02  2.2992E-01 -7.2291E-01  1.8509E-03  5.4082E-02
             1.4370E-03
 GRADIENT:   2.7958E+01 -5.8069E+00 -6.4662E+00  2.0235E-01  5.6659E+00  5.1506E+00 -2.4197E+00  8.5518E-01 -4.9962E+00  2.1024E+00
            -7.8073E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1755.59055516057        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      390
 NPARAMETR:  1.0527E+00  8.7519E-01  1.0588E+00  1.1342E+00  1.0011E+00  9.9278E-01  1.2555E+00  4.7933E-01  8.9021E-01  9.5127E-01
             9.2519E-01
 PARAMETER:  1.5140E-01 -3.3317E-02  1.5710E-01  2.2592E-01  1.0114E-01  9.2751E-02  3.2753E-01 -6.3537E-01 -1.6294E-02  5.0039E-02
             2.2239E-02
 GRADIENT:   9.1469E+00  6.5984E+00  3.0071E+00  8.8422E+00 -3.6211E+00  6.7512E-02 -3.7440E-01 -2.0360E-01 -2.7362E-01 -1.1559E+00
            -9.3260E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1755.72727182687        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      566
 NPARAMETR:  1.0480E+00  7.3834E-01  1.0645E+00  1.2083E+00  9.5034E-01  9.9095E-01  1.4400E+00  4.2636E-01  8.4181E-01  9.4272E-01
             9.2479E-01
 PARAMETER:  1.4688E-01 -2.0336E-01  1.6248E-01  2.8918E-01  4.9063E-02  9.0907E-02  4.6467E-01 -7.5246E-01 -7.2198E-02  4.1011E-02
             2.1806E-02
 GRADIENT:   2.0237E+00 -5.7544E-01  1.6562E+00 -2.6373E+00 -1.1446E+00 -1.1999E-02 -2.4987E-01 -4.1980E-01  3.0647E-01 -1.8644E-01
            -1.6252E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1755.80392932918        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      742
 NPARAMETR:  1.0448E+00  6.4097E-01  1.1616E+00  1.2746E+00  9.6225E-01  9.8862E-01  1.5596E+00  5.8044E-01  8.2231E-01  9.6377E-01
             9.2447E-01
 PARAMETER:  1.4386E-01 -3.4477E-01  2.4982E-01  3.4260E-01  6.1520E-02  8.8555E-02  5.4443E-01 -4.4397E-01 -9.5633E-02  6.3094E-02
             2.1464E-02
 GRADIENT:  -6.9424E-01  1.4840E-01  3.3680E-02 -4.7283E-02  2.6846E-01 -1.8334E-01  1.2542E-01 -4.5995E-02  1.4439E-01 -9.6241E-02
             1.3417E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1755.81375104574        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      918
 NPARAMETR:  1.0430E+00  5.3429E-01  1.2518E+00  1.3462E+00  9.6471E-01  9.8699E-01  1.6695E+00  6.8474E-01  8.0988E-01  9.8348E-01
             9.2443E-01
 PARAMETER:  1.4208E-01 -5.2681E-01  3.2460E-01  3.9732E-01  6.4073E-02  8.6903E-02  6.1253E-01 -2.7872E-01 -1.1087E-01  8.3340E-02
             2.1423E-02
 GRADIENT:  -3.5658E-01  1.2769E+00  7.0736E-01  3.1037E+00 -1.7094E+00 -1.0989E-01  1.9380E-01  3.9699E-02 -4.4748E-02  1.4278E-01
             8.8523E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1755.92548187904        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1097
 NPARAMETR:  1.0385E+00  3.3871E-01  1.3723E+00  1.4733E+00  9.5277E-01  9.8280E-01  1.9785E+00  8.0480E-01  7.8317E-01  1.0049E+00
             9.2573E-01
 PARAMETER:  1.3776E-01 -9.8261E-01  4.1649E-01  4.8751E-01  5.1620E-02  8.2646E-02  7.8233E-01 -1.1716E-01 -1.4440E-01  1.0493E-01
             2.2823E-02
 GRADIENT:  -2.1820E+00  2.1888E+00  1.2828E+00  1.0186E+01 -3.3648E+00 -5.6619E-01  3.2436E-01 -1.7715E-01 -2.2437E-01  5.2691E-01
             8.1341E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1756.06142260441        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  1.0365E+00  2.1233E-01  1.4752E+00  1.5488E+00  9.5805E-01  9.8116E-01  2.2001E+00  9.2062E-01  7.6496E-01  1.0140E+00
             9.2505E-01
 PARAMETER:  1.3584E-01 -1.4496E+00  4.8882E-01  5.3751E-01  5.7143E-02  8.0976E-02  8.8851E-01  1.7290E-02 -1.6794E-01  1.1387E-01
             2.2091E-02
 GRADIENT:  -1.1161E+00  1.2948E-01  2.8941E-01 -1.2246E+00  1.4068E-01 -4.2013E-01  4.2328E-02  3.1796E-04 -4.5961E-01 -3.2600E-01
             6.4683E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1756.12699420095        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1447
 NPARAMETR:  1.0348E+00  1.1804E-01  1.5098E+00  1.6081E+00  9.4506E-01  9.8063E-01  2.3571E+00  9.6199E-01  7.5069E-01  1.0216E+00
             9.2260E-01
 PARAMETER:  1.3418E-01 -2.0367E+00  5.1198E-01  5.7506E-01  4.3497E-02  8.0437E-02  9.5745E-01  6.1250E-02 -1.8677E-01  1.2136E-01
             1.9440E-02
 GRADIENT:  -8.4809E-01  2.3737E-01 -1.1652E+00  3.1349E+00 -1.5850E-01 -6.3989E-02  6.4612E-02  4.3688E-01 -3.7624E-01  6.8620E-01
             2.2414E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1756.16640412036        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1622
 NPARAMETR:  1.0342E+00  5.7211E-02  1.5671E+00  1.6451E+00  9.4880E-01  9.7995E-01  2.3382E+00  1.0110E+00  7.3882E-01  1.0215E+00
             9.2231E-01
 PARAMETER:  1.3362E-01 -2.7610E+00  5.4920E-01  5.9783E-01  4.7443E-02  7.9743E-02  9.4937E-01  1.1095E-01 -2.0270E-01  1.2128E-01
             1.9129E-02
 GRADIENT:   4.8169E-01 -1.6605E-02  2.8017E-01 -2.0956E+00  4.4154E-01  6.2864E-02  4.6529E-02 -5.0079E-02  8.7515E-03 -3.0921E-01
            -2.2871E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1756.17458445362        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1797
 NPARAMETR:  1.0334E+00  3.4023E-02  1.5541E+00  1.6584E+00  9.3818E-01  9.7935E-01  2.2394E+00  9.9784E-01  7.3407E-01  1.0203E+00
             9.2224E-01
 PARAMETER:  1.3289E-01 -3.2807E+00  5.4089E-01  6.0583E-01  3.6192E-02  7.9132E-02  9.0621E-01  9.7836E-02 -2.0915E-01  1.2007E-01
             1.9052E-02
 GRADIENT:  -5.5715E-02 -1.0361E-02 -1.7482E-01 -9.8402E-01  2.3095E-01 -1.1081E-02  1.9737E-02  1.4199E-02  5.9778E-02  6.0391E-02
             3.3631E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1756.17602324947        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1961
 NPARAMETR:  1.0333E+00  3.0479E-02  1.5586E+00  1.6609E+00  9.3864E-01  9.7929E-01  1.9545E+00  9.9779E-01  7.3291E-01  1.0194E+00
             9.2221E-01
 PARAMETER:  1.3279E-01 -3.3907E+00  5.4382E-01  6.0733E-01  3.6672E-02  7.9077E-02  7.7016E-01  9.7792E-02 -2.1073E-01  1.1925E-01
             1.9015E-02
 GRADIENT:  -1.3668E-01 -7.3588E-03  3.6319E-01 -6.3645E-01  1.5316E-01 -6.2397E-03  1.3367E-02 -2.5290E-01 -8.9760E-02 -2.4088E-01
            -1.5062E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1756.18543726331        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2137
 NPARAMETR:  1.0334E+00  3.2786E-02  1.5592E+00  1.6598E+00  9.3928E-01  9.7930E-01  3.3056E-01  1.0036E+00  7.3408E-01  1.0211E+00
             9.2227E-01
 PARAMETER:  1.3286E-01 -3.3178E+00  5.4419E-01  6.0668E-01  3.7358E-02  7.9086E-02 -1.0070E+00  1.0358E-01 -2.0913E-01  1.2084E-01
             1.9088E-02
 GRADIENT:  -1.3349E-01 -1.2206E-02  6.8887E-02  6.9017E-02 -3.6684E-01 -3.3308E-02  5.6957E-04  5.4101E-02  2.0063E-02  1.1984E-01
             2.7766E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1756.20417829394        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2313
 NPARAMETR:  1.0343E+00  6.8656E-02  1.5582E+00  1.6389E+00  9.4818E-01  9.7828E-01  1.0000E-02  1.0005E+00  7.4258E-01  1.0262E+00
             9.2261E-01
 PARAMETER:  1.3375E-01 -2.5787E+00  5.4352E-01  5.9400E-01  4.6785E-02  7.8037E-02 -1.0398E+01  1.0052E-01 -1.9762E-01  1.2584E-01
             1.9455E-02
 GRADIENT:   2.4416E-01  4.6739E-02  6.2060E-01  1.3401E+00 -1.2682E+00 -6.9047E-01  0.0000E+00 -2.9036E-02 -6.6842E-01  3.3081E-01
            -8.7373E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1756.22315416263        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2489
 NPARAMETR:  1.0354E+00  1.2431E-01  1.5503E+00  1.6075E+00  9.6085E-01  9.8087E-01  1.0000E-02  9.9449E-01  7.5877E-01  1.0287E+00
             9.2330E-01
 PARAMETER:  1.3482E-01 -1.9850E+00  5.3844E-01  5.7470E-01  6.0058E-02  8.0684E-02 -2.0904E+01  9.4471E-02 -1.7606E-01  1.2830E-01
             2.0197E-02
 GRADIENT:   1.4394E-01  3.7114E-01  4.1495E-01  6.9579E+00 -1.2700E+00 -4.6918E-02  0.0000E+00 -5.9954E-02 -6.6104E-01 -2.2224E-02
            -1.3826E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1756.24306217120        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2667
 NPARAMETR:  1.0366E+00  1.7923E-01  1.5383E+00  1.5724E+00  9.7315E-01  9.8127E-01  1.0000E-02  9.8306E-01  7.7573E-01  1.0351E+00
             9.2393E-01
 PARAMETER:  1.3591E-01 -1.6191E+00  5.3065E-01  5.5261E-01  7.2783E-02  8.1093E-02 -3.0500E+01  8.2914E-02 -1.5395E-01  1.3446E-01
             2.0878E-02
 GRADIENT:   2.9023E-01  2.8481E-01 -9.5054E-02  3.5449E+00 -2.9156E-01 -2.4732E-01  0.0000E+00 -2.0553E-02 -5.8060E-01  1.1448E-01
            -5.2679E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1756.24911371095        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2852             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0371E+00  2.1240E-01  1.5317E+00  1.5504E+00  9.8057E-01  9.8243E-01  1.0000E-02  9.7767E-01  7.8758E-01  1.0373E+00
             9.2433E-01
 PARAMETER:  1.3643E-01 -1.4493E+00  5.2637E-01  5.3849E-01  8.0377E-02  8.2272E-02 -3.6157E+01  7.7417E-02 -1.3879E-01  1.3661E-01
             2.1311E-02
 GRADIENT:   7.2506E+01  3.2217E+00  8.9442E-01  1.1262E+02  6.6385E-01  5.4281E+00  0.0000E+00  1.8081E-02  1.6642E+00  1.1755E-01
             8.8075E-02

0ITERATION NO.:   92    OBJECTIVE VALUE:  -1756.24911371095        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:     2912
 NPARAMETR:  1.0371E+00  2.1240E-01  1.5317E+00  1.5504E+00  9.8057E-01  9.8243E-01  1.0000E-02  9.7767E-01  7.8758E-01  1.0373E+00
             9.2433E-01
 PARAMETER:  1.3643E-01 -1.4493E+00  5.2637E-01  5.3849E-01  8.0377E-02  8.2272E-02 -3.6157E+01  7.7417E-02 -1.3879E-01  1.3661E-01
             2.1311E-02
 GRADIENT:   1.2251E-01 -1.0927E-02 -5.6373E-02 -1.9718E-01 -3.0238E-03  3.9251E-03  0.0000E+00  1.2929E-03  2.2159E-02 -1.3882E-03
             8.0756E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2912
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9461E-05 -7.8809E-05 -2.6459E-02 -5.5030E-03 -3.5784E-02
 SE:             2.9853E-02  3.7696E-05  1.4922E-02  2.9367E-02  2.2201E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9734E-01  3.6557E-02  7.6197E-02  8.5136E-01  1.0700E-01

 ETASHRINKSD(%)  1.0000E-10  9.9874E+01  5.0010E+01  1.6176E+00  2.5625E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  7.5010E+01  3.2090E+00  4.4683E+01
 EBVSHRINKSD(%)  3.7499E-01  9.9881E+01  5.3009E+01  1.9607E+00  2.2675E+01
 EBVSHRINKVR(%)  7.4857E-01  1.0000E+02  7.7918E+01  3.8830E+00  4.0209E+01
 RELATIVEINF(%)  9.7163E+01  7.6579E-06  5.0819E+00  6.4048E+00  8.4251E+00
 EPSSHRINKSD(%)  4.4149E+01
 EPSSHRINKVR(%)  6.8807E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1756.2491137109455     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1021.0982871472073     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.54
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     5.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1756.249       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  2.12E-01  1.53E+00  1.55E+00  9.81E-01  9.82E-01  1.00E-02  9.78E-01  7.88E-01  1.04E+00  9.24E-01
 


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
+        1.12E+03
 
 TH 2
+       -1.83E+01  3.94E+02
 
 TH 3
+        7.82E+00  4.03E+01  9.51E+01
 
 TH 4
+        5.47E+01  4.72E+02 -3.81E+01  7.44E+02
 
 TH 5
+       -2.98E+01 -2.12E+02 -2.42E+02 -3.55E+01  6.74E+02
 
 TH 6
+        2.61E+01 -2.95E+01  1.32E+01 -4.09E+01 -1.53E+01  1.87E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        8.65E+00 -1.88E+00 -4.98E+00  1.23E+00  6.37E+00  5.43E+00  0.00E+00  1.98E+00
 
 TH 9
+        1.12E+02 -1.50E+02  1.79E+01 -3.84E+01  1.52E-02  1.39E+01  0.00E+00 -4.29E-01  3.88E+02
 
 TH10
+        4.54E+00  1.93E+01  1.54E+01  6.59E+00 -5.98E+01 -7.81E+00  0.00E+00  2.70E+00  1.06E+00  1.32E+01
 
 TH11
+        4.14E+01  2.68E+00 -6.31E+01  5.71E+01  5.40E+01  1.86E+01  0.00E+00  2.72E+01  1.63E+01  4.83E+01  4.05E+02
 
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
+        1.07E+03
 
 TH 2
+       -2.69E+01  3.65E+02
 
 TH 3
+       -8.25E-01  5.01E+01  1.16E+02
 
 TH 4
+       -7.77E+00  4.45E+02 -2.46E+01  7.17E+02
 
 TH 5
+        2.96E+00 -2.12E+02 -2.21E+02 -4.14E+01  6.65E+02
 
 TH 6
+       -8.13E+00 -5.05E+00  4.03E-01 -1.60E+00 -5.89E+00  2.11E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        5.29E+00 -4.22E+00 -2.72E+01 -3.39E+00 -5.35E+00  1.13E+01  0.00E+00  3.24E+01
 
 TH 9
+        1.42E+00 -9.87E+01  8.66E+00 -4.80E-02  3.80E+00 -1.60E+00  0.00E+00  9.67E-01  3.06E+02
 
 TH10
+       -3.60E+00  8.20E+00 -1.69E+00 -1.46E-01 -6.92E+01  1.39E-01  0.00E+00  1.24E+01  1.02E+00  7.73E+01
 
 TH11
+       -6.35E+00 -1.26E+01 -1.77E+01 -1.05E+01  1.84E+00  1.01E+00  0.00E+00  1.05E+01  4.91E+00  2.21E+01  2.59E+02
 
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
+        1.06E+03
 
 TH 2
+       -5.49E+01  3.43E+02
 
 TH 3
+       -1.72E+01  5.43E+01  1.14E+02
 
 TH 4
+       -7.60E+01  4.37E+02 -1.69E+01  7.17E+02
 
 TH 5
+        3.67E+01 -2.07E+02 -2.14E+02 -3.81E+01  6.61E+02
 
 TH 6
+       -3.85E+01  2.79E+01  7.70E+00  4.04E+01  1.20E+01  2.28E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.92E+00 -6.00E+00 -2.68E+01 -4.26E+00  3.31E+00 -9.71E+00  0.00E+00  1.98E+01
 
 TH 9
+       -8.55E+01 -5.53E+01  9.89E+00  3.27E+01  7.92E+00 -1.29E+01  0.00E+00 -5.22E+00  2.51E+02
 
 TH10
+        1.52E+01  5.00E-01 -4.79E+00 -1.69E+01 -8.47E+01 -2.29E+01  0.00E+00  1.50E+01 -4.25E+00  8.09E+01
 
 TH11
+       -3.51E+01 -2.41E+01  6.50E+00 -5.13E+01 -3.32E+01 -5.65E+00  0.00E+00  5.09E+00  1.19E+00  1.65E+01  1.71E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       40.141
Stop Time:
Sat Sep 25 12:45:53 CDT 2021
