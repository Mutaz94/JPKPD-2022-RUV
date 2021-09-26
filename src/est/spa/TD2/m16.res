Sat Sep 25 13:21:39 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat16.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1722.60522937183        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.6363E+01 -1.0198E+02 -6.3740E+01 -3.4078E+01  1.3084E+02  4.3451E+01 -7.5416E+00  3.7883E+00  2.9552E+01 -1.9780E+01
             4.0248E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1739.51719281070        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0471E+00  1.0760E+00  1.0253E+00  9.9878E-01  9.6085E-01  8.7727E-01  1.0172E+00  9.9482E-01  8.9621E-01  1.0339E+00
             8.9586E-01
 PARAMETER:  1.4605E-01  1.7322E-01  1.2502E-01  9.8780E-02  6.0065E-02 -3.0942E-02  1.1701E-01  9.4807E-02 -9.5802E-03  1.3329E-01
            -9.9661E-03
 GRADIENT:   1.0755E+02 -5.3953E+00 -7.4079E+00  1.3265E+01  2.1713E+01 -4.2434E-01 -4.7644E+00  1.0571E+00  9.1665E+00 -5.1239E+00
            -3.2750E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1740.33268564110        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0432E+00  1.0394E+00  9.4026E-01  1.0125E+00  9.1269E-01  8.9555E-01  1.2195E+00  8.2320E-01  7.8774E-01  9.8240E-01
             8.9810E-01
 PARAMETER:  1.4226E-01  1.3866E-01  3.8406E-02  1.1244E-01  8.6390E-03 -1.0315E-02  2.9847E-01 -9.4562E-02 -1.3858E-01  8.2242E-02
            -7.4788E-03
 GRADIENT:   9.1618E+01 -1.5882E+00 -8.0325E+00  1.0076E+01  2.4157E+01  7.7446E+00  6.1753E+00  7.9968E-01  8.5297E-01 -6.0994E+00
            -2.1422E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1740.69028106301        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0288E+00  1.0588E+00  8.7032E-01  9.9380E-01  8.8713E-01  8.8633E-01  1.1753E+00  6.7065E-01  8.0087E-01  1.0055E+00
             8.9461E-01
 PARAMETER:  1.2842E-01  1.5713E-01 -3.8900E-02  9.3776E-02 -1.9760E-02 -2.0666E-02  2.6151E-01 -2.9951E-01 -1.2205E-01  1.0553E-01
            -1.1372E-02
 GRADIENT:   4.4307E+01 -3.9861E+00 -7.8616E+00  4.2655E+00  1.4065E+01  3.5527E+00  3.3991E+00  1.1745E+00  1.1096E+00  2.3060E-01
            -2.6189E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1740.69122101156        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0273E+00  1.0619E+00  8.6243E-01  9.9131E-01  8.8382E-01  8.8557E-01  1.1697E+00  6.5039E-01  8.0215E-01  1.0038E+00
             8.9533E-01
 PARAMETER:  1.2694E-01  1.6009E-01 -4.8003E-02  9.1271E-02 -2.3505E-02 -2.1523E-02  2.5673E-01 -3.3018E-01 -1.2046E-01  1.0376E-01
            -1.0563E-02
 GRADIENT:   3.8987E+01 -3.5721E+00 -7.0258E+00  3.7892E+00  1.2462E+01  3.1182E+00  3.0015E+00  1.0682E+00  9.8189E-01  2.5829E-01
            -2.2918E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1740.69133796743        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0269E+00  1.0630E+00  8.5962E-01  9.9046E-01  8.8267E-01  8.8537E-01  1.1680E+00  6.4290E-01  8.0257E-01  1.0031E+00
             8.9553E-01
 PARAMETER:  1.2654E-01  1.6108E-01 -5.1264E-02  9.0414E-02 -2.4803E-02 -2.1746E-02  2.5529E-01 -3.4176E-01 -1.1993E-01  1.0310E-01
            -1.0338E-02
 GRADIENT:   3.7493E+01 -3.4409E+00 -6.7662E+00  3.6464E+00  1.1990E+01  2.9980E+00  2.8874E+00  1.0306E+00  9.4491E-01  2.5366E-01
            -2.2028E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1740.69134425089        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      454
 NPARAMETR:  1.0267E+00  1.0635E+00  8.5831E-01  9.9007E-01  8.8214E-01  8.8529E-01  1.1672E+00  6.3935E-01  8.0277E-01  1.0028E+00
             8.9562E-01
 PARAMETER:  1.2636E-01  1.6153E-01 -5.2793E-02  9.0018E-02 -2.5407E-02 -2.1842E-02  2.5464E-01 -3.4730E-01 -1.1969E-01  1.0279E-01
            -1.0239E-02
 GRADIENT:   3.6837E+01 -3.3817E+00 -6.6493E+00  3.5828E+00  1.1781E+01  2.9456E+00  2.8369E+00  1.0131E+00  9.2848E-01  2.5014E-01
            -2.1641E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1741.36150549184        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      634
 NPARAMETR:  1.0429E+00  1.2267E+00  7.3044E-01  8.9000E-01  8.8284E-01  8.9021E-01  1.0333E+00  3.8898E-01  8.6628E-01  1.0028E+00
             9.0087E-01
 PARAMETER:  1.4205E-01  3.0432E-01 -2.1410E-01 -1.6530E-02 -2.4612E-02 -1.6303E-02  1.3276E-01 -8.4424E-01 -4.3547E-02  1.0279E-01
            -4.3968E-03
 GRADIENT:   9.8062E+00  1.5989E+00 -1.3198E-02  2.4689E+00 -3.7016E+00  4.3827E-01 -2.9871E-02  3.9092E-01  3.2337E-01  1.8275E+00
             2.2619E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1741.61605939357        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      809
 NPARAMETR:  1.0397E+00  1.3851E+00  6.3568E-01  7.8474E-01  9.1483E-01  8.8963E-01  9.4223E-01  1.6319E-01  9.3917E-01  9.9805E-01
             9.0021E-01
 PARAMETER:  1.3898E-01  4.2575E-01 -3.5306E-01 -1.4240E-01  1.0987E-02 -1.6947E-02  4.0495E-02 -1.7128E+00  3.7236E-02  9.8046E-02
            -5.1223E-03
 GRADIENT:  -2.2885E-01 -8.6482E-01 -1.6892E-01 -8.8722E-01  4.2287E-01 -6.4323E-02  5.2766E-02  7.6144E-02 -2.4319E-02  6.6850E-02
             1.2033E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1741.65334514726        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  1.0399E+00  1.3564E+00  6.4481E-01  8.0381E-01  9.0505E-01  8.8979E-01  9.5874E-01  4.6314E-02  9.2639E-01  9.9698E-01
             9.0034E-01
 PARAMETER:  1.3912E-01  4.0483E-01 -3.3880E-01 -1.1840E-01  2.3278E-04 -1.6768E-02  5.7870E-02 -2.9723E+00  2.3541E-02  9.6977E-02
            -4.9857E-03
 GRADIENT:   2.0233E-01  1.2877E-01 -8.1760E-01  1.3795E+00  1.0171E+00  1.3544E-02  2.4773E-01  6.3571E-03  3.5779E-01  3.8421E-01
             1.9137E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1741.65791681368        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  1.0399E+00  1.3616E+00  6.4143E-01  7.9992E-01  9.0534E-01  8.8979E-01  9.5519E-01  1.0000E-02  9.2753E-01  9.9447E-01
             9.0013E-01
 PARAMETER:  1.3908E-01  4.0863E-01 -3.4405E-01 -1.2325E-01  5.5510E-04 -1.6771E-02  5.4155E-02 -4.5478E+00  2.4773E-02  9.4451E-02
            -5.2155E-03
 GRADIENT:   5.0290E-02  4.6818E-02 -6.3286E-02  1.2724E-01  7.3728E-02  6.5877E-03  2.1506E-02  0.0000E+00  2.8547E-02  2.5786E-02
             1.0825E-02

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1741.65792903959        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1220
 NPARAMETR:  1.0398E+00  1.3613E+00  6.4150E-01  7.9999E-01  9.0524E-01  8.8977E-01  9.5523E-01  1.0000E-02  9.2733E-01  9.9430E-01
             9.0012E-01
 PARAMETER:  1.3906E-01  4.0848E-01 -3.4394E-01 -1.2316E-01  4.4806E-04 -1.6787E-02  5.4193E-02 -4.5388E+00  2.4559E-02  9.4285E-02
            -5.2243E-03
 GRADIENT:   3.2183E-03 -2.9714E-04 -3.4399E-03  2.3196E-03  3.0934E-03  7.8012E-04  1.4235E-03  0.0000E+00  2.4210E-03  1.9278E-03
             7.5177E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1220
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9182E-04 -1.4431E-02 -3.7691E-04  9.9909E-03 -2.0125E-02
 SE:             2.9852E-02  2.3386E-02  1.4284E-04  2.2618E-02  2.3912E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9487E-01  5.3718E-01  8.3215E-03  6.5869E-01  3.9999E-01

 ETASHRINKSD(%)  1.0000E-10  2.1655E+01  9.9521E+01  2.4227E+01  1.9892E+01
 ETASHRINKVR(%)  1.0000E-10  3.8620E+01  9.9998E+01  4.2585E+01  3.5827E+01
 EBVSHRINKSD(%)  4.2496E-01  2.1545E+01  9.9589E+01  2.5836E+01  1.7968E+01
 EBVSHRINKVR(%)  8.4811E-01  3.8448E+01  9.9998E+01  4.4996E+01  3.2708E+01
 RELATIVEINF(%)  9.8961E+01  2.8525E+00  1.4937E-04  2.3702E+00  8.4336E+00
 EPSSHRINKSD(%)  4.5143E+01
 EPSSHRINKVR(%)  6.9908E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1741.6579290395932     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1006.5071024758550     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.64
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1741.658       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.36E+00  6.42E-01  8.00E-01  9.05E-01  8.90E-01  9.55E-01  1.00E-02  9.27E-01  9.94E-01  9.00E-01
 


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
+       -6.01E+00  4.20E+02
 
 TH 3
+        1.55E+01  1.55E+02  5.14E+02
 
 TH 4
+       -1.63E+01  4.04E+02 -3.77E+02  1.12E+03
 
 TH 5
+       -3.86E+00 -2.50E+02 -5.58E+02  3.97E+02  9.22E+02
 
 TH 6
+        9.25E-01 -8.93E-01  3.04E+00 -1.49E+00 -1.65E+00  2.47E+02
 
 TH 7
+        1.70E+00  2.05E+01 -1.83E+01 -1.29E+01  1.09E+00  3.35E-01  8.74E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.39E+00 -1.97E+01 -4.09E+01  4.66E+01  8.37E+00 -6.61E-01  2.58E+01  0.00E+00  8.10E+01
 
 TH10
+       -1.55E+00 -1.58E+01 -5.38E+01 -3.55E+00 -5.76E+01 -1.83E+00  1.61E+01  0.00E+00  1.32E+01  8.87E+01
 
 TH11
+       -8.54E+00 -1.44E+01 -2.87E+01  2.48E+00  9.53E+00  2.32E+00  5.51E+00  0.00E+00  7.11E+00  1.75E+01  2.56E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.327
Stop Time:
Sat Sep 25 13:22:00 CDT 2021
