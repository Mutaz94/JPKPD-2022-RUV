Wed Sep 29 04:54:31 CDT 2021
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
$DATA ../../../../data/int/SL3/dat92.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      979
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

 TOT. NO. OF OBS RECS:      879
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -617.524636050453        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5251E+02  2.6273E+01  1.5670E+02  1.0991E+02  1.8513E+02  6.0327E+01 -1.2585E+02 -2.3632E+02 -1.2933E+02 -1.7902E+01
            -5.8234E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2680.29839726952        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0454E+00  1.3325E+00  8.8353E-01  9.4829E-01  1.0423E+00  8.4030E-01  1.3805E+00  1.0142E+00  1.2203E+00  9.4235E-01
             2.5990E+00
 PARAMETER:  1.4437E-01  3.8703E-01 -2.3830E-02  4.6901E-02  1.4142E-01 -7.3994E-02  4.2247E-01  1.1406E-01  2.9908E-01  4.0620E-02
             1.0551E+00
 GRADIENT:   2.1127E+02  1.2707E+02 -1.2936E+01  8.2003E+01 -1.3640E+01 -4.0368E+01  3.2415E+01  2.3868E+00  1.1692E+01  3.1472E+00
            -1.9754E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2684.05751167307        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      222
 NPARAMETR:  1.0546E+00  1.5002E+00  9.9270E-01  8.6849E-01  1.1965E+00  9.4264E-01  1.2402E+00  9.5319E-01  1.3128E+00  9.7209E-01
             2.6158E+00
 PARAMETER:  1.5317E-01  5.0558E-01  9.2669E-02 -4.1002E-02  2.7942E-01  4.0930E-02  3.1526E-01  5.2058E-02  3.7220E-01  7.1688E-02
             1.0616E+00
 GRADIENT:   1.1633E+02  8.7411E+01  7.7237E-01  7.7516E+01 -9.2494E+00 -7.8136E-01  1.5416E+01 -2.0507E+00  7.0382E+00 -1.2995E+01
            -2.6800E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2697.60217500281        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  1.0076E+00  1.6157E+00  1.2792E+00  7.1320E-01  1.4543E+00  9.3058E-01  9.5100E-01  1.7293E+00  1.2879E+00  1.2717E+00
             2.6250E+00
 PARAMETER:  1.0758E-01  5.7978E-01  3.4621E-01 -2.3799E-01  4.7451E-01  2.8057E-02  4.9761E-02  6.4774E-01  3.5300E-01  3.4034E-01
             1.0651E+00
 GRADIENT:   6.6927E+00 -1.1119E+00 -4.5023E-01  2.3236E+01  1.3735E+01  5.8628E-01 -7.1243E+00 -1.7575E+00 -5.1200E+00 -3.5188E-01
             2.3785E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2702.67596654697        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      574
 NPARAMETR:  1.0060E+00  2.0878E+00  1.1420E+00  4.0054E-01  1.7580E+00  9.2409E-01  8.5089E-01  2.1568E+00  1.8166E+00  1.5364E+00
             2.6194E+00
 PARAMETER:  1.0594E-01  8.3609E-01  2.3282E-01 -8.1493E-01  6.6416E-01  2.1057E-02 -6.1476E-02  8.6862E-01  6.9697E-01  5.2941E-01
             1.0629E+00
 GRADIENT:   3.5861E+00  1.2660E+01  2.5524E+00  5.9626E+00 -6.8295E+00 -1.9307E+00 -2.2527E+00  1.0476E-02 -3.0871E-01  2.0484E+00
             4.2448E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2706.14256403000        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      752
 NPARAMETR:  1.0058E+00  2.5057E+00  4.2900E-01  1.3091E-01  1.9665E+00  9.4170E-01  7.8566E-01  6.8935E-01  3.5042E+00  1.6601E+00
             2.5976E+00
 PARAMETER:  1.0576E-01  1.0186E+00 -7.4629E-01 -1.9332E+00  7.7626E-01  3.9931E-02 -1.4123E-01 -2.7201E-01  1.3540E+00  6.0690E-01
             1.0546E+00
 GRADIENT:   3.7505E+00  3.3216E+01 -3.7140E-01  7.9707E+00 -5.1870E+00  4.4799E+00 -3.6894E+00  4.5834E-01  8.9670E+00 -3.6771E+00
            -1.1917E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2707.79014993939        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      931
 NPARAMETR:  1.0051E+00  2.5335E+00  2.9258E-01  9.0128E-02  2.0070E+00  9.2780E-01  7.9469E-01  4.4174E-01  3.7116E+00  1.7043E+00
             2.6026E+00
 PARAMETER:  1.0508E-01  1.0296E+00 -1.1290E+00 -2.3065E+00  7.9663E-01  2.5064E-02 -1.2980E-01 -7.1704E-01  1.4115E+00  6.3313E-01
             1.0565E+00
 GRADIENT:   2.3802E+00 -5.0935E+00 -2.6322E+00  4.2337E+00 -1.2825E+00 -5.5345E-01 -1.1003E+00  2.3575E-01  4.7761E+00 -3.0962E+00
            -2.3931E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2707.89417604434        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1109
 NPARAMETR:  1.0031E+00  2.5541E+00  2.7291E-01  7.4210E-02  2.0331E+00  9.2932E-01  7.9494E-01  3.7553E-01  4.1394E+00  1.7319E+00
             2.6004E+00
 PARAMETER:  1.0311E-01  1.0377E+00 -1.1986E+00 -2.5009E+00  8.0954E-01  2.6702E-02 -1.2949E-01 -8.7940E-01  1.5206E+00  6.4923E-01
             1.0557E+00
 GRADIENT:  -2.7299E+00 -1.6384E+00 -1.0743E+00  6.4895E-01  3.7001E+00 -5.6473E-02  1.1744E+00  1.5097E-01 -5.9157E-02  1.3874E+00
             4.6476E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2707.93638828713        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1286            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0041E+00  2.5596E+00  2.7911E-01  6.9926E-02  2.0238E+00  9.2939E-01  7.9359E-01  3.5665E-01  4.2779E+00  1.7217E+00
             2.5995E+00
 PARAMETER:  1.0409E-01  1.0399E+00 -1.1762E+00 -2.5603E+00  8.0499E-01  2.6773E-02 -1.3119E-01 -9.3099E-01  1.5535E+00  6.4329E-01
             1.0553E+00
 GRADIENT:   6.8110E+01  4.8615E+02  3.2964E-01  4.3176E+00  3.0351E+01  6.2772E+00  5.5271E+00  1.3785E-01  6.7764E+00  7.3951E+00
             1.9122E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2707.99520916895        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1422
 NPARAMETR:  1.0049E+00  2.5586E+00  2.7817E-01  6.9853E-02  2.0261E+00  9.2959E-01  7.9131E-01  1.3422E-01  4.2752E+00  1.7203E+00
             2.6002E+00
 PARAMETER:  1.0493E-01  1.0395E+00 -1.1795E+00 -2.5614E+00  8.0610E-01  2.6992E-02 -1.3406E-01 -1.9083E+00  1.5528E+00  6.4250E-01
             1.0556E+00
 GRADIENT:   1.8949E+00 -2.3921E+00 -2.6913E-01  1.2006E-01  4.7832E-01  1.1588E-01 -1.0160E+00  1.9662E-02  4.3994E-01 -2.5906E-01
             9.0592E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2708.00786065715        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1597
 NPARAMETR:  1.0042E+00  2.5594E+00  2.8373E-01  7.0159E-02  2.0248E+00  9.2940E-01  7.9379E-01  3.3597E-02  4.2827E+00  1.7204E+00
             2.5988E+00
 PARAMETER:  1.0414E-01  1.0398E+00 -1.1597E+00 -2.5570E+00  8.0545E-01  2.6781E-02 -1.3093E-01 -3.2933E+00  1.5546E+00  6.4255E-01
             1.0551E+00
 GRADIENT:  -3.0386E-02  6.1797E-01  1.1515E-01 -3.9708E-01  2.9774E-01 -9.5932E-03  4.0489E-01  1.2469E-03 -4.2052E-01  5.0183E-02
             8.4893E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -2708.00952409180        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:     1722
 NPARAMETR:  1.0041E+00  2.5592E+00  2.8301E-01  7.0254E-02  2.0260E+00  9.2938E-01  7.9378E-01  1.0358E-02  4.2779E+00  1.7222E+00
             2.5988E+00
 PARAMETER:  1.0414E-01  1.0396E+00 -1.1623E+00 -2.5554E+00  8.0596E-01  2.6758E-02 -1.3095E-01 -4.4257E+00  1.5536E+00  6.4358E-01
             1.0550E+00
 GRADIENT:  -3.1869E-01 -3.7227E+01 -3.3993E-02  1.6242E+01 -6.1560E+01 -3.4696E-02  2.7061E-01  1.3021E-04  2.7417E+01 -6.7121E-02
             3.0535E-02
 NUMSIGDIG:         2.9         2.5         3.3         2.4         2.3         3.0         2.8         0.4         2.4         3.6
                    4.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1722
 NO. OF SIG. DIGITS IN FINAL EST.:  0.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4372E-03 -2.2489E-02 -5.1680E-06  3.1719E-02 -2.3543E-02
 SE:             2.9357E-02  2.7214E-02  1.3452E-05  1.3987E-02  2.5891E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6095E-01  4.0859E-01  7.0084E-01  2.3342E-02  3.6317E-01

 ETASHRINKSD(%)  1.6516E+00  8.8281E+00  9.9955E+01  5.3142E+01  1.3263E+01
 ETASHRINKVR(%)  3.2760E+00  1.6877E+01  1.0000E+02  7.8043E+01  2.4768E+01
 EBVSHRINKSD(%)  1.8361E+00  6.5533E+00  9.9919E+01  6.7602E+01  9.0964E+00
 EBVSHRINKVR(%)  3.6385E+00  1.2677E+01  1.0000E+02  8.9504E+01  1.7365E+01
 RELATIVEINF(%)  9.6311E+01  4.2528E+01  5.6874E-05  4.8002E+00  6.9901E+01
 EPSSHRINKSD(%)  1.6553E+01
 EPSSHRINKVR(%)  3.0366E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          879
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1615.4939413738146     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2708.0095240917958     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1092.5155827179813     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    47.93
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2708.010       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.56E+00  2.83E-01  7.03E-02  2.03E+00  9.29E-01  7.94E-01  1.08E-02  4.28E+00  1.72E+00  2.60E+00
 


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
+        1.23E+03
 
 TH 2
+       -1.09E+01  4.17E+02
 
 TH 3
+       -2.47E+01  2.86E+00  1.17E+04
 
 TH 4
+        1.59E+00  2.48E+02 -1.05E+02  4.01E+04
 
 TH 5
+       -3.08E+00 -1.18E+01  2.33E+03  4.59E+01  5.50E+02
 
 TH 6
+        6.64E+00 -2.75E+00 -1.26E+01 -7.20E+00 -1.01E+00  2.15E+02
 
 TH 7
+       -2.72E+01 -2.41E+01  3.68E+04  1.32E+02  7.38E+03 -1.77E+01  1.06E+03
 
 TH 8
+        6.97E-02 -8.90E-03 -4.76E-02  1.76E-01 -1.09E-02  1.16E-01  8.96E-02  1.48E+00
 
 TH 9
+        4.33E-01 -2.20E-01  6.06E-01  1.08E+03  4.10E-01 -2.06E-01  6.86E+00  5.26E-03  3.04E+01
 
 TH10
+       -1.10E+01 -3.69E+00  3.45E+03  3.24E+01  6.91E+02 -5.20E+00  1.10E+04 -2.86E-02  1.05E+00  1.07E+03
 
 TH11
+       -2.01E+01 -6.35E+00  1.43E+03 -9.67E+01  2.86E+02  1.65E+00  4.52E+03 -2.85E-03 -1.22E+00  4.28E+02  3.43E+02
 
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
 #CPUT: Total CPU Time in Seconds,       62.076
Stop Time:
Wed Sep 29 04:55:35 CDT 2021
