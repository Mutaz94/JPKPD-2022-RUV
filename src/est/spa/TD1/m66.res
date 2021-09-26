Sat Sep 25 13:02:05 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat66.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1641.90814780625        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4162E+02 -7.7522E+01 -5.3836E+00 -9.4745E+01 -4.6807E+00  1.2033E+01 -7.6802E+00  8.5816E+00  2.0495E+01 -1.0583E-01
            -7.5414E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1651.86954544100        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.6580E-01  1.1734E+00  1.0460E+00  9.5729E-01  1.1208E+00  9.4371E-01  1.1414E+00  8.9753E-01  7.5325E-01  1.0607E+00
             1.0512E+00
 PARAMETER:  6.5202E-02  2.5987E-01  1.4497E-01  5.6354E-02  2.1402E-01  4.2065E-02  2.3229E-01 -8.1041E-03 -1.8336E-01  1.5897E-01
             1.4990E-01
 GRADIENT:   6.2028E+01  7.2551E+00  1.0778E+01 -1.3652E+01  6.7726E+00 -5.8702E+00 -2.3117E+00 -8.4213E-01 -8.3157E+00 -7.8086E+00
             4.6935E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1653.43766773683        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.7125E-01  8.9677E-01  1.0148E+00  1.1503E+00  9.6049E-01  9.3324E-01  1.5039E+00  5.2679E-01  6.6272E-01  9.9267E-01
             1.0406E+00
 PARAMETER:  7.0825E-02 -8.9610E-03  1.1469E-01  2.4002E-01  5.9688E-02  3.0910E-02  5.0803E-01 -5.4095E-01 -3.1141E-01  9.2647E-02
             1.3982E-01
 GRADIENT:   7.9773E+01  2.9393E+01  1.3762E+01  4.9446E+01 -2.2950E+00 -1.0706E+01  3.2874E+00 -1.2224E+00 -8.5150E+00 -6.1128E+00
             8.0035E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1655.46135991708        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.4765E-01  9.4184E-01  8.2546E-01  1.0908E+00  8.8107E-01  9.4753E-01  1.3702E+00  3.0585E-01  7.4275E-01  9.1749E-01
             1.0231E+00
 PARAMETER:  4.6233E-02  4.0075E-02 -9.1815E-02  1.8688E-01 -2.6618E-02  4.6104E-02  4.1495E-01 -1.0847E+00 -1.9739E-01  1.3887E-02
             1.2283E-01
 GRADIENT:   1.4437E+01  4.2667E+00 -6.3967E+00  1.5889E+01  8.9425E+00 -4.2086E+00  1.9504E+00  5.4374E-01  1.1618E+00  5.5897E-01
             5.8362E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1655.46583678396        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.4519E-01  9.3991E-01  8.0982E-01  1.0875E+00  8.6924E-01  9.5187E-01  1.3647E+00  2.8003E-01  7.4274E-01  9.0396E-01
             1.0225E+00
 PARAMETER:  4.3629E-02  3.8029E-02 -1.1095E-01  1.8390E-01 -4.0135E-02  5.0678E-02  4.1090E-01 -1.1729E+00 -1.9741E-01 -9.7311E-04
             1.2228E-01
 GRADIENT:   7.9682E+00  2.0783E+00 -4.8650E+00  9.6497E+00  6.0626E+00 -2.4908E+00  1.2642E+00  4.8903E-01  1.0262E+00  6.6240E-01
             4.9669E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1655.47009633128        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.4342E-01  9.3491E-01  7.9670E-01  1.0869E+00  8.5778E-01  9.5520E-01  1.3648E+00  2.4380E-01  7.4112E-01  8.9188E-01
             1.0222E+00
 PARAMETER:  4.1760E-02  3.2694E-02 -1.2728E-01  1.8329E-01 -5.3412E-02  5.4169E-02  4.1100E-01 -1.3114E+00 -1.9959E-01 -1.4419E-02
             1.2192E-01
 GRADIENT:   3.3914E+00  5.4037E-01 -3.1815E+00  4.6717E+00  3.4816E+00 -1.2060E+00  7.0635E-01  3.8215E-01  7.6655E-01  5.9252E-01
             3.5716E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1655.47295304664        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  9.4216E-01  9.2912E-01  7.8608E-01  1.0874E+00  8.4778E-01  9.5771E-01  1.3674E+00  2.0300E-01  7.3912E-01  8.8202E-01
             1.0219E+00
 PARAMETER:  4.0419E-02  2.6480E-02 -1.4070E-01  1.8382E-01 -6.5137E-02  5.6792E-02  4.1288E-01 -1.4945E+00 -2.0229E-01 -2.5542E-02
             1.2170E-01
 GRADIENT:   1.3673E-01 -5.2353E-01 -1.7482E+00  9.6853E-01  1.4327E+00 -2.5616E-01  2.7188E-01  2.6951E-01  5.0929E-01  4.8667E-01
             2.2769E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1655.47898176455        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      509
 NPARAMETR:  9.4121E-01  9.2285E-01  7.7714E-01  1.0888E+00  8.3880E-01  9.5968E-01  1.3714E+00  1.5452E-01  7.3694E-01  8.7365E-01
             1.0218E+00
 PARAMETER:  3.9411E-02  1.9715E-02 -1.5214E-01  1.8507E-01 -7.5782E-02  5.8846E-02  4.1586E-01 -1.7675E+00 -2.0524E-01 -3.5071E-02
             1.2155E-01
 GRADIENT:  -2.2958E+00 -1.3042E+00 -4.8069E-01 -1.9418E+00 -2.7284E-01  4.7980E-01 -8.4801E-02  1.5790E-01  2.6128E-01  3.5917E-01
             1.0223E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1655.63697240486        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  9.4245E-01  9.1143E-01  7.7204E-01  1.0956E+00  8.3189E-01  9.5792E-01  1.3896E+00  1.0000E-02  7.3345E-01  8.7403E-01
             1.0223E+00
 PARAMETER:  4.0724E-02  7.2619E-03 -1.5872E-01  1.9128E-01 -8.4057E-02  5.7005E-02  4.2902E-01 -6.3643E+00 -2.1000E-01 -3.4643E-02
             1.2207E-01
 GRADIENT:   5.9347E-01 -1.2016E+00 -2.6757E+00  1.0288E+00  2.2863E+00 -1.9740E-01  4.6317E-01  0.0000E+00  4.6982E-01  8.9446E-01
             4.4991E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1656.56506705397        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      766
 NPARAMETR:  9.5561E-01  7.2594E-01  8.4075E-01  1.2199E+00  7.9929E-01  9.5267E-01  1.6998E+00  1.0000E-02  6.8771E-01  8.8653E-01
             1.0293E+00
 PARAMETER:  5.4596E-02 -2.2029E-01 -7.3463E-02  2.9880E-01 -1.2403E-01  5.1511E-02  6.3053E-01 -1.4806E+01 -2.7438E-01 -2.0440E-02
             1.2887E-01
 GRADIENT:   3.8255E+00  5.1409E+00 -2.5095E+00  1.5018E+01  4.9048E+00 -4.9892E+00  6.9671E-02  0.0000E+00  6.1211E-01 -1.2310E+00
             1.9278E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1656.79636538188        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      941
 NPARAMETR:  9.5303E-01  6.2150E-01  8.3474E-01  1.2696E+00  7.5836E-01  9.6381E-01  1.9144E+00  1.0000E-02  6.6250E-01  8.7365E-01
             1.0232E+00
 PARAMETER:  5.1891E-02 -3.7562E-01 -8.0633E-02  3.3869E-01 -1.7660E-01  6.3143E-02  7.4942E-01 -1.2164E+01 -3.1173E-01 -3.5076E-02
             1.2290E-01
 GRADIENT:  -1.4912E-01  2.0570E-01  3.4409E-01  7.0656E-02 -5.4689E-01  2.5499E-02 -1.0184E-02  0.0000E+00 -3.7102E-02  4.3835E-03
            -2.3699E-02

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1656.79653566016        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1033
 NPARAMETR:  9.5306E-01  6.1899E-01  8.3553E-01  1.2709E+00  7.5827E-01  9.6367E-01  1.9202E+00  1.0000E-02  6.6215E-01  8.7416E-01
             1.0232E+00
 PARAMETER:  5.1922E-02 -3.7967E-01 -7.9685E-02  3.3972E-01 -1.7672E-01  6.2992E-02  7.5245E-01 -1.2197E+01 -3.1226E-01 -3.4487E-02
             1.2296E-01
 GRADIENT:   1.5149E-02 -6.5482E-03  3.3729E-02 -6.2128E-02 -3.9263E-02 -1.9413E-02 -5.0296E-04  0.0000E+00  3.3272E-03  1.9488E-03
             8.2363E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1033
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1977E-04  2.2320E-02 -4.8509E-04 -2.5354E-02 -7.0464E-04
 SE:             2.9827E-02  2.1510E-02  2.1475E-04  2.2893E-02  2.3309E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8877E-01  2.9943E-01  2.3893E-02  2.6808E-01  9.7588E-01

 ETASHRINKSD(%)  7.5903E-02  2.7940E+01  9.9281E+01  2.3305E+01  2.1911E+01
 ETASHRINKVR(%)  1.5175E-01  4.8073E+01  9.9995E+01  4.1179E+01  3.9021E+01
 EBVSHRINKSD(%)  4.8320E-01  2.8781E+01  9.9318E+01  2.2165E+01  1.9608E+01
 EBVSHRINKVR(%)  9.6406E-01  4.9278E+01  9.9995E+01  3.9417E+01  3.5371E+01
 RELATIVEINF(%)  9.8427E+01  6.7136E+00  4.3953E-04  8.7016E+00  5.4258E+00
 EPSSHRINKSD(%)  4.2620E+01
 EPSSHRINKVR(%)  6.7076E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1656.7965356601621     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -921.64570909642396     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.95
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1656.797       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  6.19E-01  8.36E-01  1.27E+00  7.58E-01  9.64E-01  1.92E+00  1.00E-02  6.62E-01  8.74E-01  1.02E+00
 


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
+        1.31E+03
 
 TH 2
+       -1.23E+01  4.22E+02
 
 TH 3
+        2.14E+01  1.83E+02  8.14E+02
 
 TH 4
+       -9.64E+00  3.79E+02 -3.25E+02  9.78E+02
 
 TH 5
+       -6.56E+00 -3.53E+02 -1.04E+03  3.54E+02  1.64E+03
 
 TH 6
+       -3.73E+00 -1.95E+00  6.35E+00 -2.69E+00 -3.80E+00  2.14E+02
 
 TH 7
+        1.41E+00  3.61E+01 -2.35E+00 -1.39E+01  5.95E-01 -1.19E-02  2.00E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.93E-01 -2.04E+01 -4.24E+01 -1.81E+01  4.31E+01  1.23E+00  1.59E+01  0.00E+00  1.81E+02
 
 TH10
+       -2.82E+00 -4.97E+00 -6.96E+01 -2.94E+01 -5.09E+01  2.08E+00  3.26E+00  0.00E+00  1.51E+01  1.10E+02
 
 TH11
+       -5.50E+00 -1.02E+01 -4.30E+01 -6.41E+00  1.83E+01  2.03E+00  2.04E+00  0.00E+00  1.66E+01  2.57E+01  2.16E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.112
Stop Time:
Sat Sep 25 13:02:23 CDT 2021
