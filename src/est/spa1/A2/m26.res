Wed Sep 29 23:12:00 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat26.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1311.55389511598        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1924E+02  3.2607E+01  1.0327E+02 -1.3079E+00  9.5092E+01  5.6226E+01 -5.5085E+01 -2.5955E+02 -1.8517E+01 -7.9518E+01
            -1.1417E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1747.12513281885        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1438E+00  1.1756E+00  9.4739E-01  9.6208E-01  1.0284E+00  1.0364E+00  1.0678E+00  1.2467E+00  7.4590E-01  7.9292E-01
             2.4623E+00
 PARAMETER:  2.3436E-01  2.6179E-01  4.5955E-02  6.1340E-02  1.2804E-01  1.3573E-01  1.6556E-01  3.2053E-01 -1.9317E-01 -1.3204E-01
             1.0011E+00
 GRADIENT:   4.1993E+02  2.3527E+01 -3.5130E-01  2.4194E+01 -3.1866E+00  2.2791E+01 -3.0839E+00  2.0532E+00 -1.0231E+01  6.3672E+00
             1.0370E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1769.26949201551        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.1004E+00  1.0053E+00  3.5666E-01  9.8414E-01  5.4232E-01  9.9244E-01  1.1322E+00  1.7730E+00  8.9209E-01  1.5260E-01
             2.1305E+00
 PARAMETER:  1.9565E-01  1.0533E-01 -9.3096E-01  8.4011E-02 -5.1189E-01  9.2407E-02  2.2417E-01  6.7269E-01 -1.4193E-02 -1.7800E+00
             8.5637E-01
 GRADIENT:   3.5510E+02  5.1965E+01  2.1048E+01  4.4788E+01 -3.3830E+01  1.8754E+01  1.4067E+01  4.6063E+01  9.7272E+00  1.0007E+00
             8.0863E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1787.80757382271        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      276
 NPARAMETR:  1.0211E+00  8.9692E-01  3.2640E-01  1.0234E+00  4.9771E-01  8.9744E-01  1.1154E+00  1.5417E+00  9.1655E-01  1.4087E-01
             2.0067E+00
 PARAMETER:  1.2089E-01 -8.7851E-03 -1.0196E+00  1.2311E-01 -5.9775E-01 -8.2041E-03  2.0920E-01  5.3290E-01  1.2857E-02 -1.8599E+00
             7.9650E-01
 GRADIENT:   4.6062E+01  1.6103E+01 -1.3266E+01  4.4562E+01  3.8534E+00 -6.0791E+00  4.7703E+00  3.2236E+01  1.0927E+01  5.4225E-01
             4.6502E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1796.07673822107        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  9.9455E-01  8.7149E-01  3.1056E-01  1.0061E+00  4.8379E-01  8.9844E-01  1.1751E+00  1.0136E+00  8.5725E-01  1.3062E-01
             1.9251E+00
 PARAMETER:  9.4537E-02 -3.7546E-02 -1.0694E+00  1.0609E-01 -6.2611E-01 -7.0910E-03  2.6132E-01  1.1350E-01 -5.4026E-02 -1.9354E+00
             7.5499E-01
 GRADIENT:  -2.0942E+01  4.4501E-01 -1.5576E+01  1.0680E+01  2.5502E+01 -6.3348E+00  3.3248E-02 -3.7186E-01  1.6699E+00 -1.0274E-01
             2.5504E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1803.34220870068        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      628
 NPARAMETR:  1.0021E+00  5.5079E-01  2.3575E-01  1.0759E+00  3.2394E-01  9.2305E-01  1.3029E+00  1.2386E+00  9.5109E-01  4.4972E-02
             1.6529E+00
 PARAMETER:  1.0213E-01 -4.9641E-01 -1.3450E+00  1.7312E-01 -1.0272E+00  1.9927E-02  3.6458E-01  3.1399E-01  4.9854E-02 -3.0017E+00
             6.0253E-01
 GRADIENT:   5.4645E+00  1.0906E+01  4.1123E+00  1.0166E+00 -1.2745E+01  9.4124E-01 -2.4367E+00 -3.9230E-01 -1.8205E+00 -2.1164E-01
            -3.0450E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1805.75491648750        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  9.8349E-01  5.2011E-01  2.2719E-01  1.0698E+00  3.1027E-01  9.1636E-01  1.3130E+00  1.2147E+00  9.5700E-01  2.9061E-01
             1.6458E+00
 PARAMETER:  8.3357E-02 -5.5372E-01 -1.3820E+00  1.6744E-01 -1.0703E+00  1.2651E-02  3.7230E-01  2.9449E-01  5.6053E-02 -1.1358E+00
             5.9824E-01
 GRADIENT:   7.3541E+01  3.1620E+01  5.7885E+01  2.8314E+01  1.6323E+02  6.0862E+00  6.5375E+00  9.5769E+00  1.6924E+00 -7.7270E-03
             1.1076E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1807.07651076192        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      926
 NPARAMETR:  9.8356E-01  5.0467E-01  2.1951E-01  1.0695E+00  3.0238E-01  9.1347E-01  1.2324E+00  1.0995E+00  9.5324E-01  4.5688E-01
             1.6198E+00
 PARAMETER:  8.3421E-02 -5.8385E-01 -1.4164E+00  1.6721E-01 -1.0961E+00  9.4965E-03  3.0894E-01  1.9489E-01  5.2108E-02 -6.8333E-01
             5.8229E-01
 GRADIENT:  -4.1908E+01  1.2780E+00  2.2699E+00 -1.4114E+01  2.9208E+01 -3.4167E+00  1.6477E+00  4.0081E+00 -2.3051E+00  3.2208E+00
            -8.2488E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1813.71138896675        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1102
 NPARAMETR:  9.7978E-01  3.2060E-01  1.2748E-01  9.7640E-01  1.9213E-01  9.0871E-01  7.3095E-01  3.9373E-01  9.9239E-01  8.3238E-01
             1.7366E+00
 PARAMETER:  7.9573E-02 -1.0375E+00 -1.9598E+00  7.6120E-02 -1.5496E+00  4.2761E-03 -2.1341E-01 -8.3209E-01  9.2361E-02 -8.3470E-02
             6.5191E-01
 GRADIENT:  -5.6465E+01  1.2217E+01  2.2303E+01  5.4906E+01 -5.5794E+01 -4.1203E+00 -3.2836E-01 -2.4527E+00 -1.5564E+01  2.1343E+01
             1.2351E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1818.91691079348        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1280
 NPARAMETR:  9.9935E-01  2.9224E-01  1.0809E-01  8.9791E-01  1.7758E-01  9.1479E-01  7.0690E-01  3.8094E-01  1.1110E+00  7.6932E-01
             1.6763E+00
 PARAMETER:  9.9353E-02 -1.1302E+00 -2.1248E+00 -7.6811E-03 -1.6284E+00  1.0938E-02 -2.4686E-01 -8.6510E-01  2.0525E-01 -1.6224E-01
             6.1656E-01
 GRADIENT:   3.9892E-01  1.4829E+00  4.9746E+00  1.1545E+00  3.0557E+00 -5.2233E-01 -2.7767E+00 -4.4430E+00  5.8024E-01 -2.4943E+00
            -4.5974E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1819.41618582655        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1440
 NPARAMETR:  9.9915E-01  2.9018E-01  1.0826E-01  8.9939E-01  1.7777E-01  9.1591E-01  9.3676E-01  3.8120E-01  1.1108E+00  7.6859E-01
             1.6770E+00
 PARAMETER:  9.9146E-02 -1.1373E+00 -2.1233E+00 -6.0426E-03 -1.6272E+00  1.2160E-02  3.4673E-02 -8.6444E-01  2.0510E-01 -1.6320E-01
             6.1701E-01
 GRADIENT:  -5.6044E-01  1.3625E+00  4.6466E+00 -3.6243E-02  7.8838E+00  2.7368E-02  1.2018E-01 -4.3012E+00  2.0708E+00  2.4795E-02
             1.2840E+01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1819.41618582655        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     1465
 NPARAMETR:  9.9915E-01  2.9018E-01  1.0826E-01  8.9939E-01  1.7777E-01  9.1591E-01  9.3676E-01  3.8120E-01  1.1108E+00  7.6859E-01
             1.6770E+00
 PARAMETER:  9.9146E-02 -1.1373E+00 -2.1233E+00 -6.0426E-03 -1.6272E+00  1.2160E-02  3.4673E-02 -8.6444E-01  2.0510E-01 -1.6320E-01
             6.1701E-01
 GRADIENT:   3.4389E+04  3.0256E+03  1.6571E+03  3.4367E+04  2.1294E+03 -3.4386E+04 -3.4389E+04  3.9481E+03  8.3614E+03 -2.1063E+04
             2.7952E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1465
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.4859E-04  1.6094E-02  4.0248E-03 -5.5204E-03  4.6389E-03
 SE:             2.9568E-02  1.5355E-02  9.1485E-03  2.6840E-02  2.6721E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8250E-01  2.9460E-01  6.5998E-01  8.3704E-01  8.6218E-01

 ETASHRINKSD(%)  9.4364E-01  4.8558E+01  6.9351E+01  1.0081E+01  1.0480E+01
 ETASHRINKVR(%)  1.8784E+00  7.3537E+01  9.0607E+01  1.9146E+01  1.9862E+01
 EBVSHRINKSD(%)  1.1108E+00  4.9391E+01  7.3646E+01  7.9072E+00  1.0747E+01
 EBVSHRINKVR(%)  2.2093E+00  7.4388E+01  9.3055E+01  1.5189E+01  2.0340E+01
 RELATIVEINF(%)  9.7727E+01  6.5567E+00  6.6011E-01  2.5387E+01  6.5460E+00
 EPSSHRINKSD(%)  3.5836E+01
 EPSSHRINKVR(%)  5.8830E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1819.4161858265547     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -900.47765262188204     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1819.416       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  2.90E-01  1.08E-01  8.99E-01  1.78E-01  9.16E-01  9.37E-01  3.81E-01  1.11E+00  7.69E-01  1.68E+00
 


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
+        8.61E+06
 
 TH 2
+        2.61E+06  1.58E+06
 
 TH 3
+       -6.00E+01  1.11E+06  1.61E+06
 
 TH 4
+       -1.22E+04 -3.84E+03 -7.67E+03  1.06E+07
 
 TH 5
+       -5.95E+06 -1.81E+06  1.24E+06 -4.79E+03  1.07E+06
 
 TH 6
+        1.29E+03  3.68E+02  6.26E+02  1.43E+03  4.36E+02  1.02E+07
 
 TH 7
+        2.07E+01  6.76E+01 -9.27E+00  1.30E+04  3.16E+06 -1.38E+03  9.80E+06
 
 TH 8
+       -3.36E+04 -1.02E+04  1.10E+06 -3.68E+03 -9.05E+05  3.87E+02  2.80E+06  7.81E+05
 
 TH 9
+       -2.16E+04 -6.56E+03 -9.17E+03 -5.36E+03 -1.30E+06  5.69E+02  2.31E+04 -6.48E+03  1.65E+06
 
 TH10
+        5.02E+03  1.55E+03  2.37E+03  5.59E+03  1.67E+03  7.48E+06 -5.35E+03  1.54E+03  2.21E+03  5.46E+06
 
 TH11
+        4.77E+01  2.52E+05 -3.69E+05 -1.17E+03 -5.74E+05  1.26E+02 -8.88E+05  6.26E+02 -2.06E+03  4.89E+02  8.02E+04
 
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
 #CPUT: Total CPU Time in Seconds,       28.889
Stop Time:
Wed Sep 29 23:12:30 CDT 2021
