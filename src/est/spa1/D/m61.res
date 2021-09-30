Thu Sep 30 03:21:36 CDT 2021
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
$DATA ../../../../data/spa1/D/dat61.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   15579.8284765041        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.6517E+02  3.2645E+02 -3.1980E+01  2.0741E+02  4.2941E+02 -2.3866E+03 -7.3647E+02 -4.0537E+01 -1.2972E+03 -7.5982E+02
            -2.9801E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -677.973922434082        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2491E+00  9.4820E-01  8.2625E-01  2.2190E+00  1.2714E+00  2.9279E+00  1.3286E+00  9.5042E-01  2.1141E+00  1.3268E+00
             1.2745E+01
 PARAMETER:  3.2241E-01  4.6814E-02 -9.0854E-02  8.9704E-01  3.4008E-01  1.1743E+00  3.8416E-01  4.9146E-02  8.4862E-01  3.8276E-01
             2.6452E+00
 GRADIENT:  -2.6695E+01  4.0633E+01 -3.3910E+01  1.0677E+02 -2.1542E+00  7.4034E+01  4.0736E-01  7.6346E+00 -2.0523E+01  4.8691E+00
             1.7590E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -715.027355242571        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.1891E+00  9.1325E-01  3.4751E+00  2.4639E+00  4.4575E+00  3.3447E+00  6.6021E+00  3.3879E-01  3.2292E+00  7.7362E+00
             1.0027E+01
 PARAMETER:  2.7323E-01  9.2557E-03  1.3456E+00  1.0017E+00  1.5946E+00  1.3074E+00  1.9874E+00 -9.8238E-01  1.2722E+00  2.1459E+00
             2.4052E+00
 GRADIENT:  -1.1461E+01  1.5029E+01  1.1876E+01  6.5780E+01 -1.0426E+01  1.2863E+02  1.5861E+01 -1.6317E-01  6.4918E+01  2.0710E+01
             1.2979E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -753.406462623271        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.2437E+00  9.7793E-01  1.2119E+00  1.4921E+00  4.3967E+00  2.3355E+00  1.0988E+00  3.0660E-02  2.4085E+00  6.8053E+00
             1.0830E+01
 PARAMETER:  3.1813E-01  7.7681E-02  2.9218E-01  5.0018E-01  1.5809E+00  9.4822E-01  1.9424E-01 -3.3848E+00  9.7902E-01  2.0177E+00
             2.4824E+00
 GRADIENT:  -1.0680E+01  2.0196E+01  1.9920E+01  8.6759E-01 -1.3544E+01 -1.2566E+01  4.7231E+00 -3.9478E-03 -8.9056E+00  5.7091E+00
             1.7674E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -821.068259762629        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0026E+00  1.3810E-01  1.7639E-01  1.1400E+00  1.8495E+01  1.7196E+00  1.1326E-01  8.1043E-02  1.2062E+00  1.0416E+01
             7.8038E+00
 PARAMETER:  1.0258E-01 -1.8798E+00 -1.6350E+00  2.3103E-01  3.0175E+00  6.4207E-01 -2.0781E+00 -2.4128E+00  2.8744E-01  2.4433E+00
             2.1546E+00
 GRADIENT:   1.0499E+02  9.5034E+00 -3.2013E+01  1.8554E+02 -1.5133E+01 -1.2355E+02  1.2221E-01 -1.0873E-01 -3.3341E+01  3.6832E+01
            -1.4608E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -891.018039329116        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  7.5663E-01  4.7405E-02  8.6196E-02  7.2563E-01  2.6010E+01  1.9061E+00  1.0000E-02  1.0386E+00  7.5557E-01  1.0528E+01
             7.1605E+00
 PARAMETER: -1.7888E-01 -2.9490E+00 -2.3511E+00 -2.2071E-01  3.3585E+00  7.4504E-01 -5.0771E+00  1.3788E-01 -1.8028E-01  2.4541E+00
             2.0686E+00
 GRADIENT:   6.3938E+01  2.9333E+00 -4.5765E+01  1.9300E+02  8.7865E-01  1.1052E+00  0.0000E+00 -1.3073E+01 -5.1299E+00 -5.3042E+00
            -1.7811E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -898.958494366102        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      456
 NPARAMETR:  6.3591E-01  1.9536E-02  4.1767E-02  4.5311E-01  3.9145E+01  1.7768E+00  1.0000E-02  8.5803E-01  3.8239E-01  1.2083E+01
             7.5160E+00
 PARAMETER: -3.5270E-01 -3.8355E+00 -3.0756E+00 -6.9161E-01  3.7673E+00  6.7482E-01 -7.3154E+00 -5.3118E-02 -8.6132E-01  2.5918E+00
             2.1170E+00
 GRADIENT:   1.2324E+02 -1.6022E+00 -1.1020E+02  2.8967E+02  2.7597E-01 -2.5067E+01  0.0000E+00 -3.6424E+01 -1.2807E+01 -6.6001E-01
            -1.4531E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -931.674931471681        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      639
 NPARAMETR:  5.8414E-01  1.3791E-02  4.7513E-02  4.4873E-01  3.6481E+01  1.9460E+00  1.0000E-02  1.1257E+00  2.6950E-01  1.2225E+01
             8.5033E+00
 PARAMETER: -4.3762E-01 -4.1837E+00 -2.9468E+00 -7.0133E-01  3.6968E+00  7.6579E-01 -8.1317E+00  2.1842E-01 -1.2112E+00  2.6035E+00
             2.2405E+00
 GRADIENT:   1.4735E+00 -1.0110E-01  4.4525E-01  2.7439E+00  2.8329E-01 -4.6591E-01  0.0000E+00 -2.3785E+00  4.8221E-01  1.6032E+00
             1.9306E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -936.804857415260        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      824             RESET HESSIAN, TYPE I
 NPARAMETR:  5.5453E-01  1.0841E-02  4.2121E-02  4.1106E-01  3.6061E+01  1.9518E+00  1.0000E-02  1.1326E+00  2.2785E-01  1.3345E+01
             8.2539E+00
 PARAMETER: -4.8963E-01 -4.4244E+00 -3.0672E+00 -7.8902E-01  3.6852E+00  7.6874E-01 -8.3250E+00  2.2454E-01 -1.3791E+00  2.6911E+00
             2.2107E+00
 GRADIENT:   5.0828E+01 -3.0389E-02  7.5502E+01  2.5653E+01  4.1432E-02  2.9837E+01  0.0000E+00  1.0962E+00  6.7242E-01  1.5497E-02
             1.9131E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -936.869526178347        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      985
 NPARAMETR:  5.5331E-01  1.1037E-02  4.2015E-02  4.1067E-01  3.4537E+01  1.9455E+00  1.0000E-02  1.1330E+00  1.8350E-01  1.2279E+01
             8.2634E+00
 PARAMETER: -4.9184E-01 -4.4065E+00 -3.0697E+00 -7.8997E-01  3.6420E+00  7.6554E-01 -8.3250E+00  2.2491E-01 -1.5955E+00  2.6079E+00
             2.2118E+00
 GRADIENT:  -1.2283E+00 -3.2494E-02 -4.9554E-01  6.9803E-01  3.3508E-02 -2.0977E-01  0.0000E+00 -1.0798E+00  1.6380E-01 -4.6377E-03
            -2.9670E+00

0ITERATION NO.:   48    OBJECTIVE VALUE:  -936.876604703276        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1086
 NPARAMETR:  5.5619E-01  1.0754E-02  4.2802E-02  4.1044E-01  3.4921E+01  1.9568E+00  1.0000E-02  1.1332E+00  1.7713E-01  1.2598E+01
             8.3914E+00
 PARAMETER: -4.9155E-01 -4.3998E+00 -3.0695E+00 -7.8999E-01  3.6262E+00  7.6566E-01 -8.3250E+00  2.2677E-01 -1.6188E+00  2.6533E+00
             2.2120E+00
 GRADIENT:  -9.2426E-01  5.1310E+01 -7.2062E+01  7.4257E-01 -6.2176E+01 -2.9587E+02  0.0000E+00  9.9535E+02  6.9942E+01  8.5994E+01
            -1.0535E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1086
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.6121E-03 -7.5785E-07  6.8575E-03 -6.8331E-03 -1.2115E-03
 SE:             2.9300E-02  8.8212E-07  2.3436E-02  5.4000E-03  6.0418E-04
 N:                     100         100         100         100         100

 P VAL.:         8.2146E-01  3.9027E-01  7.6982E-01  2.0574E-01  4.4943E-02

 ETASHRINKSD(%)  1.8418E+00  9.9997E+01  2.1488E+01  8.1909E+01  9.7976E+01
 ETASHRINKVR(%)  3.6496E+00  1.0000E+02  3.8358E+01  9.6727E+01  9.9959E+01
 EBVSHRINKSD(%)  1.9738E+00  9.9996E+01  2.1931E+01  8.2168E+01  9.8574E+01
 EBVSHRINKVR(%)  3.9086E+00  1.0000E+02  3.9053E+01  9.6820E+01  9.9980E+01
 RELATIVEINF(%)  6.7446E+00  8.3952E-09  5.5641E-01  2.0063E-02  2.3959E-03
 EPSSHRINKSD(%)  1.1322E+01
 EPSSHRINKVR(%)  2.1362E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -936.87660470327648     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -17.938071498603790     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.90
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -936.877       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.53E-01  1.11E-02  4.20E-02  4.11E-01  3.40E+01  1.95E+00  1.00E-02  1.14E+00  1.79E-01  1.28E+01  8.26E+00
 


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
+        9.12E+02
 
 TH 2
+       -2.70E+02  2.11E+06
 
 TH 3
+       -1.63E+05 -7.93E+05  6.30E+05
 
 TH 4
+        6.37E+04  1.93E+03  9.57E+04  6.19E+03
 
 TH 5
+        2.25E-02  3.30E+00 -1.66E+00 -1.41E+02  3.71E-01
 
 TH 6
+        6.15E+00  2.02E+01  1.33E+02 -3.45E+01 -2.79E-03  2.59E+03
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -8.07E+04 -1.69E+02 -2.31E+02 -2.68E+01  5.69E-02  2.37E+00  0.00E+00  8.53E+04
 
 TH 9
+       -3.77E-01 -7.95E+01  1.16E+02 -2.62E+01  2.83E-02  3.05E+00  0.00E+00  1.45E+00  6.71E+04
 
 TH10
+       -2.12E-02 -3.53E+02 -2.43E+01  5.17E+02 -2.65E-02  2.15E-02  0.00E+00 -2.11E-01 -1.04E-01  4.96E+00
 
 TH11
+        1.13E+03 -9.70E+01  2.27E+03 -9.87E+02  3.64E-02  2.10E+02  0.00E+00 -1.20E+03  1.57E+00 -1.26E-01  2.45E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.812
Stop Time:
Thu Sep 30 03:22:06 CDT 2021
