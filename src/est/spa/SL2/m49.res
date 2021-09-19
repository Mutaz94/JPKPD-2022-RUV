Sat Sep 18 12:18:17 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat49.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1579.31732316455        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4187E+02 -5.3165E+01 -1.3198E+01 -4.5277E+01  2.5226E+01 -2.1566E+01 -3.7712E+01  1.6910E+00 -2.4556E+01 -2.0743E+01
            -1.8757E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1592.87563711016        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.5236E-01  1.1208E+00  1.0388E+00  9.4301E-01  1.0784E+00  1.0427E+00  1.3940E+00  9.6886E-01  1.1000E+00  1.1473E+00
             1.0413E+00
 PARAMETER:  5.1183E-02  2.1404E-01  1.3809E-01  4.1324E-02  1.7551E-01  1.4184E-01  4.3214E-01  6.8363E-02  1.9531E-01  2.3737E-01
             1.4047E-01
 GRADIENT:   3.1833E+01 -7.1106E+00  4.5737E-01 -1.1859E+01 -8.5834E-01  3.6013E+00  7.5113E+00  2.7296E+00  1.4081E+01  1.9832E+00
             3.6884E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1593.38970614841        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      166
 NPARAMETR:  9.6306E-01  1.0884E+00  1.0654E+00  9.9643E-01  1.0784E+00  1.0641E+00  1.4356E+00  8.3735E-01  1.0350E+00  1.1937E+00
             1.0461E+00
 PARAMETER:  6.2356E-02  1.8467E-01  1.6332E-01  9.6419E-02  1.7545E-01  1.6216E-01  4.6155E-01 -7.7515E-02  1.3437E-01  2.7704E-01
             1.4507E-01
 GRADIENT:   5.3572E+01  1.3952E+01  2.0136E+00  2.0424E+01  4.9618E+00  1.2571E+01  5.5595E+00 -6.4817E-01  7.9892E+00  3.5236E+00
             4.8054E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1594.33199681613        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      237
 NPARAMETR:  9.4371E-01  1.0848E+00  8.3642E-01  9.6826E-01  9.5760E-01  1.0477E+00  1.4588E+00  5.1493E-01  9.9218E-01  1.0457E+00
             1.0325E+00
 PARAMETER:  4.2061E-02  1.8139E-01 -7.8624E-02  6.7747E-02  5.6676E-02  1.4663E-01  4.7764E-01 -5.6372E-01  9.2149E-02  1.4465E-01
             1.3196E-01
 GRADIENT:   1.0616E+01  2.3703E+00 -4.9341E+00  7.9958E+00  3.5425E+00  4.8843E+00  4.8588E+00  1.2042E+00  4.9950E+00  2.7132E+00
             3.0524E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1594.35615007859        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  9.4043E-01  1.0971E+00  7.8868E-01  9.5295E-01  9.3493E-01  1.0414E+00  1.4318E+00  4.1615E-01  9.8010E-01  1.0105E+00
             1.0283E+00
 PARAMETER:  3.8579E-02  1.9263E-01 -1.3739E-01  5.1811E-02  3.2721E-02  1.4058E-01  4.5891E-01 -7.7671E-01  7.9898E-02  1.1048E-01
             1.2795E-01
 GRADIENT:   2.5367E+00  4.5074E-01 -3.2353E+00  3.0835E+00  1.5551E+00  1.7223E+00  2.1716E+00  9.2728E-01  2.3033E+00  1.6102E+00
             1.5486E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1594.37009046794        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.3919E-01  1.1023E+00  7.5531E-01  9.4472E-01  9.1704E-01  1.0384E+00  1.4194E+00  3.1974E-01  9.7167E-01  9.8402E-01
             1.0255E+00
 PARAMETER:  3.7258E-02  1.9737E-01 -1.8063E-01  4.3132E-02  1.3391E-02  1.3765E-01  4.5020E-01 -1.0403E+00  7.1257E-02  8.3892E-02
             1.2522E-01
 GRADIENT:  -8.5359E-01 -4.6498E-01 -1.8732E+00  4.1373E-01  4.9260E-01  1.1904E-01  5.8623E-01  6.0951E-01  6.8241E-01  7.4998E-01
             5.6053E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1594.38836104752        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  9.3857E-01  1.1039E+00  7.2936E-01  9.3992E-01  9.0160E-01  1.0365E+00  1.4133E+00  2.1427E-01  9.6475E-01  9.6201E-01
             1.0233E+00
 PARAMETER:  3.6602E-02  1.9884E-01 -2.1558E-01  3.8038E-02 -3.5873E-03  1.3584E-01  4.4594E-01 -1.4405E+00  6.4109E-02  6.1269E-02
             1.2306E-01
 GRADIENT:  -2.7192E+00 -9.9691E-01 -7.9438E-01 -1.2668E+00 -2.0255E-01 -9.0543E-01 -5.1214E-01  3.0421E-01 -4.5501E-01  6.6455E-02
            -1.6616E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1594.42804927831        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  9.3858E-01  1.1034E+00  7.1795E-01  9.3878E-01  8.9415E-01  1.0361E+00  1.4132E+00  1.3574E-01  9.6196E-01  9.5236E-01
             1.0224E+00
 PARAMETER:  3.6616E-02  1.9838E-01 -2.3135E-01  3.6827E-02 -1.1878E-02  1.3549E-01  4.4587E-01 -1.8970E+00  6.1220E-02  5.1184E-02
             1.2218E-01
 GRADIENT:  -2.8926E+00 -1.0213E+00 -3.4623E-01 -1.5565E+00 -3.6097E-01 -1.1310E+00 -8.1525E-01  1.2892E-01 -7.8731E-01 -1.8826E-01
            -4.0615E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1594.43971636858        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  9.3859E-01  1.1034E+00  7.1805E-01  9.3896E-01  8.9414E-01  1.0366E+00  1.4142E+00  1.3494E-01  9.6263E-01  9.5252E-01
             1.0224E+00
 PARAMETER:  3.6626E-02  1.9837E-01 -2.3121E-01  3.7020E-02 -1.1891E-02  1.3598E-01  4.4657E-01 -1.9029E+00  6.1915E-02  5.1356E-02
             1.2218E-01
 GRADIENT:  -3.0033E+00 -7.8379E-01 -4.1170E-01 -1.2157E+00 -3.7118E-01 -1.0717E+00 -6.9988E-01  1.2758E-01 -6.9887E-01 -1.5190E-01
            -3.8594E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1594.63423099040        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      848             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3859E-01  1.1033E+00  7.2230E-01  9.4397E-01  8.9416E-01  1.0632E+00  1.4328E+00  1.3501E-01  9.5962E-01  9.5164E-01
             1.0224E+00
 PARAMETER:  3.6623E-02  1.9833E-01 -2.2532E-01  4.2340E-02 -1.1869E-02  1.6129E-01  4.5962E-01 -1.9024E+00  5.8785E-02  5.0430E-02
             1.2216E-01
 GRADIENT:  -7.4283E-01  4.3367E+00  8.4212E-01  2.8520E+00 -1.4757E+00  1.1014E+01  1.2855E+00  1.1431E-01 -6.8007E-01 -6.3025E-01
            -4.6904E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1594.63850480381        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1031
 NPARAMETR:  9.3861E-01  1.1032E+00  7.2215E-01  9.4650E-01  8.9421E-01  1.0633E+00  1.4325E+00  1.3487E-01  9.5960E-01  9.5722E-01
             1.0223E+00
 PARAMETER:  3.6573E-02  1.9823E-01 -2.2553E-01  4.5011E-02 -1.1819E-02  1.6140E-01  4.5939E-01 -1.9015E+00  5.8738E-02  5.6278E-02
             1.2210E-01
 GRADIENT:  -3.9616E+01 -2.2008E+00  8.1237E+05 -1.0239E-01 -9.4382E-01 -2.3294E-02 -1.9653E+00  1.1532E-01 -1.1106E+00  3.3938E-02
            -3.6464E-01
 NUMSIGDIG:         0.4         1.5         3.3         2.8         1.9         3.4         1.3         0.2         0.9         2.3
                    1.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1031
 NO. OF SIG. DIGITS IN FINAL EST.:  0.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8993E-02 -3.2036E-04 -6.1828E-03 -2.3047E-03 -1.5241E-02
 SE:             2.9774E-02  2.3956E-02  2.2300E-03  2.3057E-02  2.2307E-02
 N:                     100         100         100         100         100

 P VAL.:         5.2353E-01  9.8933E-01  5.5623E-03  9.2038E-01  4.9448E-01

 ETASHRINKSD(%)  2.5280E-01  1.9745E+01  9.2529E+01  2.2755E+01  2.5268E+01
 ETASHRINKVR(%)  5.0496E-01  3.5592E+01  9.9442E+01  4.0332E+01  4.4151E+01
 EBVSHRINKSD(%)  4.0528E-01  1.9628E+01  9.3419E+01  2.4032E+01  2.3846E+01
 EBVSHRINKVR(%)  8.0892E-01  3.5403E+01  9.9567E+01  4.2288E+01  4.2005E+01
 RELATIVEINF(%)  9.9002E+01  6.4024E+00  6.4738E-02  5.5662E+00  7.5240E+00
 EPSSHRINKSD(%)  4.3986E+01
 EPSSHRINKVR(%)  6.8624E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1594.6385048038087     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -859.48767824007052     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.18
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1594.639       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.39E-01  1.10E+00  7.22E-01  9.46E-01  8.94E-01  1.06E+00  1.43E+00  1.35E-01  9.60E-01  9.57E-01  1.02E+00
 


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
+        1.10E+03
 
 TH 2
+       -4.25E+00  2.99E+02
 
 TH 3
+        2.61E+08  1.44E+09  1.28E+09
 
 TH 4
+       -5.83E+00  2.45E+02 -1.98E+09  7.62E+02
 
 TH 5
+       -5.38E+02 -2.21E+02 -1.59E+09  3.35E+08  8.73E+02
 
 TH 6
+        3.91E+01 -7.58E-01  4.97E+07 -2.57E+00 -2.42E+00  1.73E+02
 
 TH 7
+        8.39E-01  2.51E+01 -1.06E+07 -1.33E+01  4.10E+00  1.31E-01  4.34E+01
 
 TH 8
+       -2.56E+09 -2.30E+07 -1.55E+09 -3.24E+09 -1.42E+08 -3.05E-01 -2.48E+04  1.23E+09
 
 TH 9
+        1.20E+09 -1.68E+01 -3.74E+09  6.88E+09  9.98E+00 -1.56E-01  1.38E+01  1.38E+00  8.56E+01
 
 TH10
+       -7.99E-02 -1.08E+01 -1.71E+09 -1.82E+01 -5.53E+01 -2.49E-01  9.38E+00 -2.65E+08  1.30E+01  7.59E+01
 
 TH11
+       -3.67E+08 -1.20E+01 -5.45E+09  1.48E+03 -6.82E+09  1.14E+00  4.24E+00  2.06E+09  7.45E+00  1.57E+01  3.56E+09
 
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
 #CPUT: Total CPU Time in Seconds,       23.099
Stop Time:
Sat Sep 18 12:18:41 CDT 2021
