Sat Sep 18 10:03:25 CDT 2021
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
$DATA ../../../../data/spa/A2/dat76.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1017.60567858861        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4486E+02  6.7329E+01  7.9370E+01 -9.6557E+00  5.2575E+01 -3.9791E+01 -3.5358E+01 -2.0793E+01 -6.9057E+01 -9.9824E+01
            -9.8906E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1356.13629600365        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.4041E-01  8.4997E-01  9.1395E-01  1.0902E+00  8.8677E-01  1.0657E+00  1.0475E+00  9.6813E-01  1.1912E+00  1.0410E+00
             2.1724E+00
 PARAMETER:  3.8559E-02 -6.2556E-02  1.0023E-02  1.8638E-01 -2.0169E-02  1.6362E-01  1.4639E-01  6.7615E-02  2.7499E-01  1.4019E-01
             8.7583E-01
 GRADIENT:  -4.3069E+01  3.0293E+00  5.2341E-01  4.3532E-01  1.7599E+01 -2.1710E+00 -1.1794E-01  3.8248E+00  8.9863E+00 -4.2196E+00
            -3.7910E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1360.91210750552        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.4736E-01  5.3604E-01  4.0707E-01  1.2996E+00  4.2037E-01  1.1235E+00  1.1523E+00  4.0701E-01  1.0866E+00  5.4507E-01
             2.1108E+00
 PARAMETER:  4.5922E-02 -5.2354E-01 -7.9876E-01  3.6206E-01 -7.6663E-01  2.1642E-01  2.4173E-01 -7.9893E-01  1.8302E-01 -5.0685E-01
             8.4709E-01
 GRADIENT:  -4.1624E+01  7.4325E+01  5.6606E+01  1.7382E+02 -9.1342E+01  1.0712E+01 -1.1185E+01 -3.6652E+00 -3.2194E+00 -1.5257E+01
            -3.5517E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1376.44767072733        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8419E-01  4.9612E-01  3.0026E-01  1.1409E+00  3.5165E-01  1.0534E+00  1.1738E+00  1.5645E-01  1.0055E+00  6.3539E-01
             2.1893E+00
 PARAMETER:  8.4063E-02 -6.0094E-01 -1.1031E+00  2.3180E-01 -9.4511E-01  1.5205E-01  2.6027E-01 -1.7550E+00  1.0545E-01 -3.5352E-01
             8.8360E-01
 GRADIENT:   2.4223E+01  2.5533E+01  7.9887E+00  3.2911E+01 -1.2650E+01 -1.2711E+01 -3.4507E+00 -1.5475E-02 -1.2532E+01  6.1280E+00
             1.2578E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1378.83974049817        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.7293E-01  3.5927E-01  2.1318E-01  1.0955E+00  2.6015E-01  1.0953E+00  1.2358E+00  5.0532E-02  1.0889E+00  6.4195E-01
             2.0241E+00
 PARAMETER:  7.2559E-02 -9.2369E-01 -1.4456E+00  1.9121E-01 -1.2465E+00  1.9105E-01  3.1171E-01 -2.8851E+00  1.8520E-01 -3.4325E-01
             8.0512E-01
 GRADIENT:   4.9204E+00  1.8374E+01  4.9906E+00  3.2722E+01 -2.1618E+01 -1.4018E+00 -2.2683E-01 -2.1548E-02 -8.3682E+00  1.2228E-01
             8.3615E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1380.95144991273        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      481
 NPARAMETR:  9.7249E-01  2.7173E-01  2.8436E-01  1.1898E+00  2.9115E-01  1.0940E+00  1.6215E+00  2.5623E-02  1.0398E+00  7.0347E-01
             2.0709E+00
 PARAMETER:  7.2109E-02 -1.2030E+00 -1.1575E+00  2.7375E-01 -1.1339E+00  1.8988E-01  5.8333E-01 -3.5643E+00  1.3899E-01 -2.5174E-01
             8.2797E-01
 GRADIENT:  -2.9449E+00  6.0942E+00  1.6894E+01  2.4462E+00 -2.9036E+01  2.0048E+00  2.2956E-01 -3.9290E-03 -4.2140E-01 -7.6899E-01
            -1.5606E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1382.02355219088        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  9.6949E-01  1.6857E-01  3.5321E-01  1.2931E+00  3.2441E-01  1.0728E+00  2.2276E+00  1.0000E-02  9.8272E-01  7.7506E-01
             2.1072E+00
 PARAMETER:  6.9017E-02 -1.6804E+00 -9.4070E-01  3.5701E-01 -1.0258E+00  1.7031E-01  9.0093E-01 -5.9570E+00  8.2570E-02 -1.5481E-01
             8.4536E-01
 GRADIENT:   4.0884E+00  3.2775E+00  1.5444E+01  9.6400E+00 -2.1588E+01 -1.8759E+00 -1.3494E+00  0.0000E+00 -2.0782E+00  4.7491E+00
            -5.7498E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1383.39730895302        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.5922E-01  4.7173E-02  3.4690E-01  1.3298E+00  3.1088E-01  1.0698E+00  5.1160E+00  1.0000E-02  9.7358E-01  7.7120E-01
             2.1033E+00
 PARAMETER:  5.8362E-02 -2.9539E+00 -9.5873E-01  3.8506E-01 -1.0683E+00  1.6742E-01  1.7324E+00 -1.3482E+01  7.3226E-02 -1.5980E-01
             8.4351E-01
 GRADIENT:  -9.0515E-01  2.8953E-01  2.2661E+00  9.2059E+00 -5.5428E+00 -1.3924E+00 -1.1354E-01  0.0000E+00 -1.2076E-01  3.5966E-01
            -3.3879E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1383.49331282851        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1009
 NPARAMETR:  9.5686E-01  1.9697E-02  3.4428E-01  1.3314E+00  3.0717E-01  1.0731E+00  7.9811E+00  1.0000E-02  9.7342E-01  7.7554E-01
             2.1033E+00
 PARAMETER:  5.5900E-02 -3.8273E+00 -9.6629E-01  3.8624E-01 -1.0804E+00  1.7055E-01  2.1771E+00 -1.8776E+01  7.3057E-02 -1.5419E-01
             8.4350E-01
 GRADIENT:  -1.7366E+00  8.9869E-02  1.0499E+00  8.0002E-01 -2.1614E+00  1.7315E-01  6.6662E-02  0.0000E+00  4.8995E-01  4.4911E-01
             4.8908E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1383.51492898624        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1184
 NPARAMETR:  9.5710E-01  1.0000E-02  3.4780E-01  1.3376E+00  3.0886E-01  1.0721E+00  1.1121E+01  1.0000E-02  9.6888E-01  7.7588E-01
             2.1037E+00
 PARAMETER:  5.6153E-02 -4.5252E+00 -9.5614E-01  3.9086E-01 -1.0749E+00  1.6960E-01  2.5088E+00 -2.3003E+01  6.8383E-02 -1.5376E-01
             8.4371E-01
 GRADIENT:   1.9163E-01  0.0000E+00 -5.0081E-02  7.9988E-03  1.2443E-01 -3.2150E-02  1.3773E-03  0.0000E+00 -5.3944E-02 -4.8167E-02
            -6.0591E-02

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1383.51496787733        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  9.5699E-01  1.0000E-02  3.4742E-01  1.3373E+00  3.0859E-01  1.0722E+00  1.1110E+01  1.0000E-02  9.6923E-01  7.7604E-01
             2.1037E+00
 PARAMETER:  5.6034E-02 -4.5227E+00 -9.5721E-01  3.9063E-01 -1.0757E+00  1.6972E-01  2.5078E+00 -2.2988E+01  6.8751E-02 -1.5355E-01
             8.4369E-01
 GRADIENT:  -4.8346E-02  0.0000E+00  7.8459E-02  9.8574E-02 -1.6577E-01  9.5058E-03 -3.6019E-03  0.0000E+00 -3.1716E-03  1.7179E-02
             1.3419E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1241
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.6079E-04  1.8475E-03  1.6957E-05 -7.4181E-03 -3.4973E-03
 SE:             2.9353E-02  1.8802E-03  2.5427E-04  2.7859E-02  2.4561E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9019E-01  3.2581E-01  9.4683E-01  7.9003E-01  8.8677E-01

 ETASHRINKSD(%)  1.6626E+00  9.3701E+01  9.9148E+01  6.6683E+00  1.7716E+01
 ETASHRINKVR(%)  3.2977E+00  9.9603E+01  9.9993E+01  1.2892E+01  3.2294E+01
 EBVSHRINKSD(%)  1.5757E+00  9.4906E+01  9.9189E+01  5.8188E+00  1.7086E+01
 EBVSHRINKVR(%)  3.1266E+00  9.9740E+01  9.9993E+01  1.1299E+01  3.1252E+01
 RELATIVEINF(%)  8.4714E+01  3.5885E-02  3.0299E-04  1.7786E+01  2.7954E+00
 EPSSHRINKSD(%)  3.7554E+01
 EPSSHRINKVR(%)  6.1005E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1383.5149678773334     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -648.36414131359527     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1383.515       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.57E-01  1.00E-02  3.47E-01  1.34E+00  3.09E-01  1.07E+00  1.11E+01  1.00E-02  9.69E-01  7.76E-01  2.10E+00
 


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
+        1.02E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.87E+01  0.00E+00  6.07E+03
 
 TH 4
+       -1.45E+01  0.00E+00 -4.56E+02  5.88E+02
 
 TH 5
+        8.87E+01  0.00E+00 -8.38E+03 -1.72E+02  1.34E+04
 
 TH 6
+       -6.93E+00  0.00E+00  1.02E+01 -7.78E+00 -2.22E+00  1.60E+02
 
 TH 7
+        1.18E-02  0.00E+00 -1.16E-01 -7.10E-02  2.84E-01  1.02E-02  1.80E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -8.87E-01  0.00E+00  8.40E+01 -1.33E+01  1.22E+01  3.86E+00  6.80E-02  0.00E+00  1.74E+02
 
 TH10
+       -2.02E+01  0.00E+00 -4.77E+01  8.96E+00 -3.83E+00  8.34E+00 -1.11E-01  0.00E+00  1.03E+01  1.51E+02
 
 TH11
+       -1.32E+01  0.00E+00 -3.36E+01 -6.44E+00  9.42E+00  2.97E+00 -4.49E-02  0.00E+00  6.32E+00  2.49E+01  5.67E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.330
Stop Time:
Sat Sep 18 10:03:47 CDT 2021
