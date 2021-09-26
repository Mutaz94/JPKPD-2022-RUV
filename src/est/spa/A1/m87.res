Sat Sep 25 08:18:40 CDT 2021
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
$DATA ../../../../data/spa/A1/dat87.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -999.162440957190        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.2384E+01 -6.1779E+01 -2.9435E+01 -5.1719E+01  1.6600E+02 -8.9912E+00 -4.7899E+01 -1.1382E+01 -7.5931E+01 -7.4042E+01
            -1.1368E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1368.56775162211        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0849E+00  1.1106E+00  1.1581E+00  1.0684E+00  1.0455E+00  9.1264E-01  8.6547E-01  9.6834E-01  8.1663E-01  7.3964E-01
             4.0987E+00
 PARAMETER:  1.8151E-01  2.0494E-01  2.4681E-01  1.6616E-01  1.4449E-01  8.5855E-03 -4.4480E-02  6.7833E-02 -1.0256E-01 -2.0160E-01
             1.5107E+00
 GRADIENT:   7.9756E+01 -5.0104E+00 -1.5942E+01  1.0252E+01  9.2563E+00 -2.2413E+01  7.0967E+00  4.2768E+00  7.1261E+00  1.1715E+01
             1.8558E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1401.51892869795        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0501E+00  1.2415E+00  9.7992E-01  9.4764E-01  9.4678E-01  9.7505E-01  5.7539E-01  9.2709E-01  1.0854E+00  4.1696E-01
             3.2723E+00
 PARAMETER:  1.4887E-01  3.1631E-01  7.9713E-02  4.6217E-02  4.5315E-02  7.4734E-02 -4.5270E-01  2.4293E-02  1.8195E-01 -7.7478E-01
             1.2855E+00
 GRADIENT:   3.3889E+01  3.3574E+01  1.1799E+01  1.9863E+01 -4.5134E+01 -2.7667E+00 -3.3138E+00  2.8141E+00  6.3477E+00  1.6225E+00
             6.2808E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1407.74851327267        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0264E+00  8.6704E-01  1.0010E+00  1.1619E+00  8.5901E-01  9.8090E-01  8.8612E-01  4.8116E-01  8.7674E-01  3.9454E-01
             2.9747E+00
 PARAMETER:  1.2603E-01 -4.2667E-02  1.0103E-01  2.5002E-01 -5.1973E-02  8.0716E-02 -2.0906E-02 -6.3155E-01 -3.1547E-02 -8.3003E-01
             1.1901E+00
 GRADIENT:   6.4963E-02 -7.6457E+00 -1.1942E+01  9.4749E-01  1.7017E+01 -1.1179E+00 -1.4056E+00  9.7826E-01  3.3692E-02  1.5635E-01
             2.6220E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1408.66289256903        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  1.0237E+00  7.8161E-01  1.1471E+00  1.2198E+00  8.8205E-01  9.7752E-01  9.5012E-01  8.7403E-02  8.3599E-01  4.0886E-01
             3.0004E+00
 PARAMETER:  1.2343E-01 -1.4639E-01  2.3720E-01  2.9867E-01 -2.5511E-02  7.7267E-02  4.8828E-02 -2.3372E+00 -7.9137E-02 -7.9438E-01
             1.1988E+00
 GRADIENT:  -1.1166E+00 -1.6252E+00 -2.8909E+00  2.7735E-01  4.7478E+00 -3.1301E-01 -2.8617E-01  3.9504E-02 -4.4797E-02  1.2628E-01
             1.6198E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1408.71868437852        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.0237E+00  7.5753E-01  1.1727E+00  1.2345E+00  8.8032E-01  9.7748E-01  9.8577E-01  4.7072E-02  8.2544E-01  3.9370E-01
             2.9995E+00
 PARAMETER:  1.2339E-01 -1.7769E-01  2.5933E-01  3.1069E-01 -2.7474E-02  7.7220E-02  8.5670E-02 -2.9561E+00 -9.1837E-02 -8.3216E-01
             1.1984E+00
 GRADIENT:   3.9407E-02 -3.9089E-01 -1.2539E-01 -3.1159E-01  2.7467E-01 -1.2723E-02  5.0702E-03  1.1100E-02  3.1589E-02 -6.4306E-03
             1.4729E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1408.78247854127        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      536
 NPARAMETR:  1.0266E+00  6.6823E-01  1.2206E+00  1.2957E+00  8.6985E-01  9.7801E-01  1.0603E+00  1.0000E-02  7.9657E-01  3.6734E-01
             3.0159E+00
 PARAMETER:  1.2624E-01 -3.0312E-01  2.9937E-01  3.5904E-01 -3.9430E-02  7.7762E-02  1.5856E-01 -5.4007E+00 -1.2744E-01 -9.0147E-01
             1.2039E+00
 GRADIENT:   1.9537E+00  1.4070E+00  4.6844E-01  2.9436E+00 -1.1078E+00  1.4139E-02 -9.6804E-02  0.0000E+00 -2.6135E-01 -1.4998E-01
             5.8669E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1408.84373221981        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  1.0247E+00  5.4272E-01  1.2129E+00  1.3694E+00  8.2938E-01  9.7887E-01  1.1706E+00  1.0000E-02  7.7406E-01  4.6169E-01
             2.9873E+00
 PARAMETER:  1.2443E-01 -5.1115E-01  2.9305E-01  4.1435E-01 -8.7080E-02  7.8647E-02  2.5754E-01 -7.7603E+00 -1.5610E-01 -6.7287E-01
             1.1944E+00
 GRADIENT:   6.7305E-01 -4.2568E-02 -5.0027E-01  4.9843E-01  8.3417E-01  6.9686E-02 -1.3145E-01  0.0000E+00 -1.2798E-01  2.6609E-02
             1.8758E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1408.85379897622        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      888
 NPARAMETR:  1.0240E+00  5.1025E-01  1.1967E+00  1.3874E+00  8.1267E-01  9.7929E-01  1.2295E+00  1.0000E-02  7.6964E-01  4.9374E-01
             2.9657E+00
 PARAMETER:  1.2374E-01 -5.7286E-01  2.7960E-01  4.2747E-01 -1.0744E-01  7.9068E-02  3.0661E-01 -8.5825E+00 -1.6183E-01 -6.0574E-01
             1.1871E+00
 GRADIENT:   1.7724E-02  2.3518E-01  1.1474E-01  5.3763E-01 -3.4559E-01  1.3227E-02  1.4427E-02  0.0000E+00  2.5439E-03  3.8542E-02
             3.1998E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1408.85517223019        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  1.0236E+00  4.8432E-01  1.2015E+00  1.4028E+00  8.0738E-01  9.7912E-01  1.2639E+00  1.0000E-02  7.6436E-01  4.9912E-01
             2.9619E+00
 PARAMETER:  1.2333E-01 -6.2501E-01  2.8355E-01  4.3844E-01 -1.1396E-01  7.8900E-02  3.3422E-01 -9.4615E+00 -1.6872E-01 -5.9491E-01
             1.1858E+00
 GRADIENT:  -5.7119E-03  1.0568E-02  4.6254E-03  5.3617E-02 -4.1819E-03 -2.5160E-03  2.4695E-04  0.0000E+00  9.7503E-03 -2.9476E-04
            -5.6627E-03

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1408.85517223019        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1085
 NPARAMETR:  1.0236E+00  4.8432E-01  1.2015E+00  1.4028E+00  8.0738E-01  9.7912E-01  1.2639E+00  1.0000E-02  7.6436E-01  4.9912E-01
             2.9619E+00
 PARAMETER:  1.2333E-01 -6.2501E-01  2.8355E-01  4.3844E-01 -1.1396E-01  7.8900E-02  3.3422E-01 -9.4615E+00 -1.6872E-01 -5.9491E-01
             1.1858E+00
 GRADIENT:  -5.7119E-03  1.0568E-02  4.6254E-03  5.3617E-02 -4.1819E-03 -2.5160E-03  2.4695E-04  0.0000E+00  9.7503E-03 -2.9476E-04
            -5.6627E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1085
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6306E-04 -3.7928E-03  1.0460E-04 -1.1349E-02 -8.7432E-03
 SE:             2.8916E-02  1.0045E-02  1.1249E-04  2.3603E-02  1.1191E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9550E-01  7.0574E-01  3.5244E-01  6.3063E-01  4.3463E-01

 ETASHRINKSD(%)  3.1277E+00  6.6348E+01  9.9623E+01  2.0925E+01  6.2510E+01
 ETASHRINKVR(%)  6.1575E+00  8.8675E+01  9.9999E+01  3.7472E+01  8.5945E+01
 EBVSHRINKSD(%)  3.0240E+00  6.6736E+01  9.9582E+01  2.0608E+01  6.2917E+01
 EBVSHRINKVR(%)  5.9566E+00  8.8935E+01  9.9998E+01  3.6970E+01  8.6249E+01
 RELATIVEINF(%)  8.6047E+01  9.1270E-02  7.6559E-05  7.3516E-01  4.7259E-01
 EPSSHRINKSD(%)  2.4541E+01
 EPSSHRINKVR(%)  4.3060E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1408.8551722301931     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -673.70434566645497     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.78
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1408.855       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  4.84E-01  1.20E+00  1.40E+00  8.07E-01  9.79E-01  1.26E+00  1.00E-02  7.64E-01  4.99E-01  2.96E+00
 


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
+        1.04E+03
 
 TH 2
+       -5.88E+01  3.25E+02
 
 TH 3
+        9.59E+00  9.00E+01  1.25E+02
 
 TH 4
+       -7.10E+01  3.91E+02  4.68E+00  6.09E+02
 
 TH 5
+        8.75E+00 -3.11E+02 -3.18E+02 -1.27E+02  8.70E+02
 
 TH 6
+       -1.82E+00 -8.68E+00  5.38E+00 -1.86E+01 -1.05E+01  1.79E+02
 
 TH 7
+        4.42E-01  1.15E+00  9.58E-01 -2.71E+00  1.47E+00  7.42E-01  3.14E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.30E+00 -1.89E+01  1.16E+01 -6.95E+00  3.54E+00  1.62E+00  9.57E+00  0.00E+00  1.33E+02
 
 TH10
+       -7.32E+00 -6.01E+00 -4.19E+00 -1.05E+01  1.12E+01 -1.09E+00  1.18E+00  0.00E+00  1.52E+00  1.19E+01
 
 TH11
+       -1.43E+01 -6.71E+00 -1.26E+00 -1.08E+01 -6.16E+00  3.34E+00  1.59E+00  0.00E+00  1.31E+01  1.54E+01  4.33E+01
 
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
 #CPUT: Total CPU Time in Seconds,       17.866
Stop Time:
Sat Sep 25 08:18:59 CDT 2021
