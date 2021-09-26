Sat Sep 25 00:35:23 CDT 2021
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
$DATA ../../../../data/int/SL1/dat81.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2539.50261826662        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.1429E+00 -3.7648E+01  1.3933E+02  3.6230E+00  1.0871E+02 -3.8626E+01 -5.7454E+01 -3.3026E+02 -8.8803E+01 -3.4092E+01
            -2.2141E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3356.08831666807        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.9220E-01  1.1357E+00  9.8410E-01  9.5682E-01  1.0298E+00  1.0797E+00  9.9928E-01  1.1611E+00  1.0592E+00  1.0507E+00
             1.7271E+00
 PARAMETER:  9.2170E-02  2.2723E-01  8.3970E-02  5.5860E-02  1.2934E-01  1.7670E-01  9.9284E-02  2.4938E-01  1.5749E-01  1.4948E-01
             6.4645E-01
 GRADIENT:  -4.6516E+01  2.4347E+00 -1.5039E+01  3.2020E+01  1.4001E+01 -8.7409E+00  6.9086E+00  1.2367E+00  8.2575E+00 -1.7231E+00
             1.4741E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3358.67130179488        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0129E+00  1.2649E+00  1.0933E+00  8.7890E-01  1.1713E+00  1.0518E+00  9.2884E-01  1.3247E+00  1.1485E+00  1.1594E+00
             1.7206E+00
 PARAMETER:  1.1284E-01  3.3497E-01  1.8923E-01 -2.9080E-02  2.5813E-01  1.5047E-01  2.6185E-02  3.8118E-01  2.3848E-01  2.4786E-01
             6.4270E-01
 GRADIENT:  -2.5609E+01 -1.3300E+01 -1.3074E+01  1.9435E+01  2.6433E+01 -2.3891E+01  9.9045E+00  7.5133E-01  1.4231E+01  2.1968E+00
             1.3593E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3359.91585064576        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      386
 NPARAMETR:  1.0179E+00  1.2651E+00  1.0933E+00  8.5971E-01  1.1712E+00  1.1274E+00  8.8874E-01  1.3249E+00  1.1487E+00  1.1558E+00
             1.7202E+00
 PARAMETER:  1.1778E-01  3.3514E-01  1.8916E-01 -5.1157E-02  2.5800E-01  2.1988E-01 -1.7952E-02  3.8133E-01  2.3860E-01  2.4478E-01
             6.4244E-01
 GRADIENT:   6.7695E+00 -1.8145E+01 -9.8932E+00 -4.4699E+00  2.4126E+01  1.2103E+01  5.1580E+00  7.4153E-01  1.3263E+01  5.0167E-01
             1.3585E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3359.96568692538        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  1.0179E+00  1.2651E+00  1.0933E+00  8.6482E-01  1.1712E+00  1.1137E+00  8.8874E-01  1.3249E+00  1.1487E+00  1.1557E+00
             1.7202E+00
 PARAMETER:  1.1778E-01  3.3514E-01  1.8916E-01 -4.5238E-02  2.5800E-01  2.0771E-01 -1.7952E-02  3.8133E-01  2.3860E-01  2.4468E-01
             6.4244E-01
 GRADIENT:  -1.3712E+01 -2.7722E+01 -1.0715E+01  1.7280E-02  2.1709E+01 -1.2444E-02  4.9135E+00  6.5529E-01  1.2991E+01 -1.1260E-02
             1.3461E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3362.78218507269        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  1.0224E+00  1.2711E+00  1.1318E+00  8.6536E-01  1.1570E+00  1.1085E+00  8.5191E-01  1.3543E+00  1.0809E+00  1.1604E+00
             1.6912E+00
 PARAMETER:  1.2218E-01  3.3985E-01  2.2381E-01 -4.4609E-02  2.4579E-01  2.0305E-01 -6.0277E-02  4.0327E-01  1.7781E-01  2.4878E-01
             6.2542E-01
 GRADIENT:   1.5914E+01  5.5735E+00  1.7278E+00 -3.0577E+00 -2.2353E+00  5.2835E+00 -8.1968E-01 -4.3398E-01 -1.2439E-01 -1.3934E-01
             9.8017E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3362.79818926791        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      838
 NPARAMETR:  1.0253E+00  1.2710E+00  1.1246E+00  8.6906E-01  1.1570E+00  1.1129E+00  8.6098E-01  1.3541E+00  1.0842E+00  1.1637E+00
             1.6914E+00
 PARAMETER:  1.2503E-01  3.3977E-01  2.1743E-01 -4.0339E-02  2.4585E-01  2.0699E-01 -4.9685E-02  4.0317E-01  1.8081E-01  2.5156E-01
             6.2557E-01
 GRADIENT:  -5.4846E-01 -5.8235E+00 -4.9854E-01  4.3790E-01 -3.0940E+00 -1.8447E-01  2.7718E-01 -2.3930E-01  5.1766E-01  3.2577E-01
             9.7857E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3362.85946544586        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1011
 NPARAMETR:  1.0255E+00  1.2714E+00  1.1247E+00  8.6807E-01  1.1567E+00  1.1135E+00  8.6092E-01  1.3547E+00  1.0840E+00  1.1573E+00
             1.6904E+00
 PARAMETER:  1.2519E-01  3.4011E-01  2.1749E-01 -4.1481E-02  2.4560E-01  2.0748E-01 -4.9753E-02  4.0358E-01  1.8069E-01  2.4609E-01
             6.2494E-01
 GRADIENT:   2.1818E+01  8.6214E+00 -1.4001E-02  1.8503E+00  3.2660E-01  7.1027E+00  3.2047E-01 -1.5874E-01  8.4086E-01 -3.1901E-01
             9.7527E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3362.86204017668        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1192
 NPARAMETR:  1.0255E+00  1.2720E+00  1.1249E+00  8.6864E-01  1.1572E+00  1.1135E+00  8.6067E-01  1.3558E+00  1.0838E+00  1.1620E+00
             1.6904E+00
 PARAMETER:  1.2520E-01  3.4058E-01  2.1772E-01 -4.0828E-02  2.4599E-01  2.0747E-01 -5.0049E-02  4.0438E-01  1.8043E-01  2.5014E-01
             6.2498E-01
 GRADIENT:  -2.2027E-01 -4.9976E+00 -4.1276E-01  4.2841E-01 -3.6056E+00  3.2203E-04  2.0442E-01 -2.5235E-01  3.2328E-01 -4.0179E-02
             9.6588E+01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -3362.86362644117        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  1.0255E+00  1.2722E+00  1.1252E+00  8.6860E-01  1.1574E+00  1.1135E+00  8.6048E-01  1.3562E+00  1.0836E+00  1.1620E+00
             1.6904E+00
 PARAMETER:  1.2521E-01  3.4076E-01  2.1782E-01 -4.0845E-02  2.4614E-01  2.0747E-01 -5.0162E-02  4.0468E-01  1.8035E-01  2.5015E-01
             6.2498E-01
 GRADIENT:  -1.6882E-01  3.1312E+05 -4.1799E-01  5.4262E-01 -2.1675E+05 -7.2673E-03  1.9653E-01  2.6360E+05  2.8845E-01 -4.4183E-02
            -1.7077E+05
 NUMSIGDIG:         2.8         3.3         1.7         2.1         3.3         4.1         1.5         3.3         1.8         2.8
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1259
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.5462E-04 -2.8288E-02 -1.9762E-02  1.7575E-02 -2.3370E-02
 SE:             2.9785E-02  2.0763E-02  1.5023E-02  2.5558E-02  2.4959E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7711E-01  1.7307E-01  1.8838E-01  4.9168E-01  3.4910E-01

 ETASHRINKSD(%)  2.1533E-01  3.0441E+01  4.9670E+01  1.4376E+01  1.6384E+01
 ETASHRINKVR(%)  4.3021E-01  5.1616E+01  7.4669E+01  2.6685E+01  3.0084E+01
 EBVSHRINKSD(%)  5.7766E-01  3.1064E+01  5.2419E+01  1.5797E+01  1.5264E+01
 EBVSHRINKVR(%)  1.1520E+00  5.2478E+01  7.7361E+01  2.9098E+01  2.8199E+01
 RELATIVEINF(%)  9.8837E+01  1.4156E+01  1.6427E+01  2.6242E+01  2.5320E+01
 EPSSHRINKSD(%)  2.2470E+01
 EPSSHRINKVR(%)  3.9891E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3362.8636264411707     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1708.7742666727599     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.28
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3362.864       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.27E+00  1.13E+00  8.69E-01  1.16E+00  1.11E+00  8.61E-01  1.36E+00  1.08E+00  1.16E+00  1.69E+00
 


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
+        1.62E+09
 
 TH 2
+        9.88E+02  1.42E+08
 
 TH 3
+       -8.48E+08  1.95E+04  4.44E+08
 
 TH 4
+       -9.19E+00 -1.23E+04 -4.61E+01  8.78E+02
 
 TH 5
+       -1.51E+03 -5.95E+02 -2.97E+04  1.94E+04  3.29E+08
 
 TH 6
+        3.42E+00  1.71E+03  5.13E-01 -3.19E+00 -2.60E+03  1.57E+02
 
 TH 7
+       -1.87E+00 -2.02E+03 -1.26E+09 -9.29E+00  3.12E+03 -5.90E-01  7.24E+01
 
 TH 8
+        7.84E+02  1.12E+08  1.54E+04 -1.00E+04 -3.09E+02  1.35E+03 -1.62E+03  8.85E+07
 
 TH 9
+        1.06E+09 -4.50E+03 -5.57E+08  3.62E+01  6.83E+03 -1.03E+00  1.59E+09 -3.54E+03  6.98E+08
 
 TH10
+       -2.92E-01 -1.56E+01  5.73E+00  5.60E+00 -1.63E+01  2.09E-01  1.81E+01 -4.80E-01  3.07E+00  7.50E+01
 
 TH11
+       -4.16E+02  2.26E+03 -7.99E+03  5.20E+03  1.52E+02 -7.00E+02  8.52E+02  2.30E+04  1.86E+03  7.87E+00  2.39E+07
 
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
 #CPUT: Total CPU Time in Seconds,       51.105
Stop Time:
Sat Sep 25 00:36:15 CDT 2021
