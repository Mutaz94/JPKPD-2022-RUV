Thu Sep 30 05:37:51 CDT 2021
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
$DATA ../../../../data/spa2/A2/dat30.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -951.899611452882        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4587E+02  1.5264E+02  1.7254E+02  1.0721E+02  2.1011E+02  5.7468E+01 -1.3943E+02 -8.8042E+01 -2.8754E+01 -9.0907E+01
            -2.6856E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1882.38604158014        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  1.0425E+00  8.4896E-01  7.6243E-01  1.0926E+00  7.6281E-01  8.8206E-01  1.2565E+00  9.5588E-01  1.1773E+00  1.0405E+00
             2.0688E+00
 PARAMETER:  1.4162E-01 -6.3745E-02 -1.7125E-01  1.8859E-01 -1.7075E-01 -2.5500E-02  3.2833E-01  5.4874E-02  2.6319E-01  1.3973E-01
             8.2695E-01
 GRADIENT:   2.3254E+02  4.6124E+01  2.9463E+01  7.0936E+01  3.5946E+01 -2.7101E+01 -5.5379E+00  7.9891E+00  2.9726E+01  7.4680E+00
            -2.1221E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1900.06351089971        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      184
 NPARAMETR:  1.0682E+00  5.0247E-01  3.8280E-01  1.2474E+00  4.2449E-01  9.0824E-01  1.3166E+00  7.5368E-01  1.1397E+00  7.7165E-01
             1.9014E+00
 PARAMETER:  1.6594E-01 -5.8822E-01 -8.6026E-01  3.2109E-01 -7.5686E-01  3.7530E-03  3.7505E-01 -1.8279E-01  2.3078E-01 -1.5922E-01
             7.4260E-01
 GRADIENT:   1.3718E+02  2.6869E+01 -4.5365E+01  1.2852E+02  8.0063E+01 -3.6375E+01  5.4980E+00  4.0933E+00 -1.0294E+01  1.6388E+01
            -2.6299E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1924.57163491571        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      362
 NPARAMETR:  1.0503E+00  4.7419E-01  3.6494E-01  1.2314E+00  4.0260E-01  9.5418E-01  1.1996E+00  2.8926E-01  1.2717E+00  7.9854E-01
             2.0642E+00
 PARAMETER:  1.4907E-01 -6.4615E-01 -9.0801E-01  3.0817E-01 -8.0980E-01  5.3095E-02  2.8201E-01 -1.1404E+00  3.4036E-01 -1.2497E-01
             8.2475E-01
 GRADIENT:   8.1842E+01  1.0125E+01 -3.4464E+01  8.6283E+01  6.5569E+01 -1.0870E+01  3.0314E+00  4.7342E-02  2.2140E+01  8.4080E+00
            -1.4643E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1942.40454913065        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  1.0178E+00  3.8897E-01  2.9502E-01  1.1432E+00  3.2308E-01  9.7920E-01  1.0650E+00  3.9701E-01  1.1834E+00  7.1763E-01
             2.2870E+00
 PARAMETER:  1.1761E-01 -8.4425E-01 -1.1207E+00  2.3379E-01 -1.0299E+00  7.8984E-02  1.6299E-01 -8.2380E-01  2.6842E-01 -2.3180E-01
             9.2724E-01
 GRADIENT:  -2.5320E+00  6.5276E+00  6.1496E+00  1.2210E+01 -1.5002E+01  7.8672E-01 -4.9467E-01 -8.5892E-01 -4.1087E-01 -3.2218E+00
             2.3141E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1943.45029919328        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      714
 NPARAMETR:  1.0204E+00  3.6570E-01  2.6726E-01  1.1111E+00  2.9834E-01  9.7512E-01  1.1082E+00  1.1881E+00  1.1916E+00  6.0158E-01
             2.2322E+00
 PARAMETER:  1.2022E-01 -9.0595E-01 -1.2195E+00  2.0538E-01 -1.1095E+00  7.4806E-02  2.0275E-01  2.7232E-01  2.7533E-01 -4.0819E-01
             9.0300E-01
 GRADIENT:   2.4553E+00  1.7950E+01  1.2689E+01 -9.6072E-01 -1.8621E+01 -2.2900E+00  1.3415E+00  6.4934E+00 -1.7552E+00  1.3480E+01
             2.8508E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1945.94739662837        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      898             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0191E+00  3.4285E-01  2.5012E-01  1.1062E+00  2.8759E-01  9.8151E-01  1.1152E+00  1.2157E+00  1.2234E+00  4.0502E-01
             2.1932E+00
 PARAMETER:  1.1888E-01 -9.7047E-01 -1.2858E+00  2.0096E-01 -1.1462E+00  8.1337E-02  2.0901E-01  2.9532E-01  3.0166E-01 -8.0382E-01
             8.8535E-01
 GRADIENT:   8.7000E+01  2.1470E+01  4.4510E+01  4.9389E+01  2.3107E+02  6.2463E+00 -8.8780E+00 -9.7371E+00  1.0937E+01 -1.0603E+00
             3.6227E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1946.54759554406        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1053
 NPARAMETR:  1.0191E+00  3.4284E-01  2.5013E-01  1.1062E+00  2.8758E-01  9.8286E-01  1.2454E+00  1.2157E+00  1.2274E+00  4.1357E-01
             2.1931E+00
 PARAMETER:  1.1892E-01 -9.7050E-01 -1.2858E+00  2.0095E-01 -1.1463E+00  8.2710E-02  3.1944E-01  2.9531E-01  3.0493E-01 -7.8294E-01
             8.8531E-01
 GRADIENT:   1.5800E-01 -7.3372E+00 -2.3552E+00  5.8800E+00  2.8868E+00 -6.1355E-02  4.7843E-02 -7.2945E+00  4.7196E-02  5.1472E-02
             3.0352E+01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1946.54759554406        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:     1103
 NPARAMETR:  1.0191E+00  3.4243E-01  2.5052E-01  1.1040E+00  2.8798E-01  9.8291E-01  1.2452E+00  1.2152E+00  1.2275E+00  4.1343E-01
             2.1954E+00
 PARAMETER:  1.1892E-01 -9.7050E-01 -1.2858E+00  2.0095E-01 -1.1463E+00  8.2710E-02  3.1944E-01  2.9531E-01  3.0493E-01 -7.8294E-01
             8.8531E-01
 GRADIENT:  -1.1584E-01  3.7279E+03 -2.8096E+03  5.4267E+00 -3.1588E+03 -7.2983E-02  4.2459E-02  1.2275E+04 -6.3491E-02  3.4404E-02
            -4.1316E+03
 NUMSIGDIG:         3.4         2.3         2.3         1.4         2.3         2.7         2.9         2.3         3.1         2.8
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1103
 NO. OF SIG. DIGITS IN FINAL EST.:  1.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4170E-04  1.5199E-02 -7.0154E-03 -8.5492E-03  4.3217E-03
 SE:             2.9483E-02  2.0945E-02  2.3040E-02  2.8067E-02  1.4230E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7722E-01  4.6804E-01  7.6076E-01  7.6067E-01  7.6135E-01

 ETASHRINKSD(%)  1.2287E+00  2.9832E+01  2.2813E+01  5.9734E+00  5.2328E+01
 ETASHRINKVR(%)  2.4424E+00  5.0765E+01  4.0422E+01  1.1590E+01  7.7273E+01
 EBVSHRINKSD(%)  1.5033E+00  2.9763E+01  2.5073E+01  5.6369E+00  5.3454E+01
 EBVSHRINKVR(%)  2.9841E+00  5.0668E+01  4.3860E+01  1.0956E+01  7.8335E+01
 RELATIVEINF(%)  9.6966E+01  1.1916E+01  8.8111E+00  6.3989E+01  1.8359E+00
 EPSSHRINKSD(%)  3.1559E+01
 EPSSHRINKVR(%)  5.3159E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1946.5475955440634     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -843.82135569845627     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.61
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1946.548       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  3.43E-01  2.50E-01  1.11E+00  2.88E-01  9.83E-01  1.25E+00  1.22E+00  1.23E+00  4.14E-01  2.19E+00
 


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
+        1.08E+03
 
 TH 2
+       -2.25E+06  8.21E+05
 
 TH 3
+       -3.41E+02  3.51E+03  8.74E+05
 
 TH 4
+       -3.37E+06  1.22E+06  1.36E+03  1.83E+06
 
 TH 5
+       -2.82E+02 -3.11E+03  8.41E+05  1.17E+03  8.51E+05
 
 TH 6
+        1.69E+00  1.82E+02 -1.84E+02 -3.98E+00 -1.77E+02  1.95E+02
 
 TH 7
+       -2.89E-02  6.85E+05  5.26E+02  1.03E+06  3.93E+02 -3.89E-02  3.09E+01
 
 TH 8
+       -2.08E+06  7.58E+05  1.08E+03  1.14E+06  1.15E+03  1.74E+02  6.35E+05  7.03E+05
 
 TH 9
+        5.23E+00  1.00E+03 -9.94E+02 -1.09E+06 -9.31E+02  1.21E+00  3.88E+00 -6.74E+05  9.72E+01
 
 TH10
+        4.66E-01  5.80E+00  8.02E+00 -4.80E+00  2.27E+02  1.41E+00  2.36E+01  7.80E+05  9.30E+00  7.34E+01
 
 TH11
+       -6.95E+01  3.20E+02  1.47E+05  2.56E+02 -3.68E+02 -3.11E+01  8.89E+01  2.04E+02 -1.69E+02  1.04E-01  2.49E+04
 
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
 #CPUT: Total CPU Time in Seconds,       31.990
Stop Time:
Thu Sep 30 05:38:25 CDT 2021
