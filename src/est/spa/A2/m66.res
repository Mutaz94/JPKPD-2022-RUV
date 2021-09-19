Sat Sep 18 09:59:38 CDT 2021
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
$DATA ../../../../data/spa/A2/dat66.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1111.23908079258        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7173E+02 -5.0037E+01  4.6804E+01 -1.2328E+02  4.6909E+01  1.2794E+01 -1.7303E+01 -7.2706E+00 -2.4191E+01 -4.9257E+01
            -9.2183E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1410.16770652701        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.3423E-01  1.0432E+00  1.0099E+00  1.0761E+00  1.0069E+00  8.6627E-01  9.8001E-01  9.2200E-01  8.9704E-01  9.1995E-01
             2.2442E+00
 PARAMETER:  3.1965E-02  1.4234E-01  1.0989E-01  1.7337E-01  1.0690E-01 -4.3559E-02  7.9806E-02  1.8793E-02 -8.6514E-03  1.6565E-02
             9.0835E-01
 GRADIENT:  -4.6980E+01  5.2878E+00 -3.3040E+00  1.0895E+01  1.3038E+01 -2.9498E+01 -4.1164E+00  2.7985E+00 -5.4688E+00 -2.2886E-01
            -2.4023E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1416.13823480228        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.3846E-01  8.0932E-01  6.5079E-01  1.2187E+00  6.8065E-01  9.1591E-01  1.4289E+00  4.4051E-01  8.5791E-01  4.1536E-01
             2.2949E+00
 PARAMETER:  3.6486E-02 -1.1156E-01 -3.2956E-01  2.9779E-01 -2.8471E-01  1.2168E-02  4.5693E-01 -7.1982E-01 -5.3251E-02 -7.7861E-01
             9.3071E-01
 GRADIENT:  -3.9949E+01  2.6554E+01 -1.0348E+01  7.6016E+01  1.9186E+01 -9.6554E+00  2.2464E+00  5.5224E-01  6.1385E+00 -1.1599E+00
            -7.4115E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1421.33491279363        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.6072E-01  6.6216E-01  3.8272E-01  1.1902E+00  4.4468E-01  9.5139E-01  1.5000E+00  1.5300E-01  7.4330E-01  3.0781E-01
             2.2180E+00
 PARAMETER:  5.9926E-02 -3.1225E-01 -8.6044E-01  2.7410E-01 -7.1041E-01  5.0170E-02  5.0547E-01 -1.7773E+00 -1.9665E-01 -1.0783E+00
             8.9659E-01
 GRADIENT:   1.4803E+01  2.9967E+01  7.5016E+00  4.2733E+01 -2.0291E+01  2.1083E+00 -6.5444E-01  5.4300E-02 -1.1487E+01  8.5929E-02
             7.1473E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1423.27006378723        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      343
 NPARAMETR:  9.5095E-01  5.2653E-01  3.4452E-01  1.2003E+00  3.8662E-01  9.4718E-01  1.6375E+00  5.7440E-02  7.8373E-01  3.3895E-01
             2.1408E+00
 PARAMETER:  4.9709E-02 -5.4144E-01 -9.6560E-01  2.8260E-01 -8.5030E-01  4.5730E-02  5.9320E-01 -2.7570E+00 -1.4369E-01 -9.8191E-01
             8.6120E-01
 GRADIENT:  -1.1609E+01  5.0138E-01  7.1914E+00 -1.3741E+01 -9.6819E+00 -7.1010E-01 -5.2113E+00  5.2418E-03 -6.4050E-01 -7.2190E-01
            -3.9796E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1424.22278296555        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.5415E-01  4.1088E-01  3.9519E-01  1.2858E+00  3.9365E-01  9.4532E-01  2.0933E+00  2.4377E-02  7.6513E-01  4.4505E-01
             2.1578E+00
 PARAMETER:  5.3067E-02 -7.8946E-01 -8.2840E-01  3.5137E-01 -8.3229E-01  4.3771E-02  8.3875E-01 -3.6141E+00 -1.6771E-01 -7.0956E-01
             8.6909E-01
 GRADIENT:  -1.7763E+00  1.4549E-01 -1.4862E+00 -2.9083E+00  1.2316E+00  1.2214E-01  2.4558E-01  4.8359E-03  1.7217E-01  4.5678E-01
             7.8030E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1424.25557880619        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      694            RESET HESSIAN, TYPE II
 NPARAMETR:  9.5430E-01  3.8795E-01  4.0955E-01  1.3048E+00  3.9834E-01  9.4389E-01  2.1964E+00  1.1116E-02  7.6111E-01  4.6117E-01
             2.1603E+00
 PARAMETER:  5.3226E-02 -8.4689E-01 -7.9269E-01  3.6604E-01 -8.2045E-01  4.2253E-02  8.8684E-01 -4.3994E+00 -1.7298E-01 -6.7399E-01
             8.7026E-01
 GRADIENT:   6.6861E+00  1.0313E+00  1.7698E+00  6.9387E+00  4.1651E+00  6.1603E-01  6.7045E-01  1.0002E-03  2.2596E-01  5.8261E-02
             6.9409E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1424.25566581566        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      851
 NPARAMETR:  9.5426E-01  3.8800E-01  4.0942E-01  1.3047E+00  3.9827E-01  9.4391E-01  2.1958E+00  1.0000E-02  7.6117E-01  4.6126E-01
             2.1601E+00
 PARAMETER:  5.3181E-02 -8.4674E-01 -7.9301E-01  3.6598E-01 -8.2063E-01  4.2279E-02  8.8652E-01 -5.4654E+00 -1.7290E-01 -6.7379E-01
             8.7016E-01
 GRADIENT:  -7.7057E-02 -6.2112E-03 -2.3187E-02  2.7473E-02  1.7654E-02  3.9622E-03 -2.0933E-03  0.0000E+00  4.3896E-03  2.8534E-03
            -2.1481E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1424.25566850749        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  9.5429E-01  3.8802E-01  4.0941E-01  1.3047E+00  3.9826E-01  9.4390E-01  2.1957E+00  1.0000E-02  7.6115E-01  4.6120E-01
             2.1602E+00
 PARAMETER:  5.3215E-02 -8.4670E-01 -7.9303E-01  3.6597E-01 -8.2064E-01  4.2265E-02  8.8652E-01 -5.4321E+00 -1.7293E-01 -6.7392E-01
             8.7022E-01
 GRADIENT:   4.0560E-04 -3.2585E-03  1.0009E-03 -3.9232E-04 -1.2452E-02 -5.2087E-04  4.8630E-04  0.0000E+00  1.2696E-04  7.9138E-04
            -7.1570E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      943
 NO. OF SIG. DIGITS IN FINAL EST.:  4.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0076E-03  4.2219E-02 -5.0704E-04 -2.7893E-02  1.8931E-02
 SE:             2.9335E-02  2.0175E-02  2.3505E-04  2.4727E-02  1.5070E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4544E-01  3.6378E-02  3.0990E-02  2.5931E-01  2.0903E-01

 ETASHRINKSD(%)  1.7247E+00  3.2412E+01  9.9213E+01  1.7162E+01  4.9515E+01
 ETASHRINKVR(%)  3.4197E+00  5.4319E+01  9.9994E+01  3.1379E+01  7.4512E+01
 EBVSHRINKSD(%)  1.9515E+00  3.5071E+01  9.9158E+01  1.6056E+01  4.7033E+01
 EBVSHRINKVR(%)  3.8649E+00  5.7843E+01  9.9993E+01  2.9534E+01  7.1945E+01
 RELATIVEINF(%)  9.5620E+01  9.3672E+00  2.1976E-04  2.5970E+01  7.6085E-01
 EPSSHRINKSD(%)  3.4845E+01
 EPSSHRINKVR(%)  5.7549E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1424.2556685074935     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -689.10484194375533     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1424.256       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.54E-01  3.88E-01  4.09E-01  1.30E+00  3.98E-01  9.44E-01  2.20E+00  1.00E-02  7.61E-01  4.61E-01  2.16E+00
 


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
+       -4.18E+01  6.82E+02
 
 TH 3
+        5.04E+00  9.44E+02  5.25E+03
 
 TH 4
+       -3.43E+01  2.34E+02 -6.54E+02  8.08E+02
 
 TH 5
+        6.79E+01 -1.70E+03 -6.87E+03  3.67E+02  9.70E+03
 
 TH 6
+       -3.78E+00 -7.18E+00  1.76E+01 -1.10E+01  6.67E-01  2.04E+02
 
 TH 7
+        2.59E+00  5.10E+01 -2.77E+01 -6.72E+00  9.27E+00  3.31E-01  1.40E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.18E+01 -1.23E+01 -9.54E+00 -1.55E+01  5.26E+01  1.62E+00  4.52E+00  0.00E+00  1.82E+02
 
 TH10
+       -3.51E+00  6.04E+00 -2.35E+02 -2.83E+01  2.52E+02  1.71E+00  6.88E+00  0.00E+00  6.69E-01  9.96E+01
 
 TH11
+       -1.42E+01 -4.27E+00 -6.79E+01 -1.40E+01  5.32E+01  3.14E+00  2.36E+00  0.00E+00  1.43E+01  2.38E+01  5.87E+01
 
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
 #CPUT: Total CPU Time in Seconds,       17.222
Stop Time:
Sat Sep 18 09:59:56 CDT 2021
