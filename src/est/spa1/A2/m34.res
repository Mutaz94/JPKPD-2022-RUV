Wed Sep 29 23:16:53 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat34.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1668.22306302168        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8800E+02 -7.6553E+01 -1.6859E+01 -3.5553E+01  1.1735E+02 -1.1442E+01 -9.0903E+00  7.1087E-01  3.3133E+01 -4.9151E+01
            -7.4594E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1848.52452264089        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      149
 NPARAMETR:  9.8679E-01  1.0892E+00  1.0039E+00  1.0508E+00  9.6405E-01  1.0982E+00  9.9419E-01  9.6656E-01  8.5679E-01  1.0642E+00
             1.7206E+00
 PARAMETER:  8.6702E-02  1.8543E-01  1.0390E-01  1.4957E-01  6.3392E-02  1.9366E-01  9.4177E-02  6.5991E-02 -5.4564E-02  1.6224E-01
             6.4266E-01
 GRADIENT:   5.7460E+01  5.6461E+00 -8.4354E+00  2.5048E+01  7.6198E+00  7.7508E+00 -3.1307E+00  6.0857E+00  1.0341E+01 -4.0914E+00
             1.3683E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1850.82998389871        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  9.8330E-01  9.6469E-01  1.1448E+00  1.1264E+00  9.9084E-01  1.0782E+00  1.2498E+00  6.0361E-01  7.7548E-01  1.2786E+00
             1.7141E+00
 PARAMETER:  8.3161E-02  6.4050E-02  2.3524E-01  2.1902E-01  9.0802E-02  1.7525E-01  3.2301E-01 -4.0482E-01 -1.5428E-01  3.4577E-01
             6.3887E-01
 GRADIENT:   5.7238E+01  6.8238E+00 -3.8369E+00  9.4049E+00  3.0428E-01  1.8657E+00  5.6003E+00  1.2901E+00  1.0201E+01  1.1285E+01
             1.5313E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1853.70954135271        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  9.5262E-01  8.4998E-01  9.2252E-01  1.1745E+00  8.2790E-01  1.0679E+00  1.3420E+00  2.8432E-01  6.9788E-01  1.0469E+00
             1.6830E+00
 PARAMETER:  5.1465E-02 -6.2545E-02  1.9349E-02  2.6084E-01 -8.8864E-02  1.6574E-01  3.9416E-01 -1.1576E+00 -2.5971E-01  1.4580E-01
             6.2058E-01
 GRADIENT:  -3.3151E+00  4.9131E+00 -5.8833E-01  7.7117E+00 -6.2819E-01 -3.0632E-01 -9.3888E-01  3.2915E-01 -8.9594E-01 -5.3385E-01
            -1.0342E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1854.40732023331        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      676
 NPARAMETR:  9.5230E-01  5.3920E-01  1.1191E+00  1.3749E+00  8.1652E-01  1.0622E+00  1.8261E+00  2.8937E-02  6.6088E-01  1.1311E+00
             1.6955E+00
 PARAMETER:  5.1123E-02 -5.1767E-01  2.1249E-01  4.1841E-01 -1.0270E-01  1.6035E-01  7.0216E-01 -3.4427E+00 -3.1418E-01  2.2321E-01
             6.2795E-01
 GRADIENT:   5.1054E+00  6.1197E+00  1.0882E+00  1.6661E+01 -4.2772E+00 -4.2821E-01  4.1251E-01  1.6562E-03 -7.6311E-02  6.7488E-01
             1.6528E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1854.76356866549        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      854
 NPARAMETR:  9.4328E-01  3.1290E-01  1.2707E+00  1.5188E+00  8.1688E-01  1.0615E+00  2.4813E+00  1.0000E-02  6.3982E-01  1.1816E+00
             1.6915E+00
 PARAMETER:  4.1611E-02 -1.0619E+00  3.3959E-01  5.1790E-01 -1.0227E-01  1.5965E-01  1.0088E+00 -7.6182E+00 -3.4656E-01  2.6685E-01
             6.2564E-01
 GRADIENT:  -5.7649E+00  3.6354E+00  2.0862E+00  1.9945E+01 -4.9435E+00  4.7524E-01 -4.6317E-01  0.0000E+00 -1.3978E+00 -9.7045E-02
            -2.8272E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1855.11082529675        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1034
 NPARAMETR:  9.4416E-01  1.4729E-01  1.3905E+00  1.6263E+00  8.1974E-01  1.0562E+00  3.6513E+00  1.0000E-02  6.3506E-01  1.2284E+00
             1.7004E+00
 PARAMETER:  4.2540E-02 -1.8154E+00  4.2965E-01  5.8634E-01 -9.8763E-02  1.5469E-01  1.3951E+00 -1.4611E+01 -3.5403E-01  3.0574E-01
             6.3088E-01
 GRADIENT:   1.2449E+00  2.7678E+00  1.3099E+00  2.3270E+01 -5.1429E+00 -4.7179E-01  1.2741E+00  0.0000E+00  6.1006E-01  3.4727E-01
             4.1424E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1856.22920004413        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1212
 NPARAMETR:  9.3943E-01  4.0861E-02  1.4461E+00  1.6935E+00  8.1366E-01  1.0609E+00  6.1593E+00  1.0000E-02  6.2681E-01  1.2005E+00
             1.6980E+00
 PARAMETER:  3.7523E-02 -3.0976E+00  4.6888E-01  6.2678E-01 -1.0622E-01  1.5909E-01  1.9180E+00 -2.7680E+01 -3.6712E-01  2.8275E-01
             6.2943E-01
 GRADIENT:  -5.3887E+00 -1.5444E+00  2.8109E+00  4.8046E+01 -6.6982E+00  1.5230E+00 -6.0609E+00  0.0000E+00  7.1605E+00 -4.5236E+00
            -3.5583E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1857.18099815446        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1382
 NPARAMETR:  9.4178E-01  3.0574E-02  1.4618E+00  1.6877E+00  8.2455E-01  1.0556E+00  7.0143E+00  1.0000E-02  6.1040E-01  1.2481E+00
             1.6998E+00
 PARAMETER:  3.9561E-02 -3.3771E+00  4.7992E-01  6.2165E-01 -9.3452E-02  1.5374E-01  2.0533E+00 -3.0717E+01 -3.9208E-01  3.2489E-01
             6.3058E-01
 GRADIENT:  -6.3518E-01  1.7117E+01  1.2556E-01 -9.4708E+01 -3.8530E-01 -2.1131E-01  2.3918E+01  0.0000E+00  4.9542E-01  9.8936E-01
             1.1486E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1382
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8541E-03  1.5136E-02 -7.4370E-05 -1.4681E-02 -2.6217E-02
 SE:             2.9599E-02  8.9638E-03  1.2120E-04  2.7508E-02  2.3472E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5005E-01  9.1314E-02  5.3947E-01  5.9356E-01  2.6402E-01

 ETASHRINKSD(%)  8.4048E-01  6.9970E+01  9.9594E+01  7.8441E+00  2.1365E+01
 ETASHRINKVR(%)  1.6739E+00  9.0982E+01  9.9998E+01  1.5073E+01  3.8165E+01
 EBVSHRINKSD(%)  8.9981E-01  8.0702E+01  9.9581E+01  6.8154E+00  1.8381E+01
 EBVSHRINKVR(%)  1.7915E+00  9.6276E+01  9.9998E+01  1.3166E+01  3.3383E+01
 RELATIVEINF(%)  9.7993E+01  1.9588E+00  2.1446E-04  3.4409E+01  8.6579E+00
 EPSSHRINKSD(%)  2.8833E+01
 EPSSHRINKVR(%)  4.9353E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1857.1809981544643     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -938.24246494979161     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.84
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1857.181       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.41E-01  3.09E-02  1.46E+00  1.68E+00  8.24E-01  1.06E+00  7.05E+00  1.00E-02  6.11E-01  1.25E+00  1.70E+00
 


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
+        1.11E+03
 
 TH 2
+       -6.12E+01  4.44E+05
 
 TH 3
+       -6.83E+00 -2.76E+00  9.27E+01
 
 TH 4
+       -1.11E+01 -2.89E+04 -6.77E+01  2.38E+03
 
 TH 5
+        1.62E+01 -3.43E+02 -2.32E+02 -8.76E+00  7.52E+02
 
 TH 6
+        5.68E+00 -3.79E+01 -1.07E-01 -3.72E+00  1.75E+00  1.69E+02
 
 TH 7
+       -1.98E-01  2.14E+03 -3.88E-01 -2.15E+02 -2.88E-01 -2.60E-01  2.39E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -8.05E+00 -2.48E+02  2.65E+01 -4.17E+01 -1.36E+02 -1.27E+01 -5.03E+02  0.00E+00  2.94E+04
 
 TH10
+       -9.69E+00  2.60E+02  1.21E+01 -3.01E+01 -1.65E+02 -9.10E+00  1.54E+00  0.00E+00  1.70E+04  1.01E+04
 
 TH11
+       -1.38E+01 -2.90E-01 -7.26E+00 -2.32E+01 -3.57E+00  2.87E+00 -3.94E-02  0.00E+00  6.51E+03  3.83E+03  1.56E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.365
Stop Time:
Wed Sep 29 23:17:25 CDT 2021
