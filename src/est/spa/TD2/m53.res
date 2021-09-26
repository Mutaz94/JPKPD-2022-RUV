Sat Sep 25 13:36:32 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat53.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m53.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1678.37766413118        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.9676E+01 -7.4237E+01 -3.3683E+01 -4.8757E+01  3.1777E+01 -2.3613E+01  1.2251E+01  8.1433E+00  5.3128E+01  1.5568E+01
            -1.0991E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1692.79292671343        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7518E-01  1.1059E+00  1.0639E+00  9.9121E-01  1.0347E+00  1.0584E+00  9.1375E-01  9.7031E-01  7.3337E-01  9.2208E-01
             1.0220E+00
 PARAMETER:  7.4870E-02  2.0067E-01  1.6190E-01  9.1167E-02  1.3412E-01  1.5675E-01  9.7974E-03  6.9861E-02 -2.1011E-01  1.8872E-02
             1.2175E-01
 GRADIENT:   3.1974E+01  1.8150E+01  5.8452E+00  1.5956E+01  1.2815E+00  5.8103E+00 -1.4792E+00 -2.3427E+00 -1.0858E+00 -3.2154E+00
            -1.2106E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1692.98381083945        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.7307E-01  1.0647E+00  9.6681E-01  1.0117E+00  9.7372E-01  1.0634E+00  9.7934E-01  8.9032E-01  6.9835E-01  8.7354E-01
             1.0210E+00
 PARAMETER:  7.2701E-02  1.6266E-01  6.6247E-02  1.1168E-01  7.3371E-02  1.6147E-01  7.9125E-02 -1.6177E-02 -2.5903E-01 -3.5201E-02
             1.2075E-01
 GRADIENT:   2.6427E+01  1.2538E+01 -2.2511E+00  2.0814E+01  6.2677E+00  7.6800E+00  7.3838E-02  2.1388E-01 -1.3373E+00 -1.4726E+00
             5.2825E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1693.21418481365        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      340
 NPARAMETR:  9.8106E-01  1.1374E+00  7.5662E-01  9.4981E-01  8.9300E-01  1.0756E+00  9.4536E-01  6.1909E-01  7.2088E-01  8.1339E-01
             1.0152E+00
 PARAMETER:  8.0876E-02  2.2874E-01 -1.7890E-01  4.8510E-02 -1.3165E-02  1.7288E-01  4.3810E-02 -3.7950E-01 -2.2728E-01 -1.0654E-01
             1.1510E-01
 GRADIENT:  -2.5482E-01  4.3532E+00  8.0535E-01  3.8625E+00 -4.8738E+00  5.9834E-01  6.5887E-01  5.9788E-01  2.9335E-01  6.1728E-01
            -3.4239E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1693.40790194006        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  9.8189E-01  1.2735E+00  6.0897E-01  8.5259E-01  8.7931E-01  1.0756E+00  8.6348E-01  3.2163E-01  7.6585E-01  7.9147E-01
             1.0169E+00
 PARAMETER:  8.1722E-02  3.4179E-01 -3.9598E-01 -5.9477E-02 -2.8615E-02  1.7289E-01 -4.6785E-02 -1.0343E+00 -1.6677E-01 -1.3386E-01
             1.1678E-01
 GRADIENT:  -8.8470E-01 -1.7042E+00 -2.1304E+00  1.7839E+00  2.2310E+00  5.1860E-02  3.6364E-01  1.5552E-01  5.5517E-01  7.0596E-01
             5.1389E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1693.45498797861        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      690
 NPARAMETR:  9.8254E-01  1.3810E+00  5.5521E-01  7.8323E-01  9.0493E-01  1.0755E+00  8.0880E-01  1.8855E-01  8.0911E-01  7.9872E-01
             1.0170E+00
 PARAMETER:  8.2386E-02  4.2282E-01 -4.8841E-01 -1.4433E-01  9.8263E-05  1.7283E-01 -1.1220E-01 -1.5684E+00 -1.1182E-01 -1.2474E-01
             1.1686E-01
 GRADIENT:   7.0641E-02  2.9008E-01 -2.7143E-01  4.0627E-01  3.4758E-01 -5.1922E-02  6.7402E-02  3.5271E-02 -3.6825E-02  6.8343E-02
             1.0854E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1693.46621683467        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      870
 NPARAMETR:  9.8264E-01  1.3831E+00  5.4942E-01  7.8038E-01  9.0327E-01  1.0762E+00  8.0749E-01  8.5006E-02  8.1154E-01  7.9862E-01
             1.0168E+00
 PARAMETER:  8.2483E-02  4.2435E-01 -4.9888E-01 -1.4798E-01 -1.7286E-03  1.7343E-01 -1.1383E-01 -2.3650E+00 -1.0882E-01 -1.2487E-01
             1.1671E-01
 GRADIENT:   2.4606E-01 -1.9105E+00 -4.0394E-01 -9.6248E-01  1.1114E+00  1.7201E-01  8.4626E-02  6.3144E-03  1.0284E-01  1.2863E-01
             8.1057E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1693.47200068698        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1045
 NPARAMETR:  9.8258E-01  1.3608E+00  5.5292E-01  7.9453E-01  8.9152E-01  1.0760E+00  8.1834E-01  2.4750E-02  7.9965E-01  7.9062E-01
             1.0167E+00
 PARAMETER:  8.2425E-02  4.0811E-01 -4.9255E-01 -1.3001E-01 -1.4833E-02  1.7322E-01 -1.0047E-01 -3.5989E+00 -1.2358E-01 -1.3494E-01
             1.1652E-01
 GRADIENT:   8.9329E-02 -1.3945E-01 -5.1576E-02 -5.1858E-02  1.5397E-01  7.1226E-02  2.8880E-02  5.0577E-04 -2.6215E-02 -1.9614E-02
             1.3326E-04

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1693.47223219783        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  9.8254E-01  1.3614E+00  5.5265E-01  7.9420E-01  8.9158E-01  1.0758E+00  8.1791E-01  1.0000E-02  8.0012E-01  7.9082E-01
             1.0167E+00
 PARAMETER:  8.2378E-02  4.0853E-01 -4.9305E-01 -1.3042E-01 -1.4745E-02  1.7304E-01 -1.0110E-01 -4.9110E+00 -1.2305E-01 -1.3479E-01
             1.1651E-01
 GRADIENT:  -1.6960E-02  4.6842E-03 -4.2207E-03  1.4680E-02  2.0978E-02 -7.0839E-03 -6.2905E-03  0.0000E+00 -5.8065E-03 -5.5050E-03
            -2.4813E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1178
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.9443E-06 -1.1748E-02 -3.3719E-04  6.6170E-03 -1.8928E-02
 SE:             2.9866E-02  2.3629E-02  1.5530E-04  2.2407E-02  2.2327E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9981E-01  6.1905E-01  2.9912E-02  7.6776E-01  3.9656E-01

 ETASHRINKSD(%)  1.0000E-10  2.0840E+01  9.9480E+01  2.4934E+01  2.5202E+01
 ETASHRINKVR(%)  1.0000E-10  3.7337E+01  9.9997E+01  4.3651E+01  4.4053E+01
 EBVSHRINKSD(%)  3.6394E-01  2.0828E+01  9.9528E+01  2.5941E+01  2.4457E+01
 EBVSHRINKVR(%)  7.2656E-01  3.7318E+01  9.9998E+01  4.5152E+01  4.2933E+01
 RELATIVEINF(%)  9.9218E+01  2.2920E+00  1.2397E-04  1.8542E+00  4.6170E+00
 EPSSHRINKSD(%)  4.3271E+01
 EPSSHRINKVR(%)  6.7819E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1693.4722321978311     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -958.32140563409291     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.46
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1693.472       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.36E+00  5.53E-01  7.94E-01  8.92E-01  1.08E+00  8.18E-01  1.00E-02  8.00E-01  7.91E-01  1.02E+00
 


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
+        9.88E+02
 
 TH 2
+       -7.37E+00  5.75E+02
 
 TH 3
+        1.07E+01  2.58E+02  8.62E+02
 
 TH 4
+       -1.77E+01  5.07E+02 -6.15E+02  1.52E+03
 
 TH 5
+       -2.90E+00 -4.02E+02 -9.27E+02  5.73E+02  1.27E+03
 
 TH 6
+       -1.03E-01 -1.49E+00  2.35E+00 -5.08E+00 -2.26E+00  1.70E+02
 
 TH 7
+       -2.31E-01  2.31E+01 -2.77E+01 -1.59E+01  4.28E-01 -7.59E-01  1.24E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.99E-01 -2.12E+01 -5.01E+01  5.40E+01  3.85E+00  5.41E-01  3.20E+01  0.00E+00  1.10E+02
 
 TH10
+        7.61E-01 -1.63E+01 -5.73E+01 -1.31E+01 -7.41E+01  2.10E-01  2.12E+01  0.00E+00  1.47E+01  1.14E+02
 
 TH11
+       -3.16E+00 -1.75E+01 -3.16E+01 -2.42E+00 -5.64E-03  1.69E+00  9.80E+00  0.00E+00  1.62E+01  2.29E+01  2.08E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.883
Stop Time:
Sat Sep 25 13:36:54 CDT 2021
