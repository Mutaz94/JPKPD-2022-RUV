Wed Sep 29 23:40:57 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1122.09186923083        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6574E+02 -1.3903E+01  1.9589E+01 -2.3492E+01  1.2417E+02  2.9984E+01 -3.8713E+01 -2.0741E+00 -7.7447E+00 -7.1142E+01
            -1.8198E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1717.24452967462        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0071E+00  1.0695E+00  1.1511E+00  1.0329E+00  1.0511E+00  9.0604E-01  1.0750E+00  8.6384E-01  9.7377E-01  8.7217E-01
             2.4726E+00
 PARAMETER:  1.0709E-01  1.6716E-01  2.4075E-01  1.3241E-01  1.4985E-01  1.3262E-03  1.7235E-01 -4.6370E-02  7.3415E-02 -3.6768E-02
             1.0053E+00
 GRADIENT:   1.2808E+02 -1.8442E+01 -1.5827E+01 -1.0491E+01  2.2006E+01 -2.4947E+01 -2.0579E+00  5.6712E+00  2.1996E+00  3.6748E+00
             4.0920E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1723.31015847413        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.9076E-01  7.8813E-01  1.0946E+00  1.2379E+00  8.7813E-01  9.5994E-01  1.1796E+00  3.0024E-01  9.4848E-01  8.4486E-01
             2.3816E+00
 PARAMETER:  9.0713E-02 -1.3809E-01  1.9041E-01  3.1340E-01 -2.9959E-02  5.9117E-02  2.6515E-01 -1.1032E+00  4.7107E-02 -6.8582E-02
             9.6775E-01
 GRADIENT:   8.9515E+01  1.6926E+01 -4.0020E+00  7.7828E+01 -1.1669E-01 -1.1148E+00 -2.4843E+00  1.0043E+00  6.4617E+00  5.7401E-01
             2.0141E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1725.23178910881        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      314
 NPARAMETR:  9.7832E-01  7.1301E-01  1.0118E+00  1.2651E+00  8.1601E-01  9.8371E-01  1.5191E+00  1.9771E-01  8.4390E-01  8.2035E-01
             2.3182E+00
 PARAMETER:  7.8084E-02 -2.3826E-01  1.1170E-01  3.3519E-01 -1.0333E-01  8.3577E-02  5.1812E-01 -1.5210E+00 -6.9723E-02 -9.8019E-02
             9.4078E-01
 GRADIENT:  -5.2912E+00  9.5746E+00 -1.3149E+00  1.5063E+01  8.7200E-02  2.6263E+00 -2.9505E+00  4.8290E-01 -2.5889E+00  2.7002E-01
            -5.9805E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1728.05125939549        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  9.7787E-01  3.9054E-01  7.1360E-01  1.3931E+00  5.5723E-01  9.6589E-01  2.3127E+00  1.0000E-02  8.1660E-01  7.3004E-01
             2.2673E+00
 PARAMETER:  7.7618E-02 -8.4022E-01 -2.3743E-01  4.3156E-01 -4.8477E-01  6.5290E-02  9.3843E-01 -4.5077E+00 -1.0260E-01 -2.1465E-01
             9.1859E-01
 GRADIENT:   1.4760E+00  9.1629E+00  1.2651E+01  1.7430E+01 -2.4301E+01 -1.9285E+00  3.3881E-01  1.2437E-03  3.0513E-01 -2.5279E-01
            -5.8963E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1729.29790613635        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      669
 NPARAMETR:  9.7111E-01  2.2627E-01  8.7785E-01  1.5169E+00  6.0559E-01  9.6763E-01  2.8912E+00  1.0000E-02  8.2649E-01  8.3896E-01
             2.2675E+00
 PARAMETER:  7.0684E-02 -1.3860E+00 -3.0275E-02  5.1670E-01 -4.0155E-01  6.7091E-02  1.1617E+00 -7.0886E+00 -9.0563E-02 -7.5596E-02
             9.1869E-01
 GRADIENT:  -2.9332E+00  5.8185E+00  9.6583E+00  1.9972E+01 -1.8980E+01 -2.2394E-01  2.3934E+00  0.0000E+00  2.9355E+00  7.8002E-01
            -3.5610E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1730.31291855049        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      847
 NPARAMETR:  9.7008E-01  1.0157E-01  9.4548E-01  1.5904E+00  6.1961E-01  9.6515E-01  4.1448E+00  1.0000E-02  7.9010E-01  8.5530E-01
             2.2661E+00
 PARAMETER:  6.9619E-02 -2.1870E+00  4.3937E-02  5.6400E-01 -3.7866E-01  6.4524E-02  1.5219E+00 -1.1537E+01 -1.3559E-01 -5.6308E-02
             9.1805E-01
 GRADIENT:   2.0791E+00  1.5449E-01  2.8351E+00  1.9496E+01 -5.1758E+00 -4.5825E-01 -3.0342E+00  0.0000E+00 -1.0797E+00 -3.2585E-01
            -1.8490E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1730.72330720858        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  9.6525E-01  5.1238E-02  8.9805E-01  1.6030E+00  5.9022E-01  9.6671E-01  5.7698E+00  1.0000E-02  7.9846E-01  8.1516E-01
             2.2730E+00
 PARAMETER:  6.4636E-02 -2.8713E+00 -7.5307E-03  5.7190E-01 -4.2725E-01  6.6143E-02  1.8526E+00 -1.5719E+01 -1.2507E-01 -1.0437E-01
             9.2108E-01
 GRADIENT:  -5.7388E+00 -1.6966E+00  3.2873E+00  1.9135E+01 -4.3557E+00  5.8740E-01 -4.9760E+00  0.0000E+00  3.7522E+00 -6.0723E-01
             3.7383E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1731.96694271853        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1209             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6574E-01  1.9967E-02  7.0143E-01  1.5526E+00  4.9286E-01  9.6338E-01  8.8036E+00  1.0000E-02  8.0814E-01  7.5096E-01
             2.2275E+00
 PARAMETER:  6.5140E-02 -3.8137E+00 -2.5464E-01  5.3995E-01 -6.0753E-01  6.2692E-02  2.2752E+00 -2.1935E+01 -1.1301E-01 -1.8640E-01
             9.0086E-01
 GRADIENT:   7.5084E+01  4.8043E+00  1.5751E+00  1.8462E+02  2.2462E+01  5.6217E+00  4.1759E+01  0.0000E+00  5.1980E-01 -5.7596E-02
             3.8070E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1731.97459012431        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1371
 NPARAMETR:  9.6540E-01  1.9825E-02  7.0173E-01  1.5535E+00  4.9229E-01  9.6360E-01  8.7662E+00  1.0000E-02  8.0887E-01  7.4871E-01
             2.2346E+00
 PARAMETER:  6.4723E-02 -3.8143E+00 -2.5405E-01  5.3957E-01 -6.0765E-01  6.2950E-02  2.2748E+00 -2.1935E+01 -1.1150E-01 -1.8753E-01
             9.0316E-01
 GRADIENT:  -4.5440E-01  7.3488E+01  5.2936E-01 -2.7439E+02  2.4462E+02  3.0174E-02  1.5943E+02  0.0000E+00  4.2184E-01  3.9914E-01
            -2.4965E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1371
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3325E-03  1.6541E-02 -2.0459E-06 -1.1205E-02 -1.0606E-02
 SE:             2.9315E-02  8.5539E-03  1.7929E-04  2.7817E-02  2.1101E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6375E-01  5.3144E-02  9.9090E-01  6.8707E-01  6.1523E-01

 ETASHRINKSD(%)  1.7918E+00  7.1343E+01  9.9399E+01  6.8106E+00  2.9308E+01
 ETASHRINKVR(%)  3.5515E+00  9.1788E+01  9.9996E+01  1.3157E+01  5.0026E+01
 EBVSHRINKSD(%)  1.7274E+00  8.0487E+01  9.9297E+01  6.4097E+00  2.8398E+01
 EBVSHRINKVR(%)  3.4249E+00  9.6192E+01  9.9995E+01  1.2409E+01  4.8732E+01
 RELATIVEINF(%)  9.6337E+01  2.9361E+00  3.0697E-04  4.3118E+01  3.2447E+00
 EPSSHRINKSD(%)  2.6553E+01
 EPSSHRINKVR(%)  4.6055E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1731.9745901243100     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -813.03605691963730     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.85
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1731.975       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.65E-01  2.00E-02  7.02E-01  1.55E+00  4.93E-01  9.64E-01  8.80E+00  1.00E-02  8.09E-01  7.50E-01  2.23E+00
 


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
+        1.24E+03
 
 TH 2
+       -4.15E+02  1.96E+06
 
 TH 3
+       -3.21E+01  6.50E+03  1.26E+03
 
 TH 4
+       -3.98E+00 -2.92E+04 -4.10E+02  5.84E+03
 
 TH 5
+        6.88E+00  8.22E+04 -2.97E+04 -1.47E+04  4.49E+04
 
 TH 6
+        4.42E+00  1.97E+02  1.39E+01 -1.92E+01  2.57E+01  1.98E+02
 
 TH 7
+       -1.64E+00  4.99E+03  2.91E+01 -1.12E+02  3.14E+02  1.12E+00  2.87E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.76E+00  4.44E+03  1.68E+02 -2.15E+02 -5.53E+04  9.37E+00  2.07E+01  0.00E+00  3.88E+02
 
 TH10
+       -9.31E+00  1.53E+03  6.00E+01 -7.84E+01  2.13E+02  7.77E+00  6.96E+00  0.00E+00  7.37E+01  1.36E+02
 
 TH11
+       -1.83E+01  2.45E+02  1.05E+04 -2.08E+01 -6.17E+03  4.56E+00 -5.43E+01  0.00E+00  3.84E+01  4.46E+01  1.03E+03
 
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
 #CPUT: Total CPU Time in Seconds,       32.831
Stop Time:
Wed Sep 29 23:41:35 CDT 2021
