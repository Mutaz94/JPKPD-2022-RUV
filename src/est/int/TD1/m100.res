Sat Sep 25 04:30:52 CDT 2021
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
$DATA ../../../../data/int/TD1/dat100.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3741.99879392522        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3506E+02 -6.3401E+01  1.9284E+01  5.5635E+01  3.2155E+01 -1.3390E+01  7.2320E+00 -2.9839E+01 -2.6212E-01  1.4970E+00
            -1.9092E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3756.14398712070        NO. OF FUNC. EVALS.:  84
 CUMULATIVE NO. OF FUNC. EVALS.:       97
 NPARAMETR:  9.2371E-01  1.0814E+00  9.8446E-01  9.4247E-01  1.0140E+00  1.0159E+00  9.6930E-01  1.0599E+00  9.9678E-01  1.0053E+00
             1.1022E+00
 PARAMETER:  2.0641E-02  1.7825E-01  8.4335E-02  4.0748E-02  1.1391E-01  1.1579E-01  6.8818E-02  1.5819E-01  9.6771E-02  1.0525E-01
             1.9729E-01
 GRADIENT:  -4.6763E+01 -7.7946E+00 -4.4248E+00 -6.3139E+00 -1.9641E+01 -5.1257E+00  1.1429E+01 -1.5160E+01  2.7058E-01 -5.9107E+00
             2.8863E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3756.44988213926        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:      221
 NPARAMETR:  9.2296E-01  1.0817E+00  9.8212E-01  9.4383E-01  1.0172E+00  1.0167E+00  9.6416E-01  1.0726E+00  9.9535E-01  1.0164E+00
             1.1024E+00
 PARAMETER:  1.9831E-02  1.7856E-01  8.1962E-02  4.2191E-02  1.1706E-01  1.1653E-01  6.3506E-02  1.7011E-01  9.5343E-02  1.1628E-01
             1.9748E-01
 GRADIENT:  -8.4428E+01 -2.3349E+01 -8.6280E+00 -9.3805E+00 -1.9373E+01 -1.1158E+01  1.1057E+01 -1.4333E+01 -5.7017E-01 -4.6684E+00
             2.9899E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3760.56514549065        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      405
 NPARAMETR:  9.3437E-01  1.1880E+00  1.0534E+00  8.9967E-01  1.1132E+00  9.9366E-01  8.0131E-01  1.0772E+00  9.5818E-01  1.2273E+00
             1.1058E+00
 PARAMETER:  3.2122E-02  2.7229E-01  1.5200E-01 -5.7271E-03  2.0728E-01  9.3637E-02 -1.2151E-01  1.7436E-01  5.7276E-02  3.0478E-01
             2.0054E-01
 GRADIENT:  -6.1498E+01 -6.3460E+00  3.9606E+00  4.4402E+00 -1.2622E+01 -1.9042E+01  2.6173E+00 -2.2380E+01 -1.4658E+01  1.4402E+01
             2.3427E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3762.77971947771        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  9.5909E-01  1.1914E+00  1.0455E+00  8.9671E-01  1.1243E+00  1.0399E+00  7.7695E-01  1.0780E+00  1.0164E+00  1.1535E+00
             1.1058E+00
 PARAMETER:  5.8229E-02  2.7515E-01  1.4451E-01 -9.0207E-03  2.1714E-01  1.3911E-01 -1.5238E-01  1.7506E-01  1.1626E-01  2.4282E-01
             2.0054E-01
 GRADIENT:  -2.7152E+00 -8.6756E+00 -1.4817E+00  5.7067E+00  1.7684E+00  1.2886E+00 -6.5769E-01 -2.1898E+01 -1.6934E+00 -2.4885E+00
             2.5021E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3762.98740272324        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  9.6074E-01  1.2237E+00  1.0554E+00  8.7685E-01  1.1488E+00  1.0358E+00  7.6436E-01  1.0780E+00  1.0318E+00  1.1837E+00
             1.1058E+00
 PARAMETER:  5.9949E-02  3.0187E-01  1.5388E-01 -3.1425E-02  2.3874E-01  1.3513E-01 -1.6872E-01  1.7506E-01  1.3128E-01  2.6862E-01
             2.0054E-01
 GRADIENT:   7.5797E-01 -1.6149E-01  1.6097E-01 -5.5575E-01 -1.2060E-02 -2.8745E-01 -8.1157E-02 -2.3134E+01 -1.2405E-01 -1.3067E-01
             2.1730E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3763.84835209086        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      923
 NPARAMETR:  9.6032E-01  1.2231E+00  1.0549E+00  8.7735E-01  1.1483E+00  1.0365E+00  7.6475E-01  1.1198E+00  1.0321E+00  1.1839E+00
             1.1049E+00
 PARAMETER:  5.9508E-02  3.0141E-01  1.5345E-01 -3.0844E-02  2.3826E-01  1.3584E-01 -1.6820E-01  2.1319E-01  1.3156E-01  2.6885E-01
             1.9976E-01
 GRADIENT:   3.6013E+01  2.5198E+01 -1.4334E+00  5.7172E+00  8.4085E+00  6.4742E+00  8.3271E-01 -2.0935E+01  8.8652E-01  1.2702E+00
             2.4513E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3763.87417940741        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1096
 NPARAMETR:  9.6030E-01  1.2231E+00  1.0551E+00  8.7736E-01  1.1483E+00  1.0365E+00  7.6472E-01  1.1209E+00  1.0320E+00  1.1839E+00
             1.1046E+00
 PARAMETER:  5.9495E-02  3.0141E-01  1.5361E-01 -3.0836E-02  2.3827E-01  1.3583E-01 -1.6825E-01  2.1411E-01  1.3155E-01  2.6884E-01
             1.9950E-01
 GRADIENT:  -1.4204E-01 -6.4065E-01 -1.9954E+00 -2.6160E-01 -7.6335E-01 -3.5046E-02  1.0320E-01 -2.1020E+01  7.8582E-02 -2.7923E-02
             2.3704E+01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -3763.87417940741        NO. OF FUNC. EVALS.:  32
 CUMULATIVE NO. OF FUNC. EVALS.:     1128
 NPARAMETR:  9.6031E-01  1.2232E+00  1.0551E+00  8.7738E-01  1.1484E+00  1.0365E+00  7.6459E-01  1.1209E+00  1.0320E+00  1.1839E+00
             1.1046E+00
 PARAMETER:  5.9495E-02  3.0141E-01  1.5361E-01 -3.0836E-02  2.3827E-01  1.3583E-01 -1.6825E-01  2.1411E-01  1.3155E-01  2.6884E-01
             1.9950E-01
 GRADIENT:  -6.1878E+05 -5.1013E-01 -4.0280E+05 -2.2088E-01 -6.8678E-01 -2.5252E-02  9.8539E-02  5.7779E+05  7.1514E-02 -2.8406E-02
            -6.2033E+05
 NUMSIGDIG:         3.3         2.7         3.3         2.6         2.3         3.0         1.9         3.3         2.4         3.4
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1128
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.1426E-04 -3.5227E-02 -2.0015E-02  1.9607E-02 -2.4861E-02
 SE:             2.9884E-02  2.0272E-02  1.9206E-02  2.6737E-02  2.6045E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8627E-01  8.2264E-02  2.9735E-01  4.6335E-01  3.3980E-01

 ETASHRINKSD(%)  1.0000E-10  3.2086E+01  3.5658E+01  1.0428E+01  1.2747E+01
 ETASHRINKVR(%)  1.0000E-10  5.3877E+01  5.8601E+01  1.9769E+01  2.3869E+01
 EBVSHRINKSD(%)  3.0377E-01  3.2510E+01  4.6821E+01  1.1695E+01  1.1698E+01
 EBVSHRINKVR(%)  6.0662E-01  5.4451E+01  7.1719E+01  2.2023E+01  2.2027E+01
 RELATIVEINF(%)  9.9392E+01  1.7664E+01  2.3006E+01  3.8109E+01  3.4753E+01
 EPSSHRINKSD(%)  2.1307E+01
 EPSSHRINKVR(%)  3.8074E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3763.8741794074122     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2109.7848196390014     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3763.874       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.60E-01  1.22E+00  1.06E+00  8.77E-01  1.15E+00  1.04E+00  7.65E-01  1.12E+00  1.03E+00  1.18E+00  1.10E+00
 


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
+        3.35E+09
 
 TH 2
+        3.82E+04  5.73E+02
 
 TH 3
+        1.99E+09  2.26E+04  1.18E+09
 
 TH 4
+        1.61E+05  4.21E+02 -2.18E+09  4.02E+09
 
 TH 5
+        3.88E+03 -2.44E+02  2.18E+03  1.36E+02  4.69E+02
 
 TH 6
+       -3.08E+03 -1.62E+00 -1.82E+03  4.34E+00 -6.03E-01  1.84E+02
 
 TH 7
+       -3.33E+03  3.58E+00 -1.98E+03 -2.74E+09 -1.70E+00  1.26E+00  8.17E+01
 
 TH 8
+        1.30E+03 -1.53E+04  8.30E+04 -6.44E+04 -1.56E+03  1.23E+03  1.33E+03  5.37E+08
 
 TH 9
+        1.75E+05 -7.14E+00  1.04E+05 -2.60E+09  9.04E-01 -4.45E+00  2.64E+01 -7.02E+04  1.68E+09
 
 TH10
+        1.50E+02 -2.33E+01  9.41E+01 -1.11E+09 -1.65E+01  5.27E-01  2.70E+01 -6.06E+01 -3.61E-01  8.84E+01
 
 TH11
+        1.46E+09  1.66E+04 -9.05E+04  7.02E+04  1.68E+03 -1.34E+03 -1.45E+03 -3.93E+04  7.65E+04  6.99E+01  6.37E+08
 
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
 #CPUT: Total CPU Time in Seconds,       44.722
Stop Time:
Sat Sep 25 04:31:39 CDT 2021
