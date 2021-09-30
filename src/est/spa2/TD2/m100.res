Thu Sep 30 08:29:00 CDT 2021
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
$DATA ../../../../data/spa2/TD2/dat100.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2137.22187023255        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0814E+02 -1.8033E+01  8.0251E+01  2.2154E+01  5.7337E+01  3.1174E+01  2.9780E+00 -4.2977E+02 -7.4555E+01 -5.6686E+00
            -1.0386E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2335.02786604107        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      163             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3784E+00  1.0796E+00  9.7609E-01  1.0383E+00  9.9474E-01  1.5851E+00  9.9411E-01  1.7587E+00  9.0273E-01  9.9000E-01
             1.0708E+00
 PARAMETER:  4.2090E-01  1.7657E-01  7.5802E-02  1.3758E-01  9.4731E-02  5.6064E-01  9.4092E-02  6.6456E-01 -2.3288E-03  8.9945E-02
             1.6843E-01
 GRADIENT:   1.9046E+03  1.3697E+02 -3.6691E+01  2.1987E+02  9.3345E+00  2.5085E+02  6.2092E+00  3.0871E+01  2.1311E+00 -2.4582E-01
             3.2239E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2336.82819491277        NO. OF FUNC. EVALS.: 113
 CUMULATIVE NO. OF FUNC. EVALS.:      276
 NPARAMETR:  1.3769E+00  1.0795E+00  1.0442E+00  1.0382E+00  1.0073E+00  1.4746E+00  9.9407E-01  1.7582E+00  8.9511E-01  9.9651E-01
             1.0709E+00
 PARAMETER:  4.1986E-01  1.7650E-01  1.4324E-01  1.3753E-01  1.0731E-01  4.8841E-01  9.4055E-02  6.6430E-01 -1.0806E-02  9.6501E-02
             1.6849E-01
 GRADIENT:   3.9127E+02  5.7271E+01 -2.7426E+01  1.1627E+02 -2.3681E+01 -3.0196E+01  1.3787E+00  2.1641E+01 -3.4370E+00 -2.2962E+00
             2.9677E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2350.58002793172        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      453
 NPARAMETR:  1.3737E+00  1.0795E+00  1.9062E+00  1.0382E+00  1.2902E+00  1.5662E+00  9.9408E-01  1.7582E+00  8.6431E-01  1.2207E+00
             1.0709E+00
 PARAMETER:  4.1749E-01  1.7646E-01  7.4514E-01  1.3748E-01  3.5477E-01  5.4867E-01  9.4066E-02  6.6429E-01 -4.5820E-02  2.9944E-01
             1.6846E-01
 GRADIENT:   3.4976E+02  3.8595E+01  3.5301E+00  6.9756E+01  3.5495E+00 -3.5078E+00 -1.1088E+00 -8.8453E-01  5.3619E-01  1.9861E-01
             1.0942E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2410.71472678087        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  9.6819E-01  1.0738E+00  3.0911E+00  1.0316E+00  1.4803E+00  8.0861E-01  9.9503E-01  1.7644E+00  8.6053E-01  1.5384E+00
             1.0683E+00
 PARAMETER:  6.7671E-02  1.7123E-01  1.2285E+00  1.3113E-01  4.9225E-01 -1.1244E-01  9.5016E-02  6.6782E-01 -5.0211E-02  5.3075E-01
             1.6608E-01
 GRADIENT:   3.1662E+01  3.5688E+01  3.7052E+01  1.3613E+01 -1.5560E+01 -1.2248E+02 -1.1533E+00 -2.2818E+01  6.5855E-01  1.6514E+01
            -1.1189E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2435.67622586547        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  9.6505E-01  1.0740E+00  1.8840E+00  1.0317E+00  1.2782E+00  1.0313E+00  9.9519E-01  1.7664E+00  8.6199E-01  1.2188E+00
             1.0680E+00
 PARAMETER:  6.4423E-02  1.7138E-01  7.3342E-01  1.3119E-01  3.4541E-01  1.3078E-01  9.5182E-02  6.6893E-01 -4.8513E-02  2.9789E-01
             1.6578E-01
 GRADIENT:   1.2127E+01  3.1428E+01  3.1110E+00  5.5830E+01 -9.8567E-01 -4.0720E-01 -1.2979E+00  5.8159E-01  2.4481E-01  9.4675E-01
             1.1156E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2435.79645853255        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      995
 NPARAMETR:  9.6069E-01  1.0740E+00  1.8433E+00  1.0316E+00  1.2693E+00  1.0318E+00  1.0015E+00  1.7658E+00  8.6171E-01  1.2066E+00
             1.0638E+00
 PARAMETER:  5.9898E-02  1.7142E-01  7.1154E-01  1.3111E-01  3.3848E-01  1.3130E-01  1.0150E-01  6.6860E-01 -4.8831E-02  2.8779E-01
             1.6188E-01
 GRADIENT:   2.7272E+00  3.1287E+01  7.9652E-01  5.7918E+01 -3.8170E-01 -1.6722E-01 -8.6314E-01  1.7060E+00  4.2332E-01  1.7264E-01
             7.2111E+00

0ITERATION NO.:   31    OBJECTIVE VALUE:  -2435.79645853255        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     1025
 NPARAMETR:  9.6071E-01  1.0741E+00  1.8436E+00  1.0316E+00  1.2693E+00  1.0319E+00  1.0025E+00  1.7658E+00  8.6153E-01  1.2065E+00
             1.0639E+00
 PARAMETER:  5.9898E-02  1.7142E-01  7.1154E-01  1.3111E-01  3.3848E-01  1.3130E-01  1.0150E-01  6.6860E-01 -4.8831E-02  2.8779E-01
             1.6188E-01
 GRADIENT:  -6.9219E+04 -4.0357E+04 -9.7056E+03 -5.2757E+04  4.0301E-01 -3.3302E-01 -8.4797E-01 -2.4693E+02  3.7331E-01  2.4052E+04
            -4.2764E+04
 NUMSIGDIG:         2.3         2.3         2.3         2.3         7.3         1.9         0.7         4.2         1.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1025
 NO. OF SIG. DIGITS IN FINAL EST.:  0.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.0774E-05 -3.1491E-02 -3.9035E-02 -1.4204E-02 -3.3410E-02
 SE:             2.9905E-02  1.9510E-02  1.7569E-02  2.4103E-02  2.3621E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9758E-01  1.0650E-01  2.6299E-02  5.5566E-01  1.5724E-01

 ETASHRINKSD(%)  1.0000E-10  3.4641E+01  4.1141E+01  1.9253E+01  2.0867E+01
 ETASHRINKVR(%)  1.0000E-10  5.7281E+01  6.5356E+01  3.4800E+01  3.7380E+01
 EBVSHRINKSD(%)  3.4630E-01  3.4824E+01  4.4596E+01  2.0321E+01  1.7230E+01
 EBVSHRINKVR(%)  6.9139E-01  5.7520E+01  6.9304E+01  3.6513E+01  3.1491E+01
 RELATIVEINF(%)  9.9262E+01  5.2385E+00  1.0740E+01  8.1932E+00  3.2236E+01
 EPSSHRINKSD(%)  3.0345E+01
 EPSSHRINKVR(%)  5.1482E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2435.7964585325535     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1333.0702186869464     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2435.796       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  1.07E+00  1.84E+00  1.03E+00  1.27E+00  1.03E+00  1.00E+00  1.77E+00  8.62E-01  1.21E+00  1.06E+00
 


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
+        3.75E+07
 
 TH 2
+       -9.78E+06  1.02E+07
 
 TH 3
+       -1.38E+06  3.00E+01  1.00E+05
 
 TH 4
+       -5.39E+00  6.95E+06 -9.77E+05  1.89E+07
 
 TH 5
+        5.49E+01 -2.19E+06 -3.06E+05 -2.98E+06  1.88E+06
 
 TH 6
+       -2.34E+02 -1.23E+02 -1.68E+01 -1.67E+02 -6.86E-01  1.85E+02
 
 TH 7
+       -1.77E+07 -9.24E+06 -6.16E+01 -6.05E+02 -4.60E-01  1.75E-01  1.67E+07
 
 TH 8
+        1.54E+06 -6.33E+00  2.67E+03  1.05E+01 -3.45E+05 -2.82E-01 -6.51E+01  2.48E+05
 
 TH 9
+        2.14E+03  1.09E+03  1.57E+02  1.50E+03  3.96E+00  1.66E-01  4.23E+01 -1.68E+06  1.12E+02
 
 TH10
+       -5.19E+06 -1.02E+01 -3.79E+05  6.81E+00 -3.54E+01  6.47E+01  2.33E+02 -1.01E+04 -5.89E+02  1.44E+06
 
 TH11
+       -7.32E+00 -5.46E+06  3.49E+02 -8.92E+00 -6.49E+00 -1.29E+02 -4.61E+02  2.05E+04  1.20E+03 -2.89E+06  5.84E+06
 
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
 #CPUT: Total CPU Time in Seconds,       30.460
Stop Time:
Thu Sep 30 08:29:32 CDT 2021
