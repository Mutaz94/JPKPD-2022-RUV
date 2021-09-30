Wed Sep 29 18:52:19 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1701.06330242860        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1647E+02 -4.8869E+01 -1.5880E+01 -5.6346E+01 -3.5151E+01  5.9660E+01 -1.1378E+00  1.8713E+01  8.4023E+00  3.1544E+01
             7.9644E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1713.10469595106        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.9723E-01  1.1278E+00  1.2154E+00  1.0135E+00  1.1585E+00  8.7817E-01  1.0211E+00  8.8852E-01  9.5666E-01  8.4641E-01
             9.8074E-01
 PARAMETER:  9.7226E-02  2.2024E-01  2.9509E-01  1.1341E-01  2.4714E-01 -2.9911E-02  1.2088E-01 -1.8194E-02  5.5698E-02 -6.6751E-02
             8.0551E-02
 GRADIENT:   1.2502E+01  1.4050E+01  1.9413E+01 -4.8561E-01 -5.6851E+00 -2.3930E+01 -3.8843E+00 -1.3823E+00 -8.6116E+00 -1.4887E+01
            -1.7263E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1714.72287075076        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      347
 NPARAMETR:  9.9552E-01  1.1427E+00  1.1686E+00  1.0030E+00  1.1739E+00  9.2923E-01  1.0543E+00  5.7130E-01  9.5565E-01  9.5689E-01
             1.0007E+00
 PARAMETER:  9.5507E-02  2.3336E-01  2.5584E-01  1.0295E-01  2.6032E-01  2.6598E-02  1.5292E-01 -4.5984E-01  5.4641E-02  5.5938E-02
             1.0075E-01
 GRADIENT:   6.4379E+00  1.1052E+01  2.0150E+01  5.5688E-01  4.7530E+00 -2.4472E-01  1.2089E-01 -3.3946E+00 -9.1941E+00 -4.4503E+00
            -5.7983E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1716.92996164967        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  9.9528E-01  1.2223E+00  9.0241E-01  9.3930E-01  1.0672E+00  9.3644E-01  1.0078E+00  3.4621E-01  1.0298E+00  8.6936E-01
             9.9556E-01
 PARAMETER:  9.5268E-02  3.0072E-01 -2.6897E-03  3.7379E-02  1.6508E-01  3.4334E-02  1.0773E-01 -9.6071E-01  1.2935E-01 -3.9997E-02
             9.5549E-02
 GRADIENT:   9.5447E-01  4.9951E+00  2.3441E+00  2.6581E+00 -1.0925E+01  1.8618E+00 -8.6079E-02  8.7506E-01  1.2791E+00  2.3015E+00
             1.0999E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1717.38169972467        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  9.9596E-01  1.4856E+00  7.9438E-01  7.6895E-01  1.1739E+00  9.3287E-01  8.7752E-01  1.5891E-01  1.1825E+00  9.0685E-01
             9.9744E-01
 PARAMETER:  9.5950E-02  4.9584E-01 -1.3020E-01 -1.6273E-01  2.6036E-01  3.0508E-02 -3.0659E-02 -1.7394E+00  2.6764E-01  2.2189E-03
             9.7440E-02
 GRADIENT:   4.1980E-02  2.2653E+00 -1.7868E+00  3.9034E+00  1.0404E+00 -1.1122E-02 -3.8605E-01  2.0515E-01 -4.0474E-01  1.0039E-01
            -3.7928E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1717.46335857880        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.9613E-01  1.5122E+00  7.8727E-01  7.4853E-01  1.1892E+00  9.3299E-01  8.6778E-01  4.9635E-02  1.2099E+00  9.1596E-01
             9.9805E-01
 PARAMETER:  9.6119E-02  5.1359E-01 -1.3919E-01 -1.8965E-01  2.7326E-01  3.0644E-02 -4.1818E-02 -2.9031E+00  2.9057E-01  1.2217E-02
             9.8050E-02
 GRADIENT:   4.7461E-01 -1.3357E+00 -5.4270E-01  2.1578E-01  7.3663E-01  5.0472E-02  6.0411E-02  1.9779E-02  8.7547E-02 -8.0438E-02
             1.1682E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1717.46762732363        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.9580E-01  1.5159E+00  7.8762E-01  7.4655E-01  1.1909E+00  9.3281E-01  8.6586E-01  2.0574E-02  1.2121E+00  9.1894E-01
             9.9783E-01
 PARAMETER:  9.5795E-02  5.1603E-01 -1.3874E-01 -1.9229E-01  2.7470E-01  3.0443E-02 -4.4029E-02 -3.7837E+00  2.9239E-01  1.5469E-02
             9.7831E-02
 GRADIENT:  -3.3135E-01 -9.7825E-02  2.0385E-01  1.3435E-01 -5.6885E-01 -2.4495E-02 -1.1520E-02  3.3598E-03 -8.5785E-02  7.3645E-02
            -2.9639E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1717.47093445029        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  9.9637E-01  1.5147E+00  7.8792E-01  7.4710E-01  1.1907E+00  9.3301E-01  8.6627E-01  1.0000E-02  1.2125E+00  9.1833E-01
             9.9783E-01
 PARAMETER:  9.6363E-02  5.1525E-01 -1.3836E-01 -1.9155E-01  2.7457E-01  3.0659E-02 -4.3563E-02 -4.8167E+00  2.9267E-01  1.4804E-02
             9.7831E-02
 GRADIENT:   1.1156E+00 -6.7812E-01  8.4266E-03  1.1024E-01 -6.9524E-02  6.0434E-02  3.8849E-03  0.0000E+00  7.2508E-02  2.5284E-03
            -2.3190E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1717.47093445029        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1264
 NPARAMETR:  9.9637E-01  1.5147E+00  7.8792E-01  7.4710E-01  1.1907E+00  9.3301E-01  8.6627E-01  1.0000E-02  1.2125E+00  9.1833E-01
             9.9783E-01
 PARAMETER:  9.6363E-02  5.1525E-01 -1.3836E-01 -1.9155E-01  2.7457E-01  3.0659E-02 -4.3563E-02 -4.8167E+00  2.9267E-01  1.4804E-02
             9.7831E-02
 GRADIENT:   1.1156E+00 -6.7812E-01  8.4266E-03  1.1024E-01 -6.9524E-02  6.0434E-02  3.8849E-03  0.0000E+00  7.2508E-02  2.5284E-03
            -2.3190E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1264
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.4251E-05 -2.2035E-02 -3.2406E-04  1.4162E-02 -3.3056E-02
 SE:             2.9815E-02  2.2303E-02  1.2684E-04  2.3377E-02  2.1816E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9828E-01  3.2315E-01  1.0622E-02  5.4464E-01  1.2973E-01

 ETASHRINKSD(%)  1.1485E-01  2.5283E+01  9.9575E+01  2.1683E+01  2.6913E+01
 ETASHRINKVR(%)  2.2956E-01  4.4174E+01  9.9998E+01  3.8664E+01  4.6583E+01
 EBVSHRINKSD(%)  4.8586E-01  2.4362E+01  9.9595E+01  2.2674E+01  2.5910E+01
 EBVSHRINKVR(%)  9.6936E-01  4.2789E+01  9.9998E+01  4.0207E+01  4.5107E+01
 RELATIVEINF(%)  9.8850E+01  3.0622E+00  2.2423E-04  3.6469E+00  8.9979E+00
 EPSSHRINKSD(%)  4.2561E+01
 EPSSHRINKVR(%)  6.7007E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1717.4709344502874     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -982.32010788654918     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.40
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1717.471       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.96E-01  1.51E+00  7.88E-01  7.47E-01  1.19E+00  9.33E-01  8.66E-01  1.00E-02  1.21E+00  9.18E-01  9.98E-01
 


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
+        1.27E+03
 
 TH 2
+       -6.07E+00  3.87E+02
 
 TH 3
+        1.25E+01  1.55E+02  3.28E+02
 
 TH 4
+       -1.37E+01  3.24E+02 -1.99E+02  8.19E+02
 
 TH 5
+       -3.95E+00 -1.99E+02 -3.00E+02  2.06E+02  5.12E+02
 
 TH 6
+        6.39E-01 -9.89E-01  2.89E+00 -3.94E+00 -1.38E+00  2.25E+02
 
 TH 7
+        7.73E-01  1.73E+01 -2.59E+00 -1.47E+01 -9.79E+00 -1.39E-01  9.78E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.26E+00 -2.18E+01 -2.77E+01  4.49E+01 -7.64E-01 -6.12E-01  1.75E+01  0.00E+00  6.13E+01
 
 TH10
+        1.61E-01 -9.80E+00 -3.10E+01 -1.04E+01 -6.10E+01  3.75E-01  1.35E+01  0.00E+00  6.20E+00  8.23E+01
 
 TH11
+       -8.97E+00 -2.05E+01 -3.12E+01  1.62E+00  1.49E+00  2.71E+00  1.05E+01  0.00E+00  6.05E+00  2.28E+01  2.24E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.450
Stop Time:
Wed Sep 29 18:52:43 CDT 2021
