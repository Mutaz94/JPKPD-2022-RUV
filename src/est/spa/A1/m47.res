Wed Sep 29 12:09:34 CDT 2021
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
$DATA ../../../../data/spa/A1/dat47.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1222.30627780133        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4605E+02  2.8621E+01  3.0738E+01 -1.5095E+00  1.3720E+02  6.3708E+01 -5.4761E+01 -1.2861E+01 -8.0934E+01 -1.0094E+02
            -5.9546E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1408.31756480148        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0544E+00  8.0671E-01  8.4796E-01  1.2324E+00  7.2469E-01  1.0072E+00  1.3523E+00  9.1608E-01  1.2798E+00  1.2542E+00
             2.0525E+00
 PARAMETER:  1.5294E-01 -1.1480E-01 -6.4925E-02  3.0893E-01 -2.2201E-01  1.0716E-01  4.0183E-01  1.2351E-02  3.4673E-01  3.2652E-01
             8.1906E-01
 GRADIENT:   2.7658E+02  5.4801E+01  2.1555E+01  1.1627E+02 -4.8016E+01  3.5633E+01  2.8820E+00  9.8419E+00  3.9443E+01  2.7498E+01
             4.0375E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1418.18054155817        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      201
 NPARAMETR:  1.0470E+00  5.2014E-01  6.9297E-01  1.4328E+00  5.5003E-01  9.9464E-01  1.9617E+00  5.5971E-01  1.0697E+00  1.0364E+00
             1.9525E+00
 PARAMETER:  1.4594E-01 -5.5366E-01 -2.6677E-01  4.5966E-01 -4.9779E-01  9.4627E-02  7.7382E-01 -4.8033E-01  1.6738E-01  1.3576E-01
             7.6909E-01
 GRADIENT:   1.2090E+02  5.8055E+01  4.6363E+01  1.3445E+02 -8.1706E+01  2.3901E+01  9.1599E-01  4.0667E+00 -9.0763E+00  1.5973E+01
             1.5906E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1441.79992306880        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  9.8106E-01  2.7809E-01  5.6935E-01  1.4368E+00  4.6223E-01  8.9674E-01  2.4283E+00  1.3288E-01  9.6491E-01  8.1215E-01
             1.8017E+00
 PARAMETER:  8.0877E-02 -1.1798E+00 -4.6326E-01  4.6239E-01 -6.7168E-01 -8.9846E-03  9.8718E-01 -1.9183E+00  6.4277E-02 -1.0807E-01
             6.8873E-01
 GRADIENT:  -1.1559E+01  1.9422E+01  1.6342E+01  6.6233E+01 -3.0572E+00 -5.6231E+00  8.5251E+00 -3.2312E-01 -3.0008E+01 -1.8052E+01
            -2.5934E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1450.50273726920        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  9.8215E-01  1.1175E-01  4.4339E-01  1.4084E+00  3.6473E-01  9.0527E-01  3.0361E+00  1.0000E-02  1.0255E+00  7.8301E-01
             1.8825E+00
 PARAMETER:  8.1989E-02 -2.0915E+00 -7.1329E-01  4.4245E-01 -9.0859E-01  4.7873E-04  1.2106E+00 -4.6134E+00  1.2522E-01 -1.4461E-01
             7.3263E-01
 GRADIENT:   1.3191E+00  9.9932E-01  2.9383E+00  5.7626E+00 -4.4197E+00 -9.8723E-01 -4.2987E-01  0.0000E+00 -1.0874E-02  3.4497E-01
            -1.7766E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1450.68382991556        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  9.7889E-01  6.2319E-02  4.4548E-01  1.4271E+00  3.6059E-01  9.0538E-01  3.6806E+00  1.0000E-02  1.0152E+00  7.7922E-01
             1.8924E+00
 PARAMETER:  7.8661E-02 -2.6755E+00 -7.0861E-01  4.5567E-01 -9.2002E-01  5.9985E-04  1.4031E+00 -5.5360E+00  1.1510E-01 -1.4946E-01
             7.3785E-01
 GRADIENT:  -6.4679E-01  1.8477E-01  3.5678E+00  2.4562E+00 -5.2708E+00 -2.9982E-01 -5.2777E-01  0.0000E+00  4.9645E-01 -1.0350E+00
            -3.8206E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1450.78625540759        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  9.7630E-01  1.3688E-02  4.5958E-01  1.4551E+00  3.6429E-01  9.0445E-01  6.6226E+00  1.0000E-02  9.9751E-01  7.9471E-01
             1.9006E+00
 PARAMETER:  7.6017E-02 -4.1912E+00 -6.7744E-01  4.7508E-01 -9.0979E-01 -4.2589E-04  1.9905E+00 -7.9297E+00  9.7506E-02 -1.2978E-01
             7.4218E-01
 GRADIENT:  -4.9148E-01  2.6736E-02  1.1248E+00  1.5945E+00 -1.9161E+00 -7.6719E-03 -1.0230E-01  0.0000E+00 -3.2605E-01  2.0072E-01
             1.4988E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1450.81577934386        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1084            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7634E-01  1.0000E-02  4.5941E-01  1.4537E+00  3.6411E-01  9.0435E-01  9.7238E+00  1.0000E-02  9.9734E-01  7.9385E-01
             1.9009E+00
 PARAMETER:  7.6058E-02 -5.0206E+00 -6.7782E-01  4.7413E-01 -9.1029E-01 -5.3467E-04  2.3746E+00 -9.3141E+00  9.7337E-02 -1.3086E-01
             7.4230E-01
 GRADIENT:   9.1692E+01  0.0000E+00  1.8797E+01  1.7270E+02  1.0755E+02  7.4308E+00  3.0667E-01  0.0000E+00  8.3597E+00  1.1142E+00
             5.8930E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1450.81577934386        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  9.7630E-01  1.0000E-02  4.5933E-01  1.4538E+00  3.6412E-01  9.0423E-01  9.4956E+00  1.0000E-02  9.9683E-01  7.9391E-01
             1.9012E+00
 PARAMETER:  7.6058E-02 -5.0206E+00 -6.7782E-01  4.7413E-01 -9.1029E-01 -5.3467E-04  2.3746E+00 -9.3141E+00  9.7337E-02 -1.3086E-01
             7.4230E-01
 GRADIENT:   5.5014E-02  0.0000E+00  1.0593E-01 -6.9426E-02 -3.4528E-02  2.2871E-02  7.8024E-02  0.0000E+00  7.4445E-02 -6.5198E-03
            -3.8231E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1150
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0701E-05  6.3195E-04 -3.3934E-05 -6.0047E-03 -6.5163E-03
 SE:             2.9405E-02  1.6306E-03  2.5789E-04  2.8254E-02  2.4517E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9944E-01  6.9835E-01  8.9531E-01  8.3169E-01  7.9040E-01

 ETASHRINKSD(%)  1.4895E+00  9.4537E+01  9.9136E+01  5.3471E+00  1.7864E+01
 ETASHRINKVR(%)  2.9569E+00  9.9702E+01  9.9993E+01  1.0408E+01  3.2537E+01
 EBVSHRINKSD(%)  1.6520E+00  9.5464E+01  9.9190E+01  4.5193E+00  1.6649E+01
 EBVSHRINKVR(%)  3.2767E+00  9.9794E+01  9.9993E+01  8.8343E+00  3.0526E+01
 RELATIVEINF(%)  8.6937E+01  2.7283E-02  3.0356E-04  2.2921E+01  2.5085E+00
 EPSSHRINKSD(%)  3.8233E+01
 EPSSHRINKVR(%)  6.1849E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1450.8157793438647     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -715.66495278012655     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1450.816       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.00E-02  4.59E-01  1.45E+00  3.64E-01  9.04E-01  9.72E+00  1.00E-02  9.97E-01  7.94E-01  1.90E+00
 


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
+        1.38E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.66E+01  0.00E+00  3.48E+03
 
 TH 4
+       -2.00E+01  0.00E+00 -2.35E+02  4.84E+02
 
 TH 5
+        7.22E+01  0.00E+00 -5.32E+03 -1.09E+02  9.35E+03
 
 TH 6
+       -8.58E-01  0.00E+00  4.64E+00 -4.61E+00  3.32E+00  2.28E+02
 
 TH 7
+        1.06E-01  0.00E+00 -8.19E-01 -2.38E-01  3.35E-01  1.20E-01  3.86E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.30E+00  0.00E+00  3.22E+01 -4.35E+00  9.88E+00  1.54E+00  4.61E-01  0.00E+00  1.66E+02
 
 TH10
+       -3.56E+00  0.00E+00 -7.33E+01  3.15E+00 -6.07E+01 -1.32E+00 -1.17E-01  0.00E+00 -3.52E+00  1.52E+02
 
 TH11
+       -1.37E+01  0.00E+00 -6.54E+00 -6.36E+00  1.52E+01  2.93E+00 -6.21E+01  0.00E+00  5.43E+00  2.72E+01  7.02E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.178
Stop Time:
Wed Sep 29 12:09:56 CDT 2021
