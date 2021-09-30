Wed Sep 29 23:28:56 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1380.10725200683        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6920E+02 -2.1438E+01 -3.4113E+00  4.6237E+00  1.0074E+02  7.2961E+01 -7.9908E+00  1.0617E+01  2.6567E+01 -1.8011E+01
            -1.5280E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1838.55606293962        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0356E+00  1.0865E+00  1.1781E+00  1.0163E+00  1.0379E+00  7.7933E-01  8.6511E-01  8.5151E-01  7.8617E-01  7.0966E-01
             1.9956E+00
 PARAMETER:  1.3498E-01  1.8299E-01  2.6391E-01  1.1613E-01  1.3724E-01 -1.4932E-01 -4.4903E-02 -6.0749E-02 -1.4059E-01 -2.4297E-01
             7.9095E-01
 GRADIENT:   1.4155E+02  1.5292E+01  1.3391E+01 -7.2664E+00 -5.7233E-02 -4.1051E+01 -6.8457E+00  3.0776E-01 -2.2060E+01 -2.3203E+00
            -1.0701E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1845.98486001858        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0288E+00  1.0712E+00  6.6390E-01  9.9479E-01  7.9549E-01  8.1072E-01  1.0749E+00  5.9620E-01  8.9044E-01  3.1595E-01
             2.0095E+00
 PARAMETER:  1.2843E-01  1.6882E-01 -3.0962E-01  9.4777E-02 -1.2879E-01 -1.0983E-01  1.7222E-01 -4.1717E-01 -1.6043E-02 -1.0522E+00
             7.9790E-01
 GRADIENT:   9.6856E+01  4.2682E+00 -1.4292E+01  2.6671E+01  4.9076E+01 -2.9209E+01  1.0483E+01  2.9467E+00  1.0659E+01 -6.3278E-01
            -8.3256E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1849.79920811735        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0123E+00  9.4686E-01  3.6177E-01  1.0072E+00  5.3328E-01  8.7691E-01  1.0813E+00  4.6078E-01  8.0080E-01  1.6270E-01
             2.0631E+00
 PARAMETER:  1.1223E-01  4.5395E-02 -9.1675E-01  1.0721E-01 -5.2871E-01 -3.1354E-02  1.7817E-01 -6.7483E-01 -1.2214E-01 -1.7158E+00
             8.2419E-01
 GRADIENT:   1.4308E+01  1.4350E+01 -2.5700E+01  4.1632E+01  4.9297E+01 -3.3946E+00  9.2733E+00  2.0246E+00 -3.4537E+00  7.6945E-01
            -1.7628E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1855.56149988029        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      406
 NPARAMETR:  1.0450E+00  9.0514E-01  5.1466E-01  1.0625E+00  6.1795E-01  8.8044E-01  1.1111E+00  1.7899E-01  7.8558E-01  2.8931E-01
             2.1735E+00
 PARAMETER:  1.4400E-01  3.3089E-04 -5.6424E-01  1.6063E-01 -3.8135E-01 -2.7334E-02  2.0535E-01 -1.6204E+00 -1.4133E-01 -1.1402E+00
             8.7632E-01
 GRADIENT:   5.0598E+00  1.3229E+00  2.3311E+00 -2.9231E+00 -3.0916E+00 -3.7498E+00 -2.9351E+00  2.8313E-01 -1.6900E-01  9.4071E-02
             5.2419E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1856.18423852029        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  1.0414E+00  7.4750E-01  5.3070E-01  1.1446E+00  5.7188E-01  8.8878E-01  1.3275E+00  6.9031E-02  7.4229E-01  3.7228E-01
             2.1411E+00
 PARAMETER:  1.4054E-01 -1.9103E-01 -5.3356E-01  2.3504E-01 -4.5883E-01 -1.7911E-02  3.8333E-01 -2.5732E+00 -1.9802E-01 -8.8811E-01
             8.6132E-01
 GRADIENT:   4.6505E-01 -3.1914E-01 -9.2940E-01 -2.6037E+00  2.0064E+00  4.5535E-01  2.1024E-01  5.9738E-02 -2.9772E-01  1.4216E-01
             1.3638E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1856.25321163828        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      761
 NPARAMETR:  1.0400E+00  6.6656E-01  5.4520E-01  1.1887E+00  5.5437E-01  8.8634E-01  1.4572E+00  3.7380E-02  7.2459E-01  4.1245E-01
             2.1329E+00
 PARAMETER:  1.3920E-01 -3.0563E-01 -5.0660E-01  2.7286E-01 -4.8992E-01 -2.0649E-02  4.7655E-01 -3.1866E+00 -2.2215E-01 -7.8565E-01
             8.5747E-01
 GRADIENT:   3.0455E-02 -1.7603E-01 -2.6625E-01 -4.1119E-01 -1.1073E-01  1.1083E-02  6.3540E-02  1.9725E-02 -7.6247E-02  9.9908E-02
             1.1946E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1856.26259512901        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  1.0398E+00  6.6517E-01  5.5062E-01  1.1909E+00  5.5736E-01  8.8632E-01  1.4611E+00  1.0000E-02  7.2402E-01  4.1449E-01
             2.1304E+00
 PARAMETER:  1.3903E-01 -3.0772E-01 -4.9671E-01  2.7468E-01 -4.8454E-01 -2.0676E-02  4.7921E-01 -5.4974E+00 -2.2294E-01 -7.8070E-01
             8.5631E-01
 GRADIENT:  -3.3959E-02 -1.6558E-01  1.0908E-01 -1.4068E-01 -1.1075E-01  4.8602E-02 -7.1082E-02  0.0000E+00 -7.4092E-02 -6.9764E-02
            -4.9574E-01

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1856.26259512901        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      921
 NPARAMETR:  1.0398E+00  6.6517E-01  5.5062E-01  1.1909E+00  5.5736E-01  8.8632E-01  1.4611E+00  1.0000E-02  7.2402E-01  4.1449E-01
             2.1304E+00
 PARAMETER:  1.3903E-01 -3.0772E-01 -4.9671E-01  2.7468E-01 -4.8454E-01 -2.0676E-02  4.7921E-01 -5.4974E+00 -2.2294E-01 -7.8070E-01
             8.5631E-01
 GRADIENT:  -3.3959E-02 -1.6558E-01  1.0908E-01 -1.4068E-01 -1.1075E-01  4.8602E-02 -7.1082E-02  0.0000E+00 -7.4092E-02 -6.9764E-02
            -4.9574E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      921
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5137E-03  2.0203E-02 -2.8828E-04 -1.8514E-02  4.4750E-03
 SE:             2.9403E-02  2.1074E-02  2.0107E-04  2.4408E-02  1.3566E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5894E-01  3.3774E-01  1.5164E-01  4.4813E-01  7.4150E-01

 ETASHRINKSD(%)  1.4946E+00  2.9398E+01  9.9326E+01  1.8230E+01  5.4552E+01
 ETASHRINKVR(%)  2.9668E+00  5.0154E+01  9.9995E+01  3.3137E+01  7.9345E+01
 EBVSHRINKSD(%)  1.7165E+00  3.0074E+01  9.9272E+01  1.7979E+01  5.4354E+01
 EBVSHRINKVR(%)  3.4035E+00  5.1104E+01  9.9995E+01  3.2726E+01  7.9164E+01
 RELATIVEINF(%)  9.5987E+01  5.7544E+00  2.3893E-04  1.3488E+01  7.6368E-01
 EPSSHRINKSD(%)  2.6038E+01
 EPSSHRINKVR(%)  4.5296E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1856.2625951290095     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -937.32406192433677     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1856.263       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  6.65E-01  5.51E-01  1.19E+00  5.57E-01  8.86E-01  1.46E+00  1.00E-02  7.24E-01  4.14E-01  2.13E+00
 


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
+        1.26E+03
 
 TH 2
+       -3.40E+01  5.93E+02
 
 TH 3
+        2.18E+01  5.90E+02  2.04E+03
 
 TH 4
+       -3.13E+01  3.29E+02 -5.50E+02  1.02E+03
 
 TH 5
+        1.49E+01 -1.10E+03 -2.89E+03  4.57E+02  4.38E+03
 
 TH 6
+        2.45E+00 -5.81E+00  9.46E+00 -1.06E+01 -2.52E+00  2.35E+02
 
 TH 7
+        2.19E+00  4.41E+01 -1.14E+01 -7.62E+00 -1.68E+00  1.29E+00  2.79E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.38E+00 -2.11E+01 -1.12E+01 -5.38E+00  4.00E+01 -6.10E-01  1.73E+01  0.00E+00  1.72E+02
 
 TH10
+       -1.57E+00 -1.18E+00 -9.08E+01 -1.60E+01  8.89E+01  2.29E-02  8.04E+00  0.00E+00  4.25E-01  6.32E+01
 
 TH11
+       -1.51E+01 -9.96E+00 -3.70E+01 -1.89E+01  9.92E+00  3.26E+00  4.50E+00  0.00E+00  1.57E+01  2.80E+01  1.04E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.365
Stop Time:
Wed Sep 29 23:29:22 CDT 2021
