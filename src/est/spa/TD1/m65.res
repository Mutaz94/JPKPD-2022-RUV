Sat Sep 25 13:01:42 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat65.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.86104202231        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.9893E+01 -5.8630E+01 -5.9500E+01 -2.1813E+00  1.1754E+02  4.4000E+01 -6.0479E+00  6.2069E+00 -7.4366E+00 -1.7201E+01
            -3.6474E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1700.47445434357        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      133
 NPARAMETR:  1.0474E+00  1.0304E+00  1.0529E+00  1.0047E+00  9.5031E-01  8.3537E-01  1.0133E+00  9.6912E-01  1.0459E+00  1.0241E+00
             1.0139E+00
 PARAMETER:  1.4631E-01  1.2999E-01  1.5159E-01  1.0468E-01  4.9035E-02 -7.9884E-02  1.1317E-01  6.8635E-02  1.4492E-01  1.2386E-01
             1.1382E-01
 GRADIENT:   5.2122E+01 -2.6096E-01 -3.1434E+00  4.1245E+00  9.7273E-01 -2.4873E+01  5.0477E-01  3.3035E+00  9.6908E-01 -2.2116E+00
             8.5700E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1701.07268285780        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      310
 NPARAMETR:  1.0438E+00  1.0609E+00  1.0177E+00  9.8569E-01  9.5538E-01  8.4326E-01  9.2857E-01  7.7461E-01  1.0851E+00  1.0907E+00
             1.0090E+00
 PARAMETER:  1.4283E-01  1.5908E-01  1.1756E-01  8.5585E-02  5.4355E-02 -7.0479E-02  2.5895E-02 -1.5540E-01  1.8172E-01  1.8682E-01
             1.0901E-01
 GRADIENT:   3.9506E+01  6.5783E-01 -1.3620E+00  7.5672E+00  1.1959E+00 -2.0701E+01 -1.9920E+00  2.7684E-01  2.1127E+00  3.7351E+00
            -8.1576E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1702.27611703601        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      487
 NPARAMETR:  1.0314E+00  1.1895E+00  8.5081E-01  8.9577E-01  9.3586E-01  8.8900E-01  9.4679E-01  5.7252E-01  1.1167E+00  1.0223E+00
             1.0106E+00
 PARAMETER:  1.3092E-01  2.7352E-01 -6.1571E-02 -1.0074E-02  3.3708E-02 -1.7661E-02  4.5327E-02 -4.5770E-01  2.1035E-01  1.2208E-01
             1.1051E-01
 GRADIENT:  -1.4163E+00  3.4301E+00 -2.5254E-01  4.9699E+00 -2.3966E+00  7.8729E-01 -1.0739E-01  3.8552E-01 -4.1296E-01  2.5962E-01
             9.8055E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1702.66652202481        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  1.0339E+00  1.4529E+00  6.9203E-01  7.2018E-01  9.9770E-01  8.8521E-01  8.5443E-01  3.1438E-01  1.2873E+00  1.0325E+00
             1.0136E+00
 PARAMETER:  1.3329E-01  4.7353E-01 -2.6813E-01 -2.2825E-01  9.7701E-02 -2.1930E-02 -5.7321E-02 -1.0572E+00  3.5252E-01  1.3194E-01
             1.1354E-01
 GRADIENT:   2.5517E+00  2.4628E+00 -9.1071E-01  2.8666E+00  1.5607E+00 -1.2526E+00  8.6739E-01  1.4695E-01 -5.6754E-01 -1.7996E-01
             4.4680E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1702.73774754201        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      838
 NPARAMETR:  1.0330E+00  1.5599E+00  6.3755E-01  6.4664E-01  1.0309E+00  8.8806E-01  8.1286E-01  2.0986E-01  1.3901E+00  1.0482E+00
             1.0139E+00
 PARAMETER:  1.3249E-01  5.4462E-01 -3.5013E-01 -3.3597E-01  1.3041E-01 -1.8713E-02 -1.0720E-01 -1.4613E+00  4.2936E-01  1.4704E-01
             1.1379E-01
 GRADIENT:   1.9395E-02 -2.5093E-01 -1.7018E-01 -3.3835E-01 -1.8222E-01 -1.6153E-02  1.0036E-01  7.6964E-02  7.4048E-02  9.8920E-02
             1.3689E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1702.76934172471        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1016
 NPARAMETR:  1.0332E+00  1.5289E+00  6.4500E-01  6.6758E-01  1.0164E+00  8.8785E-01  8.2273E-01  6.8814E-02  1.3597E+00  1.0412E+00
             1.0134E+00
 PARAMETER:  1.3266E-01  5.2457E-01 -3.3851E-01 -3.0409E-01  1.1623E-01 -1.8952E-02 -9.5124E-02 -2.5763E+00  4.0727E-01  1.4040E-01
             1.1335E-01
 GRADIENT:   4.3224E-01  1.0022E+00 -7.4104E-02  8.4215E-01 -4.4866E-01 -1.3080E-01 -8.3877E-02  6.6471E-03  1.0434E-01  1.6573E-01
             5.5681E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1702.77422847084        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1191
 NPARAMETR:  1.0330E+00  1.5438E+00  6.3838E-01  6.5705E-01  1.0222E+00  8.8813E-01  8.1880E-01  1.2112E-02  1.3734E+00  1.0431E+00
             1.0135E+00
 PARAMETER:  1.3251E-01  5.3427E-01 -3.4882E-01 -3.2000E-01  1.2195E-01 -1.8632E-02 -9.9912E-02 -4.3136E+00  4.1727E-01  1.4219E-01
             1.1343E-01
 GRADIENT:   2.2712E-02  4.1806E-02 -3.9701E-02  6.9273E-02 -8.3983E-03 -2.7942E-03  1.0770E-02  2.1180E-04  1.2195E-02  1.8997E-02
             1.5784E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1702.77426078208        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  1.0331E+00  1.5441E+00  6.3841E-01  6.5651E-01  1.0225E+00  8.8815E-01  8.1856E-01  1.0000E-02  1.3739E+00  1.0431E+00
             1.0135E+00
 PARAMETER:  1.3251E-01  5.3464E-01 -3.4895E-01 -3.2059E-01  1.2225E-01 -1.8628E-02 -1.0020E-01 -4.5719E+00  4.1770E-01  1.4234E-01
             1.1343E-01
 GRADIENT:  -2.4400E-02  7.6253E-02 -1.1654E-02  3.8685E-02 -7.9824E-03 -9.9155E-04  2.4776E-04  0.0000E+00  2.7226E-03  5.7822E-03
             4.0041E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1254
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1424E-05 -3.0988E-02 -3.0336E-04  2.1455E-02 -3.2101E-02
 SE:             2.9799E-02  2.1655E-02  1.1791E-04  2.3886E-02  2.3580E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9969E-01  1.5243E-01  1.0089E-02  3.6906E-01  1.7340E-01

 ETASHRINKSD(%)  1.7083E-01  2.7453E+01  9.9605E+01  1.9980E+01  2.1004E+01
 ETASHRINKVR(%)  3.4136E-01  4.7369E+01  9.9998E+01  3.5968E+01  3.7596E+01
 EBVSHRINKSD(%)  5.5689E-01  2.6654E+01  9.9667E+01  2.1237E+01  1.9286E+01
 EBVSHRINKVR(%)  1.1107E+00  4.6203E+01  9.9999E+01  3.7963E+01  3.4852E+01
 RELATIVEINF(%)  9.8787E+01  3.2059E+00  1.5459E-04  4.1784E+00  1.2212E+01
 EPSSHRINKSD(%)  4.4536E+01
 EPSSHRINKVR(%)  6.9238E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1702.7742607820846     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -967.62343421834646     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.44
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1702.774       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.54E+00  6.38E-01  6.57E-01  1.02E+00  8.88E-01  8.19E-01  1.00E-02  1.37E+00  1.04E+00  1.01E+00
 


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
+        1.30E+03
 
 TH 2
+       -7.21E+00  3.98E+02
 
 TH 3
+        1.23E+01  1.56E+02  3.38E+02
 
 TH 4
+       -1.66E+01  3.39E+02 -2.11E+02  8.52E+02
 
 TH 5
+       -3.95E+00 -2.25E+02 -3.37E+02  2.20E+02  6.33E+02
 
 TH 6
+       -3.34E+00 -1.19E+00  1.96E+00 -5.22E+00  1.65E+00  2.49E+02
 
 TH 7
+        1.76E+00  7.77E+00 -5.85E+00 -1.11E+01 -1.43E+01 -4.03E-01  1.01E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.56E+00 -2.21E+01 -2.69E+01  4.71E+01  4.85E-01 -4.98E-02  1.88E+01  0.00E+00  4.94E+01
 
 TH10
+        1.02E+00 -1.51E+01 -3.97E+01 -4.99E+00 -5.31E+01 -1.85E+00  1.52E+01  0.00E+00  5.92E+00  8.31E+01
 
 TH11
+       -9.09E+00 -1.77E+01 -2.91E+01  1.53E+00  2.86E+00  4.85E+00  9.21E+00  0.00E+00  4.12E+00  1.80E+01  2.08E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.433
Stop Time:
Sat Sep 25 13:02:05 CDT 2021
