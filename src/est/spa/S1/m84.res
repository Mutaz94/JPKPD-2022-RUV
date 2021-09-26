Sat Sep 25 10:08:25 CDT 2021
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
$DATA ../../../../data/spa/S1/dat84.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1708.34115378128        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.0585E+01 -5.6862E+01 -3.3986E+01 -5.6933E+01  1.5059E+01  4.2306E+01 -9.3027E+00  1.2997E+01 -1.2382E+01  1.3654E+01
             1.9545E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1715.51223221872        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0591E+00  1.0537E+00  1.2316E+00  9.8344E-01  1.1076E+00  8.3581E-01  1.0825E+00  8.6906E-01  1.0943E+00  9.2808E-01
             1.0361E+00
 PARAMETER:  1.5737E-01  1.5229E-01  3.0831E-01  8.3304E-02  2.0224E-01 -7.9359E-02  1.7929E-01 -4.0344E-02  1.9009E-01  2.5368E-02
             1.3544E-01
 GRADIENT:   1.1463E+02 -2.1580E+01  1.6652E+01 -4.3254E+01  9.0277E+00 -2.0246E+01  2.9698E+00 -2.5430E+00  4.3363E+00 -1.6300E+01
             3.0134E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1717.67312629506        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0437E+00  8.3958E-01  1.1011E+00  1.1560E+00  9.8731E-01  8.4446E-01  1.0850E+00  3.6277E-01  9.6976E-01  1.0013E+00
             1.0369E+00
 PARAMETER:  1.4276E-01 -7.4854E-02  1.9635E-01  2.4501E-01  8.7232E-02 -6.9053E-02  1.8159E-01 -9.1400E-01  6.9297E-02  1.0127E-01
             1.3622E-01
 GRADIENT:   6.2266E+01 -6.2334E+00 -1.7808E+01  3.5565E+01  3.5977E+01 -1.4528E+01 -7.1142E+00 -4.5806E-01 -2.9250E+00 -5.8683E+00
             5.2426E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1718.94694537847        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0291E+00  1.0125E+00  1.0064E+00  1.0294E+00  1.0039E+00  8.7518E-01  1.1551E+00  4.4287E-01  1.0080E+00  9.7772E-01
             1.0058E+00
 PARAMETER:  1.2873E-01  1.1244E-01  1.0642E-01  1.2902E-01  1.0392E-01 -3.3325E-02  2.4415E-01 -7.1449E-01  1.0797E-01  7.7468E-02
             1.0577E-01
 GRADIENT:   1.5447E+01 -6.1363E+00 -6.2724E+00  1.3886E-01  1.1402E+01 -1.3360E+00  6.7248E-01  3.9079E-01 -1.1441E+00  4.6478E-01
            -2.5617E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1719.39801727891        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      409
 NPARAMETR:  1.0431E+00  1.0401E+00  9.7416E-01  1.0205E+00  9.9262E-01  8.8224E-01  1.1404E+00  2.8993E-01  1.0279E+00  9.7242E-01
             1.0164E+00
 PARAMETER:  1.4217E-01  1.3927E-01  7.3820E-02  1.2028E-01  9.2595E-02 -2.5293E-02  2.3139E-01 -1.1381E+00  1.2748E-01  7.2037E-02
             1.1629E-01
 GRADIENT:  -2.4859E+00  7.1571E-01 -1.2453E+00  2.6059E+00  5.9550E-01 -7.1860E-01  2.8736E-01  9.2308E-02  5.4132E-01  4.1313E-01
             6.6744E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1719.45551860973        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.0444E+00  1.1289E+00  9.2767E-01  9.6185E-01  1.0119E+00  8.8382E-01  1.0797E+00  7.1798E-02  1.0688E+00  9.7371E-01
             1.0169E+00
 PARAMETER:  1.4342E-01  2.2122E-01  2.4921E-02  6.1108E-02  1.1187E-01 -2.3502E-02  1.7671E-01 -2.5339E+00  1.6652E-01  7.3358E-02
             1.1674E-01
 GRADIENT:  -7.5609E-01  6.8841E-01  1.2598E-01  7.1193E-01 -6.7822E-01 -3.1242E-01 -1.8103E-02  6.3286E-03  1.1741E-01  1.7870E-01
             2.2222E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1719.45809391445        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      761
 NPARAMETR:  1.0447E+00  1.1479E+00  9.2091E-01  9.4886E-01  1.0184E+00  8.8447E-01  1.0677E+00  4.3931E-02  1.0787E+00  9.7465E-01
             1.0165E+00
 PARAMETER:  1.4371E-01  2.3795E-01  1.7610E-02  4.7504E-02  1.1822E-01 -2.2770E-02  1.6547E-01 -3.0251E+00  1.7580E-01  7.4322E-02
             1.1638E-01
 GRADIENT:  -1.9331E-01 -1.9044E-02  1.7112E-01 -3.6012E-01 -3.9998E-01 -7.6917E-02 -4.1818E-03  2.5346E-03  3.9452E-02  1.2681E-01
             2.1206E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1719.45969189444        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  1.0447E+00  1.1391E+00  9.2285E-01  9.5464E-01  1.0152E+00  8.8460E-01  1.0741E+00  1.0000E-02  1.0735E+00  9.7319E-01
             1.0165E+00
 PARAMETER:  1.4375E-01  2.3025E-01  1.9716E-02  5.3579E-02  1.1507E-01 -2.2616E-02  1.7149E-01 -4.9227E+00  1.7089E-01  7.2822E-02
             1.1634E-01
 GRADIENT:   3.8074E-03  2.0713E-02 -1.8468E-02  4.1972E-02  5.9780E-03  1.2230E-03  1.5913E-03  0.0000E+00  1.3348E-03  2.7724E-03
             7.7675E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1719.45969290825        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      994
 NPARAMETR:  1.0447E+00  1.1393E+00  9.2291E-01  9.5450E-01  1.0153E+00  8.8460E-01  1.0739E+00  1.0000E-02  1.0736E+00  9.7325E-01
             1.0165E+00
 PARAMETER:  1.4374E-01  2.3042E-01  1.9776E-02  5.3433E-02  1.1520E-01 -2.2620E-02  1.7132E-01 -4.9244E+00  1.7101E-01  7.2890E-02
             1.1634E-01
 GRADIENT:  -3.1705E-03  5.3011E-03  1.1090E-03  8.3846E-03  3.8598E-03  1.3779E-04 -8.5424E-04  0.0000E+00 -7.5692E-04 -2.6258E-03
            -2.1198E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      994
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.9109E-04 -1.2915E-02 -3.5982E-04  4.6215E-03 -2.5024E-02
 SE:             2.9800E-02  2.0072E-02  1.5610E-04  2.4547E-02  2.3217E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8953E-01  5.1996E-01  2.1163E-02  8.5066E-01  2.8110E-01

 ETASHRINKSD(%)  1.6534E-01  3.2756E+01  9.9477E+01  1.7764E+01  2.2220E+01
 ETASHRINKVR(%)  3.3041E-01  5.4782E+01  9.9997E+01  3.2372E+01  3.9503E+01
 EBVSHRINKSD(%)  5.5101E-01  3.2013E+01  9.9525E+01  1.8312E+01  2.0699E+01
 EBVSHRINKVR(%)  1.0990E+00  5.3778E+01  9.9998E+01  3.3270E+01  3.7114E+01
 RELATIVEINF(%)  9.8409E+01  2.1322E+00  3.3114E-04  3.7790E+00  7.7903E+00
 EPSSHRINKSD(%)  4.2893E+01
 EPSSHRINKVR(%)  6.7388E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1719.4596929082520     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -984.30886634451383     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.18
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1719.460       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.14E+00  9.23E-01  9.55E-01  1.02E+00  8.85E-01  1.07E+00  1.00E-02  1.07E+00  9.73E-01  1.02E+00
 


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
+        1.29E+03
 
 TH 2
+       -8.32E+00  3.62E+02
 
 TH 3
+        1.68E+01  1.49E+02  3.12E+02
 
 TH 4
+       -1.30E+01  3.21E+02 -1.53E+02  7.08E+02
 
 TH 5
+       -6.27E+00 -2.52E+02 -3.81E+02  1.81E+02  7.12E+02
 
 TH 6
+        1.63E+00 -8.14E-01  5.93E+00 -4.94E-01 -8.32E-03  2.49E+02
 
 TH 7
+        2.18E+00  1.80E+01  7.48E+00 -7.42E+00 -1.31E+01 -5.38E-02  4.25E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.07E+00 -2.17E+01 -2.27E+01  3.10E+01  4.45E+00  5.14E-01  2.11E+01  0.00E+00  8.51E+01
 
 TH10
+        8.35E-01 -3.90E+00 -3.32E+01 -1.38E+01 -5.02E+01  7.63E-01  9.43E+00  0.00E+00  5.93E+00  8.21E+01
 
 TH11
+       -9.40E+00 -1.77E+01 -3.94E+01 -6.92E+00  6.26E+00  3.23E+00  6.38E+00  0.00E+00  7.06E+00  2.22E+01  2.13E+02
 
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
 #CPUT: Total CPU Time in Seconds,       17.114
Stop Time:
Sat Sep 25 10:08:44 CDT 2021
