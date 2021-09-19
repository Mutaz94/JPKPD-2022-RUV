Sat Sep 18 13:53:23 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat17.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1650.64023845834        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.8452E+01 -6.4677E+01  2.2435E+01 -1.1099E+02 -3.0682E+01 -4.3214E-01 -1.6939E+01  8.7706E-01 -8.9050E-02 -3.3798E+01
            -3.1009E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1664.51922204998        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0542E+00  1.1867E+00  1.0496E+00  9.7608E-01  1.2124E+00  9.9596E-01  1.2281E+00  9.2088E-01  8.8938E-01  1.4746E+00
             1.0583E+00
 PARAMETER:  1.5282E-01  2.7121E-01  1.4842E-01  7.5790E-02  2.9263E-01  9.5951E-02  3.0549E-01  1.7571E-02 -1.7231E-02  4.8835E-01
             1.5662E-01
 GRADIENT:   9.2671E+01  1.8858E+01 -5.7389E+00  1.3333E+01  8.1901E+00 -2.9791E-01  8.9289E+00  3.0124E-01 -4.5026E+00  1.0709E+01
            -4.6982E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.68221846216        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0530E+00  1.1912E+00  1.0906E+00  9.7530E-01  1.2409E+00  9.9833E-01  1.2151E+00  8.9180E-01  9.2670E-01  1.4690E+00
             1.0778E+00
 PARAMETER:  1.5160E-01  2.7495E-01  1.8672E-01  7.4988E-02  3.1587E-01  9.8330E-02  2.9485E-01 -1.4512E-02  2.3879E-02  4.8461E-01
             1.7491E-01
 GRADIENT:   3.1262E+01  9.8291E+00  1.3298E+00  6.4228E+00  1.0546E+01 -2.9194E+00  8.4704E+00 -1.9950E+00 -1.6667E+00  5.6495E+00
             8.9848E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1666.03167351049        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  1.0375E+00  1.2684E+00  9.2144E-01  9.1551E-01  1.1578E+00  1.0076E+00  1.0603E+00  7.1042E-01  1.0137E+00  1.3248E+00
             1.0622E+00
 PARAMETER:  1.3685E-01  3.3775E-01  1.8184E-02  1.1725E-02  2.4653E-01  1.0754E-01  1.5860E-01 -2.4190E-01  1.1358E-01  3.8124E-01
             1.6038E-01
 GRADIENT:  -4.8722E+00  2.0899E+00 -2.7876E-01  5.1457E+00  3.0967E-01  1.3976E-01 -3.8376E-01  5.2484E-02  4.7589E-01 -5.0808E-01
            -3.0407E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1666.42954176282        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  1.0418E+00  1.5631E+00  7.1876E-01  7.2709E-01  1.2184E+00  1.0093E+00  9.0605E-01  4.1798E-01  1.1967E+00  1.3360E+00
             1.0699E+00
 PARAMETER:  1.4100E-01  5.4664E-01 -2.3023E-01 -2.1871E-01  2.9754E-01  1.0923E-01  1.3415E-03 -7.7232E-01  2.7958E-01  3.8968E-01
             1.6754E-01
 GRADIENT:   2.1544E-01  1.3657E+01  2.1651E+00  8.0841E+00 -3.8652E+00 -7.2620E-02 -6.6216E-01  1.7572E-04 -3.1269E-01 -7.6747E-01
            -1.0935E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1666.72760944182        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      756
 NPARAMETR:  1.0428E+00  1.7977E+00  5.7143E-01  5.6472E-01  1.3131E+00  1.0104E+00  8.1972E-01  2.0883E-01  1.4174E+00  1.3809E+00
             1.0729E+00
 PARAMETER:  1.4190E-01  6.8648E-01 -4.5961E-01 -4.7143E-01  3.7240E-01  1.1030E-01 -9.8795E-02 -1.4662E+00  4.4881E-01  4.2271E-01
             1.7040E-01
 GRADIENT:   9.5409E-01  3.8476E+00  2.5391E-01  1.5230E+00 -4.4090E-01  3.7357E-02 -3.9396E-01  4.6182E-02 -8.1757E-02 -2.8888E-01
             1.8390E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1666.73492424479        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      934
 NPARAMETR:  1.0426E+00  1.8172E+00  5.5551E-01  5.4861E-01  1.3223E+00  1.0106E+00  8.1517E-01  1.7864E-01  1.4418E+00  1.3859E+00
             1.0716E+00
 PARAMETER:  1.4175E-01  6.9728E-01 -4.8787E-01 -5.0036E-01  3.7935E-01  1.1052E-01 -1.0436E-01 -1.6224E+00  4.6587E-01  4.2634E-01
             1.6912E-01
 GRADIENT:   5.8280E-01 -1.6304E+00 -4.3915E-01 -3.4892E-01  4.3299E-01  9.3492E-02  5.7538E-02  4.0816E-02  1.1039E-01  1.7095E-02
            -6.2846E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1666.75547123838        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  1.0424E+00  1.8045E+00  5.6224E-01  5.5755E-01  1.3151E+00  1.0103E+00  8.1957E-01  3.1770E-02  1.4266E+00  1.3830E+00
             1.0725E+00
 PARAMETER:  1.4151E-01  6.9028E-01 -4.7582E-01 -4.8421E-01  3.7390E-01  1.1022E-01 -9.8970E-02 -3.3492E+00  4.5532E-01  4.2428E-01
             1.6997E-01
 GRADIENT:   1.0420E-01 -1.8318E-01 -1.6886E-02 -9.6884E-02 -6.5197E-02 -9.2644E-03  4.0553E-02  1.1190E-03  2.3270E-03  5.2453E-02
             1.1261E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1666.75599072431        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1274
 NPARAMETR:  1.0423E+00  1.8048E+00  5.6211E-01  5.5743E-01  1.3153E+00  1.0103E+00  8.1936E-01  1.0000E-02  1.4271E+00  1.3829E+00
             1.0722E+00
 PARAMETER:  1.4146E-01  6.9046E-01 -4.7605E-01 -4.8441E-01  3.7405E-01  1.1026E-01 -9.9234E-02 -4.6157E+00  4.5564E-01  4.2415E-01
             1.6974E-01
 GRADIENT:  -8.0916E-03 -1.8574E-02 -4.8682E-03 -5.6864E-03  8.0625E-04  3.8793E-03  3.1173E-03  0.0000E+00  4.2091E-03  1.2177E-03
             1.8298E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1274
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.1093E-05 -3.4916E-02 -3.2413E-04  2.9509E-02 -4.5574E-02
 SE:             2.9826E-02  2.3421E-02  1.0317E-04  2.2167E-02  2.2664E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9944E-01  1.3602E-01  1.6791E-03  1.8311E-01  4.4343E-02

 ETASHRINKSD(%)  7.7641E-02  2.1536E+01  9.9654E+01  2.5738E+01  2.4072E+01
 ETASHRINKVR(%)  1.5522E-01  3.8435E+01  9.9999E+01  4.4852E+01  4.2350E+01
 EBVSHRINKSD(%)  4.8535E-01  2.0315E+01  9.9719E+01  2.8790E+01  2.0615E+01
 EBVSHRINKVR(%)  9.6834E-01  3.6503E+01  9.9999E+01  4.9292E+01  3.6980E+01
 RELATIVEINF(%)  9.8910E+01  4.1886E+00  1.0406E-04  3.2522E+00  1.8349E+01
 EPSSHRINKSD(%)  4.3916E+01
 EPSSHRINKVR(%)  6.8545E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1666.7559907243060     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -931.60516416056782     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1666.756       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.80E+00  5.62E-01  5.57E-01  1.32E+00  1.01E+00  8.19E-01  1.00E-02  1.43E+00  1.38E+00  1.07E+00
 


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
+        9.90E+02
 
 TH 2
+       -5.62E+00  3.46E+02
 
 TH 3
+        1.13E+01  1.25E+02  3.21E+02
 
 TH 4
+       -1.50E+01  3.16E+02 -2.58E+02  9.10E+02
 
 TH 5
+       -3.56E-01 -8.73E+01 -1.59E+02  1.45E+02  2.15E+02
 
 TH 6
+        2.14E+00 -7.68E-01  3.60E+00 -3.73E+00 -1.77E-01  1.93E+02
 
 TH 7
+        1.17E+00  6.63E+00 -9.41E+00 -1.73E+01 -7.30E+00 -4.74E+00  1.25E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.37E-01 -1.79E+01 -3.47E+01  5.87E+01  9.54E-01 -6.39E-01  1.63E+01  0.00E+00  3.67E+01
 
 TH10
+        4.79E-01 -9.60E+00 -2.42E+01  1.96E+00 -3.50E+01 -2.51E-01  7.03E+00  0.00E+00  4.01E+00  4.72E+01
 
 TH11
+       -6.66E+00 -2.12E+01 -3.91E+01  6.03E+00 -1.05E+00  5.68E+00  1.12E+01  0.00E+00  5.62E+00  8.05E+00  1.88E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.524
Stop Time:
Sat Sep 18 13:53:47 CDT 2021
