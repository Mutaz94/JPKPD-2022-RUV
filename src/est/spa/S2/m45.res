Wed Sep 29 17:25:01 CDT 2021
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
$DATA ../../../../data/spa/S2/dat45.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1617.91565859230        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1547E+02 -2.3794E+00 -3.6314E+01  6.4919E+01  7.3617E+01  4.9937E+00  5.9171E+00  7.0765E+00  1.2161E+01  5.2886E+00
            -3.5995E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1627.91351222146        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8653E-01  1.0246E+00  1.0313E+00  1.0109E+00  9.7928E-01  1.1639E+00  9.8644E-01  9.8036E-01  9.8835E-01  9.6771E-01
             1.0862E+00
 PARAMETER:  8.6441E-02  1.2426E-01  1.3086E-01  1.1088E-01  7.9067E-02  2.5176E-01  8.6350E-02  8.0164E-02  8.8277E-02  6.7182E-02
             1.8269E-01
 GRADIENT:   2.5199E+01  3.4096E-01 -9.4305E+00  1.1696E+01  1.4412E+01 -1.6955E-01  3.0459E+00  4.3719E+00  1.7126E+00  5.1710E+00
            -7.0012E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1628.95052166847        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8411E-01  9.8699E-01  8.6905E-01  1.0226E+00  8.8794E-01  1.1555E+00  9.6762E-01  5.9989E-01  9.6114E-01  9.0088E-01
             1.1045E+00
 PARAMETER:  8.3985E-02  8.6900E-02 -4.0349E-02  1.2239E-01 -1.8852E-02  2.4450E-01  6.7081E-02 -4.1100E-01  6.0364E-02 -4.3794E-03
             1.9942E-01
 GRADIENT:   1.9226E+01 -5.8504E+00 -1.5524E+01  1.2857E+01  1.7818E+01 -3.5038E+00 -2.9587E+00  2.0263E+00 -1.6222E+00  3.3936E+00
             6.3179E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1629.85826158225        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.6879E-01  1.1249E+00  7.6182E-01  9.2591E-01  8.8799E-01  1.1688E+00  9.4854E-01  3.1345E-01  1.0247E+00  8.7212E-01
             1.0900E+00
 PARAMETER:  6.8289E-02  2.1767E-01 -1.7205E-01  2.3022E-02 -1.8792E-02  2.5598E-01  4.7169E-02 -1.0601E+00  1.2437E-01 -3.6832E-02
             1.8618E-01
 GRADIENT:  -7.9517E+00  4.1815E-01 -2.9795E+00  3.3098E+00  5.2282E-02  9.1826E-01  2.4720E-01  4.5185E-01  6.9589E-01  1.4895E+00
             3.3452E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1630.05253545809        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      710
 NPARAMETR:  9.7387E-01  1.2317E+00  7.3096E-01  8.5603E-01  9.2662E-01  1.1662E+00  8.8605E-01  1.0188E-01  1.0888E+00  8.8918E-01
             1.0917E+00
 PARAMETER:  7.3518E-02  3.0837E-01 -2.1340E-01 -5.5446E-02  2.3786E-02  2.5377E-01 -2.0983E-02 -2.1840E+00  1.8509E-01 -1.7456E-02
             1.8771E-01
 GRADIENT:   1.7645E-01 -4.8248E-01 -1.1193E-01 -3.5536E-01  4.5752E-01 -1.0439E-01  7.2954E-02  3.9027E-02  1.7407E-01 -1.0810E-01
            -1.1499E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1630.05314484384        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.7375E-01  1.2448E+00  7.2590E-01  8.4752E-01  9.3080E-01  1.1665E+00  8.7870E-01  7.9276E-02  1.0968E+00  8.9129E-01
             1.0920E+00
 PARAMETER:  7.3403E-02  3.1899E-01 -2.2034E-01 -6.5447E-02  2.8287E-02  2.5398E-01 -2.9310E-02 -2.4348E+00  1.9238E-01 -1.5087E-02
             1.8799E-01
 GRADIENT:  -6.4817E-02 -6.0607E-01 -2.3383E-01 -4.4936E-01  5.1955E-01 -2.9288E-02  2.9624E-02  2.4289E-02  1.1362E-01 -1.5297E-02
            -1.0630E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1630.05647570581        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1063
 NPARAMETR:  9.7371E-01  1.2542E+00  7.2198E-01  8.4136E-01  9.3370E-01  1.1667E+00  8.7373E-01  3.7061E-02  1.1022E+00  8.9285E-01
             1.0922E+00
 PARAMETER:  7.3359E-02  3.2650E-01 -2.2576E-01 -7.2739E-02  3.1395E-02  2.5417E-01 -3.4979E-02 -3.1952E+00  1.9729E-01 -1.3334E-02
             1.8824E-01
 GRADIENT:  -1.7696E-01 -7.7060E-01 -3.2186E-01 -6.1053E-01  6.0527E-01  3.7705E-02 -2.0223E-03  5.4673E-03  1.0163E-02  5.6738E-02
             8.3059E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1630.08480566418        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1248             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7516E-01  1.1980E+00  7.3795E-01  8.7683E-01  9.1378E-01  1.1756E+00  9.0356E-01  1.0000E-02  1.0672E+00  8.8475E-01
             1.0917E+00
 PARAMETER:  7.4851E-02  2.8066E-01 -2.0387E-01 -3.1446E-02  9.8356E-03  2.6178E-01 -1.4124E-03 -7.4867E+00  1.6500E-01 -2.2454E-02
             1.8775E-01
 GRADIENT:   3.7793E+02  1.3395E+02  2.6655E+00  5.0539E+01  8.5778E+00  1.8863E+02  3.1906E+00  0.0000E+00  8.7557E+00  5.2546E-01
             1.5970E+00

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1630.08769410753        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  9.7503E-01  1.1990E+00  7.3874E-01  8.7745E-01  9.1313E-01  1.1714E+00  9.0523E-01  1.0000E-02  1.0677E+00  8.8462E-01
             1.0914E+00
 PARAMETER:  7.3991E-02  2.8087E-01 -2.0352E-01 -3.1320E-02  9.5998E-03  2.5869E-01 -5.6145E-04 -7.4867E+00  1.6499E-01 -2.2450E-02
             1.8768E-01
 GRADIENT:  -5.3201E-01 -3.5680E-01 -1.7171E-01 -3.5878E-01  3.7815E-01  8.8308E-02 -4.3532E-02  0.0000E+00 -4.4593E-02  9.8242E-03
             4.6558E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1362
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2549E-04 -1.6567E-02 -3.3430E-04  7.3198E-03 -2.3117E-02
 SE:             2.9845E-02  2.0251E-02  1.5124E-04  2.4697E-02  2.3072E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9130E-01  4.1330E-01  2.7079E-02  7.6694E-01  3.1639E-01

 ETASHRINKSD(%)  1.4169E-02  3.2157E+01  9.9493E+01  1.7262E+01  2.2705E+01
 ETASHRINKVR(%)  2.8336E-02  5.3973E+01  9.9997E+01  3.1544E+01  4.0254E+01
 EBVSHRINKSD(%)  3.7890E-01  3.1664E+01  9.9531E+01  1.7615E+01  2.1911E+01
 EBVSHRINKVR(%)  7.5636E-01  5.3302E+01  9.9998E+01  3.2127E+01  3.9021E+01
 RELATIVEINF(%)  9.9089E+01  1.9201E+00  2.3358E-04  3.5963E+00  6.3157E+00
 EPSSHRINKSD(%)  4.3111E+01
 EPSSHRINKVR(%)  6.7636E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1630.0876941075287     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -894.93686754379053     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.42
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1630.088       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  1.20E+00  7.38E-01  8.77E-01  9.14E-01  1.17E+00  9.04E-01  1.00E-02  1.07E+00  8.85E-01  1.09E+00
 


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
+        8.46E+02
 
 TH 2
+       -5.33E+00  4.63E+02
 
 TH 3
+        9.60E+00  2.27E+02  4.80E+02
 
 TH 4
+       -1.02E+01  3.72E+02 -2.21E+02  8.65E+02
 
 TH 5
+       -1.85E+00 -3.76E+02 -5.88E+02  2.45E+02  1.03E+03
 
 TH 6
+        3.74E-01 -1.23E+00  2.17E+00 -2.86E+00 -9.47E-01  1.43E+02
 
 TH 7
+        6.46E-01  1.83E+01  2.36E+00 -8.90E+00 -1.49E+01  2.31E-02  5.79E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.69E+00 -2.58E+01 -2.62E+01  3.65E+01 -2.90E-01 -3.29E-01  2.13E+01  0.00E+00  8.77E+01
 
 TH10
+       -7.80E-01 -8.60E+00 -4.95E+01 -1.48E+01 -5.70E+01  2.27E-01  1.95E+01  0.00E+00  6.35E+00  9.41E+01
 
 TH11
+       -5.80E+00 -1.82E+01 -3.28E+01 -2.93E+00  5.25E+00  1.22E+00  8.17E+00  0.00E+00  8.31E+00  2.49E+01  1.80E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.354
Stop Time:
Wed Sep 29 17:25:26 CDT 2021
