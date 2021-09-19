Sat Sep 18 09:23:43 CDT 2021
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
$DATA ../../../../data/spa/A1/dat71.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1472.40939049553        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5342E+02 -2.3427E+01  7.5494E+01 -1.1071E+02 -1.8459E+01  1.7146E+01 -2.8717E+01 -2.0521E+01 -2.2145E+01 -2.7177E+01
            -2.1742E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1517.84432744237        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.1981E-01  1.0062E+00  8.9865E-01  1.0692E+00  9.5859E-01  8.8271E-01  1.1130E+00  1.0158E+00  1.0118E+00  1.0584E+00
             1.4454E+00
 PARAMETER:  1.6409E-02  1.0616E-01 -6.8626E-03  1.6690E-01  5.7710E-02 -2.4759E-02  2.0709E-01  1.1563E-01  1.1170E-01  1.5673E-01
             4.6835E-01
 GRADIENT:  -8.0730E+01  1.1504E+01  6.5076E+00  1.1517E+01 -9.2274E+00 -2.8547E+01 -1.0232E+01  4.5225E+00  1.5815E+00  3.1569E+00
             8.0372E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1525.99445824863        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.3782E-01  8.4591E-01  3.6925E-01  1.1076E+00  5.2052E-01  9.0204E-01  1.5133E+00  2.9897E-01  8.8374E-01  6.3593E-01
             1.3438E+00
 PARAMETER:  3.5802E-02 -6.7347E-02 -8.9628E-01  2.0217E-01 -5.5292E-01 -3.0972E-03  5.1432E-01 -1.1074E+00 -2.3592E-02 -3.5267E-01
             3.9549E-01
 GRADIENT:  -2.8702E+01  3.6824E+01 -2.6963E+01  1.0554E+02  2.6915E+01 -2.0192E+01  1.9328E+01  2.0604E+00  1.1602E+01  1.0784E+01
             1.3295E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1532.51532589852        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.7035E-01  8.2790E-01  2.9630E-01  1.0142E+00  4.7083E-01  9.6068E-01  1.3662E+00  7.3752E-02  8.2615E-01  5.4042E-01
             1.2977E+00
 PARAMETER:  6.9905E-02 -8.8859E-02 -1.1164E+00  1.1409E-01 -6.5326E-01  5.9887E-02  4.1203E-01 -2.5070E+00 -9.0980E-02 -5.1541E-01
             3.6062E-01
 GRADIENT:   6.8932E+01 -5.2726E-01 -1.4250E+00 -1.6050E+01  1.6143E+01  5.0318E+00  1.0156E+01  1.1827E-01 -6.6508E+00  8.6362E+00
             4.1451E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1533.39677967460        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.5504E-01  7.2220E-01  2.5790E-01  1.0469E+00  4.0271E-01  9.5966E-01  1.4101E+00  3.2112E-02  8.3655E-01  4.4169E-01
             1.2857E+00
 PARAMETER:  5.4003E-02 -2.2546E-01 -1.2552E+00  1.4580E-01 -8.0953E-01  5.8827E-02  4.4368E-01 -3.3385E+00 -7.8472E-02 -7.1716E-01
             3.5131E-01
 GRADIENT:   3.5183E+01  8.1402E+00 -5.6695E+00  1.1128E+01  3.3375E+00  4.4438E+00  6.5937E+00  1.7949E-02 -3.0260E+00  3.0981E+00
             3.6187E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1534.34149718516        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      457
 NPARAMETR:  9.4953E-01  7.8227E-01  2.7889E-01  1.0303E+00  4.4044E-01  9.5032E-01  1.3367E+00  3.8470E-02  8.5115E-01  4.4474E-01
             1.2883E+00
 PARAMETER:  4.8208E-02 -1.4555E-01 -1.1769E+00  1.2989E-01 -7.1997E-01  4.9049E-02  3.9020E-01 -3.1579E+00 -6.1163E-02 -7.1027E-01
             3.5336E-01
 GRADIENT:   3.3298E-01 -2.8142E-01  1.0463E+00 -1.0913E+00 -1.0559E+00 -2.6399E-01 -4.0461E-01  1.9707E-02  1.6964E-02  1.4304E-01
             1.5776E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1534.34573353863        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      634
 NPARAMETR:  9.4909E-01  7.8502E-01  2.7708E-01  1.0285E+00  4.4017E-01  9.5102E-01  1.3355E+00  3.5616E-02  8.5170E-01  4.3786E-01
             1.2886E+00
 PARAMETER:  4.7749E-02 -1.4204E-01 -1.1834E+00  1.2807E-01 -7.2060E-01  4.9784E-02  3.8934E-01 -3.2350E+00 -6.0520E-02 -7.2585E-01
             3.5356E-01
 GRADIENT:  -3.7176E-01  1.7906E-01  7.4135E-02  6.1024E-01 -3.0930E-01  3.4146E-02 -1.1732E-01  1.6516E-02 -5.6199E-02 -1.3615E-01
             8.5799E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1534.35402964310        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  9.4924E-01  7.8642E-01  2.7739E-01  1.0278E+00  4.4095E-01  9.5092E-01  1.3347E+00  1.0000E-02  8.5203E-01  4.4006E-01
             1.2882E+00
 PARAMETER:  4.7910E-02 -1.4026E-01 -1.1823E+00  1.2740E-01 -7.1883E-01  4.9675E-02  3.8872E-01 -4.6067E+00 -6.0138E-02 -7.2084E-01
             3.5328E-01
 GRADIENT:   1.1645E-02  1.1922E-02  1.7737E-02 -2.1068E-02 -2.6857E-02  6.8345E-04  4.3170E-03  0.0000E+00  2.4242E-03 -3.6648E-03
            -1.8630E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1534.35402964310        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      834
 NPARAMETR:  9.4924E-01  7.8642E-01  2.7739E-01  1.0278E+00  4.4095E-01  9.5092E-01  1.3347E+00  1.0000E-02  8.5203E-01  4.4006E-01
             1.2882E+00
 PARAMETER:  4.7910E-02 -1.4026E-01 -1.1823E+00  1.2740E-01 -7.1883E-01  4.9675E-02  3.8872E-01 -4.6067E+00 -6.0138E-02 -7.2084E-01
             3.5328E-01
 GRADIENT:   1.1645E-02  1.1922E-02  1.7737E-02 -2.1068E-02 -2.6857E-02  6.8345E-04  4.3170E-03  0.0000E+00  2.4242E-03 -3.6648E-03
            -1.8630E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      834
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.6890E-04  8.2154E-03 -4.2380E-04 -6.6684E-03  8.7313E-03
 SE:             2.9751E-02  2.5153E-02  1.9238E-04  2.7247E-02  1.4745E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8743E-01  7.4396E-01  2.7596E-02  8.0666E-01  5.5375E-01

 ETASHRINKSD(%)  3.2915E-01  1.5735E+01  9.9356E+01  8.7206E+00  5.0602E+01
 ETASHRINKVR(%)  6.5722E-01  2.8994E+01  9.9996E+01  1.6681E+01  7.5598E+01
 EBVSHRINKSD(%)  6.9575E-01  1.4130E+01  9.9363E+01  8.9656E+00  5.2292E+01
 EBVSHRINKVR(%)  1.3867E+00  2.6263E+01  9.9996E+01  1.7127E+01  7.7240E+01
 RELATIVEINF(%)  9.7554E+01  7.6071E+00  5.8684E-04  1.8307E+01  1.4215E+00
 EPSSHRINKSD(%)  4.3620E+01
 EPSSHRINKVR(%)  6.8213E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1534.3540296431004     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -799.20320307936220     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1534.354       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  7.86E-01  2.77E-01  1.03E+00  4.41E-01  9.51E-01  1.33E+00  1.00E-02  8.52E-01  4.40E-01  1.29E+00
 


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
+        1.35E+03
 
 TH 2
+       -8.65E+00  7.55E+02
 
 TH 3
+       -1.18E+02  1.15E+03  6.30E+03
 
 TH 4
+       -1.66E+01  1.44E+02 -1.55E+03  1.21E+03
 
 TH 5
+        3.85E+01 -1.76E+03 -5.71E+03  1.09E+03  6.73E+03
 
 TH 6
+        1.89E+00 -5.07E-01  2.55E+00 -4.57E+00  6.04E+00  2.14E+02
 
 TH 7
+        8.24E-01  3.40E+01 -1.32E+02 -4.36E+00  4.65E+01 -1.30E-01  5.47E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.67E+00 -1.19E+01  3.47E+01  5.75E+00 -1.15E+01 -7.87E-01  1.33E+01  0.00E+00  1.93E+02
 
 TH10
+       -4.90E-01 -1.73E+01 -2.55E+02 -2.39E+01  9.16E+01 -9.42E-01  2.37E+01  0.00E+00  4.38E+00  1.05E+02
 
 TH11
+       -9.99E+00 -1.31E+01 -7.02E+01 -5.04E+00  1.19E+01  2.47E+00  8.13E+00  0.00E+00  7.84E+00  2.33E+01  1.30E+02
 
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
 #CPUT: Total CPU Time in Seconds,       13.562
Stop Time:
Sat Sep 18 09:23:58 CDT 2021
