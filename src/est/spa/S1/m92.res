Sat Sep 25 10:11:18 CDT 2021
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
$DATA ../../../../data/spa/S1/dat92.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1669.59745602567        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.0856E+01 -8.2102E+01 -4.1483E+01 -6.7833E+01  1.6194E+01  2.0507E+01 -2.6228E+01  1.1392E+01 -1.2054E+01  1.6399E+01
             1.8856E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1678.68156997161        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      110
 NPARAMETR:  9.9281E-01  1.1798E+00  1.2391E+00  9.4018E-01  1.1466E+00  8.9808E-01  1.2650E+00  9.1800E-01  1.0762E+00  8.9128E-01
             9.2777E-01
 PARAMETER:  9.2784E-02  2.6533E-01  3.1436E-01  3.8317E-02  2.3684E-01 -7.4919E-03  3.3505E-01  1.4442E-02  1.7342E-01 -1.5094E-02
             2.5029E-02
 GRADIENT:   8.1688E+00  3.8359E+00  2.8877E+01 -1.7023E+01 -2.0064E+00 -2.6653E+01  2.2720E+00 -7.8098E+00  4.0869E+00 -1.9025E+01
            -2.3157E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1679.82370666739        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      287
 NPARAMETR:  9.8859E-01  1.3239E+00  1.1158E+00  8.7234E-01  1.1703E+00  9.1684E-01  1.1493E+00  7.7085E-01  1.1916E+00  1.0048E+00
             9.1912E-01
 PARAMETER:  8.8525E-02  3.8061E-01  2.0957E-01 -3.6576E-02  2.5723E-01  1.3176E-02  2.3911E-01 -1.6026E-01  2.7527E-01  1.0475E-01
             1.5667E-02
 GRADIENT:  -5.7738E+00  2.4928E+01  2.2326E+01  1.1898E+01 -1.1015E+01 -1.8000E+01  1.8907E+00 -5.1557E+00  9.3414E+00 -4.5810E+00
            -2.4377E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1682.88848250657        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      462
 NPARAMETR:  9.9199E-01  1.2633E+00  8.5105E-01  8.8121E-01  1.0286E+00  9.6125E-01  1.2175E+00  5.3396E-01  1.0352E+00  8.6600E-01
             9.6582E-01
 PARAMETER:  9.1953E-02  3.3372E-01 -6.1282E-02 -2.6464E-02  1.2820E-01  6.0478E-02  2.9678E-01 -5.2744E-01  1.3457E-01 -4.3875E-02
             6.5222E-02
 GRADIENT:  -1.3323E+00  1.7864E+00 -8.0037E-01  1.3190E+00 -4.6614E+00  8.7802E-01  5.8518E-01  6.6431E-01  1.5654E-01  1.4120E+00
             1.8101E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1683.06854150751        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      637
 NPARAMETR:  9.9338E-01  1.4700E+00  7.2795E-01  7.4698E-01  1.0820E+00  9.6061E-01  1.0781E+00  2.8065E-01  1.1619E+00  8.8351E-01
             9.6309E-01
 PARAMETER:  9.3353E-02  4.8529E-01 -2.1752E-01 -1.9172E-01  1.7885E-01  5.9812E-02  1.7518E-01 -1.1707E+00  2.5003E-01 -2.3857E-02
             6.2393E-02
 GRADIENT:  -1.6137E-01  2.9315E+00  1.3393E-01  1.6893E+00 -1.2496E+00  2.8375E-01 -5.2970E-02  6.5105E-02 -1.8506E-01 -6.9192E-02
             2.6128E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1683.08234210989        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      812
 NPARAMETR:  9.9355E-01  1.5317E+00  7.0006E-01  7.0459E-01  1.1077E+00  9.6013E-01  1.0427E+00  2.1366E-01  1.2129E+00  9.0089E-01
             9.6235E-01
 PARAMETER:  9.3525E-02  5.2640E-01 -2.5659E-01 -2.5013E-01  2.0232E-01  5.9315E-02  1.4185E-01 -1.4434E+00  2.9301E-01 -4.3687E-03
             6.1627E-02
 GRADIENT:  -3.3549E-02 -4.5683E-01 -1.8576E-01 -2.9110E-01  3.7505E-01  5.3345E-02  6.2052E-02  4.0995E-02 -5.4424E-03  2.8254E-02
             8.2462E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1683.09795885244        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      992
 NPARAMETR:  9.9359E-01  1.5247E+00  6.9058E-01  7.0860E-01  1.0972E+00  9.6022E-01  1.0481E+00  7.4895E-02  1.2040E+00  8.9085E-01
             9.6249E-01
 PARAMETER:  9.3572E-02  5.2177E-01 -2.7022E-01 -2.4447E-01  1.9280E-01  5.9402E-02  1.4698E-01 -2.4917E+00  2.8565E-01 -1.5580E-02
             6.1768E-02
 GRADIENT:  -7.4035E-02  2.8558E-02 -1.0758E-01  5.3076E-02  5.6019E-02  5.6375E-02  7.8631E-02  3.9163E-03 -2.2750E-02  1.2808E-02
             7.6594E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1683.09971556726        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  9.9363E-01  1.5292E+00  6.8816E-01  7.0558E-01  1.0988E+00  9.6010E-01  1.0453E+00  1.0134E-02  1.2080E+00  8.9182E-01
             9.6236E-01
 PARAMETER:  9.3610E-02  5.2476E-01 -2.7374E-01 -2.4874E-01  1.9419E-01  5.9285E-02  1.4430E-01 -4.4918E+00  2.8897E-01 -1.4487E-02
             6.1631E-02
 GRADIENT:  -6.1088E-03  1.6068E-02  7.7316E-04  4.7425E-03 -1.2397E-02  6.7087E-03  2.2663E-03  6.5698E-05  4.5434E-05  1.2238E-03
             6.7785E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1683.09971653438        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1224
 NPARAMETR:  9.9363E-01  1.5292E+00  6.8815E-01  7.0556E-01  1.0988E+00  9.6009E-01  1.0453E+00  1.0000E-02  1.2080E+00  8.9183E-01
             9.6235E-01
 PARAMETER:  9.3613E-02  5.2476E-01 -2.7375E-01 -2.4876E-01  1.9420E-01  5.9272E-02  1.4427E-01 -4.8067E+00  2.8900E-01 -1.4475E-02
             6.1619E-02
 GRADIENT:   7.8709E-04 -3.2561E-03 -4.2667E-03  1.7523E-03 -1.6560E-03  1.4293E-03 -1.9244E-03  0.0000E+00  1.4657E-03  1.3601E-03
             2.0911E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1224
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.0474E-05 -1.6591E-02 -3.5209E-04  1.4636E-02 -2.9105E-02
 SE:             2.9850E-02  2.4679E-02  1.3424E-04  2.2473E-02  2.1566E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9758E-01  5.0142E-01  8.7186E-03  5.1487E-01  1.7715E-01

 ETASHRINKSD(%)  1.0000E-10  1.7323E+01  9.9550E+01  2.4712E+01  2.7752E+01
 ETASHRINKVR(%)  1.0000E-10  3.1645E+01  9.9998E+01  4.3318E+01  4.7802E+01
 EBVSHRINKSD(%)  4.2526E-01  1.6945E+01  9.9623E+01  2.6184E+01  2.6438E+01
 EBVSHRINKVR(%)  8.4871E-01  3.1018E+01  9.9999E+01  4.5511E+01  4.5887E+01
 RELATIVEINF(%)  9.9017E+01  5.5181E+00  1.8583E-04  4.0643E+00  9.4778E+00
 EPSSHRINKSD(%)  4.3983E+01
 EPSSHRINKVR(%)  6.8621E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1683.0997165343822     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -947.94888997064402     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.30
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1683.100       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  1.53E+00  6.88E-01  7.06E-01  1.10E+00  9.60E-01  1.05E+00  1.00E-02  1.21E+00  8.92E-01  9.62E-01
 


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
+        1.21E+03
 
 TH 2
+       -3.87E+00  3.12E+02
 
 TH 3
+        1.34E+01  1.31E+02  3.85E+02
 
 TH 4
+       -1.21E+01  2.51E+02 -2.95E+02  8.57E+02
 
 TH 5
+       -5.20E+00 -1.76E+02 -3.64E+02  3.18E+02  6.30E+02
 
 TH 6
+        6.37E-01 -6.80E-01  2.70E+00 -3.36E+00 -9.34E-01  2.12E+02
 
 TH 7
+        8.09E-01  1.92E+01 -1.81E+01 -1.79E+01  2.97E+00  9.64E-01  9.29E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.17E+00 -1.69E+01 -3.54E+01  5.31E+01 -1.24E+00 -8.07E-01  1.25E+01  0.00E+00  5.36E+01
 
 TH10
+        1.37E+00 -1.37E+01 -3.57E+01 -8.53E+00 -7.15E+01 -2.45E+00  9.69E+00  0.00E+00  8.84E+00  8.45E+01
 
 TH11
+       -6.27E+00 -1.37E+01 -2.83E+01  8.92E-01 -4.14E+00  1.17E+00  7.35E+00  0.00E+00  8.15E+00  1.98E+01  2.36E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.331
Stop Time:
Sat Sep 25 10:11:46 CDT 2021
