Sat Sep 18 13:34:57 CDT 2021
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
$DATA ../../../../data/spa/S2/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1675.47568769684        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.4210E+01 -1.0040E+02 -2.5527E+00 -1.6200E+02 -1.6352E+01 -9.3240E+00 -1.3978E+01  1.0315E+01 -2.3480E+01  1.1098E+01
             5.4056E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1682.70851459845        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9620E-01  1.0489E+00  1.0204E+00  9.7057E-01  1.0874E+00  1.0136E+00  1.1184E+00  9.3680E-01  1.1046E+00  9.8393E-01
             8.6923E-01
 PARAMETER:  9.6195E-02  1.4770E-01  1.2021E-01  7.0129E-02  1.8375E-01  1.1353E-01  2.1187E-01  3.4710E-02  1.9944E-01  8.3798E-02
            -4.0150E-02
 GRADIENT:   6.9848E+01 -9.9494E+01 -1.4701E+01 -1.1403E+02  4.1358E+01 -2.0461E+00 -5.8777E+00  5.7287E+00  3.3325E+00 -3.7364E+00
             2.4622E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1690.18374655830        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.9263E-01  1.4870E+00  7.7499E-01  7.8964E-01  1.1319E+00  1.0051E+00  1.0285E+00  3.6517E-01  1.2445E+00  1.0541E+00
             9.0986E-01
 PARAMETER:  9.2601E-02  4.9675E-01 -1.5491E-01 -1.3617E-01  2.2393E-01  1.0512E-01  1.2813E-01 -9.0740E-01  3.1875E-01  1.5268E-01
             5.5394E-03
 GRADIENT:   5.1381E+01  5.8645E+01  1.9460E+01 -1.1924E+01 -2.1298E+01 -7.7208E+00  1.2650E+01  1.2549E-01  1.3703E+00  3.4613E+00
             1.7024E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1692.09434221710        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.7249E-01  1.1705E+00  7.7808E-01  9.7219E-01  9.6180E-01  1.0239E+00  1.0843E+00  2.9970E-01  1.0441E+00  8.9018E-01
             8.6288E-01
 PARAMETER:  7.2104E-02  2.5746E-01 -1.5092E-01  7.1791E-02  6.1053E-02  1.2363E-01  1.8090E-01 -1.1050E+00  1.4311E-01 -1.6331E-02
            -4.7477E-02
 GRADIENT:   1.3371E+01 -1.0380E+00  2.7169E+00 -1.1591E+01  2.9820E+00  1.5921E+00 -6.3620E+00  3.3468E-01 -2.9304E+00 -4.8485E+00
            -1.7682E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1692.24287381597        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.6556E-01  1.0672E+00  7.9142E-01  1.0416E+00  9.1718E-01  1.0217E+00  1.2073E+00  2.0988E-01  9.9651E-01  9.0240E-01
             8.6427E-01
 PARAMETER:  6.4957E-02  1.6505E-01 -1.3393E-01  1.4074E-01  1.3549E-02  1.2146E-01  2.8842E-01 -1.4612E+00  9.6506E-02 -2.7000E-03
            -4.5867E-02
 GRADIENT:  -1.7304E+00 -1.5511E+00 -1.9812E+00  1.3959E+00  2.1182E+00  1.4782E-01 -4.9725E-01  2.3195E-01 -1.5052E-01  4.0029E-01
             4.8023E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1693.07260442072        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      421
 NPARAMETR:  9.8573E-01  1.0599E+00  7.8373E-01  1.0508E+00  9.0781E-01  1.0390E+00  1.2343E+00  1.0291E-01  9.8623E-01  8.8954E-01
             8.6734E-01
 PARAMETER:  8.5628E-02  1.5821E-01 -1.4369E-01  1.4951E-01  3.2750E-03  1.3827E-01  3.1054E-01 -2.1739E+00  8.6136E-02 -1.7051E-02
            -4.2323E-02
 GRADIENT:  -1.0041E+01 -3.6023E+00 -1.0976E+00 -3.5740E+00  1.9740E+00 -6.3918E-01 -7.5227E-01  3.4977E-02 -1.6383E+00 -1.1542E+00
             1.3345E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1693.15441298170        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      597
 NPARAMETR:  9.9069E-01  1.1245E+00  7.7715E-01  1.0147E+00  9.3521E-01  1.0412E+00  1.1839E+00  5.3088E-02  1.0235E+00  9.1109E-01
             8.6527E-01
 PARAMETER:  9.0641E-02  2.1737E-01 -1.5212E-01  1.1454E-01  3.3012E-02  1.4039E-01  2.6879E-01 -2.8358E+00  1.2321E-01  6.8818E-03
            -4.4715E-02
 GRADIENT:  -2.2448E-01 -9.4765E-04  2.2525E-01 -6.5758E-01 -3.7523E-01  6.6042E-02  1.9567E-01  9.8152E-03  9.7510E-03  2.0534E-02
             1.8679E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1693.15925524566        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      773
 NPARAMETR:  9.9086E-01  1.1226E+00  7.7657E-01  1.0163E+00  9.3393E-01  1.0413E+00  1.1837E+00  1.4446E-02  1.0222E+00  9.1060E-01
             8.6489E-01
 PARAMETER:  9.0816E-02  2.1568E-01 -1.5287E-01  1.1612E-01  3.1642E-02  1.4043E-01  2.6868E-01 -4.1373E+00  1.2200E-01  6.3484E-03
            -4.5151E-02
 GRADIENT:   1.2666E-01  1.6119E-01  3.6284E-02  1.8184E-01 -2.9473E-02  7.8280E-02 -3.6496E-02  7.0911E-04 -5.1481E-02 -3.4138E-02
            -1.3776E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1693.15948058116        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  9.9078E-01  1.1218E+00  7.7672E-01  1.0166E+00  9.3368E-01  1.0410E+00  1.1847E+00  1.0000E-02  1.0221E+00  9.1070E-01
             8.6489E-01
 PARAMETER:  9.0734E-02  2.1493E-01 -1.5268E-01  1.1644E-01  3.1383E-02  1.4018E-01  2.6946E-01 -4.7923E+00  1.2190E-01  6.4590E-03
            -4.5151E-02
 GRADIENT:  -3.0842E-02 -3.3919E-02 -1.4425E-02 -3.2047E-02  1.3602E-02 -1.5825E-02  6.8912E-03  0.0000E+00  9.5630E-03  8.2328E-03
             2.5256E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      865
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7940E-05 -8.7346E-03 -5.2834E-04  3.7987E-03 -1.9258E-02
 SE:             2.9879E-02  2.1850E-02  2.0163E-04  2.4916E-02  2.2783E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9925E-01  6.8934E-01  8.7836E-03  8.7882E-01  3.9795E-01

 ETASHRINKSD(%)  1.0000E-10  2.6799E+01  9.9325E+01  1.6530E+01  2.3676E+01
 ETASHRINKVR(%)  1.0000E-10  4.6416E+01  9.9995E+01  3.0327E+01  4.1746E+01
 EBVSHRINKSD(%)  3.0259E-01  2.6148E+01  9.9422E+01  1.6893E+01  2.2659E+01
 EBVSHRINKVR(%)  6.0427E-01  4.5458E+01  9.9997E+01  3.0932E+01  4.0184E+01
 RELATIVEINF(%)  9.9211E+01  3.8080E+00  5.5766E-04  6.0285E+00  6.7222E+00
 EPSSHRINKSD(%)  4.4806E+01
 EPSSHRINKVR(%)  6.9536E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1693.1594805811617     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -958.00865401742351     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.29
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
 





 #OBJV:********************************************    -1693.159       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  1.12E+00  7.77E-01  1.02E+00  9.34E-01  1.04E+00  1.18E+00  1.00E-02  1.02E+00  9.11E-01  8.65E-01
 


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
+        1.04E+03
 
 TH 2
+       -3.99E+00  3.61E+02
 
 TH 3
+        1.28E+01  2.30E+02  6.94E+02
 
 TH 4
+       -4.25E+00  2.57E+02 -2.67E+02  7.13E+02
 
 TH 5
+       -3.21E+00 -3.08E+02 -6.34E+02  2.75E+02  9.09E+02
 
 TH 6
+        1.26E+00 -1.17E+00  2.81E+00 -1.47E+00 -1.14E+00  1.81E+02
 
 TH 7
+        7.67E-01  2.50E+01 -1.55E+01 -7.60E+00 -3.57E+00 -6.06E-01  4.77E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.31E-01 -2.32E+01 -3.24E+01  3.55E+01 -5.04E-01 -1.60E+00  1.46E+01  0.00E+00  9.98E+01
 
 TH10
+        4.53E-01 -1.04E+01 -6.53E+01 -1.67E+01 -5.29E+01 -1.32E+00  1.37E+01  0.00E+00  1.05E+01  8.69E+01
 
 TH11
+       -8.24E+00 -1.58E+01 -5.09E+01 -5.04E+00  6.67E+00  1.26E+00  6.94E+00  0.00E+00  5.60E+00  2.24E+01  2.79E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.317
Stop Time:
Sat Sep 18 13:35:14 CDT 2021
