Fri Sep 24 20:21:43 CDT 2021
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
$DATA ../../../../data/int/A1/dat16.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3478.40663439345        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.7066E+01 -4.6582E+00 -8.7583E+00 -3.4417E+01  1.5595E+02  4.7913E+01 -6.2917E+01 -1.4060E+01  4.6171E+01 -8.9245E+01
            -7.8896E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3601.24632235665        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0554E+00  9.0308E-01  9.0757E-01  1.1029E+00  8.0451E-01  7.6973E-01  1.3275E+00  8.8882E-01  6.9270E-01  1.3275E+00
             1.2911E+00
 PARAMETER:  1.5390E-01 -1.9430E-03  3.0125E-03  1.9792E-01 -1.1753E-01 -1.6172E-01  3.8329E-01 -1.7865E-02 -2.6715E-01  3.8329E-01
             3.5546E-01
 GRADIENT:   1.3071E+02  1.5153E+01 -6.1669E+00  1.0133E+02 -4.1190E+01 -5.5127E+01 -1.0033E+00  1.1356E+01 -2.5278E+01  3.9338E+01
            -6.2166E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3608.46821234801        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0512E+00  8.9757E-01  8.5763E-01  1.0952E+00  7.8623E-01  7.8665E-01  1.3875E+00  3.9216E-01  7.1461E-01  1.2821E+00
             1.3254E+00
 PARAMETER:  1.4994E-01 -8.0672E-03 -5.3581E-02  1.9095E-01 -1.4051E-01 -1.3997E-01  4.2749E-01 -8.3608E-01 -2.3602E-01  3.4847E-01
             3.8169E-01
 GRADIENT:   1.0988E+02  2.1685E+01 -7.9304E+00  5.8208E+01 -3.5251E+01 -4.3859E+01  1.0233E+01  8.7778E-01 -1.6724E+01  3.8985E+01
            -3.4175E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3617.04148258216        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0200E+00  8.6575E-01  8.3607E-01  1.0880E+00  7.8782E-01  8.6995E-01  1.2783E+00  3.0588E-01  7.8958E-01  1.0652E+00
             1.3491E+00
 PARAMETER:  1.1980E-01 -4.4158E-02 -7.9048E-02  1.8437E-01 -1.3849E-01 -3.9314E-02  3.4555E-01 -1.0846E+00 -1.3626E-01  1.6321E-01
             3.9943E-01
 GRADIENT:   7.1644E-01  1.0342E+00 -7.0898E+00  5.2493E-01  2.7860E+00  3.5409E-01  5.7425E-01  7.5780E-01 -9.2063E-01  1.4487E+00
             2.0258E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3617.11423386007        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0196E+00  8.5132E-01  8.2454E-01  1.0939E+00  7.7350E-01  8.6961E-01  1.2662E+00  1.7504E-01  7.9780E-01  1.0518E+00
             1.3508E+00
 PARAMETER:  1.1944E-01 -6.0964E-02 -9.2932E-02  1.8974E-01 -1.5683E-01 -3.9712E-02  3.3601E-01 -1.6427E+00 -1.2590E-01  1.5053E-01
             4.0067E-01
 GRADIENT:  -4.2748E-01 -2.6202E-01 -4.6505E-01 -1.6647E-01 -4.4744E-01  2.1243E-01  1.0287E-02  1.8199E-01 -2.3096E-01  1.6685E-01
             2.1420E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3617.14591390665        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0198E+00  8.4750E-01  8.1932E-01  1.0955E+00  7.6934E-01  8.6909E-01  1.2634E+00  4.7527E-02  8.0032E-01  1.0493E+00
             1.3521E+00
 PARAMETER:  1.1960E-01 -6.5465E-02 -9.9278E-02  1.9125E-01 -1.6222E-01 -4.0303E-02  3.3379E-01 -2.9465E+00 -1.2274E-01  1.4811E-01
             4.0169E-01
 GRADIENT:  -1.8743E-03 -2.2541E-02 -2.0336E-01 -9.1118E-03  1.9251E-02 -1.1874E-03  2.9918E-03  1.1992E-02  3.4223E-04  3.3733E-02
             1.9616E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3617.40831796737        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  1.0299E+00  8.6821E-01  8.3887E-01  1.0914E+00  7.9073E-01  8.7251E-01  1.2856E+00  2.8880E-02  7.9089E-01  1.0672E+00
             1.3535E+00
 PARAMETER:  1.2950E-01 -4.1327E-02 -7.5702E-02  1.8745E-01 -1.3480E-01 -3.6376E-02  3.5124E-01 -3.4446E+00 -1.3459E-01  1.6506E-01
             4.0272E-01
 GRADIENT:   2.6150E-01  1.3877E-04  1.6954E-01  2.1211E-01 -9.2533E-02  1.0569E-01 -7.4127E-02  2.9803E-04 -1.2131E-02 -1.2222E-02
            -3.0835E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3617.40854394391        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      712
 NPARAMETR:  1.0298E+00  8.6826E-01  8.3872E-01  1.0913E+00  7.9076E-01  8.7228E-01  1.2863E+00  1.0000E-02  7.9084E-01  1.0672E+00
             1.3538E+00
 PARAMETER:  1.2941E-01 -4.1268E-02 -7.5882E-02  1.8734E-01 -1.3476E-01 -3.6645E-02  3.5176E-01 -4.5462E+00 -1.3467E-01  1.6505E-01
             4.0288E-01
 GRADIENT:   9.9631E-03  8.2960E-04 -8.7391E-03  6.0100E-03  5.9951E-05  9.1569E-04 -2.5164E-03  0.0000E+00 -2.0877E-03 -2.5668E-04
             7.5509E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -3617.40854394391        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  1.0298E+00  8.6826E-01  8.3872E-01  1.0913E+00  7.9076E-01  8.7228E-01  1.2863E+00  1.0000E-02  7.9084E-01  1.0672E+00
             1.3538E+00
 PARAMETER:  1.2941E-01 -4.1268E-02 -7.5882E-02  1.8734E-01 -1.3476E-01 -3.6645E-02  3.5176E-01 -4.5462E+00 -1.3467E-01  1.6505E-01
             4.0288E-01
 GRADIENT:   9.9631E-03  8.2960E-04 -8.7391E-03  6.0100E-03  5.9951E-05  9.1569E-04 -2.5164E-03  0.0000E+00 -2.0877E-03 -2.5668E-04
             7.5509E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      734
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4863E-04 -1.5551E-02 -1.7331E-04  8.2603E-03 -1.4193E-02
 SE:             2.9791E-02  2.3734E-02  2.2746E-04  2.7132E-02  2.5915E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9066E-01  5.1234E-01  4.4609E-01  7.6079E-01  5.8391E-01

 ETASHRINKSD(%)  1.9806E-01  2.0487E+01  9.9238E+01  9.1046E+00  1.3182E+01
 ETASHRINKVR(%)  3.9572E-01  3.6776E+01  9.9994E+01  1.7380E+01  2.4627E+01
 EBVSHRINKSD(%)  5.5982E-01  1.9605E+01  9.9255E+01  1.0725E+01  1.3505E+01
 EBVSHRINKVR(%)  1.1165E+00  3.5367E+01  9.9994E+01  2.0300E+01  2.5187E+01
 RELATIVEINF(%)  9.8872E+01  2.9951E+01  2.3425E-03  5.2243E+01  2.2712E+01
 EPSSHRINKSD(%)  1.9385E+01
 EPSSHRINKVR(%)  3.5012E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3617.4085439439141     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1963.3191841755033     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.75
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3617.409       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  8.68E-01  8.39E-01  1.09E+00  7.91E-01  8.72E-01  1.29E+00  1.00E-02  7.91E-01  1.07E+00  1.35E+00
 


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
+        1.36E+03
 
 TH 2
+       -5.45E+00  5.78E+02
 
 TH 3
+       -2.25E+00 -6.36E+01  8.23E+02
 
 TH 4
+       -8.67E+00  2.64E+02 -6.53E+01  1.17E+03
 
 TH 5
+       -4.52E+00 -3.36E+02 -6.07E+02  2.77E+02  1.15E+03
 
 TH 6
+       -1.60E+00 -1.55E-01  5.83E-01 -2.72E+00 -2.55E+00  2.58E+02
 
 TH 7
+        2.73E-01  5.61E+00 -3.19E+00  1.94E+01 -1.50E+01 -1.35E-01  4.46E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.68E-01  8.70E+00  4.02E+01 -4.80E+00  2.53E+01  1.85E+00  2.53E+01  0.00E+00  1.94E+02
 
 TH10
+        8.27E-01 -1.98E+01 -4.55E+01  9.90E+00 -2.54E+01 -6.28E-01  6.83E+00  0.00E+00  2.65E+00  1.06E+02
 
 TH11
+       -1.14E+01 -2.00E+01 -2.11E+01 -1.67E+01 -4.22E+00  1.88E+00  9.60E+00  0.00E+00  1.35E+01  1.15E+01  6.06E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.938
Stop Time:
Fri Sep 24 20:22:15 CDT 2021
