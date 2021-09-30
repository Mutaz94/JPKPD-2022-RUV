Wed Sep 29 18:36:24 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat92.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.91683350375        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4719E+02 -4.5098E+01 -6.3177E+01  2.7936E+01  8.3957E+01  5.4726E+01 -1.5620E+01  9.5711E+00 -8.9735E+00  1.6218E+01
             2.3660E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1660.45136506338        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.8677E-01  1.1191E+00  1.2433E+00  9.5561E-01  1.0769E+00  8.9864E-01  1.1488E+00  9.3225E-01  1.1482E+00  8.3533E-01
             9.9806E-01
 PARAMETER:  8.6683E-02  2.1255E-01  3.1779E-01  5.4598E-02  1.7413E-01 -6.8683E-03  2.3876E-01  2.9843E-02  2.3819E-01 -7.9932E-02
             9.8061E-02
 GRADIENT:   1.3889E+01 -1.1724E+01  1.4967E+01 -1.1556E+01  2.0427E+01 -2.6068E+01  8.2689E-02 -8.7385E+00  9.0297E+00 -1.8427E+01
            -1.0289E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1661.28836108270        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  9.8670E-01  1.3260E+00  1.0777E+00  8.4886E-01  1.0803E+00  9.2949E-01  1.0922E+00  1.0034E+00  1.2707E+00  7.9247E-01
             1.0156E+00
 PARAMETER:  8.6609E-02  3.8220E-01  1.7481E-01 -6.3862E-02  1.7720E-01  2.6880E-02  1.8820E-01  1.0339E-01  3.3959E-01 -1.3260E-01
             1.1546E-01
 GRADIENT:   9.1159E+00  2.4926E+01  1.5985E+01  1.4091E+01 -9.4858E+00 -1.2323E+01  6.3348E+00 -3.5541E+00  1.1849E+01 -1.8038E+01
            -2.3383E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1663.85022443124        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  9.8516E-01  1.2646E+00  9.9647E-01  8.7315E-01  1.0418E+00  9.6395E-01  1.1037E+00  8.8369E-01  1.1312E+00  8.9016E-01
             1.0067E+00
 PARAMETER:  8.5047E-02  3.3476E-01  9.6461E-02 -3.5651E-02  1.4094E-01  6.3288E-02  1.9869E-01 -2.3653E-02  2.2326E-01 -1.6354E-02
             1.0665E-01
 GRADIENT:   4.2709E+00  6.3410E+00  2.6334E+00  4.0563E+00 -8.1672E+00  2.4520E+00  1.1050E+00  1.5542E-01  4.5735E-01  1.4029E+00
             3.1739E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.13030032146        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.8196E-01  1.5491E+00  7.5216E-01  6.9340E-01  1.0802E+00  9.4831E-01  9.5057E-01  6.6571E-01  1.3082E+00  8.7290E-01
             1.0057E+00
 PARAMETER:  8.1794E-02  5.3769E-01 -1.8480E-01 -2.6614E-01  1.7710E-01  4.6922E-02  4.9306E-02 -3.0690E-01  3.6867E-01 -3.5932E-02
             1.0564E-01
 GRADIENT:  -7.2991E+00  2.0625E+01  8.9571E-01  1.3771E+01 -5.9867E+00 -4.6517E+00 -1.7517E+00  4.1328E-01 -2.6897E+00 -1.5674E+00
            -1.3734E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.80236193189        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.8804E-01  1.7946E+00  5.8142E-01  5.2201E-01  1.1506E+00  9.6743E-01  8.7176E-01  3.7322E-01  1.6197E+00  9.2326E-01
             1.0119E+00
 PARAMETER:  8.7969E-02  6.8476E-01 -4.4227E-01 -5.5008E-01  2.4026E-01  6.6886E-02 -3.7239E-02 -8.8558E-01  5.8222E-01  2.0151E-02
             1.1182E-01
 GRADIENT:   6.4783E+00  7.3987E+00 -1.5844E+00  7.4207E+00 -7.3760E+00  3.0727E+00  2.8007E+00  3.1789E-01  1.9229E+00  2.7692E+00
             2.0498E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1665.08260756292        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  9.8560E-01  1.8577E+00  5.5029E-01  4.7144E-01  1.1846E+00  9.5986E-01  8.4096E-01  1.8880E-01  1.7248E+00  9.2801E-01
             1.0095E+00
 PARAMETER:  8.5500E-02  7.1932E-01 -4.9731E-01 -6.5197E-01  2.6941E-01  5.9032E-02 -7.3217E-02 -1.5671E+00  6.4510E-01  2.5282E-02
             1.0944E-01
 GRADIENT:   7.8736E-01 -4.9149E+00  1.4172E-01  1.5047E+00 -3.6568E-01  6.3200E-02  1.1032E-01  6.3359E-02 -2.1030E-01 -2.5286E-01
            -6.6610E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1665.11676269615        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  9.8577E-01  1.8672E+00  5.4355E-01  4.6234E-01  1.1888E+00  9.5993E-01  8.3722E-01  4.7082E-02  1.7455E+00  9.3160E-01
             1.0096E+00
 PARAMETER:  8.5672E-02  7.2442E-01 -5.0964E-01 -6.7146E-01  2.7291E-01  5.9100E-02 -7.7672E-02 -2.9559E+00  6.5701E-01  2.9144E-02
             1.0956E-01
 GRADIENT:   1.2420E+00 -8.6387E+00  5.2991E-01 -3.2879E-01 -9.2345E-01  1.0463E-01 -5.3158E-03  4.0537E-03 -3.5282E-01 -2.3492E-02
             1.3983E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1665.12090842458        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1425
 NPARAMETR:  9.8593E-01  1.8652E+00  5.4226E-01  4.6190E-01  1.1890E+00  9.5997E-01  8.3701E-01  1.0000E-02  1.7503E+00  9.3160E-01
             1.0095E+00
 PARAMETER:  8.5828E-02  7.2339E-01 -5.1200E-01 -6.7240E-01  2.7311E-01  5.9146E-02 -7.7921E-02 -4.9364E+00  6.5977E-01  2.9146E-02
             1.0944E-01
 GRADIENT:   1.6364E+00 -1.2332E+01 -1.4088E-01 -5.7512E-01  4.8004E-01  1.1906E-01  5.6846E-02  0.0000E+00  1.7913E-01  9.0572E-02
             8.3139E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1665.12090842458        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1447
 NPARAMETR:  9.8593E-01  1.8652E+00  5.4226E-01  4.6190E-01  1.1890E+00  9.5997E-01  8.3701E-01  1.0000E-02  1.7503E+00  9.3160E-01
             1.0095E+00
 PARAMETER:  8.5828E-02  7.2339E-01 -5.1200E-01 -6.7240E-01  2.7311E-01  5.9146E-02 -7.7921E-02 -4.9364E+00  6.5977E-01  2.9146E-02
             1.0944E-01
 GRADIENT:   1.6364E+00 -1.2332E+01 -1.4088E-01 -5.7512E-01  4.8004E-01  1.1906E-01  5.6846E-02  0.0000E+00  1.7913E-01  9.0572E-02
             8.3139E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1447
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.7669E-05 -3.1624E-02 -2.3366E-04  3.1935E-02 -4.0646E-02
 SE:             2.9850E-02  2.4395E-02  8.6853E-05  2.2197E-02  2.1869E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9766E-01  1.9486E-01  7.1388E-03  1.5024E-01  6.3080E-02

 ETASHRINKSD(%)  1.0000E-10  1.8274E+01  9.9709E+01  2.5636E+01  2.6737E+01
 ETASHRINKVR(%)  1.0000E-10  3.3209E+01  9.9999E+01  4.4700E+01  4.6326E+01
 EBVSHRINKSD(%)  4.6966E-01  1.7098E+01  9.9752E+01  2.8333E+01  2.5232E+01
 EBVSHRINKVR(%)  9.3712E-01  3.1273E+01  9.9999E+01  4.8638E+01  4.4097E+01
 RELATIVEINF(%)  9.9017E+01  7.0045E+00  1.0287E-04  5.0867E+00  1.5804E+01
 EPSSHRINKSD(%)  4.4223E+01
 EPSSHRINKVR(%)  6.8889E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.1209084245804     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -929.97008186084224     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.27
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.121       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.87E+00  5.42E-01  4.62E-01  1.19E+00  9.60E-01  8.37E-01  1.00E-02  1.75E+00  9.32E-01  1.01E+00
 


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
+        1.23E+03
 
 TH 2
+       -5.23E+00  3.31E+02
 
 TH 3
+        7.51E+00  1.15E+02  2.71E+02
 
 TH 4
+       -1.42E+01  2.69E+02 -2.53E+02  9.16E+02
 
 TH 5
+       -5.42E+00 -1.46E+02 -2.31E+02  2.54E+02  5.08E+02
 
 TH 6
+       -2.71E-01 -6.54E-01  1.81E+00 -3.41E+00 -7.28E-01  2.13E+02
 
 TH 7
+        7.28E-01  7.22E+00 -8.69E+00 -1.43E+01 -1.49E+01 -5.47E-01  1.41E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.27E-01 -1.46E+01 -2.83E+01  5.79E+01  8.36E-01 -4.37E-01  1.47E+01  0.00E+00  2.48E+01
 
 TH10
+        2.70E-02 -1.38E+01 -2.51E+01 -3.90E+00 -6.50E+01  6.75E-01  9.60E+00  0.00E+00  5.07E+00  8.85E+01
 
 TH11
+       -7.59E+00 -1.59E+01 -1.82E+01  1.39E+00 -7.69E+00  2.68E+00  8.22E+00  0.00E+00  3.79E+00  1.78E+01  2.11E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.780
Stop Time:
Wed Sep 29 18:36:52 CDT 2021
