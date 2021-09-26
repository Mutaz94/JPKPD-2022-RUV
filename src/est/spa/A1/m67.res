Sat Sep 25 08:11:05 CDT 2021
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
$DATA ../../../../data/spa/A1/dat67.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1271.78280702962        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.5052E+01 -7.3727E+01  9.4885E-01 -1.2167E+02  4.6125E+01  5.7171E+00 -3.0657E+01 -6.3037E+00 -6.5203E+01 -2.0012E+01
            -6.7566E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1466.98640492428        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0502E+00  1.1039E+00  1.1634E+00  1.0335E+00  1.1100E+00  9.9335E-01  1.0228E+00  9.3583E-01  9.8367E-01  8.1936E-01
             2.3626E+00
 PARAMETER:  1.4895E-01  1.9888E-01  2.5136E-01  1.3293E-01  2.0439E-01  9.3326E-02  1.2258E-01  3.3680E-02  8.3535E-02 -9.9233E-02
             9.5978E-01
 GRADIENT:   1.0977E+02 -3.0429E+01 -9.6881E+00 -3.4283E+01  6.4142E+00 -3.6780E+00  1.9682E-01  3.3092E+00 -3.3945E-01  8.6455E+00
             2.4383E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1469.45067555431        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0367E+00  1.3486E+00  8.8924E-01  9.0841E-01  1.0547E+00  9.9832E-01  9.3281E-01  7.5271E-01  1.1350E+00  6.2013E-01
             2.3244E+00
 PARAMETER:  1.3605E-01  3.9910E-01 -1.7391E-02  3.9391E-03  1.5325E-01  9.8314E-02  3.0443E-02 -1.8407E-01  2.2666E-01 -3.7782E-01
             9.4346E-01
 GRADIENT:   7.4531E+01  2.5246E+01  4.9390E+00  1.2362E+01 -2.5200E+01 -8.5433E-01  1.0718E+00  1.5656E+00  6.1436E+00  3.4697E+00
             9.4679E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1471.97281809782        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9635E-01  1.0699E+00  1.1735E+00  1.0804E+00  1.0755E+00  9.9989E-01  1.0743E+00  4.2111E-01  9.6753E-01  5.8276E-01
             2.3331E+00
 PARAMETER:  9.6347E-02  1.6758E-01  2.5997E-01  1.7735E-01  1.7282E-01  9.9885E-02  1.7166E-01 -7.6486E-01  6.6987E-02 -4.3998E-01
             9.4719E-01
 GRADIENT:  -5.8178E+00  5.7867E+00  2.4943E+00  3.8774E+00 -5.0820E+00  2.4617E+00  4.7139E-01  1.8528E-01 -9.9987E-01  3.7187E-01
            -2.4592E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1472.14917562333        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.9898E-01  1.0250E+00  1.1850E+00  1.1037E+00  1.0681E+00  9.9205E-01  1.0803E+00  2.0728E-01  9.6542E-01  5.3814E-01
             2.3625E+00
 PARAMETER:  9.8982E-02  1.2469E-01  2.6973E-01  1.9863E-01  1.6589E-01  9.2019E-02  1.7728E-01 -1.4737E+00  6.4806E-02 -5.1964E-01
             9.5974E-01
 GRADIENT:   3.7610E-01 -1.8703E+00 -7.6313E-01 -1.2720E+00  1.3506E+00 -6.6269E-02 -1.0675E-01  4.5098E-02  4.7714E-02  1.2418E-02
             3.9704E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1472.19577332437        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      418
 NPARAMETR:  1.0014E+00  1.0132E+00  1.2236E+00  1.1162E+00  1.0784E+00  9.9394E-01  1.0797E+00  9.5837E-02  9.6357E-01  5.5169E-01
             2.3698E+00
 PARAMETER:  1.0142E-01  1.1308E-01  3.0178E-01  2.0990E-01  1.7545E-01  9.3920E-02  1.7664E-01 -2.2451E+00  6.2888E-02 -4.9476E-01
             9.6281E-01
 GRADIENT:  -1.3867E+00  4.5803E-01  7.0607E-02 -7.3037E-02 -2.5637E-01 -1.1315E-01  1.0933E-01  9.6471E-03  4.7463E-02  1.4572E-02
             4.4590E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1472.24091337538        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      594
 NPARAMETR:  1.0003E+00  8.4952E-01  1.3749E+00  1.2260E+00  1.0713E+00  9.9263E-01  1.1394E+00  1.0000E-02  9.1925E-01  5.8740E-01
             2.3768E+00
 PARAMETER:  1.0033E-01 -6.3082E-02  4.1837E-01  3.0377E-01  1.6889E-01  9.2603E-02  2.3052E-01 -7.3018E+00  1.5798E-02 -4.3205E-01
             9.6576E-01
 GRADIENT:  -2.6495E-01  1.6232E+00  6.1352E-01  2.5200E+00 -1.2418E+00 -6.6691E-02 -1.1073E-01  0.0000E+00 -2.9717E-02 -4.3873E-02
             2.7983E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1472.25686250655        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      770
 NPARAMETR:  9.9930E-01  7.4209E-01  1.4254E+00  1.2941E+00  1.0508E+00  9.9186E-01  1.2288E+00  1.0000E-02  8.8709E-01  6.0824E-01
             2.3713E+00
 PARAMETER:  9.9298E-02 -1.9829E-01  4.5445E-01  3.5784E-01  1.4954E-01  9.1827E-02  3.0606E-01 -1.2828E+01 -1.9812E-02 -3.9719E-01
             9.6343E-01
 GRADIENT:  -6.8521E-02  9.8276E-01  3.7842E-01  1.7661E+00 -8.6834E-01 -5.1195E-02  1.1972E-01  0.0000E+00 -1.1166E-01 -2.2953E-02
            -1.3873E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1472.26142341918        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      945
 NPARAMETR:  9.9869E-01  6.9045E-01  1.4625E+00  1.3265E+00  1.0470E+00  9.9151E-01  1.2497E+00  1.0000E-02  8.7657E-01  6.2029E-01
             2.3712E+00
 PARAMETER:  9.8689E-02 -2.7041E-01  4.8018E-01  3.8253E-01  1.4588E-01  9.1470E-02  3.2291E-01 -1.6405E+01 -3.1742E-02 -3.7758E-01
             9.6338E-01
 GRADIENT:  -4.9397E-02  9.1517E-02  2.8549E-02  2.1080E-01 -7.2644E-02 -2.8551E-03 -3.7370E-02  0.0000E+00 -3.0554E-02 -2.6071E-03
            -4.2008E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1472.26168130364        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1113
 NPARAMETR:  9.9849E-01  6.6957E-01  1.4727E+00  1.3399E+00  1.0433E+00  9.9131E-01  1.2707E+00  1.0000E-02  8.7087E-01  6.2216E-01
             2.3712E+00
 PARAMETER:  9.8474E-02 -3.0141E-01  4.8701E-01  3.9249E-01  1.4247E-01  9.1287E-02  3.3931E-01 -1.7910E+01 -3.8209E-02 -3.7468E-01
             9.6341E-01
 GRADIENT:  -3.1162E-03 -1.0068E-02 -4.1195E-03 -1.8087E-02  9.6148E-03  5.8385E-04 -6.6226E-04  0.0000E+00  9.9416E-04 -9.8118E-05
             1.2906E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1113
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.1841E-05 -1.2541E-03  8.9944E-05 -9.8782E-03 -1.5510E-02
 SE:             2.9236E-02  1.3363E-02  1.0125E-04  2.4161E-02  1.2932E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9749E-01  9.2523E-01  3.7438E-01  6.8265E-01  2.3040E-01

 ETASHRINKSD(%)  2.0564E+00  5.5232E+01  9.9661E+01  1.9059E+01  5.6676E+01
 ETASHRINKVR(%)  4.0705E+00  7.9958E+01  9.9999E+01  3.4485E+01  8.1231E+01
 EBVSHRINKSD(%)  2.0687E+00  5.6337E+01  9.9645E+01  1.8519E+01  5.7851E+01
 EBVSHRINKVR(%)  4.0947E+00  8.0935E+01  9.9999E+01  3.3609E+01  8.2234E+01
 RELATIVEINF(%)  9.1800E+01  2.1106E-01  8.9038E-05  9.6319E-01  1.0806E+00
 EPSSHRINKSD(%)  2.8070E+01
 EPSSHRINKVR(%)  4.8261E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1472.2616813036393     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -737.11085473990113     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.96
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
 





 #OBJV:********************************************    -1472.262       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.98E-01  6.69E-01  1.47E+00  1.34E+00  1.04E+00  9.91E-01  1.27E+00  1.00E-02  8.71E-01  6.22E-01  2.37E+00
 


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
+        1.08E+03
 
 TH 2
+       -3.59E+01  2.87E+02
 
 TH 3
+        5.68E+00  5.53E+01  6.12E+01
 
 TH 4
+       -3.94E+01  3.53E+02 -1.46E+00  5.44E+02
 
 TH 5
+       -4.57E-01 -1.84E+02 -1.54E+02 -5.73E+01  4.37E+02
 
 TH 6
+       -2.73E-01 -5.47E+00  3.06E+00 -1.12E+01 -3.52E+00  1.83E+02
 
 TH 7
+        5.81E-01  1.54E+00  2.82E+00 -2.73E+00 -3.72E-01 -8.28E-03  6.46E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.55E-01 -1.60E+01  6.00E+00  1.45E-01  2.63E+00  4.66E+00  1.58E+01  0.00E+00  1.15E+02
 
 TH10
+       -3.16E+00  5.50E-02  1.55E-01 -7.09E+00 -5.64E+00 -1.47E-01  1.82E+00  0.00E+00  2.83E+00  1.12E+01
 
 TH11
+       -1.27E+01 -2.56E+00 -8.00E-01 -1.03E+01 -1.23E+01  2.88E+00  2.96E+00  0.00E+00  8.73E+00  2.04E+01  6.12E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.296
Stop Time:
Sat Sep 25 08:11:26 CDT 2021
