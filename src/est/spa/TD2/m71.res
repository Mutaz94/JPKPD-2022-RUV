Sat Sep 25 13:44:08 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat71.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1620.02488193217        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4604E+02 -4.1096E+01  5.1080E+01 -1.0693E+02 -3.7136E+01  9.6695E+00 -1.6355E+01 -1.2860E+01 -6.2474E+00 -9.9500E+00
             2.3424E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1632.14669973015        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.5743E-01  1.0788E+00  8.3614E-01  1.0075E+00  9.7106E-01  9.6869E-01  1.2292E+00  1.0777E+00  9.4628E-01  1.0644E+00
             9.0233E-01
 PARAMETER:  5.6501E-02  1.7588E-01 -7.8957E-02  1.0746E-01  7.0635E-02  6.8191E-02  3.0638E-01  1.7479E-01  4.4780E-02  1.6237E-01
            -2.7734E-03
 GRADIENT:   5.5357E+01  4.0602E+00  4.3606E-02  3.0560E+00 -1.0468E+01  2.6179E+00  2.9092E+00  5.4817E+00  4.0990E+00  1.0420E+01
            -3.9830E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1632.80961290489        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.5636E-01  1.0670E+00  6.9112E-01  1.0013E+00  8.6043E-01  9.6341E-01  1.2731E+00  9.0494E-01  8.8776E-01  8.8967E-01
             8.9132E-01
 PARAMETER:  5.5382E-02  1.6481E-01 -2.6944E-01  1.0128E-01 -5.0325E-02  6.2719E-02  3.4144E-01  1.1040E-04 -1.9052E-02 -1.6903E-02
            -1.5051E-02
 GRADIENT:   5.0280E+01  8.4393E+00  3.9595E+00  2.0361E+00 -1.4412E+01 -3.2287E-03  4.6894E+00  5.0523E+00 -1.2196E+00  4.1250E+00
            -8.3752E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1633.44412214750        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:      240
 NPARAMETR:  9.4384E-01  1.0961E+00  5.1452E-01  9.5958E-01  7.5517E-01  9.6471E-01  1.2269E+00  5.5315E-01  8.7900E-01  7.9458E-01
             8.9885E-01
 PARAMETER:  4.2202E-02  1.9173E-01 -5.6451E-01  5.8745E-02 -1.8081E-01  6.4070E-02  3.0452E-01 -4.9212E-01 -2.8973E-02 -1.2994E-01
            -6.6345E-03
 GRADIENT:   1.5070E+01  6.8170E+00 -2.5897E+00  9.9275E+00 -4.7085E+00 -1.6601E-01  4.3465E+00  3.0595E+00  5.5766E-01  5.7456E+00
            -2.1906E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1634.73648580540        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      424
 NPARAMETR:  9.6233E-01  1.2620E+00  4.3898E-01  8.5060E-01  7.8942E-01  9.7922E-01  1.0900E+00  3.5260E-01  9.3790E-01  7.5979E-01
             9.0770E-01
 PARAMETER:  6.1601E-02  3.3271E-01 -7.2331E-01 -6.1808E-02 -1.3645E-01  7.8996E-02  1.8618E-01 -9.4241E-01  3.5883E-02 -1.7472E-01
             3.1561E-03
 GRADIENT:   1.6944E+01  3.6436E-01 -1.7303E+00  1.7050E+00 -1.0951E-01  1.0552E+00  3.7468E-01  8.2573E-01 -1.6727E+00 -1.3644E-01
            -2.8514E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1635.05207901252        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      601
 NPARAMETR:  9.5460E-01  1.3584E+00  4.1810E-01  7.9081E-01  8.3245E-01  9.7630E-01  1.0299E+00  2.1498E-01  9.9708E-01  7.9976E-01
             9.1008E-01
 PARAMETER:  5.3536E-02  4.0634E-01 -7.7204E-01 -1.3470E-01 -8.3385E-02  7.6017E-02  1.2950E-01 -1.4372E+00  9.7071E-02 -1.2345E-01
             5.7770E-03
 GRADIENT:  -1.4432E+00 -8.5722E-01 -6.9944E-01 -1.3371E+00 -6.8178E-01  1.0113E-02 -1.7521E-01  2.7562E-01  6.5675E-01  7.2421E-01
             4.2521E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1635.21850226792        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      777
 NPARAMETR:  9.5521E-01  1.3053E+00  4.3751E-01  8.2531E-01  8.1655E-01  9.7682E-01  1.0672E+00  4.6655E-02  9.6874E-01  7.9873E-01
             9.0892E-01
 PARAMETER:  5.4175E-02  3.6642E-01 -7.2665E-01 -9.1994E-02 -1.0266E-01  7.6551E-02  1.6501E-01 -2.9650E+00  6.8238E-02 -1.2473E-01
             4.5034E-03
 GRADIENT:  -2.9584E-01  3.4396E-01  5.9808E-01 -4.6562E-01 -5.0553E-01  1.7714E-01 -2.0727E-01  1.1639E-02 -1.7192E-02 -2.3895E-01
             1.1658E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1635.22523774052        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      952
 NPARAMETR:  9.5530E-01  1.3011E+00  4.3650E-01  8.2754E-01  8.1342E-01  9.7638E-01  1.0700E+00  1.0000E-02  9.6605E-01  7.9741E-01
             9.0837E-01
 PARAMETER:  5.4268E-02  3.6317E-01 -7.2897E-01 -8.9292E-02 -1.0651E-01  7.6097E-02  1.6762E-01 -4.6366E+00  6.5455E-02 -1.2638E-01
             3.8930E-03
 GRADIENT:  -4.2021E-02 -1.1789E-03 -8.0450E-03  1.2437E-02  9.9122E-03  8.8121E-04  3.3215E-03  0.0000E+00  6.0344E-04  1.7382E-03
             5.7294E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1635.22523774052        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      974
 NPARAMETR:  9.5530E-01  1.3011E+00  4.3650E-01  8.2754E-01  8.1342E-01  9.7638E-01  1.0700E+00  1.0000E-02  9.6605E-01  7.9741E-01
             9.0837E-01
 PARAMETER:  5.4268E-02  3.6317E-01 -7.2897E-01 -8.9292E-02 -1.0651E-01  7.6097E-02  1.6762E-01 -4.6366E+00  6.5455E-02 -1.2638E-01
             3.8930E-03
 GRADIENT:  -4.2021E-02 -1.1789E-03 -8.0450E-03  1.2437E-02  9.9122E-03  8.8121E-04  3.3215E-03  0.0000E+00  6.0344E-04  1.7382E-03
             5.7294E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      974
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4137E-04 -1.1624E-02 -4.0005E-04  9.2148E-03 -1.9185E-02
 SE:             2.9872E-02  2.4982E-02  1.5557E-04  2.4663E-02  2.0736E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9622E-01  6.4174E-01  1.0128E-02  7.0868E-01  3.5486E-01

 ETASHRINKSD(%)  1.0000E-10  1.6306E+01  9.9479E+01  1.7375E+01  3.0533E+01
 ETASHRINKVR(%)  1.0000E-10  2.9953E+01  9.9997E+01  3.1731E+01  5.1743E+01
 EBVSHRINKSD(%)  3.5916E-01  1.6011E+01  9.9529E+01  1.8098E+01  3.0449E+01
 EBVSHRINKVR(%)  7.1703E-01  2.9459E+01  9.9998E+01  3.2921E+01  5.1627E+01
 RELATIVEINF(%)  9.9194E+01  5.4704E+00  2.0065E-04  4.7818E+00  4.4250E+00
 EPSSHRINKSD(%)  4.5954E+01
 EPSSHRINKVR(%)  7.0791E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1635.2252377405155     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -900.07441117677729     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.78
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1635.225       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.55E-01  1.30E+00  4.36E-01  8.28E-01  8.13E-01  9.76E-01  1.07E+00  1.00E-02  9.66E-01  7.97E-01  9.08E-01
 


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
+        1.27E+03
 
 TH 2
+       -5.37E+00  4.14E+02
 
 TH 3
+       -5.37E+00  3.15E+02  1.38E+03
 
 TH 4
+       -9.94E+00  2.52E+02 -7.80E+02  1.17E+03
 
 TH 5
+       -3.15E+00 -3.92E+02 -1.15E+03  6.10E+02  1.27E+03
 
 TH 6
+        1.34E+00 -1.05E+00  2.47E-01 -1.36E+00 -5.49E-01  2.05E+02
 
 TH 7
+        1.88E-01  1.85E+01 -6.65E+01 -6.14E+00  9.66E+00 -2.26E-01  8.61E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.73E+00 -1.76E+01 -5.09E+01  5.08E+01 -5.91E+00 -3.56E-01  1.93E+01  0.00E+00  1.00E+02
 
 TH10
+       -1.35E+00 -1.59E+01 -6.80E+01 -2.18E+01 -7.40E+01 -3.12E-02  1.40E+01  0.00E+00  1.88E+01  8.66E+01
 
 TH11
+       -7.41E+00 -1.24E+01 -3.33E+01 -2.25E-01 -7.67E+00  1.58E+00  7.68E+00  0.00E+00  8.06E+00  1.83E+01  2.46E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.413
Stop Time:
Sat Sep 25 13:44:27 CDT 2021
