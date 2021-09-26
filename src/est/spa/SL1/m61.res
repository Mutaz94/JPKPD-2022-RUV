Sat Sep 25 10:37:41 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat61.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1703.17818063649        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.4233E+01 -5.0868E+01  2.8860E-01 -8.3729E+01 -2.6323E+01  8.5917E+00  8.9559E+00  1.2543E+01  1.6714E+01  9.8833E+00
            -1.2315E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1709.11281857722        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      132
 NPARAMETR:  1.0418E+00  1.0649E+00  1.0725E+00  1.0239E+00  1.0788E+00  9.5789E-01  9.2823E-01  9.0706E-01  8.7216E-01  9.5996E-01
             1.0723E+00
 PARAMETER:  1.4090E-01  1.6289E-01  1.6996E-01  1.2360E-01  1.7585E-01  5.6973E-02  2.5520E-02  2.4529E-03 -3.6777E-02  5.9132E-02
             1.6980E-01
 GRADIENT:   2.2964E+01  3.0498E+00  6.2109E+00 -8.2859E+00  5.9529E+00 -1.1583E+01 -3.1672E-01  2.9617E+00 -1.2114E+01 -1.0898E+01
             6.9274E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1711.12834350295        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.0368E+00  1.0128E+00  1.0120E+00  1.0542E+00  1.0447E+00  9.5912E-01  7.2681E-01  5.2489E-01  9.5901E-01  1.0889E+00
             1.0836E+00
 PARAMETER:  1.3610E-01  1.1267E-01  1.1189E-01  1.5277E-01  1.4377E-01  5.8258E-02 -2.1908E-01 -5.4456E-01  5.8149E-02  1.8521E-01
             1.8025E-01
 GRADIENT:   1.0936E+01 -6.7332E+00 -1.0095E+01  3.4670E+00  1.1589E+01 -1.0472E+01 -2.6592E+00  1.1113E+00 -4.0571E-01 -9.2553E-02
             1.5241E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1712.22461230241        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      485
 NPARAMETR:  1.0343E+00  1.1269E+00  8.8788E-01  9.7740E-01  1.0295E+00  9.8708E-01  8.7937E-01  3.9968E-01  9.6855E-01  1.0309E+00
             1.0410E+00
 PARAMETER:  1.3373E-01  2.1947E-01 -1.8914E-02  7.7138E-02  1.2908E-01  8.6996E-02 -2.8548E-02 -8.1708E-01  6.8049E-02  1.3045E-01
             1.4017E-01
 GRADIENT:   2.3372E+00 -5.4570E+00 -7.9895E+00  2.3283E+00  9.7256E+00 -8.8818E-02  1.3298E+00  9.8506E-01  1.5448E+00  2.3567E+00
             1.1180E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1712.59033950157        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.0338E+00  1.0317E+00  8.6572E-01  1.0336E+00  9.5881E-01  9.8353E-01  9.6440E-01  1.5083E-01  9.0587E-01  9.7274E-01
             1.0427E+00
 PARAMETER:  1.3320E-01  1.3118E-01 -4.4198E-02  1.3306E-01  5.7942E-02  8.3395E-02  6.3753E-02 -1.7916E+00  1.1355E-03  7.2366E-02
             1.4178E-01
 GRADIENT:   1.3926E+00 -6.8675E-01 -1.3951E+00  5.4373E-01  1.2812E+00 -1.3904E+00 -1.3944E-01  1.0760E-01 -2.5879E-01 -2.0214E-01
             1.0077E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.64198325481        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  1.0329E+00  1.0415E+00  8.7020E-01  1.0286E+00  9.6709E-01  9.8807E-01  9.5296E-01  3.2077E-02  9.1494E-01  9.9063E-01
             1.0405E+00
 PARAMETER:  1.3237E-01  1.4065E-01 -3.9027E-02  1.2823E-01  6.6539E-02  8.7997E-02  5.1816E-02 -3.3396E+00  1.1103E-02  9.0590E-02
             1.3967E-01
 GRADIENT:  -3.9580E-01  4.4051E-01 -4.4458E-01  1.3803E+00  8.1659E-02  4.4401E-01  1.5170E-01  4.2028E-03  2.9383E-01  6.8019E-01
            -1.1079E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1712.64995610748        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1010
 NPARAMETR:  1.0334E+00  1.0872E+00  8.5393E-01  9.9911E-01  9.8017E-01  9.8733E-01  9.1886E-01  1.3364E-02  9.3589E-01  9.8800E-01
             1.0411E+00
 PARAMETER:  1.3281E-01  1.8358E-01 -5.7904E-02  9.9110E-02  7.9974E-02  8.7245E-02  1.5377E-02 -4.2152E+00  3.3745E-02  8.7923E-02
             1.4026E-01
 GRADIENT:  -9.0702E-02  5.2132E-02 -8.3677E-02  1.5165E-01  1.8367E-02 -3.1852E-02 -1.3652E-02  7.1943E-04  7.7448E-03  3.1462E-02
             2.2126E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1712.65013488930        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  1.0334E+00  1.0852E+00  8.5467E-01  1.0003E+00  9.7958E-01  9.8734E-01  9.2049E-01  1.0000E-02  9.3487E-01  9.8773E-01
             1.0411E+00
 PARAMETER:  1.3283E-01  1.8174E-01 -5.7039E-02  1.0033E-01  7.9370E-02  8.7261E-02  1.7146E-02 -4.5716E+00  3.2652E-02  8.7655E-02
             1.4025E-01
 GRADIENT:  -2.6625E-02  2.9207E-02 -1.6068E-02  6.2250E-02 -2.2682E-04 -1.7404E-02 -4.9821E-03  0.0000E+00  2.8308E-03  7.9491E-04
             4.8315E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1712.65013820466        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  1.0334E+00  1.0857E+00  8.5463E-01  9.9997E-01  9.7985E-01  9.8738E-01  9.2010E-01  1.0000E-02  9.3515E-01  9.8789E-01
             1.0411E+00
 PARAMETER:  1.3284E-01  1.8223E-01 -5.7089E-02  9.9974E-02  7.9643E-02  8.7303E-02  1.6731E-02 -4.5717E+00  3.2954E-02  8.7813E-02
             1.4025E-01
 GRADIENT:  -2.2749E-03  7.2213E-04 -2.2825E-03  3.8071E-03  1.2964E-03 -1.9150E-03 -3.0628E-04  0.0000E+00  6.9113E-04 -7.4558E-05
             7.2562E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1242
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.6947E-04 -1.4492E-02 -3.8359E-04  3.5723E-03 -2.4679E-02
 SE:             2.9831E-02  1.8623E-02  1.6848E-04  2.4969E-02  2.3672E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9279E-01  4.3646E-01  2.2801E-02  8.8623E-01  2.9715E-01

 ETASHRINKSD(%)  6.2593E-02  3.7611E+01  9.9436E+01  1.6351E+01  2.0696E+01
 ETASHRINKVR(%)  1.2515E-01  6.1076E+01  9.9997E+01  3.0028E+01  3.7109E+01
 EBVSHRINKSD(%)  4.6197E-01  3.7163E+01  9.9472E+01  1.6741E+01  1.8913E+01
 EBVSHRINKVR(%)  9.2180E-01  6.0515E+01  9.9997E+01  3.0679E+01  3.4250E+01
 RELATIVEINF(%)  9.8558E+01  1.2672E+00  3.1581E-04  2.8530E+00  6.5311E+00
 EPSSHRINKSD(%)  4.2582E+01
 EPSSHRINKVR(%)  6.7032E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1712.6501382046570     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -977.49931164091879     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.84
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1712.650       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.09E+00  8.55E-01  1.00E+00  9.80E-01  9.87E-01  9.20E-01  1.00E-02  9.35E-01  9.88E-01  1.04E+00
 


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
+        1.06E+03
 
 TH 2
+       -7.54E+00  4.63E+02
 
 TH 3
+        1.33E+01  2.14E+02  4.54E+02
 
 TH 4
+       -1.30E+01  4.25E+02 -1.72E+02  8.85E+02
 
 TH 5
+       -4.02E+00 -3.12E+02 -4.93E+02  1.69E+02  7.79E+02
 
 TH 6
+       -1.64E+00 -2.87E+00  6.85E+00 -7.04E-01 -3.08E+00  1.99E+02
 
 TH 7
+       -8.65E-01  2.03E+01  8.19E+00 -9.46E+00 -1.60E+01 -6.29E-01  4.16E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.68E+00 -2.94E+01 -2.78E+01  3.35E+01  2.83E+00 -3.39E+00  2.62E+01  0.00E+00  1.17E+02
 
 TH10
+       -1.00E+00 -3.28E+00 -4.58E+01 -1.31E+01 -5.05E+01  2.30E+00  1.48E+01  0.00E+00  1.05E+01  8.58E+01
 
 TH11
+       -8.38E+00 -2.10E+01 -4.29E+01 -7.93E+00  8.58E+00  3.09E+00  7.60E+00  0.00E+00  1.15E+01  1.94E+01  2.04E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.625
Stop Time:
Sat Sep 25 10:38:03 CDT 2021
