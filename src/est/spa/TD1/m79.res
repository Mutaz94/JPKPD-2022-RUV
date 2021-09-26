Sat Sep 25 13:07:14 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat79.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1719.38399821485        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -7.2731E+00 -6.9712E+01 -4.0561E+01 -4.4509E+01  3.0860E+01  2.5907E+01 -3.8051E+00  1.1115E+01  1.7302E+01  1.9886E+01
             4.8149E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1726.26540239457        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0296E+00  1.1437E+00  1.1391E+00  9.7007E-01  1.0816E+00  9.0465E-01  1.0690E+00  9.2232E-01  8.6748E-01  8.5381E-01
             1.0297E+00
 PARAMETER:  1.2922E-01  2.3424E-01  2.3023E-01  6.9609E-02  1.7845E-01 -2.0338E-04  1.6675E-01  1.9141E-02 -4.2162E-02 -5.8044E-02
             1.2925E-01
 GRADIENT:   6.8830E+01  2.7234E+01  1.5260E+01  1.9820E+01  9.4816E+00 -1.1228E+01 -5.5682E+00 -4.6482E+00 -5.7563E+00 -1.2604E+01
             1.1401E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1726.81071973123        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0254E+00  1.0669E+00  1.0374E+00  1.0074E+00  1.0058E+00  9.2132E-01  1.2098E+00  7.4454E-01  8.4234E-01  7.6687E-01
             9.9692E-01
 PARAMETER:  1.2513E-01  1.6473E-01  1.3670E-01  1.0736E-01  1.0576E-01  1.8047E-02  2.9045E-01 -1.9499E-01 -7.1571E-02 -1.6544E-01
             9.6918E-02
 GRADIENT:   5.9518E+01  2.0263E+01  1.8104E+01  1.5870E+01  1.2132E+01 -3.6634E+00  1.7475E+00 -4.9023E+00 -1.5164E+00 -1.8394E+01
            -1.2537E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1728.88869573371        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0081E+00  1.0528E+00  8.0957E-01  9.9204E-01  8.9858E-01  9.3477E-01  1.2157E+00  4.1636E-01  8.3539E-01  7.8497E-01
             9.9848E-01
 PARAMETER:  1.0811E-01  1.5142E-01 -1.1125E-01  9.2010E-02 -6.9388E-03  3.2542E-02  2.9535E-01 -7.7621E-01 -7.9852E-02 -1.4211E-01
             9.8475E-02
 GRADIENT:   3.5129E+00  1.1197E+00 -6.3724E+00  1.0380E+01  1.0984E+01  8.7241E-01  2.1246E+00  6.1101E-01  1.2265E+00  1.8566E-01
            -1.2005E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1729.47627302170        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      408
 NPARAMETR:  1.0322E+00  1.2011E+00  7.1332E-01  9.0187E-01  9.0573E-01  9.4279E-01  1.0844E+00  2.4366E-01  8.9246E-01  7.7374E-01
             1.0020E+00
 PARAMETER:  1.3165E-01  2.8326E-01 -2.3783E-01 -3.2813E-03  9.8399E-04  4.1087E-02  1.8100E-01 -1.3120E+00 -1.3769E-02 -1.5652E-01
             1.0204E-01
 GRADIENT:   1.2638E+01  7.6818E+00  5.3426E-01  5.3370E+00 -7.4492E+00  5.7805E-01  2.1729E-01  3.0618E-01 -6.3861E-01  1.5829E+00
            -5.2079E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1729.68551029875        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  1.0272E+00  1.3280E+00  6.7935E-01  8.1794E-01  9.5986E-01  9.4185E-01  9.9089E-01  1.1058E-01  9.7384E-01  7.9310E-01
             1.0035E+00
 PARAMETER:  1.2682E-01  3.8367E-01 -2.8662E-01 -1.0097E-01  5.9032E-02  4.0091E-02  9.0845E-02 -2.1020E+00  7.3492E-02 -1.3181E-01
             1.0345E-01
 GRADIENT:   1.5146E-01 -6.2654E-01 -2.7421E-01 -2.8534E-01  6.5209E-01 -2.5259E-02 -1.3710E-01  5.4656E-02  7.5967E-02 -5.9270E-02
             2.5666E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1729.71381562030        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      763
 NPARAMETR:  1.0264E+00  1.2976E+00  6.8649E-01  8.3705E-01  9.4707E-01  9.4146E-01  1.0109E+00  3.4023E-02  9.5680E-01  7.8681E-01
             1.0031E+00
 PARAMETER:  1.2609E-01  3.6048E-01 -2.7616E-01 -7.7871E-02  4.5618E-02  3.9674E-02  1.1089E-01 -3.2807E+00  5.5838E-02 -1.3976E-01
             1.0314E-01
 GRADIENT:  -1.5271E+00 -9.3802E-02  4.0840E-03  2.5982E-01  3.3526E-01 -1.3017E-01 -1.8828E-01  4.9903E-03  1.4596E-01 -2.1529E-01
             2.9862E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1729.71715543674        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  1.0270E+00  1.2927E+00  6.8715E-01  8.3994E-01  9.4477E-01  9.4173E-01  1.0154E+00  1.0000E-02  9.5265E-01  7.8688E-01
             1.0029E+00
 PARAMETER:  1.2666E-01  3.5675E-01 -2.7520E-01 -7.4428E-02  4.3191E-02  3.9958E-02  1.1525E-01 -4.8215E+00  5.1496E-02 -1.3969E-01
             1.0289E-01
 GRADIENT:  -9.3845E-02 -3.8635E-02 -2.1613E-02  3.8299E-03  4.4173E-02 -6.7022E-03 -5.7624E-03  0.0000E+00  2.2156E-03 -5.4033E-03
             2.0694E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1729.71715809968        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      995
 NPARAMETR:  1.0271E+00  1.2929E+00  6.8712E-01  8.3981E-01  9.4486E-01  9.4174E-01  1.0152E+00  1.0000E-02  9.5276E-01  7.8695E-01
             1.0029E+00
 PARAMETER:  1.2669E-01  3.5692E-01 -2.7524E-01 -7.4577E-02  4.3284E-02  3.9971E-02  1.1513E-01 -4.8249E+00  5.1613E-02 -1.3958E-01
             1.0289E-01
 GRADIENT:  -2.4473E-02 -1.5112E-02 -6.2745E-03 -1.8324E-03  1.4105E-02 -2.0757E-03 -1.8506E-03  0.0000E+00  7.1694E-04 -1.7846E-03
             7.6431E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      995
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.1805E-05 -9.4283E-03 -3.7646E-04  5.3454E-03 -1.9865E-02
 SE:             2.9806E-02  2.3589E-02  1.5753E-04  2.3136E-02  2.1709E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9888E-01  6.8939E-01  1.6855E-02  8.1728E-01  3.6015E-01

 ETASHRINKSD(%)  1.4547E-01  2.0973E+01  9.9472E+01  2.2493E+01  2.7273E+01
 ETASHRINKVR(%)  2.9073E-01  3.7548E+01  9.9997E+01  3.9926E+01  4.7108E+01
 EBVSHRINKSD(%)  5.0415E-01  2.0555E+01  9.9520E+01  2.3237E+01  2.6673E+01
 EBVSHRINKVR(%)  1.0058E+00  3.6885E+01  9.9998E+01  4.1075E+01  4.6231E+01
 RELATIVEINF(%)  9.8788E+01  3.8007E+00  2.3281E-04  3.6190E+00  5.6383E+00
 EPSSHRINKSD(%)  4.3289E+01
 EPSSHRINKVR(%)  6.7839E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1729.7171580996778     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -994.56633153593964     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1729.717       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.29E+00  6.87E-01  8.40E-01  9.45E-01  9.42E-01  1.02E+00  1.00E-02  9.53E-01  7.87E-01  1.00E+00
 


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
+        1.18E+03
 
 TH 2
+       -6.46E+00  4.20E+02
 
 TH 3
+        1.61E+01  2.03E+02  6.01E+02
 
 TH 4
+       -1.67E+01  3.36E+02 -3.85E+02  1.03E+03
 
 TH 5
+       -4.85E+00 -3.18E+02 -6.77E+02  4.31E+02  1.08E+03
 
 TH 6
+        1.01E+00 -1.01E+00  3.77E+00 -4.54E+00 -1.54E+00  2.19E+02
 
 TH 7
+        1.80E-01  2.75E+01 -1.91E+01 -1.42E+01  6.21E-01 -6.59E-01  8.15E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.50E-01 -2.17E+01 -3.55E+01  4.68E+01 -9.62E-01 -3.00E+00  2.10E+01  0.00E+00  8.72E+01
 
 TH10
+       -9.12E-01 -1.31E+01 -5.22E+01 -1.70E+01 -7.25E+01  1.36E+00  1.52E+01  0.00E+00  1.26E+01  9.96E+01
 
 TH11
+       -9.34E+00 -1.60E+01 -3.12E+01  3.16E+00  1.39E+00  3.96E+00  7.59E+00  0.00E+00  1.25E+01  2.73E+01  2.10E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.363
Stop Time:
Sat Sep 25 13:07:32 CDT 2021
